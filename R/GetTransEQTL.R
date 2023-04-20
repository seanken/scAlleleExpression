
#' Get trans eqtl
#' 
#' A tool that uses variation in ASE between cells to help estimate trans-eQTL for a given gene/SNP pair.
#' 
#' @param meta The meta data as produced by LoadPipeline
#' @param gene The gene whose ASE is being tested 
#' @param snp The snp used to align alleles between samples
#' @param tab The SNP level allele specific information, calculated internally if not given
#' @param dat The UMI expression data, calculated internally if not given
#' @param multIndiv Set to true if multiple individuals, false if only onesystem()
#' @param DSCells DS number of cells
#' @return A table with the results of the trans-eqtl testing
#' @export 
GetTransEQTL<-function(meta,gene,snp,tab=NULL,dat=NULL,samp_col="orig.ident",snp_col="SNP",raw_col="Raw",all_col="Allele",cbc_col="CBC",cond="PD",on_disk=T,method="nb",multIndiv=T,form=NULL,DSCells=-1,pseudo=1)
{
    print("Load ASE information")
    if(is.null(tab))
    {
        tab=GetSNPs(meta=meta,genes=c(gene),snp_list=c(snp),samp=samp_col,snp_col=snp_col,all_col=all_col,cbc_col=cbc_col,bulk=F,cond=cond)
    }
    tab=tab[tab$Gene==gene & tab$SNP==snp,]
    tab<-tab %>% unite(Name,Sample,CBC,sep="_",remove=F) %>% spread(Allele,Count,fill=0)
    if(!("alt" %in% colnames(tab))){tab["alt"]=0}
    if(!("ref" %in% colnames(tab))){tab["ref"]=0}
    tab=tab[,c("alt","ref","CBC","Sample","Name")]
    

    print("Combine ASE and meta")
    meta=data.frame(Raw=meta[,raw_col],Sample=meta[,samp_col],CBC=meta[,cbc_col])
    #if(!is.null(dat)){colnames(dat)=tab[,"Name"]}
    meta=inner_join(meta,tab)
    rownames(meta)=meta[,"Name"]
    samps=table(meta$Sample)
    print(dim(meta))
    print(samps)
    samps=names(samps)[samps>10]
    meta=meta[meta$Sample %in% samps,]
    print(dim(meta))
    if(is.null(dim(meta))){
        print("Not enough cells!")
        return(NULL)
    }
    if(sum(as.numeric(table(meta$Sample))>10)<3){
        print("Not enough samples!")
        return(NULL)
    }
    if(!is.null(dat)){dat=dat[,rownames(meta)]}

    if(is.null(dat))
    {
        print("Load expression data for remaining cells")
        dat=LoadExpression(meta,cbc_col="CBC",raw_col="Raw",samp_col="Sample")
    }
    
    print("Get meta data and data to agree")
    inter=intersect(colnames(dat),rownames(meta))
    if(length(inter)<10)
    {
        print("Less than 10 cells shared between data and meta! Killing program.")
        return(NULL);
    }
    if(DSCells>0 & DSCells<length(inter))
    {
        inter=sample(inter,DSCells)
        print("DS!")
        print(length(inter))
    }
    dat=dat[,inter]
    meta=meta[inter,]
    meta["Gene"]=dat[grep(gene,rownames(dat))[1],]

    print("Run testing")
    mrk=TestTransEQTL(meta,dat,samp_col="Sample",on_disk=on_disk,multIndiv=multIndiv,form=form,method=method,pseudo=pseudo)
    mrk["Gene_ASE"]=gene
    mrk["SNP_ASE"]=snp
    print("Done!")
    return(mrk)
}

#'Test trans-eQTL
#' 
#' Given ASE data and count data tests for transeQTL
#' 
#' @param meta Meta data frame including ASE data
#' @param dat Matrix (can be sparse) of expression counts, genes by cells
#' @param samp_col The column with sample of origin
#' @param pseudo Pseudocount to use (corresponds to prior in beta prior)
#' @param percCut Percent of cells that need to express a gene for it to be tested
TestTransEQTL<-function(meta,dat,samp_col="orig.ident",pseudo=1,percCut=.05,on_disk=T,method="nb",multIndiv=T,form=NULL)
{
    meta_loc=data.frame(Gene=meta[,"Gene"],rat=(meta[,"ref"]+pseudo)/(meta[,"alt"]+meta[,"ref"]+2*pseudo),indiv=meta[,samp_col],alt=meta[,"alt"],ref=meta[,"ref"])
    meta_loc["UMI"]=Matrix::colSums(dat)
    meta_loc["Gene"]=log(10000*meta_loc[,"Gene"]/meta_loc[,"UMI"]+1)
    rownames(meta_loc)=rownames(meta)
    #nUMI=colSums(dat)
    mn=Matrix::rowMeans(dat>0)
    dat=dat[mn>percCut,]

    mrk=NULL
    print("Run!")
    if(method=="nb")
    {
        mrk=RunNB(dat,meta_loc,on_disk,multIndiv=multIndiv,form=form)
    }
    if(method=="betabinom")
    {
        mrk=RunBetaBin(dat,meta_loc,multIndiv=multIndiv)
    }
    if(method=="betabinom_mix")
    {
        print("Mixed Beta BInomial model not yet implemented, will use the fixed effect model instead")
        mrk=RunBetaBin(dat,meta_loc,multIndiv=multIndiv)
        #mrk=RunBetaBinMix(dat,meta_loc)
    }
    print("Return results")
    return(mrk)

}


#' Run BetaBin
#' 
#' @export 
RunBetaBin=function(dat,meta,multIndiv=T)
{
    print("betabinomial!")
    genes=rownames(dat)
    #genes=c("CNR1_ENSG00000118432")
    #for(i in colnames(dat)){dat[i]=log(10000*dat[,i]/sum(dat[,i])+1)} ##normalize
    meta["nUMI"]=Matrix::colSums(dat)
    print("Run per gene")
    out=lapply(genes,function(x){
        ##print(x);
        meta["Gene"]=log(10000*dat[x,]/meta[,"nUMI"]+1)
        form=cbind(alt, ref) ~ indiv+Gene 
        if(!multIndiv)
        {
            form=cbind(alt, ref)~ Gene
        }
        
        fit=tryCatch({betabin(form,~1,data=meta)},error=function(cond){return(NULL)})
        if(is.null(fit))
        {
            return(NULL)
        }
        
        coef=tryCatch({summaryAOD(fit)@Coef},error=function(cond){return(NULL)})
        if(is.null(coef)){return(NULL)}
        coef=data.frame(coef)
        
        coef=coef[grep("Gene",rownames(coef)),]
        
        colnames(coef)[4]="pval"
        coef["logP"]=log(2)+pnorm(abs(coef[,"z.value"]),lower.tail=FALSE, log.p=TRUE)
        coef["pval"]=exp(coef[,"logP"])
        coef["Test"]=rownames(coef)
        coef["GeneTested"]=x;
        coef["NumSamp"]=length(unique(meta$indiv))
        return(coef)
    })
    mrk=do.call(rbind,out)
    mrk=mrk[order(mrk$pval),]
    mrk["FDR"]=p.adjust(mrk$pval,"fdr")
    return(mrk)
}

RunBetaBinMix=function(dat,meta_loc)
{
    return(NULL)
}

RunNB=function(dat,meta_loc,on_disk,multIndiv=T,form=NULL)
{
    print("Fit model")
    if(on_disk)
    {
        genes=rownames(dat)
        #cells=colnames(dat)
        dat <- HDF5Array::writeHDF5Array(dat)
        #rownames(dat)=genes
        #colnames(dat)=cells
    }
    if(is.null(form))
    {
        form=~rat+indiv
        if(!multIndiv)
        {
            form=rat
        }
    }
    
    fit=glm_gp(data=dat,col_data=meta_loc,design=form,on_disk = on_disk)
    
    print("Run DE")
    mrk=test_de(fit,"rat")
    mrk["gene"]=mrk[,1]

    if(on_disk)
    {
        mrk["gene"]=genes[as.numeric(sub("row_","",mrk[,1]))]
    }
    

    
    mrk=mrk[order(mrk$pval),]
    mrk["gene_name"]=as.character(lapply(mrk[,1],function(x){strsplit(x,"_")[[1]][1]}))
    mrk["gene_id"]=as.character(lapply(mrk[,1],function(x){strsplit(x,"_")[[1]][2]}))

    return(mrk)

}

#' Get gene level trans-eQTL analysis
#' 
#' Combines the results of the trans-eQTL analysis to see if there are genes that are overrepresented among the top hits.
#' 
#' @param mrk A data frame generated by combining the results from GetTransEQTL for multiple SNP/genes
#' @return A table with pvalues and FDR at the genes level
#' @export
GetGeneLevel<-function(mrk)
{
    tab<-mrk %>% group_by(gene) %>% summarise(num=sum(pval<.05),tot=length(pval)) %>% filter(tot>20) %>% as.data.frame()
    tab["pval_binom"]=apply(tab,1,function(x){num=as.numeric(x["num"]);tot=as.numeric(x["tot"]);binom.test(num,tot,p=.05,alternative="greater")$p.value})
    tab=tab[order(tab$pval_binom),]
    tab["FDR_binom"]=p.adjust(tab$pval_binom,"fdr")
    return(tab)
}
