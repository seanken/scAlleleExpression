##
##Takes in a table of SNP counts produced by LoadSNPs.R and tests for AI SNPs 
##For COndition: form=Count~Allele+Sample+Condition*Allele,coef_ret="Alleleref:Conditionyes"
##
#'Test SNP for allelic imbalance with BetaBinomial Regression
#'
#'Combines allele expression data for numerous individuals to see if there is significant allelic imbalance at that SNP.
#'Based on a betabinomial model
#'
#' @param dat A dataframe with columns corresponding to the expression of each SNP. Should be a column called Feature (the SNP of interest), Gene (the gene being tested), Sample (which individual it came from), and Allele (alt or ref). The naive model just tests for presence of allelic imbalance at the SNP, but can test for conditional as well
#' @param minCount The min number of reads over all samples for the SNP to be tested
#' @param minSamp The min number of samples with allelic information about this gene for it to be used
#' @param form The formula for the betabinomial model being fitted. The default works in most cases in a format as used by aod 
#' @param coef_ret The coefficient to be tested (in most cases the intercept)
#' @param numTest Used mostly for debugging, will only test this many SNPs and ignore the rest
#' @return A dataframe with the results of the allelic imbalance test
#' @export
TestSNP_aod_Sierra<-function(dat,minCount=50,minSamp=10,form=cbind(Count, Tot - Count) ~ Allele + Sample,form2=cbind(alt, Num - alt) ~ 0,coef_ret="Allele",numTest=-1)
{
    print("Format")
    dat<-dat %>% unite(Feature,Gene,SNP,sep="_",remove=F)
    tab<-dat %>% group_by(Feature,Sample) %>% summarise(Num=length(unique(Allele)),Count=sum(Count)) %>% group_by(Feature) %>% summarise(Count=sum(Count),NumSamp=sum(Num>1)) %>% as.data.frame()
    print(length(unique(dat$Gene)))
    feats=tab[tab$NumSamp>minSamp & tab$Count>minCount,"Feature"]
    dat=dat[dat$Feature %in% feats,]
    print("Number to Test:")
    print(length(feats))
    print(length(unique(dat$Gene)))
    dat["GeneName"]=as.character(lapply(dat[,"Gene"],function(x){strsplit(x,"_")[[1]][1]}))
    tab<-dat %>% group_by(Allele,GeneName,Sample) %>% summarise(Tot=sum(Count)) %>% as.data.frame()
    dat<-dat %>% spread(Allele,Count,fill=0) %>% gather(Allele,Count,alt,ref)
    dat=inner_join(dat,tab) ##might modify so include those with 0 in peak but some in gene
    print(head(dat))
    print("Split by SNP")
    #return(dat)########
    #bySNP=lapply(feats,function(x){dat[dat[,1]==x,]})
    #names(bySNP)=feats
    bySNP=split(dat,f=dat$Feature)
    if(numTest>0){bySNP=bySNP[sample(1:length(bySNP),numTest)]}
    print("Test")
    tic()
    nams=names(bySNP)
    out=lapply(names(bySNP),function(cur_feat){
        x=bySNP[[cur_feat]]
        fit=tryCatch({betabin(form,~1,data=x)},error=function(cond){return(NULL)})
        if(is.null(fit))
        {
            return(NULL)
        }


        coef=tryCatch({summaryAOD(fit)@Coef},error=function(cond){return(NULL)})
        if(is.null(coef))
        {
            return(NULL)
        }
        coef=data.frame(coef)
        colnames(coef)[4]="pval"
        coef["Test"]=rownames(coef)
        coef=coef[grep(coef_ret,rownames(coef)),]
        coef["SNP"]=cur_feat
        coef["NumSamp"]=length(unique(x$Sample))
        return(coef)

    })
    toc()
    bySNP=out
    names(bySNP)=nams
    bySNP[sapply(bySNP, is.null)]=NULL

    ret=do.call(rbind,bySNP)

    ret=data.frame(ret)

    ret=ret[order(ret$pval),]
    ret["logP"]=log(2)+pnorm(abs(ret[,"z.value"]),lower.tail=FALSE, log.p=TRUE)
    ret["pval"]=exp(ret[,"logP"])
    ret["padj"]=p.adjust(ret[,"pval"],"fdr")
    ret["Gene"]=as.character(lapply(ret[,"SNP"],function(x){strsplit(x,split="_")[[1]][1]}))
    rownames(ret)=NULL
    return(ret)

}
