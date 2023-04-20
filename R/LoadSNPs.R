


#'Load SNP Data
#'
#'Loads SNP data for genes of interest.
#'
#' @param seur A seurat object prepared to have the correct columns (see below). Assumes the cell names in Seurat are [samp name]_[CBC] or CBC is given as column
#' @param meta A dataframe with meta data, can be given instead of Seurat object
#' @param genes A list of genes to get SNP specific information for
#' @param snp_list List of SNPs, if given will only return snps from the list, must be the same length as genes
#' @param snp_col The name of a column in seur@meta.datapointing to the SNPLevelCounts/snps.bed output from the AI pipeline for this cell
#' @param all_col The name of a column in seur@meta.datapointing to the AlleleCounts/counts.txt output from the AI pipeline for this cell
#' @param cbc_col Location of cell barcode in meta data, if not given assums sample names are [samp name]_[CBC]
#' @param cond The column in seur@meta.data with condition ifnormation
#' @param countCells If True returns cell counts instead of read counts
#' @param cellType The name of the column in seur@meta.data with cell type information if required
#' @param bulk If true makes pseudobulk, else returns at cell level
#' @return A data frame with bulk allele expression information at the SNP level
#' @export
GetSNPs<-function(seur=NULL,meta=NULL,genes,snp_list=c(),samp="orig.ident",snp_col="SNP",all_col="Allele",cbc_col="",cond="PD",countCells=F,cellType="",bulk=T)
{
    if(!is.null(seur))
    {
        meta=seur@meta.data
    }

    if(nchar(cbc_col)>2)
    {
        rownames(meta)=apply(meta,1,function(x){paste(x[samp],x[cbc_col],sep="_")})
    }

    print("Make Nice")
    tab=meta[,c(samp,snp_col,all_col,cond)]
    tab["CellType"]="None"
    if(nchar(cellType)>1)
    {
        tab=meta[,c(samp,snp_col,all_col,cond,cellType)]
    }
    colnames(tab)=c("Samp","SNP","Allele","Condition","CellType")
    tab["Name"]=rownames(meta)
    dat<-tab %>% group_by(Samp,SNP,Allele,Condition) %>% summarise() %>% as.data.frame()

    celltypes=tab[,"CellType"]
    names(celltypes)=tab[,"Name"]

    print("Get SNP data")
    out=lapply(1:dim(dat)[1],function(i){
        print(i)

        #print(head(celltypes))
        mat=tryCatch({mat=loadSNPs2(dat[i,1],dat[i,2],dat[i,3],tab[tab[,"Samp"]==dat[i,1] & tab[,"SNP"]==dat[i,2] & tab[,"Allele"]==dat[i,3],"Name"],celltypes=celltypes,genes=genes,countCells=countCells,snp_list=snp_list,bulk=bulk)
            mat["Condition"]=dat[i,"Condition"]
            #print(head(mat))
            mat
        },error=function(cond){return(NULL)})
        #print(dim(mat))
        return(mat)
    })

    out[sapply(out, is.null)]=NULL
    dat=do.call(rbind,out)

    dat=dat[!is.na(dat[,"Condition"]),]

    return(dat)

}


##
##Gets pseudobulk allele specific and SNP specific counts
##Takes in:
##Sample name (samp)
##snp comb.bed file location (snp)
##Allele gene level fount data (allele)
##list of cell names to include (nams)
##list of genes to use (genes)
##
##Returns Psuedobulk
#' Load SNP data
#'
#' Helper function to load SNP level ASE data (new version). Either as pseudobulk or at cell level.
#'
#' @param samp Sample name
#' @param snp Location of snps.bed file
#' @param allele Gene level ASE count data
#' @param nams Cell names to use
#' @param genes List of genes to load SNP data for
#' @param celltypes Cell type information
#' @param getVals For testing
#' @param countCells For testing
#' @param bulk If true makes pseudobulk, else returns at cell level
#' @param snp_list A list of SNPs to look at, must be the same length as genes
#' @return A data from with SNP level ASE information
#' @export
loadSNPs2=function(samp,snp,allele,nams,genes,celltypes,getVals=F,countCells=F,snp_list=c(),bulk=T)
{
    if(length(snp_list)!=length(genes))
    {
        print("Need same number of SNPs and Genes!")
        return(NULL)
    }

    print("Load Allele Data")
    cnts=read.table(allele,stringsAsFactors=F)
    colnames(cnts)=c("CBC","Gene","Allele","Count")
    cnts=cnts[cnts$Allele!="Ambig",]
    if(countCells)
    {
        print("Cell Counts")
        cnts_tab=cnts %>% group_by(CBC,Gene) %>% summarise(Allele=Allele[which.max(Count)],Max=max(Count),NumMax=sum(Count==Max)) %>% as.data.frame()
        cnts_tab=cnts_tab[cnts_tab$NumMax==1,]
        cnts_tab["Count"]=1
        cnts=cnts_tab[,c("CBC","Gene","Allele","Count")]
    }

    print("Process")
    cnts["Name"]=sub("^",paste(samp,"_",sep=""),cnts[,"CBC"])
    cnts=cnts[cnts$Name %in% nams,]

    print(head(cnts))

    cnts["CellType"]=celltypes[as.character(cnts$Name)]
    print(dim(cnts))
    print(head(cnts))
    if(bulk)
    {
        print("Make PseudoBulk")
        cnts<-cnts %>% group_by(Gene,Allele,CellType) %>% summarise(Count=sum(Count)) %>% as.data.frame()
    }
    cnts["Sample"]=samp
    print(head(cnts))
    print("Load SNP info")
    bed=read.table(snp,stringsAsFactors=F)
    bed=bed[,4:5]
    colnames(bed)=c("SNP","Geno")
    bed=bed[bed$SNP %in% snp_list,]
    bed <- bed %>% group_by(SNP,Geno) %>% summarise() %>% as.data.frame()
    tmp=data.frame(cbind(snp_list,genes))
    colnames(tmp)=c("SNP","Gene")
    tmp<-tmp %>% group_by(SNP,Gene) %>% summarise() %>% as.data.frame()
    bed=inner_join(bed,tmp)
    bed=bed[,c("SNP","Gene","Geno")]
    bed=bed[bed$Gene %in% genes,]
    if(length(snp_list)>1)
    {
        bed=bed[bed$SNP %in% snp_list,]
    }
    print(head(bed))
    print(length(unique(bed$Gene)))
    print(length(unique(cnts$Gene)))
    print(length(intersect(unique(bed$Gene),unique(cnts$Gene))))
    if(getVals)
    {
        lst=list()
        lst[[1]]=cnts
        lst[[2]]=bed
        return(lst)
    }
    print("Combine with Allele data")
    cnts=inner_join(cnts,bed)
    print(length(unique(cnts$Gene)))
    lst=c("alt","ref","ref","alt")
    names(lst)=c("All1_1|0","All1_0|1","All2_1|0","All2_0|1")
    cnts<-cnts %>% unite(Comb,Allele,Geno,sep="_")
    cnts=cnts[cnts$Comb %in% names(lst),]
    cnts["Allele"]=lst[cnts[,"Comb"]]
    if(bulk)
    {
        cnts=cnts[,c("Gene","Sample","SNP","CellType","Allele","Count")]
    }
    if(!bulk)
    {
        cnts=cnts[,c("Gene","Sample","SNP","CellType","Allele","Count","CBC")]
    }

    print(dim(cnts))
    print(head(cnts))
    print(length(unique(cnts$Gene)))
    print("Return!")
    return(cnts)
}
