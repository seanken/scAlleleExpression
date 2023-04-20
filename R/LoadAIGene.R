#'Load Gene ASE Data
#'
#'Loads Gene level ASE data 
#'
#' @param seur A seurat object prepared to have the correct columns (see below). Assumes the cell names in Seurat are [samp name]_[CBC] or CBC is given as column
#' @param meta A dataframe with meta data, can be given instead of Seurat object
#' @param all_col The name of a column in seur@meta.datapointing to the AlleleCounts/counts.txt output from the AI pipeline for this cell
#' @param cbc_col Location of cell barcode in meta data, if not given assums sample names are [samp name]_[CBC]
#' @param cond The column in seur@meta.data with condition ifnormation
#' @param cellType The name of the column in seur@meta.data with cell type information if required
#' @return A data frame with allele expression information at the gene level
#' @export
loadASEGene<-function(seur=NULL,meta=NULL,samp="orig.ident",all_col="Allele",cbc_col="CBC",cond="PD",celltype="CellType")
{
    print("Prep to load data")
    dat=meta[,c(samp,all_col,cond)]
    colnames(dat)=c("Samp","Allele","Condition")
    dat <- dat %>% group_by(Samp,Allele,Condition) %>% summarise() %>% as.data.frame()
    out=lapply(1:dim(dat)[1],function(i){
        print(i)
        tab=read.table(dat[i,"Allele"])
        tab["Name"]=sub("^",paste(dat[i,"Samp"],"_",sep=""),tab[,1])
        meta2=meta[meta[,samp]==dat[i,"Samp"],]
        rownames(meta2)=meta2[,cbc_col]
        print(dim(tab))
        print(head(tab))
        #print(head(meta2))
        tab=tab[tab[,1] %in% rownames(meta2),]
        print(dim(tab))
        colnames(tab)[2]="Gene"
        tab["Sample"]=dat[i,"Samp"]
        colnames(tab)[3]="Allele"
        colnames(tab)[4]="Count"
        tab["SNP"]=tab["Gene"]
        tab[tab$Allele=="All1","Allele"]="alt"
        tab[tab$Allele=="All2","Allele"]="ref"
        tab["Cond"]=dat[i,"Condition"]
        tab["CellType"]=meta2[as.character(tab[,1]),celltype]
        return(tab)
    })
    return(out)
    tab=do.call(rbind,out)
    return(tab)

}
