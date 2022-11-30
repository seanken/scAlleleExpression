#' Loads expression data based on meta data
#' 
#' Load the expression data from the pipeline.
#' 
#' @param meta The metadata with pipeline info as produced by LoadMeta
#' @param cbc_col Name of column with cell barcode info
#' @param raw_col Which column contains raw count data
#' @return A sparse count matrix corresponding to meta
#' @export 
LoadExpression<-function(meta,cbc_col="CBC",raw_col="Raw")
{   
    samps=unique(meta[,raw_col])
    out=lapply(samps,function(x){
        dat=LoadCountMatrix(x)
        meta_tmp=meta[meta[,raw_col]==x,]
        cbcs=meta_tmp[,cbc_col]
        nams=rownames(meta_tmp)
        names(nams)=cbc
        dat=dat[,cbc]
        colnames(dat)=nams[colnames(dat)]
        return(dat)
    })
    dat=do.call(rbind,out)
    dat=dat[,rownames(meta)]
    return(dat)
}

#' Load expression
#' 
#' Loads MM data as produced by STARSolo
#' 
#' @param dirRaw The directory with expression info
#' @return Sparse count matrix
#' @export 
LoadCountMatrix<-function(dirRaw)
{
    dat=readMM(paste(dirRaw,"",sep=""))
    cbc=scan(paste(dirRaw,"",sep=""),"")
    colnames(dat)=cbc
    genes=read.table(paste(dirRaw,"",sep=""),sep="\t")
    colnames(genes)=c("Name","ID","Type")
    genes<-genes %>% unite(Row,Name,ID,sep="_")
    rownames(dat)=genes[,"Row"]
    return(dat)

}