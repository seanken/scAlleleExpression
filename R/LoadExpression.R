#' Loads expression data based on meta data
#' 
#' Load the expression data from the pipeline.
#' 
#' @param meta The metadata with pipeline info as produced by LoadMeta
#' @param cbc_col Name of column with cell barcode info
#' @param raw_col Which column contains raw count data
#' @param samp_col Name of the column with sample of origin info
#' @return A sparse count matrix corresponding to meta
#' @export 
LoadExpression<-function(meta,cbc_col="CBC",raw_col="Raw",samp_col="orig.ident",genes_include=NULL)
{   
    samps=unique(meta[,raw_col])
    rownames(meta)=apply(meta,1,function(x){paste(x[samp_col],x[cbc_col],sep="_")})
    out=lapply(samps,function(x){
        print(x);
        dat=LoadCountMatrix(x)
        meta_tmp=meta[meta[,raw_col]==x,]
        cbc=meta_tmp[,cbc_col]
        nams=rownames(meta_tmp)
        names(nams)=cbc
        print("Number in meta")
        print(dim(meta_tmp))
        print("Number in dat before:")
        print(dim(dat))
        print(head(meta_tmp))
        print(head(colnames(dat)))
        dat=dat[,cbc]
        print("Number in dat after:")
        print(dim(dat))

        colnames(dat)=as.character(nams[colnames(dat)])
        if(!is.null(genes_include))
        {
            genes=intersect(rownames(dat),genes_include)
            dat=dat[genes,]
        }
        return(dat)
    })
    dat=do.call(cbind,out)
    dat=dat[,rownames(meta)]
    return(dat)
}

#' Load expression
#' 
#' Loads MM data as produced by STARSolo. Gene names are the gene name and gene id field combined.
#' 
#' @param dirRaw The directory with expression info
#' @return Sparse count matrix
#' @export 
LoadCountMatrix<-function(dirRaw)
{
    dat=readMM(paste(dirRaw,"/matrix.mtx",sep=""))
    cbc=scan(paste(dirRaw,"/barcodes.tsv",sep=""),"")
    colnames(dat)=cbc
    genes=read.table(paste(dirRaw,"/features.tsv",sep=""),sep="\t")
    colnames(genes)=c("ID","Name","Type")
    genes<-genes %>% unite(Row,Name,ID,sep="_")
    rownames(dat)=genes[,"Row"]
    return(dat)

}