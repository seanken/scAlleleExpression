
#' Add info about pipeline output to meta data
#' 
#' Given a data frame containing meta data on a per cell basis and names list mapping sample to pipeline output directory, adds the pipeline info to the meta data
#' @param meta A data frame with cell level info
#' @param lst An named character array, one entry for each sample, named by sample (should agree with the naming in meta)
#' @param samp_col The name of the column in meta with info about sample of origin
#' @param cbc_col The name of the column in meta with info about cbc (can be excluded for now)
#' @return The meta dataframe with the added columns about the pipeline added
#' @export
LoadPipeline<-function(meta,lst,samp_col="orig.ident",cbc_col="CBC")
{
    cols=c("dir_ASE","Allele","SNP","Raw")
    if(length(intersect(colnames(meta),cols))>0)
    {
        print("Error, can not have columns named:");
        print(cols);
        print("In input meta data")
        return(meta);
    }    

    if(length(lst)!=length(unique(meta[,samp_col])) | length(intersect(names(lst),unique(meta[,samp_col])))!=length(lst) )
    {
        print("Sample names do not agree between lst and meta data")
        return(meta);
    }

    meta["dirASE"]=lst[meta[,samp_col]]
    meta["Allele"]=sub("$","/output/AlleleCounts/counts.txt",meta[,"dirASE"])
    meta["SNP"]=sub("$","/output/SNPLevelCounts/snps.bed",meta[,"dirASE"])
    meta["Raw"]=sub("$","/output/STARSolo/output/resultsSolo.out/GeneFull/raw",meta[,"dirASE"])

    return(meta)


}
