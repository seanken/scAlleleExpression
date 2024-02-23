
library(purrr)
library(tidyr)
library(dplyr)

#' Run DAESC
#'
#' A method that takes in the count data produced by this package and uses DAESC with it.
#' 
#' @param dat The ASE count data produced by the GetSNPs command
#' @param minCells The min number of cells with at least one phased UMI in a SNP required to test a SNP
#' @param minSamps The min number of samples with at least one phased UMI in a SNP required to test a SNP
#' @param form The formula for the model to use in LRT
#' @param formNull The formula for the null model to compare to
#' @param method The name of the DAESC method to use, either bb or mix (if anything else is passed will use mix)
#' @return A data frame with the results for the DAESC method
#' @export

RunDAESC<-function(dat,minCells=50,minSamps=10,form=~Condition,formNull=~1,method="bb")
{
    library(DAESC)
    dat<-dat %>% unite(Test,SNP,Gene,sep="_",remove=F) %>% spread(Allele,Count,fill=0)
    counts=split(dat,dat$Test)
    counts=counts[map_dbl(counts,function(x) dim(x)[1])>minCells]
    counts=counts[map_dbl(counts,function(x) length(unique(x$Sample)))>minSamps]
    results=map(counts,function(df){
        print(head(df))
        print(dim(df))
        print(dim(cbind(1,df$Condition)))
        print(dim(matrix(1,nrow=nrow(df),ncol=1)))
        modMat=model.matrix(form,df)
        modMatNull=model.matrix(formNull,df)
        daesc=if(method=="bb") daesc_bb else daesc_mix

        res <- tryCatch({daesc(y=df[,"alt"], n=(df[,"alt"]+df[,"ref"]), subj=df$Sample, x=modMat, xnull=modMatNull, niter=200, niter_laplace=2, num.nodes=3,optim.method="BFGS", converge_tol=1e-8)},error=function(e){print("Yuck!");return(NULL)})
        print("Ran!")
        res
    })

    results=results %>% discard(is.null)
    results=map(results,function(tab) tab[map_dbl(tab,function(x){ret=length(dim(x));if(is.null(ret)){ret=0};ret})==0])
    results=map(results,function(x) x[names(x)!="p"])

    results=map(names(results),function(x){
        df=data.frame(results[[x]])
        df["Coef"]=sub("^","beta_",rownames(df))
        df["Test"]=x
        df<-df %>% spread(Coef,b)
        return(df)
    })

    mrk=do.call(rbind,results)
    mrk=mrk[order(mrk$p.value),]
    mrk["FDR"]=p.adjust(mrk[,"p.value"],"fdr")
    mrk["FWER"]=p.adjust(mrk[,"p.value"])
    return(mrk)
}