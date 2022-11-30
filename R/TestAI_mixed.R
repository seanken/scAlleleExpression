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
TestSNP_mixed<-function(dat,minCount=50,minSamp=10,form=cbind(alt, Num - alt) ~ 1 + (1|Sample),form2=cbind(alt, Num - alt) ~ 0,coef_ret="(Intercept)",numTest=-1)
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



    print("Split by SNP")
    bySNP=split(dat,f=dat$Feature)
    if(numTest>0){bySNP=bySNP[1:numTest]}
    print("Test")
    tic()
    nams=names(bySNP)
    out=lapply(names(bySNP),function(cur_feat){
        x=bySNP[[cur_feat]]
        x=x %>% spread(Allele,Count,fill=0)
        x["Num"]=x["alt"]+x["ref"]
        fit=tryCatch({glmmTMB(form,family=list(family="betabinomial",link="logit"),data=x)},error=function(cond){return(NULL)})
        if(is.null(fit))
        {
            return(NULL)
        }


        coef=coef(summary(fit))$cond
        coef=data.frame(coef)
        coef["Test"]=rownames(coef)
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

    colnames(ret)[4]="pval"

    ret=ret[order(ret$pval),]
    ret["padj"]=p.adjust(ret[,"pval"])
    ret["Gene"]=as.character(lapply(ret[,"SNP"],function(x){strsplit(x,split="_")[[1]][1]}))
    rownames(ret)=NULL
    return(ret)

}
