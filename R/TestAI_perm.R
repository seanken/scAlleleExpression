##
#'Test SNP for allelic imbalance with BetaBinomial Regression and permuted p-values
#'
#'Combines allele expression data for numerous individuals to see if there is significant allelic imbalance at that SNP.
#'Based on a betabinomial model with permutations to get p-values
#'
#' @param dat A dataframe with columns corresponding to the expression of each SNP. Should be a column called Feature (the SNP of interest), Gene (the gene being tested), Sample (which individual it came from), SNP (if SNPLevel=T), and Allele (alt or ref). The naive model just tests for presence of allelic imbalance at the SNP, but can test for conditional as well. If SNPLevel=T suggest using pseudobulked data, if =F use cell level data from one individual.
#' @param nReps Number of reps to run to get p-value
#' @param SNPLevel A binary, true by default. If false does gene level analysis. If true minSamp is set to 0.
#' @param minCount The min number of reads over all samples for the SNP to be tested
#' @param minSamp The min number of samples with allelic information about this gene for it to be used
#' @param form The formula for the betabinomial model being fitted. The default works in most cases in a format as used by aod 
#' @param coef_ret The coefficient to be tested (in most cases the intercept)
#' @param numTest Used mostly for debugging, will only test this many SNPs and ignore the rest
#' @return A dataframe with the results of the allelic imbalance test
#' @export
TestSNP_perm<-function(dat,nRep=1000,minCount=50,minSamp=10,form=cbind(alt, Num - alt) ~ 1,coef_ret="(Intercept)",numTest=-1)
{
    print("Real test")
    mrk=TestSNP_aod(dat,minCount=minCount,minSamp=minSamp,form=form,coef_ret=coef_ret,numTest=numTest,correctMethod ="fdr")
    mrk["NumTested"]=0
    mrk["NumMoreSig"]=0
    rownames(mrk)=mrk[,"SNP"]
    print(head(mrk))
    for(x in 1:nRep)
    {
        print("Next Rep:")
        print(x)
        tab=PermDat(dat,x);
        res=TestSNP_aod(tab,minCount=minCount,minSamp=minSamp,form=form,coef_ret=coef_ret,numTest=numTest,correctMethod ="fdr")
        rownames(res)=res[,"SNP"]
        print(head(res))
        snps=intersect(rownames(res),rownames(mrk))
        res=res[snps,]
        mrk[snps,"NumTested"]=mrk[snps,"NumTested"]+1
        snps=snps[res[snps,"pval"]<=mrk[snps,"pval"]]
        snps=snps[!is.na(snps)]
        mrk[snps,"NumMoreSig"]=mrk[snps,"NumMoreSig"]+1
        print(head(mrk))
    }
    mrk["pval_perm"]=(mrk[,"NumMoreSig"]+1)/(mrk[,"NumTested"]+1)
    mrk["padj_perm"]=p.adjust(mrk[,"pval_perm"],"fdr")
    return(mrk)

}


#'Permutes ASE data
#' 
#' Randomly permutes the ref vs alt allele in ASE data 
#' 
#' @param dat A dataframe with columns corresponding to the expression of each SNP. Should be a column called Feature (the SNP of interest), Gene (the gene being tested), Sample (which individual it came from), SNP (if SNPLevel=T), and Allele (alt or ref). The naive model just tests for presence of allelic imbalance at the SNP, but can test for conditional as well. If SNPLevel=T suggest using pseudobulked data, if =F use cell level data from one individual.
#' @return A dataframe with the results of the allele permuted
#' @export
PermDat<-function(dat,seed=1)
{
    set.seed(seed)
    dat=dat[dat$Allele %in% c("ref","alt"),]
    samps=unique(dat[,"Sample"])
    samps=samps[sample(c(0,1),length(samps),replace=T)==1]
    dct=c("alt","ref")
    names(dct)=c("ref","alt")
    dat[dat[,"Sample"] %in% samps,"Allele"]=as.character(dct[dat[dat[,"Sample"] %in% samps,"Allele"]])
    return(dat)
}
