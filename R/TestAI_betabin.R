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
TestSNP_betabin<-function(dat,minCount=50,minSamp=10,form=cbind(alt, Num - alt) ~ 1,form2=cbind(alt, Num - alt) ~ 0,coef_ret="(Intercept)",numTest=-1)
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
#return(dat)########
#bySNP=lapply(feats,function(x){dat[dat[,1]==x,]})
#names(bySNP)=feats
bySNP=split(dat,f=dat$Feature)
if(numTest>0){bySNP=bySNP[1:numTest]}
print("Test")
tic()
nams=names(bySNP)
out=lapply(names(bySNP),function(cur_feat){
x=bySNP[[cur_feat]]
x=x %>% spread(Allele,Count,fill=0)
x["Num"]=x["alt"]+x["ref"]
#print(head(x))
#tab=x[,c("Sample","Allele","Count")] %>% spread(Sample,Count,fill=0) %>% gather(Sample,Count,-Allele) %>% unite(Nam,Sample,Allele,sep="_") %>% spread(Nam,Count)
#meta=data.frame(Nam=colnames(tab))
#meta["Sample"]=
#fit=glm.nb(Count~Allele+Sample+Condition*Allele,x)
##Using ado fit=tryCatch({betabin(form,~1,data=x)},error=function(cond){return(NULL)})
fit=tryCatch({vglm(form,betabinomial,data=x)},error=function(cond){return(NULL)})
#fit2=tryCatch({betabin(form2,~1,data=x)},error=function(cond){return(NULL)})
#fit=tryCatch({glm.nb(form,x)},error=function(cond){return(NULL)})
#fit=tryCatch({glm.nb(Count~Allele+Condition*Allele+Sample,x)},error=function(cond){print("Error!");return(NULL)})
#print("Yay!")
if(is.null(fit))
{
return(NULL)
}

#print("NotNull")
#if(is.null(fit2))
#{
#return(NULL)
#}
#print("NotNull2")


#if(!is.null(fit$th.warn)){fit=glmer(Count~offset(log(tot))+Allele+(1|sample),x,family=poisson)}

#anov=anova(fit,fit2,test="LRT")
#print(str(anov))

##Using AOD coef=summary(fit@Coef)
#coef=coef(fit)
coef=coef(summaryvglm(fit))
coef=data.frame(coef)
#print(head(coef))
#coef["pval"]=0
#print(rownames(coef))
coef["Test"]=rownames(coef)
#coef=coef[grep(paste("^",coef_ret,sep=""),rownames(coef)),]
#if(coef_ret=="(Intercept)")
#{
#wld=wald.test(b = coef(fit), Sigma = vcov(fit), Terms = 1,H0=c(log(2)))
#pval=wld$result$chi2[3]
#coef[4]=pval
#}
coef["SNP"]=cur_feat
coef["NumSamp"]=length(unique(x$Sample))
#print(coef)
#coef["Theta"]=theta
#coef["Count"]=sum(x[,"Count"])
#coef["Warn"]=warned
return(coef)
#return(fit)

})
print("Fitted")
toc()
#names(bySNP)=feats[1:10]
bySNP=out
names(bySNP)=nams
bySNP[sapply(bySNP, is.null)]=NULL

ret=do.call(rbind,bySNP)

ret=data.frame(ret)
#ret["SNP"]=names(bySNP)

colnames(ret)[4]="pval"
#print(head(ret))

ret=ret[order(ret$pval),]
ret["padj"]=p.adjust(ret[,"pval"])
ret["logP"]=log(2)+pnorm(abs(ret[,"z.value"]),lower.tail=FALSE, log.p=TRUE)
ret["Gene"]=as.character(lapply(ret[,"SNP"],function(x){strsplit(x,split="_")[[1]][1]}))
rownames(ret)=NULL
return(ret)

}
