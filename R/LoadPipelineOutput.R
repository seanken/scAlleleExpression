#'Make Seurat object
#'
#'Makes a Seurat object with gene level expression and allele specific expression. Meant to work with the loadIntoMatrices function below.
#'
#' @param dat_gene A matrix of gene expression, cells in the columns, genes in the rows
#' @param dat_allele A matrix of allele specific information, cells in the columns, alleles in the rows
#' @param minGenes The min number of genes below which a cell will be removed
#' @param regress The parameters to be regressed out in the ScaleData step
#' @return A Seurat object with gene level data in the RNA assay, Allele specific data in the ALLELE assay
#' @export
MakeSeurat<-function(dat_gene,dat_allele,minGenes=500,regress="nFeature_RNA")
{
print("Load into Seurat")
seur<-CreateSeuratObject(dat_gene,"Seurat",min.features=minGenes)#,normalization.method="LogNormalize",scale.factor=10000)
print("Process RNA data")
seur<-NormalizeData(seur,normalization.method="LogNormalize",scale.factor=10000)
seur<-FindVariableFeatures(seur)
seur<-ScaleData(seur,features=seur@assays$RNA@var.features,vars.to.regress=regress)
seur<-RunPCA(seur,npcs=60)
seur<-RunUMAP(seur,dims=1:20)
seur<-FindNeighbors(seur,dims=1:20)
seur<-FindClusters(seur)

print("Add Allele data")
seur[["ALLELE"]]=CreateAssayObject(dat_allele[,names(seur@active.ident)])

print("Done!")
return(seur)

}




#'Load Pipeline Into Matrix
#'
#'Loads the results of the pipeline into a list of matrices
#'
#' @param dirs A named list of directories from the AI pipeline
#' @return A list of matrices, one with Gene level Expression (named Expression) and one with allele level expression for each gene (names Allele)
#' @export
loadIntoMatrices<-function(dirs)
{
ret=list()
nams=names(dirs)
print("Load expression data")
lst=sub("$","/STARSolo/output/resultsSolo.out/GeneFull/filtered",dirs)
dat=ReadInSTARSolo(lst)

print("Load gene level allelic expression")
lst=sub("$","/AlleleCounts/counts.txt",dirs)
tab=ReadInAllele(lst)
inter=intersect(colnames(dat),colnames(tab))
dat=dat[,inter]
tab=tab[,inter]

ret[["Expression"]]=dat
ret[["Allele"]]=tab

print("Load dataframe of SNP data")
lst=sub("$","/SNPLevelCounts/comb.bed",dirs)
print(lst)
snps=ReadInSNPs(lst)
ret[["SNPs"]]=snps


return(ret)
}


#' Read in Gene level allele information
#'
#'Given a list of allele level count data from the AI pipeline reads the results into a matrix
#' @param fils A list of files with gene level allele specific expression data (in the AlleleCounts/counts.txt subdirectory produced by AI pipeline)
#' @return A sparse matrix of gene level allele specific expression data
#' @export
ReadInAllele<-function(fils)
{
nams=sub("^","samp",1:length(fils))
if(length(fils)==length(names(fils)))
{
nams=names(fils)
}
print(nams)

print("Read in each matrix")
mats=lapply(1:length(fils),function(x){
cnts=read.table(fils[x],stringsAsFactors=F)
colnames(cnts)=c("CBC","Gene","Allele","UMI")
cnts["CBC"]=sub("^",paste(nams[x],"_",sep=""),cnts[,"CBC"])
cnts=cnts[cnts$Allele!="Ambig",]
cnts<-cnts %>% unite(Feature,Gene,Allele,sep=":")
return(cnts)
})

comb=do.call(rbind,mats)

print("Change from tidy to matrix")
mat=cast_sparse(comb,"Feature","CBC","UMI")

return(mat)

}



##Read in STARSolo data, given as list of directories (dirs), giving samples names from names of dirs if there, otherwise based off order
##Can run empty dropleta utility if runFilter=T
#'Read in STARSolo data
#'
#'Reads data out of STARSolo and into a matrix
#'
#' @param dirs A (possibly named) list of directories with STARSolo output
#' @param runFilts If T performs filtering with EmptyDroplets
#' @param numCells If filtering, numCells is the expected number of cells
#' @param var_genes A list of variable genes, if non-empty will only return the data for those genes
#' @param numCellSamp If >0 subsamples to the given number of cells
#' @return Expression information in sparse matrix form
#' @export
ReadInSTARSolo<-function(dirs,runFilts=F,numCells=10000,var_genes=c(),numCellsSamp=0)
{

nams=sub("^","samp",1:length(dirs))
if(length(dirs)==length(names(dirs)))
{
nams=names(dirs)
}

print("Read in each matrix")
mats=lapply(1:length(dirs),function(x){
print("Load!")
print(x)
nam=nams[x]
print(nam)
dir=dirs[x]

mat=readMM(paste(dir,"/matrix.mtx",sep=""))


colnames(mat)=sub("^",paste(nam,"_",sep=""),scan(paste(dir,"/barcodes.tsv",sep=""),""))
tab=read.table(paste(dir,"/features.tsv",sep=""),stringsAsFactors=F,sep="\t")
tab["Name"]=tab[,2]
tab[duplicated(tab[,2]),"Name"]=tab[duplicated(tab[,2]),1]
rownames(mat)=tab[,"Name"]


if(length(var_genes)>0)
{
var_genes=intersect(var_genes,rownames(mat))
mat=mat[rownames(mat) %in% var_genes,]
}
if(numCellsSamp>0)
{
if(numCellsSamp>dim(mat)[2])
{
mat=mat[,sample(colnames(mat),numCellsSamp)]
}
}



if(runFilts)
{
print("Filter!")
e.out <- emptyDrops(mat)
keep=rownames(e.out[!is.na(e.out$FDR) & e.out$FDR<.01,])
mat=mat[,keep]
}
return(mat)
})

print("Combine Matrices")
mat=do.call(cbind,mats)
rm(mats)
print("Done!")
return(mat)

}


##Takes in a list of the SNP comb.bed files produced by the ASE pipeline and forms them into one table
##
#' Read in table of SNP info
#'
#' Takes in the comb.bed SNP file produced by the ASE Pipeline and returns as a data.frame
#'
#' @param fils A list of comb.bed files to load in
#' @return Returns a matric with 4 columns corresponding to the SNP, the gene, the Allele from the vcf (0|1 or 1|0) telling you if which value is in the ref vs alt, and the Sample name
#' @export
ReadInSNPs<-function(fils)
{
nams=sub("^","samp",1:length(fils))
if(length(fils)==length(names(fils)))
{
nams=names(fils)
}

print("Read in each matrix")
mats=lapply(1:length(fils),function(x){
snp=read.table(fils[x],stringsAsFactors=F)
colnames(snp)=c("SNP","Gene","Allele")
snp["Sample"]=nams[x]
return(snp)
})
snp_ret=do.call(rbind,mats)
return(snp_ret)
}

