#' Normalize data
#'
#' An alternative to the NormalizeData from Seurat
#'
#' @param seur The seurat object
#' @return Seurat object with normalization
#' @export
NormData<-function(seur)
{
dat=seur@assays$RNA@counts
scaleBy=seur@meta.data[,"nCount_RNA"]
if(class(dat)=="dgCMatrix")
{
dat@x<-dat@x/rep.int(scaleBy,diff(dat@p))
}
if(class(dat)=="dgTMatrix")
{
dat@x/scaleBy[dat@j+1L]
}
dat@x=log(10000*dat@x+1)
seur@assays$RNA@data=dat
print(dim(dat))
rm(dat)
return(seur)
}



#'Loads data from the pipeline as used in practice
#'
#'
#'Loads the data into a Seurat object
#' @param fils List of directories
#' @param var_gene_only If true only load variable genes
#' @param numVar Number of samples a gene must be variable in to be used
#' @param minGenes Min number of genes
#' @param regress Covariates to regress out
#' @return A Seurat object
#' @export
LoadData<-function(fils,var_gene_only=T,numVar=10,minGenes=500,regress="nFeature_RNA")
{
dirs=fils
nams=sub("^","samp",1:length(fils))
if(length(fils)==length(names(fils)))
{
nams=names(fils)
}
vargenes=c()
if(var_gene_only)
{
vargenes=sub("$","/Seurat/res.var.genes.txt",dirs)
vargenes=lapply(vargenes,function(x){scan(x,"")})
vargenes=table(do.call(c,vargenes))
vargenes=names(vargenes)[vargenes>numVar]
print(paste("Number variable genes:",length(vargenes)))
}
lst=sub("$","/STARSolo/output/resultsSolo.out/GeneFull/filtered",dirs)
print("Read in data")
dat=ReadInSTARSolo(lst,var_genes=vargenes)
print("Get cells to use and celltype labels")
cells=c()
celltype=c()
celltype_score=c()
nUMI=c()
nGene=c()
for(nam in nams)
{
print(nam)
dir=dirs[nam]
print(dir)
seur=readRDS(sub("$","/Seurat/res.seur.RDS",dir))
meta=read.table(sub("$","/Seurat/res.azimuth.cortex.txt",dir),sep="\t",header=T)
seur_nams=sub("^",paste(nam,"_",sep=""),as.character(lapply(names(seur@active.ident),function(x){s=strsplit(x,"_")[[1]];s[length(s)]})))
cells=c(cells,seur_nams)
print(dim(meta))
print(dim(seur@meta.data))
rownames(meta)=seur_nams
#rownames(meta)=sub("^",paste(nam,"_",sep=""),names(seur@active.ident))
celltype=c(celltype,meta[,"predicted.subclass"])
celltype_score=c(celltype_score,meta[,"predicted.subclass.score"])
print("UMI")
nUMI=c(nUMI,meta[,"nCount_RNA"])
nGene=c(nGene,meta[,"nFeature_RNA"])
}
names(celltype)=cells
names(celltype_score)=cells
names(nUMI)=cells
names(nGene)=cells
print(head(celltype))
print(head(celltype_score))
print(length(cells))
print(head(cells))
print(head(colnames(dat)))
cells=intersect(cells,colnames(dat))
print(length(cells))
dat=dat[,cells]
print("Create Seurat")
seur<-CreateSeuratObject(dat,"Seurat",min.features=0)
seur@meta.data["predicted.subclass"]=celltype[names(seur@active.ident)]
seur@meta.data["predicted.subclass.score"]=celltype_score[names(seur@active.ident)]
seur@meta.data["nGene"]=nGene[names(seur@active.ident)]
seur@meta.data["nUMI"]=nUMI[names(seur@active.ident)]
seur=subset(seur,nGene>minGenes)
print(length(seur@active.ident))
print(head(seur@meta.data))
print("Process RNA data")
#seur<-NormalizeData(seur,normalization.method="LogNormalize",scale.factor=10000)
seur=NormData(seur)
print(head(seur@assays$RNA@data))
#seur<-FindVariableFeatures(seur)
seur@assays$RNA@var.features=vargenes
seur<-ScaleData(seur,features=seur@assays$RNA@var.features,vars.to.regress=regress)
seur<-RunPCA(seur,npcs=60)
seur<-RunUMAP(seur,dims=1:20)
seur<-FindNeighbors(seur,dims=1:20)
seur<-FindClusters(seur)

print("Add metadata info")
seur@meta.data["SNP"]=""
seur@meta.data["Allele"]=""
seur@meta.data["CBC"]=""

#seur@meta.data["SNP"]=sub("$","/SNPLevelCounts/comb.bed",dirs)
#seur@meta.data["Allele"]=sub("$","/AlleleCounts/counts.txt",dirs)
#seur@meta.data["CBC"]=as.character(lapply(names(seur@active.ident),function(x){s=strsplit(x,"_")[[1]];s[length(s)]}))

seur@meta.data["Scrub"]=-1
#seur@meta.data["predicted.subclass.score"]=-1
#seur@meta.data["predicted.subclass"]="NA"
seur@meta.data["Sample"]=""
for(nam in nams)
{
dir=dirs[nam]
print(dir)
cbc=scan(sub("$","/STARSolo/output/resultsSolo.out/GeneFull/filtered/barcodes.tsv",dir),"")
cbc=sub("^",paste(nam,"_",sep=""),cbc)
scrub=scan(sub("$","/Scrublet/scrub.txt",dir))
names(scrub)=cbc
inter=intersect(names(seur@active.ident),cbc)
seur@meta.data[inter,"Sample"]=nam
seur@meta.data[inter,"Scrub"]=scrub[inter]
seur@meta.data[inter,"SNP"]=sub("$","/SNPLevelCounts/comb.bed",dir)
seur@meta.data[inter,"Allele"]=sub("$","/AlleleCounts/counts.txt",dir)
#meta=read.table(sub("$","/Seurat/res.azimuth.cortex.txt",dir),sep="\t",header=T)
#print(dim(meta))
#print(length(cbc))
#rownames(meta)=cbc

#seur@meta.data[inter,"predicted.subclass.score"]=meta[inter,"predicted.subclass.score"]
#seur@meta.data[inter,"predicted.subclass"]=meta[inter,"predicted.subclass"]

}
#seur@meta.data["Scrub"]=sub("$","",dirs) 

#seur@meta.data["Azimuth"]=sub("$","",dirs)

return(seur)
}
