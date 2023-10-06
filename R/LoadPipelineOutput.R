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


