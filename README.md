# scAlleleExpression


This R package is build to use single cell or nuclei allele specific expression (ASE) data to investigate Cell Type specific regulatory genetics (eQTLs, etc). It is build to work with the output of our ASE pipeline (https://github.com/seanken/ASE_pipeline), but ourdownstream analysis tools should work with most allele specific expression pipelines.

The package is built to:
1) Load the output of our ASE pipeline into a matrix or a Seurat object.
2) Load snps or gene level allele specific expression information in either a pseudobulk or cell specific count matrix.
3) Perform varying types of allele specific expression analysis.

Tasks 1 and 2 are pipeline specific, but taks 3 is more general.

## Install

Details to be added.

## Load pipeline output

The pipeline has many outputs. In most use cases the first step would be to load the expression data output by STARSolo. The simplest way to do this: need a named array `dirs` with a list of directories, each directory being the filtered STARSolo output from the pipeline, the names being the sample names. Then one can run 

```
dat=ReadInSTARSolo(dirs)
```

One can then process the data with Seurat or a similiar tool. To prepare to load the ASE data, one needs a data frame of meta data `meta` and a named character array, `lst` with one entry for each sample pointing to the output directory from the pipeline for that sample, named by sample (should agree with the naming in meta/in dirs). The meta data needs at least 2 columns: a column telling which sample a cell is from (for Seurat will usually be the orig.ident column) and a column telling the cells CBC. Then one can run

```
meta=LoadPipeline(meta, lst, samp_col = "orig.ident", cbc_col = "CBC")
```

This will return a version of `meta` with other columns used for loading ASE. Can also add the information into the meta data manually.

## Extracting ASE counts

To load the SNP level ASE data, use the command `GetSNPs`. So with the meta data produced above run:

`dat=GetSNPs(meta=meta,genes=genes,snp_list=snps,cond = "condition",cellType = "CellType",bulk = T)`

Here `genes` is a list of genes (can have duplicates) we want ASE for, `snps` is a list of snps (with the same length as `genes`) where each entry tells the program which allele to use as the reference vs alternate for the corresponding gene (if heterozygous at that SNP the reference allele for the gene will be the one with the reference allele for the snps, if homozygous returns nothing), a column with condition information, a column with cell type information (optional), and an argument (bulk) that tells the program to return a psudobulk if True and return cell level information if False. Note if no list of snps is passed will return information for all snps in the genes given. Note can also specify columns (which column has sample level infomration, etc) in GetSNPs for non-standard meta data set ups.

## Downstream testing

Once we have the allele specific expression information in the data frame `dat` generated above, there are many possibly downstream analysis one might want to run. The most basic is testing for ASE. There are two methods for doing this, the mixed model based one (requires bulk=F when loading `dat`) and the pseudobulk based one (requires bulk=T when loading `dat`). For the pseudobulk based one, if want to look at one cell type, can run

```
dat2=dat[dat$CellType==celltype,]
mrk=TestSNP_aod(dat2)
```

Where `mrk` are the results. Note this tests the intercept term, but can use other models as well by setting `form` (needs to be consistant with aod style formulas).

This work was funded by Aligning Science Across Parkinson's [grant # ASAP-000301] through the Michael J. Fox Foundation for Parkinson's Research (MJFF).
