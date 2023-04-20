#' Draw boxplot
#' 
#' @param out ASE expression info
#' @param nam Name of snp/gene pair you want to plot
#' @param minExp Min number of phased UMIs to plot
#' @return A plot!
#' @export 
drawBoxPlot=function(out,nam,minExp=5)
{
    tab<-out %>% group_by(SNP,Gene,CellType,Sample,Allele) %>% summarise(Count=sum(Count)) %>% spread(Allele,Count,fill=0) %>% as.data.frame()
    tab=tab %>% unite(Name,Gene,SNP,sep="_")
    tab["Ratio"]=tab["alt"]/(tab["alt"]+tab["ref"])
    tab["Tot"]=tab["alt"]+tab["ref"]

    tab=tab[tab$Name==nam,]
    alt=strsplit(nam,":")[[1]][4]

    p= ggplot(tab[tab$Tot>minExp,],aes(x=CellType,y=Ratio))+geom_boxplot()+coord_flip()+ylab(paste("Proportion phase UMI mapping to alt allele (",alt,")",sep=""))+geom_hline(yintercept=.5,linetype="dotted")+ggtitle(snp)+ylim(c(0,1))+geom_point(position = position_jitter(w = 0.1, h = 0))
    return(p)
}