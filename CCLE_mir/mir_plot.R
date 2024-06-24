source("../dist_group.R")
source("../MEV Diffu.R")
source("../main_fun.R") 
library(ggplot2)
library(mclust)
library(cluster)
library(aricode)
set.seed(2024) 

#######################
load("CCLE_miRNA_20181103.bin")
load("CCLE_RNAseq_genes_rpkm_20180929.bin")

mir<-log10(1+mir)
gene<-log10(1+gene)
gene<-as.matrix(gene)
mir<-as.matrix(mir)

mir<-mir[,which(colnames(mir) %in% colnames(gene))]
gene<-gene[,which(colnames(gene) %in% colnames(mir))]

mir<-mir[,order(colnames(mir))]
gene<-gene[,order(colnames(gene))]

sum(colnames(mir)==colnames(gene))
dim(mir)

mir<-t(mir)
gene<-t(gene)

cv.mir<-apply(mir,2,sd)/apply(mir,2,mean)
mir<-mir[,cv.mir>=0.1]

ze<-apply(gene==0,2,sum)
gene<-gene[,which(ze<=0.25*nrow(gene))]
cv.gene<-apply(gene,2,sd)/apply(gene,2,mean)
gene<-gene[,cv.gene>=0.5]

############################

cells<-rownames(gene)
cell.type<-cells
for(i in 1:length(cells))
{
	this<-strsplit(cells[i], "_")[[1]][-1]
	this<-paste(this, collapse="_")
	cell.type[i]<-this
}

ttt<-table(cell.type)
sel<-which(cell.type %in% names(ttt)[ttt<10])
cell.type[sel]<-NA

anno<-read.table("Cell_lines_annotations_20181226.txt",header=T,sep="\t")
anno<-anno[which(anno[,1] %in% cells),]
anno<-anno[order(anno[,1]),]

##### limit to cells with annotations
s.anno<-which(cells %in% anno[,1] & !is.na(cell.type))
cells<-cells[s.anno]
cell.type<-cell.type[s.anno]
mir<-mir[s.anno,]
gene<-gene[s.anno,]

mir2<-mir
gene2<-gene

for(i in 1:ncol(mir2)) mir2[,i]<-(mir2[,i]-mean(mir2[,i]))/sd(mir2[,i])
for(i in 1:ncol(gene2)) gene2[,i]<-(gene2[,i]-mean(gene2[,i]))/sd(gene2[,i])

all.colors<-c("grey50","green","blue","cyan", "yellow","orange","red","black", "wheat","purple","darkblue","dodgerblue4","darkred","darkorange","darkcyan","magenta","firebrick","khaki4","cornsilk3","darkgoldenrod4")

################### 
dat<-mir
info<-as.factor(cell.type)
load("all embeded dim 3.bin")
n_components<-3
dataset = "mir"


source("../plot.R")