
source("../dist_group.R")
source("../MEV Diffu.R")
source("../main_fun.R") 
library(ggplot2)

#######################
load("Quake_Smartseq2_Lung.Rdata")
dat<-Quake_Smartseq2_Lung$data
cell.type=info=as.factor(Quake_Smartseq2_Lung$data.cellType)


# candidate.out = candidate.visual(dat,dim = 2, method=c("PCA", "MDS", "iMDS", "Sammon", "HLLE", "Isomap", 
#                                                  "kPCA", "LEIM", "UMAP", "tSNE","PHATE"),tsne.perplexity = c(10, 30))
# e<-candidate.out[[1]]
# name = candidate.out[[2]]
# save(e, file="all embeded dim 2.bin")
# print(name)

dist.power.list = c(0.5)
conn.prop.list=c(0.00788888888888889)
diffu.factor.list = c(3.5)


dataset = "Quake"
load("all embeded dim 3.bin")
n_components<-3
source("../scanner.R")
