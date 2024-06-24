source("../dist_group.R")
source("../MEV Diffu.R")
source("../main_fun.R") 
library(ggplot2)

#######################
load("Micheli_mouse.Rdata")
dat<-Micheli_mouse
cell.type=info=as.factor(Micheli_mouse.cellType)

# candidate.out = candidate.visual(dat,dim = 2, method=c("PCA", "MDS", "iMDS", "Sammon", "HLLE", "Isomap", 
#                                                   "kPCA", "LEIM", "UMAP", "tSNE","PHATE"),tsne.perplexity = c(10, 30))
# e<-candidate.out[[1]]
# name = candidate.out[[2]]
# save(e, file="all embeded dim 2.bin")
# print(name)



load("all embeded dim 3.bin")
n_components<-3
dataset = "Micheli"

dist.power.list = c(0.5)
conn.prop.list=c(0.012777777)
diffu.factor.list = c(1.5)
source("../scanner.R")
