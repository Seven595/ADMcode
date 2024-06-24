
source("../dist_group.R")
source("../MEV Diffu.R")
source("../main_fun.R") #https://github.com/rongstat/meta-visualization
library(ggplot2)

#######################
load("Oihane.Rdata")
dat<-Oihane
cell.type=info=as.factor(Oihane.info.cellType)

# candidate.out = candidate.visual(dat,dim = 2, method=c("PCA", "MDS", "iMDS", "Sammon", "HLLE", "Isomap", 
#                                                  "kPCA", "LEIM", "UMAP", "tSNE","PHATE"),tsne.perplexity = c(10, 30))
# e<-candidate.out[[1]]
# name = candidate.out[[2]]
# print(name)
# save(e, file="all embeded dim 2.bin")


# dist.power.list = c(1)
# conn.prop.list=c(0.0225555555555556 )
# diffu.factor.list = c(1.5)

load("all embeded dim 3.bin")
n_components<-3
dataset= "Oihane"
source("../scanner.R")

