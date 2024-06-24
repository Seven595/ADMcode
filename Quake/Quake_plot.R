source("../dist_group.R")
source("../MEV Diffu.R")
source("../main_fun.R") 
library(ggplot2)
library(mclust)
library(cluster)
library(aricode)
set.seed(2024) 

#######################
load("Quake_Smartseq2_Lung.Rdata")
dat<-Quake_Smartseq2_Lung$data
cell.type=info=as.factor(Quake_Smartseq2_Lung$data.cellType)

load("all embeded dim 3.bin")
n_components<-3
dataset = "Quake"

dist.power.list = c(0.5)
conn.prop.list=c(0.012777777)
diffu.factor.list = c(1.5)
source("../plot.R")


# load(paste("./Data/",dataset, "viz_out.Rdata"))
# load(paste("./Data/",dataset, "Diffu_out.Rdata"))

# rec<-new("list")
# ARI_NMI = new("list")
# # ensemble.out = ensemble.viz(e, names(e))

# ensemble.data=umap(as.dist(ensemble.out$ensemble.dist.mat), n_components = n_components)
# this.rec<-dist.grp(ensemble.out$ensemble.dist.mat, info, k=c(1,2,5,10,20))
# this.rec<-cbind(this.rec, umap5(x=as.matrix(ensemble.out$ensemble.dist.mat), info=info,  n_components = n_components, do.plot=FALSE))
# this.rec<-cbind(this.rec, dist.grp(as.matrix(dist(eigen(max(ensemble.out$ensemble.dist.mat)-ensemble.out$ensemble.dist.mat)$vec[,1:n_components])), info, k=c(1,2,5,10,20))[,2])
# colnames(this.rec)<-c("n_neighbors","viz_raw","viz_umap","viz_pca")
# rec[[1]]<-this.rec

# ## ARI/NMI
# k = length(unique(info))
# set.seed(2024) 
# umap_viz <- umap(ensemble.out$ensemble.dist.mat)
# print(str(umap_viz))
# set.seed(2024)
# cluster_viz <- kmeans(umap_viz, centers = k)

# ARI_viz <- adjustedRandIndex(info, cluster_viz$cluster)
# NMI_viz <- NMI(info, cluster_viz$cluster)
# rec_ari_nmi = cbind(ARI_viz,NMI_viz)
# colnames(rec_ari_nmi)<-c("ARI_viz","NMI_viz")
# ARI_NMI[[1]] = rec_ari_nmi
# print(ARI_NMI[[1]])

# set.seed(2024) 
# umap_viz <- umap(ensemble.out$ensemble.dist.mat)
# umap_viz <- as.data.frame(umap_viz)
# #使用ggplot2绘制UMAP结果
# p1  = ggplot(umap_viz, aes(x = V1, y = V2, color = info)) +
#   geom_point(alpha = 0.5) +
#   labs(title ="meta-spec",
#        x = "",
#        y = " ") +
#     theme( panel.grid.major = element_blank(),  #移除栅格线
#     panel.grid.minor = element_blank(),
#     panel.background = element_rect(fill = "white"),  # 设置背景为白色
#     panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),  # 设置绘图区域边框
#     legend.position="none",
#     axis.text = element_blank(),  # 隐藏刻度值
#     axis.ticks = element_blank(),
#     axis.title = element_text(face = "bold"),
#     plot.title = element_text(hjust = 0.5,size = 30))
# ggsave(paste("./visulization/",dataset,"viz1.png"), plot = p1, width = 6, height = 4, dpi = 420)
# # save(ensemble.out, file=paste("./Data/",dataset, "viz_out.Rdata"))



#  ### neighbours
# this.rec<-dist.grp(mev.out$diffu.dist, info, k=c(1,2,5,10,20))
# this.rec<-cbind(this.rec, umap5(x=mev.out$diffu.dist, info=info,  n_components = n_components, do.plot=FALSE))
# this.rec<-cbind(this.rec, dist.grp(as.matrix(dist(eigen(max(mev.out$diffu.dist)-mev.out$diffu.dist)$vec[,1:n_components])), info, k=c(1,2,5,10,20))[,2])

# colnames(this.rec)<-c("n_neighbors","diffu_raw","diffu_umap","diffu_pca")
# # this.name<-paste("dist.power = ",dist.power, "conn.prop = ", conn.prop,"diffu.factor= ",diffu.factor)
# rec[[2]]<-this.rec
# #    names(rec)[i]<-this.name

#     ##  ARI/NMI
# set.seed(2024) 
# umap_diffu <- umap(mev.out$diffu.dist)
    
# set.seed(2024) 
# cluster_diffu <- kmeans(umap_diffu, centers = k)
# ARI_mev <- adjustedRandIndex(info, cluster_diffu$cluster)
# NMI_mev <- NMI(info, cluster_diffu$cluster)
# rec_ari_nmi = cbind(ARI_mev,NMI_mev)
# colnames(rec_ari_nmi)<-c("ARI_diffu","NMI_diffu")
# ARI_NMI[[2]] = rec_ari_nmi
# print(rec[[1]])
# print(ARI_NMI[[1]])
# print(rec[[2]])
# print((ARI_NMI[[2]]))
    
# umap_df <- as.data.frame(umap_diffu)   
# p  = ggplot(umap_df, aes(x = V1, y = V2, color = info)) +
#   geom_point(alpha = 0.5) +
#   labs(title = "meta-diffu",
#        x = "",
#        y = " ") +
#     theme( panel.grid.major = element_blank(),  #移除栅格线
#     panel.grid.minor = element_blank(),
#     panel.background = element_rect(fill = "white"),  # 设置背景为白色
#     panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),  # 设置绘图区域边框
#     legend.position="none",
#     axis.text = element_blank(),  # 隐藏刻度值
#     axis.ticks = element_blank(),
#     axis.title = element_text(face = "bold"),
#     plot.title = element_text(hjust = 0.5,size = 30,family = "Arial"))
# ggsave(paste("./visulization/",dataset,"Diffu_good.png"), plot = p, width = 6, height = 4, dpi = 420)