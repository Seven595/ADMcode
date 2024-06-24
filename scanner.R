library(mclust)
library(cluster)
library(aricode)
set.seed(2024) 

rec<-new("list")
ARI_NMI = new("list")

ensemble.out = ensemble.viz(e, names(e))
ensemble.data=umap(as.dist(ensemble.out$ensemble.dist.mat), n_components = n_components)
this.rec<-dist.grp(ensemble.out$ensemble.dist.mat, info, k=c(1,2,5,10,20))
this.rec<-cbind(this.rec, umap5(x=as.matrix(ensemble.out$ensemble.dist.mat), info=info,  n_components = n_components, do.plot=FALSE))
this.rec<-cbind(this.rec, dist.grp(as.matrix(dist(eigen(max(ensemble.out$ensemble.dist.mat)-ensemble.out$ensemble.dist.mat)$vec[,1:n_components])), info, k=c(1,2,5,10,20))[,2])
colnames(this.rec)<-c("n_neighbors","viz_raw","viz_umap","viz_pca")
rec[[1]]<-this.rec


## ARI/NMI
k = length(unique(info))
set.seed(2024) 
umap_viz <- umap(ensemble.out$ensemble.dist.mat)
print(str(umap_viz))
set.seed(2024)
cluster_viz <- kmeans(umap_viz, centers = k)

ARI_viz <- adjustedRandIndex(info, cluster_viz$cluster)
NMI_viz <- NMI(info, cluster_viz$cluster)
rec_ari_nmi = cbind(ARI_viz,NMI_viz)
colnames(rec_ari_nmi)<-c("ARI_viz","NMI_viz")
ARI_NMI[[1]] = rec_ari_nmi
print(ARI_NMI[[1]])

set.seed(2024) 
umap_viz <- umap(ensemble.out$ensemble.dist.mat)
umap_viz <- as.data.frame(umap_viz)
# 使用ggplot2绘制UMAP结果
# p1 = ggplot(umap_viz, aes(x = V1, y = V2, color = info)) +
# geom_point(alpha = 0.5) +
# labs(title = paste(dataset,"viz","ARI = ",ARI_viz,"NMI=", NMI_viz),
#        x = "Dimension 1",
#        y = "Dimension 2") +
#   theme_minimal()
# ggsave(paste("./visulization/",dataset,"viz.png"), plot = p1, width = 6, height = 4, dpi = 300)
# save(ensemble.out, file=paste("./Data/",dataset, "viz_out.Rdata"))



p1  = ggplot(umap_viz, aes(x = V1, y = V2, color = info)) +
  geom_point(alpha = 0.5) +
  labs(title = paste(dataset,"meta-viz"),
       x = "",
       y = " ") +
    theme( panel.grid.major = element_blank(),  #移除栅格线
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white"),  # 设置背景为白色
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),  # 设置绘图区域边框
    legend.position="none",
    axis.text = element_blank(),  # 隐藏刻度值
    axis.ticks = element_blank(),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5))
ggsave(paste("./visulization/",dataset,"viz1.png"), plot = p1, width = 6, height = 4, dpi = 300)


# dist.power.list = c(0.5,1)
# conn.prop.list=seq(0.003,0.025,length.out = 10)
# diffu.factor.list = seq(1.5,3.5,length.out = 4)


# dist.power.list = c(1)
# conn.prop.list=c(0.0225555555555556)
# diffu.factor.list = c(1.5)

combos<-expand.grid(dist.power.list,conn.prop.list,diffu.factor.list)

for(i in 1:nrow(combos))
{
   dist.power = combos[i,1]
	conn.prop=combos[i,2]
   diffu.factor= combos[i,3]
	mev.out = mev(e, dist.power=dist.power, conn.prop=conn.prop,diffu.factor = diffu.factor)
   save(mev.out, file=paste("./Data/",dataset, "Diffu_out.Rdata"))
    ### neighbours
	this.rec<-dist.grp(mev.out$diffu.dist, info, k=c(1,2,5,10,20))
	this.rec<-cbind(this.rec, umap5(x=mev.out$diffu.dist, info=info,  n_components = n_components, do.plot=FALSE))
	this.rec<-cbind(this.rec, dist.grp(as.matrix(dist(eigen(max(mev.out$diffu.dist)-mev.out$diffu.dist)$vec[,1:n_components])), info, k=c(1,2,5,10,20))[,2])
	
	colnames(this.rec)<-c("n_neighbors","diffu_raw","diffu_umap","diffu_pca")
	this.name<-paste("dist.power = ",dist.power, "conn.prop = ", conn.prop,"diffu.factor= ",diffu.factor)
	rec[[i+1]]<-this.rec
   names(rec)[i]<-this.name

    ##  ARI/NMI
   set.seed(2024) 
   umap_diffu <- umap(mev.out$diffu.dist)
    
   set.seed(2024) 
   cluster_diffu <- kmeans(umap_diffu, centers = k)
   ARI_mev <- adjustedRandIndex(info, cluster_diffu$cluster)
   NMI_mev <- NMI(info, cluster_diffu$cluster)
   rec_ari_nmi = cbind(ARI_mev,NMI_mev)
   colnames(rec_ari_nmi)<-c("ARI_diffu","NMI_diffu")
   ARI_NMI[[i+1]] = rec_ari_nmi
   names(ARI_NMI)[i+1]<-this.name
   print(rec[[1]] )
   print(ARI_NMI[[1]])
   print(i)
	print(this.name)
	print((rec[[i+1]]))
   print((ARI_NMI[[i+1]]))
    
   umap_df <- as.data.frame(umap_diffu)
#    p = ggplot(umap_df, aes(x = V1, y = V2, color = info)) +
#       geom_point(alpha = 0.5) +
#       labs(title = paste(dataset,"diffu","ARI = ",ARI_mev,"NMI=", NMI_mev),
#        x = "Dimension 1",
#        y = "Dimension 2") +
#       theme_minimal()
    
   p  = ggplot(umap_df, aes(x = V1, y = V2, color = info)) +
  geom_point(alpha = 0.5) +
  labs(title = paste(dataset,"meta-diffu"),
       x = "",
       y = " ") +
    theme( panel.grid.major = element_blank(),  #移除栅格线
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white"),  # 设置背景为白色
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),  # 设置绘图区域边框
    legend.position="none",
    axis.text = element_blank(),  # 隐藏刻度值
    axis.ticks = element_blank(),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5))
#    ggsave(paste("./visulization/",dataset,"Diffu", this.name,".png"), plot = p, width = 6, height = 4, dpi = 300)
    ggsave(paste("./visulization/",dataset,"Diffu_good.png"), plot = p, width = 6, height = 4, dpi = 300)
}
# save(rec, file=paste("search_CSRA_gamma.Rdata"))
# save(rec, file=paste("search_ari_nmi_gamma.Rdata"))