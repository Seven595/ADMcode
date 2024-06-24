#library(dynamicTreeCut)
library(DDoutlier)
library(diffudist)
library(fitdistrplus)
library(igraph)
library(MixGHD)
library(BiocParallel)

zp.quantile <- function(x,y)
		  # x and y are two vectors of equal length
		  # x serves as the quantile normalization template. Only updated y is returned
		{
		  o.x<-order(x)
		  #r.x<-rank(x, ties.method = "random")

		  o.y<-order(y)
		  r.y<-rank(y, ties.method = "average")

		  x<-x[o.x]
		  y<-y[o.y]

		  #x2<-x[x>0]
		  #y2<-y[y>0]

		  z.x<-seq(0,1,length.out=length(x))
		  z.y<-seq(0,1,length.out=length(y))

		  new.y<-stats::approx(x=z.x, y=x, xout=z.y)$y
		  #y[y>0]<-new.y2
		  y<-new.y[r.y]
		}



mev_info<-function(e, k.dim=NULL, dist.power=0.5, conn.prop=0.02, raw.d.pwr=.5, diffu.steps=NA, diffu.factor=2, distr.template="gamma", gamma.shape=3, gamma.rate=3, scale.dist=TRUE, symmetrize="mean", dist.quantile=0.25)  #扩散距离转换为信息距离
{
	########################
	if(is.null(k.dim)) k.dim=ncol(e[[1]])

	#this.e is one dimension reduction result
	
	fake.fun<-function(this.e, raw.d.pwr, diffu.steps, conn.prop, scale.dist, symmetrize, diffu.factor, dist.quantile)
	{
	
		library(DDoutlier)
		library(diffudist)
		library(MixGHD)
		library(BiocParallel)
		library(fitdistrplus)
		
		zp.quantile <- function(x,y)
		  # x and y are two vectors of equal length
		  # x serves as the quantile normalization template. Only updated y is returned
		{
		  o.x<-order(x)
		  #r.x<-rank(x, ties.method = "random")

		  o.y<-order(y)
		  r.y<-rank(y, ties.method = "average")

		  x<-x[o.x]
		  y<-y[o.y]

		  #x2<-x[x>0]
		  #y2<-y[y>0]

		  z.x<-seq(0,1,length.out=length(x))
		  z.y<-seq(0,1,length.out=length(y))

		  new.y<-stats::approx(x=z.x, y=x, xout=z.y)$y
		  #y[y>0]<-new.y2
		  y<-new.y[r.y]
		}

		estimate_mode <- function(x) {
		  d <- density(x)
		  d$x[which.max(d$y)]
		}
	
	
		move.outlier<-function(x, d=NULL, fraction=0.01)
		{
			x.is.vector<-F
			if(is.null(nrow(x))) 
			{
				x<-matrix(x, ncol=1)
				x.is.vector<-TRUE
			}
			
			for(n in 1:ncol(x))
			{
				this.x<-x[,n]
				if(is.null(d))
				{
					this.d<-diff(quantile(this.x,c(0.25,0.75)))/3
				}else{
					this.d<-d
				}
				
				this.out<-DB(matrix(this.x,ncol=1), this.d, fraction)
				this.sel<-which(this.out$classification == "Outlier")
				
				if(length(this.sel)>0)
				{
					new.x<-zp.quantile(this.x[-this.sel], this.x)
					x[,n]<-new.x
				}
			}
			x
		}
		
		this.e<-move.outlier(this.e)
		this.d<-as.matrix(dist(this.e))
     write.csv(this.d,"./Data/distance.csv")
		n<-nrow(this.d)
		
		diag(this.d)<-Inf
		
		######## global
		
		conn.cutoff<-quantile(as.dist(this.d), conn.prop)
		this.conn<-1*(this.d <= conn.cutoff)
		write.csv(this.conn,"./Data/global_connec.csv")
		########## local
		
		for(i in 1:nrow(this.d))
		{
			sel<-which(this.d[i,]<=quantile(this.d[i,], conn.prop))
			this.conn[i,sel]<-1
		}
		write.csv(this.conn,"./Data/global_local_connec.csv")
		this.graph<-graph_from_adjacency_matrix(this.conn, mode="undirected")
		#this.d[which(this.d==0)] <- min(this.d[which(this.d>0)])
		
		pmat<-1/this.d^raw.d.pwr
		diag(pmat)=0
		pmat[pmat == Inf]<-max(pmat[pmat != Inf], na.rm=TRUE)
		pmat[is.na(pmat)]<-max(pmat, na.rm=TRUE)
		
		pmat<-pmat*this.conn
     write.csv(pmat,"./Data/pmat.csv")
		for(i in 1:nrow(pmat)) pmat[i,]<-pmat[i,]/sum(pmat[i,])
     write.csv(pmat,"./Data/pmat.csv")
		if(is.na(diffu.steps[1])) 
		{
			sp.mat<-shortest.paths(this.graph)
			sp.mat[sp.mat==Inf]<-NA
			mean.steps<-apply(sp.mat,1,quantile, na.rm=TRUE, probs=dist.quantile)
			diffu.steps<-quantile(mean.steps, c(0.1, 0.3, 0.5, 0.7, 0.9), na.rm=TRUE)*diffu.factor
        write.csv(sp.mat,"./Data/shortest_path.csv")
		}
      
		
		# sp.mat is shortest distance
		# diffu.steps is all the steps we try
		
		all.d<-new("list")
		for(i in 1:length(diffu.steps))
		{
            step = diffu.steps[i]
            
            pmat_pwr = pmat^step+ 1e-7
            
            pmat_log = -1 * log(pmat_pwr)
            for(i in 1:nrow(pmat_log)) pmat_log[i,]<-pmat_log[i,]/sum(pmat_log[i,])
			this.d<-as.matrix(suppressMessages(get_distance_matrix_from_T(pmat_log, step)))
            
         write.csv(this.d,"./Data/diffusion_distance"+str(i)+".csv")
			if(scale.dist) 
			{
				this.d<- 1-cor(this.d, method="spearman")
				this.d[this.d>1]<-1
                
			}
			all.d[[i]]<-this.d
		}
		
		# merge.d is the matrix in which each point takes different diffusion time
		# the steps for each point is based on its mean distance to other points
		
		merge.d<-all.d[[1]]
		for(i in 1:nrow(merge.d))
		{
			this.dist<-abs(mean.steps[i]*diffu.factor-diffu.steps)
			this.closest<-which(this.dist==min(this.dist))[1]
			merge.d[i,]<-all.d[[this.closest]][i,]
		}
		this.d<-(merge.d+t(merge.d))/2
		
		if(symmetrize=="mean") {
			this.d<-this.d+t(this.d)
		}else if(symmetrize=="min") {
			this.d.2<-t(this.d)
			this.d[this.d>this.d.2]<-this.d.2[this.d>this.d.2]
		}else if(symmetrize=="max") {
			this.d.2<-t(this.d)
			this.d[this.d<this.d.2]<-this.d.2[this.d<this.d.2]			
		}

		to.return<-list(this.d=this.d)
		return(to.return)
	}

	# e is a list object. Each item is a dimension reduction result, an nXk matrix. K is the dimension.
	
	d<-bplapply(e, fake.fun, raw.d.pwr=raw.d.pwr, diffu.steps=diffu.steps, conn.prop=conn.prop, scale.dist=scale.dist, symmetrize=symmetrize, diffu.factor=diffu.factor, dist.quantile=dist.quantile)
	

	
	# the code below combines the diffusion distance matrices. A quantile normalization is taken
	# combine: merge all into a long vector as template
	# gamma: find gamma parameter from each, and take average of the parameters
	# parametric: specify gamma parameters by user
	
	if(distr.template != "none")
	{
		
		if(distr.template == "combine")
		{
			d.template<-NULL
			for(i in 1:length(d)) d.template<-c(d.template, as.dist(d[[i]][[1]]))
		}else if(distr.template == "gamma"){
			for(m in 1:length(d))
			{
				this.fit<-fitdist(as.vector(as.dist(d[[m]][[1]])),"gamma",method="mme")$estimate
				if(m == 1) all.fit=this.fit
				else all.fit=rbind(all.fit, this.fit)
			}
			ave.fit= apply(all.fit,2,mean)
			
			d.template<-rgamma(nrow(d[[1]][[1]])*ncol(d[[1]][[1]]), shape=ave.fit[1], rate=ave.fit[2])
		}else if(distr.template=="parametric"){
			d.template<-rgamma(nrow(d[[1]][[1]])*ncol(d[[1]][[1]]), shape=gamma.shape, rate=gamma.rate)
		}

		for(m in 1:length(d))
		{
			this.d<-as.dist(d[[m]][[1]])
			this.d2<-zp.quantile(d.template, this.d)
			attributes(this.d2)<-attributes(this.d)			
			this.d2<-as.matrix(this.d2)
			d[[m]][[1]]<-this.d2
		}
	}		
	
	dd<-d[[1]][[1]]*0
	for(m in 1:length(d)) dd<-dd+d[[m]][[1]]^dist.power

	to.return=list(diffu.dist=dd)
	return(to.return)
}





    
mev<-function(e, k.dim=NULL, dist.power=0.5, conn.prop=0.02, raw.d.pwr=.5, diffu.steps=NA, diffu.factor=2, distr.template="combine", gamma.shape=3, gamma.rate=3, scale.dist=TRUE, symmetrize="mean", dist.quantile=0.25)
{
	########################
	if(is.null(k.dim)) k.dim=ncol(e[[1]])

	#this.e is one dimension reduction result
	
	fake.fun<-function(this.e, raw.d.pwr, diffu.steps, conn.prop, scale.dist, symmetrize, diffu.factor, dist.quantile)
	{
	
		library(DDoutlier)
		library(diffudist)
		library(MixGHD)
		library(BiocParallel)
		library(fitdistrplus)
		
		zp.quantile <- function(x,y)
		  # x and y are two vectors of equal length
		  # x serves as the quantile normalization template. Only updated y is returned
		{
		  o.x<-order(x)
		  #r.x<-rank(x, ties.method = "random")

		  o.y<-order(y)
		  r.y<-rank(y, ties.method = "average")

		  x<-x[o.x]
		  y<-y[o.y]

		  #x2<-x[x>0]
		  #y2<-y[y>0]

		  z.x<-seq(0,1,length.out=length(x))
		  z.y<-seq(0,1,length.out=length(y))

		  new.y<-stats::approx(x=z.x, y=x, xout=z.y)$y
		  #y[y>0]<-new.y2
		  y<-new.y[r.y]
		}

		estimate_mode <- function(x) {
		  d <- density(x)
		  d$x[which.max(d$y)]
		}
	
	
		move.outlier<-function(x, d=NULL, fraction=0.01)
		{
			x.is.vector<-F
			if(is.null(nrow(x))) 
			{
				x<-matrix(x, ncol=1)
				x.is.vector<-TRUE
			}
			
			for(n in 1:ncol(x))
			{
				this.x<-x[,n]
				if(is.null(d))
				{
					this.d<-diff(quantile(this.x,c(0.25,0.75)))/3
				}else{
					this.d<-d
				}
				
				this.out<-DB(matrix(this.x,ncol=1), this.d, fraction)
				this.sel<-which(this.out$classification == "Outlier")
				
				if(length(this.sel)>0)
				{
					new.x<-zp.quantile(this.x[-this.sel], this.x)
					x[,n]<-new.x
				}
			}
			x
		}
		
		this.e<-move.outlier(this.e)
		this.d<-as.matrix(dist(this.e))
     write.csv(this.d,"./Data/distance.csv")
		n<-nrow(this.d)
		
		diag(this.d)<-Inf
		
		######## global
		
		conn.cutoff<-quantile(as.dist(this.d), conn.prop)
		this.conn<-1*(this.d <= conn.cutoff)
		write.csv(this.conn,"./Data/conoec_global.csv")
		########## local
		
		for(i in 1:nrow(this.d))
		{
			sel<-which(this.d[i,]<=quantile(this.d[i,], conn.prop))
			this.conn[i,sel]<-1
		}
		write.csv(this.conn,"./Data/conoec_adapt.csv")
		this.graph<-graph_from_adjacency_matrix(this.conn, mode="undirected")
		#this.d[which(this.d==0)] <- min(this.d[which(this.d>0)])
		
		pmat<-1/this.d^raw.d.pwr
		diag(pmat)=0
		pmat[pmat == Inf]<-max(pmat[pmat != Inf], na.rm=TRUE)
		pmat[is.na(pmat)]<-max(pmat, na.rm=TRUE)
		
		pmat<-pmat*this.conn
		for(i in 1:nrow(pmat)) pmat[i,]<-pmat[i,]/sum(pmat[i,])
     write.csv(pmat,"./Data/pmat.csv")
		if(is.na(diffu.steps[1])) 
		{
			sp.mat<-shortest.paths(this.graph)
			sp.mat[sp.mat==Inf]<-NA
			mean.steps<-apply(sp.mat,1,quantile, na.rm=TRUE, probs=dist.quantile)
			diffu.steps<-quantile(mean.steps, c(0.1, 0.3, 0.5, 0.7, 0.9), na.rm=TRUE)*diffu.factor
		}
		
		# sp.mat is shortest distance
		# diffu.steps is all the steps we try
		
		all.d<-new("list")
		for(i in 1:length(diffu.steps))
		{
			this.d<-as.matrix(suppressMessages(get_distance_matrix_from_T(pmat, diffu.steps[i])))
			if(scale.dist) 
			{
				this.d<- 1-cor(this.d, method="spearman")
				this.d[this.d>1]<-1
			}
       
			all.d[[i]]<-this.d
		}
		
		# merge.d is the matrix in which each point takes different diffusion time
		# the steps for each point is based on its mean distance to other points
		
		merge.d<-all.d[[1]]
     write.csv(pmat,"./Data/pmat.csv")
		for(i in 1:nrow(merge.d))
		{
			this.dist<-abs(mean.steps[i]*diffu.factor-diffu.steps)
			this.closest<-which(this.dist==min(this.dist))[1]
			merge.d[i,]<-all.d[[this.closest]][i,]
		}
		this.d<-(merge.d+t(merge.d))/2
		
		if(symmetrize=="mean") {
			this.d<-this.d+t(this.d)
		}else if(symmetrize=="min") {
			this.d.2<-t(this.d)
			this.d[this.d>this.d.2]<-this.d.2[this.d>this.d.2]
		}else if(symmetrize=="max") {
			this.d.2<-t(this.d)
			this.d[this.d<this.d.2]<-this.d.2[this.d<this.d.2]			
		}

		to.return<-list(this.d=this.d)
		return(to.return)
	}

	# e is a list object. Each item is a dimension reduction result, an nXk matrix. K is the dimension.
	
	d<-bplapply(e, fake.fun, raw.d.pwr=raw.d.pwr, diffu.steps=diffu.steps, conn.prop=conn.prop, scale.dist=scale.dist, symmetrize=symmetrize, diffu.factor=diffu.factor, dist.quantile=dist.quantile)
	
	################### for Junning to change #######################
	
	# the code below combines the diffusion distance matrices. A quantile normalization is taken
	# combine: merge all into a long vector as template
	# gamma: find gamma parameter from each, and take average of the parameters
	# parametric: specify gamma parameters by user
	
	if(distr.template != "none")
	{
		
		if(distr.template == "combine")
		{
			d.template<-NULL
			for(i in 1:length(d)) d.template<-c(d.template, as.dist(d[[i]][[1]]))
		}else if(distr.template == "gamma"){
			for(m in 1:length(d))
			{
				this.fit<-fitdist(as.vector(as.dist(d[[m]][[1]])),"gamma",method="mme")$estimate
				if(m == 1) all.fit=this.fit
				else all.fit=rbind(all.fit, this.fit)
			}
			ave.fit= apply(all.fit,2,mean)
			
			d.template<-rgamma(nrow(d[[1]][[1]])*ncol(d[[1]][[1]]), shape=ave.fit[1], rate=ave.fit[2])
		}else if(distr.template=="parametric"){
			d.template<-rgamma(nrow(d[[1]][[1]])*ncol(d[[1]][[1]]), shape=gamma.shape, rate=gamma.rate)
		}

		for(m in 1:length(d))
		{
			this.d<-as.dist(d[[m]][[1]])
			this.d2<-zp.quantile(d.template, this.d)
			attributes(this.d2)<-attributes(this.d)			
			this.d2<-as.matrix(this.d2)
			d[[m]][[1]]<-this.d2
		}
	}		
	for(m in 1:length(d)){
    write.csv(d[[m]],paste("./Data/diffu_dist",m, ".csv"))
    }
	dd<-d[[1]][[1]]*0
	for(m in 1:length(d)) dd<-dd+d[[m]][[1]]^dist.power

	to.return=list(diffu.dist=dd)
   write.csv(dd,"./Data/diffu_dist_merged.csv")
	return(to.return)
}


########### other functions


dist.grp.majority<-function(distmat, grp, k=seq(1,21,by=2))
{
	
	grpmat<-matrix(0,nrow=length(grp),ncol=length(grp))
	for(i in 1:length(grp)) for(j in 1:length(grp)) if(grp[i]==grp[j]) grpmat[i,j]<-1

	rec<-cbind(k,k)
	
	for(j in 1:length(k))
	{
		r<-NULL
		diag(distmat)<-Inf
		for(i in 1:nrow(distmat))
		{
			sel<-which(distmat[i,]<=quantile(distmat[i,], k[j]/ncol(distmat)))
			this<-grpmat[i,sel]
			this<-1*(sum(this) >= length(this)*0.5)
			r<-c(r,this)
		}
		rec[j,2]<-sum(r==1)/length(r)
	}
	rec
}
library(BiocParallel)

dist.grp<-function(distmat, grp, k=1:20)
{
	grp<-as.numeric(grp)
	grpmat<-as.matrix(dist(grp))
	grpmat<-1*(grpmat==0)
	for(i in 1:nrow(distmat)) distmat[i,]<-rank(distmat[i,])-1
	
	rec<-cbind(k,k)
	diag(distmat)<-Inf
	

	
	rec<-cbind(k,k)
	for(i in 1:length(k))
	{
		rankmat<-1*(distmat <= k[i])
		rec[i,2]<-sum(rankmat*grpmat)/sum(rankmat)
	}
	
	#rec[,2]<-unlist(bplapply(k, dist.fun, distmat=distmat, grpmat=grpmat))
	return(rec)
}

dist.grp.all.pairs<-function(distmat, grp)
### notice: this is for distance specifically. No minus 0.5. 
### samples from the same group is supposed to have smaller distances
{
	library(MLmetrics)
	library(pROC)
	
	grpmat<-matrix(0,nrow=length(grp),ncol=length(grp))
	for(i in 1:length(grp)) for(j in 1:length(grp)) if(grp[i]==grp[j]) grpmat[i,j]<-1
	
	grpmat<-as.dist(grpmat)
	distmat<-as.dist(distmat)
	
	r<-c(0,0)
	r[1]<-auc(roc(as.vector(grpmat)~as.vector(distmat), direction=">"))
	r[2]<-PRAUC(-as.vector(distmat), as.vector(grpmat))
	
	names(r)<-c("ROCAUC","PRAUC")
	r
}


umap5<-function(x, info, do.plot=TRUE, n_components, k=c(1,2,5,10,20))
{
	library(uwot)

	for(n in 1:5)
	{
		ensemble.data<-umap(as.dist(x), n_components = n_components)
		ds<-dist.grp(as.matrix(dist(ensemble.data)), info, k=k)[,2]
		if(n==1) rec<-ds
		else rec<-cbind(rec,ds)
	}
	if(do.plot) pairs(ensemble.data, col=rainbow(length(unique(info)))[as.numeric(info)], pch=(1:length(unique(info)))[as.numeric(info)])
	to.return<-apply(rec,1,mean)
	return(to.return)
}


umap5.parallel<-function(x, info, do.plot=TRUE, n_components, k=c(1,2,5,10,20))
{
	library(uwot)
	d<-matrix(0, ncol=5, nrow=5)
	fun<-function(x, info, n_components, k)
	{
		library(uwot)
		ensemble.data<-umap(as.dist(x), n_components = n_components)
		ds<-dist.grp(as.matrix(dist(ensemble.data)), info, k=k)[,2]
		en<-ensemble.data
		to.return<-list(ds=ds, en=en)
		return(to.return)
	}
	info.repeat<-list(i1=info,i2=info,i3=info,i4=info,i5=info)
	l<-bplapply(info.repeat, fun, x=x, n_components=n_components, k=k)
	to.return<-apply(cbind(l[[1]][[1]],l[[2]][[1]],l[[3]][[1]],l[[4]][[1]],l[[5]][[1]]),1,mean)
	if(do.plot) pairs(l[[3]][[2]], col=rainbow(length(unique(info)))[as.numeric(info)], pch=(1:length(unique(info)))[as.numeric(info)])

	return(to.return)
}

grp.to.d<-function(info)
{
	d<-matrix(0, ncol=length(info), nrow=length(info))
	for(i in 1:nrow(d))
	{
		d[i,]<-1*(info==info[i])
	}
	return(d)
}

