source("../dist_group.R")
source("../MEV Diffu.R")
source("../main_fun.R") 
library(ggplot2)
library(mclust)
library(cluster)
library(aricode)
set.seed(2024) 

#######################
load("Micheli_mouse.Rdata")
dat<-Micheli_mouse
cell.type=info=as.factor(Micheli_mouse.cellType)

load("all embeded dim 3.bin")
n_components<-3
dataset = "Micheli"

dist.power.list = c(0.5)
conn.prop.list=c(0.012777777)
diffu.factor.list = c(1.5)

source("../plot.R")