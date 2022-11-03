#缺失比例按照总数据设置

rm(list = ls())
library(parallel)
#Calculate the number of cores检查电脑当前可用核数
no_cores<-detectCores() - 2		#F-物理CPU核心数/T-逻辑CPU核心数 logical=F

#Initiate cluster发起集群,同时创建数个R进行并行计算
#只是创建待用的核,而不是并行运算环境
cl<-makeCluster(no_cores)

rept = 50
r = c(1:rept)

#现只需要使用并行化版本的lapply,parLapply就可以
simul = parallel::parLapply(cl, r,function(r){
  
  Fvalue = function(CPreal,CPfind,error=5){
    l1 = length(CPreal)
    l2 = length(CPfind)
    TPvec = sapply(c(1:l1), function(i){
      ii = which(CPfind%in%c((CPreal[i]-5):(CPreal[i]+5)))
      length(ii)
    })
    TPvec[which(TPvec>0)]=1
    TP = sum(TPvec)
    P = TP/l2
    R = TP/l1
    if(l2==0 | (P+R)==0){P=10^(-5)}
    return(2*(P*R)/(P+R))
  }
  
  
  library(Rfast)
  library(TDAstats)
  library(corrplot)
  library(ecp)
  
  source('leader_follower_time_delay.R')
  source('netconstruction.R')
  
  
  Ntpool = c(20,40,60)   #数量
  kpool = c(2,5,6)
  clusall = matrix(0,nrow = 6,ncol = length(Ntpool))
  id = 1
  for (Nt in Ntpool) {
    uav = c(1:Nt)   #无人机索引
    kind = 3        # 网络类型
    k1 = kpool[id]           #子群数目
    k2 = 8
    m = 3
    disturbnum = 400
    
    N1 = 150
    N2 = 100
    N = N1+N2
    sequence = c(1:(N1+N2))
    td = 1
    
    # 生成参考网络
    g1 = netconstruction(Nt,kind,disturbnum,k1,m)
    g1[which(g1>1)] = 1
    diag(g1)=0
    # G = graph.adjacency(g1,mode="undirected")
    # plot(G)
    
    # 生成轨迹信息
    R1 = 1
    x0 <- runif(Nt*2*(1+td),min = 0,max = 10)
    v0 <- runif(Nt*2*(1+td),min = 0,max = 1)
    dim(x0) = c(Nt,(td+1),2)
    dim(v0) = c(Nt,(1+td),2)
    xvg1 = traj(k1,td,Nt,N1,uav,g1,x0,v0,R1)
    
    
    ### second stage
    # 生成参考网络
    g2 = netconstruction(Nt,4,disturbnum,k2,m)
    g2[which(g2>1)] = 1
    diag(g2)=0
    # G = graph.adjacency(g2,mode="undirected")
    # plot(G)
    
    ## 轨迹
    R2 = 1
    x0 = xvg1[[1]][,((N1-td):N1),]
    v0 = xvg1[[2]][,((N1-td):N1),]
    xvg2 = traj(k2,td,Nt,N2+td+1,uav,g2,x0,v0,R2)
    
    xg = rep(0,Nt*N*2)
    dim(xg) = c(Nt,N,2)
    xg[,1:N1,] = xvg1[[1]]
    xg[,(N1-td):N,] = xvg2[[1]]
    
    
    vg = rep(0,Nt*N*2)
    dim(vg) = c(Nt,N,2)
    vg[,1:N1,] = xvg1[[2]]
    vg[,(N1-td):N,] = xvg2[[2]]
    
    start = 50         #去掉开头不稳定的数据
    xg = xg[,(start+1):N,]
    vg = vg[,(start+1):N,]
    
    CPreal = N1        #原本
    CPreal = CPreal-start    #去掉一截数据后
    N = N-start
    include = c(2:(N-2))    #排除边界点
    
    
    i = 0
    for (i in 1:2) {
      xg[,,i] = (xg[,,i]-mean(xg[,,i]))/sqrt(sum((xg[,,i]-mean(xg[,,i]))^2)/(length(xg[,,i])-1))
      vg[,,i] = (vg[,,i]-mean(vg[,,i]))/sqrt(sum((vg[,,i]-mean(vg[,,i]))^2)/(length(vg[,,i])-1))
    }
    
    # layout(1)
    plot(xg[1,,1],xg[1,,2],type = "l")
    for (i in 2:Nt) {
      lines(xg[i,,1],xg[i,,2])
    }
    
    #### prepare data xv
    xv = rep(0,Nt*N*4)
    dim(xv) = c(Nt,N,4)
    xv[,,1:2] = xg
    xv[,,3:4] = vg
    
    distmat =  lapply(c(1:N), function(i){
      piece = xv[,i,]
      as.matrix(dist(piece))
    })
    
    library(Rfast)
    em = mean(sapply(c(1:N),function(i){
      temp = distmat[[i]]
      diag(temp)= max(temp)
      max(colMins(temp,value = TRUE))
    }))
    
    res = c(Nt/2,Nt/2,c(1,3,5,7)*Nt)
    clus = vector()
    
    for (rr in 1:length(res)) {
      e = seq(0,em,em/res[rr])
      
      bettimat = sapply(c(1:N), function(i){
        phom <- calculate_homology(distmat[[i]],format = "distmat")
        phom = phom[which(phom[,1]==0),]
        betti = sapply(e, function(e){
          length(which(phom[,3]>e))
        })
        betti
      })
      
      y = e.agglo(t(bettimat))  #2013快   agglomerative hierarchical 
      #, alpha=1, penalty=function(cps){0}
      pool = y$estimates
      CPfind =pool[which(pool%in%include)]
      clus[rr] = Fvalue(CPreal,CPfind)
    }
    
    clusall[,id] = clus
    id = id+1
  }
  
  
  final = as.vector(clusall)
  
  return(final)
})
stopCluster(cl)

D = do.call(rbind, simul)

# write.table(D,file = 'PtestB.txt', sep = ' ', row.names = FALSE, col.names = FALSE)
D = as.matrix(read.table('PtestB.txt',sep = ' ', header = FALSE))

E = apply(D, 2, mean)

library(RColorBrewer)
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
color1 <- colorRampPalette(rev(brewer.pal(11,'Spectral')))(28)
color2 = c("#dc745c","#efaa6e","#f2c368","#b0bcd6","#6a84b5","#3b5899")

# 画图AUC-Netdens
type = as.factor(c('N/4','N/2','N','3N','5N','7N'))


ggcover = data.frame(F1 = E, Resolution = rep(type,3),N = rep(as.factor(c('N = 20','N = 40','N = 60')),each = 6))
ggcover$Resolution <- factor(ggcover$Resolution, levels=type)

library(ggplot2)
p3 = ggplot(data = ggcover, mapping = aes(x = Resolution, y = F1,group = N, color = N,shape = N))+
  geom_point(size = 2.5) + geom_line(size = 0.8)+
  theme_classic()+
  labs(title="(a)",x = expression(paste('Resolution (',italic(M),')'))) +
  theme(axis.title.x = element_text(size = 12), axis.text.x = element_text(angle = 0, hjust=1, size = 12),
        axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 12),
        plot.title = element_text(size=12,hjust=0.5),legend.position = c(0.3,0.3),
        legend.text = element_text(size = 12),legend.title = element_blank())+
  scale_color_manual(values=c("#dc745c","black",'#a9a9a9'))+
  ylim(limits = c(0,1))
# geom_vline(data = data.frame(x = c(2,35)), aes(xintercept = x), linetype = 'dashed', size = 0.5)
p3


library(cowplot)
gg <- ggdraw() +
  draw_plot(p3, 0,0, 0.48,1) + draw_plot(p1, 0.52,0,0.48,1)
print(gg)



####################################################################画图######################################################
rm(list = ls())
Nt <- 40   #数量
uav = c(1:Nt)   #无人机索引
kind = 3        # 网络类型
k1 = 5           #子群数目
k2 = 8
m = 3
disturbnum = 400

N1 = 150
N2 = 100
N = N1+N2
sequence = c(1:(N1+N2))
td = 1

source('leader_follower_time_delay.R')
source('netconstruction.R')

library(TDAstats)

# 生成参考网络
g1 = netconstruction(Nt,kind,disturbnum,k1,m)
g1[which(g1>1)] = 1
diag(g1)=0
# G = graph.adjacency(g1,mode="undirected")
# plot(G)

# 生成轨迹信息
R1 = 1
x0 <- runif(Nt*2*(1+td),min = 0,max = 10)
v0 <- runif(Nt*2*(1+td),min = 0,max = 1)
dim(x0) = c(Nt,(td+1),2)
dim(v0) = c(Nt,(1+td),2)
xvg1 = traj(k1,td,Nt,N1,uav,g1,x0,v0,R1)


### second stage
# 生成参考网络
g2 = netconstruction(Nt,4,disturbnum,k2,m)
g2[which(g2>1)] = 1
diag(g2)=0
# G = graph.adjacency(g2,mode="undirected")
# plot(G)

## 轨迹
R2 = 1
x0 = xvg1[[1]][,((N1-td):N1),]
v0 = xvg1[[2]][,((N1-td):N1),]
xvg2 = traj(k2,td,Nt,N2+td+1,uav,g2,x0,v0,R2)

xg = rep(0,Nt*N*2)
dim(xg) = c(Nt,N,2)
xg[,1:N1,] = xvg1[[1]]
xg[,(N1-td):N,] = xvg2[[1]]


vg = rep(0,Nt*N*2)
dim(vg) = c(Nt,N,2)
vg[,1:N1,] = xvg1[[2]]
vg[,(N1-td):N,] = xvg2[[2]]

start = 50         #去掉开头不稳定的数据
xg = xg[,(start+1):N,]
vg = vg[,(start+1):N,]

CPreal = N1        #原本
CPreal = CPreal-start    #去掉一截数据后
N = N-start
include = c(2:(N-2))    #排除边界点


i = 0
for (i in 1:2) {
  xg[,,i] = (xg[,,i]-mean(xg[,,i]))/sqrt(sum((xg[,,i]-mean(xg[,,i]))^2)/(length(xg[,,i])-1))
  vg[,,i] = (vg[,,i]-mean(vg[,,i]))/sqrt(sum((vg[,,i]-mean(vg[,,i]))^2)/(length(vg[,,i])-1))
}

# layout(1)
plot(xg[1,,1],xg[1,,2],type = "l")
for (i in 2:Nt) {
  lines(xg[i,,1],xg[i,,2])
}

#### prepare data xv
xv = rep(0,Nt*N*4)
dim(xv) = c(Nt,N,4)
xv[,,1:2] = xg
xv[,,3:4] = vg

distmat =  lapply(c(1:N), function(i){
  piece = xv[,i,]
  as.matrix(dist(piece))
})

library(Rfast)
em = max(sapply(c(1:N),function(i){
  temp = distmat[[i]]
  diag(temp)= max(temp)
  max(colMins(temp,value = TRUE))
}))

e = seq(0,em,em/80)

bettimat = sapply(c(1:N), function(i){
  phom <- calculate_homology(distmat[[i]],format = "distmat")
  phom = phom[which(phom[,1]==0),]
  betti = sapply(e, function(e){
    length(which(phom[,3]>e))
  })
  betti
})

library(reshape2)
melted_cormat <- melt(bettimat[,76:125], na.rm = TRUE)
# Heatmap
library(ggplot2)
p4 = ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  geom_vline(xintercept=25.5, size = 0.5,color = 'black')+
  labs(x=expression(italic(t)),y=expression(italic(m))) +
  scale_fill_gradient2(low = "grey", high = "black",
                       limit = c(0,max(bettimat)), 
                       name=expression(italic(b^epsilon[m]))) +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        title=element_text(size=12),
        legend.text = element_text(size = 12),
        legend.position = 'top')
p4


ggsave(p4, file='plot\\betti1.svg', width=3.5, height=3) 

