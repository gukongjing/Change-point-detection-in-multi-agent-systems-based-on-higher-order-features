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
  
  topo = function(xv){
    clus = vector()
    
    life = sapply(c(1:N), function(i){
      piece = xv[,i,]
      phom <- calculate_homology(piece)
      phom0 = as.data.frame(phom[which(phom[,1]==0),])
      phom0$death
    })
    
    # plot(life[1,],type = "l",ylim= c(0,2))
    # for (i in 2:(Nt-1)) {
    #   lines(life[i,])
    # }
    
    for (i in 1:length(choice)) {
      
      selife = life[choice[[i]],]
      
      # plot(selife[1,],type = "l",ylim= c(0,2))
      # for (i in 2:lim) {
      #   lines(selife[i,])
      # }
      
      y = e.agglo(t(selife))  #2013快   agglomerative hierarchical 
      pool = y$estimates
      CPfind =pool[which(pool%in%include)]
      clus[i] = Fvalue(CPreal,CPfind)
    }
    
    return(clus)
  }
  
  
  source('leader_follower_time_delay.R')
  source('netconstruction.R')
  
  library(TDAstats)
  library(ecp)
  
  Ntpool = c(20,40,60)   #数量
  kpool = c(2,5,6)
  clusall = matrix(0,nrow = 5,ncol = length(Ntpool))
  id = 1
  for (Nt in Ntpool) {
    gap = floor(Nt/10)
    choice = list(c((Nt/2-gap+1):(Nt/2+gap-1)),c((Nt/2-2*gap+1):(Nt/2+2*gap-1)),c((Nt/2-3*gap+1):(Nt/2+3*gap-1)),
                  c((Nt/2-4*gap+1):(Nt/2+4*gap-1)),c((Nt/2-5*gap+1):(Nt/2+5*gap-1)))
    
    uav = c(1:Nt)   #无人机索引
    kind1 = 3        # 网络类型
    kind2 = 3
    k1 = kpool[id]           #子群数目
    k2 = 8
    m = 3
    disturbnum = 400
    
    N1 = 150
    N2 = 100
    N = N1+N2
    sequence = c(1:(N1+N2))
    td = 1
    CPreal = N1        #原本
    
    # 生成参考网络
    g1 = netconstruction(Nt,kind1,disturbnum,k1,m)
    g1[which(g1>1)] = 1
    diag(g1)=0
    G = graph.adjacency(g1,mode="undirected")
    plot(G)
    
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
    G = graph.adjacency(g2,mode="undirected")
    plot(G)
    
    ## 轨迹
    R2 = 2
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
    
    clusall[,id] = topo(xv)
    id = id+1
  }
  
  final = as.vector(clusall)
  return(final)
})

stopCluster(cl)

D = do.call(rbind, simul)

# write.table(D,file = 'PtestP.txt', sep = ' ', row.names = FALSE, col.names = FALSE)   #需要重新跑了
D = as.matrix(read.table('PtestP.txt',sep = ' ', header = FALSE))

E = apply(D, 2, mean)

library(RColorBrewer)
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
color1 <- colorRampPalette(rev(brewer.pal(11,'Spectral')))(28)
color2 = c("#dc745c","#efaa6e","#f2c368","#b0bcd6","#6a84b5","#3b5899")

Variable = c('1/5','2/5','3/5','4/5','1')
# choice = list(c(1:2),c(1:10),c(1:39),c(11:20),c(20:29),c(30:39),c(11:30),c(11:39),c(38:39))
# labe = c('v1-2','v1-10','v1-39','v11-20','v20-29','v30-39','v11-30','v11-39','v38-39')
# Variable = as.factor(labe)
ggcover = data.frame(F1 = E, Proportion = rep(Variable,3),N = rep(as.factor(c('N = 20','N = 40','N = 60')),each = 5))
ggcover$Proportion <- factor(ggcover$Proportion, levels=Variable)

library(ggplot2)
p1 = ggplot(data = ggcover, mapping = aes(x = Proportion, y = F1,group = N,color= N,shape=N))+
  geom_point(size = 2.5) + geom_line(size = 0.8)+
  labs(title = "(b)",x = expression(paste('Proportion (',italic(l),')'))) +
  theme_classic()+
  theme(axis.title.x = element_text(size = 12), axis.text.x = element_text(angle = 0, hjust=1, size = 12),
        axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 12),
        plot.title = element_text(size=12,hjust=0.5),legend.position = 'none',
        legend.text = element_text(size = 12),legend.title = element_blank())+
  scale_color_manual(values=c("#dc745c","black",'#a9a9a9'))+   ##  #efaa6e
  ylim(limits = c(0,1))
# geom_vline(data = data.frame(x = c(2,35)), aes(xintercept = x), linetype = 'dashed', size = 0.5)
p1


#################################################画图################################################################
rm(list = ls())
source('leader_follower_time_delay.R')
source('netconstruction.R')

library(TDAstats)

Nt <- 40   #数量

uav = c(1:Nt)   #无人机索引
kind1 = 3        # 网络类型
kind2 = 3
k1 = 5           #子群数目
k2 = 8
m = 3
disturbnum = 400

N1 = 150
N2 = 100
N = N1+N2
sequence = c(1:(N1+N2))
td = 1
CPreal = N1        #原本

# 生成参考网络
g1 = netconstruction(Nt,kind1,disturbnum,k1,m)
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
R2 = 2
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

life = sapply(c(1:N), function(i){
  piece = xv[,i,]
  phom <- calculate_homology(piece)
  phom0 = as.data.frame(phom[which(phom[,1]==0),])
  phom0$death
})

life = life[c((Nt/2-11):(Nt/2+11)),]
d = dim(life)[1]
### 画Persistence
i=0
ggxg = life[1,]
for (i in 2:d) {
  ggxg = rbind(ggxg,life[i,])
}
# 
# # ggplot traj
ggdata <- data.frame(
  Homology_Class = rep(c(1:d),each = N),    #homology class
  Persistence = as.vector(t(ggxg)),
  Timestamp = rep(c(1:N),d)
)

p2 = ggplot() +  geom_point(data = ggdata, 
                            mapping = aes(x = Timestamp, y=Persistence,color = Homology_Class, group = Homology_Class),
                            size = 0.2) +
  labs(x=expression(italic(t)),y=expression(italic(p)^paste('[',italic(z),']'))) +
  geom_vline(xintercept=100, size = 0.5,color = 'black')+
  theme_classic()+
  theme(axis.title.x=element_text(size=12), axis.text.x=element_text(size=12),
        axis.title.y=element_text(size=12), axis.text.y=element_text(size=12), 
        title=element_text(size=12),legend.text = element_text(size = 12),legend.position = 'top')+
  # scale_colour_manual(values=c('black','grey75','dimgray'))+
  scale_colour_gradient(low = "#dc745c",high = "black",name=expression(paste('[',z,']')))
p2

ggsave(p2, file='plot\\persistence1.svg', width=3.5, height=3) 

