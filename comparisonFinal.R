#去掉KcpRS,persistence取3/5

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
  maxK = 5
  
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
  
  
  covering = function(N,CPreal,CPfind){
    l1 = length(CPreal)+1
    l2 = length(CPfind)+1
    CPreal0 = c(1,CPreal,N)
    CPfind0 = c(1,CPfind,N)
    
    Jaccard = sapply(c(1:l1),function(i){
      temp = sapply(c(1:l2), function(j){
        intersection = intersect(c(CPreal0[i]:CPreal0[i+1]),c(CPfind0[j]:CPfind0[j+1]))
        U = union(c(CPreal0[i]:CPreal0[i+1]),c(CPfind0[j]:CPfind0[j+1]))
        length(intersection)/length(U)
      })
      max(temp)*abs(CPreal0[i+1]-CPreal0[i]+1)
    })
    return(sum(Jaccard)/N)
  }
  
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
  
  datalist = list()
  datalist[[1]] = vg[,,1]
  datalist[[2]] = vg
  datalist[[3]] = xv

  clusall = matrix(0, nrow = 18, ncol = 3)
  
  # cons = 1
  
  for (cons in 1:3) {

    clus = vector()
    
    data1 = datalist[[cons]]
    if(cons>1){
      x = data1[,,1]
      for (i in 2:dim(data1)[3]) {
        x = rbind(x,data1[,,i])
      }
    }else{
      x = data1
    }
    
    x = t(x)
    
    library(bcp)
    #################1111111111111111111111111111
    res2 <- bcp(x)     #2015 
    temp = order(res2$posterior.prob,decreasing=T)[1:maxK]    #限定最大K
    pool = temp[which(res2$posterior.prob[temp]>=0.1)]
    CPfind = pool[which(pool%in%include)]
    clus[1] = Fvalue(CPreal,CPfind)
    
    
    library(ecp)
    ################222222222222222222222222222222222
    y = e.agglo(x)  #2013快   agglomerative hierarchical 
    pool = y$estimates
    CPfind = pool[which(pool%in%include)]
    clus[2] = Fvalue(CPreal,CPfind)
    
    
    ##################333333333333333333333333333333333
    y1 = e.divisive(X=x,k = maxK)   #2013中快 divisive hierarchical
    pool = y1$estimates
    CPfind = pool[which(pool%in%include)]
    clus[3] = Fvalue(CPreal,CPfind)
    
    
    ###################4444444444444444444444444444
    y1 = e.cp3o_delta(Z=x,K = 5)     #2017快   dynamic pruning
    pool = y1$estimates
    CPfind = pool[which(pool%in%include)]
    clus[4] = Fvalue(CPreal,CPfind)
    
    
    #####################55555555555555555555555555
    y = ks.cp3o_delta(Z=x,K=5)   #2017快   pruning   , minsize=30, verbose=FALSE
    pool = y$estimates
    CPfind = pool[which(pool%in%include)]
    clus[5] = Fvalue(CPreal,CPfind)
    
    #####################666666666666666666666666666
    y2 = kcpa(x,L=5,C=5)   #2019中   kernel    
    CPfind = y2[which(y2%in%include)]
    clus[6] = Fvalue(CPreal,CPfind)
    
    
    ###########################################################
    if(cons==1){
      distmat =  lapply(c(1:N), function(i){
        piece = rbind(data1[,i],rep(0,Nt))
        as.matrix(dist(t(piece)))
      })
    }else{
      distmat =  lapply(c(1:N), function(i){
        piece = data1[,i,]
        as.matrix(dist(piece))
      })
    }
    library(Rfast)
    # em = max(sapply(c(1:N),function(i){
    #   temp = distmat[[i]]
    #   diag(temp)= max(temp)
    #   max(colMins(temp,value = TRUE))
    # }))
    em = mean(sapply(c(1:N),function(i){
      temp = distmat[[i]]
      diag(temp)= max(temp)
      max(colMins(temp,value = TRUE))
    }))
    e = seq(0,em,em/(Nt))
    
    
    bettimat = sapply(c(1:N), function(i){
      phom <- calculate_homology(distmat[[i]],format = "distmat")
      phom = phom[which(phom[,1]==0),]
      betti = sapply(e, function(e){
        length(which(phom[,3]>e))
      })
      betti
    })
    
    selife = bettimat    #[1:Nt,]
    
    library(bcp)
    ##############################0000000000000000000000000
    res2 <- bcp(t(selife))     #2015
    temp = order(res2$posterior.prob,decreasing=T)[1:maxK]    #限定最大K
    pool = temp[which(res2$posterior.prob[temp]>=0.5)]
    CPfind =pool[which(pool%in%include)]
    clus[7] = Fvalue(CPreal,CPfind)
    
    #########################111111111111111111111111
    y = e.agglo(t(selife))  #2013快   agglomerative hierarchical 
    #, alpha=1, penalty=function(cps){0}
    pool = y$estimates
    CPfind =pool[which(pool%in%include)]
    clus[8] = Fvalue(CPreal,CPfind)
    
    ##########################2222222222222222222222222
    y1 = e.divisive(X=t(selife),k = maxK)   #2013中快 divisive hierarchical
    pool = y1$estimates
    CPfind = pool[which(pool%in%include)]
    clus[9] = Fvalue(CPreal,CPfind)
    
    ########################33333333333333333333333333
    y1 = e.cp3o_delta(Z=t(selife),K = 5)     #2017快   dynamic pruning
    pool = y1$estimates
    CPfind = pool[which(pool%in%include)]
    clus[10] = Fvalue(CPreal,CPfind)
    
    ########################444444444444444444444444444
    y = ks.cp3o_delta(Z=t(selife),K = 5)   #2017快   pruning
    pool = y$estimates
    pool = pool[which(pool%in%include)]
    CPfind =pool[which(pool%in%include)]
    clus[11] = Fvalue(CPreal,CPfind)
    
    #######################5555555555555555555555555555
    y2 = kcpa(t(selife), L=5, C = 5)   #2019中   kernel
    pool = y2[which(y2%in%include)]
    CPfind =pool[which(pool%in%include)]
    clus[12] = Fvalue(CPreal,CPfind)

    
    ##########################################################
    # life = sapply(c(1:N), function(i){
    #   phom <- calculate_homology(distmat[[i]],dim = 0,format = "distmat")
    #   phom[,3]
    # })
    
    life = sapply(c(1:N), function(i){
      phom = calculate_homology(distmat[[i]],format = "distmat")
      phom0 = as.data.frame(phom[which(phom[,1]==0),])
      phom0$death
    })
    # plot(life[1,],type = "l",ylim= c(0,4))
    # for (i in 2:(Nt-1)) {
    #   lines(life[i,])
    # }
    
    selife = life[(Nt/2-11):(Nt/2+11),]
    
    
    
    library(bcp)
    ##############################999999999999999999999999999
    res2 <- bcp(t(selife))     #2015
    temp = order(res2$posterior.prob,decreasing=T)[1:maxK]    #限定最大K
    pool = temp[which(res2$posterior.prob[temp]>=0.5)]
    CPfind =pool[which(pool%in%include)]
    clus[13] = Fvalue(CPreal,CPfind)
    
    #########################0000000000000000000000000000000
    y = e.agglo(t(selife))  #2013快   agglomerative hierarchical 
    #, alpha=1, penalty=function(cps){0}
    pool = y$estimates
    CPfind =pool[which(pool%in%include)]
    clus[14] = Fvalue(CPreal,CPfind)
    
    ##########################1111111111111111111111111111111
    y1 = e.divisive(X=t(selife),k = maxK)   #2013中快 divisive hierarchical
    pool = y1$estimates
    CPfind = pool[which(pool%in%include)]
    clus[15] = Fvalue(CPreal,CPfind)
    
    ########################222222222222222222222222222222222222
    y1 = e.cp3o_delta(Z=t(selife),K = 5)     #2017快   dynamic pruning
    pool = y1$estimates
    CPfind = pool[which(pool%in%include)]
    clus[16] = Fvalue(CPreal,CPfind)
    
    ########################333333333333333333333333333333
    y = ks.cp3o_delta(Z=t(selife),K = 5)   #2017快   pruning
    pool = y$estimates
    pool = pool[which(pool%in%include)]
    CPfind =pool[which(pool%in%include)]
    clus[17] = Fvalue(CPreal,CPfind)
    
    #######################44444444444444444444444444444444
    y2 = kcpa(t(selife), L=5, C = 5)   #2019中   kernel
    pool = y2[which(y2%in%include)]
    CPfind =pool[which(pool%in%include)]
    clus[18] = Fvalue(CPreal,CPfind)
    
    clusall[,cons] = clus
  }

  final = as.vector(clusall)
  
  return(final)
})

stopCluster(cl)

D = do.call(rbind, simul)
# write.table(D,file = 'simulComparison1.txt', sep = ' ', row.names = FALSE, col.names = FALSE)   #simulComparison0
D = as.matrix(read.table('simulComparison1.txt',sep = ' ', header = FALSE))

E = apply(D, 2, mean)
ta = E
dim(ta) = c(6,9)
# ta = t(ta)
ta = round(ta,3)
library(xlsx)
write.xlsx(x = ta, file = "1.xlsx", row.names = FALSE)

clus = E

library(ggplot2)
library(RColorBrewer)
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
color1 <- colorRampPalette(rev(brewer.pal(11,'Spectral')))(28)
color2 = c("#dc745c","#efaa6e","#f2c368","#b0bcd6","#6a84b5","#3b5899")

# 画图AUC-Netdens
method = c('BCP','E-Agglo','E-divisive','E-CP3O','KS-CP3O','KCPA')

ggclus = data.frame(F1 = clus, 
                    Method = rep(method,9),
                    Datatype = rep(c('v1','v','xv'),each = 18),
                    face = rep(rep(c('CPD Method','Betti Number','Persistence'),each = 6),3))
ggclus$Method <- factor(ggclus$Method, levels=method)
ggclus$Datatype <- factor(ggclus$Datatype, levels=c('v1','v','xv'))
ggclus$face <- factor(ggclus$face, levels=c('CPD Method','Betti Number','Persistence'))

p2 = ggplot(ggclus,aes(x = Method, y = F1,fill = Datatype))+
  geom_bar(stat = 'identity', 
           #柱状图位置并排:
           position = 'dodge', #使用position=position_dodge(width=0.9),可使组内柱子间隔,自行试一下。
           width = 0.8,      #设置柱子宽度,使变量之间分开
           color='black')+        
  labs(x=NULL)+ 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 65, hjust=1, size = 12),
        axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 12),
        legend.position = 'none')+ 
  scale_fill_manual(values=color2)+
  facet_grid(Datatype~face)+   #,scales= "free",space= "free"  
  theme(strip.text.x = element_text(size = 12),strip.text.y = element_text(size = 12))
p2





