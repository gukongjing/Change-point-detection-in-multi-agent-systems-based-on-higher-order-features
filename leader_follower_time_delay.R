# this content aim to explore collective behavior under different connection matrix
# rm(list = ls())
# 
# library(network)
# library(ggplot2)
# library(RColorBrewer)
# # select the observed data
# # layout(matrix(c(1,2),nrow = 1, byrow = TRUE))
# 
# Nt = 150
# kind = 6
# disturbnum = 2000
# k=5
# m = 3
# td=1
# N = 200


# 壳
traj = function(k,td,Nt,N,uav,A,x0,v0,R){
  library(ggplot2)
  library(network)
  library(GGally)
  library(sna)
  library(RColorBrewer)
  # td =1 
  t <- 0.1
  gnum <- Nt/k
  # R <- 1
  td = td     # time delay
  sdg = 0    # 高斯噪声的方差   
  
  
  
  # known & unknown parameters
  f = 2
  alpha = 2
  beta = 1
  sigma = 0.5
  
  c <- rep(0,Nt*(N-1)*2)
  w <- rep(0,Nt*(N-1)*2)
  u <- rep(0,Nt*(N-1)*2)
  
  dim(c) <- c(Nt,N-1,2)
  dim(w) <- c(Nt,N-1,2)
  dim(u) <- c(Nt,N-1,2)
  
  
  # initial values
  x <- rep(0,Nt*N*2)
  v <- rep(0,Nt*N*2)
  xg <- rep(0,Nt*N*2)    # 加高斯
  vg <- rep(0,Nt*N*2)    # 加高斯
  
  dim(x) <- c(Nt,N,2)
  dim(v) <- c(Nt,N,2)
  dim(xg) <- c(Nt,N,2)    # 高斯版
  dim(vg) <- c(Nt,N,2)    # 高斯版
  
  x[,(1:(1+td)),] <- x0
  v[,(1:(1+td)),] <- v0
  
  # model
  for(m in (1+td):(N-1)){
    for(j in 1:Nt){
      temp0 <- (x[,m-td,1]-x[j,m,1])^2+(x[,m-td,2]-x[j,m,2])^2
      temp1 <- 1+temp0
      temp2 <- 1/(temp1^f)*A[j,]
      temp31 <- v[,m-td,1]-v[j,m,1]
      temp32 <- v[,m-td,2]-v[j,m,2]
      temp3 <- cbind(temp31,temp32)
      c[j,m,] <- temp2%*%temp3
      
      xm <- colMeans(x[,m-td,],1)
      temp4 <- xm-x[j,m,]
      temp4_ab <- sqrt(sum(temp4^2))
      w[j,m,] <- (1-R/temp4_ab)*temp4
      
      for(i in 1:Nt){
        temp5 <- x[j,m,]-x[i,m-td,]
        if(sqrt(sum(temp5^2)) < 2*sin(pi/Nt)&&sqrt(sum(temp5^2)) != 0){
          temp5_ab <- sqrt(sum(temp5^2))
          temp6 <- temp5/temp5_ab-temp5
          u[j,m,] <- temp6*sin(pi/Nt)/2+u[j,m,]
        }
      }
    }
    
    v[,m+1,] <- (alpha*c[,m,]+beta*u[,m,]+sigma*w[,m,])*t+v[,m,]   #
    x[,m+1,] <- x[,m,]+(v[,m+1,]+v[,m,])*t/2
    
  }
  
  # layout(matrix(c(1,2,3,4),nr=2,byrow = TRUE))
  # layout(1)
  # 
  # plot(x[1,25:N,1],x[1,25:N,2],xlim=c(0,22),ylim=c(0,25),xlab = "x-coordinate",ylab = "y-coordinate",type = 'l',
  #      cex.lab=2,cex.axis=2)   # cex.main=1.5,
  # for (z in 2:Nt) {
  #   lines(x[z,25:N,1],x[z,25:N,2])
  # }
  
  # 添加高斯噪声
  i=0
  for (i in 1:N) {
    temp = sapply(uav, function(uav){
      x[uav,i,] + rnorm(2, mean = 0, sd = sdg)
    })
    xg[,i,] = t(temp)
  }
  
  i=0
  for (i in 1:N) {
    temp = sapply(uav, function(uav){
      v[uav,i,] + rnorm(2, mean = 0, sd = sdg)
    })
    vg[,i,] = t(temp)
  }
  
  # 加速度
  a <- rep(0,Nt*N*2)
  ag <- rep(0,Nt*N*2)    # 加高斯
  dim(a) <- c(Nt,N,2)
  dim(ag) <- c(Nt,N,2)
  
  a[,1:(N-1),] = v[,2:N,]-v[,1:(N-1),]
  a[,N,] = a[,N-1,]
  i = 0
  for (i in 1:N) {
    temp = sapply(uav, function(uav){
      a[uav,i,] + rnorm(2, mean = 0, sd = sdg)
    })
    ag[,i,] = t(temp)
  }
  
  
  # # 归一化
  # i = 1
  # for (i in 1:2) {
  #   xg[,,i] = (xg[,,i]-mean(xg[,,i]))/sqrt(sum((xg[,,i]-mean(xg[,,i]))^2)/(length(xg[,,i])-1))  # z-score 标准化
  # }
  # i = 1
  # for (i in 1:2) {
  #   vg[,,i] = (vg[,,i]-mean(vg[,,i]))/sqrt(sum((vg[,,i]-mean(vg[,,i]))^2)/(length(vg[,,i])-1))  # z-score 标准化
  # }
  # i = 1
  # for (i in 1:2) {
  #   ag[,,i] = (ag[,,i]-mean(ag[,,i]))/sqrt(sum((ag[,,i]-mean(ag[,,i]))^2)/(length(ag[,,i])-1))  # z-score 标准化
  # }

  resu= list(xg,vg,ag)
  return(resu)
}

# i=0
# ggxg = xg[1,25:(N-25),]
# for (i in 2:Nt) {
#   ggxg = rbind(ggxg,xg[i,25:(N-25),])
# }
# # 
# # # ggplot traj
# ggdata <- data.frame(
#   Trajectory=paste('Module',rep(c(1:k), each = (N-24-25)*gnum)),
#   poin = factor(rep(c(1:Nt),each = N-24-25)),
#   traj.x1=ggxg[,1],
#   traj.x2=ggxg[,2]
# )
# 
# seg = 50
# segpoint = xg[,seq(25,N,seg),]
# 
# ggseg = data.frame(
#   traj.x1 = as.vector(segpoint[,,1]),    #-c(1,2,3)
#   traj.x2 = as.vector(segpoint[,,2]),
#   piece = rep(c(1:((N-24)/seg)),each = Nt)
# )
# 
# meany = apply(segpoint[,,2], 2, mean)    #每组数据的平均y值
# 
# xseg = sapply(1:((N-24)/seg),function(u){
#   pickup = which(segpoint[,u,2] > meany[u])
#   pickbottom = which(segpoint[,u,2] < meany[u])
#   xup = mean(segpoint[pickup,u,1])
#   xbottom = mean(segpoint[pickbottom,u,1])
#   c(xup,xbottom)
# })
# yb = apply(segpoint[,,2], 2, min)
# yu = apply(segpoint[,,2], 2, max)
# 
# ggline = data.frame(
#   traj.xb = xseg[2,],    #下方节点x值
#   traj.xu = xseg[1,],    #上方节点x值
#   traj.yb = apply(segpoint[,,2], 2, min),
#   traj.yu = apply(segpoint[,,2], 2, max),
#   piece = c(1:((N-24)/seg))
# )
# display.brewer.pal(9,"Greys")
# grey = brewer.pal(9,"Greys")
# mycolors<- grey[c(5,6,7,8,9)]
# ggplot() +  geom_path(data = ggdata, mapping = aes(x = traj.x1, y=traj.x2,color = Trajectory,group = poin),size = 1) +
#   annotate('segment', x=xseg[1,], xend=xseg[2,], y=yu+0.5, yend=yb-0.5,fill = 'grey',color = 'black',alpha = 0.6,size = 15)+
#   # stat_ellipse(data = ggseg,aes(x = traj.x1, y=traj.x2, group  = piece),geom = "polygon",
#   #              color = 'grey',alpha = 0.2,level = 0.90, show.legend = F) +
#   # geom_point(data = ggseg,aes(x = traj.x1, y=traj.x2),color = 'grey',alpha = 0.2,size = 10) +
#   geom_point(data = ggseg,aes(x = traj.x1, y=traj.x2),color = '#FFC83F',alpha = 1,size = 1) +
#   theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
#         axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(),
#         title=element_text(size=15),legend.text = element_text(size = 15),legend.position = c(0.85,0.25))+
#   # scale_colour_manual(values=c('black','grey75','dimgray'))+
#   scale_color_manual(values = mycolors)+
#   theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank()) +  # 去掉背景 
#   guides(color=FALSE)
# 
# plot(xg[,1,])
# plot(xg[,300,])
# plot(xg[1,25:324,1],xg[1,25:324,2],type = "l")
# for (i in 2:Nt) {
#   lines(xg[i,25:324,1],xg[i,25:324,2])
# }

# #保存数据
# i = 0
# for (i in 1:Nt) {
#   name = paste("data/",Nt,"uavs/uav",i,".txt", sep = "")
#   write.table(cbind(xg[i,,],vg[i,,]), file = name, sep = " ",
#               col.names = c('location_x','location_y','velocity_x','velocity_y'))
# }
# write.table(A,file = "data/150uavs/network.txt",sep = " ")
# # 画图
# 
# i=0
# ggvg = vg[1,,]
# for (i in 2:Nt) {
#   ggvg = rbind(ggvg,vg[i,,])
# }
# # 
# # # ggplot
# ggdatav <- data.frame(
#   cluster=factor(rep(c(1:k), each = N*gnum)),
#   poin = factor(rep(c(1:Nt),each = N)),
#   color = c("red"),
#   traj.x=ggvg[,1],
#   traj.y=ggvg[,2]
# )
# library(ggplot2)
# ggplot(data = ggdatav, mapping = aes(x = traj.x, y=traj.y,color = cluster,group = cluster)) +
#   ggtitle('Trajectories')+geom_path(aes(group = poin),size = 1) + scale_colour_brewer(palette="Blues")
# #+ scale_colour_manual(values=c("#84CEFC","#ACDAE4","#5C9EA4","#04FEFC","#4CD2CC"))
# 
# i=0
# ggxg = xg[1,,]
# for (i in 2:Nt) {
#   ggxg = rbind(ggxg,xg[i,,])
# }
# ggdata <- data.frame(
#   cluster=factor(rep(c(1:k), each = N*gnum)),
#   poin = factor(rep(c(1:Nt),each = N)),
#   color = c("red"),
#   traj.x=ggxg[,1],
#   traj.y=ggxg[,2]
# )
# library(ggplot2)
# ggplot(data = ggdata, mapping = aes(x = traj.x, y=traj.y,color = cluster,group = cluster)) +
#   ggtitle('Trajectories')+geom_path(aes(group = poin),size = 1) + scale_colour_brewer(palette="Blues")
# 


