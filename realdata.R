##################################################三维空间 高维数据################################################
rm(list = ls())
bow = read.csv("data/MSRC12/P1_1_7_p21.csv", header = FALSE,sep = ' ')
# throw = read.csv("data/MSRC12/P1_1_8_p21.csv", header = FALSE,sep = ' ')

first = bow[,-1]    #931:1030
first = first[-which(rowSums(first)==0),]
# second = throw[,-1]
# second = second[-which(rowSums(second)==0),]

# data = rbind(first,second[,])    #1:100
data = first
data = data[,-c(1:20)*4]

N = dim(data)[1]
Nt = dim(data)[2]/3

CPreal = c(1,55,97,160,187,250,277,346,367,430,451,520,538,607,625,688,709,769,796,853,883,943,N) #manual annotation

library(ecp)
y2 <- e.agglo(as.matrix(data))     #原方法差一点点
pool1 = y2$estimates
pool1 = intersect(union(pool1,c(1,N)),c(1:N))

N = N-1
distmat =  lapply(c(1:N), function(i){
  piece = as.matrix(data[i,])
  dim(piece) = c(3,20)
  as.matrix(dist(t(piece)))
})

library(Rfast)
em = mean(sapply(c(1:N),function(i){
  temp = distmat[[i]]
  diag(temp)= max(temp)
  max(colMins(temp,value = TRUE))
}))
e = seq(0,em,em/15)

library(TDAstats)
bettimat = sapply(c(1:N), function(i){
  phom <- calculate_homology(distmat[[i]],dim = 0,format = "distmat")
  betti = sapply(e, function(e){
    length(which(phom[,3]>e))
  })
  betti
})

library(corrplot)
# corrplot(bettimat[,1:50]/max(bettimat[,1:50]))
y <- e.agglo(t(bettimat))     
pool2 = y$estimates
pool2 = intersect(union(pool2,c(1,N)),c(1:N))

life = sapply(c(1:N), function(i){
  phom <- calculate_homology(distmat[[i]],dim = 0,format = "distmat")
  phom[,3]
})
plot(life[1,],type = "l",ylim= c(0,1))
for (i in 2:(Nt-1)) {
  lines(life[i,])
}

lim = c((Nt/2-5):(Nt/2+5))
selife = life[lim,]
plot(selife[1,],type = "l",ylim= c(0,1))
for (i in 2:length(lim)) {
  lines(selife[i,])
}
selife = life[1:(Nt-1),]
y <- e.agglo(t(selife))     
pool3 = y$estimates
pool3 = intersect(union(pool3,c(1,N)),c(1:N))

# write.table(pool1,file = 'pool1.txt', sep = ' ', row.names = FALSE, col.names = FALSE)
# write.table(pool2,file = 'pool2.txt', sep = ' ', row.names = FALSE, col.names = FALSE)
# write.table(pool3,file = 'pool3.txt', sep = ' ', row.names = FALSE, col.names = FALSE)
pool1 = as.vector(as.matrix(read.table('pool1.txt',sep = ' ', header = FALSE)))
pool2 = as.vector(as.matrix(read.table('pool2.txt',sep = ' ', header = FALSE)))
pool3 = as.vector(as.matrix(read.table('pool3.txt',sep = ' ', header = FALSE)))


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

clus1 = Fvalue(CPreal,pool1)
clus2 = Fvalue(CPreal,pool2)
clus3 = Fvalue(CPreal,pool3)


ggseg = data.frame(CP = c(pool1,pool2,pool3),
                   Method = c(rep('E.Agglo',length(pool1)),rep('The Betti Number',length(pool2)),rep('The Persistence',length(pool3))),
                   y = c(rep(0.6,length(pool1)),rep(1.8,length(pool2)),rep(1.2,length(pool3))))
ggseg$Method <- factor(ggseg$Method, levels=c('The Betti Number','The Persistence','E.Agglo'))

dots = CPreal
df <- data.frame(
  xmin = dots[1:(length(dots)-1)],
  xmax = dots[2:length(dots)],
  ymin = rep(0.47,length(dots)-1),
  ymax = rep(1.93,length(dots)-1),
  Action = rep(c('Bow','Transition'),11)
)

library(ggsci)
library(ggplot2)
library(RColorBrewer)

display.brewer.pal(9,"Oranges")

c1 = brewer.pal(9,"Oranges")[7]
c2 = brewer.pal(9,"Greys")[8]
c3 = brewer.pal(9,"Greys")[9]
c4 = brewer.pal(9,"Oranges")[2]
c5 = brewer.pal(9,"Oranges")[4]
ggplot()+
  geom_rect(data = df, mapping = aes(xmin = xmin, xmax = xmax,ymin = ymin,ymax=ymax, fill = Action), alpha = 0.6)+
  geom_point(data=ggseg,mapping=aes(x=CP,y=y,group=Method,color = Method),shape = 3,size = 7)+
  geom_line(data=ggseg,mapping=aes(x=CP,y=y,group=Method,color=Method),size = 0.8)+
  annotate('text',x=-145,y=0.6,label='F1 = 0.21', size=4.63)+
  annotate('text',x=-145,y=1.2,label='F1 = 0.43', size=4.63)+
  annotate('text',x=-145,y=1.8,label='F1 = 0.64', size=4.63)+
  annotate("segment", x = -170, xend = -170, y = 1, yend = 1.2,
           size = 1.5, color = 'white',alpha = 0)+       #透明线  加画图区域用
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=13), axis.ticks.x=element_blank(),
        axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        legend.title = element_blank(),legend.text = element_text(size = 13),
        legend.position = 'right')+
  scale_color_manual(values=c(c1,'#4472C4','black'))+
  scale_fill_manual(values=c(c4,c5))+
  scale_y_continuous(limits = c(0.47, 1.93))+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank())   # 去掉背景 

##################################################画图#########################
gnet = matrix(0, nrow = Nt+2, ncol = Nt+2)
gnet[1,2] = 1
gnet[2,3] = 1
gnet[3,4] = 1
gnet[1,13] = 1
gnet[13,14] = 1
gnet[14,15] = 1
gnet[15,16] = 1
gnet[1,17] = 1
gnet[17,18] = 1
gnet[18,19] = 1
gnet[19,20] = 1
gnet[3,5] = 1
gnet[5,6] = 1
gnet[6,7] = 1
gnet[7,8] = 1
gnet[3,9] = 1
gnet[9,10] = 1
gnet[10,11] = 1
gnet[11,12] = 1
gnet[lower.tri(gnet)] <- t(gnet)[lower.tri(gnet)]
diag(gnet) = 0

library(GGally)
library(sna)

g = as.network.matrix(gnet,matrix.type ="adjacency",directed = FALSE)
g %v% 'name' = c(1:(Nt+2))

t = 97
g %v% "x" = c(data[t,c(1:20)*3-2],-0.5,0.2)
g %v% "y" = c(data[t,c(1:20)*3-1],-0.9,0.9)

p1 = ggnet2(g,vjust = -0.6,mode = c("x", "y"),color = "black",
            size = 1, edge.size = 1,edge.alpha = 1)+
  guides(color = 'none',size = 'none')+
  # geom_point( size = 12, color = "white") +
  geom_point( size = 5, color = "black",alpha = 0.7) +
  geom_point( size = 2.5,color = "white") +
  # geom_text(aes(label = c(1:(Nt+2))), color = "white",size =3,fontface = "bold") +
  theme(title=element_text(size=12),legend.text = element_text(size = 12))
p1

t = 117
g %v% "x" = c(data[t,c(1:20)*3-2],-0.5,0.2)
g %v% "y" = c(data[t,c(1:20)*3-1],-0.9,0.9)

p2 = ggnet2(g,vjust = -0.6,mode = c("x", "y"),color = "black",
            size = 1, edge.size = 1,edge.alpha = 1)+
  guides(color = 'none',size = 'none')+
  # geom_point( size = 12, color = "white") +
  geom_point( size = 5, color = "black",alpha = 0.7) +
  geom_point( size = 2.5,color = "white") +
  # geom_text(aes(label = c(1:(Nt+2))), color = "white",size =3,fontface = "bold") +
  theme(title=element_text(size=12),legend.text = element_text(size = 12))
p2

t = 137
g %v% "x" = c(data[t,c(1:20)*3-2],-0.5,0.2)
g %v% "y" = c(data[t,c(1:20)*3-1],-0.9,0.9)

p3 = ggnet2(g,vjust = -0.6,mode = c("x", "y"),color = "black",
            size = 1, edge.size = 1,edge.alpha = 1)+
  guides(color = 'none',size = 'none')+
  # geom_point( size = 12, color = "white") +
  geom_point( size = 5, color = "black",alpha = 0.7) +
  geom_point( size = 2.5,color = "white") +
  # geom_text(aes(label = c(1:(Nt+2))), color = "white",size =3,fontface = "bold") +
  theme(title=element_text(size=12),legend.text = element_text(size = 12))
p3

t = 147
g %v% "x" = c(data[t,c(1:20)*3-2],-0.5,0.2)
g %v% "y" = c(data[t,c(1:20)*3-1],-0.9,0.9)

p4 = ggnet2(g,vjust = -0.6,mode = c("x", "y"),color = "black",
            size = 1, edge.size = 1,edge.alpha = 1)+
  guides(color = 'none',size = 'none')+
  # geom_point( size = 12, color = "white") +
  geom_point( size = 5, color = "black",alpha = 0.7) +
  geom_point( size = 2.5,color = "white") +
  # geom_text(aes(label = c(1:(Nt+2))), color = "white",size =3,fontface = "bold") +
  theme(title=element_text(size=12),legend.text = element_text(size = 12))
p4

t = 157
g %v% "x" = c(data[t,c(1:20)*3-2],-0.5,0.2)
g %v% "y" = c(data[t,c(1:20)*3-1],-0.9,0.9)

p5 = ggnet2(g,vjust = -0.6,mode = c("x", "y"),
            size = 1, edge.size = 1,edge.alpha = 1)+
  guides(color = 'none',size = 'none')+
  # geom_point( size = 12, color = "white") +
  geom_point( size = 5, color = "black",alpha = 0.7) +
  geom_point( size = 2.5,color = "white") +
  # geom_text(aes(label = c(1:(Nt+2))), color = "white",size =3,fontface = "bold") +
  theme(title=element_text(size=12),legend.text = element_text(size = 12))
p5


library(cowplot)
gg <- ggdraw() + 
  draw_plot(p1, 0,0,0.18, 1) + draw_plot(p2, 0.2,0,0.18,1) +draw_plot(p3, 0.4,0,0.18,1) +
  draw_plot(p4, 0.6,0,0.18,1) +draw_plot(p5, 0.8,0,0.18,1) 
print(gg)
#########################################################################################例5
rm(list = ls())

library(jsonlite)
fish <- fromJSON("data/GoldenShinerfish/schooling_frames.json")
sq = seq(1,5000,3)       #采样 将30Hz的频率降成3Hz
data = lapply(sq, function(i){
  temp = fish[[i]]
  if(length(temp[[1]])>=1){
    temp
  }else{NA}
})
data=data[!is.na(data)]

N = length(data)
include = c(2:(N-2))    #排除边界点

Ntpool = sapply(c(1:N), function(i){
  length(data[[i]][[1]])
})

1-sum(Ntpool)/(300*length(Ntpool))

###归一化
x1 = unlist(sapply(1:N, function(i){
  data[[i]][[1]]
}))
x2 = unlist(sapply(1:N, function(i){
  data[[i]][[2]]
}))
v1 = unlist(sapply(1:N, function(i){
  data[[i]][[3]]
}))
v2 = unlist(sapply(1:N, function(i){
  data[[i]][[4]]
}))
del = union(which(is.na(v1)),which(is.na(v2)))  #去掉速度中的NA值
x = rbind(x1[-del],x2[-del])
v = rbind(v1[-del],v2[-del])
xv0 = rbind(x,v)   #不归一化 用来计算极化和旋转
x = (x-mean(x))/sqrt(sum((x-mean(x))^2)/(length(x)-1))
# v = (v-mean(v))/sqrt(sum((v-mean(v))^2)/(length(v)-1))
Ntsum = vector()
Ntsum[1] = 0
for (i in 1:N) {
  Ntsum[i+1] = Ntsum[i]+Ntpool[i]
}
for (i in 1:length(del)) {
  temp = which(Ntsum>=del[i])
  Ntsum[temp] = Ntsum[temp]-1
}
data = lapply(1:N, function(i){
  a = x[,(Ntsum[i]+1):Ntsum[i+1]]
  b = v[,(Ntsum[i]+1):Ntsum[i+1]]
  rbind(a,b)
})
data0 = lapply(1:N, function(i){
  xv0[,(Ntsum[i]+1):Ntsum[i+1]]
})

####计算断点
Nt = floor(mean(Ntpool))
distmat =  lapply(c(1:N), function(i){
  piece = data[[i]]
  n = dim(piece)[2]
  if(n < Nt){
    if(i>1){
      addpool = data[[i-1]]
    }else{
      addpool = data[[i+1]]
    }
    nn = dim(addpool)[2]
    ind = sample(nn,Nt-n,replace = TRUE)
    piece = cbind(piece,addpool[,ind])
  }else if(n > Nt){
    ind = sample(n,n-Nt)
    piece = piece[,-ind]
  }
  as.matrix(dist(t(piece)))
})

library(Rfast)
em = mean(sapply(c(1:N),function(i){
  temp = distmat[[i]]
  diag(temp)= max(temp)
  max(colMins(temp,value = TRUE))
}))
e = seq(0,em,em/(Nt/2))

library(TDAstats)
bettimat = sapply(c(1:N), function(i){
  phom <- calculate_homology(distmat[[i]],dim = 0,format = "distmat")
  betti = sapply(e, function(e){
    length(which(phom[,3]>e))
  })
  betti
})

# library(corrplot)
# library(RColorBrewer)
# corrplot(a, is.corr = FALSE, col.lim = c(0, 1), method = "color",tl.pos = 'n',
#          cl.pos = 'n',tl.col = 'black',col = colorRampPalette(brewer.pal(9, "Greys"))(200),
#          bg = 'white')
# corrplot(a, is.corr = FALSE, col.lim = c(0, 1), method = "color",tl.pos = 'n',
#          cl.pos = 'n',tl.col = 'black',col = colorRampPalette(brewer.pal(9, "Greys"))(200),
#          bg = 'white')

betti = bettimat    #[1:40,]

library(ecp)
y <- e.agglo(t(betti))     #betti有效果
poolB = y$estimates
poolB = poolB[which(poolB%in%include)]

# y <- e.divisive(t(betti))     #betti有效果
# poolB = y$estimates
# poolB = poolB[which(poolB%in%include)]

life = lapply(c(1:N), function(i){
  phom <- calculate_homology(distmat[[i]],dim = 0,format = 'distmat')[,3]
  phom
})
l = sapply(1:N, function(j){length(life[[j]])})
se = min(l)
mid = floor(se/10*3)
selife = sapply(1:N, function(j){
  temp = length(life[[j]])
  life[[j]][(floor(temp/2)-mid+1):(floor(temp/2)+mid)]
})

plot(selife[1,],type = "l",ylim= c(0,4))
for (i in 2:dim(selife)[1]) {
  lines(selife[i,])
}

y <- e.agglo(t(selife))    
poolL = y$estimates
poolL = poolL[which(poolL%in%include)]

######################################################################################################
# https://github.com/mgsosna/Rotation_v_Polarization
# Calculate the global rotation and polarization, and label the group state. Based on Couzin et al.
# 2002 (Journal of Theoretical Biology) and Tunstrom et al. 2013 (PLoS Computational Biology).
# Matt Grobis | Jan 2018
#
# Global rotation, polarization functions:
# - Inputs: matrices, where each row is an individual's time series for x- and y-coordinates and
#           the unit vectors of their orientation.
# - Outputs: vectors of group rotation or polarization  
#   
# Group state function:
# - Input: global rotation and polarization vectors
# - Output: vector of group labels
#
#--------------------------------------------------------------------------------
# calculate_polarization
# - Definition: how aligned is the group? 0 = no alignment (everyone facing diff 
#   directions), 1 = perfectly aligned (everyone facing the same direction)
#
# - Inputs: the x- and y-unit matrices of the individuals' orientations. 
#    o Rows = individuals
#    o Columns = time points
#
#-----------------------------------------------------------------------------------------
# calculate_rotation
# - Definition: the group's degree of rotation around its center of mass. 
#   0 = no rotation, 1 = everyone swimming in a big torus
#
# - Inputs: the x- and y-position matrices, and the x- and y-unit vectors of orientation.
#    o Rows = individuals
#    o Columns = time points
#
#-----------------------------------------------------------------------------------------
# identify_state
# - Definition: assign the group state
#    o Rotation < 0.35, Polarization < 0.35: swarm
#    o Rotation < 0.35, Polarization > 0.65: polarized
#    o Rotation > 0.65, Polarization < 0.35: milling
#    o All others: transition
#
######################################################################################################

calculate_polarization <- function(BO.x, BO.y){   
  # Step 1: Get mean heading in x direction, mean in y direction
  # - You want to take all the negative (left/down) components, all the right (right/up)
  # components, let them accumulate and cancel each other out, and divide what remains by N
  head.x <- colMeans(BO.x, na.rm = T)
  head.x[head.x == "NaN"] <- NA
  
  head.y <- colMeans(BO.y, na.rm = T)
  head.y[head.y == "NaN"] <- NA
  
  # Step 2: find the magnitude of the resulting vector. 
  # - This magnitude will be between zero and 1. 0 = all the individuals' headings canceled 
  #   each other out. 1 = all the individuals' headings were in the same exact direction.
  
  return(sqrt(head.x^2 + head.y^2))
  
}

#################################################################################################

calculate_rotation <- function(xs, ys, BO.x, BO.y){
  # 1. Calculate the centroid of the group
  centX <- colMeans(xs, na.rm = T)
  centX[centX == "NaN"] <- NA
  
  centY <- colMeans(ys, na.rm = T)
  centY[centY == "NaN"] <- NA
  
  #--------------------------------------------------------------------------------
  # 2. Find the distance to centroid
  dis.x <- sweep(xs, 2, centX, "-")  # Position relative to centroid in x domain
  dis.y <- sweep(ys, 2, centY, "-")  # Position relative to centroid in y domain
  
  tot.dist <- sqrt(dis.x^2 + dis.y^2)
  
  #--------------------------------------------------------------------------------
  # 3. Find relative unit vector from centroid towards each fish
  rel.u.x <- dis.x / tot.dist
  rel.u.y <- dis.y / tot.dist
  
  #--------------------------------------------------------------------------------------
  # Group rotation
  rot <- BO.x*rel.u.y - BO.y*rel.u.x
  
  return(abs(colMeans(rot, na.rm = T)))
  
}

######################################################################################################
# Classify the group state based on the global rotation and polarization.
# - 1 Swarm: rotation < 0.35 & polarization < 0.35
# - 2 Milling: rotation > 0.65 & polarization < 0.35
# - 3 Polarized: rotation < 0.35 & polarization > 0.65
# - 4 Transition: all other combinations

identify_state <- function(rot, pol){
  
  ifelse(rot < 0.35 & pol < 0.35, 'Swarm', 
         
         ifelse(rot < 0.35 & pol > 0.65, 'Polarized',
                
                ifelse(rot > 0.65 & pol < 0.35, 'Milling', 'Transition')))
}

poro = sapply(c(1:N),function(i){
  piece = data0[[i]]
  xs = piece[1,]
  ys = piece[2,]
  BO.x = piece[3,]
  BO.y = piece[4,]
  
  po = sqrt(mean(BO.x)^2+mean(BO.y)^2)     #polarization
  
  dis.x = mean(xs)-xs
  dis.y = mean(ys)-ys
  tot.dist <- sqrt(dis.x^2 + dis.y^2)
  rel.u.x <- dis.x / tot.dist
  rel.u.y <- dis.y / tot.dist
  rot <- BO.x*rel.u.y - BO.y*rel.u.x        
  ro = abs(mean(rot))    # Group rotation
  
  state = identify_state(ro, po)      #state
  
  c(po,ro,state)
})

# write.table(poolB,file = 'poolB.txt', sep = ' ', row.names = FALSE, col.names = FALSE)
# write.table(poolL,file = 'poolL.txt', sep = ' ', row.names = FALSE, col.names = FALSE)
# write.table(poro,file = 'poro.txt', sep = ' ', row.names = FALSE, col.names = FALSE)
poolB = as.vector(as.matrix(read.table('poolB.txt',sep = ' ', header = FALSE)))
poolL = as.vector(as.matrix(read.table('poolL.txt',sep = ' ', header = FALSE)))
poro = read.table('poro.txt',sep = ' ', header = FALSE)


po = round(as.numeric(poro[1,]),2)
ro = round(as.numeric(poro[2,]),2)
state = poro[3,]


######################################去掉标签
library(RColorBrewer)
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
color1 <- colorRampPalette(rev(brewer.pal(11,'Spectral')))(28)
color2 = c("#dc745c","#ffe5ce","#FFE6CF","#b0bcd6","#6a84b5")
c1 =  brewer.pal(9,"Greys")[5]   #'#DC302E' 
c2 = brewer.pal(9,"Greys")[7]      #'#F08080'
c3 = 'black'

ggo = data.frame(Timepoint = rep(c(1:length(po)),2),
                 value = c(po,ro),
                 Parameter = rep(c('P','R'),each = length(po)))
ggcp1 = data.frame(x = poolB)

library(ggplot2)
p1 = ggplot(data = ggo,mapping = aes(x=Timepoint,y=value,group=Parameter,color=Parameter))+
  geom_line(size= 0.3)+labs(title='(a) CPD based on the Betti Number')+
  theme_set(theme_bw())+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
        axis.title.y = element_blank(), axis.text.y = element_text(size = 12),
        plot.title=element_text(size=12),legend.text = element_text(size = 12),legend.position = 'none',
        panel.grid=element_blank())+
  geom_vline(data = ggcp1, aes(xintercept = x), size = 1,color = c3)+
  scale_color_manual(values=c(c1,c2))+
  scale_y_continuous(limits = c(0,1), breaks = c(0, 0.5, 1))
  
p1

df <- data.frame(
  x = c(1:length(state)),
  z = factor(state),
  w = rep(1,length(state))
)
df$z = factor(df$z, levels=c('Swarm','Milling','Polarized','Transition'))
ggdf = data.frame(x = rep(c(1:length(po)),2),
                 y = c(po,ro),
                 Parameter = rep(c('P','R'),each = length(po)))
p2= ggplot() +
  geom_rect(data = df, mapping = aes(xmin = x - w / 2, xmax = x + w / 2, ymin = 0, ymax = 1,fill = z),
            alpha = 0.5)+labs(title='(c) Reference Change Points')+
              scale_fill_manual(values=color2[c(1,3,4,5)])+guides(fill=guide_legend(title=NULL))+
  theme_set(theme_bw())+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
        axis.title.y = element_blank(), axis.text.y = element_text(size = 12),
        plot.title=element_text(size=12),legend.text = element_text(size = 12),legend.position = 'bottom',
        panel.grid=element_blank())+
  geom_line(data = ggdf,mapping = aes(x=x,y=y,group=Parameter,color=Parameter),size = 0.3)+
  guides(color='none',fill='none')+
  scale_color_manual(values=c(c1,c2))+
  scale_y_continuous(limits = c(0,1), breaks = c(0, 0.5, 1))
  
p2

ggcp2 = data.frame(x = poolL)
p3 = ggplot(data = ggo,mapping = aes(x=Timepoint,y=value,group=Parameter,color=Parameter))+
  geom_line(size = 0.3)+labs(title='(b) CPD based on the Persistence')+
  theme_set(theme_bw())+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
        axis.title.y = element_blank(), axis.text.y = element_text(size = 12),
        plot.title=element_text(size=12),legend.text = element_text(size = 12),legend.position = 'bottom',
        panel.grid=element_blank())+
  geom_vline(data = ggcp2, aes(xintercept = x), size = 1,color = c3)+
  guides(colour='none')+
  scale_color_manual(values=c(c1,c2))+
  scale_y_continuous(limits = c(0,1), breaks = c(0, 0.5, 1))
  
p3

library(cowplot)
gg1 <- ggdraw() + 
  draw_plot(p1, 0,0.7,1, 0.3) + draw_plot(p3, 0,0.4,1,0.3) + draw_plot(p2, 0,0.1,1, 0.3) 
  # draw_plot_label(c("a", "b", "c"), c(0, 0,0), c(1,0.72,0.38), size = 15,colour = "black")
print(gg1)


#################################################只取标签
ggo = data.frame(Timepoint = rep(c(1:length(po)),2),
                 value = c(po,ro),
                 Parameter = rep(c('Op','Or'),each = length(po)))
ggcp1 = data.frame(x = poolB)

library(ggplot2)
p1 = ggplot(data = ggo,mapping = aes(x=Timepoint,y=value,group=Parameter,color=Parameter))+
  geom_line()+labs(title='Betti Number')+
  theme_classic()+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
        axis.title.y = element_blank(), axis.text.y = element_text(size = 12),
        plot.title=element_text(size=12),legend.text = element_text(size = 12),legend.position = 'none')+
  geom_vline(data = ggcp1, aes(xintercept = x), size = 1)+
  scale_color_manual(values=c(c1,c2))+
  scale_y_continuous(limits = c(0,1), breaks = c(0, 0.5, 1))

p1

df = data.frame(
  x = c(1:length(state)),
  z = factor(state),
  w = rep(1,length(state))
)
df$z = factor(df$z, levels=c('Swarm','Milling','Polarized','Transition'))
ggdf = data.frame(x = rep(c(1:length(po)),2),
                  y = c(po,ro),
                  Parameter = rep(c('Op','Or'),each = length(po)))
p2= ggplot() +
  geom_rect(data = df, mapping = aes(xmin = x - w / 2, xmax = x + w / 2, ymin = 0, ymax = 1,fill = z),
            alpha = 0.5)+labs(title='Pattern Switching based on Order Parameters')+
  scale_fill_manual(values=color2[c(1,3,4,5)])+guides(fill=guide_legend(title=NULL))+
  theme_classic()+
  theme(axis.title.x = element_text(size=12), axis.text.x = element_text(size=12),
        axis.title.y = element_blank(), axis.text.y = element_text(size = 12),
        plot.title=element_text(size=12),legend.title = element_blank(),
        legend.text = element_text(size = 12),legend.position = 'bottom')+
  geom_line(data = ggdf,mapping = aes(x=x,y=y,group=Parameter,color=Parameter),size = 1.5)+
  labs(x = expression(italic(t)))+
  scale_color_manual(values=c(c1,c2),labels = c(expression(italic(O[p]),italic(O[r]))))+
  scale_y_continuous(limits = c(0,1), breaks = c(0, 0.5, 1))

p2


library(cowplot)
gg2 <- ggdraw() + 
  draw_plot(p1, 0,0.5,1, 0.5) + draw_plot(p2, 0,0,1,0.5)
# draw_plot_label(c("a", "b", "c"), c(0, 0,0), c(1,0.72,0.38), size = 15,colour = "black")
print(gg2)

library(cowplot)
gg <- ggdraw() + 
  draw_plot(gg1, 0,0.5,1, 0.5) + draw_plot(gg2, 0,0,1,0.5)
# draw_plot_label(c("a", "b", "c"), c(0, 0,0), c(1,0.72,0.38), size = 15,colour = "black")
print(gg)


# ##############对比画图
# data = lapply(1:N, function(i){
#   a = x[1,(Ntsum[i]+1):Ntsum[i+1]]
#   b = x[2,(Ntsum[i]+1):Ntsum[i+1]]
#   c = v[1,(Ntsum[i]+1):Ntsum[i+1]]
#   d = v[2,(Ntsum[i]+1):Ntsum[i+1]]
#   list(a,b,c,d)
# })
# distmat =  lapply(c(1:N), function(i){
#   piece = data[[i]]
#   p = rbind(unlist(piece[[1]]),unlist(piece[[2]]))
#   v = rbind(unlist(piece[[3]]),unlist(piece[[4]]))
#   n = dim(p)[2]
#   if(n < Nt){
#     if(i>1){
#       addpool = data[[i-1]]
#     }else{
#       addpool = data[[i+1]]
#     }
#     nn = length(addpool[[1]])
#     ind = sample(nn,Nt-n,replace = TRUE)
#     p = cbind(p,rbind(addpool[[1]][ind],addpool[[2]][ind]))
#     v = cbind(v,rbind(addpool[[3]][ind],addpool[[4]][ind]))
#   }else if(n > Nt){
#     ind = sample(n,n-Nt)
#     p = p[,-ind]
#     v = v[,-ind]
#   }
#   piece = rbind(p,v)
#   as.matrix(dist(t(piece)))
# })
# 
# library(Rfast)
# em = max(sapply(c(1:N),function(i){
#   temp = distmat[[i]]
#   diag(temp)= max(temp)
#   max(colMins(temp,value = TRUE))
# }))
# e = seq(0,em,em/200)
# 
# library(TDAstats)
# bettimat = sapply(c(1:N), function(i){
#   phom <- calculate_homology(distmat[[i]],dim = 0,format = "distmat")
#   betti = sapply(e, function(e){
#     length(which(phom[,3]>e))
#   })
#   betti
# })
# # bettimat = bettimat[-which(rowSums(bettimat)==0),] #一般没有
# 
# a = bettimat[1:40,1:100]
# a = a/max(a)
# library(corrplot)
# library(RColorBrewer)
# corrplot(a, is.corr = FALSE, col.lim = c(0, 1), method = "color",addgrid.col = 'black',tl.pos = 'n',
#          cl.pos = 'n',tl.col = 'black',col = colorRampPalette(brewer.pal(9, "Greys"))(200),
#          bg = 'white',title = 'M=200, Adjucted=TRUE, No.factor=40')
# 
# ind = c(1:dim(a)[1])
# ind = sort(ind,decreasing=T) 
# temp = a[ind,]
# corrplot(temp, is.corr = FALSE, col.lim = c(0, 1), method = "color",addgrid.col = 'black',tl.pos = 'n',
#          cl.pos = 'n',tl.col = 'black',col = colorRampPalette(brewer.pal(9, "Greys"))(200),
#          bg = 'white',title = 'M=200, Adjucted=TRUE, No.factor=40')
# 
# ############
# e = seq(0,em,em/100)
# 
# library(TDAstats)
# bettimat = sapply(c(1:N), function(i){
#   phom <- calculate_homology(distmat[[i]],dim = 0,format = "distmat")
#   betti = sapply(e, function(e){
#     length(which(phom[,3]>e))
#   })
#   betti
# })
# # bettimat = bettimat[-which(rowSums(bettimat)==0),] #一般没有
# 
# b = bettimat[1:40,1:100]
# b = b/max(b)
# library(corrplot)
# corrplot(b, is.corr = FALSE, col.lim = c(0, 1), method = "color",addgrid.col = 'black',tl.pos = 'n',
#          cl.pos = 'n',tl.col = 'black',col = colorRampPalette(brewer.pal(9, "Greys"))(200),
#          bg = 'white',title = 'M=100, Adjucted=TRUE, No.factor=40')
# 
# ind = c(1:dim(b)[1])
# ind = sort(ind,decreasing=T) 
# temp = b[ind,]
# corrplot(temp, is.corr = FALSE, col.lim = c(0, 1), method = "color",addgrid.col = 'black',tl.pos = 'n',
#          cl.pos = 'n',tl.col = 'black',col = colorRampPalette(brewer.pal(9, "Greys"))(200),
#          bg = 'white',title = 'M=100, Adjucted=TRUE, No.factor=40')
# 
# ############
# e = seq(0,em,em/40)
# 
# library(TDAstats)
# bettimat = sapply(c(1:N), function(i){
#   phom <- calculate_homology(distmat[[i]],dim = 0,format = "distmat")
#   betti = sapply(e, function(e){
#     length(which(phom[,3]>e))
#   })
#   betti
# })
# # bettimat = bettimat[-which(rowSums(bettimat)==0),] #一般没有
# 
# c = bettimat[,1:100]
# c = c/max(c)
# library(corrplot)
# corrplot(c, is.corr = FALSE, col.lim = c(0, 1), method = "color",addgrid.col = 'black',tl.pos = 'n',
#          cl.pos = 'n',tl.col = 'black',col = colorRampPalette(brewer.pal(9, "Greys"))(200),
#          bg = 'white',title = 'M=40, Adjucted=TRUE, No.factor=40')
# 
# ind = c(1:dim(c)[1])
# ind = sort(ind,decreasing=T) 
# temp = c[ind,]
# corrplot(temp, is.corr = FALSE, col.lim = c(0, 1), method = "color",addgrid.col = 'black',tl.pos = 'n',
#          cl.pos = 'n',tl.col = 'black',col = colorRampPalette(brewer.pal(9, "Greys"))(200),
#          bg = 'white',title = 'M=40, Adjucted=TRUE, No.factor=40')
# 
# ################################
# distmat =  lapply(c(1:N), function(i){
#   piece = data[[i]]
#   p = rbind(unlist(piece[[1]]),unlist(piece[[2]]))
#   v = rbind(unlist(piece[[3]]),unlist(piece[[4]]))
#   p = (p-mean(p))/sqrt(sum((p-mean(p))^2)/(length(p)-1))
#   v = (v-mean(v))/sqrt(sum((v-mean(v))^2)/(length(v)-1))
#   piece = rbind(p,v)
#   as.matrix(dist(t(piece)))
# })
# 
# library(Rfast)
# em = max(sapply(c(1:N),function(i){
#   temp = distmat[[i]]
#   diag(temp)= max(temp)
#   max(colMins(temp,value = TRUE))
# }))
# e = seq(0,em,em/100)
# 
# library(TDAstats)
# bettimat = sapply(c(1:N), function(i){
#   phom <- calculate_homology(distmat[[i]],dim = 0,format = "distmat")
#   betti = sapply(e, function(e){
#     length(which(phom[,3]>e))
#   })
#   betti
# })
# # bettimat = bettimat[-which(rowSums(bettimat)==0),]
# 
# d = bettimat[1:40,1:100]
# d = d/max(d)
# library(corrplot)
# corrplot(d, is.corr = FALSE, col.lim = c(0, 1), method = "color",addgrid.col = 'black',tl.pos = 'n',
#          cl.pos = 'n',tl.col = 'black',col = colorRampPalette(brewer.pal(9, "Greys"))(200),
#          bg = 'white',title = 'M=100, Adjucted=FALSE, No.factor=40')
# 
# ind = c(1:dim(d)[1])
# ind = sort(ind,decreasing=T) 
# temp = d[ind,]
# corrplot(temp, is.corr = FALSE, col.lim = c(0, 1), method = "color",addgrid.col = 'black',tl.pos = 'n',
#          cl.pos = 'n',tl.col = 'black',col = colorRampPalette(brewer.pal(9, "Greys"))(200),
#          bg = 'white',title = 'M=100, Adjucted=FALSE, No.factor=40')
# 
# corrplot(d, is.corr = FALSE, col.lim = c(0, 1), method = "color",addgrid.col = 'black',tl.pos = 'n',
#          tl.col = 'black',col = colorRampPalette(brewer.pal(9, "Greys"))(200),
#          bg = 'white',title = 'M=100, Adjucted=FALSE, No.factor=40')


##################################################无效果数据，可以探索原因############################################################
# rm(list = ls())
# library(TDAstats)
# library(ecp)
# source('values.R')
# 
# ############################################例1
# #Retrospective multivariate Bayesian change-point analysis: 
# #A simultaneous single change in the mean of several hydrological sequences
# library(bcp)
# data("QuebecRivers")
# bcpr.rivers <- bcp(QuebecRivers)
# plot(bcpr.rivers, main="Quebec River Streamflow Change Point Analysis",
#      xlab="Year", xaxlab = 1972:1994)
# 
# QuebecRivers = (QuebecRivers-mean(QuebecRivers))/
#   sqrt(sum((QuebecRivers-mean(QuebecRivers))^2)/(length(QuebecRivers)-1))
# 
# N = dim(QuebecRivers)[1]
# Nt = dim(QuebecRivers)[2]
# distmat =  lapply(c(1:N), function(i){
#   piece = rbind(QuebecRivers[i,],rep(0,Nt))
#   as.matrix(dist(t(piece)))
# })
# 
# library(Rfast)
# em = max(sapply(c(1:N),function(i){
#   temp = distmat[[i]]
#   diag(temp)= max(temp)
#   max(colMins(temp,value = TRUE))
# }))
# e = seq(0,em,em/100)
# 
# 
# bettimat = sapply(c(1:N), function(i){
#   phom <- calculate_homology(distmat[[i]],dim = 0,format = "distmat")
#   betti = sapply(e, function(e){
#     length(which(phom[,3]>e))
#   })
#   betti
# })
# 
# library(corrplot)
# corrplot(bettimat/max(bettimat))
# y <- e.divisive(t(bettimat))     #betti无效果
# pool = y$estimates
# 
# life = sapply(c(1:N), function(i){
#   phom <- calculate_homology(distmat[[i]],dim = 0,format = "distmat")
#   phom[,3]
# })
# plot(life[1,],type = "l",ylim= c(0,6))
# for (i in 2:(Nt-1)) {
#   lines(life[i,])
# }
# 
# lim = Nt-1
# selife = life[1:lim,]
# y <- e.agglo(t(selife))     #life无效果
# pool = y$estimates
# 
# y2 <- e.divisive(QuebecRivers)     #原方法一致
# pool2 = y2$estimates
# 
# y2 <- e.agglo(QuebecRivers)     #原方法一致
# pool2 = y2$estimates
# 
# #######################################例2
# rm(list = ls())
# source('values.R')
# library(TDAstats)
# library(kcpRS)
# library(ecp)
# data(MentalLoad)
# cover = vector()
# clus = vector()
# 
# x = MentalLoad[,2:4]
# N = dim(x)[1]
# include = c(2:(N-2))    #排除边界点
# Nt = dim(x)[2]
# CPreal = sapply(c(1:max(MentalLoad[,1])),function(i){
#   length(which(MentalLoad[,1]==i))
# })
# for (i in 2:4) {
#   CPreal[i] = CPreal[i]+CPreal[i-1]
# }
# 
# kcp_mean <- kcpRS(data = MentalLoad[, 2:4], RS_fun = runMean, 
#                   RS_name = 'Mean', wsize = 20, nperm = 1000, Kmax = 10, 
#                   alpha = (0.05 / 4), varTest = FALSE, ncpu = detectCores())  #有一点效果
# pool = unlist(kcp_mean$CPs_given_K[kcp_mean$BestK+1,3:(kcp_mean$BestK+3)])
# pool = pool[-length(pool)]
# CPfind = pool[which(pool%in%include)]
# cover[1] = covering(N,CPreal,CPfind)
# clus[1] = Fvalue(CPreal,CPfind)
# 
# kcp_var <- kcpRS(data = MentalLoad[, 2:4],RS_fun = runVar,
#                  RS_name = 'Variance',  wsize = 20, nperm = 1000, Kmax = 10, 
#                  alpha = (0.05 / 4), varTest = FALSE, ncpu = detectCores())
# pool = unlist(kcp_var$CPs_given_K[kcp_var$BestK+1,3:(kcp_var$BestK+3)])
# pool = pool[-length(pool)]
# CPfind = pool[which(pool%in%include)]
# cover[2] = covering(N,CPreal,CPfind)
# clus[2] = Fvalue(CPreal,CPfind)
# 
# kcp_corr <- kcpRS(data = MentalLoad[, 2:4], RS_fun = runCorr,
#                   RS_name = 'Correlation', wsize = 20, nperm = 1000, Kmax = 10, 
#                   alpha = (0.05 / 4), varTest = FALSE, ncpu = detectCores())
# pool = unlist(kcp_corr$CPs_given_K[kcp_corr$BestK+1,3:(kcp_corr$BestK+3)])
# pool = pool[-length(pool)]
# CPfind = pool[which(pool%in%include)]
# cover[3] = covering(N,CPreal,CPfind)
# clus[3] = Fvalue(CPreal,CPfind)
# 
# kcp_AR <- kcpRS(data = MentalLoad[, 2:4], RS_fun = runAR, 
#                 RS_name = 'Autocorrelation', wsize = 20, nperm = 1000, Kmax = 10, 
#                 alpha = (0.05 / 4), varTest = FALSE, ncpu = detectCores())
# pool = unlist(kcp_AR$CPs_given_K[kcp_AR$BestK+1,3:(kcp_AR$BestK+3)])
# pool = pool[-length(pool)]
# CPfind = pool[which(pool%in%include)]
# cover[4] = covering(N,CPreal,CPfind)
# clus[4] = Fvalue(CPreal,CPfind)
# 
# 
# distmat =  lapply(c(1:N), function(i){
#   piece = rbind(x[i,],rep(0,Nt))
#   as.matrix(dist(t(piece)))
# })
# 
# library(Rfast)
# em = max(sapply(c(1:N),function(i){
#   temp = distmat[[i]]
#   diag(temp)= max(temp)
#   max(colMins(temp,value = TRUE))
# }))
# e = seq(0,em,em/80)
# 
# 
# bettimat = sapply(c(1:N), function(i){
#   phom <- calculate_homology(distmat[[i]],dim = 0,format = "distmat")
#   betti = sapply(e, function(e){
#     length(which(phom[,3]>e))
#   })
#   betti
# })
# 
# library(corrplot)
# corrplot(bettimat[,1:100]/max(bettimat[,1:100]))
# y <- e.divisive(t(bettimat))     #betti无效果
# pool = y$estimates
# CPfind = pool[which(pool%in%include)]
# cover[5] = covering(N,CPreal,CPfind)
# clus[5] = Fvalue(CPreal,CPfind)
# 
# life = sapply(c(1:N), function(i){
#   phom <- calculate_homology(distmat[[i]],dim = 0,format = "distmat")
#   phom[,3]
# })
# plot(life[1,],type = "l",ylim= c(0,106))
# for (i in 2:(Nt-1)) {
#   lines(life[i,])
# }
# 
# lim = Nt-1
# selife = life[1:lim,]
# y <- e.agglo(t(selife))     #life有一点效果 找的点太多
# pool = y$estimates
# CPfind = pool[which(pool%in%include)]
# cover[6] = covering(N,CPreal,CPfind)
# clus[6] = Fvalue(CPreal,CPfind)
# 
# y2 <- e.divisive(x)     #原方法更好一点点
# pool = y2$estimates
# CPfind = pool[which(pool%in%include)]
# cover[7] = covering(N,CPreal,CPfind)
# clus[7] = Fvalue(CPreal,CPfind)
# 
# y <- e.agglo(x)     #life有一点效果 找的点太多
# pool = y$estimates
# CPfind = pool[which(pool%in%include)]
# cover[8] = covering(N,CPreal,CPfind)
# clus[8] = Fvalue(CPreal,CPfind)
# 
# cover
# clus
########################################################其他数据勉强有效果###########################################
############################################例3########################################################
# source('values.R')
# library(kcpRS)
# data(CO2Inhalation)
# x = CO2Inhalation[,2:10]
# N = dim(x)[1]
# Nt = dim(x)[2]
# include = c(2:(N-2))    #排除边界点
# CPreal = sapply(c(1:max(CO2Inhalation[,1])),function(i){
#   length(which(CO2Inhalation[,1]==i))
# })
# for (i in 2:3) {
#   CPreal[i] = CPreal[i]+CPreal[i-1]
# }
# 
# cover = vector()
# clus = vector()
# 
# kcp_mean <- kcpRS(data = CO2Inhalation[, 2:10], RS_fun = runMean,
#                   RS_name = 'Mean', wsize = 25, nperm = 1000, Kmax = 10,
#                   alpha = (0.05 / 3), varTest = FALSE, ncpu = detectCores())
# pool = unlist(kcp_mean$CPs_given_K[kcp_mean$BestK+1,3:(kcp_mean$BestK+3)])
# pool = pool[-length(pool)]
# CPfind = pool[which(pool%in%include)]
# cover[1] = covering(N,CPreal,CPfind)
# clus[1] = Fvalue(CPreal,CPfind)
# 
# kcp_var <- kcpRS(data = CO2Inhalation[, 2:10], RS_fun = runVar,
#                  RS_name = 'Variance', wsize = 25, nperm = 1000, Kmax = 10,
#                  alpha = (0.05 / 3), varTest = FALSE, ncpu = detectCores())
# pool = unlist(kcp_mean$CPs_given_K[kcp_mean$BestK+1,3:(kcp_mean$BestK+3)])
# pool = pool[-length(pool)]
# CPfind = pool[which(pool%in%include)]
# cover[2] = covering(N,CPreal,CPfind)
# clus[2] = Fvalue(CPreal,CPfind)
# 
# kcp_corr <- kcpRS(data = CO2Inhalation[, 2:10], RS_fun = runCorr, 
#                   RS_name = 'Correlation', wsize = 20, nperm = 1000, Kmax = 10,
#                   alpha = (0.05 / 3), varTest = FALSE, ncpu = detectCores())
# pool = unlist(kcp_mean$CPs_given_K[kcp_mean$BestK+1,3:(kcp_mean$BestK+3)])
# pool = pool[-length(pool)]
# CPfind = pool[which(pool%in%include)]
# cover[3] = covering(N,CPreal,CPfind)
# clus[3] = Fvalue(CPreal,CPfind)
# 
# distmat =  lapply(c(1:N), function(i){
#   piece = rbind(x[i,],rep(0,Nt))
#   as.matrix(dist(t(piece)))
# })
# 
# library(Rfast)
# em = max(sapply(c(1:N),function(i){
#   temp = distmat[[i]]
#   diag(temp)= max(temp)
#   max(colMins(temp,value = TRUE))
# }))
# e = seq(0,em,em/100)
# 
# 
# bettimat = sapply(c(1:N), function(i){
#   phom <- calculate_homology(distmat[[i]],dim = 0,format = "distmat")
#   betti = sapply(e, function(e){
#     length(which(phom[,3]>e))
#   })
#   betti
# })
# 
# library(corrplot)
# corrplot(bettimat/max(bettimat))
# y <- e.divisive(t(bettimat))     #betti有效果！
# pool = y$estimates
# pool = pool[-length(pool)]
# CPfind = pool[which(pool%in%include)]
# cover[4] = covering(N,CPreal,CPfind)
# clus[4] = Fvalue(CPreal,CPfind)
# 
# y <- e.divisive(x)     #betti一般般
# pool = y$estimates
# pool = pool[-length(pool)]
# CPfind = pool[which(pool%in%include)]
# cover[5] = covering(N,CPreal,CPfind)
# clus[5] = Fvalue(CPreal,CPfind)
# 
# life = sapply(c(1:N), function(i){
#   phom <- calculate_homology(distmat[[i]],dim = 0,format = "distmat")
#   phom[,3]
# })
# plot(life[1,],type = "l",ylim= c(0,106))
# for (i in 2:(Nt-1)) {
#   lines(life[i,])
# }
# 
# lim = Nt-1
# selife = life[1:lim,]
# y <- e.agglo(t(selife))     #life无效果 找的点太多
# pool = y$estimates
# pool = pool[-length(pool)]
# CPfind = pool[which(pool%in%include)]
# cover[6] = covering(N,CPreal,CPfind)
# clus[6] = Fvalue(CPreal,CPfind)
# 
# y <- e.agglo(x)     #life无效果 找的点太多
# pool = y$estimates
# pool = pool[-length(pool)]
# CPfind = pool[which(pool%in%include)]
# cover[7] = covering(N,CPreal,CPfind)
# clus[7] = Fvalue(CPreal,CPfind)
# cover
# clus








