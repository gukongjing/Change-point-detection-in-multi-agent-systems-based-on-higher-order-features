####所有得矩阵都是完整可用得，尽量选择偶数代表矩阵（对称矩阵），方便后续辨识处理
library(igraph)
library(network)

netconstruction  = function(Nt,kind,disturbnum,k,m){
  ###构建10种不同的邻接矩阵
  
  ###输入
  #【1】Nt: 集群规模,integer，范围5-500.
  #【2】kind：指定矩阵类型，integer，范围1-10.其中1代表BCJ不对称层级矩阵,2代表CJ对称层级矩阵,3代表BLDJ不对称类对角,4代表LDJ对称类对角,
  #5代表BDC不对称矩阵,6代表DC对称矩阵,7代表BFK分块不对称矩阵,8代表FK分块对称矩阵,9代表SSJ对称上三角矩阵,10代表FSSJ分块上三角矩阵
  #【3】disturbnum: 矩阵随机连边数，integer，范围0-Nt，矩阵5 6 7 8 9 10都需要该参数。默认值设为0.1*Nt
  #【4】k：分块数目，integer，范围3-10，矩阵7 8 10需要该参数。默认值设为5
  #【5】m：分支数，integer，范围3-8，矩阵1 2需要该参数。默认值为5
  
  ###输出
  #【1】A：邻接矩阵，Nt*Nt，元素为0或1
  
  #默认设置
  # Nt = 150     
  # disturbnum = 0.1*Nt    
  # k = 5         
  # m = 5         
  
  #其他需要计算的参数
  node = c(1:Nt)     #点集
  gnum <- Nt/k   #每块数目
  n = round(log((Nt*(m-1)+1),m))  #计算层级结构矩阵的层级数
  
  ##1 层级结构
  scale = (1-m^n)/(1-m)
  A = matrix(0, nrow = max(scale,Nt), ncol = max(scale,Nt))
  for (i in 1:((1-m^(n-1))/(1-m))){
    A[i,((i-1)*m+2):(i*m+1)] = 1
  }
  if(Nt>scale){
    for (i in (scale+1):Nt) {
      temp = sample(node[-i],1)
      A[i,temp] = 1
    }
  }else{
    A = A[1:Nt,1:Nt]
  }
  CJ = A
  
  
  ##2  基础矩阵:类对角矩阵
  temp = Nt-1
  Atemp = diag(temp)
  A = matrix(0,nrow = Nt, ncol = Nt)     #类对角
  A[2:Nt,1:(Nt-1)] = Atemp
  A[upper.tri(A)] = t(A)[upper.tri(t(A))]
  LDJ = A     #类对角
  
  ##3  分块矩阵
  C = matrix(0,nrow=Nt,ncol=Nt)
  for (a in 0:(k-1)) {
    C[(a*gnum+1):((a+1)*gnum),(a*gnum+1):((a+1)*gnum)] = 1
    C[(a*gnum+1),(gnum*(a+1)+1)%%Nt] = 1
    C[(gnum*(a+1)+1)%%Nt,(a*gnum+1)] = 1
  }
  diag(C)=0
  FK = C     #分块对称矩阵
  
  
  ##4  ER network
  ernet = erdos.renyi.game(Nt, 0.01, type=c("gnp"), directed = FALSE, loops = FALSE)
  # plot(ernet)
  ernet = as.matrix(get.adjacency(ernet))
  ernet = ernet + LDJ
  
  ## ba network
  banet = sample_pa(Nt, power = 0.9, m = NULL, out.dist = NULL, out.seq = NULL,
                     out.pref = FALSE, zero.appeal = 1, directed = FALSE,
                     algorithm = c("psumtree", "psumtree-multiple", "bag"),
                     start.graph = NULL)
  # plot(banet)
  banet = as.matrix(get.adjacency(banet))
  
  # ws<-watts.strogatz.game(1,Nt,2,0.5)
  # plot(ws)
  
  #输出
  if(kind == 1){
    return(CJ)
  }else if(kind == 2){
    return(LDJ)
  }else if(kind == 3){
    return(FK)
  }else if(kind == 4){
    return(ernet)
  }else if(kind == 5){
    return(banet)
  }
  
}
