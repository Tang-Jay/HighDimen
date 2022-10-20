#===============================================
#                 SARAR经验似然统计量
#===============================================
rm(list = ls()) 
source('GlambdaChen.R')
library('MASS')
library('sp')
library('sf')
library('spData')
library('terra')
library('spdep')
SARAR_Data <- function(m,index,error=1,nsim = 1000){
  # ==================修改参数====================
  # m = 20
  # index = 0.14
  # error: 1->N(0,1) 2->N(0,0.75) 3->t(5) 4->chi2(4)
  # nsim = 1000
  
  c = c(3,4,5)
  rou1 = 0.85 
  rou2 = 0.15
  n = m*m
  ps = round(c*n^index)
  # ==================函数准备====================
  f<-function(Matrix,Vector){
    irow = nrow(Matrix)
    v =c(0)
    for(i in 2:irow){
      v[i] =Matrix[i,][1:(i-1)]%*%Vector[1:i-1]
    }
    return(v)
  }
  # ==================数据准备====================
  Wnb = cell2nb(m,m,type='queen')
  Ws = nb2listw(Wnb)
  Wn = listw2mat(Ws)
  Mn = Wn
  In = diag(n)
  An = In - rou1*Wn
  Bn = In - rou2*Mn
  Ani = solve(An)
  Bni = solve(Bn)
  # 估计方程准备
  Gn = Bn%*%Wn%*%Ani%*%Bni
  Gnn = 1/2*(Gn + t(Gn))
  Hn = Mn%*%Bni
  Hnn = 1/2*(Hn + t(Hn))
  g = diag(Gnn)
  h = diag(Hnn)
  # ==================开始模拟====================
  EL = matrix(NA,nsim,length(ps))
  colnames(EL)=c(paste0('p=',ps[1]),paste0('p=',ps[2]),paste0('p=',ps[3]))
  
  MEL = matrix(NA,nsim,length(ps))
  colnames(MEL)=c(paste0('p=',ps[1]),paste0('p=',ps[2]),paste0('p=',ps[3]))
  
  cat('nsim =',nsim,'error =',error,'index =',index,'\n')
  cat('n  ','p ','EL ',' MEL ','\n')
  j = 0
  for(p in ps){
    j  = j + 1
    mu = rep(0,p)
    Sigma = diag(p)
    # for(k in 1:p){
    #   for(l in 1:p){
    #     if(k != l){Sigma[k,l] = 0.5}
    #   }
    # }
    beta = matrix(rep(1,p),p,1)
    
    f1 = 0
    f2 = 0
    for(i in 1:nsim){
      Xn = mvrnorm(n,mu,Sigma)
      
      if(error == 1){En = rnorm(n,0,1);sigma2 = 1}
      if(error == 2){En = rnorm(n,0,sqrt(0.75));sigma2 = 0.75}
      if(error == 3){En = rt(n,5);sigma2=5/(5-2)}
      if(error == 4){En = rchisq(n,4)-4;sigma2=8}
      
      e = En
      b = t(Xn)%*%t(Bn)
      s = Bn%*%Wn%*%Ani%*%Xn%*%beta
      # 模拟Yi(程序运行不需要Yi值)
      # Yn = Ani%*%Xn%*%beta + Ani%*%Bni%*%En
      # 估计方程赋值
      z = matrix(NA,nrow=n,ncol=p+3)
      z[,1:p] = b*e
      z[,p+1] = g*(e^2-sigma2) + 2*e*f(Gnn,e) + s*e
      z[,p+2] = h*(e^2-sigma2) + 2*e*f(Hnn,e)
      z[,p+3] = e*e - rep(sigma2, n)
      
      # 计算EL值	
      lam = lambdaChen(z)
      el = 2*sum( log(1+t(lam)%*%t(z) ) )
      mel = (el-p-3)/sqrt(2*(p+3))
      if(  el<=qchisq(0.95,p+3)) f1=f1+1
      if(abs(mel)<=qnorm(0.975)) f2=f2+1
      
      EL[i,j]  =  el
      MEL[i,j] = mel
    }
    cat(n,p,f1/nsim,f2/nsim,'\n')
  }
  write.csv( EL, file = paste0('EL' ,error,n,index,'.csv'), row.names = FALSE) # 保存数据
  write.csv(MEL, file = paste0('MEL',error,n,index,'.csv'), row.names = FALSE) # 保存数据
}
#===============================================
#                   赋值运算            
#=============================================== 
index = 0.5
ms = c(15,20,30)
for(m in ms){
  for(error in 1:4){
    SARAR_Data(m,index,error)
  }
}


