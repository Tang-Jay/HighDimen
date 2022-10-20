#===============================================
#           High-dim linear model             
#===============================================
rm(list = ls()) 
library('MASS')
source('GlambdaChen.R')
LM_Data <- function(n,index,nsim=1000){
  # 给定参数
  # n = 800
  # index = 0.24
  # nsim = 1000  
  c = c(3,4,5)
  ps = round(c*n^index)    
  # 开始模拟
  EL = matrix(NA,nsim,length(ps))
  colnames(EL)=c(ps[1],ps[2],ps[3])
  
  MEL = matrix(NA,nsim,length(ps))
  colnames(MEL)=c(ps[1],ps[2],ps[3])
  
  cat('nsim =',nsim,'n =',n,'index =',index,'\n')
  j = 0
  for(p in ps){
    j  = j + 1
    f1 = 0
    f2 = 0
    
    mu = rep(0,p)
    Sigma = diag(p)
    for(k in 1:p){
      for(l in 1:p){
        if(k != l){Sigma[k,l] = 0.5}
      }
    }
    
    for(i in 1:nsim){
      # 估计方程赋值
      X = mvrnorm(n,mu,Sigma)
      e = rnorm(n)
      # beta0 = rep(1,p)
      # Y = X%*%beta0 + e
      # z = X*(Y-X%*%beta0)
      z = X*e
      
      # 计算MEL值
      lam = lambdaChen(z)
      el = 2*sum(log(1+t(lam)%*%t(z)))
      if(el<=qchisq(0.95,p)) f1 = f1 + 1
      mel = (el-p)/sqrt(2*p)
      if(abs(mel)<=qnorm(0.975)) f2 = f2 + 1
      
      # 检验lam
      # aa=1+t(lam)%*%t(z)
      # glam=rowSums(t(z)/matrix(aa,p,n,byrow = TRUE))
      # cat(m,max(abs(glam)),'\n')
      
      EL[i,j]  =  el
      MEL[i,j] = mel
    }
    cat('n =',n,'p =',p,'EL =',f1/nsim,'MEL =',f2/nsim,'\n')
  }
  write.csv( EL, file = paste0('EL' ,n,index,'.csv'), row.names = FALSE) 
  write.csv(MEL, file = paste0('MEL',n,index,'.csv'), row.names = FALSE) 
  cat('\n')
}
#===============================================
#                   赋值运算            
#=============================================== 
index = 0.5
ns = c(800)
for(n in ns){
  LM_Data(n=n,index=index)
}
