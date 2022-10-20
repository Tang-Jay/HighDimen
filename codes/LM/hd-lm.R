#===============================================
#           High-dim linear model             
#===============================================
rm(list = ls()) 
library('MASS')
source('GlambdaChen.R')
# ==================修改参数====================
nsim = 100  											  
n = 400
size = c(21,25,30)    # p
cut = qnorm(0.95)  
# ==================开始模拟====================
j = 0
C1 = matrix(NA,nsim,length(size));colnames(C1)=c('c1','c2','c3')
C2 = matrix(NA,nsim,length(size));colnames(C2)=c('c1','c2','c3')
for(p in size){
  j  = j + 1
  f1 = 0 
  f2 = 0
  EL = c()
  MEL = c()
  AEL = c()
  if(1>log(n)/2) an=1 else an=log(n)/2
  
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
    # Y = X%*%theta + e
    z = X*e
    
    # 计算MEL值
    lam = lambdaChen(z)                						
    el = 2*sum(log(1+t(lam)%*%t(z))) 
    mel = (el-p)/sqrt(2*p)
    if(mel<=cut) f1 = f1 + 1
    
    # 计算AEL值
    az=rbind(z,-an*colMeans(z))
    alam=lambdaChen(az)
    ael=2*sum( log(1+t(alam)%*%t(az)))
    ael = (ael-p)/sqrt(2*p)
    if(ael<=cut) f2 = f2 + 1
    
    EL = c(EL,el)
    MEL = c(MEL,mel)
    AEL = c(AEL,ael)
    
    C1[i,j] = mel 
    C2[i,j] = ael 
    
    # 检验lam
    # aa=1+t(lam)%*%t(z)
    # glam=rowSums(t(z)/matrix(aa,p,n,byrow = TRUE))
    # cat(m,max(abs(glam)),'\n')
  }
  cat('nsim =',nsim,'n =',n,'p =',p,'覆盖率',f1/nsim,f2/nsim,'\n')
}
write.csv(C1, file = paste0('MEL',n,p,'.csv'), row.names = FALSE) # 保存数据 
write.csv(C2, file = paste0('AEL',n,p,'.csv'), row.names = FALSE) # 保存数据 
# ===================绘制QQ图====================
m = 50                              # 分位数个数
as = (1:m-0.5)/m
qas <- qnorm(as)                    # 正态分布的分位数点
q <- function(C,j){
  c = sort(C[,j])[ceiling(as*nsim)] # mel的分位数点
  return(c)
}

par(mfrow = c(1, 2))
plot(qas,qas,xaxs = 'i', yaxs = 'i',
     yaxt = 'n', ann = F, type = 'l')
axis(2, las = 1)
title(main = paste0('n = ',n,' nsim = ',nsim),
      xlab= 'Normal quantile', ylab = 'EL quantile')
lines(qas,q(C1,1),type='l',lty=1,lwd=1.5,col='blue1')
lines(qas,q(C1,2),type='l',lty=3,lwd=1.5,col='aquamarine4')
lines(qas,q(C1,3),type='l',lty=5,lwd=1.5,col='brown2')

plot(qas,qas,xaxs = 'i', yaxs = 'i',
     yaxt = 'n', ann = F, type = 'l')
axis(2, las = 1)
title(main = paste0('n = ',n,' nsim = ',nsim),
      xlab= 'Normal quantile', ylab = 'EL quantile')
lines(qas,q(C2,1),type='l',lty=1,lwd=1.5,col='blue1')
lines(qas,q(C2,2),type='l',lty=3,lwd=1.5,col='aquamarine4')
lines(qas,q(C2,3),type='l',lty=5,lwd=1.5,col='brown2')

