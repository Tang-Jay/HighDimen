#===============================================
#             High dimensional Data             
#===============================================
rm(list = ls()) 
source('GlambdaChen.R')
source('GlambdaQin.R')
modified<-function(statistic){
  modified_statistic = (statistic-p)/sqrt(2*p)
  return(modified_statistic)
}
#==================修改参数====================#
nsim = 100

n = 800
size = c(42,55,72)

rho = 0.5
a = 0.95
#==================开始模拟====================#
cat('模型进行',nsim,'次','\n')
cat('n  ','p ',' EL ','MEL1 ','MEL2 ','\n')
for (p in size){
  mu = rep(0,p)
  Sigma = diag(p) 
  for(i in 1:p){
    for(j in 1:p){
      if(abs(i-j)==1){Sigma[i,j] = rho}
    }
  }
  
  f1 = 0
  f2 = 0
  f3 = 0
  for(i in 1:nsim){
    # 估计方程赋值
    Z = matrix(rnorm(n*(p+1),0,1), p+1, n)
    Z1 = Z[-(p+1),]
    Z2 = Z[-1,]
    X = Z1 + rho*Z2
    z = t(X)
    
    # 计算EL值
    lam = lambdaChen(z)
    el = 2*sum( log(1+t(lam)%*%t(z)) )
    mel = modified(el)
    if(  el <= qchisq(a,p) )              f1=f1+1
    if( mel <= qnorm(a,0,1))              f2=f2+1
    if( abs(mel) <= qnorm(a+(1-a)/2,0,1)) f3=f3+1
    
    # aa=1+t(lam)%*%t(z)
    # glam=rowSums(t(z)/matrix(aa,p,n,byrow = TRUE))
    # print(max(abs(glam)))
    
    # cat(el,mel,'\n')
  }
  # cat('n =',n,'p =',p,'EL =',f1/nsim,'MEL1 =',f2/nsim,'MEL2 =',f3/nsim,'\n')
  cat(n,p,f1/nsim,f2/nsim,f3/nsim,'\n')
}

