#===============================================
#             function for high-dim el             
#===============================================
# nsim = 100
# n = 800
# size = c(42,55,72)
# a = 0.95
#==================函数主体====================#
HighDim<-function(a, n, size, nsim = 200){
  # 函数准备
  source('GlambdaChen.R')
  modified<-function(statistic){
    modified_statistic = (statistic-p)/sqrt(2*p)
    return(modified_statistic)
  }
  # 参数赋值
  rho = 0.5
  re = c() # 保存函数返回的结果
  j = 0 # 记录循环次数
  EL = matrix(NA,nrow=length(size),ncol=5) # 记录每次实验数据
  colnames(EL)=c('n','p','el','mel','mel2')
  
  # 开始模拟
  cat('模型进行',nsim,'次','\n')
  cat('n  ','p ',' EL ','MEL1 ','MEL2 ','\n')
  for (p in size){
    j = j+1
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
      
    }
    cat(n,p,f1/nsim,f2/nsim,f3/nsim,'\n')
    EL[j,] = c(n,p,f1/nsim,f2/nsim,f3/nsim)
    re = c(re,f3/nsim)
    
  }
  write.csv(EL, file = paste0('EL',n,p,'.csv'), row.names = FALSE)
  return(re)
}
    

