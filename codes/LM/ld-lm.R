#===============================================
#           Low-dim linear model             
#===============================================
rm(list = ls()) 
source('GlambdaChen.R')
#==================修改参数====================#
nsim = 2000  											#模拟次数
a = 0.95    											#置信度
k = 2       											#卡方自由度
n = 500      											#样本数量
cut = qchisq(a, k)         	      #卡方分位数
#==================开始模拟====================#
f = 0 
EL = c()
for(m in 1:nsim){ 
  # 估计方程赋值
  x0 = rnorm(n)                 					
  x1 = c(x0,x0^2)               					
  x = matrix(x1,nrow=n, ncol=2)  			    
  e = runif(n,-0.5,0.5)          				
  e = t(t(e))
  e1=matrix(c(e,e),nrow=n, ncol=2) 		    
  z=x*e1  
  # 计算EL值
  lam=lambdaChen(z)                						
  el=sum(2*log(1+x%*%lam*e))   			      
  if(el<=cut) f=f+1
  EL = c(EL,el)
}
ff=f/nsim
cat('覆盖率',ff,'\n')
# ===================绘制QQ图====================
m = 50                             # 分位数个数
as = (1:m-0.5)/m
qas <- qchisq(as,k)                # 正态分布的分位数点
qel = sort(EL)[ceiling(as*nsim)]   # mel的分位数点

plot(qas,qas,
     xaxs = 'i', yaxs = 'i',
     yaxt = 'n',
     ann = F, type = 'l')
axis(2, las = 1)
title(main = paste0('ChiSquare n = ',n),
      xlab= 'Chi_square quantile', ylab = 'EL quantile',
      line = 2)
points(cut,cut,type = "p",pch=5,col='brown2',cex=0.9)
points(qas,qel,type = "p",pch=20,col='blue1',cex=0.5)

