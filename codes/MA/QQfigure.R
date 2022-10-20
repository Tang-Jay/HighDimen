#===============================================
#                  figure Q-Q             
#===============================================
source('GlambdaChen.R')
# #===================修改参数====================
nsim = 500            # 模拟次数
n = 200               # 样本个数
# size = c(42,55,72)    # p
m = 50               # 分位数个数
as = (1:m-0.5)/m
# #===================收集数据====================
# j = 0
# C = matrix(NA,nsim,length(size))
# colnames(C)=c('c1','c2','c3')
# for(p in size){
#   j  = j + 1
#   for(i in 1:nsim){
#     # 估计方程赋值
#     Z = matrix(rnorm(n*(p+1),0,1), p+1, n)
#     Z1 = Z[-(p+1),]
#     Z2 = Z[-1,]
#     X = Z1 + 0.5*Z2
#     z = t(X)
#     # 计算EL值
#     lam = lambdaChen(z)
#     el = 2*sum( log(1+t(lam)%*%t(z)) )
#     mel = (el - p) / sqrt(2*p)
#     C[i,j] = mel 
#   }
#   cat('n =',n,'p =',p,'完成模拟',nsim,'次\n')
# }
# write.csv(C, file = paste0('EL',n,p,'.csv'), row.names = FALSE) # 保存数据 
C = read.csv("EL20020.csv")
qas <- qnorm(as)                   # 正态分布的分位数点
c1 = sort(C[,1])[ceiling(as*nsim)] # mel的分位数点
c2 = sort(C[,2])[ceiling(as*nsim)]
c3 = sort(C[,3])[ceiling(as*nsim)]
# ===================绘制QQ图====================
x = seq(-3,3,0.01)
plot(x,x,
     xlim = c(-3,3), ylim = c(-3,3),
     xaxs = 'i', yaxs = 'i',
     yaxt = 'n',
     ann = F, type = 'l')
axis(2, las = 1)
title(main = paste0('Normal n = ',n),
      xlab= 'Normal quantile', ylab = 'EL quantile',
      line = 2)
lines(qas,c1,type='l',lty=1,lwd=1.5,col='blue1')
lines(qas,c2,type='l',lty=3,lwd=1.5,col='aquamarine4')
lines(qas,c3,type='l',lty=5,lwd=1.5,col='brown2')

