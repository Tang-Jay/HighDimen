#===============================================
#                  figure.Q-Q             
#===============================================
rm(list = ls()) 
source('HighDim.R')
#===================修改参数====================
n = 200     # 产生n个样本数
size = c(33,44,58)
nsim = 500   # 模拟n次
m = 50       # 产生m个正态分位数
#===================收集数据====================
c1 = c()
c2 = c()
c3 = c()
as = (1:m-0.5)/m
for(a in as){
  mel = HighDim(a,n,size,nsim)
  c1 = c(c1,mel[1])
  c2 = c(c2,mel[2])
  c3 = c(c3,mel[3])
}
#===================绘制QQ图====================
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
lines(qnorm(as),c1,type='l',lty=5,col='brown2',lwd=1.5)
lines(qnorm(as),c2,type='l',lty=3,col='aquamarine4',lwd=1.5)
lines(qnorm(as),c3,type='l',lty=1,col='blue1',lwd=1.5)

