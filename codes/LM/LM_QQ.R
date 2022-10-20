rm(list = ls()) 
LM_QQ <- function(EL, MEL){
  par(mfrow = c(1, 2))
  as = (1:50-0.5)/50                    # 名义水平
  
  # 卡方分位数图
  x <- seq(0,max(EL),5)                 # 卡方的分位数点
  c1 = sort(EL[,1])[ceiling(as*nsim)]   # mel的分位数点
  c2 = sort(EL[,2])[ceiling(as*nsim)]
  c3 = sort(EL[,3])[ceiling(as*nsim)]
  plot(x,x,xaxs = 'i', yaxs = 'i',
       xlim =c(min(c1), max(c3)),
       ylim =c(min(c1), max(c3)),
       yaxt = 'n', ann = F, type = 'l')
  axis(2, las = 1)
  title(main = paste0('ChiS n = ',n,' index = ',index),
        xlab = 'ChiSquare quantile', ylab = 'EL quantile')
  lines(qchisq(as,ps[1]),c1,type='l',lty=1,lwd=1.5,col='blue1')
  lines(qchisq(as,ps[2]),c2,type='l',lty=3,lwd=1.5,col='aquamarine4')
  lines(qchisq(as,ps[3]),c3,type='l',lty=5,lwd=1.5,col='brown2')
  points(qchisq(0.95,ps[1]),qchisq(0.95,ps[1]),pch=5,col='blue1')
  points(qchisq(0.95,ps[2]),qchisq(0.95,ps[2]),pch=5,col='aquamarine4')
  points(qchisq(0.95,ps[3]),qchisq(0.95,ps[3]),pch=5,col='brown2')
  
  # 正态分位数图
  x <- seq(-max(MEL),max(MEL),2)        # 正态的分位数点
  qas <- qnorm(as)                      # 正态分布的分位数点
  c1 = sort(MEL[,1])[ceiling(as*nsim)]  # mel的分位数点
  c2 = sort(MEL[,2])[ceiling(as*nsim)]
  c3 = sort(MEL[,3])[ceiling(as*nsim)]
  plot(x,x,xaxs = 'i', yaxs = 'i',
       xlim =c(min(c1), max(c3)),
       ylim =c(min(c1), max(c3)),
       yaxt = 'n', ann = F, type = 'l')
  axis(2, las = 1)
  title(main = paste0('Norm n = ',n,' index = ',index),
        xlab = 'Normal quantile', ylab = 'MEL quantile')
  lines(qas,c1,type='l',lty=1,lwd=1.5,col='blue1')
  lines(qas,c2,type='l',lty=3,lwd=1.5,col='aquamarine4')
  lines(qas,c3,type='l',lty=5,lwd=1.5,col='brown2')
  points(qnorm(0.95),qnorm(0.95),pch=5)
}
#===============================================
#                   赋值画图            
#=============================================== 
for(n in c(200,400,800)){
  for(index in c(0.4)){
    c = c(3,4,5)
    ps = round(c*n^index)
    MEL = read.csv(paste0('MEL',n,index,'.csv'))
    EL = read.csv(paste0('EL',n,index,'.csv'))
    nsim = dim(MEL)[1]
    LM_QQ(EL=EL,MEL=MEL)
  }
}