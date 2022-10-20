#===============================================
#             SARAR model N(0,1)
#===============================================
rm(list = ls()) 
source('GlambdaChen.R')
library('MASS')
library('sp')
library('sf')
library('spData')
library('terra')
library('spdep')

nsim = 100
rou1 = 0.85 
rou2 = 0.15
ms = c(10,13)

tol = 1e-4
a = 0.95

zero = 0
azero = 0
ff_el = c()
ff_ael = c()
for (m in ms){
  # SARAR模型
  n = m*m
  p = round(3*n^(0.4))
  beta = matrix(rep(1,p),p,1)
  
  mu = rep(0,p)
  Sigma = diag(p)
  for(k in 1:p){
    for(l in 1:p){
      if(k != l){Sigma[k,l] = 0.5}
    }
  }
  Xn = mvrnorm(n,mu,Sigma)
  cut = qchisq(a,p+3)
  
  # i = 1:n;Xn = i/(n+1)
  # Xn = rnorm(n)^2
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
  # 简化符号
  b = t(Xn)%*%t(Bn)
  g = diag(Gnn)
  h = diag(Hnn)
  s = Bn%*%Wn%*%Ani%*%Xn%*%beta
  f<-function(Matrix,Vector){
    irow = nrow(Matrix)
    v =c(0)
    for(i in 2:irow){
      v[i] =Matrix[i,][1:(i-1)]%*%Vector[1:i-1]
    }
    return(v)
  }
  if(1>log(n)/2) an=1 else an=log(n)/2
  # 启动模拟
  f1 = 0
  f2 = 0
  for(i in 1:nsim){
    # cat('样本个数为',n,'正在模拟第 ',i,'次','\n')
    En = rnorm(n,0,1)
    sigma2 = 1
    e = En
    # 模拟Yi(程序运行不需要Yi值)
    # Yn = Ani%*%Xn%*%beta + Ani%*%Bni%*%En
    # 估计方程赋值
    z = matrix(NA,nrow=n,ncol=p+3)
    z[,1:p] =  b*e
    z[,p+1] = g*(e^2-sigma2) + 2*e*f(Gnn,e) + s*e
    z[,p+2] = h*(e^2-sigma2) + 2*e*f(Hnn,e)
    z[,p+3] = e*e - rep(sigma2, n)
    
    # 计算EL值	
    lam = lambdaChen(z)
    el = 2*sum( log(1+t(lam)%*%t(z) ) )
    mel = (el-p-3)/sqrt(2*(p+3))
    if(mel<=qnorm(a)) f1=f1+1
    # aa=1+t(lam)%*%t(z)
    # glam=rowSums(t(z)/t(matrix(rep(aa,2),n,k+3)))
    # if(max(abs(glam))>tol) zero=zero+1
    
    
    # 计算AEL值
    az=rbind(z,-an*colMeans(z))	 		
    alam=lambdaChen(az)
    ael=2*sum( log(1+t(alam)%*%t(az) ) )  		
    if(ael<=cut) f2=f2+1
    # aa=1+t(alam)%*%t(az)
    # glam=rowSums(t(az)/t(matrix(rep(aa,2),n+1,k+3)))
    # if(max(abs(glam))>tol) azero=azero+1
  }
  # cat('样本个数为',n,'完成模拟',i,'次',zero,azero,'\n')
  # cat(paste0('N(0,',sqrt(sigma2),') ',n),' ',f1/nsim,f2/nsim,'\n')
  ff_el=append(ff_el,f1/nsim)
  ff_ael=append(ff_ael,f2/nsim)
  zero = 0
  azero = 0
  # Sys.sleep(10)
}
ff = matrix(NA,nrow=length(ms),ncol=2)
ff[,1]=ff_el
ff[,2]=ff_ael
rownames(ff) <- ms^2
colnames(ff) <- c("EL","AEL")
print(ff)

