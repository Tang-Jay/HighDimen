#================================================
#                     LM覆盖率
#================================================
rm(list = ls()) 
LM_CP <- function(ns=ns, indexs=indexs, a=0.95){
  # ns = c(200,400,800)
  # indexs = c(0,0.16,0.24,0.4)
  # 预存空间
  c = c(3,4,5)
  CPs = matrix(NA,length(indexs)*3,length(ns)*3)
  colnames(CPs) = rep(ns,each=3)
  rownames(CPs) = rep(indexs,each=3)
  # 计算覆盖率
  i = 0
  for(n in ns){
    i = i + 3
    j = 0
    for(index in indexs){
      j = j + 3
      # 计算p
      ps = round(c*n^index)
      # 分位数
      q1 = qchisq(a,ps)
      q2 = qnorm(0.5+0.5*a)
      # 读取数据
      EL = read.csv(paste0('EL',n,index,'.csv'))
      MEL = read.csv(paste0('MEL',n,index,'.csv'))
      nsim = dim(EL)[1]
      # 编制分位数矩阵
      Q1 = matrix(q1,nsim,3,byrow=TRUE)
      Q2 = matrix(q2,nsim,3)
      # 计算覆盖率
      CP1 = colSums(EL<=Q1)/nsim
      CP2 = colSums(abs(MEL)<=Q2)/nsim
      
      CPs[(j-2):j,i-2] = ps
      CPs[(j-2):j,i-1] = CP1
      CPs[(j-2):j,i  ] = CP2
    }
  }
  # 输出覆盖率
  cat('名义水平',a,'\n')
  print(CPs)
  
  
  # latex版输出
  library(xtable)
  print(xtable(CPs,auto=TRUE,
               label=paste0('tab-',a)
               )
        )
  
  
  # 本地版输出
  # write.csv( CPs, file = paste0('SARAR_CPs_',error,'_',a,'.csv'))
  
  # 网页版输出
  # library(kableExtra)
  # kable(CPs, "html", align = 'c') %>% 
  #   kable_styling(full_width = F,
  #     bootstrap_options = c("striped", "hover", "condensed", "responsive")) 
  cat('\n')
}
#===============================================
#                   赋值运算            
#===============================================   
ns = c(200,400,800)
indexs = c(0,0.14,0.16,0.24,0.4,0.5)

for(a in c(0.9)){
  LM_CP(ns=ns, indexs=indexs, a=a)
}
