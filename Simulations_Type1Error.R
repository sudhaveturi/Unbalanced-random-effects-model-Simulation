
rm(list=ls())

### LOADING ARGUMENTS, SETTING VARIABLES ############
if(TRUE){
  args=(commandArgs(TRUE))
  if(length(args)==0){
    print('no args')
    i=1; full=TRUE
  }else{
    print('eval agrs')
    for(i in 1:length(args)){
      eval(parse(text=args[[i]]))
    }
  }}


library(lme4)
library(lmerTest)

#### Parameters
m = c(2,5,10)
n = c(100,500)
bal = c("bal","mod","ext")
mu = 5
sigma2 = 1
#nIter = 100
num.reps = 100
pval = matrix(,1,23)
blup = list(list())
pval_bal = pval_mod = pval_ext = rep(list())
pval_m1 = pval_m2 = pval_m3 = list(pval_bal,pval_mod,pval_ext)
pval_n1 = pval_n2 = list(pval_m1, pval_m2, pval_m3)

pval_final = list(pval_n1,pval_n2)

  for(j in 1:length(n))
  {
    names(pval_final[[w]]) = c("N=100","N=500")
    for(k in 1:length(m))
    {
      names(pval_final[[w]][[j]]) = c("m=2","m=5","m=10")
      for(l in 1:length(bal))
      {   
        names(pval_final[[w]][[j]][[k]]) = c("Balanced","Moderate","Extreme")
      }
    }
  }

#######################

# levels = rep(list(list(list())), length(m))
# for(w in 1:length(m))
# {
#   for(j in 1: length(n))
#   {
#     a.bal <- sample(rep(1:m[w],each=(n[j]/m[w])))
#     a.mod <- sample(c(rep(1:round((m[w]/2),0), times=2*n[j]), sample((round((m[w]/2),0)+1):m[w], size=4*n[j], replace=TRUE)), size=n[j], replace=TRUE)
#     a.ext <- sample(c(rep(1:round((m[w]/2),0), times=2*n[j]), sample((round((m[w]/2),0)+1):m[w], size=n[j], replace=TRUE)), size=n[j], replace=TRUE)
#     print(as.numeric(table(a.bal)))
#     print(as.numeric(table(a.mod)))
#     print(as.numeric(table(a.ext)))
#     levels[[w]][[j]] = list(a.bal,a.mod,a.ext)
#   }
# }
#save(levels, file="Levels.rda")

load("Levels.rda")

# w = 1
# j = 2
# k = 1
# r = 1


for(j in 1: length(n))
{
  for(w in 1:length(m))
  {
    for(k in 1:length(bal))
    {
      
        a<-as.factor(levels[[w]][[j]][[k]])
       
        ################## Start rep loop
        
        #        for(i in 1:nIter)
        #        {      
        e0<-rnorm(n[j],0,sqrt(sigma2))
        
        yA <-mu+e0
        
        pv_rand = pv_fixed = matrix(,num.reps,2)
        
        for (m = 1:num.reps)
        {
          yA = sample(yA)
        ########## Model fitting
        
        #rand_alt <- lmer(yA ~ 1+(1|a), REML=F)
 #       fixed_alt <- lm(yA ~ a)
        dat = cbind(yA, a, e0)
        colnames(dat) = c("yA" , "a", "e0")
        
        ybari. = aggregate(yA~a,dat,mean)$yA
        ybari._m = as.numeric(a)
        for(z in 1:m[w])
          ybari._m[ybari._m == levels(a)[z]] <- ybari.[z]
        SSE = sum((yA-ybari._m)^2)
        
               
        ################### P-value distribution
        
        mu.hat = (sum(table(a)*ybari.))/n[j]
        SSA = sum(((ybari._m-mu.hat)^2))
        SST = sum((yA-mu.hat)^2)
        MSE = SSE/(n[j]-m[w])
        MSA = SSA/(m[w]-1)
        n0 = (1/(m[w]-1))*(n[j]-((sum(table(a)^2))/n[j]))          
        f = pf((MSA/MSE),m[w]-1,n[j]-m[w])
        pv_rand = 2*min(f,1-f)
        
        pv_fixed = as.numeric(unlist(summary(aov(yA ~ a)))[9])
 
    }
        
             ################# Generate p val matrix
            
        #          print(paste("N(j)=",n[j],"m(w)=",m[w],"k=",bal[k],"varA=",sigmaa2[r],"iter=",i,sep=" "))
        #        }
        pval_final[[j]][[w]][[k]] = pval
      }      
    }
  }
}

save(pval_final,blup_final, file=paste("/home/sveturi/BST760/Type1Error/Pval_Type1Error_",i,".rda",sep=""))

