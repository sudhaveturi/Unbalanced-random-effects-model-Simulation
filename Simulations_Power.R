
rm(list=ls())
#setwd("/Users/sudhaveturi/Dropbox/UAB/Coursework/Spring2014/BST760/BST760_Project/Code")
setwd("C:/Users/Guha/Dropbox/UAB/Coursework/Spring2014/BST760/BST760_Project/Code")

library(lme4)
library(lmerTest)

#### Parameters
m = c(2,5,10)
n = c(100,500)
bal = c("bal","mod","ext")
mu = 5
sigmaa2  = c(0.01,0.1,1,3,5)
sigma2 = 1
nIter = 100
pval = matrix(,nIter,23)
blup = list(list())
pval_bal = pval_mod = pval_ext = rep(list())
pval_m1 = pval_m2 = pval_m3 = list(pval_bal,pval_mod,pval_ext)
pval_n1 = pval_n2 = list(pval_m1, pval_m2, pval_m3)

pval_sig1 = pval_sig2 = pval_sig3 = pval_sig4 = pval_sig5 = list(pval_n1,pval_n2)
pval_final = list(pval_sig1,pval_sig2,pval_sig3,pval_sig4,pval_sig5)

for(i in 1:length(sigmaa2))
{
  names(pval_final) = c("VarA=0.25","VarA=0.5","VarA=1","VarA=3","VarA=5")
  for(j in 1:length(n))
  {
    names(pval_final[[i]]) = c("N=100","N=500")
    for(k in 1:length(m))
    {
      names(pval_final[[i]][[j]]) = c("m=2","m=5","m=10")
      for(l in 1:length(bal))
      {   
        names(pval_final[[i]][[j]][[k]]) = c("Balanced","Moderate","Extreme")
      }
    }
  }
}
blup_final = pval_final
#######################

# levels = rep(list(list(list())), length(m))
# for(i in 1:length(m))
# {
#   for(j in 1: length(n))
#   {
#     a.bal <- sample(rep(1:m[i],each=(n[j]/m[i])))
#     a.mod <- sample(c(rep(1:round((m[i]/2),0), times=2*n[j]), sample((round((m[i]/2),0)+1):m[i], size=4*n[j], replace=TRUE)), size=n[j], replace=TRUE)
#     a.ext <- sample(c(rep(1:round((m[i]/2),0), times=2*n[j]), sample((round((m[i]/2),0)+1):m[i], size=n[j], replace=TRUE)), size=n[j], replace=TRUE)
#     print(as.numeric(table(a.bal)))
#     print(as.numeric(table(a.mod)))
#     print(as.numeric(table(a.ext)))
#     levels[[i]][[j]] = list(a.bal,a.mod,a.ext)
#   }
# }
#save(levels, file="Levels.rda")

load("Levels.rda")

# i = 1
# j = 1
# k = 1
# r = 1


for(j in 1: length(n))
{
for(i in 1:length(m))
{
    for(k in 1:length(bal))
    {
      for(r in 1:length(sigmaa2))
      {
        
        a<-as.factor(levels[[i]][[j]][[k]])
        ################# Expected Variances of parameters
        var_mu.hat = 1/(sum((table(a)/(sigma2+(table(a)*sigmaa2[r])))))
        lambdai = sigma2+(table(a)*sigmaa2[r])
        D = ((n[j]-m[i])/sigma2^2)*(sum(((table(a)^2)/(lambdai^2))))+sum(1/(lambdai^2))*sum(((table(a))^2/(lambdai^2))) - (sum(table(a)/lambdai^2))^2
        var_sigma2.hat = (2/D)*sum((table(a)^2)/(lambdai^2))
        var_sigmaa2.hat = (2/D)*(((n[j]-m[i])/sigma2^2) + sum(1/lambdai^2))
                        
        ################## Start rep loop
        
        for(p in 1:nIter)
        {      
         e0<-rnorm(n[j],0,sqrt(sigma2))
         ef_a<-rnorm(m[i],0,sqrt(sigmaa2[r]))
      
      

         X<-matrix(nc=m[i], nr=n[j])
      
        for(q in 1:ncol(X)){X[,q]<-as.numeric(a==q)}
 #       X[1:5, 1:5]
      
 #        colSums(X)
         yA <-mu+X%*%ef_a+e0
         
         ########## Model fitting
         
         rand_alt <- lmer(yA ~ 1+(1|a), REML=F)
         fixed_alt <- lm(yA ~ a)
         dat = cbind(yA, a, e0)
         colnames(dat) = c("yA" , "a", "e0")
 
        ybari. = aggregate(yA~a,dat,mean)$yA
        ybari._m = as.numeric(a)
        for(z in 1:m[i])
          ybari._m[ybari._m == levels(a)[z]] <- ybari.[z]
        SSE = sum((yA-ybari._m)^2)
 
       fn <- function(logVar) { 
         varE<-exp(logVar[1])
         varA<-exp(logVar[2])
         lambdai. = as.numeric(varE+(varA*table(a)))
         mu. = sum(aggregate(yA~a,dat,sum)$yA/lambdai.)/(sum(table(a)/lambdai.))
         loglik = - 0.5*sum(log(lambdai.)) - 0.5*(n[j]-m[i])*log(varE)-(SSE/(2*varE))-sum((table(a)*(ybari.-mu.)^2)/(2*lambdai.))
         out = -2*loglik 
         return(out)       
       }
       fmML = optim(fn, par=log(c(0.2,0.8)*var(yA)),hessian=T)
       sigma2. = exp(fmML$par)[1]
       sigmaa2. = exp(fmML$par)[2]
       lambdai. = as.numeric(sigma2.+(sigmaa2.*table(a)))
       mu. = sum(aggregate(yA~a,dat,sum)$yA/lambdai.)/(sum(table(a)/lambdai.))
          
 
       if(sum(fmML$hessian[,2]) > 10^-8)
      {
        var = diag(solve(fmML$hessian))
        var.sigma2 = var[1] 
        var.sigmaa2 = var[2] 
      } else {
        var.sigma2 = NA 
        var.sigmaa2 = NA
      }
 
 
        if(sigmaa2.<0)
        {
          sigma2.hat = sum((yA-mean(yA))^2)/n[j]
          sigmaa2.hat = 0
          lambdai.hat = as.numeric(sigma2.hat+(sigmaa2.hat*table(a)))
          mu.hat = mean(yA)

        } else {
          sigma2.hat = sigma2.
          sigmaa2.hat = sigmaa2.
          lambdai.hat = lambdai.
          mu.hat = mu.
        }
          
   #      summary(rand_alt)  
          
      ################### P-value distribution
          

          SSA = sum(((ybari._m-mu.hat)^2))
          SST = sum((yA-mu.hat)^2)
          MSE = SSE/(n[j]-m[i])
          MSA = SSA/(m[i]-1)
          n0 = (1/(m[i]-1))*(n[j]-((sum(table(a)^2))/n[j]))          
          f = pf((MSA/MSE),m[i]-1,n[j]-m[i])
          
          
          ################### BLUP generation
          BP_rand = (sigmaa2.hat/(sigmaa2.hat+(sigma2.hat/n[j])))*(ybari.-mu.hat)
          BP_fixed = ybari.- summary(fixed_alt)[[4]][1]
          pred = cbind(BP_rand,BP_fixed)
          colnames(pred) = c("BLUP","BLUE")
          blup[[p]] = pred
 
          ################# Generate p val matrix
          
          pval[p,1] = mu
          pval[p,2] = m[i]
          pval[p,3] = n[j]
          pval[p,4] = bal[k]
          pval[p,5] = sigma2
          pval[p,6] = sigmaa2[r]
          pval[p,7] = MSA
          pval[p,8] = MSE
          pval[p,9] = mu.hat
          pval[p,10] = summary(fixed_alt)[[4]][1]
          pval[p,11] = sigma2.
          pval[p,12] = sigmaa2.
          pval[p,13] = sigma2.hat
          pval[p,14] = sigmaa2.hat
          pval[p,15] = f
          pval[p,16] = 2*min(f,1-f)
          pval[p,17] = as.numeric(unlist(summary(aov(yA ~ a)))[9])
          pval[p,18] = var.sigma2
          pval[p,19] = var.sigmaa2
          pval[p,20] = n0
          pval[p,21] = var_mu.hat
          pval[p,22] = var_sigma2.hat
          pval[p,23] = var_sigmaa2.hat
         colnames(pval) = c("mu","m","N","Type","sig2E","sig2A","MSA","MSE","mu.hat_r","mu.hat_f","sig2E.","sig2A.",
                            "sig2E.hat","sig2A.hat","f","pval_r","pval_f","optim_var(varE)","optim_var(varA)","n0",
                            "var(mu.hat)","var(sig2.hat)","var(siga2.hat)")
          
         print(paste("N(j)=",n[j],"m(i)=",m[i],"k=",bal[k],"varA=",sigmaa2[r],"iter=",p,sep=" "))
        }
        pval_final[[r]][[j]][[i]][[k]] = pval
        blup_final[[r]][[j]][[i]][[k]] = blup
      }      
    }
  }
}
save(pval_final,blup_final, file="Pval+BLUP+Estimates_Power.rda")