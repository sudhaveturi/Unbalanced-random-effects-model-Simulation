
rm(list=ls())


#### Parameters
m = c(2,5,10)
n = c(100,500)
bal = c("bal","mod","ext")
mu = 5
sigmaa2  = c(0.1,1,2.5,4,5)
sigma2 = 2.5
nIter = 1200
blup = list(list())
Est.mat = matrix(,nIter,23)

############################

pval_bal = pval_mod = pval_ext = rep(list())
pval_m1 = pval_m2 = pval_m3 = list(pval_bal,pval_mod,pval_ext)
pval_n1 = pval_n2 = list(pval_m1, pval_m2, pval_m3)

pval_sig1 = pval_sig2 = pval_sig3 = pval_sig4 = pval_sig5 = list(pval_n1,pval_n2)
pval_end = list(pval_sig1,pval_sig2,pval_sig3,pval_sig4,pval_sig5)

for(i in 1:length(sigmaa2))
{
  names(pval_end) = c("VarA=0.1","VarA=1","VarA=2.5","VarA=3","VarA=5")
  for(j in 1:length(n))
  {
    names(pval_end[[i]]) = c("N=100","N=500")
    for(k in 1:length(m))
    {
      names(pval_end[[i]][[j]]) = c("m=2","m=5","m=10")
      for(l in 1:length(bal))
      {   
        names(pval_end[[i]][[j]][[k]]) = c("Balanced","Moderate","Extreme")
      }
    }
  }
}
blup_end = pval_end
#####################################################

for(j in 1: length(n))
{
  for(w in 1:length(m))
  {
    for(k in 1:length(bal))
    {
      for(r in 1:length(sigmaa2))
      {       
        for(p in 1:nIter)
        {
          load(paste("/home/sveturi/BST760/PowerII/Pval+BLUP+Estimates_Power_",p,".rda",sep=""))
          Est.mat[p,]= pval_final[[r]][[j]][[w]][[k]]
          blup[[p]] = blup_final[[r]][[j]][[w]][[k]] 
          colnames(Est.mat) = c("mu","m","N","Type","sig2E","sig2A","MSA","MSE","mu.hat_r","mu.hat_f","sig2E.","sig2A.",
                             "sig2E.hat","sig2A.hat","f","pval_r","pval_f","optim_var(varE)","optim_var(varA)","n0",
                             "var(mu.hat)","var(sig2.hat)","var(siga2.hat)")
          
          print(paste("N(j)=",n[j],"m(w)=",m[w],"k=",bal[k],"varA=",sigmaa2[r],"iter=",p,sep=" "))
        }  
        pval_end[[r]][[j]][[w]][[k]] = Est.mat
        blup_end[[r]][[j]][[w]][[k]] = blup
      }
    }
  }
}


save(pval_end,blup_end, file="/home/sveturi/BST760/Pval+BLUP+Estimates_PowerII.rda")

###############################

load("/home/sveturi/BST760/Pval+BLUP+Estimates_PowerII.rda")
blup_edit = blup_end
nIter = 1000
for(j in 1: length(n))
{
  for(w in 1:length(m))
  {
    for(k in 1:length(bal))
    {
      for(r in 1:length(sigmaa2))
      {   
        t=which(is.na(pval_end[[r]][[j]][[w]][[k]][,19]))
        if(length(t)>0)
        {
          pval_end[[r]][[j]][[w]][[k]] = pval_end[[r]][[j]][[w]][[k]][-t,]
          for(z in 1:length(t))
          {
            blup_edit[[r]][[j]][[w]][[k]][[z]] = NULL
          }        
        }
        pval_end[[r]][[j]][[w]][[k]] = pval_end[[r]][[j]][[w]][[k]][1:nIter,]
        print(paste("N(j)=",n[j],"m(w)=",m[w],"k=",bal[k],"varA=",sigmaa2[r],sep=" "))
        
      }
    }
  }
}
pval_edit = pval_end
save(pval_edit, blup_edit,file="/home/sveturi/BST760/Pval+BLUP+Estimates_Edited_PowerII.rda")

###############################
m = c(2,5,10)
n = c(100,500)
bal = c("bal","mod","ext")
mu = 5
sigmaa2  = c(0.1,1,2.5,4,5)
sigma2 = 2.5
nIter = 1000

load("/home/sveturi/BST760/Pval+BLUP+Estimates_Edited_PowerII.rda")

compile.mat = matrix(,(length(n)*length(m)*length(bal)*length(sigmaa2)),28)

for(j in 1: length(n))
{
  for(w in 1:length(m))
  {
    for(k in 1:length(bal))
    {
      for(r in 1:length(sigmaa2))
      {   
        tmp = pval_edit[[r]][[j]][[w]][[k]][,c(1:6,20:23,16,17,7:10,13,14,18,19)]
        compile.mat[(j-1)*length(m)*length(bal)*length(sigmaa2)+(w-1)*length(sigmaa2)*length(bal)+(k-1)*length(sigmaa2)+r,13:20] = apply(matrix(as.numeric(tmp[,13:20]),nIter,8),2,mean)
        compile.mat[(j-1)*length(m)*length(bal)*length(sigmaa2)+(w-1)*length(sigmaa2)*length(bal)+(k-1)*length(sigmaa2)+r,21:28] = apply(matrix(as.numeric(tmp[,13:20]),nIter,8),2,var)
        compile.mat[(j-1)*length(m)*length(bal)*length(sigmaa2)+(w-1)*length(sigmaa2)*length(bal)+(k-1)*length(sigmaa2)+r,1:10] = tmp[1,1:10]
        compile.mat[(j-1)*length(m)*length(bal)*length(sigmaa2)+(w-1)*length(sigmaa2)*length(bal)+(k-1)*length(sigmaa2)+r,11] = (length(which(as.numeric(tmp[,11])<0.05))/nIter)*100
        compile.mat[(j-1)*length(m)*length(bal)*length(sigmaa2)+(w-1)*length(sigmaa2)*length(bal)+(k-1)*length(sigmaa2)+r,12] = (length(which(as.numeric(tmp[,12])<0.05))/nIter)*100
        colnames(compile.mat) = c("mu","m","N","Type","sig2E","sig2A","n0","var(mu.hat)","var(sig2.hat)",
                                  "var(siga2.hat)","power_r","power_f","mean(MSA)","mean(MSE)","mean(mu.hat_r)","mean(mu.hat_f)",
                              "mean(sig2E.hat)","mean(sig2A.hat)","mean(optim_var(varE))","mean(optim_var(varA))",
                              "var(MSA)","var(MSE)","var(mu.hat_r)","var(mu.hat_f)","var(sig2E.hat)","var(sig2A.hat)","var(optim_var(varE))",
                              "var(optim_var(varA))")                                                                                                                                                                                                                                                                          
      }
    }
  }
}

write.table(compile.mat, "/home/sveturi/BST760/CompiledEstimates_PowerII.txt",row.names=F,col.names=T,sep="\t")
#write.table(compile.mat, "/Users/Guha/Dropbox/UAB/Coursework/Spring2014/BST760/BST760_Project/Code/CompiledEstimates_PowerII.txt",row.names=F,col.names=T,sep="\t")
