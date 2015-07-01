
#### Parameters
m = c(2,5,10)
n = c(100,500)
bal = c("bal","mod","ext")
mu = 5
sigma2 = 1
nIter = 1000
num.reps = 100
pval = rep(list(list()),nIter)
######################################

pval_bal = pval_mod = pval_ext = rep(list())
pval_m1 = pval_m2 = pval_m3 = list(pval_bal,pval_mod,pval_ext)
pval_n1 = pval_n2 = list(pval_m1, pval_m2, pval_m3)

pval_end = list(pval_n1,pval_n2)

for(j in 1:length(n))
{
  names(pval_end) = c("N=100","N=500")
  for(k in 1:length(m))
  {
    names(pval_end[[j]]) = c("m=2","m=5","m=10")
    for(l in 1:length(bal))
    {   
      names(pval_end[[j]][[k]]) = c("Balanced","Moderate","Extreme")
    }
  }
}

#######################

for(j in 1: length(n))
{
  for(w in 1:length(m))
  {
    for(k in 1:length(bal))
    {      
      for(p in 1:nIter)
        {

          load(paste("/home/sveturi/BST760/Type1Error/Pval_Type1Error_",p,".rda",sep=""))
          pval[[p]]= pval_final[[j]][[w]][[k]][[p]]
          print(paste("N(j)=",n[j],"m(w)=",m[w],"k=",bal[k],"iter=",p,sep=" "))
        }  
      pval_end[[j]][[w]][[k]] = pval
      }
    }
  }


save(pval_end, file="/home/sveturi/BST760/Pval_Type1Error_full.rda")


##########################

#### Parameters
m = c(2,5,10)
n = c(100,500)
bal = c("bal","mod","ext")
mu = 5
sigma2 = 1
nIter = 1000
num.reps = 100

load("/Users/sudhaveturi/Dropbox/UAB/Coursework/Spring2014/BST760/BST760_Project/Code/Pval_Type1Error_full.rda")

pval.mat_r = pval.mat_f = matrix(,(length(n)*length(m)*length(bal)),num.reps+3)


colnames(pval.mat_r) = colnames(pval.mat_f) = c("N","m","Type",1:num.reps)

for(j in 1: length(n))
{
  for(w in 1:length(m))
  {
    for(k in 1:length(bal))
    {
      pval.mat_r[(j-1)*length(m)*length(bal)+(w-1)*length(bal)+k,1] = n[j]
      pval.mat_r[(j-1)*length(m)*length(bal)+(w-1)*length(bal)+k,2] = m[w]
      pval.mat_r[(j-1)*length(m)*length(bal)+(w-1)*length(bal)+k,3] = bal[k]
      
      pval.mat_f[(j-1)*length(m)*length(bal)+(w-1)*length(bal)+k,1] =n[j]
      pval.mat_f[(j-1)*length(m)*length(bal)+(w-1)*length(bal)+k,2] =m[w]
      pval.mat_f[(j-1)*length(m)*length(bal)+(w-1)*length(bal)+k,3] =bal[k]
      
      for(i in 1:num.reps)
      {
        tmp_r = tmp_f = 0
        for(t in 1:nIter)
        {
          tmp = pval_end[[j]][[w]][[k]]
          tmp_rand = length(which(tmp[[t]][i,1]<0.05))
          tmp_r = tmp_r+tmp_rand
          tmp_fixed = length(which(tmp[[t]][i,2]<0.05))
          tmp_f = tmp_f+tmp_fixed
        }
        pval.mat_r[(j-1)*length(m)*length(bal)+(w-1)*length(bal)+k,3+i] = (tmp_r/nIter)*100
        pval.mat_f[(j-1)*length(m)*length(bal)+(w-1)*length(bal)+k,3+i] = (tmp_f/nIter)*100
        
      }
    }      
  }
}
save(pval.mat_r,pval.mat_f,file="/Users/sudhaveturi/Dropbox/UAB/Coursework/Spring2014/BST760/BST760_Project/Code/Pval_rand+fix_Type1Error.rda")
#write.table(pval.mat,"/home/sveturi/BST760/Type1Error/Pval_Type1Error.txt",row.names=F,col.names=T,sep="\t")