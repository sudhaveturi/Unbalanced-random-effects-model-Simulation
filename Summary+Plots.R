
rm(list=ls())


#### Parameters
m = c(2,5,10)
n = c(100,500)
bal = c("bal","mod","ext")
mu = 5
sigmaa2  = c(0.1,1,2.5,4,5)
sigma2 = 2.5
nIter = 1000
num.reps = 100

############################

Est.mat = read.table("C:/Users/Guha/Dropbox/UAB/Coursework/Spring2014/BST760/BST760_Project/Code/CompiledEstimates_PowerII.txt",header=T,sep="\t")
load("C:/Users/Guha/Dropbox/UAB/Coursework/Spring2014/BST760/BST760_Project/Code/Pval_rand+fix_Type1Error.rda")

#Est.mat = read.table("/Users/sudhaveturi/Dropbox/UAB/Coursework/Spring2014/BST760/BST760_Project/Code/CompiledEstimates_PowerII.txt",header=T,sep="\t")
load("C:/Users/Guha/Dropbox/UAB/Coursework/Spring2014/BST760/BST760_Project/Code/Pval+BLUP+Estimates_Edited_PowerII.rda")
#load("/Users/sudhaveturi/Dropbox/UAB/Coursework/Spring2014/BST760/BST760_Project/Code/Pval+BLUP+Estimates_Edited_PowerII.rda")
type1_long_r = matrix(,(nrow(pval.mat_r)*num.reps),5)
for(i in 1:num.reps)
{
  type1_long_r[(((i-1)*nrow(pval.mat_r))+1):(i*nrow(pval.mat_r)),5] = pval.mat_r[,i+3]
  type1_long_r[(((i-1)*nrow(pval.mat_r))+1):(i*nrow(pval.mat_r)),4] = i
  for(j in 1:3)
  type1_long_r[(((i-1)*nrow(pval.mat_r))+1):(i*nrow(pval.mat_r)),j] = as.character(pval.mat_r[,j])
}

type1_long_r[,3] = sub("bal","1=Balanced",type1_long_r[,3])
type1_long_r[,3] = sub("ext","3=Extreme",type1_long_r[,3])
type1_long_r[,3] = sub("mod","2=Moderate",type1_long_r[,3])
colnames(type1_long_r) = c("N","Groups","Type","Reps","Type1_Error")


type1_long_f = matrix(,(nrow(pval.mat_f)*num.reps),5)
for(i in 1:num.reps)
{
  type1_long_f[(((i-1)*nrow(pval.mat_f))+1):(i*nrow(pval.mat_f)),5] = pval.mat_f[,i+3]
  type1_long_f[(((i-1)*nrow(pval.mat_f))+1):(i*nrow(pval.mat_f)),4] = i
  for(j in 1:3)
    type1_long_f[(((i-1)*nrow(pval.mat_f))+1):(i*nrow(pval.mat_f)),j] = as.character(pval.mat_f[,j])
}



type1_long_f[,3] = sub("bal","1=Balanced",type1_long_f[,3])
type1_long_f[,3] = sub("ext","3=Extreme",type1_long_f[,3])
type1_long_f[,3] = sub("mod","2=Moderate",type1_long_f[,3])
colnames(type1_long_f) = c("N","Groups","Type","Reps","Type1_Error")

type1_long_r = as.data.frame(type1_long_r)
type1_long_f = as.data.frame(type1_long_f)



head(type1_long_r)
head(type1_long_f)

Est.mat[,4] = sub("bal","1=Balanced",Est.mat[,4])
Est.mat[,4] = sub("mod","2=Moderate",Est.mat[,4])
Est.mat[,4] = sub("ext","3=Extreme",Est.mat[,4])
Est.mat[,2] = as.factor(Est.mat[,2])

library(lattice)
library(ggplot2)

Bias.mat = matrix(,nrow(Est.mat),5)
Bias.mat[,1] = (Est.mat$mu-Est.mat$mean.mu.hat_r.)
Bias.mat[,2] = (Est.mat$mu-Est.mat$mean.mu.hat_f.)
Bias.mat[,3] = (Est.mat$sig2A-Est.mat$mean.sig2A.hat.)
Bias.mat[,4] = (Est.mat$var.mu.hat.-Est.mat$var.mu.hat_r.)
Bias.mat[,5] = (Est.mat$var.siga2.hat.-Est.mat$var.sig2A.hat.)
colnames(Bias.mat) = c("Bias_mean_mu_r","Bias_mean_mu_f","Bias_mean_sig2A","Bias_var_mu","Bias_var_sig2A")

Est.mat.final = cbind(Est.mat,Bias.mat)
tmp1 = c("mu","Groups")
tmp2 = colnames(Est.mat.final[3:33])
colnames(Est.mat.final) = c(tmp1,tmp2)
head(Est.mat.final)

######################
load("Levels.rda")
shrink_1 =   matrix(,(length(n)*length(bal)*length(sigmaa2)*m[1]),8)
shrink_2 =   matrix(,(length(n)*length(bal)*length(sigmaa2)*m[2]),8)
shrink_3 =   matrix(,(length(n)*length(bal)*length(sigmaa2)*m[3]),8)

shrink = list(shrink_1,shrink_2,shrink_3)
for(j in 1:length(n))
 {
    for(k in 1:length(bal))
     {
       for(r in 1:length(sigmaa2))
       {
         for(w in 1:length(m))
         {

            a<-as.factor(levels[[w]][[j]][[k]])
            sig2A.hat = Est.mat$mean.sig2A.hat.[(j-1)*length(m)*length(bal)*length(sigmaa2)+(w-1)*length(sigmaa2)*length(bal)+(k-1)*length(sigmaa2)+r]
            sig2E.hat = Est.mat$mean.sig2E.hat.[(j-1)*length(m)*length(bal)*length(sigmaa2)+(w-1)*length(sigmaa2)*length(bal)+(k-1)*length(sigmaa2)+r]
            N = Est.mat$N[(j-1)*length(m)*length(bal)*length(sigmaa2)+(w-1)*length(sigmaa2)*length(bal)+(k-1)*length(sigmaa2)+r]
            sigma2A = Est.mat$sig2A[(j-1)*length(m)*length(bal)*length(sigmaa2)+(w-1)*length(sigmaa2)*length(bal)+(k-1)*length(sigmaa2)+r]
            sigma2 = Est.mat$sig2E[(j-1)*length(m)*length(bal)*length(sigmaa2)+(w-1)*length(sigmaa2)*length(bal)+(k-1)*length(sigmaa2)+r]
            Groups = Est.mat$m[(j-1)*length(m)*length(bal)*length(sigmaa2)+(w-1)*length(sigmaa2)*length(bal)+(k-1)*length(sigmaa2)+r]
            Type= Est.mat$Type[(j-1)*length(m)*length(bal)*length(sigmaa2)+(w-1)*length(sigmaa2)*length(bal)+(k-1)*length(sigmaa2)+r]

            for(z in 1:length(table(a)))
            {
            shrink[[w]][(j-1)*length(bal)*length(sigmaa2)*m[w]+(k-1)*length(sigmaa2)*m[w]+(r-1)*m[w]+z,8] =sig2A.hat/(sig2A.hat+(sig2E.hat/as.numeric(sort(table(a))[z])))
            shrink[[w]][(j-1)*length(bal)*length(sigmaa2)*m[w]+(k-1)*length(sigmaa2)*m[w]+(r-1)*m[w]+z,1] =N
            shrink[[w]][(j-1)*length(bal)*length(sigmaa2)*m[w]+(k-1)*length(sigmaa2)*m[w]+(r-1)*m[w]+z,2] =as.numeric(as.character(Groups))
            shrink[[w]][(j-1)*length(bal)*length(sigmaa2)*m[w]+(k-1)*length(sigmaa2)*m[w]+(r-1)*m[w]+z,3] =sigma2A
            shrink[[w]][(j-1)*length(bal)*length(sigmaa2)*m[w]+(k-1)*length(sigmaa2)*m[w]+(r-1)*m[w]+z,4] =sigma2
            shrink[[w]][(j-1)*length(bal)*length(sigmaa2)*m[w]+(k-1)*length(sigmaa2)*m[w]+(r-1)*m[w]+z,5] =Type
            shrink[[w]][(j-1)*length(bal)*length(sigmaa2)*m[w]+(k-1)*length(sigmaa2)*m[w]+(r-1)*m[w]+z,6] =as.numeric(sort(table(a)))[z]
            shrink[[w]][(j-1)*length(bal)*length(sigmaa2)*m[w]+(k-1)*length(sigmaa2)*m[w]+(r-1)*m[w]+z,7] =z
            colnames(shrink[[w]]) = c("N","Groups","Sigma2A","Sigma2","Type","n_i","Level","Shrinkage")
            }
         }
       }
   }
}

shrink.mat = as.data.frame(rbind(shrink[[1]],shrink[[2]],shrink[[3]]))

############## Power

pdf("Power_random.pdf")
power_r <- ggplot(Est.mat.final, aes(sig2A, power_r,group=Groups,shape=Groups))+ geom_line(colour="darkgrey")+geom_point() + scale_shape(solid=FALSE)+
       facet_wrap(~N+Type,ncol=3) + scale_x_continuous(breaks = seq(0,6, by = 0.5), labels = seq(0, 6, by = 0.5)) +
     opts(title = "Total sample size, Degrees of unbalancedness") +labs(x = expression(sigma[A]^2),y = "Power")+theme_bw()+
  annotate("text", x = 2.5, y = 50, label = "sigma^2==sigma[A]^2",parse=T)+ geom_vline(xintercept = 2.5,colour="grey",linetype="dashed")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(power_r)
dev.off()

pdf("Power_fixed.pdf")
power_f <- ggplot(Est.mat.final, aes(sig2A, power_f,group=Groups,shape=Groups))+ geom_line(colour="darkgrey")+geom_point() + scale_shape(solid=FALSE)+
  facet_wrap(~N+Type,ncol=3) + scale_x_continuous(breaks = seq(0,6, by = 0.5), labels = seq(0, 6, by = 0.5)) +
  opts(title = "Total sample size, Degrees of unbalancedness") +labs(x = expression(sigma[A]^2),y = "Power")+theme_bw()+
  annotate("text", x = 2.5, y = 50, label = "sigma^2==sigma[A]^2",parse=T)+ geom_vline(xintercept = 2.5,colour="grey",linetype="dashed")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(power_f)
dev.off()

############## Biases

pdf("bias_mu_random.pdf")
bias.mu.mean <- ggplot(Est.mat.final, aes(sig2A, Bias_mean_mu_r,group=Groups,shape=Groups))+ geom_line(colour="darkgrey")+geom_point() + scale_shape(solid=FALSE)+
  facet_wrap(~N+Type,ncol=3) + scale_x_continuous(breaks = seq(0,6, by = 0.5), labels = seq(0, 6, by = 0.5)) + scale_y_continuous(limits = c(-0.5,0.5))+
  opts(title = "Total sample size, Degrees of unbalancedness") +labs(x = expression(sigma[A]^2),y = expression(Bias(hat(mu))))+theme_bw()+
  annotate("text", x = 2.5, y = 0.5, label = "sigma^2==sigma[A]^2",parse=T)+ geom_vline(xintercept = 2.5,colour="grey",linetype="dashed")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(bias.mu.mean)
dev.off()

pdf("bias_mu_fixed.pdf")
bias.mu.mean <- ggplot(Est.mat.final, aes(sig2A, Bias_mean_mu_f,group=Groups,shape=Groups))+ geom_line(colour="darkgrey")+geom_point() + scale_shape(solid=FALSE)+
  facet_wrap(~N+Type,ncol=3) + scale_x_continuous(breaks = seq(0,6, by = 0.5), labels = seq(0, 6, by = 0.5)) + scale_y_continuous(limits = c(-0.5,0.5))+
  opts(title = "Total sample size, Degrees of unbalancedness") +labs(x = expression(sigma[A]^2),y = expression(Bias(hat(mu))))+theme_bw()+
  annotate("text", x = 2.5, y = 0.5, label = "sigma^2==sigma[A]^2",parse=T)+ geom_vline(xintercept = 2.5,colour="grey",linetype="dashed")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(bias.mu.mean)
dev.off()

pdf("bias_sig2A.pdf")
bias.sig2A.mean <- ggplot(Est.mat.final, aes(sig2A, Bias_mean_sig2A,group=Groups,shape=Groups))+ geom_line(colour="darkgrey")+geom_point() + scale_shape(solid=FALSE)+
  facet_wrap(~N+Type,ncol=3) + scale_x_continuous(breaks = seq(0,6, by = 0.5), labels = seq(0, 6, by = 0.5)) + 
  opts(title = "Total sample size, Degrees of unbalancedness") +labs(x = expression(sigma[A]^2),y = expression(Bias(hat(sigma[A]^2))))+theme_bw()+
  annotate("text", x = 2.5, y = 2.5, label = "sigma^2==sigma[A]^2",parse=T)+ geom_vline(xintercept = 2.5,colour="grey",linetype="dashed")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(bias.sig2A.mean)
dev.off()

pdf("bias_Var_mu.pdf")
bias.mu.var <- ggplot(Est.mat.final, aes(sig2A, Bias_var_mu,group=Groups,shape=Groups))+ geom_line(colour="darkgrey")+geom_point() + scale_shape(solid=FALSE)+
  facet_wrap(~N+Type,ncol=3) + scale_x_continuous(breaks = seq(0,6, by = 0.5), labels = seq(0, 6, by = 0.5)) + scale_y_continuous(limits = c(-0.5, 0.5))+
  opts(title = "Total sample size, Degrees of unbalancedness") +labs(x = expression(sigma[A]^2),y = expression(Bias(hat(mu))))+theme_bw()+
  annotate("text", x = 2.5, y = 0.4, label = "sigma^2==sigma[A]^2",parse=T)+ geom_vline(xintercept = 2.5,colour="grey",linetype="dashed")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(bias.mu.var)
dev.off()

pdf("bias_Var_sig2A.pdf")
bias.sig2A.var <- ggplot(Est.mat.final, aes(sig2A, Bias_var_sig2A,group=Groups,shape=Groups))+ geom_line(colour="darkgrey")+geom_point() + scale_shape(solid=FALSE)+
  facet_wrap(~N+Type,ncol=3) + scale_x_continuous(breaks = seq(0,6, by = 0.5), labels = seq(0, 6, by = 0.5)) +
  opts(title = "Total sample size, Degrees of unbalancedness") +labs(x = expression(sigma[A]^2),y = expression(Bias(hat(var(sigma[A]^2)))))+theme_bw()+
  annotate("text", x = 2.5, y = 12, label = "sigma^2==sigma[A]^2",parse=T)+ geom_vline(xintercept = 2.5,colour="grey",linetype="dashed")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(bias.sig2A.var)
dev.off()

############## Type 1 Errors
type1_long_r$Groups = as.numeric(as.character(type1_long_r$Groups))
pdf("TypeIError_random.pdf")
type1_r.plot <- ggplot(type1_long_r, aes(Groups, 0.01*(as.numeric(as.character(Type1_Error)))))+ 
  geom_boxplot() + facet_wrap(~N+Type) +   labs(y = expression(alpha),x = "No. of groups") + scale_y_continuous(limits = c(0, 0.1))+
  opts(title = "Total sample size, Degrees of unbalancedness")+scale_x_continuous(label=as.numeric(levels(type1_long_r$Groups)))
print(type1_r.plot)
dev.off()

pdf("TypeIError_fixed.pdf")
type1_f.plot <- ggplot(type1_long_f, aes(Groups, 0.01*(as.numeric(as.character(Type1_Error)))))+ 
  geom_boxplot() + facet_wrap(~N+Type) +   labs(y = expression(alpha),x = "No. of groups") + scale_y_continuous(limits = c(0, 0.1))+
  opts(title = "Total sample size, Degrees of unbalancedness")+scale_x_discrete(label=levels(type1_long_f$Groups))
print(type1_f.plot)
dev.off()
############## Shrinkage

shrink.mat$Shrinkage = as.numeric(as.character(shrink.mat$Shrinkage))
shrink.mat$Sigma2A = as.numeric(as.character(shrink.mat$Sigma2A))
shrink.mat$Level = as.numeric(as.character(shrink.mat$Level))
shrink.mat = na.omit(shrink.mat)
shrink.mat_G10 = subset(shrink.mat,Groups==10)
#shrink.plot <- ggplot(data=shrink.mat_G10,aes(Level,Shrinkage,group=Type,shape=Type))+scale_x_discrete(label=unique(shrink.mat_G10$Level))+

pdf("ShrinkagePlot.pdf")
shrink.plot <- ggplot(data=shrink.mat_G10,aes(Level,Shrinkage,group=Type,colour=Type))+scale_x_discrete(unique(shrink.mat_G10$Level))+
  geom_point() + scale_colour_manual(values=c("black", "blue","red"))+scale_shape(solid=FALSE)+ geom_line(colour="darkgrey")+
  facet_wrap(~N+Sigma2A,ncol=5) + 
  opts(title = bquote("Total sample size,"~ sigma[A]^2)) +labs(x = "Level",y="Shrinkage factor")+theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(shrink.plot)
dev.off()




