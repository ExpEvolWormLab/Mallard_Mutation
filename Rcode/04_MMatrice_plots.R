rm(list=ls())
library(data.table)
library(lme4)
library(plyr)
library(dplyr)
library(ggpubr)
library(RColorBrewer)
library(matrixStats)
library(boot) 
library(MCMCglmm)

#### Code that produce the plots for Figure 2 and some exports for tables

load("Output_files/RData/MCMC_glmm_output.RData")
v_col = c("chartreuse","cadetblue1", "cornflowerblue", "slateblue2", "violetred1","darkorange","firebrick")

null_matrix=new.env()
load("Output_files/RData/Mmatrix_Null.RData",envir=null_matrix)

PB_null = null_matrix$PB_G1_null
N2_null = null_matrix$N2_G1_null

rm(null_matrix);gc()

## We need to reformat the NULL
PB_null_v2 = NULL
N2_null_v2 = NULL
for(i in 1:1000){
  
  N2_null_v2=rbind(N2_null_v2,c(N2_null[(i*6-5):(i*6),]/2))
  PB_null_v2=rbind(PB_null_v2,c(PB_null[(i*6-5):(i*6),]/2))
  
}


pdf("plots/Genetic_Variances.pdf",w=5)

temp_95_N2 <- HPDinterval(VCV_mat[[1]]$VCV_Mat[,1:36]/2,prob=.95)
temp_95_PB <- HPDinterval(VCV_mat[[2]]$VCV_Mat[,1:36]/2,prob=.95)
temp_95_NULL_N2 <- HPDinterval(as.mcmc(N2_null_v2),prob=.95)
temp_95_NULL_PB <- HPDinterval(as.mcmc(PB_null_v2),prob=.95)

  
vect_CoVar <- c(2:6,9:12,16:18,23,24,30)
vect_Diag <- c(1,8,15,22,29,36)

plot(1:6,c(VCV_mat[[1]]$G1_mat/2)[vect_Diag],pch=21,col="black",las=1,bty="n",ylab=c("Genetic Variance"),bg=v_col[2],ylim=c(0,.18),type="n",xlim=c(0.7,7),xaxt="n",xlab="Transition rates")
axis(side=1,at=1:6,labels=c("SF","SB","FS","FB","BS","BF"))

arrows(1:6-.15,temp_95_NULL_N2[vect_Diag,1],1:6-.15,temp_95_NULL_N2[vect_Diag,2],code=3,length=.05,col="orange",angle=90,lwd=2)
points(1:6-.15,c(VCV_mat[[1]]$G1_mat/2)[vect_Diag],pch=21,col="black",bg=v_col[2])

arrows(1:6+.15,temp_95_NULL_PB[vect_Diag,1],1:6+.15,temp_95_NULL_PB[vect_Diag,2],code=3,length=.05,col="orange",angle=90,lwd=2)
points(1:6+.15,c(VCV_mat[[2]]$G1_mat/2)[vect_Diag],pch=21,col="black",bg=v_col[3])

points(c(1:6)-.15,rep(0,6),pch=8)
points(c(1,3,4)+.15,rep(0,3),pch=8)

legend(.9,.16,"Posterior mean",bty="n")
legend(1.15,.152,pch=16,c("N2","PB306"),bty="n",ncol=2,col=v_col[2:3])
legend(.9,.14,lwd=2,"95% CI of null posterior\nmeans",bty="n",col="orange")
dev.off()

pdf("plots/Genetic_CoVariances.pdf",w=7)
par(mar=c(8,4,4,2))
plot(1:15,c(VCV_mat[[1]]$G1_mat/2)[vect_CoVar],pch=21,col="black",las=1,bty="n",ylab=c("Genetic covariances"),bg=v_col[2],ylim=c(-.12,.18),xlim=c(0.7,16),xaxt="n",xlab="",type="n")
axis(side=1,at=1:15,labels=c("SF*SB","SF*FS","SF*FB","SF*BS","SF*BF",
                             "SB*FS","SB*FB","SB*BS","SB*BF","FS*FB",
                             "FS*BS","FS*BF","FB*BS","FB*BF","BS*BF"),las=2)
abline(h=0,lty=2)
arrows(1:15-.2,temp_95_N2[vect_CoVar,1],1:15-.2,temp_95_N2[vect_CoVar,2],code=3,length=.05,col="black",angle=90)
arrows(1:15+.2,temp_95_PB[vect_CoVar,1],1:15+.2,temp_95_PB[vect_CoVar,2],code=3,length=.05,col="black",angle=90)

points(1:15+.2,c(VCV_mat[[2]]$G1_mat/2)[vect_CoVar],bg=v_col[3],pch=21)
points(1:15-.2,c(VCV_mat[[1]]$G1_mat/2)[vect_CoVar],bg=v_col[2],pch=21)

points(c(1,2,8,11,14,15,11)+c(rep(-.2,6),.2),rep(-.1,7),pch=8)

dev.off()


### Export tables - Null

t1=round(temp_95_NULL_N2[vect_Diag,],4)
as.data.frame(paste0("[",t1[,1]," - ",t1[,2],"]"))
t1=round(temp_95_NULL_PB[vect_Diag,],4)
as.data.frame(paste0("[",t1[,1]," - ",t1[,2],"]"))

#### Export the eigen decomposition of the M-matrices

# Mode
data.frame(Ev=round(eigen(VCV_mat[[1]]$G1_mat)$values,4))
data.frame(Ev=round(eigen(VCV_mat[[2]]$G1_mat)$values,4))

data.frame(Ev=round(eigen(VCV_mat[[1]]$G1_mat)$vectors,4))

# Proportion
data.frame(cbind(round(eigen(VCV_mat[[1]]$G1_mat)$values/sum(eigen(VCV_mat[[1]]$G1_mat)$values),4),round(eigen(VCV_mat[[2]]$G1_mat)$values/sum(eigen(VCV_mat[[2]]$G1_mat)$values),4)))


all_EV_N2 = NULL
all_EV_PB = NULL

for(i in 1:nrow(VCV_mat[[1]]$VCV_Mat)){
  all_EV_N2=rbind(all_EV_N2,eigen(matrix(VCV_mat[[1]]$VCV_Mat[i,1:36],6,6))$values)
  all_EV_PB=rbind(all_EV_PB,eigen(matrix(VCV_mat[[2]]$VCV_Mat[i,1:36],6,6))$values)
}
t(round(apply(all_EV_N2,2,function(x){HPDinterval(as.mcmc(x))}),4))
t(round(apply(all_EV_PB,2,function(x){HPDinterval(as.mcmc(x))}),4))

#### Export the eigen decomposition of the M-matrices

# Mode
data.frame(Ev=round(eigen(VCV_mat[[1]]$G1_mat)$values,4))
data.frame(Ev=round(eigen(VCV_mat[[2]]$G1_mat)$values,4))

data.frame(Ev=round(eigen(VCV_mat[[1]]$G1_mat)$vectors,4))

# Proportion
data.frame(cbind(round(eigen(VCV_mat[[1]]$G1_mat)$values/sum(eigen(VCV_mat[[1]]$G1_mat)$values),4),round(eigen(VCV_mat[[2]]$G1_mat)$values/sum(eigen(VCV_mat[[2]]$G1_mat)$values),4)))


all_EV_N2 = NULL
all_EV_PB = NULL

for(i in 1:nrow(VCV_mat[[1]]$VCV_Mat)){
  all_EV_N2=rbind(all_EV_N2,eigen(matrix(VCV_mat[[1]]$VCV_Mat[i,1:36],6,6))$values)
  all_EV_PB=rbind(all_EV_PB,eigen(matrix(VCV_mat[[2]]$VCV_Mat[i,1:36],6,6))$values)
}
t(round(apply(all_EV_N2,2,function(x){HPDinterval(as.mcmc(x))}),4))
t(round(apply(all_EV_PB,2,function(x){HPDinterval(as.mcmc(x))}),4))


### The same tables with some important scaled vectors

load("Output_files/RData/MCMC_glmm_output_SCALED.RData")
VCV_mat_with_MA=VCV_mat_A6140
VCV_mat_with_MA[2]=VCV_mat[1]
VCV_mat_with_MA[3]=VCV_mat[2]
rm(VCV_mat_A6140);rm(VCV_mat);gc()

data.frame(A6140=round(eigen(VCV_mat_with_MA[[1]]$G1_mat)$values,4),
           N2=round(eigen(VCV_mat_with_MA[[2]]$G1_mat)$values,4),
           PB306=round(eigen(VCV_mat_with_MA[[3]]$G1_mat)$values,4))

# Proportion
data.frame(cbind(
  round(eigen(VCV_mat_with_MA[[1]]$G1_mat)$values/sum(eigen(VCV_mat_with_MA[[1]]$G1_mat)$values),4),
  round(eigen(VCV_mat_with_MA[[2]]$G1_mat)$values/sum(eigen(VCV_mat_with_MA[[2]]$G1_mat)$values),4),
  round(eigen(VCV_mat_with_MA[[3]]$G1_mat)$values/sum(eigen(VCV_mat_with_MA[[3]]$G1_mat)$values),4)))


all_EV_A6140 = NULL
for(i in 1:nrow(VCV_mat_with_MA[[1]]$VCV_Mat)){
  all_EV_A6140=rbind(all_EV_A6140,eigen(matrix(VCV_mat_with_MA[[1]]$VCV_Mat[i,1:36],6,6))$values)
  
}
t(round(apply(all_EV_A6140,2,function(x){HPDinterval(as.mcmc(x))}),4))

# Loadings
data.frame(Ev=round(eigen(VCV_mat_with_MA[[1]]$G1_mat)$vectors,4))




