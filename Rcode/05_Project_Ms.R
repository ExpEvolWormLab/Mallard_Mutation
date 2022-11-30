rm(list = ls())
gc()

library(MCMCglmm)
library(psych)
library(gplots)
library(data.table)
library(matrixStats)
library(boot)
library(dplyr)
v_col = c("chartreuse","cadetblue1", "cornflowerblue", "slateblue2", "violetred1","darkorange","firebrick")


#Figure 3C and 5D

load("Output_files/RData/M_matrices_estimates.RData")

## Compute the projections of mmax / gmax - posterior modes
N2mmax = eigen(VCV_mat[[1]]$G1_mat/2)$vectors
PBmmax = eigen(VCV_mat[[2]]$G1_mat/2)$vectors

N2_proj_on_PBmax <- (t(PBmmax[,1])%*%(VCV_mat[[1]]$G1_mat/2)%*%PBmmax[,1])/sum(PBmmax[,1]^2)
pi_E_N2 <- N2_proj_on_PBmax/eigen(VCV_mat[[1]]$G1_mat/2)$values[1] # 47%
pi_0_N2 <- mean(eigen(VCV_mat[[1]]$G1_mat/2)$values)/eigen(VCV_mat[[1]]$G1_mat/2)$values[1] # 26%

PB_proj_on_N2max<- (t(N2mmax[,1])%*%(VCV_mat[[2]]$G1_mat/2)%*%N2mmax[,1])/sum(N2mmax[,1]^2)
pi_E_PB <- PB_proj_on_N2max/eigen(VCV_mat[[2]]$G1_mat/2)$values[1] # 61%
pi_0_PB <- mean(eigen(VCV_mat[[2]]$G1_mat/2)$values)/eigen(VCV_mat[[2]]$G1_mat/2)$values[1] # 36 %


## Compute the projections of mmax / gmax - CI and null
# Posterior distribution of Pi
vect_rand_piE_N2 <- NULL
vect_rand_piE_PB <- NULL

# Posterior distribution of Pi0
vect_rand_pi_0_N2 <- NULL
vect_rand_pi_0_PB <- NULL

n_iter=1000
spld_idx <- sample(1:nrow(VCV_mat[[1]]$VCV_Mat),n_iter)

for(i in spld_idx){

  temp_N2 <- matrix(VCV_mat[[1]]$VCV_Mat[i,1:36],6,6)
  temp_PB <- matrix(VCV_mat[[2]]$VCV_Mat[i,1:36],6,6)
  
  temp_gmax_proj <- (t(PBmmax[,1])%*%(temp_N2/2)%*%PBmmax[,1])/sum(PBmmax[,1]^2)
  vect_rand_piE_N2 <- c(vect_rand_piE_N2 ,temp_gmax_proj/eigen(temp_N2/2)$values[1])
  
  temp_gmax_proj <- (t(N2mmax[,1])%*%(temp_PB/2)%*%N2mmax[,1])/sum(N2mmax[,1]^2)
  vect_rand_piE_PB <- c(vect_rand_piE_PB , temp_gmax_proj/eigen(temp_PB/2)$values[1])	
  
  vect_rand_pi_0_N2 <- c(vect_rand_pi_0_N2,mean(eigen(temp_N2/2)$values)/eigen(temp_N2/2)$values[1])
  vect_rand_pi_0_PB <- c(vect_rand_pi_0_PB,mean(eigen(temp_PB/2)$values)/eigen(temp_PB/2)$values[1])
  
}
rm(temp_N2, temp_PB, temp_gmax_proj)


# Figure 3C

pdf('plots/Pi_adapted_Mmatrices.pdf',w=3.5)
plot(1,1,type="n",xlim=c(0,1),ylim=c(0.5,2.5),bty="n",las=1,yaxt="n",xlab=expression(Pi),ylab="",xaxt="n")
axis(1,at=c(0,0.5,1))

temp_95 <- HPDinterval(as.mcmc(vect_rand_piE_N2))
temp_80 <- HPDinterval(as.mcmc(vect_rand_piE_N2),prob=0.83)

arrows(temp_95[1,1],1.8,temp_95[1,2],1.8,code=3,angle=90,length=0.05)
arrows(temp_80[1,1],1.8,temp_80[1,2],1.8,code=3,angle=90,length=0,col=v_col[2],lwd=2)
points(mean(vect_rand_piE_N2),1.8,pch=16)

temp_95 <- HPDinterval(as.mcmc(vect_rand_pi_0_N2))
temp_80 <- HPDinterval(as.mcmc(vect_rand_pi_0_N2),prob=0.83)
arrows(temp_95[1,1],1.55,temp_95[1,2],1.55,code=3,angle=90,length=0.05)
arrows(temp_80[1,1],1.55,temp_80[1,2],1.55,code=3,angle=90,length=0,col=v_col[2],lwd=2)
points(mean(vect_rand_pi_0_N2),1.55,pch=8)


temp_95 <- HPDinterval(as.mcmc(vect_rand_piE_PB))
temp_80 <- HPDinterval(as.mcmc(vect_rand_piE_PB),prob=0.83)

arrows(temp_95[1,1],1,temp_95[1,2],1,code=3,angle=90,length=0.05)
arrows(temp_80[1,1],1,temp_80[1,2],1,code=3,angle=90,length=0,col=v_col[3],lwd=2)
points(mean(vect_rand_piE_PB),1,pch=16)

temp_95 <- HPDinterval(as.mcmc(vect_rand_pi_0_PB));temp_80 <- HPDinterval(as.mcmc(vect_rand_pi_0_PB),prob=0.83)
arrows(temp_95[1,1],.75,temp_95[1,2],.75,code=3,angle=90,length=0.05)
arrows(temp_80[1,1],.75,temp_80[1,2],.75,code=3,angle=90,length=0,col=v_col[3],lwd=2)
points(mean(vect_rand_pi_0_PB),.75,pch=8)

mtext(side=2,"N2 proj. on\n PB306mmax",padj=0,adj=0.6,cex=.9)
mtext(side=2,"PBPB306 proj.\n on N2mmax",padj=0,adj=0.1,cex=.9)

legend(0.6,2.5,c(expression(Pi),expression(Pi[0])),pch=c(16,8),lwd=1)#,bty="n")

dev.off()

### The same thing with the G

A6140_env <- new.env()
load(file='data/VCV_A6140.RData',envir=A6140_env)
VCV_mat_A6140 = A6140_env$VCV_mat_A6140
rm(A6140_env);gc()

####################################################################################
A6gmax <- eigen(VCV_mat_A6140[[1]]$G1_mat/2)$vectors

N2_proj_on_A6max <- (t(A6gmax[,1])%*%(VCV_mat[[1]]$G1_mat/2)%*%PBmmax[,1])/sum(A6gmax[,1]^2)
pi_E_N2 <- N2_proj_on_A6max/eigen(VCV_mat[[1]]$G1_mat/2)$values[1] # 39 %
pi_0_N2 <- mean(eigen(VCV_mat[[1]]$G1_mat/2)$values)/eigen(VCV_mat[[1]]$G1_mat/2)$values[1] # 26%

PB_proj_on_A6max <- (t(A6gmax[,1])%*%(VCV_mat[[2]]$G1_mat/2)%*%A6gmax[,1])/sum(A6gmax[,1]^2)
pi_E_PB <- PB_proj_on_A6max/eigen(VCV_mat[[2]]$G1_mat/2)$values[1] # 61%
pi_0_PB <- mean(eigen(VCV_mat[[2]]$G1_mat/2)$values)/eigen(VCV_mat[[2]]$G1_mat/2)$values[1] # 36 %

vect_rand_piE_N2_A6 <- NULL
vect_rand_piE_PB_A6 <- NULL

vect_rand_pi_0_N2 <- NULL
vect_rand_pi_0_PB <- NULL

n_iter=1000
spld_idx <- sample(1:nrow(VCV_mat[[1]]$VCV_Mat),n_iter)

for(i in spld_idx){
  
  temp_N2 <- matrix(VCV_mat[[1]]$VCV_Mat[i,1:36],6,6)
  temp_PB <- matrix(VCV_mat[[2]]$VCV_Mat[i,1:36],6,6)
  
  temp_gmax_proj <- (t(A6gmax[,1])%*%(temp_N2/2)%*%A6gmax[,1])/sum(A6gmax[,1]^2)
  vect_rand_piE_N2_A6 <- c(vect_rand_piE_N2_A6 ,temp_gmax_proj/eigen(temp_N2/2)$values[1])
  
  temp_gmax_proj <- (t(A6gmax[,1])%*%(temp_PB/2)%*%A6gmax[,1])/sum(A6gmax[,1]^2)
  vect_rand_piE_PB_A6 <- c(vect_rand_piE_PB_A6 ,temp_gmax_proj/eigen(temp_PB/2)$values[1])
  
  vect_rand_pi_0_N2 <- c(vect_rand_pi_0_N2,mean(eigen(temp_N2/2)$values)/eigen(temp_N2/2)$values[1])
  vect_rand_pi_0_PB <- c(vect_rand_pi_0_PB,mean(eigen(temp_PB/2)$values)/eigen(temp_PB/2)$values[1])
  
}
rm(temp_N2, temp_PB, temp_gmax_proj)

#Figure 5D
pdf('plots/Pi_adapted_Mmatrices_onA6_gmax.pdf',w=6,h=7)
plot(1,1,type="n",xlim=c(0,1),ylim=c(0.5,3),bty="n",las=1,yaxt="n",xlab=expression(Pi),ylab="",xaxt="n")
axis(1,at=c(0,0.5,1))

temp_95 <- HPDinterval(as.mcmc(vect_rand_piE_N2_A6))
temp_80 <- HPDinterval(as.mcmc(vect_rand_piE_N2_A6),prob=0.83)

arrows(temp_95[1,1],2.3,temp_95[1,2],2.3,code=3,angle=90,length=0.05)
arrows(temp_80[1,1],2.3,temp_80[1,2],2.3,code=3,angle=90,length=0,col=v_col[2],lwd=2)
points(mean(vect_rand_piE_N2_A6),2.3,pch=16)

temp_95 <- HPDinterval(as.mcmc(vect_rand_pi_0_N2))
temp_80 <- HPDinterval(as.mcmc(vect_rand_pi_0_N2),prob=0.83)
arrows(temp_95[1,1],2.05,temp_95[1,2],2.05,code=3,angle=90,length=0.05)
arrows(temp_80[1,1],2.05,temp_80[1,2],2.05,code=3,angle=90,length=0,col=v_col[2],lwd=2)
points(mean(vect_rand_pi_0_N2),2.05,pch=8)

temp_95 <- HPDinterval(as.mcmc(vect_rand_piE_PB_A6))
temp_80 <- HPDinterval(as.mcmc(vect_rand_piE_PB_A6),prob=0.83)

arrows(temp_95[1,1],1,temp_95[1,2],1,code=3,angle=90,length=0.05)
arrows(temp_80[1,1],1,temp_80[1,2],1,code=3,angle=90,length=0,col=v_col[3],lwd=2)
points(mean(vect_rand_piE_PB_A6),1,pch=16)

temp_95 <- HPDinterval(as.mcmc(vect_rand_pi_0_PB));temp_80 <- HPDinterval(as.mcmc(vect_rand_pi_0_PB),prob=0.83)
arrows(temp_95[1,1],.75,temp_95[1,2],.75,code=3,angle=90,length=0.05)
arrows(temp_80[1,1],.75,temp_80[1,2],.75,code=3,angle=90,length=0,col=v_col[3],lwd=2)
points(mean(vect_rand_pi_0_PB),.75,pch=8)

mtext(side=2,"N2 proj. on\n gmax",padj=0,adj=0.67,cex=.9)
mtext(side=2,"PB306 proj. on\n gmax",padj=0,adj=0.15,cex=.9)

legend(0.6,3,c(expression(Pi),expression(Pi[0])),pch=c(16,8),lwd=1)#,bty="n")

dev.off()



