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

load("Output_files/RData/Mmatrix_Null.RData")
load("Output_files/RData/M_matrices_estimates.RData")
v_col = c("chartreuse","cadetblue1", "cornflowerblue", "slateblue2", "violetred1","darkorange","firebrick")

###Reformat N2_G1_null/PB_G1_null
df_G1_N2=NULL
df_G1_PB=NULL
for(i in 1:1000){
 idx = (6*(i-1)+1) : (6*(i))
df_G1_N2 = rbind(df_G1_N2, as.numeric(N2_G1_null[idx,]))
df_G1_PB = rbind(df_G1_PB, as.numeric(PB_G1_null[idx,]))
}

pdf(file='plots/M_mat_with_NULLs.pdf',h=8,w=5.5)

true_G1=VCV_mat[[1]]
true_G2=VCV_mat[[2]]
df_G1=df_G1_N2
df_G2=df_G1_PB

    par(mar=c(5,7,4,2))
  vect_Var <- c(2:6,9:12,16:18,23,24,30,1,8,15,22,29,36)
  vProb <- .95
  
  
  vectX <- c(true_G1$G1_mat/2)[vect_Var]
  vectX2 <- c(true_G2$G1_mat/2)[vect_Var]

  plot(vectX,c(24:10,6:1),yaxt="n",bty="n",xlim=c(-.1,.30),xlab="Genetic (co)variances",xaxt="n",type='n',ylab="",cex.lab=1.2)
  lines(c(0,0),c(24.5,8.5))
  lines(c(0,0),c(.5,5.5),col="red")
  
  axis(side=1,pos=0)
  axis(side=2,at=c(24:10,6:1),labels=c("SF*SB","SF*FS","SF*FB","SF*BS","SF*BF",
                                       "SB*FS","SB*FB","SB*BS","SB*BF","FS*FB",
                                       "FS*BS","FS*BF","FB*BS","FB*BF","BS*BF",
                                       "SF","SB","FS","FB","BS","BF"),las=1)
  
  #### G1
  temp_95 <- HPDinterval(true_G1$VCV_Mat[,1:36]/2,prob=.95)
  temp_80 <- HPDinterval(true_G1$VCV_Mat[,1:36]/2,prob=.8)
  
  arrows(temp_95[vect_Var,1],c(24:10,6:1)+.3,temp_95[vect_Var,2],c(24:10,6:1)+.3,code=3,length=.02,angle=90)
  arrows(temp_80[vect_Var,1],c(24:10,6:1)+.3,
         temp_80[vect_Var,2],c(24:10,6:1)+.3,code=3,length=0,angle=90,lwd=2,col=v_col[2])
  
  points(vectX,c(24:10,6:1)+.3,pch=21,bg="black",cex=.6)
  
  #### G2
  temp_95 <- HPDinterval(true_G2$VCV_Mat[,1:36]/2,prob=.95)
  temp_80 <- HPDinterval(true_G2$VCV_Mat[,1:36]/2,prob=.8)

  arrows(temp_95[vect_Var,1],c(24:10,6:1)+.6,temp_95[vect_Var,2],c(24:10,6:1)+.6,code=3,length=.02,angle=90)
  arrows(temp_80[vect_Var,1],c(24:10,6:1)+.6,
         temp_80[vect_Var,2],c(24:10,6:1)+.6,code=3,length=0,angle=90,lwd=2,col=v_col[3])
  
  points(vectX2,c(24:10,6:1)+.6,pch=21,bg="black",cex=.6)
  
  ##Random 1  
  temp_95_rand <- HPDinterval(as.mcmc(df_G1[,1:36]/2),prob=.95)
  i=1;arrows(temp_95_rand[vect_Var,1],c(24:10,6:1)+.3+(.3*(i-1)) - .1,temp_95_rand[vect_Var,2],c(24:10,6:1)+.3+(.3*(i-1))- .1,code=3,length=.02,angle=90,col="orange",lwd=1.3)
  
  ##Random 2
  temp_95_rand <- HPDinterval(as.mcmc(df_G2[,1:36]/2),prob=.95)
  i=2;arrows(temp_95_rand[vect_Var,1],c(24:10,6:1)+.3+(.3*(i-1)) - .1,temp_95_rand[vect_Var,2],c(24:10,6:1)+.3+(.3*(i-1))- .1,code=3,length=.02,angle=90,col="orange",lwd=1.3)  
  legend(.1,25,c("N2","PB306","Null distribution"),lwd=2,col=c(v_col[2:3],"orange"),bty="n")
  

dev.off()





