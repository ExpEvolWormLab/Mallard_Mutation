rm(list = ls())
gc()

library(MCMCglmm)
library(psych)
library(gplots)
library(data.table)
library(matrixStats)
library(boot)
library(dplyr)

## This code computes the two M-matrices using MCMCglmm - Figure 3, Figure S1 and S2

# A function to compute the angle between two vectors x and y
angle_theta <- function(x, y) {
  dot.prod <- x %*% y
  norm.x <- norm(x, type = "2")
  norm.y <- norm(y, type = "2")
  theta <- 180/pi * as.numeric(acos(dot.prod/(norm.x * norm.y)))
  as.numeric(theta)
}

# We load the data
final_merged=read.table("data/Final_merged_data_MA_Lines.txt",sep="\t",h=TRUE)

# Trait definition
vect_P_traits = c("T12", "T13", "T21", "T23", "T31","T32")
nb_trait = length(vect_P_traits)

# Ancestral line of each MA line: either N2 or PB (short for PB306)
final_merged$anc_line=substring(final_merged$pop_label,1,2)

### Center the environmental covariates
for(i in c("temperature","rel_humidity","logD")){
           final_merged[,i] <- (final_merged[,i]-mean(final_merged[,i]))/sd(final_merged[,i])
}

# A function to plot the autocorrelation from the MCMCglmm results
plot.acfs <- function(x) {
	n <- dim(x)[2]
	par(mfrow = c(ceiling(n/2), 2), mar = c(3, 2, 3, 0.5))
	for (i in 1:n) {
		acf(x[, i], lag.max = 100, main = colnames(x)[i], ci.type = "ma", ylim = c(0, 0.3))
		grid()
	}
}

# This will be a list containing the results of the MCMCglmm
# It can be loaded from the Dryad directory
VCV_mat = NULL

# Split the data set by ancestral lines and only retain the MA lines
N2lines= subset(final_merged, anc_line =="N2" & pop_label!="N2ancL0")
PBlines= subset(final_merged, anc_line =="PB" & pop_label!="PBancL0")

k=0
## We will now estimate the two M matrix using MCMCglmm for each set of MA lines separately
for (i in c("N2", "PB306")) {

# First N2 and then PB306
	if(i=="N2") temp_df = N2lines
	if(i=="PB306") temp_df = PBlines
	
# Transform the block column from a factor to a character string format
	temp_df$data_group_name = as.character(temp_df$data_group_name)
	
	# Create the prior
	phen.var = diag(nb_trait) * diag(var(subset(temp_df, select = vect_P_traits)))
	prior_mod <- list(G = list(G1 = list(V = phen.var/3, n = nb_trait), G2 = list(V = phen.var/3, n = nb_trait)), 
		R = list(V = phen.var/3, n = nb_trait))
  # Compute the M matrix
	model_MCMC <- MCMCglmm(cbind(c(T12, T13, T21, T23, T31, T32)) ~  (temperature+rel_humidity+logD) + trait - 1, random = ~us(trait):pop_label + us(trait):data_group_name, rcov = ~us(trait):units, family = rep("gaussian", nb_trait), data = temp_df, prior = prior_mod, verbose = TRUE,nitt=500000, burnin=50000,thin=10)

	# Plot the autocorrelation results
	pdf(file = paste0("Output_files/auto_corr_plots_MCMCglmm/M_for_Model_pdf_MCMC_", 
		i, ".pdf"), height = 20)
	par(mfrow = c(10, 2), mar = c(2, 2, 1, 0))
	plot(model_MCMC$Sol, auto.layout = F)
	dev.off()

	pdf(file = paste0("Output_files/auto_corr_plots_MCMCglmm/M_for_Model_pdf_MCMC_autocorr_", 
		i, ".pdf"), height = 10)
	plot.acfs(model_MCMC$Sol)
	dev.off()

	# Store the MCMCglmm results
	post_dist = posterior.mode(model_MCMC$VCV)
	k=k+1
	VCV_mat[[k]]=list(Population = i, N_measurement = nrow(temp_df), G1_mat = matrix(post_dist[1:nb_trait^2], 
		nb_trait, nb_trait), G2_mat = matrix(post_dist[(nb_trait^2 + 1):(2 * nb_trait^2)], nb_trait, nb_trait), 
		R_mat = matrix(post_dist[(2 * nb_trait^2 + 1):(3 * nb_trait^2)], nb_trait, nb_trait), VCV_Mat = model_MCMC$VCV)

}

# Rename the M matrices
dimnames(VCV_mat[[1]]$G1_mat) = list(vect_P_traits, vect_P_traits)
dimnames(VCV_mat[[2]]$G1_mat) = list(vect_P_traits, vect_P_traits)

# Export the matrices as tables
write.table(VCV_mat[[1]]$G1_mat/2,file="Output_files/txt/VM_estimates_N2.txt",sep="\t",quote=FALSE)
write.table(VCV_mat[[2]]$G1_mat/2,file="Output_files/txt/VM_estimates_PB.txt",sep="\t",quote=FALSE)

# Get the distribution of total genetic variance in the posterior distribution
sampled_variance_N2=NULL
sampled_variance_PB=NULL
for(i in 1:nrow(VCV_mat[[1]]$VCV_Mat)){
	sampled_variance_N2=rbind(sampled_variance_N2,c(sum(eigen(matrix(VCV_mat[[1]]$VCV_Mat[i,1:36],6,6)/2)$values),eigen(matrix(VCV_mat[[1]]$VCV_Mat[i,1:36],6,6)/2)$values))
	sampled_variance_PB=rbind(sampled_variance_PB,c(sum(eigen(matrix(VCV_mat[[2]]$VCV_Mat[i,1:36],6,6)/2)$values),eigen(matrix(VCV_mat[[2]]$VCV_Mat[i,1:36],6,6)/2)$values))
}

# Output the MCMCglmm results for data archiving
save(list="VCV_mat",file="Output_files/RData/MCMC_glmm_output.RData")

# Work can be saved here
#save(list=ls(),file="Output_files/RData/M_matrices_estimates.RData")
### Restart from saved work - Produce the plots for Figure 3 and Figure S2
#rm(list=ls())
#load("Output_files/RData/M_matrices_estimates.RData")
# Vector for colors in plot

v_col = c("chartreuse","cadetblue1", "cornflowerblue", "slateblue2", "violetred1","darkorange","firebrick")

# Load the 1000 random matrices generated with data permutations
null_matrix=new.env()
load("Output_files/RData/Mmatrix_Null.RData",envir=null_matrix)
PB_null = null_matrix$PB_G1_null
N2_null = null_matrix$N2_G1_null

# Extract the null distribution for the total genetic variance
sampled_variance_N2_null=NULL
sampled_variance_PB_null=NULL

for(i in 1:1000){
  sampled_variance_N2_null=rbind(sampled_variance_N2_null,c(sum(eigen(matrix(N2_null[(i*6-5):(i*6),],6,6)/2)$values),eigen(matrix(N2_null[(i*6-5):(i*6),],6,6)/2)$values))
  sampled_variance_PB_null=rbind(sampled_variance_PB_null,c(sum(eigen(matrix(PB_null[(i*6-5):(i*6),],6,6)/2)$values),eigen(matrix(PB_null[(i*6-5):(i*6),],6,6)/2)$values))
}


# Load the A6140 matrix
A6140_env <- new.env()
load(file='data/VCV_A6140.RData',envir=A6140_env)
VCV_mat_A6140 = A6140_env$VCV_mat_A6140
rm(A6140_env);gc()


# Extract more distribution for experimental and null matrices
# Here angles between mmax of both genotype and mmax and gmax

all_gmax_angles=NULL

all_gmax_angles_PB_A6140=NULL
all_gmax_angles_N2_A6140=NULL

# Random shuffling of the posterior indexes (we will only sample 1000)
idx_N2=sample(1:nrow(VCV_mat[[1]]$VCV_Mat))
idx_PB=sample(1:nrow(VCV_mat[[1]]$VCV_Mat))
idx_A6140=sample(1:nrow(VCV_mat_A6140[[1]]$VCV_Mat))

for(i in 1:1000){
  all_gmax_angles=c(all_gmax_angles,angle_theta(
    eigen(matrix(VCV_mat[[1]]$VCV_Mat[idx_N2[i],1:36],6,6))$vectors[,1],
    eigen(matrix(VCV_mat[[2]]$VCV_Mat[idx_PB[i],1:36],6,6))$vectors[,1]
  ))
  
  all_gmax_angles_N2_A6140=c(all_gmax_angles_N2_A6140,angle_theta(
    eigen(matrix(VCV_mat[[1]]$VCV_Mat[idx_N2[i],1:36],6,6))$vectors[,1],
    eigen(matrix(VCV_mat_A6140[[1]]$VCV_Mat[idx_A6140[i],1:36],6,6))$vectors[,1]
  ))

  all_gmax_angles_PB_A6140=c(all_gmax_angles_PB_A6140,angle_theta(
    eigen(matrix(VCV_mat[[2]]$VCV_Mat[idx_PB[i],1:36],6,6))$vectors[,1],
    eigen(matrix(VCV_mat_A6140[[1]]$VCV_Mat[idx_A6140[i],1:36],6,6))$vectors[,1]
  ))
  
  }

all_gmax_angles[all_gmax_angles>90]= 180-all_gmax_angles[all_gmax_angles>90]
all_gmax_angles_N2_A6140[all_gmax_angles_N2_A6140>90]= 180-all_gmax_angles_N2_A6140[all_gmax_angles_N2_A6140>90]
all_gmax_angles_PB_A6140[all_gmax_angles_PB_A6140>90]= 180-all_gmax_angles_PB_A6140[all_gmax_angles_PB_A6140>90]



### PLOT FIGURE 3 ###

pdf(file='plots/Eigen_decomposition_with_NULL.pdf',h=7,w=8)

layout(mat=matrix(c(1:2),nrow=1), widths = c(.3,.7))
plot(c(0,0.2),c(
  sum(eigen(VCV_mat[[1]]$G1_mat)$values),
  sum(eigen(VCV_mat[[2]]$G1_mat)$values)),bg=c("black","black"),pch=21,ylab='Total genetic variance (trace)',xlim=c(-.2,.4),ylim=c(0,.75),type="n",bty="n",xlab="",xaxt="n")

axis(side=1,at=c(0,.2),labels=c("N2","PB306"),las=2)
temp_int <- HPDinterval(mcmc(sampled_variance_N2[,1]),prob=.95)
arrows(0,temp_int[1],0,temp_int[2],length=.05,code=3,angle=90)
temp_int <- HPDinterval(mcmc(sampled_variance_PB[,1]),prob=.95)
arrows(0.2,temp_int[1],0.2,temp_int[2],length=.05,code=3,angle=90)

temp_int <- HPDinterval(mcmc(sampled_variance_N2[,1]),prob=.8)
arrows(0,temp_int[1],0,temp_int[2],length=0,code=3,angle=90,lwd=2,col=v_col[2])
temp_int <- HPDinterval(mcmc(sampled_variance_PB[,1]),prob=.8)
arrows(0.2,temp_int[1],0.2,temp_int[2],length=0,code=3,angle=90,lwd=2,col=v_col[3])

points(c(0,0.2),c(
  posterior.mode(mcmc(sampled_variance_N2[,1])),
  posterior.mode(mcmc(sampled_variance_PB[,1]))),bg=c("black","black"),pch=21)

########### NULL

temp_N2 = HPDinterval(mcmc(sampled_variance_N2_null[,1]),prob=.95)
arrows(0+.06,temp_N2[1],0+.06,temp_N2[2],length=.05,code=3,angle=90,col="grey")

temp_PB = HPDinterval(mcmc(sampled_variance_PB_null[,1]),prob=.95)
arrows(.2+.06,temp_PB[1],.2+.06,temp_PB[2],length=.05,code=3,angle=90,col="grey")


#############

plot(eigen(VCV_mat[[1]]$G1_mat)$values~c(1:6),pch=16,xlim=c(.5,6.5),bty="n",ylim=c(0,.5),ylab="Genetic variance",type="n",xlab="Eigen vectors",xaxt="n")
axis(1,at=1:6,c(expression(m[max]),expression(m[2]),expression(m[3]),expression(m[4]),expression(m[5]),expression(m[6])))
for(i in 1:6){
  temp_int <- HPDinterval(mcmc(sampled_variance_N2[,(i+1)]),prob=.95)
  arrows(i,temp_int[1],i,temp_int[2],length=.05,code=3,angle=90)
  temp_int <- HPDinterval(mcmc(sampled_variance_PB[,(i+1)]),prob=.95)
  arrows(i+0.2,temp_int[1],i+0.2,temp_int[2],length=.05,code=3,angle=90)
  
  temp_int <- HPDinterval(mcmc(sampled_variance_N2[,(i+1)]),prob=.8)
  arrows(i,temp_int[1],i,temp_int[2],length=0,code=3,angle=90,lwd=2,col=v_col[2])
  temp_int <- HPDinterval(mcmc(sampled_variance_PB[,(i+1)]),prob=.8)
  arrows(i+0.2,temp_int[1],i+0.2,temp_int[2],length=0,code=3,angle=90,lwd=2,col=v_col[3])
  
}

points( (c(1:6)),posterior.mode(mcmc(sampled_variance_N2[,2:7])) ,pch=16,col="black")
points( (c(1:6)+0.2),posterior.mode(mcmc(sampled_variance_PB[,2:7])) ,pch=16,col="black")

legend(4,1,c("N2","PB306"),lwd=2,col=v_col[2:3],bty="n")
legend(3.5,.8,bty="n",c("Null\ndistribution"),col='grey',lwd=1)

########### NULL

temp_N2 = HPDinterval(mcmc(sampled_variance_N2_null[,2:7]),prob=.95)
arrows((1:6)+.06,temp_N2[1],(1:6)+.06,temp_N2[2],length=.05,code=3,angle=90,col="grey")

temp_PB = HPDinterval(mcmc(sampled_variance_PB_null[,2:7]),prob=.95)
arrows((1:6)+.2+.06,temp_PB[1],(1:6)+.2+.06,temp_PB[2],length=.05,code=3,angle=90,col="grey")

dev.off()


### We need to compute a null expectation for the angle between two eigenvectors
### of dimension 6

random_angles=NULL
for(i in 1:10000){
  v1 = runif(6,min=c(-1),max=c(1))
  v2 = runif(6,min=c(-1),max=c(1))
  random_angles <- c(random_angles,min(angle_theta(v1,v2),angle_theta(-v1,v2),angle_theta(v1,-v2),angle_theta(-v1,-v2)))
}
hist(random_angles)

### Figure S1
pdf(file='plots/Angles_between_mmax.pdf',w=4)
plot(all_gmax_angles~rep(1,length(all_gmax_angles)),type="n",las=1,bty="n",xlab="",xaxt="n",ylab=
       expression(paste("Angles between ",m[max])),yaxt="s",xlim=c(0.7,1.5),ylim=c(0,90))

temp95 = HPDinterval(as.mcmc(all_gmax_angles))

arrows(1,temp95[,1 ],1,temp95[,2],code=3,angle=90,length=.1)

points(1,180-angle_theta(
  eigen(VCV_mat[[1]]$G1_mat)$vectors[,1],
  eigen(VCV_mat[[2]]$G1_mat)$vectors[,1]),pch=16)

### Random expectation
temp95 = HPDinterval(as.mcmc(random_angles))
arrows(1.2,temp95[,1 ],1.2,temp95[,2],code=3,angle=90,length=.1)
points(1.2,mean((random_angles)),pch=8)
legend(0.8,15,pch=c(16,8),lwd=1,c("Experimental data","Null expectation"),bty="n")

dev.off()


#### Figure S2

### Angles with A6140
pdf(file='plots/A6140_gmax_with_mmax.pdf',w=4)
plot(all_gmax_angles~rep(1,length(all_gmax_angles)),type="n",las=1,bty="n",xlab="",xaxt="n",ylab=
       expression(paste("Angles between ",g[max]," and ",m[max])),yaxt="s",xlim=c(0.5,2.5),ylim=c(0,90))

axis(side=1,at=1:2,labels=c("N2","PB306"))
temp95 = HPDinterval(as.mcmc(all_gmax_angles_N2_A6140))
temp80 = HPDinterval(as.mcmc(all_gmax_angles_N2_A6140),prob=.8)
arrows(1,temp95[,1 ],1,temp95[,2],code=3,angle=90,length=.1)
arrows(1,temp80[,1],1,temp80,length=0,col=v_col[2],lwd=1.5)
points(1,180-angle_theta(
  eigen(VCV_mat[[1]]$G1_mat)$vectors[,1],
  eigen(VCV_mat_A6140[[1]]$G1_mat)$vectors[,1]),pch=16)


temp95 = HPDinterval(as.mcmc(all_gmax_angles_PB_A6140))
temp80 = HPDinterval(as.mcmc(all_gmax_angles_PB_A6140),prob=.8)
arrows(2,temp95[,1],2,temp95[,2],code=3,angle=90,length=.1)
arrows(2,temp80[,1],2,temp80,length=0,col=v_col[3],lwd=1.5)

points(2,angle_theta(
  eigen(VCV_mat[[2]]$G1_mat)$vectors[,1],
  eigen(VCV_mat_A6140[[1]]$G1_mat)$vectors[,1]),pch=16)

### Random expectation
temp95 = HPDinterval(as.mcmc(random_angles))
arrows(1.5,temp95[,1 ],1.5,temp95[,2],code=3,angle=90,length=.1)
points(1.5,mean((random_angles)),pch=8)
legend(0.8,15,pch=c(16,8),lwd=1,c("Experimental data","Null expectation"),bty="n")


dev.off()
  

