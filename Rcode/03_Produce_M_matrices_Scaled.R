rm(list = ls())
gc()

library(MCMCglmm)
library(psych)
library(gplots)
library(data.table)
library(matrixStats)
library(boot)
library(dplyr)

# Compute the scaled G and M matrices - plot for Figure 4

# No prior correction
final_merged=read.table("data/Final_merged_data_MA_Lines.txt",sep="\t",h=TRUE)
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

VCV_mat = NULL

# Split the data set by ancestral lines and only retain the MA lines
N2lines= subset(final_merged, anc_line =="N2" & pop_label!="N2ancL0")
PBlines= subset(final_merged, anc_line =="PB" & pop_label!="PBancL0")

## Scale with the mean P variance

meanSd_N2 <- mean(colSds(as.matrix(N2lines[,vect_P_traits])))
meanSd_PB <- mean(colSds(as.matrix(PBlines[,vect_P_traits])))

for(i in 1:6){
  N2lines[,vect_P_traits[i]] <- (N2lines[,vect_P_traits[i]]-mean(N2lines[,vect_P_traits[i]]))/meanSd_N2
  PBlines[,vect_P_traits[i]] <- (PBlines[,vect_P_traits[i]]-mean(PBlines[,vect_P_traits[i]]))/meanSd_PB
}
plot(colSds(as.matrix(N2lines[,vect_P_traits])),
colSds(as.matrix(PBlines[,vect_P_traits])),asp=1)
abline(a=0,b=1)


k=0
## We will now estimate the two M matrix using MCMCglmm for each set of MA lines separately
for (i in c("N2", "PB")) {

# First N2 and then PB306
	if(i=="N2") temp_df = N2lines
	if(i=="PB") temp_df = PBlines

# Transform the block column from a factor to a character string format
	temp_df$data_group_name = as.character(temp_df$data_group_name)
	
	phen.var = diag(nb_trait) * diag(var(subset(temp_df, select = vect_P_traits)))

	prior_mod <- list(G = list(G1 = list(V = phen.var/3, n = nb_trait), G2 = list(V = phen.var/3, n = nb_trait)), 
		R = list(V = phen.var/3, n = nb_trait))

	model_MCMC <- MCMCglmm(cbind(c(T12, T13, T21, T23, T31, T32)) ~  (temperature+rel_humidity+logD) + trait - 1, random = ~us(trait):pop_label + us(trait):data_group_name, rcov = ~us(trait):units, family = rep("gaussian", nb_trait), data = temp_df, prior = prior_mod, verbose = TRUE,nitt=500000, burnin=50000,thin=10)

	pdf(file = paste0("Output_files/auto_corr_plots_MCMCglmm_scaled/M_for_Model_pdf_MCMC_", 
		i, ".pdf"), height = 20)
	par(mfrow = c(10, 2), mar = c(2, 2, 1, 0))
	plot(model_MCMC$Sol, auto.layout = F)
	dev.off()

	pdf(file = paste0("Output_files/auto_corr_plots_MCMCglmm_scaled/M_for_Model_pdf_MCMC_autocorr_", 
		i, ".pdf"), height = 10)
	plot.acfs(model_MCMC$Sol)
	dev.off()

	post_dist = posterior.mode(model_MCMC$VCV)
	k=k+1
	VCV_mat[[k]]=list(Population = i, N_measurement = nrow(temp_df), G1_mat = matrix(post_dist[1:nb_trait^2], 
		nb_trait, nb_trait), G2_mat = matrix(post_dist[(nb_trait^2 + 1):(2 * nb_trait^2)], nb_trait, nb_trait), 
		R_mat = matrix(post_dist[(2 * nb_trait^2 + 1):(3 * nb_trait^2)], nb_trait, nb_trait), VCV_Mat = model_MCMC$VCV)

}

dimnames(VCV_mat[[1]]$G1_mat) = list(vect_P_traits, vect_P_traits)
dimnames(VCV_mat[[2]]$G1_mat) = list(vect_P_traits, vect_P_traits)

write.table(VCV_mat[[1]]$G1_mat/2,file="Output_files/txt/VM_estimates_N2_scaled.txt",sep="\t",quote=FALSE)
write.table(VCV_mat[[2]]$G1_mat/2,file="Output_files/txt/VM_estimates_PB_scaled.txt",sep="\t",quote=FALSE)

sampled_variance_N2=NULL
sampled_variance_PB=NULL
for(i in 1:nrow(VCV_mat[[1]]$VCV_Mat)){
	sampled_variance_N2=rbind(sampled_variance_N2,c(sum(eigen(matrix(VCV_mat[[1]]$VCV_Mat[i,1:36],6,6))$values),eigen(matrix(VCV_mat[[1]]$VCV_Mat[i,1:36],6,6))$values))
	sampled_variance_PB=rbind(sampled_variance_PB,c(sum(eigen(matrix(VCV_mat[[2]]$VCV_Mat[i,1:36],6,6))$values),eigen(matrix(VCV_mat[[2]]$VCV_Mat[i,1:36],6,6))$values))
}


#########################################
######### Comparison with A6140 #########
#########################################

# Load the A6140 matrix

A6140_env <- new.env()
load(file='data/VCV_A6140.RData',envir=A6140_env)

final_A6140 <-  subset(A6140_env$final_merged,population=='A6140' & location_label=='Lisbon')
final_CA50 <-  subset(A6140_env$final_merged,population%in%paste0("CA",c(1:3),"50") & location_label=='Paris')
final_CA100 <-  subset(A6140_env$final_merged,population%in%paste0("CA",c(1:3),"100"))

meanSd_A6140 <- mean(colSds(as.matrix(final_A6140[,vect_P_traits])))
meanSd_CA50 <- mean(colSds(as.matrix(final_CA50[,vect_P_traits])))
meanSd_CA100 <- mean(colSds(as.matrix(final_CA100[,vect_P_traits])))
for(i in 1:6){
  final_A6140[,vect_P_traits[i]] <- (final_A6140[,vect_P_traits[i]]-mean(final_A6140[,vect_P_traits[i]]))/meanSd_A6140
  final_CA50[,vect_P_traits[i]] <- (final_CA50[,vect_P_traits[i]]-mean(final_CA50[,vect_P_traits[i]]))/meanSd_CA50
  final_CA100[,vect_P_traits[i]] <- (final_CA100[,vect_P_traits[i]]-mean(final_CA100[,vect_P_traits[i]]))/meanSd_CA100  
}

phen.var = diag(nb_trait) * diag(var(subset(final_A6140, select = vect_P_traits)))
prior_mod <- list(G = list(G1 = list(V = phen.var/3, n = nb_trait), G2 = list(V = phen.var/3, n = nb_trait)), 
                  R = list(V = phen.var/3, n = nb_trait))

model_MCMC <- MCMCglmm(cbind(c(T12, T13, T21, T23, T31, T32)) ~  (temperature+rel_humidity+logD) + is_2012 + trait - 1, random = ~us(trait):pop_label + us(trait):date_str, rcov = ~us(trait):units, family = rep("gaussian", nb_trait), data = final_A6140, prior = prior_mod, verbose = TRUE,nitt=500000, burnin=50000,thin=10)

i="A6140"
pdf(file = paste0("Output_files/auto_corr_plots_MCMCglmm_Scaled/Model_pdf_MCMC_", 
                  i, ".pdf"), height = 20)
par(mfrow = c(10, 2), mar = c(2, 2, 1, 0))
plot(model_MCMC$Sol, auto.layout = F)
dev.off()

pdf(file = paste0("Output_files/auto_corr_plots_MCMCglmm_Scaled/Model_pdf_MCMC_autocorr_", 
                  i, ".pdf"), height = 10)
plot.acfs(model_MCMC$Sol)
dev.off()


post_dist = posterior.mode(model_MCMC$VCV)
VCV_mat_A6140=list()
VCV_mat_A6140[[1]]=list(Population = i, N_measurement = nrow(temp_df), G1_mat = matrix(post_dist[1:nb_trait^2], 
                                                                                 nb_trait, nb_trait), G2_mat = matrix(post_dist[(nb_trait^2 + 1):(2 * nb_trait^2)], nb_trait, nb_trait), 
                  R_mat = matrix(post_dist[(2 * nb_trait^2 + 1):(3 * nb_trait^2)], nb_trait, nb_trait), VCV_Mat = model_MCMC$VCV)

rm(A6140_env);gc()

# Work can be saved here
#save(list=ls(),file="Output_files/RData/M_matrices_estimates_scaled.RData")
#rm(list=ls());gc()
#load("Output_files/RData/M_matrices_estimates_scaled.RData")

# Output the MCMCglmm lists for Dryad
save(list=c("VCV_mat_A6140","VCV_mat"),file="Output_files/RData/MCMC_glmm_output_SCALED.RData")

v_col = c("chartreuse","cadetblue1", "cornflowerblue", "slateblue2", "violetred1","darkorange","firebrick")

#  Figure 4

pdf(file='plots/M_matrices_with_A6140_scaled.pdf',h=8,w=6)

par(mar=c(5,7,4,2))
vect_Var <- c(2:6,9:12,16:18,23,24,30,1,8,15,22,29,36)
vProb <- .95

plot(c(VCV_mat_A6140[[1]]$G1_mat/2)[vect_Var],c(24:10,6:1),yaxt="n",bty="n",xlim=c(-.3,1),xlab="Genetic (co)variances",xaxt="n",type='n',ylab="",cex.lab=1.2)


lines(c(0,0),c(24.5,8.5))
lines(c(0,0),c(.5,5.5),col="red")

axis(side=1,pos=0)

axis(side=2,at=c(24:10,6:1),labels=c("SF*SB","SF*FS","SF*FB","SF*BS","SF*BF",
                                     "SB*FS","SB*FB","SB*BS","SB*BF","FS*FB",
                                     "FS*BS","FS*BF","FB*BS","FB*BF","BS*BF",
                                     "SF","SB","FS","FB","BS","BF"),las=1)

i=1
temp_95 <- HPDinterval(VCV_mat_A6140[[i]]$VCV_Mat[,1:36]/2,prob=.95)
temp_80 <- HPDinterval(VCV_mat_A6140[[i]]$VCV_Mat[,1:36]/2,prob=.83)

arrows(temp_95[vect_Var,1],c(24:10,6:1)+(.25*(i-1)),temp_95[vect_Var,2],c(24:10,6:1)+(.25*(i-1)),code=3,length=.02,angle=90)
arrows(temp_80[vect_Var,1],c(24:10,6:1)+(.25*(i-1)),
       temp_80[vect_Var,2],c(24:10,6:1)+(.25*(i-1)),code=3,length=0,angle=90,lwd=2,col="chartreuse")

points(c(VCV_mat_A6140[[i]]$G1_mat/2)[vect_Var],c(24:10,6:1)+(.25*(i-1)),pch=21,bg="black",cex=.6)


for(i in 1:2){
  
  temp_95 <- HPDinterval(VCV_mat[[i]]$VCV_Mat[,1:36]/2,prob=.95)
  
  temp_80 <- HPDinterval(VCV_mat[[i]]$VCV_Mat[,1:36]/2,prob=.83)
  
  arrows(temp_95[vect_Var,1],c(24:10,6:1)+.3+(.3*(i-1)),temp_95[vect_Var,2],c(24:10,6:1)+.3+(.3*(i-1)),code=3,length=.02,angle=90)
  arrows(temp_80[vect_Var,1],c(24:10,6:1)+.3+(.3*(i-1)),
         temp_80[vect_Var,2],c(24:10,6:1)+.3+(.3*(i-1)),code=3,length=0,angle=90,lwd=2,col=v_col[i+1])
  
  points(c(VCV_mat[[i]]$G1_mat/2)[vect_Var],c(24:10,6:1)+.3+(.3*(i-1)),pch=21,bg="black",cex=.6)
  
}

legend(.35,23,c("A6140","N2","PB306"),ncol=1,v_col,bty="n")

dev.off()
