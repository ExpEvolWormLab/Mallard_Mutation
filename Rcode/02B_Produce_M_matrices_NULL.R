rm(list = ls())
gc()

library(MCMCglmm)
library(psych)
library(gplots)
library(data.table)
library(matrixStats)
library(boot)
library(dplyr)

# No prior correction
final_merged=read.table("data/Final_merged_data_MA_Lines.txt",sep="\t",h=TRUE)
head(final_merged)

vect_P_traits = c("T12", "T13", "T21", "T23", "T31","T32")
nb_trait = length(vect_P_traits)

# Ancestral line of each MA line: either N2 or PB (short for PB306)
final_merged$anc_line=substring(final_merged$pop_label,1,2)


### Center the environmental covariates
for(i in c("temperature","rel_humidity","logD")){
  final_merged[,i] <- (final_merged[,i]-mean(final_merged[,i]))/sd(final_merged[,i])
}

N2_G1_null = NULL; N2_G2_null = NULL; N2_R_null  = NULL
PB_G1_null = NULL; PB_G2_null = NULL; PB_R_null  = NULL

# Split the data set by ancestral lines and only retain the MA lines
N2lines= subset(final_merged, anc_line =="N2" & pop_label!="N2ancL0")
PBlines= subset(final_merged, anc_line =="PB" & pop_label!="PBancL0")

## We will now estimate the two M matrix using MCMCglmm for each set of MA lines separately

for(iii in 1:1000){
  print(iii)
  for (i in c("N2", "PB")) {
    
    # First N2 and then PB306
    if(i=="N2") temp_df = N2lines
    if(i=="PB") temp_df = PBlines
    
    # Transform the block column from a factor to a character string format
    temp_df$data_group_name = as.character(temp_df$data_group_name)
    phen.var = diag(nb_trait) * diag(var(subset(temp_df, select = vect_P_traits)))
    
    prior_mod <- list(G = list(G1 = list(V = phen.var/3, n = nb_trait), G2 = list(V = phen.var/3, n = nb_trait)), 
                      R = list(V = phen.var/3, n = nb_trait))
    
    
    ## Randomize
    temp_df$data_group_name = sample(temp_df$data_group_name)
    temp_df$pop_label = sample(temp_df$pop_label)
    
    model_MCMC <- MCMCglmm(cbind(c(T12, T13, T21, T23, T31, T32)) ~  (temperature+rel_humidity+logD) + trait - 1,
                           random = ~us(trait):pop_label + us(trait):data_group_name, rcov = ~us(trait):units,
                           family = rep("gaussian", nb_trait), data = temp_df, prior = prior_mod, verbose = FALSE,nitt=50000,
                           burnin=20000,thin=10)
    
    
    post_dist = posterior.mode(model_MCMC$VCV)
    
    tempG1_post =	matrix(post_dist[1:nb_trait^2],nb_trait, nb_trait)
    tempG2_post =	matrix(post_dist[(nb_trait^2 + 1):(2 * nb_trait^2)], nb_trait, nb_trait)
    tempR_post =	matrix(post_dist[(2 * nb_trait^2 + 1):(3 * nb_trait^2)], nb_trait, nb_trait)
    if(i == "N2"){
      N2_G1_null = rbind(N2_G1_null,tempG1_post)
      N2_G2_null = rbind(N2_G2_null,tempG2_post)
      N2_R_null = rbind(N2_R_null,tempR_post)
    }else{
      PB_G1_null = rbind(PB_G1_null,tempG1_post)
      PB_G2_null = rbind(PB_G2_null,tempG2_post)
      PB_R_null = rbind(PB_R_null,tempR_post)
    }
    
  }
}

save(list=ls(),file="Output_files/RData/Mmatrix_Null.RData")