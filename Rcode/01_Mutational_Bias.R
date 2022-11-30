rm(list = ls())
gc()
library(emmeans)
library(MCMCglmm)
library(psych)
library(gplots)
library(data.table)
library(matrixStats)
library(boot)
library(dplyr)
library(lme4)
library(performance)
library(car)

## This code produce the phenotypic analysis of the MA lines from the 2 founder genotype produce figure 1 of the manuscript

# Load the phenotypic data 
final_merged=read.table("data/Final_merged_data_MA_Lines.txt",sep="\t",h=TRUE)
# Define the 6 traits names
vect_P_traits = c("T12", "T13", "T21", "T23", "T31","T32")

# Rename some columns 
final_merged$pop_label2 <- final_merged$pop_label
final_merged$pop_label2[final_merged$pop_label=="PBancL0"] <- paste0(final_merged$pop_label2[final_merged$pop_label=="PBancL0"],1:16)
final_merged$pop_label2[final_merged$pop_label=="N2ancL0"] <- paste0(final_merged$pop_label2[final_merged$pop_label=="N2ancL0"],1:15)

# Create a pop coulmn - N2anc/N2lines/PBanc/PBlines
final_merged$pop <- "N2anc"
final_merged$pop[substring(final_merged$pop_label,1,2)=="N2" & final_merged$pop_label!="N2ancL0"] <- "N2lines"
final_merged$pop[substring(final_merged$pop_label,1,2)=="PB" & final_merged$pop_label!="N2ancL0"] <- "PBlines"
final_merged$pop[final_merged$pop_label=="PBancL0"] <- "PBanc"

# Create two columns - is_line: is the data an ancestor (logical) and ancestor (factor: PB/N2)
final_merged$is_line = ifelse(final_merged$pop%in%c("N2anc","PBanc"),TRUE,FALSE)
final_merged$ancestor = ifelse(final_merged$pop%in%c("N2anc","N2lines"),"N2","PB")

### Center the environmental covariates
for(i in c("temperature","rel_humidity","logD")){
  final_merged[,i] <- (final_merged[,i]-mean(final_merged[,i]))#/sd(final_merged[,i])
}

### I will do the models separately and manually for each trait
i= "T12" # SF
mod_lm1 <- lmer(final_merged[,i] ~ (temperature+rel_humidity+logD) + ancestor*is_line + (1|date_str) + (1|pop_label2), data= final_merged)
mod_lm2 <- lmer(final_merged[,i] ~ (temperature+rel_humidity+logD) + ancestor+is_line + (1|date_str) + (1|pop_label2), data= final_merged)
mod_lm3 <- lmer(final_merged[,i] ~ (temperature+rel_humidity+logD) + ancestor + (1|date_str) + (1|pop_label2), data= final_merged)
mod_lm4 <- lmer(final_merged[,i] ~ (temperature+rel_humidity+logD)  + (1|date_str) + (1|pop_label2), data= final_merged)
anova(mod_lm1,mod_lm2) # not significant
anova(mod_lm2,mod_lm3) # not significant
anova(mod_lm3,mod_lm4) # sign.

leveneTest(resid(mod_lm2)~final_merged$pop) # OK
#leveneTest(resid(mod_lm2)~final_merged$ancestor) # OK


i= "T13" # SB
mod_lm1 <- lmer(final_merged[,i] ~ (temperature+rel_humidity+logD) + ancestor*is_line + (1|date_str) + (1|pop_label2), data= final_merged)
mod_lm2 <- lmer(final_merged[,i] ~ (temperature+rel_humidity+logD) + ancestor+is_line + (1|date_str) + (1|pop_label2), data= final_merged)
mod_lm3 <- lmer(final_merged[,i] ~ (temperature+rel_humidity+logD) + ancestor + (1|date_str) + (1|pop_label2), data= final_merged)
mod_lm4 <- lmer(final_merged[,i] ~ (temperature+rel_humidity+logD)  + (1|date_str) + (1|pop_label2), data= final_merged)
anova(mod_lm1,mod_lm2) # not significant
anova(mod_lm2,mod_lm3) # not significant
anova(mod_lm3,mod_lm4) # sign.

leveneTest(resid(mod_lm2)~final_merged$pop) # OK
#leveneTest(resid(mod_lm2)~final_merged$ancestor) # OK

i= "T21" # FS
mod_lm1 <- lmer(final_merged[,i] ~ (temperature+rel_humidity+logD) + ancestor*is_line + (1|date_str) + (1|pop_label2), data= final_merged)
mod_lm2 <- lmer(final_merged[,i] ~ (temperature+rel_humidity+logD) + ancestor+is_line + (1|date_str) + (1|pop_label2), data= final_merged)
mod_lm3 <- lmer(final_merged[,i] ~ (temperature+rel_humidity+logD) + ancestor + (1|date_str) + (1|pop_label2), data= final_merged)
mod_lm3b <- lmer(final_merged[,i] ~ (temperature+rel_humidity+logD) + is_line + (1|date_str) + (1|pop_label2), data= final_merged)
anova(mod_lm1,mod_lm2) # not significant
anova(mod_lm2,mod_lm3) #significant 0.034
anova(mod_lm2,mod_lm3b) #significant

leveneTest(resid(mod_lm2)~final_merged$pop) # OK
leveneTest(resid(mod_lm2)~final_merged$ancestor) # OK

i= "T23" # FB
mod_lm1 <- lmer(final_merged[,i] ~ (temperature+rel_humidity+logD) + ancestor*is_line + (1|date_str) + (1|pop_label2), data= final_merged)
mod_lm2 <- lmer(final_merged[,i] ~ (temperature+rel_humidity+logD) + ancestor+is_line + (1|date_str) + (1|pop_label2), data= final_merged)
mod_lm3 <- lmer(final_merged[,i] ~ (temperature+rel_humidity+logD) + ancestor + (1|date_str) + (1|pop_label2), data= final_merged)
mod_lm4 <- lmer(final_merged[,i] ~ (temperature+rel_humidity+logD)  + (1|date_str) + (1|pop_label2), data= final_merged)
anova(mod_lm1,mod_lm2) # not significant
anova(mod_lm2,mod_lm3) # not significant
anova(mod_lm3,mod_lm4) # not sign.

leveneTest(resid(mod_lm2)~final_merged$pop) # OK
leveneTest(resid(mod_lm2)~final_merged$ancestor) # OK


i= "T31" # BS
mod_lm1 <- lmer(final_merged[,i] ~ (temperature+rel_humidity+logD) + ancestor*is_line + (1|date_str) + (1|pop_label2), data= final_merged)
mod_lm2 <- lmer(final_merged[,i] ~ (temperature+rel_humidity+logD) + ancestor+is_line + (1|date_str) + (1|pop_label2), data= final_merged)
mod_lm3 <- lmer(final_merged[,i] ~ (temperature+rel_humidity+logD) + ancestor + (1|date_str) + (1|pop_label2), data= final_merged)
mod_lm3b <- lmer(final_merged[,i] ~ (temperature+rel_humidity+logD) + is_line + (1|date_str) + (1|pop_label2), data= final_merged)
anova(mod_lm1,mod_lm2) # not significant
anova(mod_lm2,mod_lm3) #significant 0.038
anova(mod_lm2,mod_lm3b) # not sign.

leveneTest(resid(mod_lm2)~final_merged$pop) # OK
leveneTest(resid(mod_lm2)~final_merged$ancestor) # OK

i= "T32" # BF
mod_lm1 <- lmer(final_merged[,i] ~ (temperature+rel_humidity+logD) + ancestor*is_line + (1|date_str) + (1|pop_label2), data= final_merged)
mod_lm2 <- lmer(final_merged[,i] ~ (temperature+rel_humidity+logD) + ancestor+is_line + (1|date_str) + (1|pop_label2), data= final_merged)
mod_lm3 <- lmer(final_merged[,i] ~ (temperature+rel_humidity+logD) + ancestor + (1|date_str) + (1|pop_label2), data= final_merged)
mod_lm3b <- lmer(final_merged[,i] ~ (temperature+rel_humidity+logD) + is_line + (1|date_str) + (1|pop_label2), data= final_merged)
mod_lm4 <- lmer(final_merged[,i] ~ (temperature+rel_humidity+logD)  + (1|date_str) + (1|pop_label2)-1, data= final_merged)
anova(mod_lm1,mod_lm2) # not significant
anova(mod_lm2,mod_lm3) # not significant 0.099 (marginally)
anova(mod_lm3,mod_lm4) # sign ++

leveneTest(resid(mod_lm2)~final_merged$pop) # OK
leveneTest(resid(mod_lm2)~final_merged$ancestor) # OK


## Produce Figure 1


###########
pdf("plots/Mutational_Bias_All_Traits_Single_model_with_covariates.pdf")
kkk=0
par(mfrow=c(2,3))
for(i in vect_P_traits){
  kkk=kkk+1
  
  par(mar=c(6,4,4,2))
  
  #########
  mod_lm1 <- lmer(final_merged[final_merged$logD>-1,i] ~ (temperature+rel_humidity+logD) + ancestor*is_line-1 + (1|date_str) + (1|pop_label2), data= subset(final_merged,logD>-1))
  conf_int_anc = confint(mod_lm1)[7:8,]
  new_final = data.frame(pop=sort(unique(final_merged$pop)),is_line=c(FALSE,TRUE,FALSE,TRUE),temperature=0,rel_humidity=0,logD=0,ancestor=rep(c("N2","PB"),each=2),date_str=NA,pop_label2=NA)
  predict_mean=predict(mod_lm1,new_final,re.form=NA)
  new_final$mean=predict_mean;rm(predict_mean)
  is_N2_temp <- substring(rownames(ranef(mod_lm1)$pop_label2),1,2)=="N2"
  predict_lines= ifelse(is_N2_temp,new_final$mean[new_final$pop=="N2lines"],new_final$mean[new_final$pop=="PBlines"]) + ranef(mod_lm1)$pop_label2
  
  plot(new_final$mean~c(0,1,2,3),xlim=c(-0.5,3.5),bty="n",ylab=expression(paste("Transition rate - ",log(q[i]))),
       xlab="",xaxt="n",ylim=cbind(c(-4,0),c(-3,0),c(-3,1),c(-6.5,-3),c(-1,1),c(-5,0))[,kkk],pch=21,bg="red",type="n",main=c("S->F","S->B","F->S","F->B","B->S","B->F")[kkk])
  axis(1,at=0:3,c("N2","N2 lines","PB306","PB306 lines"),las=2)
  
  temp_df <- subset(final_merged,pop%in%c('N2anc'))[,i]
  points(temp_df~I(jitter(rep(4,length(temp_df)))-4),pch=16,col="grey")
  temp_df <- subset(final_merged,pop%in%c('N2lines'))[,i]
  points(temp_df~I(jitter(rep(2,length(temp_df)))-1),pch=16,col="grey")
  points(predict_lines[is_N2_temp,1]~jitter(rep(1,sum(is_N2_temp))),pch=16)
  
  
  arrows(0, conf_int_anc[1,1], 0 ,conf_int_anc[1,2],code=3,length=.05,angle=90)
  
  temp_df <- subset(final_merged,pop%in%c('PBanc'))[,i]
  points(temp_df~jitter(rep(2,length(temp_df))),pch=16,col="grey")
  temp_df <- subset(final_merged,pop%in%c('PBlines'))[,i]
  points(temp_df~jitter(rep(3,length(temp_df))),pch=16,col="grey")
  points(predict_lines[!is_N2_temp,1]~jitter(rep(3,sum(!is_N2_temp))),pch=16)
  
  arrows(2, conf_int_anc[2,1], 2 ,conf_int_anc[2,2],code=3,length=.05,angle=90)
  
  
  points(new_final$mean~c(0:3) ,pch=21,bg="red")
}
dev.off()




