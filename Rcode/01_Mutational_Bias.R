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

final_merged=read.table("data/Final_merged_data_MA_Lines.txt",sep="\t",h=TRUE)
vect_P_traits = c("T12", "T13", "T21", "T23", "T31","T32")

final_merged$pop_label2 <- final_merged$pop_label
final_merged$pop_label2[final_merged$pop_label=="PBancL0"] <- paste0(final_merged$pop_label2[final_merged$pop_label=="PBancL0"],1:16)
final_merged$pop_label2[final_merged$pop_label=="N2ancL0"] <- paste0(final_merged$pop_label2[final_merged$pop_label=="N2ancL0"],1:15)

final_merged$pop <- "N2anc"
final_merged$pop[substring(final_merged$pop_label,1,2)=="N2" & final_merged$pop_label!="N2ancL0"] <- "N2lines"
final_merged$pop[substring(final_merged$pop_label,1,2)=="PB" & final_merged$pop_label!="N2ancL0"] <- "PBlines"
final_merged$pop[final_merged$pop_label=="PBancL0"] <- "PBanc"

###

mutational_bias=NULL
for(i in vect_P_traits){

print(i)

temp_df <- subset(final_merged,pop%in%c('N2anc','N2lines'))
mod_lm1 <- lmer(temp_df[,i] ~   pop +  (1|date_str) + (1|pop_label2), data= temp_df)
mod_lm1_1 <- lmer(temp_df[,i] ~ 1  +(1|date_str) + (1|pop_label2), data= temp_df)

mutational_bias <- rbind(mutational_bias,
data.frame(anc="N2",trait=i,p_val1=anova(mod_lm1,mod_lm1_1)$Pr[2]))

temp_df <- subset(final_merged,pop%in%c('PBanc','PBlines'))
mod_lm1 <- lmer(temp_df[,i] ~   pop   + (1|date_str) + (1|pop_label2), data= temp_df)
mod_lm1_1 <- lmer(temp_df[,i] ~ 1  + (1|date_str) + (1|pop_label2), data= temp_df)

mutational_bias <- rbind(mutational_bias,
data.frame(anc="PB",trait=i,p_val1=anova(mod_lm1,mod_lm1_1)$Pr[2]))

}
subset(mutational_bias,p_val1<0.05)
mutational_bias$p_val_adj <- p.adjust(mutational_bias$p_val1,method="BH")
subset(mutational_bias,p_val_adj<0.05)


pdf("plots/Mutational_Bias_All_Traits.pdf")
kkk=0
par(mfrow=c(2,3))
for(i in unique(mutational_bias$trait)[1:6]){
kkk=kkk+1

par(mar=c(6,4,4,2))

#########
temp_df <- subset(final_merged,pop%in%c('N2anc','N2lines'))
mod_lm1 <- lmer(temp_df[,i] ~   pop + (1|date_str) + (1|pop_label2), data= temp_df)
new_final = data.frame(pop=c("N2anc","N2lines"),ancestor="N2",date_str=NA,pop_label2=NA)
predict_N2_mean=predict(mod_lm1,new_final,re.form=NA)
predict_N2_lines=predict_N2_mean[2] + ranef(mod_lm1)$pop_label2

plot(predict_N2_mean~c(0,1),xlim=c(-0.5,3.5),bty="n",ylab="",
      xlab="",xaxt="n",ylim=cbind(c(-3,1),c(-3,0),c(-3,1),c(-6.5,-3),c(-1,1),c(-5,0))[,kkk],pch=21,bg="red",type="n",main=c("SF","SB","FS","FB","BS","BF")[kkk])

if(i %in% unique(mutational_bias$trait)[c(1:3)]){
axis(1,at=0:3,rep("",4),las=2)
}else{
  axis(1,at=0:3,c("N2","N2 lines","PB306","PB306 lines"),las=2)
}
#expression(paste("Transition rate - ",log(q[i])))

temp_df <- subset(final_merged,pop%in%c('N2anc'))[,i]
points(temp_df~I(jitter(rep(4,length(temp_df)))-4),pch=16,col="grey")
temp_df <- subset(final_merged,pop%in%c('N2lines'))[,i]
points(temp_df~I(jitter(rep(2,length(temp_df)))-1),pch=16,col="grey")
points(predict_N2_lines[,1]~jitter(rep(1,nrow(predict_N2_lines))),pch=16)
arrows(0, predict_N2_mean[1]- summary(mod_lm1)$coef[1,2], 0 ,predict_N2_mean[1] + summary(mod_lm1)$coef[1,2],code=3,length=.05,angle=90)
points(predict_N2_mean~c(0,1) ,pch=21,bg="red")

temp_df <- subset(final_merged,pop%in%c('PBanc','PBlines'))
mod_lm1 <- lmer(temp_df[,i] ~   pop + (1|date_str) + (1|pop_label2), data= temp_df)
new_final = data.frame(pop=c("PBanc","PBlines"),ancestor="PB",date_str=NA,pop_label2=NA)
predict_PB_mean=predict(mod_lm1,new_final,re.form=NA)
predict_PB_lines=predict_PB_mean[2] + ranef(mod_lm1)$pop_label2

temp_df <- subset(final_merged,pop%in%c('PBanc'))[,i]
points(temp_df~jitter(rep(2,length(temp_df))),pch=16,col="grey")
temp_df <- subset(final_merged,pop%in%c('PBlines'))[,i]
points(temp_df~jitter(rep(3,length(temp_df))),pch=16,col="grey")
points(predict_PB_lines[,1]~jitter(rep(3,nrow(predict_PB_lines))),pch=16)

arrows(2, predict_PB_mean[1]- summary(mod_lm1)$coef[1,2], 2 ,predict_PB_mean[1] + summary(mod_lm1)$coef[1,2],code=3,length=.05,angle=90)

if(i %in% unique(mutational_bias$trait)[c(3,5)]){
points(predict_PB_mean~c(2,3) ,pch=21,bg="red",type="b")
}else{
  points(predict_PB_mean~c(2,3) ,pch=21,bg="red")
}

}
dev.off()




