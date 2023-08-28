#evaluate how elevation influences moisture content and wood density

install.packages("ggplot2")
install.packages("lifecycle", dependencies = T)
install.packages("Matrix")
install.packages("rlang")
install.packages("plyr")


#load packages
library(lifecycle)
library(ggplot2)
library(extrafont)
loadfonts(device="win")
library(remotes)
library(Matrix)
library(rlang)
library(lme4)
library(readr)
library(dplyr)
library(plyr)

#clear workspace
rm(list=ls())

#create theme basis for ggplot
theme_basis<-theme(axis.title = element_text(family = "Arial", color="black", face="bold", size=14))+
  theme(axis.text=element_text(family = "Arial",colour="black", face ="bold", size = 12)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"),axis.ticks = element_line(colour="black"))

#change working directory to import metadata 
# setwd('G:/My Drive/Tupper_fellowship/Interns/Angela_Barrera_Bello/Ecosystems_revision')
setwd('/cloud/project/')
metadata <- read.csv("/cloud/project/db_respiration_complete_20220909.csv", sep= ";")
metadata_df <- as.data.frame(metadata)
metadata$wood_ID <- as.character(metadata$wood_ID)
str(metadata)

#look at the distribution of samples
table(metadata$wood_sample)
table(metadata$wood_ID)
table(metadata$treatment)

#create separate downup and elevateddowned variables
metadata$downup <- ifelse(grepl("d",metadata$treatment),"bottom","top")
metadata$elevated_status <- ifelse(grepl("g",metadata$treatment),"downed","elevated")

#organize factors
metadata$elevated_status <- factor(metadata$elevated_status,levels=c("elevated","downed"))

#remove rows lacking data
metadata<- metadata[!is.na(metadata$wood_ID),]

# charact to num (weight)
metadata$initial_weight <- as.numeric(as.character(gsub(",",".",metadata$initial_weight)))
metadata$final_weight <- as.numeric(as.character(gsub(",",".",metadata$final_weight)))

#add colum moisture content as (wet-dry)/dry * 100

metadata$moist_cont_dry <- ((metadata$initial_weight- metadata$final_weight)/ metadata$final_weight)*100


#fix , to .

metadata$wood_density <- as.numeric(as.character(gsub(",",".",metadata$wood_density)))
metadata$moisture_content <- as.numeric(as.character(gsub(",",".",metadata$moisture_content)))
metadata$moist_cont_dry <- as.numeric(as.character(gsub(",",".",metadata$moist_cont_dry)))



#remove moisture_content = NA 
metadata_moist_no_0 <- filter(metadata, !is.na(moisture_content))

#remove moist_cont_dry = NA 
metadata_moistd_no_0 <- filter(metadata, !is.na(moist_cont_dry))


#run models
#-- -- -- -- -- -- -- -- -- 
#### Moisture content #### 
#-- -- -- -- -- -- -- -- --


metadata_moist_no_0 %>%
  group_by(elevated_status) %>%
  summarise_at(vars(moisture_content), list(name = mean))

#elevated         39.4
#downed           44.8



mod_moist1 <- lmer(moisture_content ~ elevated_status + downup + elevated_status:downup + (1|wood_ID),
                   data = metadata_moist_no_0)
mod_moist2 <- lmer(moisture_content ~ elevated_status + downup + (1|wood_ID), # AIC 629.43
                   data = metadata_moist_no_0)
mod_moist3 <- lmer(moisture_content ~ elevated_status + (1|wood_ID), # AIC 627.98
                   data = metadata_moist_no_0)
mod_moist4 <- lmer(moisture_content ~  downup + (1|wood_ID),
                   data = metadata_moist_no_0)
mod_moist5 <- lmer(moisture_content ~  1 + (1|wood_ID),
                   data = metadata_moist_no_0)

#the analyses below test whether both models similarly explain the data
anova(mod_moist1,mod_moist2) 
#if the above analysis is significant, do not move forward
anova(mod_moist2,mod_moist3)
anova(mod_moist2,mod_moist4) #*significant  
#if both of the above are significant, then do not move forward
#if one of the two above are significant, then compare the model with that variable to the null model
anova(mod_moist3,mod_moist5) 
anova(mod_moist4,mod_moist5)

summary(mod_moist3)

#examine the residuals of the best fit model 
plot(mod_moist3) # 
qqnorm(resid(mod_moist3));qqline(resid(mod_moist3))
hist(resid(mod_moist3)) # 

#Interpretation: Being downed in the forest floor affects wood moisture content (mod_moist3)


#plot the moisture content 
moisture_plot_trim_all <- ggplot(metadata_moist_no_0, aes(y = moisture_content, x = elevated_status)) +
  geom_boxplot(aes(fill = downup),color = "black") +
  theme_basis  + theme(legend.title=element_blank(),legend.background = element_blank(),
                       legend.key = element_blank(), legend.position = c(.15,.95),
                       legend.text = element_text(family = "Arial",colour="black", face ="bold", size = 10))+
  scale_x_discrete(name = element_blank(),labels=c("Suspended","Downed")) +
  scale_y_continuous(name = "Moisture content (%)",limits = c(10,72),breaks = seq(10,70,10),expand = c(0,0),
                     labels = c("10","","30","","50","","70")) +
  scale_fill_manual(values=c("grey40","white"))

moisture_plot_trim_all

#set directory to save figure 
setwd('C:/Users/angel/OneDrive/2_respiration_activity_20220816/figures/respiration_experiment')
ggsave("moisture_plot_trim_all_20220909.tiff",moisture_plot_trim_29_45,dpi = 600, compression = "lzw",width = 3.25,height = 2,scale = 1.8)


#summary sd se mean ALL ("observer")
mc_se_sd_mean_all <- ddply(metadata_moist_no_0, .(observer), summarise, 
                              M = mean(moisture_content), SE = sd(moisture_content) / sqrt((length(moisture_content))), 
                              SD = sd(moisture_content), Median=median(moisture_content))

#summary sd se mean by elevated_status
mc_se_sd_mean_updown <- ddply(metadata_moist_no_0, .(elevated_status), summarise, 
                                    M = mean(moisture_content), SE = sd(moisture_content) / sqrt((length(moisture_content))), 
                                    SD = sd(moisture_content), Median=median(moisture_content))

#summary sd se mean by elevated_status and downup
mc_se_sd_mean_elev_updown <- ddply(metadata_moist_no_0, .(elevated_status, downup), summarise, 
                                    M = mean(moisture_content), SE = sd(moisture_content) / sqrt((length(moisture_content))), 
                                    SD = sd(moisture_content), Median=median(moisture_content))


#-- -- -- -- -- -- -- -- -- 
#### Moisture content DRY #### 
#-- -- -- -- -- -- -- -- --

mod_moistd1 <- lmer(log(moist_cont_dry) ~ elevated_status + downup + elevated_status:downup + (1|wood_ID),
                   data = metadata_moistd_no_0)
mod_moistd2 <- lmer(log(moist_cont_dry) ~ elevated_status + downup + (1|wood_ID), 
                   data = metadata_moistd_no_0)
mod_moistd3 <- lmer(log(moist_cont_dry) ~ elevated_status + (1|wood_ID), 
                   data = metadata_moistd_no_0)
mod_moistd4 <- lmer(log(moist_cont_dry) ~  downup + (1|wood_ID),
                   data = metadata_moistd_no_0)
mod_moistd5 <- lmer(log(moist_cont_dry) ~  1 + (1|wood_ID),
                   data = metadata_moistd_no_0)

#the analyses below test whether both models similarly explain the data
anova(mod_moistd1,mod_moistd2) 
#if the above analysis is significant, do not move forward
anova(mod_moistd2,mod_moistd3)
anova(mod_moistd2,mod_moistd4) # Significant  
#if both of the above are significant, then do not move forward
#if one of the two above are significant, then compare the model with that variable to the null model
anova(mod_moistd3,mod_moistd5) 
anova(mod_moistd4,mod_moistd5)

summary(mod_moistd3)

#examine the residuals of the best fit model 
plot(mod_moistd3) # 
qqnorm(resid(mod_moistd3));qqline(resid(mod_moistd3))
hist(resid(mod_moistd3)) # 

#Interpretation: Being downed in the forest floor affects wood moisture content (mod_moist3)


#plot the moisture content 
moistured_plot_trim_all <- ggplot(metadata_moistd_no_0, aes(y = moist_cont_dry, x = elevated_status)) +
  geom_boxplot(aes(fill = downup),color = "black") +
  theme_basis  + theme(legend.title=element_blank(),legend.background = element_blank(),
                       legend.key = element_blank(), legend.position = c(.5,.95),
                       legend.text = element_text(family = "Arial",colour="black", face ="bold", size = 10))+
  scale_x_discrete(name = element_blank(),labels=c("Suspended","Downed")) +
  scale_y_log10(name = "Moisture content (%)",limits = c(10,240),breaks = c(seq(10,100,10),200),expand = c(0,0),
                     labels = c("10",rep("",8),"100","200")) +
  scale_fill_manual(values=c("grey40","white"))

moistured_plot_trim_all

#set directory to save figure 
setwd('C:/Users/angel/OneDrive/2_respiration_activity_20220816/figures/respiration_experiment')
ggsave("moistured_plot_trim_all_20220909.tiff",moistured_plot_trim_29_45,dpi = 600, compression = "lzw",width = 3.25,height = 2,scale = 1.8)


#summary sd se mean ALL ("observer")
mcd_se_sd_mean_all <- ddply(metadata_moistd_no_0, .(observer), summarise, 
                           M = mean(moist_cont_dry), SE = sd(moist_cont_dry) / sqrt((length(moist_cont_dry))), 
                           SD = sd(moist_cont_dry), Median=median(moist_cont_dry))



#summary sd se mean by elevated_status
mcd_se_sd_mean_elev <- ddply(metadata_moistd_no_0, .(elevated_status), summarise, 
                              M = mean(moist_cont_dry), SE = sd(moist_cont_dry) / sqrt((length(moist_cont_dry))), 
                              SD = sd(moist_cont_dry), Median=median(moist_cont_dry))

#summary sd se mean by elevated_status and downup
mcd_se_sd_mean_elev_updown <- ddply(metadata_moistd_no_0, .(elevated_status, downup), summarise, 
                               M = mean(moist_cont_dry), SE = sd(moist_cont_dry) / sqrt((length(moist_cont_dry))), 
                               SD = sd(moist_cont_dry), Median=median(moist_cont_dry))
#-- -- -- -- -- -- -- 
#### MC wet and dry: 2 panel plot ####  
#-- -- -- -- -- -- -- 

install.packages("cowplot")
library("cowplot")


plot_grid(moisture_plot_trim_all, moistured_plot_trim_all, align = "v", ncol = 1, rel_heights = c(.5, .5))


#-- -- -- -- -- -- -- 
#### Wood density ####  
#-- -- -- -- -- -- -- 
str(metadata$wood_density)

mod_dens1 <- lmer(wood_density ~ elevated_status + downup + elevated_status:downup + (1|wood_ID),
                  data = metadata)
mod_dens2 <- lmer(wood_density ~ elevated_status + downup + (1|wood_ID),
                  data = metadata)
mod_dens3 <- lmer(wood_density ~ elevated_status + (1|wood_ID),
                  data = metadata)
mod_dens4 <- lmer(wood_density ~  downup + (1|wood_ID),
                  data = metadata)
mod_dens5 <- lmer(wood_density ~  1 + (1|wood_ID),
                  data = metadata)

#the analyses below test whether both models similarly explain the data
anova(mod_dens1,mod_dens2)
#if the above analysis is significant, do not move forward
anova(mod_dens2,mod_dens3)
anova(mod_dens2,mod_dens4) 
#if both of the above are significant, then do not move forward
#if one of the two above are significant, then compare the model with that variable to the null model
anova(mod_dens3,mod_dens5) 
anova(mod_dens4,mod_dens5)

#examine the residuals of the best fit model (none with significant differences, but lower p value is
#mod_dens1, p > 0.3714)
plot(mod_dens3)
hist(resid(mod_dens3))

#plot the moisture content 
density_plot_all <- ggplot(metadata, aes(y = wood_density, x = elevated_status)) +
  geom_boxplot(aes(fill = downup),color = "black") +
  theme_basis  + theme(legend.title=element_blank(),legend.background = element_blank(),
                       legend.key = element_blank(), legend.position = c(.15,.95),
                       legend.text = element_text(family = "Arial",colour="black", face ="bold", size = 10))+ 
  scale_x_discrete(name = element_blank(),labels=c("Elevated","Downed")) +
  scale_y_continuous(name = expression(bold("Wood density"~(g~mm^{-3}))),
                     limits = c(.15,.65),breaks = seq(0.2,0.6,.1),expand = c(0,0),
                     labels = c("0.2","","0.4","","0.6")) + scale_fill_manual(values=c("grey40","white"))

density_plot_all

#save figure
setwd('C:/Users/angel/OneDrive/2_respiration_activity_20220816/figures/respiration_experiment')
ggsave("density_plot_all_20220909.tiff",density_plot,dpi = 600, compression = "lzw",width = 3.25,height = 2,scale = 1.8)

#-- -- -- -- -- -- -- 
#### Termites #### 
#-- -- -- -- -- -- -- 
str(metadata$termites)

#Fisher test to test differences among treatments 

#create a numeric termite dataframe
metadata$termite_num <- ifelse(metadata$termites=="yes",1,0)

#Fisher test to test differences among treatments 

#use logistic regression to include the random effect
termite_mod1 <- glmer(termite_num ~ elevated_status + downup + elevated_status:downup + (1|wood_ID),
                   family=binomial(link="logit"), data = metadata)
termite_mod2 <- glmer(termite_num ~ elevated_status + downup + (1|wood_ID),
                   family=binomial(link="logit"), data = metadata)
termite_mod3 <- glmer(termite_num ~ elevated_status +  (1|wood_ID),
                   family=binomial(link="logit"), data = metadata)
termite_mod4 <- glmer(termite_num ~ downup + (1|wood_ID),
                   family=binomial(link="logit"), data = metadata)

anova(termite_mod1,termite_mod2,test="LRT")
anova(termite_mod3,termite_mod2,test="LRT") #(!)
anova(termite_mod4,termite_mod2,test="LRT")

#create full contingency table for fisher.test()
table(metadata$treatment,metadata$termite_num)
fisher.test(table(metadata$treatment,metadata$termite_num))

#create contingency table for fisher.test() of elevated versus downed
table(metadata$elevated_status,metadata$termite_num)
fisher.test(table(metadata$elevated_status,metadata$termite_num))

#create contingency table for fisher.test() of top versus bottom
table(metadata$downup,metadata$termite_num)
fisher.test(table(metadata$downup,metadata$termite_num))

#calculate 95% confidence intervals for the termite data
library(DescTools)
termites_pres <- table(metadata$treatment,metadata$termite_num)[,2]
termites_total <- table(metadata$treatment)
termite_df<- data.frame(BinomCI(termites_pres,termites_total,method = "clopper-pearson"))

#add unique values for each treatment
termite_df$treatment <- names(termites_pres)
termite_df$downup <- ifelse(grepl("d",termite_df$treatment),"bottom","top")
termite_df$elevated_status <- ifelse(grepl("g",termite_df$treatment),"downed","elevated")

#organize factors
termite_df$elevated_status <- factor(termite_df$elevated_status,levels=c("elevated","downed"))

#plot  
termite_plot_all <- ggplot(termite_df, aes(y = est, x = elevated_status)) +
  geom_point(aes(shape = downup),position=position_dodge(width=.3),size = 3) +
  geom_errorbar(aes(x = elevated_status, y = est, ymin = lwr.ci, ymax = upr.ci,
                    group = downup), width = .2,position=position_dodge(width=.3)) +
  scale_shape_manual(values = c(1,16)) +
  theme_basis  + theme(legend.title=element_blank(),legend.background = element_blank(),
                       legend.key = element_blank(), legend.position = c(.13,.95),
                       legend.text = element_text(family = "Arial",colour="black", face ="bold", size = 10))+ 
  scale_x_discrete(name = element_blank(),labels=c("Suspended","Downed")) +
  scale_y_continuous(name = "Termites present (%)",
                     limits = c(.3,.9),breaks = seq(0.2,0.9,.1),expand = c(0,0),
                     labels = c("20","","40","","60","","80",""))

termite_plot_all

#save figure
setwd('C:/Users/angel/OneDrive/2_respiration_activity_20220816/figures/respiration_experiment')
ggsave("termite_plot_all_20220909.tiff",termite_plot,dpi = 600, compression = "lzw",width = 3.25,height = 2,scale = 1.8)

#mean

metadata %>%
  group_by(elevated_status) %>%
  summarise_at(vars(termite_num), list(name = mean))

#elevated        0.122
#downed          0.056

#se, sd, mean
termite_se <- ddply(metadata, .(elevated_status), summarise, 
                  M = mean(termite_num), SE = sd(termite_num) / sqrt((length(termite_num))), 
                  SD = sd(termite_num))

#Across all samples
sum(metadata$termite_num) 
nrow(metadata)

all_ter <- (sum(metadata$termite_num)/nrow(metadata))*100 
# 58

#-- -- -- -- -- -- -- 
#### Hyphae #### 
#-- -- -- -- -- -- -- 
str(metadata$visible_hyphae)

#create a numeric termite dataframe
metadata$hyphae_num <- ifelse(metadata$visible_hyphae=="yes",1,0)

#use logistic regression to include the random effect
hyph_mod1 <- glmer(hyphae_num ~ elevated_status + downup + elevated_status:downup + (1|wood_ID),
                   family=binomial(link="logit"), data = metadata)
hyph_mod2 <- glmer(hyphae_num ~ elevated_status + downup + (1|wood_ID),
                   family=binomial(link="logit"), data = metadata)
hyph_mod3 <- glmer(hyphae_num ~ elevated_status +  (1|wood_ID),
                   family=binomial(link="logit"), data = metadata)
hyph_mod4 <- glmer(hyphae_num ~ downup + (1|wood_ID),
                   family=binomial(link="logit"), data = metadata)

anova(hyph_mod1,hyph_mod2,test="LRT")
anova(hyph_mod3,hyph_mod2,test="LRT") #(!)
anova(hyph_mod4,hyph_mod2,test="LRT") #(!)

# #create full contingency table for fisher.test()
# table(metadata$treatment,metadata$hyphae_num)
# fisher.test(table(metadata$treatment,metadata$hyphae_num))
# 
# #### -- TEST THE TWO MAJOR TERMS
# #compare each pair of treatments
# # elevated versus downed (!)
# table(metadata$elevated_status,metadata$hyphae_num)
# fisher.test(table(metadata$elevated_status,metadata$hyphae_num))
# 
# #compare each pair of treatments
# # elevated versus downed (!)
# table(metadata$downup,metadata$hyphae_num)
# fisher.test(table(metadata$downup,metadata$hyphae_num))

# #### -- PAIRWISE COMPARISONS
# #compare each pair of treatments
# # ed versus eu
# table(metadata$treatment,metadata$hyphae_num)[c(1,2),]
# fisher.test(table(metadata$treatment,metadata$hyphae_num)[c(1,2),])
# 
# # ed versus gd
# table(metadata$treatment,metadata$hyphae_num)[c(1,3),]
# fisher.test(table(metadata$treatment,metadata$hyphae_num)[c(1,3),])
# 
# # ed versus gu
# table(metadata$treatment,metadata$hyphae_num)[c(1,4),]
# fisher.test(table(metadata$treatment,metadata$hyphae_num)[c(1,4),])
# 
# # eu versus gd (!)
# table(metadata$treatment,metadata$hyphae_num)[c(2,3),]
# fisher.test(table(metadata$treatment,metadata$hyphae_num)[c(2,3),])
# 
# # eu versus gu
# table(metadata$treatment,metadata$hyphae_num)[c(2,4),]
# fisher.test(table(metadata$treatment,metadata$hyphae_num)[c(2,4),])
# 
# # gd versus gu
# table(metadata$treatment,metadata$hyphae_num)[c(3,4),]
# fisher.test(table(metadata$treatment,metadata$hyphae_num)[c(3,4),])

#plot data

#calculate 95% confidence intervals for the termite data
library(DescTools)
hyphae_pres <- table(metadata$treatment,metadata$hyphae_num)[,2]
hyphae_total <- table(metadata$treatment)
hyphae_df<- data.frame(BinomCI(hyphae_pres,hyphae_total,method = "clopper-pearson"))

#add unique values for each treatment
hyphae_df$treatment <- names(hyphae_pres)
hyphae_df$downup <- ifelse(grepl("d",hyphae_df$treatment),"bottom","top")
hyphae_df$elevated_status <- ifelse(grepl("g",hyphae_df$treatment),"downed","elevated")

#organize factors
hyphae_df$elevated_status <- factor(hyphae_df$elevated_status,levels=c("elevated","downed"))

#plot 
hyphae_plot <- ggplot(hyphae_df, aes(y = est, x = elevated_status)) +
  geom_point(aes(shape = downup),position=position_dodge(width=.3),size = 3) +
  geom_errorbar(aes(x = elevated_status, y = est, ymin = lwr.ci, ymax = upr.ci,
                    group = downup), width = .2,position=position_dodge(width=.3)) +
  scale_shape_manual(values = c(1,16)) +
  theme_basis  + theme(legend.title=element_blank(),legend.background = element_blank(),
                       legend.key = element_blank(), legend.position = c(.13,.95),
                       legend.text = element_text(family = "Arial",colour="black", face ="bold", size = 10))+ 
  scale_x_discrete(name = element_blank(),labels=c("Suspended","Downed")) +
  scale_y_continuous(name = "Hyphae present (%)",
                     limits = c(0,.8),breaks = seq(0,0.9,.1),expand = c(0,0),
                     labels = c("0","","20","","40","","60","","80",""))

hyphae_plot

#save figure
setwd('D:/indep_project_wood_currentlyworking/2_respiration_activity/figures/respiration_experiment')
ggsave("hyphae_plot_29_45.tiff",hyphae_plot,dpi = 600, compression = "lzw",width = 3.25,height = 2,scale = 1.8)

#mean
metadata %>%
  group_by(elevated_status) %>%
  summarise_at(vars(hyphae_num), list(name = mean))

          #elevated        0.322
          #downed          0.478

#Across all samples
sum(metadata$hyphae_num) 
nrow(metadata)

all_hyp <- (sum(metadata$hyphae_num)/nrow(metadata))*100 
# 40


#se (standar error)

      #se, sd, mean
hyp_se <- ddply(metadata, .(elevated_status), summarise, 
              M = mean(hyphae_num), SE = sd(hyphae_num) / sqrt((length(hyphae_num))), 
              SD = sd(hyphae_num))

hyp_sum_upd <- ddply(metadata, .(downup), summarise, 
                M = mean(hyphae_num), SE = sd(hyphae_num) / sqrt((length(hyphae_num))), 
                SD = sd(hyphae_num))


#-- -- -- -- -- -- -- 
#### Algae #### 
#-- -- -- -- -- -- -- 

#create a numeric termite dataframe
metadata$algae_num <- ifelse(metadata$visible_algae_lichen=="yes",1,0)

summary(metadata$algae_num)#no NA values
metadata %>% group_by(elevated_status,downup) %>% summarise(N = length(algae_num),algae_count = sum(algae_num))

#use logistic regression to include the random effect
algae_mod1 <- glmer(algae_num ~ elevated_status + downup + elevated_status:downup + (1|wood_ID),
                   family=binomial(link="logit"), data = metadata)
algae_mod2 <- glmer(algae_num ~ elevated_status + downup + (1|wood_ID),
                   family=binomial(link="logit"), data = metadata)
algae_mod3 <- glmer(algae_num ~ elevated_status +  (1|wood_ID),
                   family=binomial(link="logit"), data = metadata)
algae_mod4 <- glmer(algae_num ~ downup + (1|wood_ID),
                   family=binomial(link="logit"), data = metadata)

anova(algae_mod1,algae_mod2,test="LRT")
anova(algae_mod3,algae_mod2,test="LRT") #(!)
anova(algae_mod4,algae_mod2,test="LRT") #(!)

# #Fisher test to test differences among treatments 
# 
# #create full contingency table for fisher.test()
# table(metadata$treatment,metadata$algae_num)
# fisher.test(table(metadata$treatment,metadata$algae_num))
# 
# #compare each pair of treatments
# # ed versus eu
# table(metadata$treatment,metadata$algae_num)[c(1,2),]
# fisher.test(table(metadata$treatment,metadata$algae_num)[c(1,2),])
# 
# # ed versus gd
# table(metadata$treatment,metadata$algae_num)[c(1,3),]
# fisher.test(table(metadata$treatment,metadata$algae_num)[c(1,3),])
# 
# # ed versus gu
# table(metadata$treatment,metadata$algae_num)[c(1,4),]
# fisher.test(table(metadata$treatment,metadata$algae_num)[c(1,4),])
# 
# # eu versus gd
# table(metadata$treatment,metadata$algae_num)[c(2,3),]
# fisher.test(table(metadata$treatment,metadata$algae_num)[c(2,3),])
# 
# # eu versus gu
# table(metadata$treatment,metadata$algae_num)[c(2,4),]
# fisher.test(table(metadata$treatment,metadata$algae_num)[c(2,4),])
# 
# # gd versus gu
# table(metadata$treatment,metadata$algae_num)[c(3,4),]
# fisher.test(table(metadata$treatment,metadata$algae_num)[c(3,4),])
# 
# str(metadata$visible_algae_lichen)

#plot data

#calculate 95% confidence intervals for the termite data
library(DescTools)
algae_pres <- table(metadata$treatment,metadata$algae_num)[,2]
algae_total <- table(metadata$treatment)
algae_df<- data.frame(BinomCI(algae_pres,algae_total,method = "clopper-pearson"))

#add unique values for each treatment
algae_df$treatment <- names(algae_pres)
algae_df$downup <- ifelse(grepl("d",algae_df$treatment),"bottom","top")
algae_df$elevated_status <- ifelse(grepl("g",algae_df$treatment),"downed","elevated")

#organize factors
algae_df$elevated_status <- factor(algae_df$elevated_status,levels=c("elevated","downed"))

#plot the moisture content 
algae_plot <- ggplot(algae_df, aes(y = est, x = elevated_status)) +
  geom_point(aes(shape = downup),position=position_dodge(width=.3),size = 3) +
  geom_errorbar(aes(x = elevated_status, y = est, ymin = lwr.ci, ymax = upr.ci,
                    group = downup), width = .2,position=position_dodge(width=.3)) +
  scale_shape_manual(values = c(1,16)) +
  theme_basis  + theme(legend.title=element_blank(),legend.background = element_blank(),
                       legend.key = element_blank(), legend.position = c(.13,.95),
                       legend.text = element_text(family = "Arial",colour="black", face ="bold", size = 10))+ 
  scale_x_discrete(name = element_blank(),labels=c("Suspended","Downed")) +
  scale_y_continuous(name = "Photosynthetic growth present (%)",
                     limits = c(-.02,.35),breaks = seq(-.05,0.4,.05),expand = c(0,0),
                     labels = c("","0","","10","","20","","30","","40"))

algae_plot

#save figure
setwd('D:/indep_project_wood_currentlyworking/2_respiration_activity/figures/respiration_experiment')
ggsave("algae_plot_29_45.tiff",algae_plot,dpi = 600, compression = "lzw",width = 3.25,height = 2,scale = 1.8)

#mean

metadata %>%
  group_by(elevated_status) %>%
  summarise_at(vars(algae_num), list(name = mean))

#elevated        0.122
#downed          0.056

      #se, sd, mean
algae_se <- ddply(metadata, .(elevated_status), summarise, 
                M = mean(algae_num), SE = sd(algae_num) / sqrt((length(algae_num))), 
                SD = sd(algae_num))

algae_sum_upd <- ddply(metadata, .(downup), summarise, 
                     M = mean(algae_num), SE = sd(algae_num) / sqrt((length(algae_num))), 
                     SD = sd(algae_num))



#Across all samples
sum(metadata$algae_num) 
nrow(metadata)

all_alg <- (sum(metadata$algae_num)/nrow(metadata))*100 
# 58

