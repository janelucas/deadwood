#Create a script for analyzing respiration rates and metadata
#first import all datasets and packages
#Run models for moisture content, wood density and respiration
#Date: February 17 2022
#By: Evan Gora & Angela Barrera

install.packages("plyr")

#load packages
library(lme4)
library(ggplot2)
library(dplyr)
library(extrafont)
library(mctoolsr)
library(plyr)


loadfonts(device="win")

#create theme basis for ggplot
theme_basis<-theme(axis.title = element_text(family = "Arial", color="black", face="bold", size=14))+
  theme(axis.text=element_text(family = "Arial",colour="black", face ="bold", size = 12)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"),axis.ticks = element_line(colour="black"))

#set directory to import data
setwd('C:\\Users\\angel\\OneDrive\\Independent_project_wood_20230823\\2_respiration_activity_20220816\\results\\respiration_experiment')

#import the data
respiration_exp <- read.csv("respiration_exp.csv")
str(respiration_exp)
unique(respiration_exp$wood_ID) #15 to 45 

#Trim data (plumbers putty)  (wood_ID >= 29)
respiration_exp <- respiration_exp[respiration_exp$wood_ID >=29,]
unique(respiration_exp$wood_ID)

#### Organize data ####
#add final respiration value to each row
last_resp_df <- respiration_exp %>%
  group_by(wood_sample) %>%
  arrange(desc(hours_since_collection)) %>%
  slice(1)
last_resp_df_trime <- last_resp_df[,c("wood_sample","slope")]
colnames(last_resp_df_trime)[2]<-"final_respiration_rate"
respiration_exp <- left_join(respiration_exp,last_resp_df_trime) 
 
#plot respiration over time, as a percent of the final value
respiration_exp$resp_prcnt_of_final <- respiration_exp$slope/respiration_exp$final_respiration_rate
 
summary(respiration_exp$resp_prcnt_of_final==1);length(unique(respiration_exp$wood_sample))
ggplot(respiration_exp[respiration_exp$resp_prcnt_of_final!=1,], 
       aes(x = hours_since_collection,y = resp_prcnt_of_final)) + geom_point() +
  geom_smooth(method="loess") + geom_vline(xintercept=3)

#what is the minimum hours since collection?
min(respiration_exp$hours_since_collection)

#convert yes-no variables to 0-1
respiration_exp$visible_algae_lichen_num <- ifelse(respiration_exp$visible_algae_lichen=="yes",1,0)
respiration_exp$visible_hyphae_num <- ifelse(respiration_exp$visible_hyphae=="yes",1,0)
respiration_exp$termites_num <- ifelse(respiration_exp$termites=="yes",1,0)

#convert downup to bottom top
respiration_exp$downup <- gsub("up","top",respiration_exp$downup)
respiration_exp$downup <- gsub("down","bottom",respiration_exp$downup)

#fix commas in diameter values
str(respiration_exp)

respiration_exp$diameter <- as.numeric(as.character(gsub(",",".",respiration_exp$diameter)))
respiration_exp$wood_density <- as.numeric(as.character(gsub(",",".",respiration_exp$wood_density)))
respiration_exp$moisture_content <- as.numeric(as.character(gsub(",",".",respiration_exp$moisture_content)))
respiration_exp$distance_elevated <- as.numeric(as.character(gsub(",",".",respiration_exp$distance_elevated)))
respiration_exp$initial_weight <- as.numeric(as.character(gsub(",",".",respiration_exp$initial_weight)))
respiration_exp$final_weight <- as.numeric(as.character(gsub(",",".",respiration_exp$final_weight)))
                      
#add colum moisture content as (wet-dry)/dry * 100

respiration_exp$moist_cont_dry <- ((respiration_exp$initial_weight- respiration_exp$final_weight)/ respiration_exp$final_weight)*100

########## SAVED DATA FOR PAPER
#change working directory to save data
setwd('C:\\Users\\angel\\Desktop\\Docs_paper\\2_Wood_respiration_macroorgan_measurements')

#keep certain columns
colnames(respiration_exp)
respiration_exp_a <- respiration_exp[ , names(respiration_exp) %in% c("wood_ID",
                                                                      "treatment",
                                                                      "slope",                                                                        "decay_status",
                                                                      "termites_num",
                                                                      "visible_hyphae_num",
                                                                      "visible_algae_lichen_num",
                                                                      "moist_cont_dry",
                                                                      "moisture_content",
                                                                      "wood_density"
                                                                      )] 
#write .csv file
write.csv(respiration_exp_a,"respiration_model_results.csv")


##########

#average all variables for analyses

elevation_exp_avg <- respiration_exp %>%
  group_by(wood_ID,wood_sample,elevated_status,downup) %>%
  summarise(mean_resp_rate = mean(slope),
            wood_density = mean(wood_density),
            mean_temp = mean(mean_temp),
            moisture_content = mean(moisture_content),
            moist_cont_dry = mean(moist_cont_dry),
            distance_ground_connected = mean(distance_ground_connected,na.rm=T),
            distance_elevated = mean(distance_elevated,na.rm=T),
            termites = mean(termites_num),
            decay_status = mean(decay_status),
            diameter = mean(diameter),
            visible_algae_lichen = mean(visible_algae_lichen_num),
            hyphae = mean(visible_hyphae_num),
            initial_weight =  mean(initial_weight),
            final_weight =  mean(final_weight))
str(elevation_exp_avg)

#run models
#-- -- -- -- -- -- -- -- -- 
#### Moisture content #### 
#-- -- -- -- -- -- -- -- --

mod_moist1 <- lmer(moisture_content ~ elevated_status + downup + elevated_status:downup + (1|wood_ID),
                  data = elevation_exp_avg)
mod_moist2 <- lmer(moisture_content ~ elevated_status + downup + (1|wood_ID), # AIC 629.43
                  data = elevation_exp_avg)
mod_moist3 <- lmer(moisture_content ~ elevated_status + (1|wood_ID), # AIC 627.98
                   data = elevation_exp_avg)
mod_moist4 <- lmer(moisture_content ~  downup + (1|wood_ID),
                   data = elevation_exp_avg)
mod_moist5 <- lmer(moisture_content ~  1 + (1|wood_ID),
                   data = elevation_exp_avg)

#the analyses below test whether both models similarly explain the data
anova(mod_moist1,mod_moist2) 
#if the above analysis is significant, do not move forward
anova(mod_moist2,mod_moist3)
anova(mod_moist2,mod_moist4)  
#if both of the above are significant, then do not move forward
#if one of the two above are significant, then compare the model with that variable to the null model
anova(mod_moist3,mod_moist5) 
anova(mod_moist4,mod_moist5)

summary(mod_moist3)

#examine the residuals of the best fit model 
plot(mod_moist3) # trim_29_45
hist(resid(mod_moist3)) # trim_29_45

#Interpretation: Being elevated from the forest floor affects wood moisture content (mod_moist3)
# Being elevated from the forest floor and at the top/bottom sections affect wood moisture content (mod_moist3)

#plot the moisture content 
moisture_plot_trim_29_45 <- ggplot(elevation_exp_avg, aes(y = moisture_content, x = elevated_status)) +
  geom_boxplot(aes(fill = downup),color = "black") +
  theme_basis  + theme(legend.title=element_blank(),legend.background = element_blank(),
                       legend.key = element_blank(), legend.position = c(.15,.95),
                       legend.text = element_text(family = "Arial",colour="black", face ="bold", size = 10))+
  scale_x_discrete(name = element_blank(),labels=c("Elevated","Downed")) +
  scale_y_continuous(name = "Moisture content (%)",limits = c(10,72),breaks = seq(10,70,10),expand = c(0,0),
                     labels = c("10","","30","","50","","70")) +
  scale_fill_manual(values=c("grey40","white"))
  
moisture_plot_trim_29_45

#set directory to save figure 
setwd('C:/Users/angel/OneDrive/2_respiration_activity_20220816/figures/respiration_experiment')
ggsave("moisture_plot_trim_29_45.tiff",moisture_plot_trim_29_45,dpi = 600, compression = "lzw",width = 3.25,height = 2,scale = 1.8)

#-- -- -- -- -- -- -- -- -- 
#### Moisture content DRY #### 
#-- -- -- -- -- -- -- -- --

mod_moist1 <- lmer(moist_cont_dry ~ elevated_status + downup + elevated_status:downup + (1|wood_ID),
                   data = elevation_exp_avg)
mod_moist2 <- lmer(moist_cont_dry ~ elevated_status + downup + (1|wood_ID), # AIC 629.43
                   data = elevation_exp_avg)
mod_moist3 <- lmer(moist_cont_dry ~ elevated_status + (1|wood_ID), # AIC 627.98
                   data = elevation_exp_avg)
mod_moist4 <- lmer(moist_cont_dry ~  downup + (1|wood_ID),
                   data = elevation_exp_avg)
mod_moist5 <- lmer(moist_cont_dry ~  1 + (1|wood_ID),
                   data = elevation_exp_avg)

#the analyses below test whether both models similarly explain the data
anova(mod_moist1,mod_moist2) 
#if the above analysis is significant, do not move forward
anova(mod_moist2,mod_moist3)
anova(mod_moist2,mod_moist4)  
#if both of the above are significant, then do not move forward
#if one of the two above are significant, then compare the model with that variable to the null model
anova(mod_moist3,mod_moist5) 
anova(mod_moist4,mod_moist5)

summary(mod_moist3)

#examine the residuals of the best fit model 
plot(mod_moist3) # trim_29_45
hist(resid(mod_moist3)) # trim_29_45

#Interpretation: Being elevated from the forest floor affects wood moisture content (mod_moist3)
# Being elevated from the forest floor and at the top/bottom sections affect wood moisture content (mod_moist3)

#plot the moisture content 
moisture_plot_trim_29_45 <- ggplot(elevation_exp_avg, aes(y = moisture_content, x = elevated_status)) +
  geom_boxplot(aes(fill = downup),color = "black") +
  theme_basis  + theme(legend.title=element_blank(),legend.background = element_blank(),
                       legend.key = element_blank(), legend.position = c(.15,.95),
                       legend.text = element_text(family = "Arial",colour="black", face ="bold", size = 10))+
  scale_x_discrete(name = element_blank(),labels=c("Elevated","Downed")) +
  scale_y_continuous(name = "Moisture content (%)",limits = c(10,72),breaks = seq(10,70,10),expand = c(0,0),
                     labels = c("10","","30","","50","","70")) +
  scale_fill_manual(values=c("grey40","white"))

moisture_plot_trim_29_45

#set directory to save figure 
setwd('C:/Users/angel/OneDrive/2_respiration_activity_20220816/figures/respiration_experiment')
ggsave("moisture_plot_trim_29_45.tiff",moisture_plot_trim_29_45,dpi = 600, compression = "lzw",width = 3.25,height = 2,scale = 1.8)





#-- -- -- -- -- -- -- 
#### Wood density #### 
#-- -- -- -- -- -- -- 
str(elevation_exp_avg$wood_density)

mod_dens1 <- lmer(wood_density ~ elevated_status + downup + elevated_status:downup + (1|wood_ID),
                   data = elevation_exp_avg)
mod_dens2 <- lmer(wood_density ~ elevated_status + downup + (1|wood_ID),
                   data = elevation_exp_avg)
mod_dens3 <- lmer(wood_density ~ elevated_status + (1|wood_ID),
                   data = elevation_exp_avg)
mod_dens4 <- lmer(wood_density ~  downup + (1|wood_ID),
                   data = elevation_exp_avg)
mod_dens5 <- lmer(wood_density ~  1 + (1|wood_ID),
                   data = elevation_exp_avg)

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
plot(mod_dens3) # trim_29_45
hist(resid(mod_dens3)) # trim_29_45

#plot the moisture content 
density_plot_trim_29_45 <- ggplot(elevation_exp_avg, aes(y = wood_density, x = elevated_status)) +
  geom_boxplot(aes(fill = downup),color = "black") +
  theme_basis  + theme(legend.title=element_blank(),legend.background = element_blank(),
                       legend.key = element_blank(), legend.position = c(.15,.95),
                       legend.text = element_text(family = "Arial",colour="black", face ="bold", size = 10))+ 
  scale_x_discrete(name = element_blank(),labels=c("Elevated","Downed")) +
 scale_y_continuous(name = expression(bold("Wood density"~(g~mm^{-3}))),
                    limits = c(.15,.65),breaks = seq(0.2,0.6,.1),expand = c(0,0),
                     labels = c("0.2","","0.4","","0.6")) + scale_fill_manual(values=c("grey40","white"))

density_plot_trim_29_45  # trim_29_45

#save figure
setwd('C:/Users/angel/OneDrive/2_respiration_activity_20220816/figures/respiration_experiment')
ggsave("density_plot_trim_29_45.tiff",density_plot_trim_29_45,dpi = 600, compression = "lzw",width = 3.25,height = 2,scale = 1.8)

#-- -- -- -- -- -- -- -- -- -- -- 
#### Respiration rate per mass #### 
#-- -- -- -- -- -- -- -- -- -- -- 

#remove non significant respiration trends
respiration_exp_sigonly <- respiration_exp[respiration_exp$pval<0.05,]

#average all variables for analyses
respiration_exp_avg <- respiration_exp_sigonly %>%
  group_by(wood_ID,wood_sample,elevated_status,downup) %>%
  summarise(mean_resp_rate = mean(slope),
            wood_density = mean(wood_density),
            mean_temp = mean(mean_temp),
            moisture_content = mean(moisture_content),
            distance_ground_connected = mean(distance_ground_connected,na.rm=T),
            distance_elevated = mean(distance_elevated,na.rm=T),
            termites = mean(termites_num),
            decay_status = mean(decay_status),
            diameter = mean(diameter),
            visible_algae_lichen = mean(visible_algae_lichen_num),
            hyphae = mean(visible_hyphae_num),
            initial_weight =  mean(initial_weight),
            final_weight =  mean(final_weight))
str(respiration_exp_avg)
# create new column for respiration rate per mass 
respiration_exp_avg$resp_mass <- respiration_exp_avg$mean_resp_rate/respiration_exp_avg$final_weight

#Do all variables have multiple samples?
table(respiration_exp_avg$wood_ID)

#Everything  
#add temperature as predictor 
mod_resp_mass1 <- lmer(resp_mass ~ elevated_status + downup + elevated_status:downup + mean_temp + (1|wood_ID),
                  data = respiration_exp_avg)
mod_resp_mass2 <- lmer(resp_mass ~ elevated_status + downup + mean_temp +(1|wood_ID),
                  data = respiration_exp_avg)
mod_resp_mass3 <- lmer(resp_mass ~ elevated_status + mean_temp + (1|wood_ID),
                  data = respiration_exp_avg)
mod_resp_mass4 <- lmer(resp_mass ~  downup + mean_temp +(1|wood_ID),
                  data = respiration_exp_avg)
mod_resp_mass5 <- lmer(resp_mass ~  mean_temp + (1|wood_ID),
                  data = respiration_exp_avg)

#the analyses below test whether both models similarly explain the data
anova(mod_resp_mass1,mod_resp_mass2) 
#if the above analysis is significant, do not move forward
anova(mod_resp_mass2,mod_resp_mass3) 
anova(mod_resp_mass2,mod_resp_mass4)

#if both of the above are significant, then do not move forward
#if one of the two above are significant, then compare the model with that variable to the null model
anova(mod_resp_mass3,mod_resp_mass5) 
anova(mod_resp_mass4,mod_resp_mass5) 

#examine the residuals of the best fit model (#mod_resp_mass1, 0.09989 .)
plot(mod_resp_mass1) # trim_29_45
qqnorm(resid(mod_resp_mass1));qqline(resid(mod_resp_mass1)) # trim_29_35
hist(resid(mod_resp_mass1)) # trim_29_35
# respiration_exp_avg[resid(mod_resp_mass1)>3e-05,c("wood_sample")]
# View(respiration_exp_avg)
summary(mod_resp_mass1)

#NO SIGNIFICANT DIFFERENCES

#plot

resp_plot_trim_29_45 <- ggplot(respiration_exp_avg, aes(y = resp_mass, x = elevated_status)) +
  geom_boxplot(aes(fill = downup),color = "black") +
  theme_basis  + theme(legend.title=element_blank(),legend.background = element_blank(),
                       legend.key = element_blank(), legend.position = c(.15,.95),
                       legend.text = element_text(family = "Arial",colour="black", face ="bold", size = 10))+ 
  scale_x_discrete(name = element_blank(),labels=c("Elevated","Downed")) +
  scale_fill_manual(values=c("grey40","white")) +
  scale_y_continuous(name = expression(bold('Respiration rate'~(mg[CO[2]]~'s'^-1~'g'[wood]^-1))))
# limits = c(.15,.65),breaks = seq(0.2,0.6,.1),expand = c(0,0),
# labels = c("0.2","","0.4","","0.6"))

resp_plot_trim_29_45

#save figure
setwd('C:/Users/angel/OneDrive/2_respiration_activity_20220816/figures/respiration_experiment')
ggsave("resp_plot_trim_29_45_20220909.tiff",resp_plot_trim_29_45,dpi = 600, compression = "lzw",width = 3.25,height = 2,scale = 1.8)

#______________________________
####  *Log transformation of resp_mass (WE USED THIS ONE)####

respiration_exp_avg$resp_masslog <- log(respiration_exp_avg$resp_mass)
colnames(respiration_exp_avg)

mod_resp_masslog1 <- lmer(resp_masslog ~ elevated_status + downup + elevated_status:downup + mean_temp + (1|wood_ID),
                       data = respiration_exp_avg)
mod_resp_masslog_notemp <- lmer(resp_masslog ~ elevated_status + downup + elevated_status:downup + (1|wood_ID),
                          data = respiration_exp_avg)
mod_resp_masslog2 <- lmer(resp_masslog ~ elevated_status + downup + mean_temp +(1|wood_ID),
                       data = respiration_exp_avg)
mod_resp_masslog3 <- lmer(resp_masslog ~ elevated_status + mean_temp + (1|wood_ID),
                       data = respiration_exp_avg)
mod_resp_masslog4 <- lmer(resp_masslog ~  downup + mean_temp +(1|wood_ID),
                       data = respiration_exp_avg)
mod_resp_masslog5 <- lmer(resp_masslog ~  mean_temp + (1|wood_ID),
                       data = respiration_exp_avg)

#the analyses below test whether both models similarly explain the data
anova(mod_resp_masslog1, mod_resp_masslog_notemp) 
anova(mod_resp_masslog1,mod_resp_masslog2) 
#if the above analysis is significant, do not move forward
anova(mod_resp_masslog2,mod_resp_masslog3) 
anova(mod_resp_masslog2,mod_resp_masslog4)

#if both of the above are significant, then do not move forward
#if one of the two above are significant, then compare the model with that variable to the null model
anova(mod_resp_masslog4,mod_resp_masslog5) 

#examine the residuals of the best fit model (#mod_resp_masslog1, 0.09989 .)
plot(mod_resp_masslog4) # _trim_29_45
qqnorm(resid(mod_resp_masslog4));qqline(resid(mod_resp_masslog4)) # _trim_29_45
hist(resid(mod_resp_masslog4)) # _trim_29_45

# respiration_exp_avg[resid(mod_resp_masslog1)>3e-05,c("wood_sample")]
# View(respiration_exp_avg)
summary(mod_resp_masslog4)# _trim_29_45

#plot (!!!!!!!!!!!!)

resp_plot_log_sus_trim_29_45 <- ggplot(respiration_exp_avg, aes(y = resp_mass, x = elevated_status)) +
  geom_boxplot(aes(fill = downup),color = "black") +
  theme_basis  + theme(legend.title=element_blank(),legend.background = element_blank(),
                       legend.key = element_blank(), legend.position = c(.15,.95),
                       legend.text = element_text(family = "Arial",colour="black", face ="bold", size = 10))+ 
  scale_x_discrete(name = element_blank(),labels=c("Suspended","Downed")) +
  scale_fill_manual(values=c("grey40","white"), labels=c("Bottom","Top")) +
  scale_y_log10(name = expression(bold('Respiration rate'~(mg[CO[2]]~'s'^-1~'g'[wood]^-1))),
                breaks = c(seq(0.000001,0.00001,0.000001),seq(0.00002,0.0001,0.00001)), 
                labels = c("0.000001",rep("",8),"0.00001",rep("",8),"0.0001"),
                limits = c(0.000001,0.00015),
                expand = c(0,0))


resp_plot_log_sus_trim_29_45

#save figure
setwd('C:/Users/angel/OneDrive/2_respiration_activity_20220816/figures/respiration_experiment')
ggsave("resp_mass_plot_log_sus_trim_29_45_20221025.tiff",resp_plot_log_trim_29_45,dpi = 600, compression = "lzw",width = 3.25,height = 2,scale = 1.8)

resp_plot_logg_sus_trim_29_45 <- ggplot(respiration_exp_avg, aes(y = resp_masslog, x = elevated_status)) +
  geom_boxplot(aes(fill = downup),color = "black") +
  theme_basis  + theme(legend.title=element_blank(),legend.background = element_blank(),
                       legend.key = element_blank(), legend.position = c(.15,.95),
                       legend.text = element_text(family = "Arial",colour="black", face ="bold", size = 10))+ 
  scale_x_discrete(name = element_blank(),labels=c("Suspended","Downed")) +
  scale_fill_manual(values=c("grey40","white"), labels=c("Bottom","Top")) +
  scale_y_log10(name = expression(bold('Respiration rate'~(mg[CO[2]]~'s'^-1~'g'[wood]^-1))),
                breaks = c(seq(0.000001,0.00001,0.000001),seq(0.00002,0.0001,0.00001)), 
                labels = c("0.000001",rep("",8),"0.00001",rep("",8),"0.0001"),
                limits = c(0.000001,0.00015),
                expand = c(0,0))


resp_plot_logg_sus_trim_29_45





#-- -- -- -- -- -- -- -- -- -- -- -- 
#summary sd se mean ALL ("observer")
rr_se_sd_mean_all <- ddply(respiration_exp_avg, .(observer), summarise, 
                            M = mean(moist_cont_dry), SE = sd(moist_cont_dry) / sqrt((length(moist_cont_dry))), 
                            SD = sd(moist_cont_dry), Median=median(moist_cont_dry))



#summary sd se mean by elevated_status
rr_se_sd_mean_elev <- ddply(respiration_exp_avg, .(elevated_status), summarise, 
                             M = mean(resp_mass), SE = sd(resp_mass) / sqrt((length(resp_mass))), 
                             SD = sd(resp_mass), Median=median(resp_mass))

#summary sd se mean by elevated_status and downup
rr_se_sd_mean_elev_updown <- ddply(respiration_exp_avg, .(elevated_status, downup), summarise, 
                                    M = mean(resp_mass), SE = sd(resp_mass) / sqrt((length(resp_mass))), 
                                    SD = sd(resp_mass), Median=median(resp_mass))





#-- -- -- -- -- -- -- -- -- -- -- --
####  *Highest r2 values          ####
#-- -- -- -- -- -- -- -- -- -- -- --
#Run the respiration models only including the respiration measurement with the 
#highest r2 rather than the average slope

#highest R2 values (we use respiration_exp_sigonly because it has all the respiration
# measurements with the same wood sample, not the mean of them)
respiration_exp_highest_r2 <- respiration_exp_avg %>% 
  group_by(wood_sample) %>% 
  arrange(desc(r2), .by_group =T) %>% 
  slice(1)

#run model Respiration rate per mass
mod_resp_mass1_r2 <- lmer(resp_mass ~ elevated_status + downup + elevated_status:downup + mean_temp + (1|wood_ID),
                          data = respiration_exp_highest_r2)
mod_resp_mass2_r2 <- lmer(resp_mass ~ elevated_status + downup + mean_temp +(1|wood_ID),
                          data = respiration_exp_highest_r2)
mod_resp_mass3_r2 <- lmer(resp_mass ~ elevated_status + mean_temp + (1|wood_ID),
                          data = respiration_exp_highest_r2)
mod_resp_mass4_r2 <- lmer(resp_mass ~  downup + mean_temp +(1|wood_ID),
                          data = respiration_exp_highest_r2)
mod_resp_mass5_r2 <- lmer(resp_mass ~  mean_temp + (1|wood_ID),
                          data = respiration_exp_highest_r2)

#the analyses below test whether both models similarly explain the data
anova(mod_resp_mass1_r2,mod_resp_mass2_r2) # 
#if the above analysis is significant, do not move forward
anova(mod_resp_mass1_r2,mod_resp_mass3_r2) #
anova(mod_resp_mass1_r2,mod_resp_mass4_r2) #

#if both of the above are significant, then do not move forward
#if one of the two above are significant, then compare the model with that variable to the null model
anova(mod_resp_mass3_r2,mod_resp_mass5_r2)
anova(mod_resp_mass4_r2,mod_resp_mass5_r2) #


#examine the residuals of the best fit model (#mod_resp_mass1, 0.09989 .)
plot(mod_resp_mass1_r2) # _trim_29_45 
qqnorm(resid(mod_resp_mass1_r2));qqline(resid(mod_resp_mass1_r2)) # _trim_29_45
hist(resid(mod_resp_mass1_r2))# _trim_29_45
# respiration_exp_avg[resid(mod_resp_mass1)>3e-05,c("wood_sample")]
# View(respiration_exp_avg)
summary(mod_resp_mass1_r2)# _trim_29_45

#plot

resp_plot_r2_trim_29_45 <- ggplot(respiration_exp_highest_r2, aes(y = resp_mass, x = elevated_status)) +
  geom_boxplot(aes(fill = downup),color = "black") +
  theme_basis  + theme(legend.title=element_blank(),legend.background = element_blank(),
                       legend.key = element_blank(), legend.position = c(.15,.95),
                       legend.text = element_text(family = "Arial",colour="black", face ="bold", size = 10))+ 
  scale_x_discrete(name = element_blank(),labels=c("Elevated","Downed")) +
  scale_fill_manual(values=c("grey40","white")) +
  scale_y_continuous(name = expression(bold('Respiration rate'~(mg[CO[2]]~'s'^-1~'g'[wood]^-1))))
# limits = c(.15,.65),breaks = seq(0.2,0.6,.1),expand = c(0,0),
# labels = c("0.2","","0.4","","0.6"))

resp_plot_r2_trim_29_45

#save figure
setwd('C:/Users/angel/OneDrive/2_respiration_activity_20220816/figures/respiration_experiment')
ggsave("resp_mass_plot_highestr2_trim_29_45.tiff",resp_plot_r2_trim_29_45,dpi = 600, compression = "lzw",width = 3.25,height = 2,scale = 1.8)


#______________________________
####  *Log transformation of respiration_exp_highest_r2$resp_mass####

mod_resp_masslog1_r2 <- lmer(resp_masslog ~ elevated_status + downup + elevated_status:downup + mean_temp + (1|wood_ID),
                          data = respiration_exp_highest_r2)
mod_resp_masslog_notemp_r2 <- lmer(resp_masslog ~ elevated_status + downup + elevated_status:downup + (1|wood_ID),
                                data = respiration_exp_highest_r2)
mod_resp_masslog2_r2 <- lmer(resp_masslog ~ elevated_status + downup + mean_temp +(1|wood_ID),
                          data = respiration_exp_highest_r2)
mod_resp_masslog3_r2 <- lmer(resp_masslog ~ elevated_status + mean_temp + (1|wood_ID),
                          data = respiration_exp_highest_r2)
mod_resp_masslog4_r2 <- lmer(resp_masslog ~  downup + mean_temp +(1|wood_ID),
                          data = respiration_exp_highest_r2)
mod_resp_masslog5_r2 <- lmer(resp_masslog ~  mean_temp + (1|wood_ID),
                          data = respiration_exp_highest_r2)

#the analyses below test whether both models similarly explain the data
anova(mod_resp_masslog1_r2, mod_resp_masslog_notemp_r2) #
anova(mod_resp_masslog_notemp_r2,mod_resp_masslog2_r2) #
#if the above analysis is significant, do not move forward
anova(mod_resp_masslog2_r2,mod_resp_masslog3_r2) #
anova(mod_resp_masslog2_r2,mod_resp_masslog4_r2) #

#if both of the above are significant, then do not move forward
#if one of the two above are significant, then compare the model with that variable to the null model
anova(mod_resp_masslog3_r2,mod_resp_masslog5_r2) #
anova(mod_resp_masslog4_r2,mod_resp_masslog5_r2) #

#examine the residuals of the best fit model (#mod_resp_masslog1, 0.09989 .)
plot(mod_resp_masslog3_r2) # _trim_29_45
qqnorm(resid(mod_resp_masslog3_r2));qqline(resid(mod_resp_masslog3_r2)) # _trim_29_35
hist(resid(mod_resp_masslog3_r2))# _trim_29_45
# respiration_exp_avg[resid(mod_resp_masslog1)>3e-05,c("wood_sample")]
# View(respiration_exp_avg)
summary(mod_resp_masslog3_r2)

#plot

resp_plot__log_r2_trim_29_45 <- ggplot(respiration_exp_highest_r2, aes(y = resp_masslog, x = elevated_status)) +
  geom_boxplot(aes(fill = downup),color = "black") +
  theme_basis  + theme(legend.title=element_blank(),legend.background = element_blank(),
                       legend.key = element_blank(), legend.position = c(.15,.95),
                       legend.text = element_text(family = "Arial",colour="black", face ="bold", size = 10))+ 
  scale_x_discrete(name = element_blank(),labels=c("Elevated","Downed")) +
  scale_fill_manual(values=c("grey40","white")) +
  scale_y_continuous(name = expression(bold('log(Respiration rate)'~(mg[CO[2]]~'s'^-1~'g'[wood]^-1))))
# limits = c(.15,.65),breaks = seq(0.2,0.6,.1),expand = c(0,0),
# labels = c("0.2","","0.4","","0.6"))

resp_plot__log_r2_trim_29_45

#save figure
setwd('C:/Users/angel/OneDrive/2_respiration_activity_20220816/figures/respiration_experiment')
ggsave("resp_plot__log_r2_trim_29_45.tiff",resp_plot__log_r2_trim_29_45,dpi = 600, compression = "lzw",width = 3.25,height = 2,scale = 1.8)


#-- -- -- -- -- -- -- -- -- -- -- --
####  *Respiration rate on diff. decomp. status ####
#-- -- -- -- -- -- -- -- -- -- -- --
# Test how respiration rate differs among different decomposition status
#data = respiration_exp_highest_r2
table(respiration_exp_highest_r2$wood_ID,respiration_exp_highest_r2$decay_status)

#plot 
resp_decomp_r2_tag_trim_29_45 <- ggplot(respiration_exp_highest_r2, aes(y = resp_mass, x = elevated_status)) +
  geom_boxplot(aes(fill = downup),color = "black") +
  theme_basis  + theme(legend.title=element_blank(),legend.background = element_blank(),
                       legend.key = element_blank(), legend.position = c(.15,.95),
                       legend.text = element_text(family = "Arial",colour="black", face ="bold", size = 10))+ 
  scale_x_discrete(name = element_blank(),labels=c("Elevated","Downed")) +
  scale_fill_manual(values=c("grey40","white")) +
  scale_y_continuous(name = expression(bold('Respiration rate'~(mg[CO[2]]~'s'^-1~'g'[wood]^-1))))+
  facet_wrap(~decay_status,ncol=3)

resp_decomp_r2_tag_trim_29_45

#save figure
setwd('C:/Users/angel/OneDrive/2_respiration_activity_20220816/figures/respiration_experiment')
ggsave("resp_decomp_r2_tag_trim_29_45.tiff",resp_decomp_r2_tag_trim_29_45,dpi = 600, compression = "lzw",width = 3.25,height = 2,scale = 1.8)


#____
####  *Run model Respiration rate per mass driven by decomp ####
mod_decomp1_r2 <- lmer(resp_masslog ~ decay_status + elevated_status + downup + elevated_status:downup + mean_temp + (1|wood_ID),
                          data = respiration_exp_highest_r2)
mod_decomp2_r2 <- lmer(resp_masslog ~ decay_status + elevated_status + downup + mean_temp +(1|wood_ID),
                          data = respiration_exp_highest_r2)
mod_decomp3_r2 <- lmer(resp_masslog ~ decay_status + elevated_status + mean_temp + (1|wood_ID),
                          data = respiration_exp_highest_r2)
mod_decomp4_r2 <- lmer(resp_masslog ~  decay_status + downup + mean_temp +(1|wood_ID),
                          data = respiration_exp_highest_r2)
mod_decomp5_r2 <- lmer(resp_masslog ~  decay_status + mean_temp + (1|wood_ID),
                          data = respiration_exp_highest_r2)
mod_decomp6_r2 <- lmer(resp_masslog ~  mean_temp + (1|wood_ID),
                       data = respiration_exp_highest_r2)

#the analyses below test whether both models similarly explain the data
anova(mod_decomp1_r2,mod_decomp2_r2) # 
#if the above analysis is significant, do not move forward
anova(mod_decomp1_r2,mod_decomp3_r2) # 
anova(mod_decomp1_r2,mod_decomp4_r2) # 

#if both of the above are significant, then do not move forward
#if one of the two above are significant, then compare the model with that variable to the null model
anova(mod_decomp1_r2,mod_decomp5_r2) # 

#examine the residuals of the best fit model 
plot(mod_decomp1_r2) # _trim_29_45
qqnorm(resid(mod_decomp1_r2));qqline(resid(mod_decomp1_r2)) # _trim_29_45
hist(resid(mod_decomp1_r2))# _trim_29_45
# respiration_exp_avg[resid(mod_resp_masslog1)>3e-05,c("wood_sample")]
# View(respiration_exp_avg)
summary(mod_decomp1_r2)



####  Evaluate whether a single factor fits the data better  ####

respiration_exp_avg$treatment <- paste(respiration_exp_avg$elevated_status,respiration_exp_avg$downup,sep="_")

mod_resp_masslog1 <- lmer(resp_masslog ~ treatment + mean_temp + decay_status +
                            diameter + (1|wood_ID),
                          data = respiration_exp_avg)
mod_resp_masslog_notemp <- lmer(resp_masslog ~ treatment + decay_status +
                                  diameter + (1|wood_ID),
                                data = respiration_exp_avg)
mod_resp_masslog2 <- lmer(resp_masslog ~ 1  + decay_status +
                            diameter + mean_temp +(1|wood_ID),
                          data = respiration_exp_avg)


#the analyses below test whether both models similarly explain the data
anova(mod_resp_masslog1, mod_resp_masslog_notemp) 
anova(mod_resp_masslog1,mod_resp_masslog2) 


#calculate corrected distance elevated
respiration_exp_avg$corrected_distance_elevated <- ifelse(respiration_exp_avg$downup=="bottom",
                                                          respiration_exp_avg$distance_elevated,
                                                          respiration_exp_avg$distance_elevated+respiration_exp_avg$diameter/2)


mod_resp_masslog1 <- lmer(resp_masslog ~ corrected_distance_elevated + mean_temp + decay_status +
                            diameter + (1|wood_ID),
                          data = respiration_exp_avg)
mod_resp_masslog_notemp <- lmer(resp_masslog ~ corrected_distance_elevated + decay_status +
                                  diameter + (1|wood_ID),
                                data = respiration_exp_avg)
mod_resp_masslog2 <- lmer(resp_masslog ~ 1  + decay_status +
                            diameter + mean_temp +(1|wood_ID),
                          data = respiration_exp_avg)

#the analyses below test whether both models similarly explain the data
anova(mod_resp_masslog1, mod_resp_masslog_notemp) 
anova(mod_resp_masslog1,mod_resp_masslog2) 

summary(mod_resp_masslog1)

#look at residuals
plot(mod_resp_masslog1)
qqnorm(resid(mod_resp_masslog1));qqline(resid(mod_resp_masslog1))
