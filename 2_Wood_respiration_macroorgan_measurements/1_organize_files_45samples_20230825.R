#Feb_02_2022 
#Edited by Angela Marcela Barrera Bello and Evan Gora

#Create a script for processing Vaisala GMP434 data to calculate respiration rate
#first import all datasets and convert to CO2 concetrations (mg per m3)
#the create figures for each sample and calculate respiration rate metrics


#load packages
library(ggplot2)
library(extrafont)
library(remotes)
library(lme4)
library(readr)
library(dplyr)


#remotes::install_version("Rttf2pt1", version = "1.3.8")
extrafont::font_import()
# font_import()
loadfonts(device="win")

#create theme basis for ggplot
theme_base<-theme(axis.title = element_text(family = "Arial", color="black", face="bold", size=14))+
  theme(axis.text=element_text(family = "Arial",colour="black", face ="plain", size = 12)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"),axis.ticks = element_line(colour="black"))


#-- -- -- -- -- -- -- -- -- -- -- --
#### Import data #### 
#-- -- -- -- -- -- -- -- -- -- -- --

#set directory
setwd('C:\\Users\\angel\\OneDrive\\Independent_project_wood_20230823\\2_respiration_activity_20220816\\vaisala_files_txt_new')

#identify file names
temp = list.files(pattern="*.txt")
myfiles = lapply(temp, FUN=read.table,skip =1)
str(myfiles)#data added without column names

#create a site-names object
site_names <- gsub(".txt","",temp)

#fix columns and combine data in a single for loop
for(i in 1:length(temp)) {
  temp_df <- as.data.frame(myfiles[i])
  if(ncol(temp_df)==6){
    temp_df$V3<- paste(temp_df$V3,temp_df$V4,sep="")
    temp_df <- temp_df[,-4]
  }
  temp_df$time_series <- seq(5,300,5)
  temp_df$sample <- site_names[i]
  colnames(temp_df)[1:5] <- c("date","hour", "am_pm", "co2","temp_c")
  temp_df$co2 <- gsub(",",".",temp_df$co2)
  temp_df$temp_c <- gsub(",",".",temp_df$temp_c)
  if(i == 1){ 
    comb_df <- temp_df}else{
      comb_df <- rbind.data.frame(comb_df,temp_df)}  
}

summary(comb_df)

comb_df$co2 <- as.numeric(comb_df$co2)
comb_df$temp_c <- as.numeric(comb_df$temp_c)

comb_df_trim <- comb_df[!is.na(comb_df$co2),]
#are there major NA issues?
summary(is.na(comb_df_trim))# no NAs 


#-- -- -- -- -- -- -- -- -- -- -- --
####  Calculate respiration rate for each file ####
#-- -- -- -- -- -- -- -- -- -- -- --

#FORMULA TO CONVERT PPM TO MG/M3
# concentration Co2 (mg/m3) = 0.0409 x concentration x molecular_weight
#m3=vol_chamber_m3; 0.0409=gas_constant; concentration=ppm; molecular_weight=co2_mol_weight

#We know: gas_constant, ppm,co2_mol_weight, vol_chamber_m3

gas_constant <- 0.0409 #constant
ppm <- comb_df_trim$co2 #ppm of c02
co2_mol_weight <- 44.01 # g/mol of co2
vol_chamber_m3 <- 0.004164 #m3

#Find Co2 gr. using formula above:
comb_df_trim$mg_co2 <- gas_constant*ppm*co2_mol_weight*vol_chamber_m3

#create a vector for our sample
sample_names <- unique(comb_df_trim$sample)

#change working directory to save data
setwd('C:\\Users\\angel\\Desktop\\Docs_paper\\2_Wood_respiration_macroorgan_measurements')

#now create for loop for plotting the data and calculating respiration rate
for(i in 1:length(sample_names)){
  working_sample <- sample_names[i]
  sample_df <- comb_df_trim[comb_df_trim$sample==working_sample,]
  sample_df_trim <- sample_df[sample_df$time_series>45,]
  
#plot the data
ggobj <- ggplot(sample_df_trim,aes(x = time_series, y = mg_co2)) +
  geom_point() + geom_smooth(method = "lm") + theme_base 

ggobj + 
  xlab("Time (seconds)") + 
  ylab(expression(bold('CO'[2]~(mg))))


filename <- paste(working_sample,"respiration.png",sep ="_")

ggsave(filename,ggobj,dpi = 150,width = 3.25, height = 2.5,scale = 1.8)

#model the data
mod1 <- lm(mg_co2 ~ time_series,data = sample_df_trim)
model_sum <- summary(mod1)
slope <- model_sum$coefficients[2,1]
slope_stder <- model_sum$coefficients[2,2]
r2 <- model_sum$adj.r.squared
Fstat <- model_sum$fstatistic[1]
numdf <- model_sum$fstatistic[2]
dendf <- model_sum$fstatistic[3]
pval <- pf(Fstat,numdf,dendf,lower.tail=F)
RMSE <- sqrt(mean(mod1$residuals^2))

#add average temperature
mean_temp <- round(mean(sample_df_trim$temp_c),digits = 2)

#add time and date
am_pm <- sample_df$am_pm[1]
hour <- sample_df$hour[1]

model_vec <- c(working_sample,slope,slope_stder,r2,Fstat,
               numdf,dendf,pval,RMSE,mean_temp,am_pm,hour)

#save model metrics (slope, RMSE, R2, mean temperature during the period)
if(i ==1){
  model_results <- data.frame("sample" = NA,
                              "slope" = NA,
                              "slope_stder" = NA,
                              "r2" = NA,
                              "Fstat" = NA,
                              "numdf" = NA,
                              "dendf" = NA,
                              "pval" = NA,
                              "RMSE" = NA,
                              "mean_temp"= NA,
                              "am_pm" = NA,
                              "hour" = NA)
  model_results[1,]<- model_vec
}else{
  model_results[i,]<- model_vec
}

}

#create grouping variables
model_results$elevated_status <- ifelse(grepl("e",model_results$sample),"elevated","ground")
model_results$downup <- ifelse(grepl("d",model_results$sample),"down","up")

#separate sample name from re-sampling iteration
model_results$wood_sample<-NA
model_results$repetition <-NA
library(stringr)
for(i in 1:nrow(model_results)){
  model_results$wood_sample[i] <- str_split(model_results$sample,"_")[[i]][1]
  if(grepl("_",model_results$sample[i])){
    model_results$repetition[i] <- str_split(model_results$sample,"_")[[i]][2]
  }else{
    model_results$repetition[i] <- 0
  }

}

#examine whether the lines represent the data
model_results$pval<-as.numeric(model_results$pval)
hist(model_results$pval)
summary(model_results$pval<0.05)

#does each sample have a good measurement?
model_results$wood_sample <- gsub("22 ","22",model_results$wood_sample)
model_results$wood_ID <- gsub("ed|eu|gd|gu","",model_results$wood_sample)
model_results_onlysignificant <- model_results[model_results$pval<0.05,]
model_results_onlysignificant$treatment <- ifelse(grepl("ed",model_results_onlysignificant$wood_sample),"ed",
                                                  ifelse(grepl("eu",model_results_onlysignificant$wood_sample),"eu",
                                                                        ifelse(grepl("gd",model_results_onlysignificant$wood_sample),"gd","gu")))

#Run the following 2 lines to check which measurements have p value > 0.05
#If p value > 0.05 and figures (output) are not ok then remeasure 
table(model_results_onlysignificant$wood_ID,model_results_onlysignificant$treatment)

#look at the bad data
model_results$sample[model_results$pval>=0.05]
#Run until here to check significant values 


#drop flawed respiration measurements
#If p value < 0.05, but figure (output) is not ok, then remove and remeasure
model_results <- model_results[model_results$sample!="ed18_1",]
model_results <- model_results[model_results$sample!="ed19_1",]
model_results <- model_results[model_results$sample!="ed19_3",]
model_results <- model_results[model_results$sample!="ed20",]
model_results <- model_results[model_results$sample!="ed20_1",]
model_results <- model_results[model_results$sample!="ed21_1",]
model_results <- model_results[model_results$sample!="ed21_2",]
model_results <- model_results[model_results$sample!="ed21_3",]
model_results <- model_results[model_results$sample!="ed26_2",]
# model_results <- model_results[model_results$sample!="ed29",]
model_results <- model_results[model_results$sample!="eu13_1",]
model_results <- model_results[model_results$sample!="eu14_2",]
model_results <- model_results[model_results$sample!="eu15_2",]
model_results <- model_results[model_results$sample!="eu15_3",]
model_results <- model_results[model_results$sample!="eu19_3",]
model_results <- model_results[model_results$sample!="eu20_2",]
model_results <- model_results[model_results$sample!="eu21_1",]#the first 100 seconds were bad
model_results <- model_results[model_results$sample!="eu24_1",]
model_results <- model_results[model_results$sample!="gd16_3",]
model_results <- model_results[model_results$sample!="gd18_3",]
model_results <- model_results[model_results$sample!="gd20",]
#sample gd15_2 should be analyzed above 90 seconds only
model_results <- model_results[model_results$sample!="gd32",]
model_results <- model_results[model_results$sample!="gu16_3",]
model_results <- model_results[model_results$sample!="gu28_2",]
#All of the above is checked with Evan...

# #### select lowest p value per $woodsample ####
# 
# model_results_onlysignificant_low_pvalue <- model_results_onlysignificant %>%
#   group_by(wood_sample) %>%
#   slice_min(n = 1, pval)

#change working directory to save data
 setwd('C:\\Users\\angel\\Desktop\\Docs_paper\\2_Wood_respiration_macroorgan_measurements')

#keep certain columns
colnames(model_results_onlysignificant)
model_results_onlysignificant_a <- model_results_onlysignificant[ , names(model_results_onlysignificant) %in% c("sample","slope","elevated_status","downup")] 
names(model_results_onlysignificant_a)[names(model_results_onlysignificant_a) == 'elevated_status'] <- 'suspended_status'

#write .csv file
write.csv(model_results_onlysignificant_a,"respiration_model_results.csv")
view(model_results)

test <- gsub("p.<a0>m.","PM",model_results$am_pm)
table(test)

test2 <- gsub("a.<a0>m.","AM",model_results$am_pm)
table(test2)

model_results$am_pm <- gsub("p.<a0>m.","PM",model_results$am_pm)

model_results$am_pm <- gsub("a.<a0>m.","AM",model_results$am_pm)

##
#check which samples are pm
table(model_results$am_pm,model_results$repetition)
#View(model_results[grepl("p|P",model_results$am_pm),c("wood_sample","repetition","hour")])
#View(model_results[grepl("p|P",model_results$am_pm)&model_results$repetition==0,])
# View(model_results[grepl("a|A",model_results$am_pm)&model_results$repetition==1,])
#View(model_results[grepl("p|P",model_results$am_pm)&model_results$repetition==1,])
# View(model_results[grepl("a|A",model_results$am_pm)&model_results$repetition==2,])
#View(model_results[grepl("p|P",model_results$am_pm)&model_results$repetition==2,])
#View(model_results[grepl("p|P",model_results$am_pm)&model_results$repetition==3,])

#-- -- -- -- -- -- -- -- -- -- -- --
#### Add metadata to the results ####
#-- -- -- -- -- -- -- -- -- -- -- --

#change working directory to import metadata 
 setwd('C:\\Users\\angel\\OneDrive\\Independent_project_wood_20230823\\2_respiration_activity_20220816\\vaisala_files_txt_new') 
 metadata <- read.csv("db_respiration_complete_20220909.csv", sep= ";")
 metadata_df <- as.data.frame(metadata)
 metadata$wood_ID <- as.character(metadata$wood_ID)
 
metadata_model <- left_join(metadata,model_results, by = c("wood_ID","wood_sample"))
metadata_model_df <- as.data.frame(metadata_model)
str(metadata_model_df)

#create combined column for time since collection
table(metadata_model_df$repetition);summary(is.na(metadata_model_df$repetition))

str(metadata_model_df$wood_ID)
metadata_model_df$wood_ID <- as.numeric(metadata_model_df$wood_ID)

metadata_model_df$hours_since_collection <- ifelse(metadata_model_df$repetition==0|metadata_model_df$wood_ID>=36,metadata_model_df$time_resp0,
                                                ifelse(metadata_model_df$repetition==1,metadata_model_df$hours_resp1,
                                                       ifelse(metadata_model_df$repetition==2,metadata_model_df$hours_resp2,
                                                              metadata_model_df$hours_resp3)))
summary(as.numeric(metadata_model_df$hours_since_collection))
#the 8 NA values (as of Feb. 17 2022) were not yet weighed, producing the NA values

#drop NA values
metadata_model_df<-metadata_model_df[!is.na(metadata_model_df$hours_since_collection),]

#remove the dates, times, and "time_respX" columns from the dataframe
colnames(metadata_model_df)
colnames(metadata_model_df)[c(5,53)]<-c("field_sample_time","resp_measure_time")
metadata_model_df_trim <- select(metadata_model_df, -c(am_pm,
                                                       beaker_volume,
                                                       added_volume,
                                                       respiration_file_date,
                                                       respiration_file,
                                                       respiration_file_1_date,
                                                       respiration_file_1,
                                                       respiration_file_2_date,
                                                       respiration_file_2,
                                                       respiration_file_3_date,
                                                       respiration_file_3,
                                                       tin_number_oven))

colnames(metadata_model_df_trim)

#change working directory to save data
 setwd('C:\\Users\\angel\\Desktop\\Docs_paper\\2_Wood_respiration_macroorgan_measurements')

  #write .csv file
write.csv(metadata_model_df_trim,"respiration_exp.csv")

#-- -- -- -- -- -- -- -- -- -- -- --
#### Look at the results ####
#-- -- -- -- -- -- -- -- -- -- -- --

model_results$slope <- as.numeric(as.character(model_results$slope))

#overall ground versus elevated
ggplot(model_results,aes(x = elevated_status,y = slope)) +
  geom_boxplot(notch = TRUE) + theme_base +
  scale_x_discrete(name=element_blank(), labels =c("Elevated","Ground")) +
  scale_y_continuous(name=expression(bold('Respiration rate'~(mg[CO[2]]~'s'^-1~'g'[wood]^-1))))+ 
  scale_fill_grey()

#ANOVA for respiration rate
anova_slope <- with(model_results, aov(slope ~ elevated_status + downup + elevated_status*downup))
summary(anova_slope)
plot(TukeyHSD(anova_slope))

#before sample 29
#                               Df   Sum Sq   Mean Sq F value Pr(>F)  
#elevated_status          1 1.60e-08 1.577e-08   0.265 0.6070  
#downup                   1 2.37e-07 2.372e-07   3.992 0.0472 * (SIGNIFICANT!)
#  elevated_status:downup   1 5.30e-08 5.305e-08   0.893 0.3460  
#Residuals              180 1.07e-05 5.942e-08                 
#---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#With samples 29-35
#
#Df    Sum Sq   Mean Sq F value Pr(>F)  
#elevated_status          1 1.260e-07 1.258e-07   1.009 0.3162  
#downup                   1 3.860e-07 3.862e-07   3.098 0.0798 . (SIGNIFICANT)
#elevated_status:downup   1 1.320e-07 1.325e-07   1.063 0.3037  
#Residuals              218 2.717e-05 1.246e-07                 
#---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


#ground versus elevated, separated by up or down
resp_rate <- ggplot(model_results,aes(x = elevated_status,y = slope)) +
  geom_boxplot(aes(fill = downup), notch = TRUE)+ theme_base +
  scale_x_discrete(name=element_blank(), labels =c("Elevated","Ground")) +
  scale_y_continuous(name=expression(bold('Respiration rate'~(mg[CO[2]]~'s'^-1~'g'[wood]^-1))))+ 
  scale_fill_grey()

resp_rate

#save figure
ggsave("resp_rate_plot.tiff",resp_rate,dpi = 600, compression = "lzw",width = 3.25,height = 2,scale = 1.8)

#we need to evaluate the 
table(model_results$sample)

inter1 <- gsub("ed","",model_results$sample)
inter2 <- gsub("eu","",inter1)
inter3 <- gsub("gd","",inter2)
model_results$tree_trunk <- gsub("gu","",inter3)

