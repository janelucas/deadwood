#Feb_02_2022 
#Preliminary visualization of wood survey on 50-ha plot on BCI
#csv file used -> woodecomp_data_bci_02022022.csv
#Edited by Angela Marcela Barrera Bello and Evan Gora

#Tasks:
#1. Number of trees recorded (ok)
#2. Number of species and families (ok)

#a. How elevation changes over time (geom_point) (ok)
#b. How decomp (mass remaining) changes over time (year_dead) (ok)
#c. Species in each level of decomp. (ok)

#9. Histogram: doubtful vs. certain data (ok)

#10.Map of all the dead trees (show dot size according to dbh and dot intensity according to mass remaining) 


#Packages
#install.packages("tidyverse")
#install.packages("tidyr") 
#install.packages("gridExtra")
#For more palettes of colors (bars ggplot)
#install.packages("devtools")
#devtools::install_github("jaredhuling/jcolors", force = TRUE)
#install.packages("ggplot2")


library(extrafont)
library(dplyr)
library(tidyverse)
library(tidyr)
library(readr)
library(ggplot2)
library(mctoolsr)
library(readr)
library(jcolors)
library(gridExtra)
library(ggplot2)

rm(list=ls())

#remotes::install_version("Rttf2pt1", version = "1.3.8")
# extrafont::font_import()
# font_import()
loadfonts(device="win")

#create theme basis for ggplot
theme_basis<-theme(axis.title = element_text(family = "Arial", color="black", face="bold", size=14))+
  theme(axis.text=element_text(family = "Arial",colour="black", face ="bold", size = 12)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"),axis.ticks = element_line(colour="black"))

load("C:\\Users\\angel\\OneDrive\\Independent_project_wood_20230823\\1_wood_survey_20220906\\analysis\\BCI_dead_trees.RData")

#drop repeated trees
table(table(BCI_dead_trees$tag))

census.identifier <- BCI_dead_trees[,c("tag","source")]

bci.spptable <- read_csv("C:\\Users\\angel\\OneDrive\\Independent_project_wood_20230823\\1_wood_survey_20220906\\data\\bci.spptable.csv")
str(bci.spptable)

#set the working directory
setwd("C:\\Users\\angel\\OneDrive\\Independent_project_wood_20230823\\1_wood_survey_20220906\\data")

#.csv files needed: 
#Latest version of wood survey (woodecomp_data_BCI.csv)
#Table of species of the 50-ha plot (bci.spptable.csv)

#import .csv files 
ws_50ha_bci <- read_csv("woodecomp_data_BCI_20220508_cleaned.csv")
ws_50ha_bci$sp <- lapply(ws_50ha_bci$sp, tolower) # we changed uppercase to lowercase to unified names
ws_50ha_bci$sp <- as.character(ws_50ha_bci$sp)
View(ws_50ha_bci)
ws_50ha_bci <- as.data.frame(ws_50ha_bci)

bci_spptable <- as.data.frame(bci.spptable)
str(bci_spptable)
str(ws_50ha_bci)

#Combining these two data frames by "sp"
ws_50ha <- left_join(ws_50ha_bci, bci_spptable[,c(1,2,5)], by = "sp")

#----------------------~
##1. Trees recorded ####
#----------------------~

#lets drop NA of spp_ws_50ha$Number to get just trees recorded in the DB and
# remove dbh less tha 10 cm (100 mm)

ws_50ha_nona <- ws_50ha %>% 
  drop_na(Number) %>% 
  filter(!dbh.mm < 100)
  
length(ws_50ha_nona$Number) #15652 trees

#modify names in ws_50ha_nona

colnames(ws_50ha_nona)[9]   <- "dbh"
colnames(ws_50ha_nona)[10]   <- "year_dead"
colnames(ws_50ha_nona)[13]   <- "elev"

#[1] 16110 trees were recorded as dead between 1985 and 2019 in the 50-ha plot

#>2018 to 2019
ws_50ha_nona$year_dead <- as.numeric(as.character(gsub(">2018", 2019, ws_50ha_nona$year_dead)))
str(ws_50ha_nona)

#-----------------------------~
##2. No. spp and families ####
#-----------------------------~

names(ws_50ha_nona)
#lets make a table with data from ws_50ha_nona arranged by spp, families, abundance, DBH, Yr.dead
spp_ws_50ha <- ws_50ha_nona %>%
  group_by(sp,Latin, Family) %>%
  summarise(abundance = length(Latin),
            minDBH = min(dbh),
            maxDBH = max(dbh),
            minyear = min(year_dead),
            maxyear = max(year_dead)) %>%
  arrange(desc(abundance))

spp_ws_50ha_df <- as.data.frame(spp_ws_50ha)

#lets drop NA of spp_ws_50ha$Latin to get just identified species
spp_ws_50ha_df_nona <- spp_ws_50ha %>% drop_na(Latin)
length(spp_ws_50ha_df_nona$Latin)
#[1] 236 species of trees recorded as dead between 1985 and 2019 in the 50-ha plot

#look at the most common species still remaining in the plot
spp_ws_50ha_df_nona_persist <- ws_50ha_nona %>% 
  filter(as.numeric(decomp) > 0&!is.na(as.numeric(decomp)))

spp_ws_50ha_persist <- spp_ws_50ha_df_nona_persist %>%
  group_by(sp,Latin, Family) %>%
  summarise(abundance = length(Latin),
            minDBH = min(dbh),
            maxDBH = max(dbh),
            minyear = min(year_dead),
            maxyear = max(year_dead)) %>%
  arrange(desc(abundance))
head(spp_ws_50ha_persist,10)

fam_ws_50ha_df_nona <- as.data.frame(table(spp_ws_50ha_df_nona$Family))
colnames(fam_ws_50ha_df_nona) <- c("family", "frequence")
length(fam_ws_50ha_df_nona$family)
#[1] 57 families of trees recorded as dead between 1985 and 2019 in the 50-ha plot

#--------------------------------------~
##3. How elevation changes over time####
#--------------------------------------~

#drop NA and NA_ in ws_50ha_nona$elev
ws_50ha_nona_elev <- ws_50ha_nona[!grepl('_|ex', ws_50ha_nona$elev),]
ws_50ha_nona_elev1 <- ws_50ha_nona_elev  %>% drop_na(elev)
str(ws_50ha_nona_elev1)
#ws_50ha_nona_elev1$elev as numeric
ws_50ha_nona_elev1$elev <- as.numeric(ws_50ha_nona_elev1$elev)
str(ws_50ha_nona_elev1)
#>2018 to 2019
ws_50ha_nona_elev1$year_dead <- as.numeric(as.character(gsub(">2018", 2019, ws_50ha_nona_elev1$year_dead)))
str(ws_50ha_nona_elev1)

#plot
ggplot(ws_50ha_nona_elev1, aes(y = elev, x = year_dead)) + 
  geom_point(size = .4, alpha=.3) + geom_smooth(method="loess") + theme_basis +
  xlab("Year") + 
  ylab("% Elevation")

#remove those with uncertain decomp values
ws_50ha_nona_elev1$decomp_num <- as.numeric(ws_50ha_nona_elev1$decomp)
str(ws_50ha_nona_elev1)
ws_50ha_nona_elev1$mass_remaining <- ifelse(ws_50ha_nona_elev1$decomp_num==0,2.5,
                                             ifelse(ws_50ha_nona_elev1$decomp_num==1,12.5,
                                                    ws_50ha_nona_elev1$decomp_num*20-10))
length(unique(ws_50ha_nona_elev1$tag)) #826 trees with certain decomp values

#look at All trees remaning (decomp 1-5):
ws_50ha_nona_elev1_decomp_1to5 <- filter(ws_50ha_nona_elev1, decomp >= 1)
length(unique(ws_50ha_nona_elev1_decomp_1to5$tag))

#keep certain column_to_rownames
ws_50ha_nona_elev1_decomp_1to5_a <- ws_50ha_nona_elev1_decomp_1to5[ , names(ws_50ha_nona_elev1_decomp_1to5) %in% c("elev", "mass_remaining")] 
names(ws_50ha_nona_elev1_decomp_1to5_a)[names(ws_50ha_nona_elev1_decomp_1to5_a) == 'mass_remaining'] <- 'volume_remaining'
names(ws_50ha_nona_elev1_decomp_1to5_a)[names(ws_50ha_nona_elev1_decomp_1to5_a) == 'elev'] <- 'percentage_suspended'

#save as csv
write.csv(ws_50ha_nona_elev1_decomp_1to5_a, "C:\\Users\\angel\\Desktop\\Docs_paper\\1_wood_decomp_survey\\ws_50ha_susp_decomp_1to5.csv", row.names=FALSE)

#create variable defining 5-year census intervals versus 1 year intervals
census.identifier$tag <- as.numeric(census.identifier$tag)
ws_50ha_nona_elev1 <- left_join(ws_50ha_nona_elev1,census.identifier, by = "tag")
ws_50ha_nona_elev1$years_since_death <- 2021 - ws_50ha_nona_elev1$year_dead
summary(is.na(ws_50ha_nona_elev1$source));table(ws_50ha_nona_elev1$source)

#Now we want to plot the standard error of the mean of each year
data_se_elev <- ws_50ha_nona_elev1 %>%
  group_by(source,year_dead,years_since_death) %>%
  summarise(mean_elev = mean(elev),
            sd_elev = sd(elev),
            se_elev = sd(elev)/sqrt(sum(!is.na(elev))),  
            upper_elev = mean_elev + se_elev,
            lower_elev = mean_elev - se_elev,
            mean_decomp = mean(decomp_num,na.rm=T),
            se_decomp = sd(decomp_num)/sqrt(sum(!is.na(decomp_num))))


ggplot(data_se_elev, aes(y = mean_elev, x = years_since_death)) + 
  geom_point(size = 2, alpha=.3, aes(shape = source)) +
  scale_shape_manual(values = c(15,16))+
  # geom_line()+
  geom_errorbar(aes(ymin = lower_elev, ymax = upper_elev), width = 0.2) + theme_basis+
  xlab("Year") + 
  ylab("% Elevation")

##Plot Elev Vs. Remaining mass####
#here we arranged info of elevation by decomp

data_se_elev_decomp <- ws_50ha_nona_elev1 %>%
  group_by(decomp_num) %>%
  summarise(mean_elev = mean(elev),
            sd_elev = sd(elev),
            se_elev = sd(elev)/sqrt(sum(!is.na(elev))),  
            upper_elev = mean_elev + se_elev,
            lower_elev = mean_elev - se_elev)

#create inverse decomp number to plot from highest to lowest
data_se_elev_decomp$inv_decomp_num <- 5-data_se_elev_decomp$decomp_num

ggplot(data_se_elev_decomp, aes(y = mean_elev, x = inv_decomp_num)) + 
  geom_point(size = 2) + 
  geom_errorbar(aes(ymin = lower_elev, ymax = upper_elev), width = 0.2) + theme_basis+
  scale_x_continuous("Remaining mass (%)",breaks=seq(0,5,1),
                     labels =c("100-80","80-60","60-40","40-20","20-5","<5")) +
  scale_y_continuous("Percent suspended",breaks = seq(20,80,10),limits = c(30,80),
                     labels = c("20","","40","","60","","80"),expand=c(0,0)) 

elevation_by_massloss <- ggplot(data_se_elev_decomp, aes(y = mean_elev, x = inv_decomp_num)) + 
  geom_point(size = 2) + 
  geom_errorbar(aes(ymin = lower_elev, ymax = upper_elev), width = 0.2) + theme_basis+
  scale_x_continuous("Volume loss (%)",breaks=seq(0,5,1),
                   labels =c("0-20","20-40","40-60","60-80","80-95",">95")) +
  scale_y_continuous("Percent suspended",breaks = seq(20,80,10),limits = c(30,80),
                     labels = c("20","","40","","60","","80"),expand=c(0,0)) 

elevation_by_massloss

setwd("C:/Users/angel/OneDrive/1_wood_survey_20220906/figures")
ggsave("elevation_by_massloss_bigger_equal_10dbh_20220906.tiff",elevation_by_massloss,dpi=600,compression="lzw",
       width= 3.25, height = 1.8,scale= 1.8)
getwd()

###   Plot for the 3 most common species####
#1) Here we are seeing individuals with >0 remaining mass... so we are missing all
#those trees = 0 remaining mass (biggest part, O NA NA). If we want to add them we 
#need to create a condiction to filer just O NA NA from NA NA NA 
#2) Pablo recommeded me to remove the tendency line, better right?

#faraoc
ws_50ha_nona_elev1_faraoc <- filter(ws_50ha_nona_elev1, sp == "faraoc")

data_se_elev_faraoc <- ws_50ha_nona_elev1_faraoc %>%
  group_by(source,year_dead,years_since_death) %>%
  summarise(mean_elev = mean(elev),
            sd_elev = sd(elev),
            se_elev = sd(elev)/sqrt(sum(!is.na(elev))),  
            upper_elev = mean_elev + se_elev,
            lower_elev = mean_elev - se_elev)

ggplot(data_se_elev_faraoc, aes(y = mean_elev, x = years_since_death)) + 
  geom_point(size = .4, alpha=.3) + 
  # geom_line()+
  geom_errorbar(aes(ymin = lower_elev, ymax = upper_elev), width = 0.2) + theme_basis+
  xlab("Year") + 
  ylab("% Elevation")+ggtitle("faraoc")

#tri2tu

ws_50ha_nona_elev1_tri2tu <- filter(ws_50ha_nona_elev1, sp == "tri2tu")

sort(table(ws_50ha_nona_elev1_tri2tu$tag))
#remove bad trees
keep1 <- ws_50ha_nona_elev1_tri2tu$tag=="831"&ws_50ha_nona_elev1_tri2tu$source=="dendro"&ws_50ha_nona_elev1_tri2tu$years_since_death==7
summary(keep1)
keep2 <- ws_50ha_nona_elev1_tri2tu$tag=="1876"&ws_50ha_nona_elev1_tri2tu$source=="dendro"&ws_50ha_nona_elev1_tri2tu$years_since_death==9
summary(keep2)
keep3 <- ws_50ha_nona_elev1_tri2tu$tag=="1964"&ws_50ha_nona_elev1_tri2tu$years_since_death==8&ws_50ha_nona_elev1_tri2tu$source=="dendro"
summary(keep3)
keep4 <- ws_50ha_nona_elev1_tri2tu$tag=="4165"&ws_50ha_nona_elev1_tri2tu$years_since_death==6&ws_50ha_nona_elev1_tri2tu$source=="dendro"
summary(keep4)
keep5 <- ws_50ha_nona_elev1_tri2tu$tag=="4167"&ws_50ha_nona_elev1_tri2tu$years_since_death==6&ws_50ha_nona_elev1_tri2tu$source=="dendro"
summary(keep5)
keep6 <- ws_50ha_nona_elev1_tri2tu$tag=="5152"&ws_50ha_nona_elev1_tri2tu$years_since_death==13&ws_50ha_nona_elev1_tri2tu$source=="dendro"
summary(keep6)
keep7 <- ws_50ha_nona_elev1_tri2tu$tag=="5451"&ws_50ha_nona_elev1_tri2tu$years_since_death==12&ws_50ha_nona_elev1_tri2tu$source=="dendro"
summary(keep7)
keep8 <- ws_50ha_nona_elev1_tri2tu$tag=="7295"&ws_50ha_nona_elev1_tri2tu$years_since_death==8&ws_50ha_nona_elev1_tri2tu$source=="dendro"
summary(keep8)

#select good trees from the erronious groups
rows_to_keep <- ws_50ha_nona_elev1_tri2tu[keep1|keep2|keep3|keep4|keep5|keep6|keep7|keep8,]

#drop all bad tags
bad_tags <- names(table(ws_50ha_nona_elev1_tri2tu$tag)[table(ws_50ha_nona_elev1_tri2tu$tag)==4])
ws_50ha_nona_elev1_tri2tu_dropbad <- ws_50ha_nona_elev1_tri2tu[!ws_50ha_nona_elev1_tri2tu$tag %in% bad_tags,]

#add back the good values
ws_50ha_nona_elev1_tri2tu_trim <- bind_rows(ws_50ha_nona_elev1_tri2tu_dropbad,rows_to_keep)
length(unique(ws_50ha_nona_elev1_tri2tu_trim$tag)) #104

# tri tu still remaining (decomp 1-5)
ws_50ha_nona_tri2tu_decomp_1to5 <- filter(ws_50ha_nona_elev1_tri2tu_trim, decomp >= 1)
length(unique(ws_50ha_nona_tri2tu_decomp_1to5$tag)) # 75

# keep certain column_to_rownames
colnames(ws_50ha_nona_tri2tu_decomp_1to5)

ws_50ha_nona_tri2tu_decomp_1to5_a <- ws_50ha_nona_tri2tu_decomp_1to5[ , names(ws_50ha_nona_tri2tu_decomp_1to5) %in% c("source.x", "years_since_death", "elev","mass_remaining")] 
names(ws_50ha_nona_tri2tu_decomp_1to5_a)[names(ws_50ha_nona_tri2tu_decomp_1to5_a) == 'source.x'] <- 'source'
names(ws_50ha_nona_tri2tu_decomp_1to5_a)[names(ws_50ha_nona_tri2tu_decomp_1to5_a) == 'elev'] <- 'percentage_suspended'

# mass_remaining to Volume loss (%)
ws_50ha_nona_tri2tu_decomp_1to5_a$Volume_loss <- 100-ws_50ha_nona_tri2tu_decomp_1to5_a$mass_remaining

names(ws_50ha_nona_tri2tu_decomp_1to5_a)[names(ws_50ha_nona_tri2tu_decomp_1to5_a) == 'Volume_loss'] <- 'Volume loss (%)'
ws_50ha_nona_tri2tu_decomp_1to5_b <- ws_50ha_nona_tri2tu_decomp_1to5_a[ , names(ws_50ha_nona_tri2tu_decomp_1to5_a) %in% c("source", "years_since_death", "percentage_suspended","Volume loss (%)")] 


#save as csv
write.csv(ws_50ha_nona_tri2tu_decomp_1to5_b, "C:\\Users\\angel\\Desktop\\Docs_paper\\1_wood_decomp_survey\\tri2tu_all_remaining.csv", row.names=FALSE)

# tri tu still remaining (decomp 1-5) and downed
ws_50ha_nona_tri2tu_decomp_1to5_d <- filter(ws_50ha_nona_tri2tu_decomp_1to5, Position == "d")
length(unique(ws_50ha_nona_tri2tu_decomp_1to5_d$tag)) # 33

data_se_elev_tri2tu <- ws_50ha_nona_elev1_tri2tu_trim %>%
  group_by(source,years_since_death) %>%
  summarise(mean_elev = mean(elev),
            sd_elev = sd(elev),
            se_elev = sd(elev)/sqrt(sum(!is.na(elev))),  
            upper_elev = mean_elev + se_elev,
            lower_elev = mean_elev - se_elev,
            sample_size = sum(!is.na(elev)),
            mean_decomp = mean(mass_remaining,na.rm=T),
            se_decomp = sd(mass_remaining)/sqrt(sum(!is.na(mass_remaining))))
length(data_se_elev_tri2tu)
#11

ggplot(data_se_elev_tri2tu, aes(y = mean_elev, x = years_since_death)) + 
  geom_point(size = 2, aes(shape = source)) +
  scale_shape_manual(values = c(15,16), labels = c("Annual Surveys","Full surveys"))+
  geom_errorbar(aes(ymin = lower_elev, ymax = upper_elev), width = 0.2) + theme_basis+
  xlab("Year") + 
  facet_wrap(~source,scales="free",ncol=1) + theme(legend.position = c(.9,.4), 
  # theme(legend.position = c(.8,.4), 
        legend.title = element_blank(),
        legend.background = element_blank(),legend.key = element_blank()) +
  # geom_smooth(method = "loess", se=F) + 
  ylab("% Elevation")+ggtitle("tri2tu")


ggplot(data_se_elev_tri2tu, aes(y = mean_elev, x = years_since_death)) + 
  geom_point(aes(size = mean_decomp)) + 
  # geom_line()+
  geom_errorbar(aes(ymin = lower_elev, ymax = upper_elev), width = 0.2) + theme_basis+
  xlab("Year of death") +
  theme(legend.position = c(.8,.4), #legend.title = element_blank(),
        legend.background = element_blank(),legend.key = element_blank()) +
  ylab("% Elevation")+ggtitle("tri2tu")

#separate by color and size
data_se_elev_tri2tu$Massloss <- 100-data_se_elev_tri2tu$mean_decomp
colnames(data_se_elev_tri2tu)[11] <- "Volume loss (%)"

data_sources <- c(dendro = "Annual surveys",
                  plot.census = "5-year surveys")

install.packages("remotes")
remotes::install_github("csdaw/ggbeeswarm2")
library(ggbeeswarm2)

elevation_over_time <- ggplot(data_se_elev_tri2tu, aes(y = mean_elev, x = years_since_death)) + 
  geom_point(aes(size = sample_size, color = `Volume loss (%)`)) +
  facet_wrap(~source, ncol=1,labeller = as_labeller(data_sources),scales="free_x") +
  scale_color_gradient(low = "black", high = "lightgrey",limits = c(0,100)) +
  guides(size=guide_legend(title="Sample size"),
         shape=guide_legend(title="Data source")) +
  scale_shape_manual(values = c(15,16), labels = c("Annual Surveys","Full surveys"))+
  geom_errorbar(aes(ymin = lower_elev, ymax = upper_elev), width = 0.2) +
  scale_x_continuous("Years since death",limits = c(0,35), breaks = seq(0,35,5),
                     expand=c(0,0),labels=c("0","","10","","20","","30","")) +
  scale_y_continuous("Percent suspended") +
  theme(legend.position=c(.8,.8),legend.background = element_blank(),
        legend.key = element_blank(),legend.box = "horizontal",
        legend.text = element_text(family="Arial",size=10,face="bold"),
        legend.title = element_text(family="Arial",size=11,face="bold"),
        legend.key.height = unit(.4,"cm"),legend.spacing.x = unit(.01,"cm")) + 
  theme_basis +theme(strip.background = element_blank(),
                     strip.text = element_text(family="Arial",size=13,face="bold"))

elevation_over_time

ggsave("elevation_over_time_bigger_equal_10dbh_20230401.tiff",elevation_over_time,dpi=600,compression="lzw",
       width= 3.25, height = 3,scale=1.8)
getwd()

#sample size
data_se_elev_tri2tu %>% 
  group_by(source) %>%
  summarize(total_trees = sum(sample_size))

data_se_elev_tri2tu$source


#keep certain column_to_rownames
colnames(data_se_elev_tri2tu)

data_se_elev_tri2tu_a <- data_se_elev_tri2tu[ , names(data_se_elev_tri2tu) %in% c("source", "years_since_death", "mean_elev","Volume loss (%)","sample_size")] 
names(data_se_elev_tri2tu_a)[names(data_se_elev_tri2tu_a) == 'mean_elev'] <- 'mean_percentage_suspended'


#save as csv
write.csv(data_se_elev_tri2tu_a, "C:\\Users\\angel\\Desktop\\Docs_paper\\1_wood_decomp_survey\\tri2tu_by_years_since_death.csv", row.names=FALSE)


#look at the same patterns, but only for downed wood
data_se_elev_tri2tu_downed <- ws_50ha_nona_elev1_tri2tu_trim %>%
  filter(grepl("d",ws_50ha_nona_elev1_tri2tu_trim$Position)) %>%
  group_by(source,year_dead,years_since_death) %>%
  summarise(mean_elev = mean(elev),
            sd_elev = sd(elev),
            se_elev = sd(elev)/sqrt(sum(!is.na(elev))),  
            upper_elev = mean_elev + se_elev,
            lower_elev = mean_elev - se_elev,
            sample_size = sum(!is.na(elev)),
            mean_decomp = mean(mass_remaining,na.rm=T),
            se_decomp = sd(mass_remaining)/sqrt(sum(!is.na(mass_remaining))))

length(data_se_elev_tri2tu_downed)
# 11

#sample size
data_se_elev_tri2tu_downed %>% 
  group_by(source) %>%
  summarize(total_trees = sum(sample_size))

#separate by color and size
data_se_elev_tri2tu_downed$Massloss <- 100-data_se_elev_tri2tu_downed$mean_decomp
colnames(data_se_elev_tri2tu_downed)[12] <- "Volume loss (%)"

data_sources <- c(dendro = "Annual surveys",
                  plot.census = "5-year surveys")

elevation_over_time_downed <- ggplot(data_se_elev_tri2tu_downed, aes(y = mean_elev, x = years_since_death)) + 
  geom_point(aes(size = sample_size, color = `Volume loss (%)`)) +
  facet_wrap(~source, ncol=1,labeller = as_labeller(data_sources),scales="free_x") +
  scale_color_gradient(low = "black", high = "lightgrey",limits = c(0,100)) +
  guides(size=guide_legend(title="Sample size"),
         shape=guide_legend(title="Data source")) +
  scale_shape_manual(values = c(15,16), labels = c("Annual Surveys","Full surveys"))+
  geom_errorbar(aes(ymin = lower_elev, ymax = upper_elev), width = 0.2) +
  scale_x_continuous("Years since death",limits = c(0,35), breaks = seq(0,35,5),
                     expand=c(0,0),labels=c("0","","10","","20","","30","")) +
  scale_y_continuous("Percent suspended") +
  theme(legend.position=c(.8,.8),legend.background = element_blank(),
        legend.key = element_blank(),legend.box = "horizontal",
        legend.text = element_text(family="Arial",size=10,face="bold"),
        legend.title = element_text(family="Arial",size=11,face="bold"),
        legend.key.height = unit(.4,"cm"),legend.spacing.x = unit(.01,"cm")) + 
  theme_basis +theme(strip.background = element_blank(),
                     strip.text = element_text(family="Arial",size=13,face="bold"))

elevation_over_time_downed

ggsave("elevation_over_time_downed_bigger_equal_10dbh_20230401.tiff",elevation_over_time_downed,dpi=600,compression="lzw",
       width= 3.25, height = 3,scale=1.8)
getwd()

####t-test for elevation> dead trees over 10 years (everything)#### USE THIS ONE!

sort(unique(ws_50ha_nona_elev1_tri2tu_trim$years_since_death))

#dead trees over 5 years
ws_50ha_nona_elev1_tri2tu_trim_over5y <- filter(ws_50ha_nona_elev1_tri2tu_trim, years_since_death > 5)
length(unique(ws_50ha_nona_elev1_tri2tu_trim_over5y$tag)) #72

#dead trees over 10 years
ws_50ha_nona_elev1_tri2tu_trim_over10y <- filter(ws_50ha_nona_elev1_tri2tu_trim, years_since_death > 10)
length(unique(ws_50ha_nona_elev1_tri2tu_trim_over10y$tag)) #25

#dead trees 5-10 years dead
ws_50ha_nona_elev1_tri2tu_trim_5to10 <- filter(ws_50ha_nona_elev1_tri2tu_trim, years_since_death == 5|
                                                 years_since_death == 6|
                                                 years_since_death == 7|
                                                 years_since_death == 8|
                                                 years_since_death == 8.5|
                                                 years_since_death == 9)

t.test(ws_50ha_nona_elev1_tri2tu_trim_over10y$elev, ws_50ha_nona_elev1_tri2tu_trim_5to10$elev, paired=FALSE)
# data:  ws_50ha_nona_elev1_tri2tu_trim_over10y$elev and ws_50ha_nona_elev1_tri2tu_trim_5to10$elev
# t = 5.0285, df = 43.701, p-value = 8.905e-06 
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   28.03483 65.54979
# sample estimates:
#   mean of x mean of y 
# 71.60000  24.80769 

####t-test for elevation> dead trees over 10 years (only d)#### 

#dead trees over 10 years
ws_50ha_nona_elev1_tri2tu_trim_over10y_d <- filter(ws_50ha_nona_elev1_tri2tu_trim_over10y, Position == 'd')
length(unique(ws_50ha_nona_elev1_tri2tu_trim_over10y_d$tag)) #5

#dead trees 5-10 years dead
ws_50ha_nona_elev1_tri2tu_trim_5to10_d <- filter(ws_50ha_nona_elev1_tri2tu_trim_5to10, Position == 'd')
length(ws_50ha_nona_elev1_tri2tu_trim_5to10_d) #26

#t-test
t.test(ws_50ha_nona_elev1_tri2tu_trim_over10y_d$elev, ws_50ha_nona_elev1_tri2tu_trim_5to10_d$elev, paired=FALSE)

# data:  ws_50ha_nona_elev1_tri2tu_trim_over10y_d$elev and ws_50ha_nona_elev1_tri2tu_trim_5to10_d$elev
# t = -1.6091, df = 32.603, p-value = 0.1172
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -14.974184   1.751962
# sample estimates:
#   mean of x mean of y 
# 2.000000  8.611111 

####t-test for elevation> dead trees over 10 years (something with d (d, d_st, d_st_su))#### USE THIS ONE!

#dead trees over 10 years
ws_50ha_nona_elev1_tri2tu_trim_over10y_withd <- ws_50ha_nona_elev1_tri2tu_trim_over10y[grepl("d", ws_50ha_nona_elev1_tri2tu_trim_over10y$Position),]
length(ws_50ha_nona_elev1_tri2tu_trim_over10y_withd)

#dead trees 5-10 years dead
ws_50ha_nona_elev1_tri2tu_trim_5to10_withd <- ws_50ha_nona_elev1_tri2tu_trim_5to10[grepl("d", ws_50ha_nona_elev1_tri2tu_trim_5to10$Position),]

#t-test
t.test(ws_50ha_nona_elev1_tri2tu_trim_over10y_withd$elev, ws_50ha_nona_elev1_tri2tu_trim_5to10_withd$elev, paired=FALSE)

# data:  ws_50ha_nona_elev1_tri2tu_trim_over10y_withd$elev and ws_50ha_nona_elev1_tri2tu_trim_5to10_withd$elev
# t = 1.4869, df = 11.534, p-value = 0.1639
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -8.351377 43.742682
# sample estimates:
#   mean of x mean of y 
# 34.00000  16.30435 


#poular

ws_50ha_nona_elev1_poular <- filter(ws_50ha_nona_elev1, sp == "poular")

data_se_elev_poular <- ws_50ha_nona_elev1_poular %>%
  group_by(year_dead) %>%
  summarise(mean_elev = mean(elev),
            sd_elev = sd(elev),
            se_elev = sd(elev)/sqrt(sum(!is.na(elev))),  
            upper_elev = mean_elev + se_elev,
            lower_elev = mean_elev - se_elev)

ggplot(data_se_elev_poular, aes(y = mean_elev, x = year_dead)) + 
  geom_point(size = .4, alpha=.3) + 
  geom_line()+
  geom_errorbar(aes(ymin = lower_elev, ymax = upper_elev), width = 0.2) + theme_basis+
  xlab("Year") + 
  ylab("% Elevation")+ggtitle("poular")


  
#--------------------------------------~
##4. How decomp (mass remaining) changes over time (year_dead) ####
#--------------------------------------~

#drop NA and NA_ in ws_50ha_nona$decomp
ws_50ha_nona_decomp <- ws_50ha_nona[!grepl('_|ex|d|c', ws_50ha_nona$decomp),]
ws_50ha_nona_decomp1 <- ws_50ha_nona_decomp  %>% drop_na(decomp)
str(ws_50ha_nona_decomp1)
#ws_50ha_nona_elev1$elev as numeric
ws_50ha_nona_decomp1$decomp <- as.numeric(ws_50ha_nona_decomp1$decomp)
str(ws_50ha_nona_decomp1$decomp)

str(ws_50ha_nona_decomp1$decomp) #This one doesnt have NA, NA_0, and others

#>2018 to 2019
ws_50ha_nona_decomp1$year_dead <- as.numeric(as.character(gsub(">2018", 2019, ws_50ha_nona_decomp1$year_dead)))
str(ws_50ha_nona_decomp1)


#plot
ggplot(ws_50ha_nona_decomp1, aes(y = decomp, x = year_dead)) + 
  geom_point(size = .4, alpha=.3) + geom_smooth(method="loess") + theme_basis + 
  xlab("Year dead") + 
  ylab("Mass remaining")
  
#Now a plot with the standard error 

data_se_mass <- ws_50ha_nona_decomp1 %>%
  group_by(year_dead) %>%
  summarise(mean_decomp = mean(decomp),
            sd_decomp = sd(decomp),
            se_decomp = sd(decomp)/sqrt(sum(!is.na(decomp))),  
            upper_decomp = mean_decomp + se_decomp,
            lower_decomp = mean_decomp - se_decomp)

ggplot(data_se_mass, aes(y = mean_decomp, x = year_dead)) + 
  geom_point(size = .4, alpha=.3) + 
  geom_line()+
  geom_errorbar(aes(ymin = lower_decomp, ymax = upper_decomp), width = 0.2) + theme_basis+
  xlab("Year") + 
  ylab("Mass remaining")

#--------------------------------------~
##5. Species in each level of decomp####
#--------------------------------------~
#which spp is most common in each level of decomp

#lets make the table
species_decomp <- table(ws_50ha_nona_decomp1$decomp,ws_50ha_nona_decomp1$sp)
species_decomp_df <- as.data.frame(species_decomp)
#drop species that have 0 frequence
colnames(species_decomp_df)[1] <- "decomp"
colnames(species_decomp_df)[2] <- "sp_code"

max_sp_mass <- species_decomp_df %>%
  group_by(decomp) %>%
  arrange(desc(Freq)) %>%
  slice_head(n=3) # change here to see more species

max_sp_mass$sp_code <- factor(max_sp_mass$sp_code, levels=unique(max_sp_mass$sp_code))
 
 
#table with names of species
species_codes_mass <- as.data.frame(unique(max_sp_mass$sp_code))
colnames(species_codes_mass)[1] <- "sp"
species_codes_mass_names <- left_join(species_codes_mass, bci_spptable[,c(1,2)], by = "sp")
grid.table(species_codes_mass_names)

 
#plot

ggplot(data = max_sp_mass, aes(x = decomp, y = Freq, fill = as.factor(sp_code))) +
  geom_bar(stat = "identity", colour = "black", position = position_dodge())+
  labs(x = "Mass remaining", y = "Frequency", fill = "Species") +
  geom_text(aes(label=sp_code), position=position_dodge(width=0.9),
            vjust=-0.0005, hjust=0.000010,angle=90) + 
  theme_basis +
  theme(legend.position = "none") + 
  scale_fill_grey()+
  ylim(0, 1200)


#----------------------~
##6. Rank abundance: spp ####
#----------------------~

spp_ws_50ha_df_nona_10 <- spp_ws_50ha_df_nona %>%
  arrange(desc(abundance)) %>%
  head(10)
spp_ws_50ha_df_nona_10$sp <- factor(spp_ws_50ha_df_nona_10$sp, levels=unique(spp_ws_50ha_df_nona_10$sp))
str(spp_ws_50ha_df_nona_10)

ggplot(spp_ws_50ha_df_nona_10,aes(x = sp, y = abundance, fill= sp))+
  geom_bar(stat = "identity")+ 
  labs(x = "Species", y = "Abundance") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme_basis +
  theme(legend.position = "none") + 
  scale_fill_grey()

#table with names of species
grid.table(species_codes_10_names) #fix >2018 (run all again?)

##7. Histogram: Data certain vs. uncertain ####

#remaining mass:
cert_data_dec <- with(ws_50ha_nona, tapply(sp, decomp, length))
cert_data_dec_df <- as.data.frame(t(cert_data_dec))
# Certain data (remaining mass) = 12173 (0= 11496)
# Uncertain data (remaining mass) = 3937

(12173*100)/16110 # working with 75% (certain data remaining mass)
(11496*100)/12173 # 94% of certain data remaining mass is 0!
(11496*100)/16110  # 71% of all remaining mass is 0!
(677*100)/16110 # 8% of all remaining mass is >0 (exists)!
(3690*100)/16110 # 22% of all remaining mass is NA
  
#position:
cert_data_pos <- with(ws_50ha_nona, tapply(sp, Position, length))
cert_data_pos_df <- as.data.frame(t(cert_data_pos))
# Certain data (position) =  888 (any position)
# Uncertain data (position) or gone (just NA) = 380 (NA_)+ 14842 (NA)=15222

(888*100)/16110 #working with 5% (certain data position)
(439*100)/888 # 49%  of certain data position are completly downed(d) 
(185*100)/888 # 20% of certain data position are  completly standing (st) "snags"


#elevation:
cert_data_elev <- with(ws_50ha_nona, tapply(sp, elev, length))
cert_data_elev_df <- as.data.frame(t(cert_data_elev))
# Certain data (position) =  885 (any elevation)
# Uncertain data (position) or gone (just NA) = 15225 

(885*100)/16110 #working with 6% (certain data position)

#plot for each one!!!!
 
## 10.Map of all the dead trees #### (same general code of the maps/ remove the lines!! 5x5 m)
#Look in ws_bci_50ha_fullmap_bef_after_02122022 script

##11. mass ~ elev ? #### 
ws_50ha_nona_elev1$elev 
str(ws_50ha_nona_decomp1)

levels

aov(decomp ~ elev, data = ws_50ha_nona)

