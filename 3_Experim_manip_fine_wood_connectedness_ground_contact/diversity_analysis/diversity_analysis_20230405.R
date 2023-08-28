#Analysis of bacterial and fungal data
#Made from previous R script "microbes_input" <- 'C:/Users/angel/OneDrive/3_microbial_diversity_20220816/Analysis'
#Edited by EM Gora 2021-6-30

#install packages
#install.packages("devtools")
#devtools::install_github("leffj/mctoolsr", force = TRUE)
#install.packages("vegan")
#install.packages("lme4")
#install.packages("mctoolsr")

#load packages
library(mctoolsr)
library(vegan)
library(dplyr)
library(ggplot2)
library(grid)
library(lme4)
library(car)

#remove
rm(list=ls())

setwd('C:/Users/angel/OneDrive/3_microbial_diversity_20220816/Exploration_results')

#Import lists 
input_bact_rarfy <- readRDS("input_bact_rarfy.RData")
input_fungi_rarfy <- readRDS("input_fungi_rarfy.RData")

#////////////////////////
## 1. VISUALIZE DATA ####
#///////////////////////
# 
# #BACTERIA: look at initial ordination
 dm = calc_dm(input_bact_rarfy$data_loaded)
 ord = calc_ordination(dm, 'nmds')
 plot_ordination(input_bact_rarfy,ord,'groundcontact_treatment','groundcontact_treatment',hulls=TRUE)
# 
# #this graphic shows there is a sample that is very different from the others and that doesnt 
# # ...let us analyze the rest of the samples. We must filter it out.
# 
 #FUNGI: look at initial ordination
 dm = calc_dm(input_fungi_rarfy$data_loaded)
 ord = calc_ordination(dm, 'nmds')
 plot_ordination(input_fungi_rarfy,ord,'groundcontact_treatment','groundcontact_treatment',hulls=TRUE)
# # other plot
 dm = calc_dm(input_fungi_rarfy$data_loaded)
 ord = calc_ordination(dm, 'nmds')
 plot_ordination(input_fungi_rarfy,ord,'groundcontact_groupcat','groundcontact_treatment',hulls=TRUE)
# this plot shows that all replicates are different and replicate d just has bottom_half treatment
# 
# #BACTERIA: one sample is driving these patterns - what if we dropped it?
 dm = calc_dm(input_bact_rarfy$data_loaded)
 ord = calc_ordination(dm, 'nmds')
 plot_ordination(input_bact_rarfy,ord,'groundcontact_groupcat','groundcontact_treatment',hulls=TRUE)
# #group c, suspended is the block we need to drop to test this
 input_bact_rarfy$map_loaded$treatment_by_group <- paste(input_bact_rarfy$map_loaded$groundcontact_treatment,
#                                                         input_bact_rarfy$map_loaded$groundcontact_groupcat,sep="_")
 input_bact_rarfy$map_loaded$treatment_by_group,
# 
#filter out this value

 trim_for_dm <- filter_data(input_bact_rarfy,filter_cat = 'treatment_by_group',filter_vals = "suspended_c"),
 dm = calc_dm(trim_for_dm$data_loaded))
 ord = calc_ordination(dm, 'nmds')
# # "The goal of NMDS is to represent the original position of communities in
# #...multidimensional space as accurately as possible using a reduced number of 
# #... dimensions that can be easily plotted and visualized"
 plot_ordination(trim_for_dm,ord,'groundcontact_treatment','groundcontact_treatment',hulls=TRUE)
# #now we see that downed and bottom_half bacterial communities are more similar; and top_half and suspended
# # ...are more similar between them
# plot_ordination(trim_for_dm,ord,'groundcontact_groupcat','groundcontact_treatment',hulls=TRUE)
# #in this way we can see that groupcat has an effect on the bacterial community composition since 
# #... the circles are not completely overlapped. I wonder if the similarity they show in 
# #...the graphic would be due to their nearness during the experiment.

#BACTERIA:
input_bact_rarfy$taxonomy_loaded
#now add the taxonomy ("taxonomy_loaded") to the bottom of the OTU table ("data_loaded")
#create a taxonomy file
TAX_table<-as.data.frame(input_bact_rarfy$taxonomy_loaded)
taxonomy<-paste(TAX_table$taxonomy1,TAX_table$taxonomy2,TAX_table$taxonomy3,TAX_table$taxonomy4,TAX_table$taxonomy5,TAX_table$taxonomy6,TAX_table$taxonomy7, sep = "; ")
#add the taxonomy file to the OTU table
rarefied3199_bacteria<-data.frame(input_bact_rarfy$data_loaded,taxonomy)
str(rarefied3199_bacteria)

#FUNGI:
input_fungi_rarfy$taxonomy_loaded
#now add the taxonomy ("taxonomy_loaded") to the bottom of the OTU table ("data_loaded")
#create a taxonomy file
TAX_tablefungi<-as.data.frame(input_fungi_rarfy$taxonomy_loaded)
taxonomyfungi<-paste(TAX_tablefungi$taxonomy1,TAX_tablefungi$taxonomy2,TAX_tablefungi$taxonomy3,TAX_tablefungi$taxonomy4,TAX_tablefungi$taxonomy5,TAX_tablefungi$taxonomy6,TAX_tablefungi$taxonomy7, sep = "; ")
#add the taxonomy file to the OTU table
rarefied4704_fungi<-data.frame(input_fungi_rarfy$data_loaded,taxonomyfungi)
str(rarefied4704_fungi)

# 
# # BACTERIA: 
# write.csv(rarefied3199_bacteria,"rarefied_bacteria.csv")
# 
# # getwd()
# #utilizo getwd() cuando quiero saber el directorio d?nde est? guardado los datos
# 
# # FUNGI: 
# write.csv(rarefied4704_fungi,"rarefied_fungi.csv")
# # getwd()
# 
# #BACTERIA:
# rarefied3199_samplesonly<-rarefied3199_bacteria[,1:20]
# 
# rarefied3199_transposed<-t(rarefied3199_samplesonly)
# str(rarefied3199_transposed)
# 
# #FUNGI:
# rarefied4704_fungi_samplesonly<-rarefied4704_fungi[,1:11]
# 
# rarefied4704_fungi_transposed<-t(rarefied4704_fungi_samplesonly)
# str(rarefied4704_fungi_transposed)
# 
# # #save and export the transposed data for analysis in PRIMER
# # write.csv(rarefied4704_fungi_transposed,"transplant_fungi_OTUmatrix.csv")
# getwd()
# 
 #### PHYLA ####
# #BACTERIA: summarize data by phyla for all data (phyla= level 2)
 taxa_sum_phyla<-summarize_taxonomy(input_bact_rarfy,level=2,report_higher_tax = FALSE)
plot_taxa_bars(taxa_sum_phyla,input_bact_rarfy$map_loaded,'groundcontact_treatment',10)
 plot_ts_heatmap(taxa_sum_phyla,input_bact_rarfy$map_loaded,.005,'groundcontact_treatment')
 plot_taxa_bars(taxa_sum_phyla,input_bact_rarfy$map_loaded,'gc_connection', 10)

# taxa_sum_phyla<-t(taxa_sum_phyla)
# write.csv(taxa_sum_phyla,"summarized_phyla_bacteria.csv")
# 
# # heatmap helps me to see things more clearly, with numbers 
# # ...of the abundance, but just to visualize bars are fine
# 
# #FUNGI: summarize data by phyla for all data (phyla= level 2)
 taxa_sum_phyla_fungi<-summarize_taxonomy(input_fungi_rarfy,level=2,report_higher_tax = FALSE)
 plot_taxa_bars(taxa_sum_phyla_fungi,input_fungi_rarfy$map_loaded,'groundcontact_treatment',10)
# 
 plot_taxa_bars(taxa_sum_phyla_fungi,input_fungi_rarfy$map_loaded,'gc_connection', 10)
 plot_taxa_bars(taxa_sum_phyla_fungi,input_fungi_rarfy$map_loaded,'gc_suspended', 10)
# 
# taxa_sum_phyla_fungi<-t(taxa_sum_phyla_fungi)
# write.csv(taxa_sum_phyla_fungi,"summarized_phyla _fungi.csv")
# 
# #Now we plot the same information in a heatmap, the number .00005 shows the min_rel_abund 
# # that I want to see in the graphic, below that doesnt show:
# plot_ts_heatmap(taxa_sum_phyla_fungi,input_fungi_rarfy$map_loaded,.00005,'groundcontact_treatment')
# 
# #BACTERIA. We see:
# #Relative abundance
# # //Acidobacteria// are more abundant in downed and bottom half:
# #//Actinobacteria// more abundant in top half and suspended
# #//Bacteriodetes// more abundant in top half and suspended
# #//Chloroflexi// common in abundance to all
# #//Cyanobacteria// more abundant in top half and suspended
# #//Firmicutes// more abundant in downed and bottom half
# #//Gemmatimonadetes// quiet more abundant in downed and bottom half
# #//Planctomycetes// apparently common in abundance to all
# #//Proteobacteria// apparently very abundant in all groups 
# #//Verrucomicrobia// apparently common in abundance to all
# 
#### CLASS ####
# #BACTERIA: summarize data by class for all data (class= level 3)
 taxa_sum_class<-summarize_taxonomy(input_bact_rarfy,level=3,report_higher_tax = FALSE)
 plot_taxa_bars(taxa_sum_class,input_bact_rarfy$map_loaded,'groundcontact_treatment',10)
 plot_ts_heatmap(taxa_sum_class,input_bact_rarfy$map_loaded,.01,'groundcontact_treatment')
# 
# #FUNGI: summarize data by class for all data (class= level 3)
 taxa_sum_class_fungi<-summarize_taxonomy(input_fungi_rarfy,level=3,report_higher_tax = FALSE)
  plot_taxa_bars(taxa_sum_class_fungi,input_fungi_rarfy$map_loaded,'groundcontact_treatment',10)
 plot_ts_heatmap(taxa_sum_class_fungi,input_fungi_rarfy$map_loaded,.01,'groundcontact_treatment')
 plot_ts_heatmap(taxa_sum_class_fungi,input_fungi_rarfy$map_loaded,.01,'gc_suspended')
 plot_ts_heatmap(taxa_sum_class_fungi,input_fungi_rarfy$map_loaded,.01,'gc_connection')
# 
 #### ORDER ####
# #BACTERIA:summarize data by order for all data
 taxa_sum_order<-summarize_taxonomy(input_bact_rarfy,level=4,report_higher_tax = FALSE)
 plot_taxa_bars(taxa_sum_order,input_bact_rarfy$map_loaded,'groundcontact_treatment',10)
 plot_ts_heatmap(taxa_sum_order,input_bact_rarfy$map_loaded,.01,'groundcontact_treatment')
 plot_ts_heatmap(taxa_sum_order,input_bact_rarfy$map_loaded,.01,'gc_suspended')
 plot_ts_heatmap(taxa_sum_order,input_bact_rarfy$map_loaded,.01,'gc_connection')
 taxa_sum_order<-t(taxa_sum_order)
write.csv(taxa_sum_order,"summarized_orders_bacteria.csv")

# 
# #FUNGI:summarize data by order for all data
 taxa_sum_order_fungi<-summarize_taxonomy(input_fungi_rarfy,level=4,report_higher_tax = FALSE)
 plot_taxa_bars(taxa_sum_order_fungi,input_fungi_rarfy$map_loaded,'groundcontact_treatment',10)
 plot_ts_heatmap(taxa_sum_order_fungi,input_fungi_rarfy$map_loaded,.01,'groundcontact_treatment')
 #now for different treatments: 

# #plot_taxa_bars(taxa_sum_order_fungi,input_fungi_rarfy$map_loaded,.01,'gc_suspended')
plot_ts_heatmap(taxa_sum_order_fungi,input_fungi_rarfy$map_loaded,.03,'gc_connection')
all_taxa_sum_order_fungi<-t(taxa_sum_order_fungi)
write.csv(all_taxa_sum_order_fungi,"summarized_orders_fungi.csv")
 
# 
 #### FAMILY ####
# #BACTERIA: summarize data by family for all data
# taxa_sum_family<-summarize_taxonomy(input_bact_rarfy,level=5,report_higher_tax = FALSE)
# plot_taxa_bars(taxa_sum_family,input_bact_rarfy$map_loaded,'groundcontact_treatment',10)
# plot_ts_heatmap(taxa_sum_family,input_bact_rarfy$map_loaded,.02,'groundcontact_treatment')
# 
# #FUNGI: summarize data by family for all data
# taxa_sum_family_fungi<-summarize_taxonomy(input_fungi_rarfy,level=5,report_higher_tax = FALSE)
# plot_taxa_bars(taxa_sum_family_fungi,input_fungi_rarfy$map_loaded,'groundcontact_treatment',10)
# plot_ts_heatmap(taxa_sum_family_fungi,input_fungi_rarfy$map_loaded,.02,'groundcontact_treatment')
# 
# #BACTERIA: you can also examine subsets of the blocks alone (connected and not connected)
# #subset the data to get
separated<-filter_data(input_bact_rarfy, filter_cat = "gc_connection", filter_vals = c("connected"))
connected<-filter_data(input_bact_rarfy, filter_cat = "gc_connection", filter_vals = c("separated"))
# 
# #FUNGI: you can also examine subsets of the blocks alone
# #subset the data to get
separated_fungi<-filter_data(input_fungi_rarfy, filter_cat = "gc_connection", filter_vals = c("connected"))
connected_fungi<-filter_data(input_fungi_rarfy, filter_cat = "gc_connection", filter_vals = c("separated"))
# 
# #BACTERIA: Canopy to ground
# dm = calc_dm(separated$data_loaded)
# ord = calc_ordination(dm, 'nmds')
# plot_ordination(separated,ord,'groundcontact_treatment',hulls=TRUE)
# 
# #FUNGI: Canopy to ground
# dm = calc_dm(separated_fungi$data_loaded)
# ord = calc_ordination(dm, 'nmds')
# plot_ordination(separated_fungi,ord,'groundcontact_treatment',hulls=TRUE)
# 
# #BACTERIA: one sample is driving these patterns - what if we dropped it?
# 
# #group c, suspended is the block we need to drop to test this
# separated$map_loaded$treatment_by_group <- paste(separated$map_loaded$groundcontact_treatment,
#                                                  separated$map_loaded$groundcontact_groupcat,sep="_")
# separated$map_loaded$treatment_by_group
# 
#filter out this value
# trim_for_dm <- filter_data(separated,filter_cat = 'treatment_by_group',filter_vals = "suspended_c")
# dm = calc_dm(trim_for_dm$data_loaded)
# ord = calc_ordination(dm, 'nmds')
# # "The goal of NMDS is to represent the original position of communities in
# #...multidimensional space as accurately as possible using a reduced number of 
# #... dimensions that can be easily plotted and visualized"
# plot_ordination(trim_for_dm,ord,'groundcontact_treatment','groundcontact_treatment',hulls=TRUE)
# #now we see that downed and bottom_half bacterial communities are more similar; and top_half and suspended
# # ...are more similar between them
# plot_ordination(trim_for_dm,ord,'groundcontact_groupcat','groundcontact_treatment',hulls=TRUE)
# #in this way we can see that groupcat has an effect on the bacterial community composition since 
# #... the circles are not completely overlapped. I wonder if the similarity they show in 
# #...the graphic would be due to their nearness during the experiment.
# 
# 
# #BACTERIA: ground to canopy - the connection really does facilitate colonization
# dm = calc_dm(connected$data_loaded)
# ord = calc_ordination(dm, 'nmds')
# plot_ordination(connected,ord,'groundcontact_treatment',hulls=TRUE)
# plot_ordination(connected,ord,'groundcontact_groupcat',hulls=TRUE)
# 
# #FUNGI: ground to canopy - the CONNECTION really does facilitate colonization
# dm = calc_dm(connected_fungi$data_loaded)
# ord = calc_ordination(dm, 'nmds')
# plot_ordination(connected_fungi,ord,'groundcontact_treatment',hulls=TRUE)
# plot_ordination(connected_fungi,ord,'groundcontact_groupcat',hulls=TRUE)

# #In vegan::metaMDS(dm, k = 2) :
# #  stress is (nearly) zero: you may have insufficient data
# 
# # BACTERIA
#### Look at the main taxa for the source environment ####
# #now do look at the main taxa for the source environment
# #summarize data by phyla for source environment
# taxa_sum_phyla<-summarize_taxonomy(separated,level=2,report_higher_tax = FALSE)
# plot_taxa_bars(taxa_sum_phyla,separated$map_loaded,'groundcontact_treatment',10)
# plot_ts_heatmap(taxa_sum_phyla,separated$map_loaded,.01,'groundcontact_treatment')
# 
# #summarize data by class for source environment
# taxa_sum_class<-summarize_taxonomy(separated,level=3,report_higher_tax = FALSE)
# plot_taxa_bars(taxa_sum_class,separated$map_loaded,'groundcontact_treatment',10)
# plot_ts_heatmap(taxa_sum_class,separated$map_loaded,.01,'groundcontact_treatment')
# 
# #summarize data by order for source environment
# taxa_sum_order<-summarize_taxonomy(separated,level=4,report_higher_tax = FALSE)
# plot_taxa_bars(taxa_sum_order,separated$map_loaded,'groundcontact_treatment',10)
# plot_ts_heatmap(taxa_sum_order,separated$map_loaded,.02,'groundcontact_treatment')
# pretreat_taxa_sum_order<-t(taxa_sum_order)
# write.csv(pretreat_taxa_sum_order,"bact_order_summary_pretreat.csv")
# 
# 
# #summarize data by family for source environment
# taxa_sum_family<-summarize_taxonomy(separated,level=5,report_higher_tax = FALSE)
# plot_taxa_bars(taxa_sum_family,separated$map_loaded,'groundcontact_treatment',10)
# plot_ts_heatmap(taxa_sum_family,separated$map_loaded,.02,'groundcontact_treatment')



# 
# ### FUNGI (!!!!!!!)
# #now do look at the main taxa for the source environment
# #summarize data by phyla for source environment
# taxa_sum_phyla_fungi<-summarize_taxonomy(input_fungi_rarfy,level=2,report_higher_tax = FALSE)
# plot_taxa_bars(taxa_sum_phyla_fungi,separated$map_loaded,'groundcontact_treatment',10)

# plot_ts_heatmap(taxa_sum_phyla_fungi,separated$map_loaded,.01,'groundcontact_treatment')
# #Error in tapply(as.numeric(x), metadata_map[, type_header], mean) : 
# #  arguments must have same length
# 
# #summarize data by class for source environment
# taxa_sum_class_fungi<-summarize_taxonomy(input_fungi_rarfy,level=3,report_higher_tax = FALSE)
# plot_taxa_bars(taxa_sum_class_fungi,separated$map_loaded,'groundcontact_treatment',10)

# plot_ts_heatmap(taxa_sum_class_fungi,separated$map_loaded,.01,'groundcontact_treatment')
# #Error in tapply(as.numeric(x), metadata_map[, type_header], mean) : 
# #  arguments must have same length
# 
# 
#summarize data by order for source environment
taxa_sum_order_fungi<-summarize_taxonomy(input_fungi_rarfy,level=4,report_higher_tax = FALSE)
plot_taxa_bars(taxa_sum_order_fungi,separated$map_loaded,'groundcontact_treatment',10)

plot_ts_heatmap(taxa_sum_order_fungi,separated$map_loaded,.02,'groundcontact_treatment')
#Error in tapply(as.numeric(x), metadata_map[, type_header], mean) :
#  arguments must have same length
pretreat_taxa_sum_order_fungi<-t(taxa_sum_order_fungi)
# write.csv(pretreat_taxa_sum_order_fungi,"fungi_order_summary_pretreat.csv")
# 
# 
# #summarize data by family for source environment
# taxa_sum_family<-summarize_taxonomy(separated,level=5,report_higher_tax = FALSE)
# plot_taxa_bars(taxa_sum_family,separated$map_loaded,'groundcontact_treatment',10)
# plot_ts_heatmap(taxa_sum_family,separated$map_loaded,.02,'groundcontact_treatment')
# 
# #BACTERIA: now do look at the main taxa for the source environment
# #summarize data by phyla for source environment
# taxa_sum_phyla<-summarize_taxonomy(connected,level=2,report_higher_tax = FALSE)
# plot_taxa_bars(taxa_sum_phyla,connected$map_loaded,'groundcontact_treatment',10)
# plot_ts_heatmap(taxa_sum_phyla,connected$map_loaded,.01,'groundcontact_treatment')
# 
# #FUNGI: now do look at the main taxa for the source environment
# #summarize data by phyla for source environment
# taxa_sum_phyla_fungi1<-summarize_taxonomy(connected_fungi,level=2,report_higher_tax = FALSE)
# plot_taxa_bars(taxa_sum_phyla_fungi1,connected_fungi$map_loaded,'groundcontact_treatment',10)
# plot_ts_heatmap(taxa_sum_phyla_fungi1,connected_fungi$map_loaded,.00005,'groundcontact_treatment')
# 
# #BACTERIA: summarize data by class for source environment
# taxa_sum_class<-summarize_taxonomy(connected,level=3,report_higher_tax = FALSE)
# plot_taxa_bars(taxa_sum_class,connected$map_loaded,'groundcontact_treatment',10)
# plot_ts_heatmap(taxa_sum_class,connected$map_loaded,.01,'groundcontact_treatment')
# 
# #FUNGI: summarize data by class for source environment
# taxa_sum_class_fungi1<-summarize_taxonomy(connected_fungi,level=3,report_higher_tax = FALSE)
# plot_taxa_bars(taxa_sum_class_fungi1,connected_fungi$map_loaded,'groundcontact_treatment',10)
# plot_ts_heatmap(taxa_sum_class_fungi1,connected_fungi$map_loaded,.0005,'groundcontact_treatment')
# 
# #BACTERIA: summarize data by order for source environment
 taxa_sum_order<-summarize_taxonomy(connected,level=4,report_higher_tax = FALSE)
 plot_taxa_bars(taxa_sum_order,connected$map_loaded,'groundcontact_treatment',10)
 plot_ts_heatmap(taxa_sum_order,connected$map_loaded,.02,'groundcontact_treatment')
 posttreat_taxa_sum_order<-t(taxa_sum_order)
# write.csv(posttreat_taxa_sum_order,"bact_order_summary_posttreat.csv")
# 
# #FUNGI: summarize data by order for source environment
# taxa_sum_order_fungi1<-summarize_taxonomy(connected_fungi,level=4,report_higher_tax = FALSE)
# plot_taxa_bars(taxa_sum_order_fungi1,connected_fungi$map_loaded,'groundcontact_treatment',10)
# plot_ts_heatmap(taxa_sum_order_fungi1,connected_fungi$map_loaded,.02,'groundcontact_treatment')
# 
# posttreat_taxa_sum_order_fungi<-t(taxa_sum_order_fungi1)
# write.csv(posttreat_taxa_sum_order_fungi,"fungi_order_summary_posttreat.csv")
# 
# 
# #BACTERIA: summarize data by family for source environment
# taxa_sum_family<-summarize_taxonomy(connected,level=5,report_higher_tax = FALSE)
# plot_taxa_bars(taxa_sum_family,connected$map_loaded,'groundcontact_treatment',15)
# plot_ts_heatmap(taxa_sum_family,connected$map_loaded,.02,'groundcontact_treatment')
# 
# #FUNGI: summarize data by family for source environment
# taxa_sum_family_fungi<-summarize_taxonomy(connected_fungi,level=5,report_higher_tax = FALSE)
# plot_taxa_bars(taxa_sum_family_fungi,connected_fungi$map_loaded,'groundcontact_treatment',10)
# plot_ts_heatmap(taxa_sum_family_fungi,connected_fungi$map_loaded,.02,'groundcontact_treatment')
# 

#///////////////////////////////////
### 2. ANALYZE OVERALL DIVERSITY: Remember to check if the bad sample is already droped for this analysis #### 
#///////////////////////////////////

###BACTERIA: extract data for additional analyses
#estimate richness and shannon diversity

sum_rich_nosingles_bact<-as.numeric(calc_diversity(input_bact_rarfy$data_loaded,"richness"))
list(sum_rich_nosingles_bact)
sort(sum_rich_nosingles_bact, decreasing = FALSE)

sum_shan_nosingles_bact<-as.numeric(calc_diversity(input_bact_rarfy$data_loaded,"shannon"))
list(sum_shan_nosingles_bact)
sort(sum_shan_nosingles_bact, decreasing = FALSE)

sum_simp_nosingles_bact<-as.numeric(calc_diversity(input_bact_rarfy$data_loaded,"simpson"))
list(sum_simp_nosingles_bact)
sort(sum_simp_nosingles_bact, decreasing = FALSE)

###FUNGI: extract data for additional analyses
#estimate richness and shannon diversity

sum_rich_nosingles_fungi<-as.numeric(calc_diversity(input_fungi_rarfy$data_loaded,"richness"))
list(sum_rich_nosingles_fungi)
sort(sum_rich_nosingles_fungi, decreasing = FALSE)

sum_shan_nosingles_fungi<-as.numeric(calc_diversity(input_fungi_rarfy$data_loaded,"shannon"))
list(sum_shan_nosingles_fungi)
sort(sum_shan_nosingles_fungi, decreasing = FALSE)

sum_simp_nosingles_fungi<-as.numeric(calc_diversity(input_fungi_rarfy$data_loaded,"simpson"))
list(sum_simp_nosingles_fungi)
sort(sum_simp_nosingles_fungi, decreasing = FALSE)

#BACTERIA
mapping_bactdiversity <- input_bact_rarfy$map_loaded

mapping_bactdiversity$richness <- sum_rich_nosingles_bact
mapping_bactdiversity$shannon <- sum_shan_nosingles_bact
mapping_bactdiversity$simpson <- sum_simp_nosingles_bact

write.csv(mapping_bactdiversity, "mapping_bactdiversity.csv")
getwd()

#FUNGI
mapping_fungidiversity <- input_fungi_rarfy$map_loaded

mapping_fungidiversity$richness <- sum_rich_nosingles_fungi
mapping_fungidiversity$shannon <- sum_shan_nosingles_fungi
mapping_fungidiversity$simpson <- sum_simp_nosingles_fungi

write.csv(mapping_fungidiversity, "mapping_fungidiversity.csv")
getwd()


### Mean richness, shannon, simpson: Bact and fungi ####

#bact each treatment richness (4)

bact_richness_se_sd_mean <- ddply(mapping_bactdiversity, 
                                       .(gc_connection, gc_suspended), 
                                       summarise, 
                                       M = mean(richness), 
                                       SE = sd(richness) / sqrt((length(richness))), 
                                       SD = sd(richness), 
                                       Median=median(richness))

#each treatment shannon (4)

bact_shannon_se_sd_mean <- ddply(mapping_bactdiversity, 
                                  .(gc_connection, gc_suspended), 
                                  summarise, 
                                  M = mean(shannon), 
                                  SE = sd(shannon) / sqrt((length(shannon))), 
                                  SD = sd(shannon), 
                                  Median=median(shannon))

#each treatment simpson (4)

bact_simpson_se_sd_mean <- ddply(mapping_bactdiversity, 
                                 .(gc_connection, gc_suspended), 
                                 summarise, 
                                 M = mean(simpson), 
                                 SE = sd(simpson) / sqrt((length(simpson))), 
                                 SD = sd(simpson), 
                                 Median=median(simpson))
-----------
#FUNGI each treatment richness (4)

fung_richness_se_sd_mean <- ddply(mapping_fungidiversity, 
                                  .(gc_connection, gc_suspended), 
                                  summarise, 
                                  M = mean(richness), 
                                  SE = sd(richness) / sqrt((length(richness))), 
                                  SD = sd(richness), 
                                  Median=median(richness))

#each treatment shannon (4)

fung_shannon_se_sd_mean <- ddply(mapping_fungidiversity, 
                                 .(gc_connection, gc_suspended), 
                                 summarise, 
                                 M = mean(shannon), 
                                 SE = sd(shannon) / sqrt((length(shannon))), 
                                 SD = sd(shannon), 
                                 Median=median(shannon))

#each treatment simpson (4)

fung_simpson_se_sd_mean <- ddply(mapping_fungidiversity, 
                                 .(gc_connection, gc_suspended), 
                                 summarise, 
                                 M = mean(simpson), 
                                 SE = sd(simpson) / sqrt((length(simpson))), 
                                 SD = sd(simpson), 
                                 Median=median(simpson))


#now extract the corresponding treatments and run anovas comparing the diversity metrics

###BACTERIA: 

names(input_bact_rarfy$data_loaded)
rep_bact <- as.data.frame(names(input_bact_rarfy$data_loaded))
rep_bact$ground_contact_treatment <- as.character(input_bact_rarfy$map_loaded$groundcontact_treatment)
rep_bact$sum_rich_nosingles_bact <- sum_rich_nosingles_bact
rep_bact$gc_connection <- as.character(input_bact_rarfy$map_loaded$gc_connection)
rep_bact$gc_suspended <- as.character(input_bact_rarfy$map_loaded$gc_suspended)

#we need to change ground_contact_treatment from character to factor to run the ANOVA

rep_bact$ground_contact_treatment <- as.factor(rep_bact$ground_contact_treatment)


#ANOVA RICHNESS BACT
rich_bact_mod <- lm(sum_rich_nosingles_bact ~ gc_connection + gc_suspended + gc_connection:gc_suspended,
                    data = rep_bact)
Anova(rich_bact_mod)

rich_bact_mod2 <- lm(sum_rich_nosingles_bact ~ gc_connection + gc_suspended,
                    data = rep_bact)
Anova(rich_bact_mod2)


#ANOVA SHANNON BACT
shan_bact_mod <- lm(sum_shan_nosingles_bact ~ gc_connection + gc_suspended + gc_connection:gc_suspended,
                    data = rep_bact)
Anova(shan_bact_mod)

shan_bact_mod2 <- lm(sum_shan_nosingles_bact ~ gc_connection + gc_suspended,
                     data = rep_bact)
Anova(shan_bact_mod2)


#ANOVA SIMPSON BACT
simp_bact_mod <- lm(sum_simp_nosingles_bact ~ gc_connection + gc_suspended + gc_connection:gc_suspended,
                    data = rep_bact)
Anova(simp_bact_mod)

simp_bact_mod2 <- lm(sum_simp_nosingles_bact ~ gc_connection + gc_suspended,
                     data = rep_bact)
Anova(simp_bact_mod2)


###FUNGI: 
names(input_fungi_rarfy$data_loaded)
rep_fungi <- as.data.frame(names(input_fungi_rarfy$data_loaded))
rep_fungi$ground_contact_treatment <- as.character(input_fungi_rarfy$map_loaded$groundcontact_treatment)
rep_fungi$sum_rich_nosingles_fungi <- sum_rich_nosingles_fungi
rep_fungi$gc_connection <- as.character(input_fungi_rarfy$map_loaded$gc_connection)
rep_fungi$gc_suspended <- as.character(input_fungi_rarfy$map_loaded$gc_suspended)

#we need to change ground_contact_treatment from character to factor to run the ANOVA

rep_fungi$ground_contact_treatment <- as.factor(rep_fungi$ground_contact_treatment)


#ANOVA RICHNESS fungi
rich_fungi_mod <- lm(sum_rich_nosingles_fungi ~ gc_connection + gc_suspended + gc_connection:gc_suspended,
                    data = rep_fungi)
Anova(rich_fungi_mod)

rich_fungi_mod2 <- lm(sum_rich_nosingles_fungi ~ gc_connection + gc_suspended,
                     data = rep_fungi)
Anova(rich_fungi_mod2)


#ANOVA SHANNON fungi
shan_fungi_mod <- lm(sum_shan_nosingles_fungi ~ gc_connection + gc_suspended + gc_connection:gc_suspended,
                    data = rep_fungi)
Anova(shan_fungi_mod)

shan_fungi_mod2 <- lm(sum_shan_nosingles_fungi ~ gc_connection + gc_suspended,
                     data = rep_fungi)
Anova(shan_fungi_mod2)


#ANOVA SIMPSON fungi
simp_fungi_mod <- lm(sum_simp_nosingles_fungi ~ gc_connection + gc_suspended + gc_connection:gc_suspended,
                    data = rep_fungi)
Anova(simp_fungi_mod)

simp_fungi_mod2 <- lm(sum_simp_nosingles_fungi ~ gc_connection + gc_suspended,
                     data = rep_fungi)
Anova(simp_fungi_mod2)


# # # # BACT # # #
# #transponer pero me generara una matriz
 data_loadedbact <- t(input_bact_rarfy$data_loaded)
 data_loadedbact <- as.data.frame(data_loadedbact)
 #repeticiones no es una columna sino que es el nombre de las filas
 data_loadedbact$samples <- as.factor(rownames(data_loadedbact))
 data_loadedbact <- data_loadedbact[,c(8282, 1:8281)]
# 
 mapdiversitybact <- cbind(data_loadedbact, 
                           input_bact_rarfy$map_loaded$groundcontact_treatment,
                           input_bact_rarfy$map_loaded$gc_connection,
                           input_bact_rarfy$map_loaded$gc_suspended)
 mapdiversitybact <- mapdiversitybact[,c(1,8283:8285, 2:8282)]
 mapdiversitybact <- mapdiversitybact[,-c(5:8292)]
 mapdiversitybact$richness <- sum_rich_nosingles_bact
 mapdiversitybact$shannon <- sum_shan_nosingles_bact
 mapdiversitybact$simpson <- sum_simp_nosingles_bact
 #Vamos a cambiar los nombres
 colnames(mapdiversitybact) <- c("samples", "treatment", "gc_connection", "gc_suspended",
                                 "richness", "shannon", "simpson")
 mapdiversitybact$treatment <- as.factor(mapdiversitybact$treatment)
 mapdiversitybact$gc_connection <- as.factor(mapdiversitybact$gc_connection)
 mapdiversitybact$gc_suspended <- as.factor(mapdiversitybact$gc_suspended)
# 
# #ANOVAS de BACT para treatment, gc_connection, gc_suspended x riqueza, 
# # shannon, simpson (9!!)
# 
# #para ANOVAS, recordar poner primero lo numerico (richness) y luego la 
# #categoria (treatment)
# 
# anova_bact_gct_rich <- aov(richness ~ treatment, mapdiversitybact)
# summary(anova_bact_gct_rich) # p value 0.466
# 
# anova_bact_gct_shan <- aov(shannon ~ treatment, mapdiversitybact)
# summary(anova_bact_gct_shan) # p value 0.564
# 
# anova_bact_gct_simp <- aov(simpson ~ treatment, mapdiversitybact)
# summary(anova_bact_gct_simp) # p value 0.73
# 
# ...
# 
# anova_bact_conn_rich <- aov(richness ~ gc_connection, mapdiversitybact)
# summary(anova_bact_conn_rich) # p value 0.915
# 
# anova_bact_conn_shan <- aov(shannon ~ gc_connection, mapdiversitybact)
# summary(anova_bact_conn_shan) # p value 0.814
# 
# anova_bact_conn_simp <- aov(simpson ~ gc_connection, mapdiversitybact)
# summary(anova_bact_conn_simp) # p value 0.77
# 
# #RESULTS (NO SIGNIFICATIVE DIFERENCES) 
# ...
# 
# anova_bact_sus_rich <- aov(richness ~ gc_suspended, mapdiversitybact)
# summary(anova_bact_sus_rich) # p value 0.937
# 
# anova_bact_sus_shan <- aov(shannon ~ gc_suspended, mapdiversitybact)
# summary(anova_bact_sus_shan) # p value 0.982
# 
# anova_bact_sus_simp <- aov(simpson ~ gc_suspended, mapdiversitybact)
# summary(anova_bact_sus_simp) # p value 0.886

# #RESULTS (NO SIGNIFICATIVE DIFERENCES) 
# 
# # # # FUNGI # # #
# #transponer pero me generara una matriz
# data_loadedfungi <- t(input_fungi_rarfy$data_loaded)
# data_loadedfungi <- as.data.frame(data_loadedfungi)
# #repeticiones no es una columna sino que es el nombre de las filas
# data_loadedfungi$samples <- as.factor(rownames(data_loadedfungi))
# data_loadedfungi <- data_loadedfungi[,c(1807, 1:1806)]
# 
# 
# mapdiversityfungi <- cbind(data_loadedfungi, 
#                            input_fungi_rarfy$map_loaded$groundcontact_treatment,
#                            input_fungi_rarfy$map_loaded$gc_connection,
#                            input_fungi_rarfy$map_loaded$gc_suspended)
# 
# mapdiversityfungi <- mapdiversityfungi[,c(1,1808:1810,4:1807)]
# mapdiversityfungi <- mapdiversityfungi[,-c(5:1810)]
# mapdiversityfungi$richness <- sum_rich_nosingles_fungi
# mapdiversityfungi$shannon <- sum_shan_nosingles_fungi
# mapdiversityfungi$simpson <- sum_simp_nosingles_fungi
# 
# 
# #Vamos a cambiar los nombres
# colnames(mapdiversityfungi) <- c("samples", "treatment", "gc_connection", "gc_suspended",
#                                  "richness", "shannon", "simpson")
# mapdiversityfungi$treatment <- as.factor(mapdiversityfungi$treatment)
# mapdiversityfungi$gc_connection <- as.factor(mapdiversityfungi$gc_connection)
# mapdiversityfungi$gc_suspended <- as.factor(mapdiversityfungi$gc_suspended)
# 
# #ANOVAS de FUNG para treatment, gc_connection, gc_suspended x riqueza, 
# # shannon, simpson (9!!)
# 
# #para ANOVAS, recordar poner primero lo numerico (richness) y luego la 
# #categoria (treatment)
# 
# anova_fungi_gct_rich <- aov(richness ~ treatment, mapdiversityfungi)
# summary(anova_fungi_gct_rich) # p value 0.293
# 
# anova_fungi_gct_shan <- aov(shannon ~ treatment, mapdiversityfungi)
# summary(anova_fungi_gct_shan) # p value 0.394
# 
# anova_fungi_gct_simp <- aov(simpson ~ treatment, mapdiversityfungi)
# summary(anova_fungi_gct_simp) # p value 0.334
# 
# #RESULTS (NO SIGNIFICATIVE DIFERENCES)
# ...
# 
# anova_fungi_conn_rich <- aov(richness ~ gc_connection, mapdiversityfungi)
# summary(anova_fungi_conn_rich) # p value 0.743
# 
# anova_fungi_conn_shan <- aov(shannon ~ gc_connection, mapdiversityfungi)
# summary(anova_fungi_conn_shan) # p value 0.325
# 
# anova_fungi_conn_simp<- aov(simpson ~ gc_connection, mapdiversityfungi)
# summary(anova_fungi_conn_simp) # p value 0.284
# 
# ...
# 
# anova_fungi_sus_rich <- aov(richness ~ gc_suspended, mapdiversityfungi)
# summary(anova_fungi_sus_rich) # p value 0.21
# 
# anova_fungi_sus_shan <- aov(shannon ~ gc_suspended, mapdiversityfungi)
# summary(anova_fungi_sus_shan) # p value 0.0836
# 
# anova_fungi_sus_simp <- aov(simpson ~ gc_suspended, mapdiversityfungi)
# summary(anova_fungi_sus_simp) # p value 0.088
# 
# #RESULTS (NO SIGNIFICATIVE DIFERENCES)

#///////////////////////////////////////
### 3. ANALYZE COMMUNITY COMPOSITION ####
#///////////////////////////////////////




# We will perform an Adonis test. Multilevel pairwise comparisons to compare
#specific pairs of block subsets

#Are the distances (relative abundances of orders in each sample?)between samples from the same group samaller than the distance 
#between samples from different groups?

#Use this website as a reference  https://chrischizinski.github.io/rstats/adonis/
## adonis(formula = all.sites ~ trt, permutations = 999, method = "bray") 

#~~~~~~
#FUNGI#         ORDERS
#~~~~~~

# #1. 
# pretreat_taxa_sum_order_fungi_df  <-  as.data.frame(pretreat_taxa_sum_order_fungi)
# pretreat_taxa_sum_order_fungi_df$samples <- rownames(pretreat_taxa_sum_order_fungi) 
# 
# rep_fungi_adonis <- rep_fungi
# colnames(rep_fungi_adonis)[1] <- "samples"
# rep_fungi_adonis_2c <- rep_fungi_adonis[-3]
# 
# taxa_sum_order_fungi_adonis <- left_join(pretreat_taxa_sum_order_fungi_df, 
#                                          rep_fungi_adonis_2c)
# 
# #2. explanatory variable = "taxa_sum_order_fungi_adonis$groundcontact_treatment"
# 
# # adonis test
# 
# #an easy way to call columns (prelast and last one)
# ncol(taxa_sum_order_fungi_adonis)
# c((ncol(taxa_sum_order_fungi_adonis)-1),ncol(taxa_sum_order_fungi_adonis))
# 
# adonis_order_fungi <- adonis(taxa_sum_order_fungi_adonis[,-c((ncol(taxa_sum_order_fungi_adonis)-1),ncol(taxa_sum_order_fungi_adonis))] ~ 
#                                taxa_sum_order_fungi_adonis$ground_contact_treatment,
#                              permutations = 10000, method = "bray") 
# 
# #RESULTS: NO SIGNIFICANT DIFERENCES (0.2941 > 0.05)
#  
# #~~~~~~
# #FUNGI#         OTUS
# #~~~~~~
# #Use data_loadedfungi instead of pretreat_taxa_sum_order_fungi!!! (for finest scale OTUs)
# 
# data_loadedfungi_df  <-  as.data.frame(data_loadedfungi)
# str(data_loadedfungi_df) #numeric! ok!
# data_loadedfungi_df$samples <- rownames(data_loadedfungi) 
# 
# rep_fungi_adonis_2c 
# 
# taxa_sum_otus_fungi_adonis <- left_join(data_loadedfungi_df, 
#                                         rep_fungi_adonis_2c)
# 
# adonis_otus_fungi <- adonis(taxa_sum_otus_fungi_adonis[,-c((ncol(taxa_sum_otus_fungi_adonis)-1),ncol(taxa_sum_otus_fungi_adonis))] ~ 
#                               taxa_sum_otus_fungi_adonis$ground_contact_treatment,
#                             permutations = 10000, method = "bray") 
# 
# #RESULTS: NO SIGNIFICANT DIFERENCES (0.4555 > 0.05) 
# 
# # Adonis for gc_suspended:
# 
# colnames(rep_fungi_adonis_2c)[2] <- "groundcontact_treatment"
# 
# rep_fungi_adonis_2c_gc <- left_join(rep_fungi_adonis_2c, input_fungi_rarfy$map_loaded , by= "groundcontact_treatment")
# 
# colnames(rep_fungi_adonis_2c_gc)
# 
# rep_fungi_adonis_2c_gc_sc <- rep_fungi_adonis_2c_gc[c("samples","groundcontact_treatment",
#                                                     "gc_connection","gc_suspended","groundcontact_groupcat")]
# 
# taxa_sum_otus_fungi_adonis_gc_sc <- left_join(data_loadedfungi_df, 
#                                              rep_fungi_adonis_2c_gc_sc , by= "samples")
# str(taxa_sum_otus_fungi_adonis_gc_sc)
# 
# adonis_otus_fungi_gc_sus <- adonis(taxa_sum_otus_fungi_adonis_gc_sc[,-c((ncol(taxa_sum_otus_fungi_adonis_gc_sc)-4),
#                                                                       (ncol(taxa_sum_otus_fungi_adonis_gc_sc)-3), 
#                                                                       (ncol(taxa_sum_otus_fungi_adonis_gc_sc)-2),
#                                                                       (ncol(taxa_sum_otus_fungi_adonis_gc_sc)-1),
#                                                                       ncol(taxa_sum_otus_fungi_adonis_gc_sc))] ~ 
#                                     taxa_sum_otus_fungi_adonis_gc_sc$gc_suspended,
#                                   permutations = 10000, method = "bray") 
# 
# #RESULTS: SIGNIFICANT DIFERENCES ( 9.999e-05 < 0.05) !
# #Interpretation: Diversity patterns are being driven by elevation from the ground (up, down)
# 
# # Adonis for gc_connection:
# 
# adonis_otus_fungi_gc_con <- adonis(taxa_sum_otus_fungi_adonis_gc_sc[,-c((ncol(taxa_sum_otus_fungi_adonis_gc_sc)-4),
#                                                                       (ncol(taxa_sum_otus_fungi_adonis_gc_sc)-3), 
#                                                                       (ncol(taxa_sum_otus_fungi_adonis_gc_sc)-2),
#                                                                       (ncol(taxa_sum_otus_fungi_adonis_gc_sc)-1),
#                                                                       ncol(taxa_sum_otus_fungi_adonis_gc_sc))] ~ 
#                                     taxa_sum_otus_fungi_adonis_gc_sc$gc_connection,
#                                   permutations = 10000, method = "bray") 
# 
# 
# #RESULTS: SIGNIFICANT DIFERENCES ( 0.0049 ** < 0.05) !
# #Interpretation: Diversity patterns are being driven by connection with the ground (connected, separated)
# 
# # ~~~~~ Interaction (interaction(gc_suspended + gc_connection + gc_suspended:gc_connection))
# 
# adonis_otus_fungi_gc_consus <- adonis(taxa_sum_otus_fungi_adonis_gc_sc[,-c((ncol(taxa_sum_otus_fungi_adonis_gc_sc)-4),
#                                                                          (ncol(taxa_sum_otus_fungi_adonis_gc_sc)-3), 
#                                                                          (ncol(taxa_sum_otus_fungi_adonis_gc_sc)-2),
#                                                                          (ncol(taxa_sum_otus_fungi_adonis_gc_sc)-1),
#                                                                          ncol(taxa_sum_otus_fungi_adonis_gc_sc))] ~ 
#                                        (taxa_sum_otus_fungi_adonis_gc_sc$gc_connection +
#                                           taxa_sum_otus_fungi_adonis_gc_sc$gc_suspended +
#                                           taxa_sum_otus_fungi_adonis_gc_sc$gc_connection:taxa_sum_otus_fungi_adonis_gc_sc$gc_suspended),
#                                      permutations = 10000, method = "bray") 

#RESULTS: SIGNIFICANT DIFERENCES FOR ALL ( Pr(>F)=0.009899 ** < 0.05) !
#Interpretation: Diversity patterns are being driven by connection and elevation with the ground! Being elevated matters,
# being connected matters and being connected changes the effect of being elevated. 

#~~~~~~~~~
#BACTERIA#           ORDERS
#~~~~~~~~~

#1. 


###THESE MIGHT NOT BE THE DATA WE WANT TO USE

# pretreat_taxa_sum_order_bact_df <-  as.data.frame(pretreat_taxa_sum_order)
# pretreat_taxa_sum_order_bact_df$samples <- rownames(pretreat_taxa_sum_order)
# 
# rep_bact_adonis <- rep_bact
# colnames(rep_bact_adonis)[1] <- "samples"
# rep_bact_adonis_2c <- rep_bact_adonis[-3]
# 
# taxa_sum_order_bact_adonis <- left_join(pretreat_taxa_sum_order_bact_df, 
#                                         rep_bact_adonis_2c)
# 
# str(taxa_sum_order_bact_adonis) #all numeric!
# 
# #2. explanatory variable = "taxa_sum_order_bact_adonis$ground_contact_treatment"
# 
# # adonis test
# 
# adonis_order_bact <- adonis(taxa_sum_order_bact_adonis[,-c((ncol(taxa_sum_order_bact_adonis)-1),ncol(taxa_sum_order_bact_adonis))] ~ 
#                               taxa_sum_order_bact_adonis$ground_contact_treatment,
#                             permutations = 10000, method = "bray") 
# 
# #RESULTS: SIGNIFICANT DIFERENCES (0.0386 > 0.05) !
# 
# #~~~~~~
# #BACTERIA#         OTUS
# #~~~~~~
# #Use data_loaded instead of pretreat_taxa_sum_order_bact!!! (for finest scale OTUs)
# 
# data_loadedbact_df  <-  as.data.frame(data_loadedbact)
# str(data_loadedbact_df) #numeric! ok!
# data_loadedbact_df$samples <- rownames(data_loadedbact) 
# 
# rep_bact_adonis_2c 
# 
# taxa_sum_otus_bact_adonis <- left_join(data_loadedbact_df, 
#                                        rep_bact_adonis_2c)
# 
# adonis_otus_bact <- adonis(taxa_sum_otus_bact_adonis[,-c((ncol(taxa_sum_otus_bact_adonis)-1),ncol(taxa_sum_otus_bact_adonis))] ~ 
#                              taxa_sum_otus_bact_adonis$ground_contact_treatment,
#                            permutations = 10000, method = "bray") 
# 
# #RESULTS: Poor SIGNIFICANT DIFERENCES (0.08719 > 0.05) !
# 
# #//////////////////////
# #Now we want to see the interaction so we run more adonis for gc_suspended and gc_connected 
# #and the interaction(gc_suspended + gc_connection + gc_suspended:gc_connection)
# #specify the df of the above 
# 
# 
# # Adonis for gc_suspended:
# 
# colnames(rep_bact_adonis_2c)[2] <- "groundcontact_treatment"
# 
# rep_bact_adonis_2c_gc <- left_join(rep_bact_adonis_2c, input_bact_rarfy$map_loaded , by= "groundcontact_treatment")
# 
# colnames(rep_bact_adonis_2c_gc)
# 
# rep_bact_adonis_2c_gc_sc <- rep_bact_adonis_2c_gc[c("samples","groundcontact_treatment",
#                                                     "gc_connection","gc_suspended","groundcontact_groupcat")]
# 
# taxa_sum_otus_bact_adonis_gc_sc <- left_join(data_loadedbact_df, 
#                                              rep_bact_adonis_2c_gc_sc , by= "samples")
# str(taxa_sum_otus_bact_adonis_gc_sc)
# 
# adonis_otus_bact_gc_sus <- adonis(taxa_sum_otus_bact_adonis_gc_sc[,-c((ncol(taxa_sum_otus_bact_adonis_gc_sc)-4),
#                                                                       (ncol(taxa_sum_otus_bact_adonis_gc_sc)-3), 
#                                                                       (ncol(taxa_sum_otus_bact_adonis_gc_sc)-2),
#                                                                       (ncol(taxa_sum_otus_bact_adonis_gc_sc)-1),
#                                                                       ncol(taxa_sum_otus_bact_adonis_gc_sc))] ~ 
#                                     taxa_sum_otus_bact_adonis_gc_sc$gc_suspended,
#                                   permutations = 10000, method = "bray") 
# 
# 
# #RESULTS: SIGNIFICANT DIFERENCES ( 9.999e-05 *** < 0.05) !
# #Interpretation: Diversity patterns are being driven by elevation from the ground (up, down)
# 
# # Adonis for gc_connection:
# 
# adonis_otus_bact_gc_con <- adonis(taxa_sum_otus_bact_adonis_gc_sc[,-c((ncol(taxa_sum_otus_bact_adonis_gc_sc)-4),
#                                                                       (ncol(taxa_sum_otus_bact_adonis_gc_sc)-3), 
#                                                                       (ncol(taxa_sum_otus_bact_adonis_gc_sc)-2),
#                                                                       (ncol(taxa_sum_otus_bact_adonis_gc_sc)-1),
#                                                                       ncol(taxa_sum_otus_bact_adonis_gc_sc))] ~ 
#                                     taxa_sum_otus_bact_adonis_gc_sc$gc_connection,
#                                   permutations = 10000, method = "bray") 
# 
# 
# #RESULTS: SIGNIFICANT DIFERENCES ( 9.999e-05 *** < 0.05) !
# #Interpretation: Diversity patterns are being driven by connection with the ground (connected, separated)
# 
# # ~~~~~ Interaction (interaction(gc_suspended + gc_connection + gc_suspended:gc_connection))
# 
# adonis_otus_bact_gc_consus <- adonis(taxa_sum_otus_bact_adonis_gc_sc[,-c((ncol(taxa_sum_otus_bact_adonis_gc_sc)-4),
#                                                                          (ncol(taxa_sum_otus_bact_adonis_gc_sc)-3), 
#                                                                          (ncol(taxa_sum_otus_bact_adonis_gc_sc)-2),
#                                                                          (ncol(taxa_sum_otus_bact_adonis_gc_sc)-1),
#                                                                          ncol(taxa_sum_otus_bact_adonis_gc_sc))] ~ 
#                                        (taxa_sum_otus_bact_adonis_gc_sc$gc_connection +
#                                           taxa_sum_otus_bact_adonis_gc_sc$gc_suspended +
#                                           taxa_sum_otus_bact_adonis_gc_sc$gc_connection:taxa_sum_otus_bact_adonis_gc_sc$gc_suspended),
#                                      permutations = 10000, method = "bray") 
# 
# #RESULTS: SIGNIFICANT DIFERENCES FOR ALL ( Pr(>F)=9.999e-05 *** < 0.05) !
# #Interpretation: Diversity patterns are being driven by connection and elevation with the ground! Being elevated matters,
# # being connected matters and being connected changes the effect of being elevated. 
# 
# #See PERMANOVA results (C:/Users/angel/OneDrive/3_microbial_diversity_20220816/Primer_results)

#///////////////////////////////////////
### 4. MOST COMMMON ORDERS ####
#///////////////////////////////////////

setwd("C:/Users/angel/OneDrive/3_microbial_diversity_20220816/Exploration_results")

#import data
summarized_orders_bacteria <- read.csv("summarized_orders_bacteria.csv")
summarized_orders_fungi<- read.csv("summarized_orders_fungi.csv")
mapping_bactdiversity <- read.csv("mapping_bactdiversity.csv")
mapping_fungidiversity <- read.csv("mapping_fungidiversity.csv")

#remove the contaminated sample from bacteria
logical_drop_Csusp <- mapping_bactdiversity$groundcontact_groupcat=="c"&mapping_bactdiversity$gc_suspended=="up"&mapping_bactdiversity$gc_connection=="separated"
summarized_orders_bacteria <- summarized_orders_bacteria[!logical_drop_Csusp,]

#identify most common orders
#bact:
bact_comm_ord <- sort(colSums(summarized_orders_bacteria[,-c(1:3)]), decreasing=T) #-c(1:3) columns as characters, not orders
bact_comm_ord_df <- as.data.frame(bact_comm_ord)
head(rownames(bact_comm_ord_df), n=10) #10 orders
head(bact_comm_ord_df, n=10) #10 orders and abundances

#o__Rhizobiales              2.5973742
#o__Actinomycetales          1.3894967
#o__Rhodospirillales         1.3301032
#o__Xanthomonadales          1.2269459
#o__Myxococcales             1.0768990
#o__.Saprospirales.          0.9149734
#o__Gemmatales               0.9112223
#o__Burkholderiales          0.6398875
#o__Cytophagales             0.6123789
#o__.Chthoniobacterales.     0.5970616

#fungi:
fungi_comm_ord <- sort(colSums(summarized_orders_fungi[,-c(1)]), decreasing=T) #-c(1) column as character, not order
fungi_comm_ord_df <- as.data.frame(fungi_comm_ord)
head(rownames(fungi_comm_ord_df), n=10) #10 orders
head(fungi_comm_ord_df, n=10) #10 orders and abundances

#o__unidentified                           3.9719388
#o__Agaricales                             1.1653912
#o__Xylariales                             0.7166241
#o__Hypocreales                            0.6707058
#o__Chaetothyriales                        0.6611395
#o__Pezizomycotina_ord_Incertae_sedis      0.5801446
#o__Magnaporthales                         0.4936224
#o__Capnodiales                            0.4455782
#o__Pleosporales                           0.3843537
#o__Polyporales                            0.2869898

#add mapping data to order files 
#bact:
colnames(summarized_orders_bacteria)[1] <- "sample"
colnames(mapping_bactdiversity)[1]<- "sample"
summarized_orders_bacteria_trim <- summarized_orders_bacteria[,colnames(summarized_orders_bacteria) %in% c("sample",head(rownames(bact_comm_ord_df), n=10))]#select only the 10 most abundant orders
bact_ord_map <- left_join(summarized_orders_bacteria_trim, mapping_bactdiversity, by= "sample" ) 
str(bact_ord_map)

#fungi:
colnames(summarized_orders_fungi)[1] <- "sample"
colnames(mapping_fungidiversity)[1]<- "sample"
summarized_orders_fungi_trim <- summarized_orders_fungi[,colnames(summarized_orders_fungi) %in% c("sample",head(rownames(fungi_comm_ord_df), n=10))]#select only the 10 most abundant orders
fungi_ord_map <- left_join(summarized_orders_fungi_trim, mapping_fungidiversity, by= "sample" ) 
str(fungi_ord_map)



#///////////////////////////////////////
### 5. MODELS TO TEST ORDERS    ####
#///////////////////////////////////////

#bact: 
  #1.Rhizobiales
mod_rhiz1 <- lmer(o__Rhizobiales ~ gc_connection + gc_suspended + gc_connection:gc_suspended +
                         (1|groundcontact_groupcat),bact_ord_map)
mod_rhiz2 <- lmer(o__Rhizobiales ~ gc_connection + gc_suspended + 
                         (1|groundcontact_groupcat),bact_ord_map)
mod_rhiz3 <- lmer(o__Rhizobiales ~ gc_connection + 
                         (1|groundcontact_groupcat),bact_ord_map)
mod_rhiz4 <- lmer(o__Rhizobiales ~ gc_suspended + 
                         (1|groundcontact_groupcat),bact_ord_map)
mod_rhiz5 <- lmer(o__Rhizobiales ~ 1 + 
                         (1|groundcontact_groupcat),bact_ord_map)

anova(mod_rhiz1,mod_rhiz2,test="LRT")
anova(mod_rhiz2,mod_rhiz3,test="LRT") # 0.03535 *mod_rhiz3 
anova(mod_rhiz2,mod_rhiz4,test="LRT")
# anova(mod_rhiz4,mod_rhiz5,test="LRT") # 0.03617 * mod_rhiz4 

summary(mod_rhiz4)#run summary of full model to see effect direction
#Updated interpretation: Being suspended increases Rhizobiales relative abundance

#evaluate model fit of best-fit model
plot(mod_rhiz4)
qqnorm(resid(mod_rhiz4));qqline(resid(mod_rhiz4))
#these residuals look great!

#look at the data
ggplot(bact_ord_map,aes(x = gc_connection,y = o__Rhizobiales)) +
  geom_boxplot(aes(color=gc_suspended),shape=1)

  #2.Actinomycetales
mod_actin1 <- lmer(o__Actinomycetales ~ gc_connection + gc_suspended + gc_connection:gc_suspended +
                    (1|groundcontact_groupcat),bact_ord_map)
mod_actin2 <- lmer(o__Actinomycetales ~ gc_connection + gc_suspended + 
                    (1|groundcontact_groupcat),bact_ord_map)
mod_actin3 <- lmer(o__Actinomycetales ~ gc_connection + 
                    (1|groundcontact_groupcat),bact_ord_map)
mod_actin4 <- lmer(o__Actinomycetales ~ gc_suspended + 
                    (1|groundcontact_groupcat),bact_ord_map)
mod_actin5 <- lmer(o__Actinomycetales ~ 1 + 
                    (1|groundcontact_groupcat),bact_ord_map)

anova(mod_actin1,mod_actin2,test="LRT")
anova(mod_actin2,mod_actin3,test="LRT") # 0.03697  
anova(mod_actin2,mod_actin4,test="LRT")
anova(mod_actin4,mod_actin5,test="LRT") # 0.03737 
#Interpretation: Being suspended affects the relative abundance of Actinomycetales (+suspended)

summary(mod_actin4)#run summary of full model to see effect direction
#Updated interpretation: Being suspended increases Actinomycetales relative abundance

#evaluate model fit of best-fit model
plot(mod_actin4)
qqnorm(resid(mod_actin4));qqline(resid(mod_actin4))
#these residuals look decent

#look at the data
ggplot(bact_ord_map,aes(x = gc_connection,y = o__Actinomycetales)) +
  geom_boxplot(aes(color=gc_suspended),shape=1)

#3.Rhodospirillales 
mod_rhod1 <- lmer(o__Rhodospirillales ~ gc_connection + gc_suspended + gc_connection:gc_suspended +
                     (1|groundcontact_groupcat),bact_ord_map)
mod_rhod2 <- lmer(o__Rhodospirillales ~ gc_connection + gc_suspended + 
                     (1|groundcontact_groupcat),bact_ord_map)
mod_rhod3 <- lmer(o__Rhodospirillales ~ gc_connection + 
                     (1|groundcontact_groupcat),bact_ord_map)
mod_rhod4 <- lmer(o__Rhodospirillales ~ gc_suspended + 
                     (1|groundcontact_groupcat),bact_ord_map)
mod_rhod5 <- lmer(o__Rhodospirillales ~ 1 + 
                     (1|groundcontact_groupcat),bact_ord_map)

anova(mod_rhod1,mod_rhod2,test="LRT")
anova(mod_rhod2,mod_rhod3,test="LRT") 
anova(mod_rhod2,mod_rhod4,test="LRT")
anova(mod_rhod4,mod_rhod5,test="LRT") 
#Interpretation: No effect

#confirm appropriate residuals with full model
plot(mod_rhod1)
qqnorm(resid(mod_rhod1));qqline(resid(mod_rhod1))
#these residuals look decent

#look at the data
ggplot(bact_ord_map,aes(x = gc_connection,y = o__Rhodospirillales)) +
  geom_boxplot(aes(color=gc_suspended),shape=1)

#4.Xanthomonadales 

mod_xan1 <- lmer(o__Xanthomonadales  ~ gc_connection + gc_suspended + gc_connection:gc_suspended +
                    (1|groundcontact_groupcat),bact_ord_map)
mod_xan2 <- lmer(o__Xanthomonadales  ~ gc_connection + gc_suspended + 
                    (1|groundcontact_groupcat),bact_ord_map)
mod_xan3 <- lmer(o__Xanthomonadales  ~ gc_connection + 
                    (1|groundcontact_groupcat),bact_ord_map)
mod_xan4 <- lmer(o__Xanthomonadales  ~ gc_suspended + 
                    (1|groundcontact_groupcat),bact_ord_map)
mod_xan5 <- lmer(o__Xanthomonadales  ~ 1 + 
                    (1|groundcontact_groupcat),bact_ord_map)

anova(mod_xan1,mod_xan2,test="LRT")
anova(mod_xan2,mod_xan3,test="LRT") 
anova(mod_xan2,mod_xan4,test="LRT")
anova(mod_xan4,mod_xan5,test="LRT") 
#Interpretation: No effect

#confirm appropriate residuals with full model
plot(mod_xan1)
qqnorm(resid(mod_xan1));qqline(resid(mod_xan1))
#these residuals look decent

#look at the data
ggplot(bact_ord_map,aes(x = gc_connection,y = o__Xanthomonadales)) +
  geom_boxplot(aes(color=gc_suspended),shape=1)

#5.Myxococcales

mod_myx1 <- lmer(o__Myxococcales  ~ gc_connection + gc_suspended + gc_connection:gc_suspended +
                   (1|groundcontact_groupcat),bact_ord_map)
mod_myx2 <- lmer(o__Myxococcales  ~ gc_connection + gc_suspended + 
                   (1|groundcontact_groupcat),bact_ord_map)
mod_myx3 <- lmer(o__Myxococcales  ~ gc_connection + 
                   (1|groundcontact_groupcat),bact_ord_map)
mod_myx4 <- lmer(o__Myxococcales  ~ gc_suspended + 
                   (1|groundcontact_groupcat),bact_ord_map)
mod_myx5 <- lmer(o__Myxococcales  ~ 1 + 
                   (1|groundcontact_groupcat),bact_ord_map)

anova(mod_myx1,mod_myx2,test="LRT")
anova(mod_myx2,mod_myx3,test="LRT") 
anova(mod_myx2,mod_myx4,test="LRT")
anova(mod_myx4,mod_myx5,test="LRT") 
#Interpretation: No effect

#confirm appropriate residuals with full model
plot(mod_myx1)
qqnorm(resid(mod_myx1));qqline(resid(mod_myx1))
#these residuals look ok - not great though
#### log transforming abundance slightly improves model fit, but result is unchanged

#look at the data
ggplot(bact_ord_map,aes(x = gc_connection,y = o__Myxococcales)) +
  geom_boxplot(aes(color=gc_suspended),shape=1)

#6.Saprospirales

mod_sapr1 <- lmer(o__.Saprospirales.  ~ gc_connection + gc_suspended + gc_connection:gc_suspended +
                   (1|groundcontact_groupcat),bact_ord_map)
mod_sapr2 <- lmer(o__.Saprospirales.  ~ gc_connection + gc_suspended + 
                   (1|groundcontact_groupcat),bact_ord_map)
mod_sapr3 <- lmer(o__.Saprospirales.  ~ gc_connection + 
                   (1|groundcontact_groupcat),bact_ord_map)
mod_sapr4 <- lmer(o__.Saprospirales.  ~ gc_suspended + 
                   (1|groundcontact_groupcat),bact_ord_map)
mod_sapr5 <- lmer(o__.Saprospirales.  ~ 1 + 
                   (1|groundcontact_groupcat),bact_ord_map)

anova(mod_sapr1,mod_sapr2,test="LRT")
anova(mod_sapr2,mod_sapr3,test="LRT") # p value 0.01591 
anova(mod_sapr2,mod_sapr4,test="LRT") 
anova(mod_sapr4,mod_sapr5,test="LRT") # p value 0.01597 
#Interpretation: Being suspended affects the relative abundance of Saprospirales (+suspended)

summary(mod_sapr4)#run summary of full model to see effect direction
#Updated interpretation: Being suspended increases Saprospirales relative abundance

#evaluate model fit of best-fit model
plot(mod_sapr4)
qqnorm(resid(mod_sapr4));qqline(resid(mod_sapr4))
#these residuals look ok, but not great - poor fit at the extremes
#slight improvement in fit with log-transformation, but the result is the same

#look at the data
ggplot(bact_ord_map,aes(x = gc_connection,y = o__.Saprospirales.)) +
  geom_boxplot(aes(color=gc_suspended),shape=1)

#7. o__Gemmatales

mod_gem1 <- lmer(o__Gemmatales  ~ gc_connection + gc_suspended + gc_connection:gc_suspended +
                   (1|groundcontact_groupcat),bact_ord_map)
mod_gem2 <- lmer(o__Gemmatales  ~ gc_connection + gc_suspended + 
                   (1|groundcontact_groupcat),bact_ord_map)
mod_gem3 <- lmer(o__Gemmatales  ~ gc_connection + 
                   (1|groundcontact_groupcat),bact_ord_map)
mod_gem4 <- lmer(o__Gemmatales  ~ gc_suspended + 
                   (1|groundcontact_groupcat),bact_ord_map)
mod_gem5 <- lmer(o__Gemmatales  ~ 1 + 
                   (1|groundcontact_groupcat),bact_ord_map)

anova(mod_gem1,mod_gem2,test="LRT")
anova(mod_gem2,mod_gem3,test="LRT") 
anova(mod_gem2,mod_gem4,test="LRT")
anova(mod_gem4,mod_gem5,test="LRT") 
#Interpretation: No effect

#confirm appropriate residuals with full model
plot(mod_gem1)
qqnorm(resid(mod_gem1));qqline(resid(mod_gem1))
#these residuals look great

#look at the data
ggplot(bact_ord_map,aes(x = gc_connection,y = o__Gemmatales)) +
  geom_boxplot(aes(color=gc_suspended),shape=1)

#8. o__Burkholderiales

mod_burk1 <- lmer(o__Burkholderiales  ~ gc_connection + gc_suspended + gc_connection:gc_suspended +
                   (1|groundcontact_groupcat),bact_ord_map)
mod_burk2 <- lmer(o__Burkholderiales  ~ gc_connection + gc_suspended + 
                   (1|groundcontact_groupcat),bact_ord_map)
mod_burk3 <- lmer(o__Burkholderiales  ~ gc_connection + 
                   (1|groundcontact_groupcat),bact_ord_map)
mod_burk4 <- lmer(o__Burkholderiales  ~ gc_suspended + 
                   (1|groundcontact_groupcat),bact_ord_map)
mod_burk5 <- lmer(o__Burkholderiales  ~ 1 + 
                   (1|groundcontact_groupcat),bact_ord_map)

anova(mod_burk1,mod_burk2,test="LRT")
anova(mod_burk2,mod_burk3,test="LRT") 
anova(mod_burk2,mod_burk4,test="LRT")
anova(mod_burk4,mod_burk5,test="LRT") 
#Interpretation: No effect

#confirm appropriate residuals with full model
plot(mod_burk1)
qqnorm(resid(mod_burk1));qqline(resid(mod_burk1))
#these residuals look ok, but not great
#no meaningful change with log transformation

#look at the data
ggplot(bact_ord_map,aes(x = gc_connection,y = o__Burkholderiales)) +
  geom_boxplot(aes(color=gc_suspended),shape=1)

#9. o__Cytophagales             

mod_cyto1 <- lmer(o__Cytophagales  ~ gc_connection + gc_suspended + gc_connection:gc_suspended +
                   (1|groundcontact_groupcat),bact_ord_map)
mod_cyto2 <- lmer(o__Cytophagales  ~ gc_connection + gc_suspended + 
                   (1|groundcontact_groupcat),bact_ord_map)
mod_cyto3 <- lmer(o__Cytophagales  ~ gc_connection + 
                   (1|groundcontact_groupcat),bact_ord_map)
mod_cyto4 <- lmer(o__Cytophagales  ~ gc_suspended + 
                   (1|groundcontact_groupcat),bact_ord_map)
mod_cyto5 <- lmer(o__Cytophagales  ~ 1 + 
                   (1|groundcontact_groupcat),bact_ord_map)

anova(mod_cyto1,mod_cyto2,test="LRT")
anova(mod_cyto2,mod_cyto3,test="LRT") 
anova(mod_cyto2,mod_cyto4,test="LRT")
anova(mod_cyto4,mod_cyto5,test="LRT") 
#Interpretation: No effect

#confirm appropriate residuals with full model
plot(mod_cyto1)
qqnorm(resid(mod_cyto1));qqline(resid(mod_cyto1))
#great residuals

#look at the data
ggplot(bact_ord_map,aes(x = gc_connection,y = o__Cytophagales)) +
  geom_boxplot(aes(color=gc_suspended),shape=1)

#10. o__.Chthoniobacterales.

mod_chth1 <- lmer(o__.Chthoniobacterales.  ~ gc_connection + gc_suspended + gc_connection:gc_suspended +
                   (1|groundcontact_groupcat),bact_ord_map)
mod_chth2 <- lmer(o__.Chthoniobacterales.  ~ gc_connection + gc_suspended + 
                   (1|groundcontact_groupcat),bact_ord_map)
mod_chth3 <- lmer(o__.Chthoniobacterales.  ~ gc_connection + 
                   (1|groundcontact_groupcat),bact_ord_map)
mod_chth4 <- lmer(o__.Chthoniobacterales.  ~ gc_suspended + 
                   (1|groundcontact_groupcat),bact_ord_map)
mod_chth5 <- lmer(o__.Chthoniobacterales.  ~ 1 + 
                   (1|groundcontact_groupcat),bact_ord_map)

anova(mod_chth1,mod_chth2,test="LRT")
anova(mod_chth2,mod_chth3,test="LRT") 
anova(mod_chth2,mod_chth4,test="LRT")
anova(mod_chth4,mod_chth5,test="LRT") 
#Interpretation: No effect

#confirm appropriate residuals with full model
plot(mod_chth1)
qqnorm(resid(mod_chth1));qqline(resid(mod_chth1))
#good residuals

#look at the data
ggplot(bact_ord_map,aes(x = gc_connection,y = o__.Chthoniobacterales.)) +
  geom_boxplot(aes(color=gc_suspended),shape=1)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#fungi:
#2. Agaricales (1st order was unidentified)  

mod_aga1 <- lmer(o__Agaricales ~ gc_connection + gc_suspended + gc_connection:gc_suspended +
                         (1|groundcontact_groupcat),fungi_ord_map)
mod_aga2 <- lmer(o__Agaricales ~ gc_connection + gc_suspended + 
                         (1|groundcontact_groupcat),fungi_ord_map)
mod_aga3 <- lmer(o__Agaricales ~ gc_connection + 
                         (1|groundcontact_groupcat),fungi_ord_map)
mod_aga4 <- lmer(o__Agaricales ~ gc_suspended + 
                         (1|groundcontact_groupcat),fungi_ord_map)
mod_aga5 <- lmer(o__Agaricales ~ 1 + 
                         (1|groundcontact_groupcat),fungi_ord_map)

anova(mod_aga1,mod_aga2,test="LRT")
anova(mod_aga2,mod_aga3,test="LRT") 
anova(mod_aga2,mod_aga4,test="LRT")
anova(mod_aga4,mod_aga5,test="LRT") 
#Interpretation: No effect (marginal interaction)

#confirm appropriate residuals with full model
plot(mod_aga1)
qqnorm(resid(mod_aga1));qqline(resid(mod_aga1))
#ok residuals (very few points)

#look at the data
ggplot(fungi_ord_map,aes(x = gc_connection,y = o__Agaricales)) +
  geom_point(aes(color=gc_suspended),shape=1)

#3. Xylariales  

mod_xyl1 <- lmer(o__Xylariales ~ gc_connection + gc_suspended + gc_connection:gc_suspended +
                   (1|groundcontact_groupcat),fungi_ord_map)
mod_xyl2 <- lmer(o__Xylariales ~ gc_connection + gc_suspended + 
                   (1|groundcontact_groupcat),fungi_ord_map)
mod_xyl3 <- lmer(o__Xylariales ~ gc_connection + 
                   (1|groundcontact_groupcat),fungi_ord_map) # boundary (singular) fit: see help('isSingular') ?
mod_xyl4 <- lmer(o__Xylariales ~ gc_suspended + 
                   (1|groundcontact_groupcat),fungi_ord_map) # boundary (singular) fit: see help('isSingular') ?
mod_xyl5 <- lmer(o__Xylariales ~ 1 + 
                   (1|groundcontact_groupcat),fungi_ord_map) # boundary (singular) fit: see help('isSingular') ?

anova(mod_xyl1,mod_xyl2,test="LRT")
anova(mod_xyl2,mod_xyl3,test="LRT") # p value 0.0374 * (mod_xyl2)
anova(mod_xyl2,mod_xyl4,test="LRT")
# anova(mod_xyl4,mod_xyl5,test="LRT") 
#Interpretation: Being suspended affects the relative abundance of Xylariales (+downed)

summary(mod_xyl4)#run summary of full model to see effect direction

#confirm appropriate residuals with full model
plot(mod_xyl1)
qqnorm(resid(mod_xyl1));qqline(resid(mod_xyl1))
#ok residuals (very few points)

#look at the data
ggplot(fungi_ord_map,aes(x = gc_connection,y = o__Xylariales)) +
  geom_point(aes(color=gc_suspended),shape=1)

#4. Hypocreales  

mod_hyp1 <- lmer(o__Hypocreales ~ gc_connection + gc_suspended + gc_connection:gc_suspended +
                   (1|groundcontact_groupcat),fungi_ord_map)
mod_hyp2 <- lmer(o__Hypocreales ~ gc_connection + gc_suspended + 
                   (1|groundcontact_groupcat),fungi_ord_map)
mod_hyp3 <- lmer(o__Hypocreales ~ gc_connection + 
                   (1|groundcontact_groupcat),fungi_ord_map)
mod_hyp4 <- lmer(o__Hypocreales ~ gc_suspended + 
                   (1|groundcontact_groupcat),fungi_ord_map)
mod_hyp5 <- lmer(o__Hypocreales ~ 1 + 
                   (1|groundcontact_groupcat),fungi_ord_map)

anova(mod_hyp1,mod_hyp2,test="LRT") # p value 0.03117 * (mod_hyp1)
### Stop because the interaction is significant
# anova(mod_hyp2,mod_hyp3,test="LRT") 
# anova(mod_hyp2,mod_hyp4,test="LRT")
# anova(mod_hyp4,mod_hyp5,test="LRT") 
#Interpretation: The interaction of being suspended and connected 
#affect the relative abundance of Hypocreales

summary(mod_hyp1)#run summary of full model to see effect direction

#confirm appropriate residuals with full model
plot(mod_hyp1)
qqnorm(resid(mod_hyp1));qqline(resid(mod_hyp1))
#ok residuals (very few points)

#look at the data
ggplot(fungi_ord_map,aes(x = gc_connection,y = o__Hypocreales)) +
  geom_point(aes(color=gc_suspended),shape=1)

#5. o__Chaetothyriales 

mod_chae1 <- lmer(o__Chaetothyriales ~ gc_connection + gc_suspended + gc_connection:gc_suspended +
                   (1|groundcontact_groupcat),fungi_ord_map)
mod_chae2 <- lmer(o__Chaetothyriales ~ gc_connection + gc_suspended + 
                   (1|groundcontact_groupcat),fungi_ord_map)
mod_chae3 <- lmer(o__Chaetothyriales ~ gc_connection + 
                   (1|groundcontact_groupcat),fungi_ord_map)
mod_chae4 <- lmer(o__Chaetothyriales ~ gc_suspended + 
                   (1|groundcontact_groupcat),fungi_ord_map)
mod_chae5 <- lmer(o__Chaetothyriales ~ 1 + 
                   (1|groundcontact_groupcat),fungi_ord_map)

anova(mod_chae1,mod_chae2,test="LRT") 
anova(mod_chae2,mod_chae3,test="LRT") 
anova(mod_chae2,mod_chae4,test="LRT") # p value 0.003648 ** (mod_chae2)
# anova(mod_chae3,mod_chae5,test="LRT")  
#Interpretation: Being connected affect the relative abundance of Chaetothyriales

summary(mod_chae3)

#confirm appropriate residuals with full model
plot(mod_chae3)
qqnorm(resid(mod_chae3));qqline(resid(mod_chae3))
#ok residuals (very few points)

#look at the data
ggplot(fungi_ord_map,aes(x = gc_connection,y = o__Chaetothyriales)) +
  geom_point(aes(color=gc_suspended),shape=1)

#6. o__Pezizomycotina_ord_Incertae_sedis

mod_pez1 <- lmer(o__Pezizomycotina_ord_Incertae_sedis ~ gc_connection + gc_suspended + gc_connection:gc_suspended +
                    (1|groundcontact_groupcat),fungi_ord_map)
mod_pez2 <- lmer(o__Pezizomycotina_ord_Incertae_sedis ~ gc_connection + gc_suspended + 
                    (1|groundcontact_groupcat),fungi_ord_map)
mod_pez3 <- lmer(o__Pezizomycotina_ord_Incertae_sedis ~ gc_connection + 
                    (1|groundcontact_groupcat),fungi_ord_map)
mod_pez4 <- lmer(o__Pezizomycotina_ord_Incertae_sedis ~ gc_suspended + 
                    (1|groundcontact_groupcat),fungi_ord_map)
mod_pez5 <- lmer(o__Pezizomycotina_ord_Incertae_sedis ~ 1 + 
                    (1|groundcontact_groupcat),fungi_ord_map)

anova(mod_pez1,mod_pez2,test="LRT") 
anova(mod_pez2,mod_pez3,test="LRT") #marginally significant
anova(mod_pez2,mod_pez4,test="LRT") 
anova(mod_pez4,mod_pez5,test="LRT")  
#Interpretation: No effect (marginally significant effect of suspension)

#confirm appropriate residuals with full model
plot(mod_pez4)
qqnorm(resid(mod_pez4));qqline(resid(mod_pez4))
#ok residuals (very few points)

#look at the data
ggplot(fungi_ord_map,aes(x = gc_suspended,y = o__Pezizomycotina_ord_Incertae_sedis)) +
  geom_point(aes(color=gc_connection),shape=1)

#7.o__Magnaporthales                         

mod_magn1 <- lmer(o__Magnaporthales   ~ gc_connection + gc_suspended + gc_connection:gc_suspended +
                   (1|groundcontact_groupcat),fungi_ord_map)
mod_magn2 <- lmer(o__Magnaporthales   ~ gc_connection + gc_suspended + 
                   (1|groundcontact_groupcat),fungi_ord_map)
mod_magn3 <- lmer(o__Magnaporthales   ~ gc_connection + 
                   (1|groundcontact_groupcat),fungi_ord_map)
mod_magn4 <- lmer(o__Magnaporthales   ~ gc_suspended + 
                   (1|groundcontact_groupcat),fungi_ord_map)
mod_magn5 <- lmer(o__Magnaporthales   ~ 1 + 
                   (1|groundcontact_groupcat),fungi_ord_map)

anova(mod_magn1,mod_magn2,test="LRT") 
anova(mod_magn2,mod_magn3,test="LRT") # p value 0.04146 * (mod_magn2)
anova(mod_magn2,mod_magn4,test="LRT") 
anova(mod_magn4,mod_magn5,test="LRT") # p value 0.04261 * (mod_magn4)
#Interpretation: Being suspended affect the relative abundance of Magnaporthales (+down)

summary(mod_magn4)#run summary of full model to see effect direction

#confirm appropriate residuals with full model
plot(mod_magn4)
qqnorm(resid(mod_magn4));qqline(resid(mod_magn4))
#ok residuals (very few points)

#look at the data
ggplot(fungi_ord_map,aes(x = gc_suspended,y = o__Magnaporthales )) +
  geom_point(aes(color=gc_connection),shape=1) 


#8. o__Capnodiales  

mod_capn1 <- lmer(o__Capnodiales  ~ gc_connection + gc_suspended + gc_connection:gc_suspended +
                   (1|groundcontact_groupcat),fungi_ord_map)
mod_capn2 <- lmer(o__Capnodiales  ~ gc_connection + gc_suspended + 
                   (1|groundcontact_groupcat),fungi_ord_map)
mod_capn3 <- lmer(o__Capnodiales  ~ gc_connection + 
                   (1|groundcontact_groupcat),fungi_ord_map)
mod_capn4 <- lmer(o__Capnodiales  ~ gc_suspended + 
                   (1|groundcontact_groupcat),fungi_ord_map)
mod_capn5 <- lmer(o__Capnodiales  ~ 1 + 
                   (1|groundcontact_groupcat),fungi_ord_map)

anova(mod_capn1,mod_capn2,test="LRT") # p value  0.05071 (mod_capn1)
anova(mod_capn2,mod_capn3,test="LRT") 
anova(mod_capn2,mod_capn4,test="LRT") 
anova(mod_capn4,mod_capn5,test="LRT")  
#Interpretation: The interaction of being suspended and connected 
#affect the relative abundance of Capnodiales

summary(mod_capn1)#run summary of full model to see effect direction

#confirm appropriate residuals with full model
plot(mod_capn1)
qqnorm(resid(mod_capn1));qqline(resid(mod_capn1))
#ok residuals (very few points)

#look at the data
ggplot(fungi_ord_map,aes(x = gc_suspended,y = o__Capnodiales )) +
  geom_point(aes(color=gc_connection),shape=1) 

#9. o__Pleosporales   

mod_ple1 <- lmer(o__Pleosporales ~ gc_connection + gc_suspended + gc_connection:gc_suspended +
                   (1|groundcontact_groupcat),fungi_ord_map)
mod_ple2 <- lmer(o__Pleosporales ~ gc_connection + gc_suspended + 
                   (1|groundcontact_groupcat),fungi_ord_map)
mod_ple3 <- lmer(o__Pleosporales ~ gc_connection + 
                   (1|groundcontact_groupcat),fungi_ord_map)
mod_ple4 <- lmer(o__Pleosporales ~ gc_suspended + 
                   (1|groundcontact_groupcat),fungi_ord_map)
mod_ple5 <- lmer(o__Pleosporales ~ 1 + 
                   (1|groundcontact_groupcat),fungi_ord_map)

anova(mod_ple1,mod_ple2,test="LRT") 
anova(mod_ple2,mod_ple3,test="LRT") 
anova(mod_ple2,mod_ple4,test="LRT") # p value 0.03367 * (mod_ple2)
anova(mod_ple4,mod_ple5,test="LRT")  
#Interpretation: Being connected affect the relative abundance of Pleosporales (+suspended)

summary(mod_ple3)#run summary of full model to see effect direction

#confirm appropriate residuals with full model
plot(mod_ple4)
qqnorm(resid(mod_ple4));qqline(resid(mod_ple4))
#ok residuals (very few points)

#look at the data
ggplot(fungi_ord_map,aes(x = gc_suspended,y = o__Pleosporales )) +
  geom_point(aes(color=gc_connection),shape=1) 

#10. o__Polyporales

mod_poly1 <- lmer(o__Polyporales ~ gc_connection + gc_suspended + gc_connection:gc_suspended +
                   (1|groundcontact_groupcat),fungi_ord_map)
mod_poly2 <- lmer(o__Polyporales ~ gc_connection + gc_suspended + 
                   (1|groundcontact_groupcat),fungi_ord_map)
mod_poly3 <- lmer(o__Polyporales ~ gc_connection + 
                   (1|groundcontact_groupcat),fungi_ord_map)
mod_poly4 <- lmer(o__Polyporales ~ gc_suspended + 
                   (1|groundcontact_groupcat),fungi_ord_map)
mod_poly5 <- lmer(o__Polyporales ~ 1 + 
                   (1|groundcontact_groupcat),fungi_ord_map)

anova(mod_poly1,mod_poly2,test="LRT") 
anova(mod_poly2,mod_poly3,test="LRT") # p value 0.04792 * (mod_poly2)
anova(mod_poly2,mod_poly4,test="LRT") 
anova(mod_poly4,mod_poly5,test="LRT") # p value 0.04763 (mod_poly4)
#Interpretation: Being suspended affect the relative abundance of Polyporales (+suspended)

summary(mod_poly4)#run summary of full model to see effect direction

#confirm appropriate residuals with full model
plot(mod_poly4)
qqnorm(resid(mod_poly4));qqline(resid(mod_poly4))
#ok residuals (very few points)

#look at the data
ggplot(fungi_ord_map,aes(x = gc_suspended,y = o__Polyporales )) +
  geom_point(aes(color=gc_connection),shape=1) 

#/////////////////////////
#### Table for presenting order-level results ####
#/////////////////////////
# create tables for presenting data and results
# bacteria 
figure_divesity_bact <- bact_ord_map %>% group_by(gc_connection,gc_suspended) %>%
  summarise(Chthoniobacterales_mean = mean(o__.Chthoniobacterales.),
            Chthoniobacterales_sd = sd(o__.Chthoniobacterales.),
            Saprospirales_mean = mean(o__.Saprospirales.),
            Saprospirales_sd = sd(o__.Saprospirales.),
            Actinomycetales_mean = mean(o__Actinomycetales),
            Actinomycetales_sd = sd(o__Actinomycetales),
            Burkholderiales_mean = mean(o__Burkholderiales),
            Burkholderiales_sd = sd(o__Burkholderiales),
            Cytophagales_mean = mean(o__Cytophagales),
            Cytophagales_sd = sd(o__Cytophagales),
            Gemmatales_mean = mean(o__Gemmatales),
            Gemmatales_sd = sd(o__Gemmatales),
            Myxococcales_mean = mean(o__Myxococcales),
            Myxococcales_sd = sd(o__Myxococcales),
            Rhizobiales_mean = mean(o__Rhizobiales),
            Rhizobiales_sd = sd(o__Rhizobiales),
            Rhodospirillales_mean = mean(o__Rhodospirillales),
            Rhodospirillales_sd = sd(o__Rhodospirillales),
            Xanthomonadales_mean = mean(o__Xanthomonadales),
            Xanthomonadales_sd = sd(o__Xanthomonadales),
            N = length(o__.Chthoniobacterales.)) %>%
  mutate(Chthoniobacterales_se = Chthoniobacterales_sd/sqrt(N),
         Saprospirales_se = Saprospirales_sd/sqrt(N),
         Actinomycetales_se = Actinomycetales_sd/sqrt(N),
         Burkholderiales_se = Burkholderiales_sd/sqrt(N),
         Cytophagales_se = Cytophagales_sd/sqrt(N),
         Gemmatales_se = Gemmatales_sd/sqrt(N),
         Myxococcales_se =  Myxococcales_sd/sqrt(N),
         Rhizobiales_se = Rhizobiales_sd/sqrt(N),
         Rhodospirillales_se =  Rhodospirillales_sd/sqrt(N),
         Xanthomonadales_se =  Xanthomonadales_sd/sqrt(N),)

# fungi 

figure_divesity_fungi <- fungi_ord_map %>% group_by(gc_connection,gc_suspended) %>%
  summarise(Agaricales_mean = mean(o__Agaricales),
            Agaricales_sd = sd(o__Agaricales),
            Xylariales_mean = mean(o__Xylariales),
            Xylariales_sd = sd(o__Xylariales),
            Hypocreales_mean = mean(o__Hypocreales),
            Hypocreales_sd = sd(o__Hypocreales),
            Chaetothyriales_mean = mean(o__Chaetothyriales),
            Chaetothyriales_sd = sd(o__Chaetothyriales),
            Pezizomycotina_incertae_sedis_mean = mean(o__Pezizomycotina_ord_Incertae_sedis),
            Pezizomycotina_incertae_sedis_sd = sd(o__Pezizomycotina_ord_Incertae_sedis),
            Magnaporthales_mean = mean(o__Magnaporthales),
            Magnaporthales_sd = sd(o__Magnaporthales),
            Capnodiales_mean = mean(o__Capnodiales),
            Capnodiales_sd = sd(o__Capnodiales),
            Pleosporales_mean = mean(o__Pleosporales),
            Pleosporales_sd = sd(o__Pleosporales),
            Polyporales_mean = mean(o__Polyporales),
            Polyporales_sd = sd(o__Polyporales),
            N = length(o__Polyporales)) %>%
  mutate(Agaricales_se = Agaricales_sd/sqrt(N),
         Xylariales_se = Xylariales_sd/sqrt(N),
         Hypocreales_se = Hypocreales_sd/sqrt(N),
         Chaetothyriales_se = Chaetothyriales_sd/sqrt(N),
         Pezizomycotina_incertae_sedis_se = Pezizomycotina_incertae_sedis_sd/sqrt(N),
         Magnaporthales_se = Magnaporthales_sd/sqrt(N),
         Capnodiales_se =  Capnodiales_sd/sqrt(N),
         Pleosporales_se = Pleosporales_sd/sqrt(N),
         Polyporales_se =  Polyporales_sd/sqrt(N))

#Save the figure as csv to modify them in excel 
setwd('C:/Users/angel/OneDrive/3_microbial_diversity_20220816/figures')
write.csv(figure_divesity_bact,"figure_divesity_bact.csv")
write.csv(figure_divesity_fungi,"figure_divesity_fungi.csv")
getwd()

#///////////////////////////////
#### Create NMDS figures: ####
#///////////////////////////////
library(vegan)

#check data
str(rarefied3199_bacteria)
str(mapping_bactdiversity)
str(rarefied4704_fungi)
str(mapping_fungidiversity)

#remove the taxonomy data
ncol(rarefied4704_fungi)
rarfy_fung_trim <- rarefied4704_fungi[,1:11]
str(rarfy_fung_trim)
ncol(rarefied3199_bacteria)
rarfy_bact_trim <- rarefied3199_bacteria[,1:20]
str(rarfy_bact_trim)


## Bacteria NMDS ####

#create square-root transformed data
bact.sqrt. <- sqrt(rarfy_bact_trim) # take the square root of the whole bacterial community dataframe. 
#This decreases the influence of super abundant taxa

taxa.matrix.b. <- t(bact.sqrt.) #now transpose the data so that samples are rows not columns 
matrix.df.b. <- as.data.frame(taxa.matrix.b.) # create a data frame

# #BACTERIA: one sample is driving odd patterns - what if we dropped it?
# #group c, suspended is the block we need to drop to visualize the data (mapping_bactdiversity we drop it from all analyses?)
logical_drop_Csusp <- mapping_bactdiversity$groundcontact_groupcat=="c"&mapping_bactdiversity$gc_suspended=="up"&mapping_bactdiversity$gc_connection=="separated"

#drop this sample
matrix.df.b_trim <- matrix.df.b.[!logical_drop_Csusp,]
# bact_ord_map_trim <- bact_ord_map[!logical_drop_Csusp,]

# Non-metric multidimensional scaling analysis = NMS
bact_mds <-metaMDS(matrix.df.b_trim,distance = "bray", k=2,try=1000,autotransform = FALSE,maxit=100)
distance <-vegdist(matrix.df.b_trim,method="bray")

#create vectors with the NMS attributes
NMS_coordinates<-scores(bact_mds,display="sites")
Bact_NMS_axes<-as.data.frame(NMS_coordinates)
#write.csv(Bact_NMS_axes,"Bact_NMS_axes.csv")
NMS_OTUscores<-scores(bact_mds,display="species")

#add NMS coordinates to the mapping file
#create dataframe with NMS Coordinates and Mapping File information
#add the proportion of total sequences of each vector to the "for_plotting" object below
### use only the 10 most common orders
for_ploting<-as.data.frame(cbind(NMS_coordinates,bact_ord_map))
str(for_ploting) #20 obs. of 27 variables

#now create the vectors
fit1<-envfit(bact_mds, bact_ord_map[,2:11], perm = 999, na.rm = TRUE)
fit1
str(for_ploting)

#create the figure -> NMDS
tiff(file="Bact_NMS_experiment.tiff", width = 1200*3.25, height = 1200*2, res = 1200)
par(mar=c(2,2,1,1))
plot(for_ploting$NMDS1 ~ for_ploting$NMDS2,
     xlab = "Bacterial NMS1",
     ylab = "Bacterial NMS2",
     font=2,
     font.lab=2,
     cex.axis=.7,
     cex.lab=.7,
     xlim = c(-.55,.35),
     mgp=c(1.15,.4,0),
     col = c("blue3","light blue","red4","lightcoral")[as.factor(for_ploting$groundcontact_treatment)],
     # bg = c("blue3","light blue","red4","lightcoral")[as.factor(unique(for_ploting$groundcontact_treatment))],
     #set the point types (plotting 'character'=pch) and colors based on your data (4 groups)
     pch = c(0,15,1,16)[as.factor(for_ploting$groundcontact_treatment)], cex=.8, 
     data = for_ploting)

#HULLS
ordiellipse(bact_mds, group=for_ploting$groundcontact_treatment,kind = "se", 
            conf=0.95, lwd=1.9, 
            #set the linetypes and points based on your data (not the same as here)
            lty = c("solid","solid","solid","solid"), 
            col = c("blue3","light blue","red4","lightcoral")[as.factor(unique(for_ploting$groundcontact_treatment))])

# names_inter<-gsub("o__","",names(fit1$vectors$arrows[,1]))
# bact_order_labels<-gsub("o__.","",names_inter)
# #Plot labels 
# plot(fit1,font=1, p.max = 0.05, col = "black", cex=0.6,
#      labels = bact_order_labels)
legend(
  x ="bottomleft",
  legend = c("down-connected","down-separated","susp.-separated","susp.-connected"), # for readability of legend
  pch = c(0,15,1,16),
  col =c("blue3","light blue","red4","lightcoral"),
  cex = 0.49 # scale the legend to look attractively sized
)

dev.off()
getwd()


##   Fungi NMDS  ####

#create square-root transformed data
fung.sqrt. <- sqrt(rarfy_fung_trim) # take the square root of the whole fungal community dataframe. 
#This decreases the influence of super abundant taxa

taxa.matrix.b.fung <- t(fung.sqrt.) #now transpose the data so that samples are rows not columns 
matrix.df.b.fung <- as.data.frame(taxa.matrix.b.fung) # create a data frame

# Non-metric multidimensional scaling analysis = NMS
fung_mds <-metaMDS(matrix.df.b.fung,distance = "bray", k=2,try=1000,autotransform = FALSE,maxit=100)
distance <-vegdist(matrix.df.b.fung,method="bray")

#create vectors with the NMS attributes
NMS_coordinates_fung<-scores(fung_mds,display="sites")
fung_NMS_axes<-as.data.frame(NMS_coordinates_fung)
#write.csv(Bact_NMS_axes,"Bact_NMS_axes.csv")
NMS_OTUscores_fung<-scores(fung_mds,display="species")

#add NMS coordinates to the mapping file
#create dataframe with NMS Coordinates and Mapping File information
#add the proportion of total sequences of each vector to the "for_plotting" object below
### use only the 10 most common orders
for_ploting_fung<-as.data.frame(cbind(NMS_coordinates_fung,fungi_ord_map))
str(for_ploting_fung) #20 obs. of 27 variables

#now create the vectors
fit1_fung<-envfit(fung_mds, fungi_ord_map[,2:9,11], perm = 999, na.rm = TRUE)
fit1_fung
str(for_ploting_fung)

#create the figure -> NMDS
tiff(file="Fungi_NMS_experiment_20220922.tiff", width = 1200*3.25, height = 1200*2, res = 1200)
par(mar=c(2,2,1,1))
plot(for_ploting_fung$NMDS1 ~ for_ploting_fung$NMDS2,
     xlab = "Fungal NMS1",
     ylab = "Fungal NMS2",
     font=2,
     font.lab=2,
     cex.axis=.7,
     cex.lab=.7,
     xlim = c(-1.5,1.5),
     mgp=c(1.15,.4,0),
     col = c("blue3","light blue","red4","lightcoral")[as.factor(for_ploting_fung$groundcontact_treatment)],
     # bg = c("blue3","light blue","red4","lightcoral")[as.factor(unique(for_ploting$groundcontact_treatment))],
     #set the point types (plotting 'character'=pch) and colors based on your data (4 groups)
     pch = c(0,15,1,16)[as.factor(for_ploting_fung$groundcontact_treatment)], cex=.8, 
     data = for_ploting_fung)

#HULLS
#ordiellipse(fung_mds, group=for_ploting_fung$groundcontact_treatment,kind = "se", 
            conf=0.95, lwd=1.9, 
            #set the linetypes and points based on your data (not the same as here)
            lty = c("solid","solid","solid","solid"), 
            col = c("blue3","light blue","red4","lightcoral")[as.factor(unique(for_ploting_fung$groundcontact_treatment))])

# names_inter<-gsub("o__","",names(fit1$vectors$arrows[,1]))
# bact_order_labels<-gsub("o__.","",names_inter)
# #Plot labels 
# plot(fit1,font=1, p.max = 0.05, col = "black", cex=0.6,
#      labels = bact_order_labels)
legend(
  x ="bottomleft",
  legend = c("down-connected","down-separated","susp.-separated","susp.-connected"), # for readability of legend
  pch = c(0,15,1,16),
  col =c("blue3","light blue","red4","lightcoral"),
  cex = 0.49 # scale the legend to look attractively sized
)

dev.off()
getwd()

###################### 
#### PLOT RAW ABUNDANCES ####
#plot the data raw values of relative abundance within each group
## individual points for each sample rather than a boxplot (so replace "geom_boxplot" with "geom_point")

#data frames to be used
summarized_orders_bacteria <- as.data.frame(summarized_orders_bacteria)
summarized_orders_fungi <- as.data.frame(summarized_orders_fungi)

#filter 10 most common orders 
#bact

#o__Rhizobiales              2.5973742
#o__Actinomycetales          1.3894967
#o__Rhodospirillales         1.3301032
#o__Xanthomonadales          1.2269459
#o__Myxococcales             1.0768990
#o__.Saprospirales.          0.9149734 # not found in taxa_sum_order_bact_t
#o__Gemmatales               0.9112223
#o__Burkholderiales          0.6398875
#o__Cytophagales             0.6123789
#o__.Chthoniobacterales.     0.5970616 # not found in taxa_sum_order_bact_t

taxa_sum_order_bact_10 <- as.data.frame(summarized_orders_bacteria[,c("X","o__Rhizobiales",
                               "o__Actinomycetales","o__Rhodospirillales",
                               "o__Xanthomonadales","o__Myxococcales",
                               "o__.Saprospirales.","o__Gemmatales",
                               "o__Burkholderiales", "o__Cytophagales",
                               "o__.Chthoniobacterales.")])
#names 

colnames(taxa_sum_order_bact_10) <- c("sample","Rhizobiales",
                                        "Actinomycetales","Rhodospirillales",
                                        "Xanthomonadales","Myxococcales","Saprospirales",
                                        "Gemmatales","Burkholderiales", "Cytophagales",
                                        "Chthoniobacterales")

#bring treatment column 

mapping_bactdiversity_treat <- mapping_bactdiversity[,c("X","groundcontact_treatment")]

colnames(mapping_bactdiversity_treat) <- c("sample", "treatment")

taxa_sum_order_bact_10_t <- left_join(mapping_bactdiversity_treat, 
                                      taxa_sum_order_bact_10, by='sample')

# fungi 

#10 most common orders fungi 

#o__Agaricales                             1.1653912
#o__Xylariales                             0.7166241
#o__Hypocreales                            0.6707058
#o__Chaetothyriales                        0.6611395
#o__Pezizomycotina_ord_Incertae_sedis      0.5801446
#o__Magnaporthales                         0.4936224
#o__Capnodiales                            0.4455782
#o__Pleosporales                           0.3843537
#o__Polyporales                            0.2869898

taxa_sum_order_fungi_9 <- as.data.frame(summarized_orders_fungi[,c("X","o__Agaricales",
                                                                      "o__Xylariales",
                                                                    "o__Hypocreales",
                                                                    "o__Chaetothyriales",
                                                                    "o__Pezizomycotina_ord_Incertae_sedis",
                                                                    "o__Magnaporthales",
                                                                    "o__Capnodiales",
                                                                    "o__Pleosporales",
                                                                    "o__Polyporales")])
#names 

colnames(taxa_sum_order_fungi_9) <- c("sample","Agaricales",
                                      "Xylariales","Hypocreales",
                                      "Chaetothyriales","Pezizomycotina_incertae_sedis",
                                      "Magnaporthales","Capnodiales",
                                      "Pleosporales", "Polyporales")

#bring treatment column 

mapping_fungidiversity_treat <- mapping_fungidiversity[,c("X","groundcontact_treatment")]

colnames(mapping_fungidiversity_treat) <- c("sample", "treatment")

taxa_sum_order_fungi_9_t <- left_join(mapping_fungidiversity_treat, 
                                      taxa_sum_order_fungi_9, by='sample')


####plot raw abundances orders bacteria####
         #"Rhizobiales"

bact_order_raw_abun_rhizob <- ggplot(data = taxa_sum_order_bact_10_t, aes(x =treatment,
                        y = Rhizobiales)) +
  geom_point(position =position_dodge(width =.1)) +
  theme(legend.position=c(.88,.92),legend.title = element_blank(),
        legend.background= element_blank(),legend.key=element_blank(),
        legend.text=element_text(family = "Arial",face="bold",size=10)) +
  scale_fill_manual(values = c("white","grey50")) +
  scale_x_discrete(name=element_blank(),labels=c("bottom half","downed","suspended","top half")) +
  scale_y_continuous(name = "Rhizobiales abundance")

bact_order_raw_abun_rhizob

#save figure
setwd('C:/Users/angel/OneDrive/3_microbial_diversity_20220816/Exploration_results/orders_abundance/bacteria')
ggsave("1_bact_order_raw_abun_rhizob.tiff",bact_order_raw_abun_rhizob,dpi=600,
       width=3.25,height=2,scale=1.8,compression = "lzw")

          #"Actinomycetales"

bact_order_raw_abun_actin <- ggplot(data = taxa_sum_order_bact_10_t, aes(x =treatment,
                                                                          y = Actinomycetales)) +
  geom_point(position =position_dodge(width =.1)) +
  theme(legend.position=c(.88,.92),legend.title = element_blank(),
        legend.background= element_blank(),legend.key=element_blank(),
        legend.text=element_text(family = "Arial",face="bold",size=10)) +
  scale_fill_manual(values = c("white","grey50")) +
  scale_x_discrete(name=element_blank(),labels=c("bottom half","downed","suspended","top half")) +
  scale_y_continuous(name = "Actinomycetales abundance")

bact_order_raw_abun_actin

#save figure
setwd('C:/Users/angel/OneDrive/3_microbial_diversity_20220816/Exploration_results/orders_abundance/bacteria')
ggsave("2_bact_order_raw_abun_actin.tiff",bact_order_raw_abun_actin,dpi=600,
       width=3.25,height=2,scale=1.8,compression = "lzw")

          #"Rhodospirillales"

bact_order_raw_abun_rhod <- ggplot(data = taxa_sum_order_bact_10_t, aes(x =treatment,
                                                                         y = Rhodospirillales)) +
  geom_point(position =position_dodge(width =.1)) +
  theme(legend.position=c(.88,.92),legend.title = element_blank(),
        legend.background= element_blank(),legend.key=element_blank(),
        legend.text=element_text(family = "Arial",face="bold",size=10)) +
  scale_fill_manual(values = c("white","grey50")) +
  scale_x_discrete(name=element_blank(),labels=c("bottom half","downed","suspended","top half")) +
  scale_y_continuous(name = "Rhodospirillales abundance")

bact_order_raw_abun_rhod

#save figure
setwd('C:/Users/angel/OneDrive/3_microbial_diversity_20220816/Exploration_results/orders_abundance/bacteria')
ggsave("3_bact_order_raw_abun_rhod.tiff",bact_order_raw_abun_rhod,dpi=600,
       width=3.25,height=2,scale=1.8,compression = "lzw")

          #"Xanthomonadales"

bact_order_raw_abun_xan <- ggplot(data = taxa_sum_order_bact_10_t, aes(x =treatment,
                                                                        y = Xanthomonadales)) +
  geom_point(position =position_dodge(width =.1)) +
  theme(legend.position=c(.88,.92),legend.title = element_blank(),
        legend.background= element_blank(),legend.key=element_blank(),
        legend.text=element_text(family = "Arial",face="bold",size=10)) +
  scale_fill_manual(values = c("white","grey50")) +
  scale_x_discrete(name=element_blank(),labels=c("bottom half","downed","suspended","top half")) +
  scale_y_continuous(name = "Xanthomonadales abundance")

bact_order_raw_abun_xan

#save figure
setwd('C:/Users/angel/OneDrive/3_microbial_diversity_20220816/Exploration_results/orders_abundance/bacteria')
ggsave("4_bact_order_raw_abun_xan.tiff",bact_order_raw_abun_xan,dpi=600,
       width=3.25,height=2,scale=1.8,compression = "lzw")

          #"Myxococcales"

bact_order_raw_abun_myx <- ggplot(data = taxa_sum_order_bact_10_t, aes(x =treatment,
                                                                       y = Myxococcales)) +
  geom_point(position =position_dodge(width =.1)) +
  theme(legend.position=c(.88,.92),legend.title = element_blank(),
        legend.background= element_blank(),legend.key=element_blank(),
        legend.text=element_text(family = "Arial",face="bold",size=10)) +
  scale_fill_manual(values = c("white","grey50")) +
  scale_x_discrete(name=element_blank(),labels=c("bottom half","downed","suspended","top half")) +
  scale_y_continuous(name = "Myxococcales abundance")

bact_order_raw_abun_myx

#save figure
setwd('C:/Users/angel/OneDrive/3_microbial_diversity_20220816/Exploration_results/orders_abundance/bacteria')
ggsave("5_bact_order_raw_abun_myx.tiff",bact_order_raw_abun_myx,dpi=600,
       width=3.25,height=2,scale=1.8,compression = "lzw")

          #"Saprospirales"

bact_order_raw_abun_sap <- ggplot(data = taxa_sum_order_bact_10_t, aes(x =treatment,
                                                                       y = Saprospirales)) +
  geom_point(position =position_dodge(width =.1)) +
  theme(legend.position=c(.88,.92),legend.title = element_blank(),
        legend.background= element_blank(),legend.key=element_blank(),
        legend.text=element_text(family = "Arial",face="bold",size=10)) +
  scale_fill_manual(values = c("white","grey50")) +
  scale_x_discrete(name=element_blank(),labels=c("bottom half","downed","suspended","top half")) +
  scale_y_continuous(name = "Saprospirales abundance")

bact_order_raw_abun_xan

#save figure
setwd('C:/Users/angel/OneDrive/3_microbial_diversity_20220816/Exploration_results/orders_abundance/bacteria')
ggsave("6_bact_order_raw_abun_sap.tiff",bact_order_raw_abun_sap,dpi=600,
       width=3.25,height=2,scale=1.8,compression = "lzw")

          #"Gemmatales"

bact_order_raw_abun_gem <- ggplot(data = taxa_sum_order_bact_10_t, aes(x =treatment,
                                                                       y = Gemmatales)) +
  geom_point(position =position_dodge(width =.1)) +
  theme(legend.position=c(.88,.92),legend.title = element_blank(),
        legend.background= element_blank(),legend.key=element_blank(),
        legend.text=element_text(family = "Arial",face="bold",size=10)) +
  scale_fill_manual(values = c("white","grey50")) +
  scale_x_discrete(name=element_blank(),labels=c("bottom half","downed","suspended","top half")) +
  scale_y_continuous(name = "Gemmatales abundance")

bact_order_raw_abun_gem

#save figure
setwd('C:/Users/angel/OneDrive/3_microbial_diversity_20220816/Exploration_results/orders_abundance/bacteria')
ggsave("7_bact_order_raw_abun_gem.tiff",bact_order_raw_abun_gem,dpi=600,
       width=3.25,height=2,scale=1.8,compression = "lzw")

          #"Burkholderiales" 

bact_order_raw_abun_bur <- ggplot(data = taxa_sum_order_bact_10_t, aes(x =treatment,
                                                                       y = Burkholderiales)) +
  geom_point(position =position_dodge(width =.1)) +
  theme(legend.position=c(.88,.92),legend.title = element_blank(),
        legend.background= element_blank(),legend.key=element_blank(),
        legend.text=element_text(family = "Arial",face="bold",size=10)) +
  scale_fill_manual(values = c("white","grey50")) +
  scale_x_discrete(name=element_blank(),labels=c("bottom half","downed","suspended","top half")) +
  scale_y_continuous(name = "Burkholderiales abundance")

bact_order_raw_abun_bur

#save figure
setwd('C:/Users/angel/OneDrive/3_microbial_diversity_20220816/Exploration_results/orders_abundance/bacteria')
ggsave("8_bact_order_raw_abun_bur.tiff",bact_order_raw_abun_bur,dpi=600,
       width=3.25,height=2,scale=1.8,compression = "lzw")

          #"Cytophagales"

bact_order_raw_abun_cyt <- ggplot(data = taxa_sum_order_bact_10_t, aes(x =treatment,
                                                                       y = Cytophagales)) +
  geom_point(position =position_dodge(width =.1)) +
  theme(legend.position=c(.88,.92),legend.title = element_blank(),
        legend.background= element_blank(),legend.key=element_blank(),
        legend.text=element_text(family = "Arial",face="bold",size=10)) +
  scale_fill_manual(values = c("white","grey50")) +
  scale_x_discrete(name=element_blank(),labels=c("bottom half","downed","suspended","top half")) +
  scale_y_continuous(name = "Cytophagales abundance")

bact_order_raw_abun_cyt

#save figure
setwd('C:/Users/angel/OneDrive/3_microbial_diversity_20220816/Exploration_results/orders_abundance/bacteria')
ggsave("9_bact_order_raw_abun_cyt.tiff",bact_order_raw_abun_cyt,dpi=600,
       width=3.25,height=2,scale=1.8,compression = "lzw")

          #"Chthoniobacterales"

bact_order_raw_abun_cht <- ggplot(data = taxa_sum_order_bact_10_t, aes(x =treatment,
                                                                       y = Chthoniobacterales)) +
  geom_point(position =position_dodge(width =.1)) +
  theme(legend.position=c(.88,.92),legend.title = element_blank(),
        legend.background= element_blank(),legend.key=element_blank(),
        legend.text=element_text(family = "Arial",face="bold",size=10)) +
  scale_fill_manual(values = c("white","grey50")) +
  scale_x_discrete(name=element_blank(),labels=c("bottom half","downed","suspended","top half")) +
  scale_y_continuous(name = "Chthoniobacterales abundance")

bact_order_raw_abun_cht

#save figure
setwd('C:/Users/angel/OneDrive/3_microbial_diversity_20220816/Exploration_results/orders_abundance/bacteria')
ggsave("10_bact_order_raw_abun_cht.tiff",bact_order_raw_abun_cht,dpi=600,
       width=3.25,height=2,scale=1.8,compression = "lzw")
  
####plot raw abundances orders fungi####

           #"Agaricales"

fungi_order_raw_abun_aga <- ggplot(data = taxa_sum_order_fungi_9_t, aes(x =treatment,
                                                                       y = Agaricales)) +
  geom_point(position =position_dodge(width =.1)) +
  theme(legend.position=c(.88,.92),legend.title = element_blank(),
        legend.background= element_blank(),legend.key=element_blank(),
        legend.text=element_text(family = "Arial",face="bold",size=10)) +
  scale_fill_manual(values = c("white","grey50")) +
  scale_x_discrete(name=element_blank(),labels=c("bottom half","downed","suspended","top half")) +
  scale_y_continuous(name = "Agaricales abundance")

fungi_order_raw_abun_aga

#save figure
setwd('C:/Users/angel/OneDrive/3_microbial_diversity_20220816/Exploration_results/orders_abundance/fungi')
ggsave("1_fungi_order_raw_abun_cht.tiff",fungi_order_raw_abun_aga,dpi=600,
       width=3.25,height=2,scale=1.8,compression = "lzw")

          #"Xylariales"

fungi_order_raw_abun_xyl <- ggplot(data = taxa_sum_order_fungi_9_t, aes(x =treatment,
                                                                        y = Xylariales)) +
  geom_point(position =position_dodge(width =.1)) +
  theme(legend.position=c(.88,.92),legend.title = element_blank(),
        legend.background= element_blank(),legend.key=element_blank(),
        legend.text=element_text(family = "Arial",face="bold",size=10)) +
  scale_fill_manual(values = c("white","grey50")) +
  scale_x_discrete(name=element_blank(),labels=c("bottom half","downed","suspended","top half")) +
  scale_y_continuous(name = "Xylariales abundance")

fungi_order_raw_abun_xyl

#save figure
setwd('C:/Users/angel/OneDrive/3_microbial_diversity_20220816/Exploration_results/orders_abundance/fungi')
ggsave("2_fungi_order_raw_abun_xyl.tiff",fungi_order_raw_abun_xyl,dpi=600,
       width=3.25,height=2,scale=1.8,compression = "lzw")

          #"Hypocreales"

fungi_order_raw_abun_hyp <- ggplot(data = taxa_sum_order_fungi_9_t, aes(x =treatment,
                                                                        y = Hypocreales)) +
  geom_point(position =position_dodge(width =.1)) +
  theme(legend.position=c(.88,.92),legend.title = element_blank(),
        legend.background= element_blank(),legend.key=element_blank(),
        legend.text=element_text(family = "Arial",face="bold",size=10)) +
  scale_fill_manual(values = c("white","grey50")) +
  scale_x_discrete(name=element_blank(),labels=c("bottom half","downed","suspended","top half")) +
  scale_y_continuous(name = "Hypocreales abundance")

fungi_order_raw_abun_hyp

#save figure
setwd('C:/Users/angel/OneDrive/3_microbial_diversity_20220816/Exploration_results/orders_abundance/fungi')
ggsave("3_fungi_order_raw_abun_hyp.tiff",fungi_order_raw_abun_hyp,dpi=600,
       width=3.25,height=2,scale=1.8,compression = "lzw")

          #"Chaetothyriales"

fungi_order_raw_abun_cha <- ggplot(data = taxa_sum_order_fungi_9_t, aes(x =treatment,
                                                                        y = Chaetothyriales)) +
  geom_point(position =position_dodge(width =.1)) +
  theme(legend.position=c(.88,.92),legend.title = element_blank(),
        legend.background= element_blank(),legend.key=element_blank(),
        legend.text=element_text(family = "Arial",face="bold",size=10)) +
  scale_fill_manual(values = c("white","grey50")) +
  scale_x_discrete(name=element_blank(),labels=c("bottom half","downed","suspended","top half")) +
  scale_y_continuous(name = "Chaetothyriales abundance")

fungi_order_raw_abun_cha

#save figure
setwd('C:/Users/angel/OneDrive/3_microbial_diversity_20220816/Exploration_results/orders_abundance/fungi')
ggsave("4_fungi_order_raw_abun_cha.tiff",fungi_order_raw_abun_cha,dpi=600,
       width=3.25,height=2,scale=1.8,compression = "lzw")

          #"Pezizomycotina_incertae_sedis"
 
fungi_order_raw_abun_pez <- ggplot(data = taxa_sum_order_fungi_9_t, aes(x =treatment,
                                                                        y = Pezizomycotina_incertae_sedis)) +
  geom_point(position =position_dodge(width =.1)) +
  theme(legend.position=c(.88,.92),legend.title = element_blank(),
        legend.background= element_blank(),legend.key=element_blank(),
        legend.text=element_text(family = "Arial",face="bold",size=10)) +
  scale_fill_manual(values = c("white","grey50")) +
  scale_x_discrete(name=element_blank(),labels=c("bottom half","downed","suspended","top half")) +
  scale_y_continuous(name = "Pezizomycotina Incertae sedis abundance")

fungi_order_raw_abun_pez

#save figure
setwd('C:/Users/angel/OneDrive/3_microbial_diversity_20220816/Exploration_results/orders_abundance/fungi')
ggsave("5_fungi_order_raw_abun_pez.tiff",fungi_order_raw_abun_pez,dpi=600,
       width=3.25,height=2,scale=1.8,compression = "lzw")

          #"Magnaporthales"

fungi_order_raw_abun_mag <- ggplot(data = taxa_sum_order_fungi_9_t, aes(x =treatment,
                                                                        y = Magnaporthales)) +
  geom_point(position =position_dodge(width =.1)) +
  theme(legend.position=c(.88,.92),legend.title = element_blank(),
        legend.background= element_blank(),legend.key=element_blank(),
        legend.text=element_text(family = "Arial",face="bold",size=10)) +
  scale_fill_manual(values = c("white","grey50")) +
  scale_x_discrete(name=element_blank(),labels=c("bottom half","downed","suspended","top half")) +
  scale_y_continuous(name = "Magnaporthales abundance")

fungi_order_raw_abun_mag

#save figure
setwd('C:/Users/angel/OneDrive/3_microbial_diversity_20220816/Exploration_results/orders_abundance/fungi')
ggsave("6_fungi_order_raw_abun_mag.tiff",fungi_order_raw_abun_mag,dpi=600,
       width=3.25,height=2,scale=1.8,compression = "lzw")

          #"Capnodiales"

fungi_order_raw_abun_cap <- ggplot(data = taxa_sum_order_fungi_9_t, aes(x =treatment,
                                                                        y = Capnodiales)) +
  geom_point(position =position_dodge(width =.1)) +
  theme(legend.position=c(.88,.92),legend.title = element_blank(),
        legend.background= element_blank(),legend.key=element_blank(),
        legend.text=element_text(family = "Arial",face="bold",size=10)) +
  scale_fill_manual(values = c("white","grey50")) +
  scale_x_discrete(name=element_blank(),labels=c("bottom half","downed","suspended","top half")) +
  scale_y_continuous(name = "Capnodiales abundance")

fungi_order_raw_abun_cap

#save figure
setwd('C:/Users/angel/OneDrive/3_microbial_diversity_20220816/Exploration_results/orders_abundance/fungi')
ggsave("7_fungi_order_raw_abun_cap.tiff",fungi_order_raw_abun_cap,dpi=600,
       width=3.25,height=2,scale=1.8,compression = "lzw")

          #"Pleosporales"

fungi_order_raw_abun_ple <- ggplot(data = taxa_sum_order_fungi_9_t, aes(x =treatment,
                                                                        y = Pleosporales)) +
  geom_point(position =position_dodge(width =.1)) +
  theme(legend.position=c(.88,.92),legend.title = element_blank(),
        legend.background= element_blank(),legend.key=element_blank(),
        legend.text=element_text(family = "Arial",face="bold",size=10)) +
  scale_fill_manual(values = c("white","grey50")) +
  scale_x_discrete(name=element_blank(),labels=c("bottom half","downed","suspended","top half")) +
  scale_y_continuous(name = "Pleosporales abundance")

fungi_order_raw_abun_ple

#save figure
setwd('C:/Users/angel/OneDrive/3_microbial_diversity_20220816/Exploration_results/orders_abundance/fungi')
ggsave("8_fungi_order_raw_abun_ple.tiff",fungi_order_raw_abun_ple,dpi=600,
       width=3.25,height=2,scale=1.8,compression = "lzw")

          #"Polyporales"

fungi_order_raw_abun_pol <- ggplot(data = taxa_sum_order_fungi_9_t, aes(x =treatment,
                                                                        y = Polyporales)) +
  geom_point(position =position_dodge(width =.1)) +
  theme(legend.position=c(.88,.92),legend.title = element_blank(),
        legend.background= element_blank(),legend.key=element_blank(),
        legend.text=element_text(family = "Arial",face="bold",size=10)) +
  scale_fill_manual(values = c("white","grey50")) +
  scale_x_discrete(name=element_blank(),labels=c("bottom half","downed","suspended","top half")) +
  scale_y_continuous(name = "Polyporales abundance")

fungi_order_raw_abun_pol

#save figure
setwd('C:/Users/angel/OneDrive/3_microbial_diversity_20220816/Exploration_results/orders_abundance/fungi')
ggsave("9_fungi_order_raw_abun_pol.tiff",fungi_order_raw_abun_pol,dpi=600,
       width=3.25,height=2,scale=1.8,compression = "lzw")

