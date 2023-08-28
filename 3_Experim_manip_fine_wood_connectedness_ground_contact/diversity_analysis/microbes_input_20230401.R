#input bacterial and fungal data
#perform exploratory analyses
#Edited by EM Gora 2021-6-30

#install packages
install.packages("devtools")
devtools::install_github("leffj/mctoolsr")
install.packages("vegan")
if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')
BiocManager::install("AnnotationForge", force = TRUE)
install.packages("biom")
install.packages("devtools")
install.Rtools(check = TRUE, check_r_update = TRUE, GUI = TRUE)
install_github("biom", "joey711")
BiocManager::install("biomformat")
install.packages("biomformat", force = TRUE)
devtools::install_github("joey711/biomformat")
install_github("biom", "joey711")


install.packages('githubinstall')
library(githubinstall)
        githubinstall("biom")

#load packages
library(mctoolsr)
library(vegan)
library(dplyr)
library(ggplot2)
library(grid)
library("devtools")

library(biom)
        
        if (!require("BiocManager", quietly = TRUE))
          install.packages("BiocManager")
        BiocManager::install(version = "3.16")
        
library("BiocManager")
library("biomformat")

#remove
rm(list=ls())

##First figure out how to import data using mctoolsr
#data location for fungi
fung_table_imp<-'C:\\Users\\angel\\OneDrive\\Independent_project_wood\\3_microbial_diversity_20220816\\Analysis\\Fungal_OTUtable.biom'
#data location for bacteria
bact_table_imp<-'C:\\Users\\angel\\OneDrive\\Independent_project_wood\\3_microbial_diversity_20220816\\Analysis\\otu_table_mc2_w_tax_no_pynast_failures.biom'
# update the mapping file to only have the priority effects
### make sure that all of covariates are present for the covariate effects
map_file<-'C:\\Users\\angel\\OneDrive\\Independent_project_wood\\3_microbial_diversity_20220816\\Analysis\\GC_mapping_file_2021-8-31.txt'

fung_input<-load_taxa_table(fung_table_imp,map_file)
bact_input<-load_taxa_table(bact_table_imp,map_file)

#what samples are missing?
table(fung_input$map_loaded$groundcontact_treatment)
table(fung_input$map_loaded$groundcontact_groupcat)
table(bact_input$map_loaded$groundcontact_treatment)
table(bact_input$map_loaded$groundcontact_groupcat)

table(fung_input$map_loaded$gc_suspended,
      fung_input$map_loaded$gc_connection)
table(bact_input$map_loaded$gc_suspended,
      bact_input$map_loaded$gc_connection)

#the middle timeframe was not imported well
###I might want to check for errors in the mapping file

str(bact_input)
str(fung_input)
#these data have 3 components: 
#1. abundance data ($data_loaded);
#2. mapping file ($map_loaded); 
#3. taxonomy file ($taxonomy_loaded)

#how many samples were imported?
length(fung_input$data_loaded)#14 samples
length(bact_input$data_loaded)#20 samples

#separate bacteria and archaea
#drop anything not assigned to either domain

#//////////////////////
## 1. ORGANIZE DATA ####
#//////////////////////
#set directory
setwd('C:\\Users\\angel\\OneDrive\\Independent_project_wood\\3_microbial_diversity_20220816\\Exploration_results')
#Filter values to just have the Kingdom Archaea
all_Archaea<-filter_taxa_from_input(bact_input, taxa_to_keep = 'k__Archaea')
taxa_sum_Arch_gen<-summarize_taxonomy(all_Archaea,level=6,relative = FALSE,report_higher_tax = FALSE)
plot_taxa_bars(taxa_sum_Arch_gen,all_Archaea$map_loaded,c('groundcontact_treatment'),10)
#write.csv(as.matrix(taxa_sum_Arch_gen),"Arch_gen_summary.csv")
taxa_sum_Arch_ord<-summarize_taxonomy(all_Archaea,level=4,relative = FALSE,report_higher_tax = FALSE)
plot_taxa_bars(taxa_sum_Arch_ord,all_Archaea$map_loaded,'groundcontact_treatment',10)
#write.csv(as.matrix(taxa_sum_Arch_ord),"Arch_ord_summary.csv")
taxa_sum_Arch_class<-summarize_taxonomy(all_Archaea,level=3,relative = FALSE,report_higher_tax = FALSE)
plot_taxa_bars(taxa_sum_Arch_class,all_Archaea$map_loaded,'groundcontact_treatment',10)
#write.csv(as.matrix(taxa_sum_Arch_class),"Arch_class_summary.csv")

#Now filter out Archaea so that I can run Bacterial analyses
all_Bacteria_inter<-filter_taxa_from_input(bact_input, taxa_to_keep = 'k__Bacteria')
all_Bacteria<-filter_taxa_from_input(all_Bacteria_inter, taxa_to_remove = 'c__Chloroplast')

#look at data summary - the number of sequences per sample
sort(colSums(all_Bacteria$data_loaded))
#all samples have at least 3265 sequences

#look at data summary for FUNGI - the number of sequences per sample
sort(colSums(fung_input$data_loaded))
#all samples have at least 33 sequences

#create a file from the mapping file that is 3 columns of relevant groupings
# all_dataloaded<-as.data.frame(all_Bacteria$map_loaded)
#View(all_dataloaded)

#create a file from the mapping file for FUNGI that is 3 columns of relevant groupings
# all_dataloadedfungi<-as.data.frame(fung_input$map_loaded)
#View(all_dataloadedfungi)

#filter singletons from the Bacteria (debemos poner 1/19 porque asi eliminamos los organismos 
# que estan solamente en una de las 20 muestras de bacterias existentes, de ahi para abjo elimina)
input_nosingles<-filter_taxa_from_input(all_Bacteria,1/19)
#10,369 taxa were removed

#filter singletons from the FUNGI 
input_nosinglesfungi<-filter_taxa_from_input(fung_input,1/13)
#4170 taxa were removed

sort(colSums(input_nosingles$data_loaded))
#after removing singletons, all samples have at least 3199 sequences

sort(colSums(input_nosinglesfungi$data_loaded))
#FOR FUNGI TOO: after removing singletons, all samples have at least 32 sequences

# FUNGI: Samples GC08  GC15  GC07 were to tiny (few OTUS) -> we remove them like this:
# we create another data frame called noweirdreplicatesfungi extracting dataloaded 
# from input_nosinglesfungi and we remove tiny samples from colums 9,13,14 which are GC08  GC15  GC07
noweirdreplicatesfungi <- as.data.frame(input_nosinglesfungi$data_loaded)
noweirdreplicatesfungi <- noweirdreplicatesfungi[,-c(9, 13:14)]
#we replace the all data frame dataloaded to noweirdreplicatesfungi. It has the same name dataloaded
input_nosinglesfungi$data_loaded <- noweirdreplicatesfungi
#Now we have 11 samples for fungi

#now rarefy data - only include the smallest number of sequences per sample
#input_bact_rarfy<-single_rarefy(input_nosingles,3199)
#It drops the samples that do not have sufficient sequences
#this means we reduced each sample to a random selection of 3301 sequences

#FUNGI: now rarefy data - only include the smallest number of sequences per sample
#input_fungi_rarfy<-single_rarefy(input_nosinglesfungi,4704)
#It drops the samples that do not have sufficient sequences
#this means we reduced each sample to a random selection of 4704 sequences. THIS IS ALMOST THE SAME
# WE DID WITH ABOVE WITH noweirdreplicatesfungi (="rarefaccion a pie")

# confirm that all samples have the same number of sequences (=rarefied)
sort(colSums(input_bact_rarfy$data_loaded))
sort(colSums(input_fungi_rarfy$data_loaded))

#Save lists of bacteria and fungi for next scripts
#saveRDS(input_bact_rarfy, file="input_bact_rarfy.RData")
#saveRDS(input_fungi_rarfy, file="input_fungi_rarfy.RData")

#Continues on R script "diversity_analysis" 
#'D:/indep_project_wood_currentlyworking/3_microbial_diversity/Analysis'