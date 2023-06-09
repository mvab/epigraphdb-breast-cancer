# build metadata
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(TwoSampleMR)
ao <- available_outcomes()
source('01_MR_related/scripts/app1_MR-EvE_app/functions.R')


all_mr_tests <- read_tsv("01_MR_related/results/mr_evidence_outputs/tidy_traits_by_catV3.tsv") %>% 
  mutate(exposure_cat = ifelse(exposure.trait == "Albumin", "Metabolites", exposure_cat))  # 202



meta <- ao %>% filter(id %in% all_mr_tests$exposure.id) %>% 
               left_join(all_mr_tests %>% select(exposure.id, exposure_cat) %>% distinct(), 
                        by = c("id" = "exposure.id")) %>% 
              select(exposure_cat, everything()) %>% 
              arrange(exposure_cat) %>% 
              select('exposure_cat', "id" , "trait",
                     "year" , "author" , "consortium" , 
                     "population" ,"sex",  "sample_size" , "ncase", "ncontrol" ,
                     "unit"   ,   "sd" , "nsnp" , 
                     "ontology" ,   "note" ,"pmid" ) %>% 
              mutate(sample_size = ifelse(is.na(sample_size) & author == "Neale lab", 361194, sample_size))

length(unique(meta$id)) # 202

sums_mol <- meta  %>% filter(exposure_cat %in% c("Proteins", "Metabolites", "Lipids")) %>% 
  select( "exposure_cat", "year" , "author" , "consortium" , "sample_size", "sex", 'nsnp') %>% 
  group_by(exposure_cat, year , author , consortium, sex) %>% 
  summarise(ss_mean = mean(sample_size), nsnp_mean = mean(nsnp),  n=n()) %>% mutate(ss_mean=round(ss_mean,0))
# 147

sums_life <- meta  %>% filter(exposure_cat %in% c("Alcohol" , "Antrophometric" ,
                                                   "Diet and supplements", "Drugs",  "Physical activity"   , 
                                                   "Reproductive", "Sleep" , "Smoking",
                                                   "Other biomarkers")) %>% 
  select( "exposure_cat", "year" , "author" , "consortium" , "sample_size", "sex", 'nsnp') %>% 
  group_by(exposure_cat,  year , author , consortium, sex) %>% 
  summarise(ss_mean = mean(sample_size), nsnp_mean = mean(nsnp),  n=n()) %>% mutate(ss_mean=round(ss_mean,0))
# 55

sums_life_sum <- meta  %>% filter(exposure_cat %in% c("Alcohol" , "Antrophometric" ,
                                                  "Diet and supplements", "Drugs",  "Physical activity"   , 
                                                  "Reproductive", "Sleep" , "Smoking",
                                                  "Other biomarkers")) %>% 
  select( "exposure_cat", "year" , "author" , "consortium" , "sample_size", "sex", 'nsnp') %>% 
  group_by(  year , author , consortium, sex) %>% 
  summarise(ss_mean = mean(sample_size), nsnp_mean = mean(nsnp),  n=n()) %>% mutate(ss_mean=round(ss_mean,0))

