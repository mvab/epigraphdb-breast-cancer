library(tidyverse)
library(epigraphdb)
source("helper_functions.R")
source("01_MR_related/scripts/mr_related_functions.R")
source("01_MR_related/scripts/app1_MR-EvE_app/functions.R")

# testing example
exposure = 'ukb-d-30760_irnt'
outcome = "ieu-a-1126"
pval_threshold = 1e-01
#test<-query_epigraphdb_as_table(mediator_query)
#xx<-tidy_conf_query_output(test, type = "mediator")



#### 

# read effect direction matrix from 03..

traits_df <- read_tsv("01_MR_related/results/mr_evidence_outputs/trait_manual_ivw_subtypes_mergedV3.tsv") %>% 
  select(id.exposure, contains('ieu')) %>% distinct()# only for BCAC 2017 cols, because mr-eve is available only for them

# extract pairs of exp-out with effect
traits_all <- bind_rows(
  traits_df %>%filter(`ieu-a-1126` !=0) %>%  select(id.exposure) %>% mutate(id.outcome = "ieu-a-1126"),
  traits_df %>%filter(`ieu-a-1127` !=0) %>%  select(id.exposure) %>% mutate(id.outcome = "ieu-a-1127"),
  traits_df %>%filter(`ieu-a-1128` !=0) %>%  select(id.exposure) %>% mutate(id.outcome = "ieu-a-1128")) %>% distinct()
traits_all %>% pull(id.exposure) %>%  unique() %>% length() # 105


# for each pair collect intermediates (all 4 scenarios)
all_results <- tibble()
for (i in 1:length(traits_all$id.exposure)){
  
  out <- query_and_tidy_conf(exposure = traits_all$id.exposure[i], 
                             outcome = traits_all$id.outcome[i], 
                             pval_threshold =  1, ##### keeping all; going to use FDR on both steps to identify reliable results
                             mediator = T)  # only mediators
  #add FDR corrected p-val to results
  out <- out %>% 
  arrange(r1.pval) %>% mutate(r1.qval = p.adjust(r1.pval, method = "BH")) %>% 
  arrange(r3.pval) %>% mutate(r3.qval = p.adjust(r3.pval, method = "BH")) 
  
  all_results<- bind_rows(all_results, out)
}
length(unique(all_results$exposure.id)) # 103 it's ok : 96 will do too :saved 103 as med_extracted_all_r3V3_backup_with_103 

# keep only those pairs, where each step q-value passed 0.05
all_results_fdr <- all_results %>% filter (r1.qval <0.05 & r3.qval <0.05)


# add trait category to intermediate items
all_results_trait_cats <- 
  all_results_fdr %>%  
  select(exposure.trait = med.trait, exposure.id = med.id) %>% 
  create_exposure_categories() %>% 
  select(med_cat = exposure_cat, med.trait = exposure.trait, med.id = exposure.id ) %>% distinct() %>% 
  filter(!grepl("arm|leg", med.trait, ignore.case = T)) %>% 
  filter(med_cat %in% c("Antrophometric" ,"Physical activity" ,   "Diet and supplements", "Reproductive"  , "Alcohol" , "Smoking" ,
                       "Metabolites" ,   "Proteins", "Other biomarkers" ,  "Sleep"  ,  "Drugs" ,   "Lipids"   )) 

# join with cats and restrubcture
results_subset <- all_results_fdr %>% right_join(all_results_trait_cats) %>% 
                                  select(exposure.trait, med.trait, outcome.id,
                                         r1.beta_CI, r2.OR_CI, r3.OR_CI, type, med_cat, r1.b, r2.b, r3.b, r1.pval, r2.pval, r3.pval, r1.qval, r3.qval, 
                                         exposure.id, med.id,med_cat) %>% 
                                  distinct()

length(unique(results_subset$exposure.id)) # 103 this is ok 96 will do too

# adhoc save of intesting results - dieatary - PA - BC
tmp <- results_subset %>% filter(exposure_cat == "Diet and supplements", med_cat == "Physical activity" ) 
write_csv(tmp, "01_MR_related/results/mr_evidence_outputs/dietary_PA_BC_mediation.csv")



dim(results_subset) # 9961

counts_raw <- results_subset %>%
  select(exposure.trait,exposure.id, med.id) %>% distinct() %>% count(exposure.trait,exposure.id) %>% rename(med_count_preval=n) 
mean(counts_raw$med_count_preval) # 75.32039 prevalidation mean meds per traits

#write_csv(results_subset, "01_MR_related/mr_evidence_outputs/conf_med_extracted.csv") # p<10e4
write_csv(results_subset, "01_MR_related/results/mr_evidence_outputs/med_extracted_all_r3V3.csv") # qval<0.05 


protein_names <- read_csv("01_MR_related/results/mr_evidence_outputs/protein_names_w_ids.csv", col_names = T) %>% select(name=exposure, gene) %>% distinct() %>% drop_na()
results_subset<- read_csv("01_MR_related/results/mr_evidence_outputs/med_extracted_all_r3V3.csv") %>% 
  filter(type %in% c('mediator')) %>% # i.e. not conf
  create_exposure_categories() %>%
  left_join(protein_names, by = c("med.trait" = 'name')) %>% rename(med.gene = gene) %>% 
  left_join(protein_names, by = c("exposure.trait" = 'name')) %>% rename(exp.gene = gene) %>% 
  select(exposure_cat,exposure.trait, exposure.id, exp.gene, med.trait, med.gene, med.id, med_cat, everything()) 
dim(results_subset)# 3435

results_subset %>% filter(med_cat == "Proteins") %>% filter(is.na(med.gene)) %>% select(med.id, med.trait, med.gene) %>% distinct() %>% View()

length(unique(results_subset$exposure.id)) # 103
# next need to validate E->M for mediators (and M->E for confounders)

#
#
#
#### -- validation step happens here

#  in 06sub_mreve_mediators_validation.R





#### for exporting supplementary

cols<- colnames(results_subset)

cols<- gsub("r1", "step1", cols)
cols<- gsub("r3", "step2", cols)
cols<- gsub("r2", "total", cols)
 
colnames(results_subset)<- cols

# rearrage cols:

results_subset <- results_subset %>% 
  mutate(outcome = case_when(outcome.id == "ieu-a-1126"~ 'Brest cancer overall',
                             outcome.id == "ieu-a-1127"~ 'Brest cancer ER+',
                              outcome.id == "ieu-a-1128"~ 'Brest cancer ER-')) %>% 
  select(exposure_cat:outcome.id,outcome, step1.beta_CI, step2.OR_CI, total.OR_CI, step1.b, step2.b, total.b, step1.pval, step2.pval,total.pval, step1.qval, step2.qval)


write_csv(results_subset, "01_MR_related/results/mr_evidence_outputs/med_extracted_all_r3V3_as_mreveraw_supl.csv") # supl data 5
