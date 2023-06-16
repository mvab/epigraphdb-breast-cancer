library(tidyverse)
library(TwoSampleMR)
source("helper_functions.R")
source("01_MR_related/scripts/mr_related_functions.R")
source("01_MR_related/scripts/app1_MR-EvE_app/functions.R")

results_subset<- read_csv("01_MR_related/results/mr_evidence_outputs/med_extracted_all_r3V3.csv") %>% 
  filter(type %in% c('mediator')) # conf coud be here too if it was extracted


ids_to_ignore <- c("prot-a-550", "prot-a-3000") # their SNPs are not available in outcomes
# protein location data, for cis instruments extraction
protein_regions<- read_csv("01_MR_related/results/mr_evidence_outputs/protein_gene_regions_ids.csv") %>% 
  filter(chr != "X") %>%  # non 1-22 will be labeles as 'not a protein'
  filter(!exposure.id %in% ids_to_ignore) # these protein's cis insruments are not available in BC output; so treating them not as proteins, we get their cis+trans SNPs


ao <- available_outcomes()
trait_ss <-ao %>%
          filter(id %in% unique(results_subset$med.id, results_subset$exposure.id)) %>% 
          mutate(sample_size = ifelse(is.na(sample_size) & author == "Neale lab", 361194, sample_size)) %>%
          select(exposure.id =id, trait, sex, exposure.sample_size =sample_size, population)
exclude_meds<- trait_ss %>% filter(population != "European" | sex == "Males")

### select pairs for validation

x<- results_subset %>%
  filter(type == 'mediator') %>% 
  create_exposure_categories() %>% 
  filter(!grepl("raw", med.id)) %>% 
  filter(!med.id %in% exclude_meds$exposure.id) %>% 
  select(exposure.trait, exposure.id, exposure_cat, med.trait, med.id, med_cat, type) %>% distinct() # there are dups because of multiple BC 


# need to do 2 types of MR:

# step1:
#   normal Exp -> Med
#   protein Exp -> Med
# step2:
#   normal Med -> Out
#   protein Med -> Out

# step 1
res_step1<-tibble()
completed = 0
for (trait_exp in unique(x$exposure.id) ){ # 103
  
  # collected all outcomes to test
  outcomes_to_test <- x %>% filter(exposure.id == !!trait_exp) %>% pull(med.id)
  
  print(paste0(trait_exp, " has ", length(outcomes_to_test) , " mediators"))
  
  # pre-extract instruments for exposure
  instruments_selected <-  instrument_selection(trait_exp, protein_regions )
  
  res_list_step1_oneexposure <- list()
  
  for (i in 1:length(outcomes_to_test )){
    print(paste0(i, " / ", length(outcomes_to_test)))
          
    if (!outcomes_to_test[i] %in% c("ieu-a-90", "ieu-a-93")){ # outcomes that fail
      
      res_list_step1_oneexposure[[i]]<- do_MR_not_BC(trait_exp_instruments = instruments_selected, 
                                          trait_out = outcomes_to_test[i], 
                                          exposure_ss_df = trait_ss)
    }
  }
  res_step1_oneexposure <- bind_rows(res_list_step1_oneexposure) %>% distinct()
  completed = completed +1
  res_step1 <- bind_rows(res_step1, res_step1_oneexposure)
  
}
res_step1 <- res_step1 %>% 
  # add FDR corrected p-val
  arrange(mr.pval) %>% mutate(qval = p.adjust(mr.pval, method = "BH")) 

write_tsv(res_step1,      "01_MR_related/results/mr_evidence_outputs/redone_MRmeds_fulloutput_fulltable_step1_V3.tsv")


# step 2
# using do_MR <- function(trait, bc_type, exposure_ss_df, protein_regions)
med_ids <- x %>% pull(med.id) %>% unique()

# test 
#a <- do_MR(trait = med_ids[2], bc_type = "all", exposure_ss_df =trait_ss,  protein_regions = protein_regions)

res_df_all <- bind_rows(lapply(med_ids, do_MR, 'all', trait_ss, protein_regions)) 
res_df_pos <- bind_rows(lapply(med_ids, do_MR, 'ER+', trait_ss, protein_regions)) 
res_df_neg <- bind_rows(lapply(med_ids, do_MR, 'ER-', trait_ss, protein_regions))

res_step2 <- bind_rows( res_df_all,res_df_pos,res_df_neg) %>% 
  # add FDR corrected p-val
  arrange(mr.pval) %>% mutate(qval = p.adjust(mr.pval, method = "BH")) %>% 
  arrange(outcome.id)

write_tsv(res_step2,      "01_MR_related/results/mr_evidence_outputs/redone_MRmeds_fulloutput_fulltable_step2_V3.tsv")




## split into MR and sens and save as is

# step 1
res<-read_tsv("01_MR_related/results/mr_evidence_outputs/redone_MRmeds_fulloutput_fulltable_step1_V3.tsv")
redone_MR <- res %>% 
  select(-starts_with('egger_intercept'), -starts_with('Q')) # TODO check which really to keep
redone_MR_sens <- res %>% 
  select(any_of(c('id.exposure', 'exposure', 'outcome', 'method', 'nsnp',
                  "egger_intercept", "egger_intercept_se","egger_intercept_pval"  , 
                  "Q" , "Q_df" ,"Q_pval",
                  'Fst', 'total_r2', 'correct_causal_direction', 'steiger_pval',  "snp_r2.exposure" ,"snp_r2.outcome" ))) %>% 
  filter(method != 'Weighted median')

write_tsv(redone_MR,      "01_MR_related/results/mr_evidence_outputs/redone_MRmeds_fulloutput_step1_V3.tsv")
write_tsv(redone_MR_sens, "01_MR_related/results/mr_evidence_outputs/redone_MRmeds_fulloutput_sens_step1_V3.tsv") 


# step 2
res<-read_tsv("01_MR_related/results/mr_evidence_outputs/redone_MRmeds_fulloutput_fulltable_step2_V3.tsv")
redone_MR <- res %>% 
  select(-starts_with('egger_intercept'), -starts_with('Q'))
redone_MR_sens <- res %>% 
  select(any_of(c('id.exposure', 'exposure', 'outcome', 'method', 'nsnp',
                  "egger_intercept", "egger_intercept_se","egger_intercept_pval"  , 
                  "Q" , "Q_df" ,"Q_pval",
                  'Fst', 'total_r2', 'correct_causal_direction', 'steiger_pval',  "snp_r2.exposure" ,"snp_r2.outcome" ))) %>% 
  filter(method != 'Weighted median')

write_tsv(redone_MR,      "01_MR_related/results/mr_evidence_outputs/redone_MRmeds_fulloutput_step2_V3.tsv")
write_tsv(redone_MR_sens, "01_MR_related/results/mr_evidence_outputs/redone_MRmeds_fulloutput_sens_step2_V3.tsv") 



## now need to find traits that are plausibe mediators when filterign by qvalue in both steps

step1_meds <-read_tsv("01_MR_related/results/mr_evidence_outputs/redone_MRmeds_fulloutput_fulltable_step1_V3.tsv") %>% 
  filter(method %in% c('Inverse variance weighted', 'Wald ratio')) %>% 
  filter(qval< 0.05) %>% pull(outcome.id) %>% unique()

step2_meds <-read_tsv("01_MR_related/results/mr_evidence_outputs/redone_MRmeds_fulloutput_fulltable_step2_V3.tsv") %>% 
  filter(method %in% c('Inverse variance weighted', 'Wald ratio')) %>% 
  filter(qval< 0.05) %>% pull(exposure.id) %>% unique()

plausible_mediators<- intersect(step1_meds, step2_meds)

# extract validated only
step1_validated <- read_tsv("01_MR_related/results/mr_evidence_outputs/redone_MRmeds_fulloutput_fulltable_step1_V3.tsv") %>% 
              select(exposure, id.exposure , id.outcome, beta_CI, pval,qval, effect_direction , nsnp, method) %>% 
              filter(method %in% c('Inverse variance weighted', 'Wald ratio')) %>% 
              filter(outcome.id %in% plausible_mediators)

step2_validated <- read_tsv("01_MR_related/results/mr_evidence_outputs/redone_MRmeds_fulloutput_fulltable_step2_V3.tsv") %>% 
             select(exposure, id.exposure , id.outcome, beta_CI, pval,qval, effect_direction , nsnp, method) %>% 
              filter(method %in% c('Inverse variance weighted', 'Wald ratio')) %>% 
              filter(exposure.id %in% plausible_mediators)

write_tsv(step1_validated, "01_MR_related/results/mr_evidence_outputs/redone_MRmeds_subsetoutput_ivw_step1_validated_V3.tsv")
write_tsv(step2_validated, "01_MR_related/results/mr_evidence_outputs/redone_MRmeds_subsetoutput_ivw_step2_validated_V3.tsv")


