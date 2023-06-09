library(tidyverse)
library(TwoSampleMR)
source("helper_functions.R")
source("01_MR_related/scripts/mr_related_functions.R")
source("01_MR_related/scripts/app1_MR-EvE_app/functions.R")

results_subset<- read_csv("01_MR_related/results/mr_evidence_outputs/med_extracted_all_r3V3.csv") %>% 
  filter(type %in% c('mediator')) # conf coud be here too if it was extracted


### select pairs for validation

x<- results_subset %>% select(exposure.trait, exposure.id, med.trait, med.id, type) %>% distinct()
dim(x)
x1 <- x %>% filter(type == 'mediator') %>% 
  select(exposure.trait, exposure.id, outcome.trait = med.trait,outcome.id = med.id, type)
x2 <- x %>% filter(type == 'confounder') %>% 
  select(outcome.trait = exposure.trait, outcome.id = exposure.id, exposure.trait = med.trait,exposure.id = med.id, type) # reverse

x2 <- tibble()
y <- bind_rows(x1, x2) %>% distinct() ## as we are only using mediators, x2 will be empty
dim(y)


res_list<-list()
for (i in 1:dim(y)[1]){
  print(paste(i, y$exposure.id[i], y$outcome.id[i]))
  
  if (!y$outcome.id[i] %in% c("ieu-a-90", "ieu-a-93")){ # outcomes that fail
    
    res_list[[i]]<- do_MR_trait_pairs( y$exposure.id[i], y$outcome.id[i] )
  }
}
res <- bind_rows(res_list) %>% distinct()
length(unique(res$id.exposure)) # 103 this is ok


write_tsv(res,      "01_MR_related/results/mr_evidence_outputs/redone_MRmeds_fulloutput_fulltableV3.tsv")

# split into MR and sens and save
redone_MR <- res %>% 
  select(-starts_with('egger_intercept'), -starts_with('Q'))
redone_MR_sens <- res %>% 
  select(id.exposure, id.outcome, exposure, outcome, method, nsnp,
         starts_with('egger_intercept'), starts_with('Q')) %>% 
  filter(method != 'Weighted median')

write_tsv(redone_MR,      "01_MR_related/results/mr_evidence_outputs/redone_MRmeds_fulloutputV3.tsv")
write_tsv(redone_MR_sens, "01_MR_related/results/mr_evidence_outputs/redone_MRmeds_fulloutput_sensV3.tsv") 



redoneMR_tidy <- redone_MR %>%  
  select(exposure, id.exposure , id.outcome, beta_CI, pval, effect_direction , nsnp, method) %>% 
  filter(effect_direction != "overlaps null") %>% 
  filter(method %in% c('Inverse variance weighted', 'Wald ratio'))

write_tsv(redoneMR_tidy, "01_MR_related/results/mr_evidence_outputs/redone_MRmeds_subsetoutput_ivwV3.tsv")

