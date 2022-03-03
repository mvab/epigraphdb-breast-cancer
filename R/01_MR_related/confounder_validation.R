library(tidyverse)
library(TwoSampleMR)
source("helper_functions.R")
source("01_MR_related/explore_MR-EvE_app/functions.R")


results_subset<- read_csv("01_MR_related/mr_evidence_outputs/conf_med_extracted_all.csv") %>% 
  filter(type %in% c('confounder', 'mediator'))


### select pairs for validation

x<- results_subset %>% select(exposure.trait, exposure.id, med.trait, med.id, type) %>% distinct()
dim(x)
x1 <- x %>% filter(type == 'mediator') %>% select(exposure.trait, exposure.id, outcome.trait = med.trait,outcome.id = med.id)
x2 <- x %>% filter(type == 'confounder') %>% select(outcome.trait = exposure.trait, outcome.id = exposure.id, exposure.trait = med.trait,exposure.id = med.id) # reverse
y <- bind_rows(x1, x2) %>% distinct()
dim(y)




do_MR_trait_pairs <- function(exp_trait, out_trait){
  
  instruments <- extract_instruments(exp_trait)
  out <- extract_outcome_data(snps = instruments$SNP,
                              outcome = out_trait)
  harmonised<- harmonise_data(exposure_dat = instruments, 
                              outcome_dat = out)
  #mr
  res <- mr(harmonised, method_list=c('mr_ivw','mr_wald_ratio','mr_egger_regression','mr_weighted_median')) %>% 
    split_outcome() %>% 
    split_exposure() %>% 
    generate_odds_ratios() %>% 
    mutate(OR_CI = paste0(round(or,3), " [",round(or_lci95,3) ,":",round(or_uci95,3), "]")) %>% 
    mutate(effect_direction = ifelse(or_lci95 > 1 & or_uci95 >= 1, 'positive',
                                     ifelse(or_lci95 < 1 & or_uci95 <= 1, 'negative', 'overlaps null'))) 
  # sensitivity
  if (dim(harmonised)[1]>1){
    res_sens <-
      full_join(mr_pleiotropy_test(harmonised),
                mr_heterogeneity(harmonised, method_list=c("mr_egger_regression", "mr_ivw"))) %>% 
      split_outcome() %>%
      split_exposure() %>% 
      rename(egger_intercept_pval = pval,
             egger_intercept_se = se)
  } else {
    # making dummy sens anlysis table as a placeholder
    res_sens <- res %>% select(id.exposure, id.outcome, exposure, outcome, nsnp) %>% distinct() %>% mutate(method = "NO SENSITIVITY")
  }
  
  # join mr and sens in one table
  out <- full_join(res, res_sens)
  
  return( out) 
} 


#res_list<-list()
for (i in 1765:dim(y)[1]){
  print(paste(i, y$exposure.id[i], y$outcome.id[i]))
  
  if (!y$outcome.id[i] %in% c("ieu-a-90", "ieu-a-93")){ # outcomes that fail
  
    res_list[[i]]<- do_MR_trait_pairs( y$exposure.id[i], y$outcome.id[i] )
  }
}


res <- bind_rows(res_list) %>% distinct()



# split into MR and sens and save
redone_MR <- res %>% 
  select(-starts_with('egger_intercept'), -starts_with('Q'))
redone_MR_sens <- res %>% 
  select(id.exposure, id.outcome, exposure, outcome, method, nsnp,
         starts_with('egger_intercept'), starts_with('Q')) %>% 
  filter(method != 'Weighted median')

write_tsv(redone_MR,      "01_MR_related/mr_evidence_outputs/redone_MRconf_fulloutput.tsv")
write_tsv(redone_MR_sens, "01_MR_related/mr_evidence_outputs/redone_MRconf_fulloutput_sens.tsv")



redoneMR_tidy <- redone_MR %>%  
  select(exposure, id.exposure , id.outcome, OR_CI, effect_direction , nsnp, method) %>% 
  filter(effect_direction != "overlaps null") %>% 
  filter(method %in% c('Inverse variance weighted', 'Wald ratio'))

write_tsv(redoneMR_tidy, "01_MR_related/mr_evidence_outputs/redone_MRconf_subsetoutput_ivw.tsv")

