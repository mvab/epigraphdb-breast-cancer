library(TwoSampleMR)
library(tidyverse)
instruments <- extract_instruments('ukb-b-5495')
instruments<-read_tsv("/Users/ny19205/OneDrive - University of Bristol/Documents - OneDrive/Mini-project2/01_Data/GWAS_tophits/igf_tophits.tsv")


out <- extract_outcome_data(
  snps = instruments$SNP,
  outcome = "ieu-a-1128")

harmonised<- harmonise_data(exposure_dat = instruments, 
                            outcome_dat = out)


res_early <- mr(harmonised) %>% 
  split_outcome() %>% 
  split_exposure() %>% 
  generate_odds_ratios()
  
betaine <- res_early %>% select(id.exposure, id.outcome, starts_with("or"))
igf <-  res_early %>% select(id.exposure, id.outcome, starts_with("or"))

instruments <- extract_instruments('ukb-b-4650')


out <- extract_outcome_data(
  snps = instruments$SNP,
  outcome = 'ieu-a-1')

harmonised<- harmonise_data(exposure_dat = instruments, 
                            outcome_dat = out)


res_early <- mr(harmonised, method_list=c('mr_ivw','mr_egger_regression','mr_weighted_median', 'mr_wald_ratio')) %>% 
  split_outcome() %>% 
  split_exposure() %>% 
  generate_odds_ratios()



#-----
  
  
instruments <- extract_instruments('met-a-362')
instruments<-read_tsv("/Users/ny19205/OneDrive - University of Bristol/Documents - OneDrive/Mini-project2/01_Data/GWAS_tophits/igf_tophits.tsv")

out <- extract_outcome_data(
  snps = instruments$SNP,
  outcome = 'met-a-362')

harmonised<- harmonise_data(exposure_dat = instruments, 
                            outcome_dat = out)


res_early <- mr(harmonised, method_list=c('mr_ivw','mr_egger_regression','mr_weighted_median', 'mr_wald_ratio')) %>% 
  split_outcome() %>% 
  split_exposure() %>% 
  generate_odds_ratios()

  
----
  # test in BC subtypes
  
instruments <- extract_instruments('met-a-362')



# MVMR: 1 mrbase, 1 text

instruments1 <- extract_instruments('met-a-362')
instruments2 <-read_tsv("/Users/ny19205/OneDrive - University of Bristol/Documents - OneDrive/Mini-project2/01_Data/GWAS_tophits/igf_tophits.tsv")

exposure_list <- list(instruments1, instruments2)

gwas1 <- extract_outcome_data(snps = exposure_list %>% purrr::reduce(bind_rows) %>% pull(SNP), 
                              outcomes = 'met-a-362')

gwas2 <- vroom::vroom("/Users/ny19205/OneDrive - University of Bristol/Documents - OneDrive/Mini-project2/01_Data/GWAS_results_tidy/igf_GWAS_tidy_outcome.txt.gz")


full_gwas_list <- list(gwas1, gwas2)



# MVMR: 2 mrbase
instruments1 <- extract_instruments('met-a-362')
instruments2 <- extract_instruments('met-a-558')
exposure_list <- list(instruments1, instruments2)

gwas1 <- extract_outcome_data(snps = exposure_list %>% purrr::reduce(bind_rows) %>% pull(SNP), 
                              outcomes = 'met-a-362')
gwas2 <- extract_outcome_data(snps = exposure_list %>% purrr::reduce(bind_rows) %>% pull(SNP), 
                              outcomes = 'met-a-558')

full_gwas_list <- list(gwas1, gwas2)



source("/Users/ny19205/OneDrive - University of Bristol/Documents - OneDrive/Mini-project2/early-bmi-breast-cancer-mr/functions_mvmr.R")
# create exposure_dat format
exposure_dat <- get_mv_exposures(exposure_list, full_gwas_list, clump_exposures = F) 

#Next, also extract those SNPs from the outcome.
outcome_dat <- extract_outcome_data(exposure_dat$SNP, 'ieu-a-1128')

#Once the data has been obtained, harmonise so that all are on the same reference allele.
mvdat <- mv_harmonise_data(exposure_dat, outcome_dat)

#Finally, perform the multivariable MR analysis
res <- mv_multiple(mvdat)

mv_res<- res$result %>%
  split_outcome() %>% 
  separate(outcome, "outcome", sep="[(]") %>% 
  generate_odds_ratios() %>% 
  select(-id.exposure, -id.outcome)

