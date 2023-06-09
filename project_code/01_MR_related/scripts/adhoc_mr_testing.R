library(TwoSampleMR)
library(tidyverse)
instruments <- extract_instruments("prot-b-20")
instruments<-read_tsv("/Users/ny19205/OneDrive - University of Bristol/Documents - OneDrive/Mini-project4/01_Data/GWAS_tophits/hdl_cholesterol_tophits.tsv")

out <- extract_outcome_data(
  snps = instruments$SNP,
  outcome = "ukb-b-17400")

#out in local data 
#out <- vroom("/Users/ny19205/OneDrive - University of Bristol/Documents - OneDrive/Mini-project2/01_Data/GWAS_results_tidy/igf_GWAS_tidy_outcome.txt.gz") %>% 
#  filter(SNP %in%instruments$SNP)

harmonised<- harmonise_data(exposure_dat = instruments, 
                            outcome_dat = out)


res_early <- TwoSampleMR::mr(harmonised, method_list = c('mr_ivw')) %>% 
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
  
  
instruments <- extract_instruments('prot-b-20')
instruments<-read_tsv("/Users/ny19205/OneDrive - University of Bristol/Documents - OneDrive/Mini-project2/01_Data/GWAS_tophits/igf_tophits.tsv")

out <- extract_outcome_data(
  snps = instruments$SNP,
  outcome = 'ieu-a-1126')

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

instruments1 <- extract_instruments('prot-b-20')
instruments2 <-read_tsv("/Users/ny19205/OneDrive - University of Bristol/Documents - OneDrive/Mini-project2/01_Data/GWAS_tophits/igf_tophits.tsv")

exposure_list <- list(instruments1, instruments2)

gwas1 <- extract_outcome_data(snps = exposure_list %>% purrr::reduce(bind_rows) %>% pull(SNP), 
                              outcomes = 'prot-b-20')

gwas2 <- vroom::vroom("/Users/ny19205/OneDrive - University of Bristol/Documents - OneDrive/Mini-project2/01_Data/GWAS_results_tidy/igf_GWAS_tidy_outcome.txt.gz")


full_gwas_list <- list(gwas1, gwas2)



# MVMR: 2 mrbase
instruments1 <- extract_instruments('ukb-d-30770_irnt')
instruments2 <- extract_instruments('prot-b-20')
exposure_list <- list(instruments1, instruments2)

gwas1 <- extract_outcome_data(snps = exposure_list %>% purrr::reduce(bind_rows) %>% pull(SNP), 
                              outcomes = 'ukb-d-30770_irnt')
gwas2 <- extract_outcome_data(snps = exposure_list %>% purrr::reduce(bind_rows) %>% pull(SNP), 
                              outcomes = 'prot-b-20')

full_gwas_list <- list(gwas1, gwas2)



source("/Users/ny19205/OneDrive - University of Bristol/Documents - OneDrive/Mini-project2/early-bmi-breast-cancer-mr/functions_mvmr.R")
# create exposure_dat format
exposure_dat <- get_mv_exposures(exposure_list, full_gwas_list, clump_exposures = T) 

#Next, also extract those SNPs from the outcome.
outcome_dat <- extract_outcome_data(exposure_dat$SNP, 'ieu-a-1126')

#Once the data has been obtained, harmonise so that all are on the same reference allele.
mvdat <- mv_harmonise_data(exposure_dat, outcome_dat)

#Finally, perform the multivariable MR analysis
res <- mv_multiple(mvdat)

mv_res<- res$result %>%
  split_outcome() %>% 
  separate(outcome, "outcome", sep="[(]") %>% 
  generate_odds_ratios() %>% 
  select(-id.exposure, -id.outcome)



res1<-mv_res


### mvmr sensitivity

library(MVMR)


mvmr_input <- make_mvmr_input(exposure_dat, outcome.id.mrbase='ieu-a-1126')

# format data to be in MVMR package-compatiable df
mvmr_out <- format_mvmr(BXGs = mvmr_input$XGs %>% select(contains("beta")),  # exposure betas
                        BYG = mvmr_input$YG$beta.outcome,                        # outcome beta
                        seBXGs = mvmr_input$XGs %>% select(contains("se")),      # exposure SEs
                        seBYG = mvmr_input$YG$se.outcome,                        # outcome SEs
                        RSID = mvmr_input$XGs$SNP)                               # SNPs

#  estimate causal effects using method in MVMR package
mvmr_res <-ivw_mvmr(r_input=mvmr_out) %>% 
  tidy_mvmr_output() %>% 
  mutate(exposure = mvmr_input$exposures,
         outcome = 'ieu-a-1126')


gen_cov <- 0


#Test for weak instruments
sres <- strength_mvmr(r_input=mvmr_out, gencov=gen_cov)
colnames(sres) = paste(c("Childhood BMI", mediator_name), "(Fst)")
print(sres)

#Test for horizontal pleiotropy
pres <- pleiotropy_mvmr(r_input=mvmr_out, gencov=gen_cov)

mvmr_sens_df <- sres
mvmr_sens_df$Qstat <- pres$Qstat
mvmr_sens_df$Qpval <- pres$Qpval



### other test
### 
### 

exp <-read_tsv("/Users/ny19205/OneDrive - University of Bristol/Documents - OneDrive/Mini-project4/01_Data/GWAS_tophits/hdl_cholesterol_tophits.tsv")
out <- vroom::vroom("/Users/ny19205/OneDrive - University of Bristol/Documents - OneDrive/Mini-project4/01_Data/GWAS_tidy/snoring_GWAS_tidy_outcome.txt.gz") %>% 
  filter(SNP %in%exp$SNP)


harmonised<- harmonise_data(exposure_dat = exp, 
                            outcome_dat = out)


res_early <- TwoSampleMR::mr(harmonised, method_list = c('mr_ivw')) %>% 
  split_outcome() %>% 
  split_exposure() %>% 
  generate_odds_ratios()
