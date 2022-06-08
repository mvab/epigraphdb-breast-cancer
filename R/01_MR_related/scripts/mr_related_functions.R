
#### Method for finding interesting things

traits_for_follow_up <- function(input) {
  # this is a method to select traits that  have effect in 2/3 of main dataset in a subtype
  
  
  ids <- input %>% pull(exposure.id) %>% unique()
  res <- list()
  
  for (i in ids){
    # pre-set of results to use for selection
    input_subset_array <- input %>% 
      filter(chip %in% c('Meta', "OncArray",  "iCOG2017")) %>% 
      select(exposure.id, exposure.trait, exposure_cat, outcome.id, outcome, mr.pval, starts_with('or'),
             effect_direction, mr.method, mr.moescore) %>%
      filter(exposure.id == i) 
    
    if (dim(input_subset_array)[1] != 0) {
      # continue if these outcomes (arrays) are present
      
      # drop results that 'overlap null' for all selected outcomes
      effect <-input_subset_array %>% pull(effect_direction) %>% unique()
      
      if (length(effect) == 1 & "overlaps null" %in% effect){
        # Skipping this trait
        i_name <- input %>% filter(exposure.id == i) %>% pull(exposure.trait) %>% unique()
        print(paste0(i," / ", i_name))
        
      } else {
        #get counts per BC type
        counts_per_type <- 
          input_subset_array %>% 
          group_by(outcome) %>% 
          janitor::tabyl(outcome, effect_direction) %>% 
          as_tibble()
        
        
        
        # identify cases when traits have >=2 effects in positive or negative direction
        if ('positive' %in% colnames(counts_per_type) & 'negative' %in% colnames(counts_per_type)){
          dir <- 'both'
          tmp <- bind_rows(
            counts_per_type %>% filter(positive >= 2),
            counts_per_type %>% filter(negative >= 2)) %>% 
            distinct()
          
          if (sum(tmp$negative) == 0){dir <- 'positive'
          } else if (sum(tmp$positive) == 0){dir <- 'negative'
          } else if (sum(tmp$negative) > sum(tmp$positive) ){dir <- 'negative'
          } else if (sum(tmp$negative) < sum(tmp$positive)){dir <- 'positive'}
          
          
        } else if ('positive' %in% colnames(counts_per_type)){
          tmp <- counts_per_type %>% filter(positive >=  2)
          dir <- 'positive'
          
        } else if ('negative' %in% colnames(counts_per_type)){
          tmp <-  counts_per_type %>% filter(negative >=2 )
          dir <- 'negative'
        }
        
        if (dim(tmp)[1] != 0) {
          # collect results for trait in a vector and keep in a list
          trait_tmp<-c(unique(input_subset_array$exposure_cat), unique(input_subset_array$exposure.id),  unique(input_subset_array$exposure.trait), 0,0,0, dir)
          names(trait_tmp)<-c('exposure_cat', 'id', 'trait', 'all', "ER+", "ER-", 'dir')
          
          if ("Breast cancer (all)" %in% tmp$outcome){ trait_tmp[["all"]] <- 1 }
          if ("ER+ postmeno"        %in% tmp$outcome){ trait_tmp[["ER+"]] <- 1 }
          if ("ER- premeno"         %in% tmp$outcome){ trait_tmp[["ER-"]] <- 1 }
          res[[unique(input_subset_array$exposure.id)]] <- trait_tmp
        }
        
      }
      if (exists("tmp")){ rm(tmp) }
      
    }
  }
  
  
  #list to df
  out<- bind_rows(res) 
  return(out)
}



### validation MR
do_MR <- function(trait, bc_type){
  
  
  if       (bc_type == 'all'){bc.id <- 'ieu-a-1126'
  }else if (bc_type == 'ER+'){bc.id <- 'ieu-a-1127'
  }else if (bc_type == 'ER-'){bc.id <- 'ieu-a-1128'}
  
  instruments <- extract_instruments(trait)
  out <- extract_outcome_data(snps = instruments$SNP,
                              outcome = bc.id)
  harmonised<- harmonise_data(exposure_dat = instruments, 
                              outcome_dat = out)
  #mr
  res <- TwoSampleMR::mr(harmonised, method_list=c('mr_ivw','mr_wald_ratio','mr_egger_regression','mr_weighted_median')) %>% 
    split_outcome() %>% 
    split_exposure() %>% 
    generate_odds_ratios() %>% 
    mutate(OR_CI = paste0(round(or,2), " [",round(or_lci95,2) ,":",round(or_uci95,2), "]")) %>% 
    mutate(effect_direction = ifelse(or_lci95 > 1 & or_uci95 >= 1, 'positive',
                                     ifelse(or_lci95 < 1 & or_uci95 <= 1, 'negative', 'overlaps null'))) 
  # sensitivity
  if (dim(harmonised)[1]>1){
    
    if (dim(mr_pleiotropy_test(harmonised))[1] ==0){
      res_sens <- res %>% select(id.exposure, id.outcome, exposure, outcome, nsnp) %>% distinct() %>% mutate(method = "NO SENSITIVITY")
    } else{
    
    res_sens <-
      full_join(mr_pleiotropy_test(harmonised),
                mr_heterogeneity(harmonised, method_list=c("mr_egger_regression", "mr_ivw"))) %>% 
      split_outcome() %>%
      split_exposure() %>% 
      rename(egger_intercept_pval = pval,
             egger_intercept_se = se)
    }
  } else {
    # making dummy sens anlysis table as a placeholder
    res_sens <- res %>% select(id.exposure, id.outcome, exposure, outcome, nsnp) %>% distinct() %>% mutate(method = "NO SENSITIVITY")
  }
  
  # join mr and sens in one table
  out <- full_join(res, res_sens)
  
  return( out) 
} 



# add sensitivity analysis cols 

add_sensitivity_cols <- function(main_results, sens_results, outcome_name, mtc){
  
  
  # main results
  all <- main_results %>% mutate(low_SNP_num = ifelse(nsnp < 5 , T, F))
  type <- all %>% filter(id.outcome == outcome_name)
  type_filtered<- type  %>% 
    filter(method %in%  c('Inverse variance weighted', "Wald ratio")  & effect_direction != "overlaps null") 
  print(paste( outcome_name, ":::", length(unique(type_filtered$id.exposure)), "/",  length(unique(type$id.exposure)) ) )
  
  ####type <- type %>% filter(id.exposure %in% type_filtered$id.exposure) 
  
  ##count alt MR success
  effect_count <- type %>% 
    filter(effect_direction != "overlaps null") %>% 
    group_by( exposure, id.exposure) %>% count( exposure, id.exposure) %>% 
    rename(num_mr_effect = n) %>% ungroup()
  
  type <- type %>%
    left_join(effect_count) %>%
    filter(method %in%  c('Inverse variance weighted', "Wald ratio"))
  
  # sens results 
  sens_type <- sens_results %>% 
    filter(id.outcome == outcome_name) %>% 
    filter(method ==  'Inverse variance weighted')
  
  sens_type <- sens_type %>%   
    mutate(`egger_intercept_less_than_0.05` = ifelse(abs(egger_intercept) < 0.05, T, F) ) %>%  # pleio 1?
    mutate(`egger_intercept_pval_less_than_0.05` = ifelse(egger_intercept_pval < 0.05, T, F) ) %>%  # pleio 2?
    mutate(`heterogeneity_Q_pval_less_than_0.05` = ifelse(Q_pval < 0.05, T, F) ) %>%  # heterogeniety 
    rename("heterogeneity_Q"="Q",
           "heterogeneity_Q_df"="Q_df",
           "heterogeneity_Q_pval"="Q_pval")
  
  
  # merge into one df
  merged<- left_join(type, sens_type) %>% left_join(mtc) %>% 
    select(everything(), -low_SNP_num, low_SNP_num, -num_mr_effect, -pval_FDR_adj,  num_mr_effect, pval_FDR_adj)
  
}



### QUICK MR

quick_mr <- function(exp, out){
  
  instruments <- extract_instruments(exp)
  outcome <- extract_outcome_data(snps = instruments$SNP,
                                  outcome = out)
  harmonised<- harmonise_data(exposure_dat = instruments, 
                              outcome_dat = outcome)
  res <- TwoSampleMR::mr(harmonised, method_list=c('mr_ivw','mr_wald_ratio')) %>% 
    split_outcome() %>% 
    split_exposure() %>% 
    generate_odds_ratios() %>% 
    mutate(beta_CI = paste0(round(b,2), " [",round(lo_ci,2) ,":",round(up_ci,2), "]")) %>% 
    mutate(OR_CI = paste0(round(or,2), " [",round(or_lci95,2) ,":",round(or_uci95,2), "]")) %>% 
    mutate(effect_direction = ifelse(or_lci95 > 1 & or_uci95 >= 1, 'positive',
                                     ifelse(or_lci95 < 1 & or_uci95 <= 1, 'negative', 'overlaps null'))) 
  
  res
  
}

quick_mvmr <- function(exp1, exp2, out){
  
  source("/Users/ny19205/OneDrive - University of Bristol/Documents - OneDrive/Mini-project2/early-bmi-breast-cancer-mr/functions_mvmr.R")
  
  
  # MVMR: 2 mrbase
  instruments1 <- extract_instruments(exp1)
  instruments2 <- extract_instruments(exp2)
  exposure_list <- list(instruments1, instruments2)
  
  gwas1 <- extract_outcome_data(snps = exposure_list %>% purrr::reduce(bind_rows) %>% pull(SNP), 
                                outcomes = exp1)
  gwas2 <- extract_outcome_data(snps = exposure_list %>% purrr::reduce(bind_rows) %>% pull(SNP), 
                                outcomes = exp2)
  
  full_gwas_list <- list(gwas1, gwas2)
  
  # create exposure_dat format
  exposure_dat <- get_mv_exposures(exposure_list, full_gwas_list, clump_exposures = F) 
  
  #Next, also extract those SNPs from the outcome.
  outcome_dat <- extract_outcome_data(exposure_dat$SNP, out)
  
  #Once the data has been obtained, harmonise so that all are on the same reference allele.
  mvdat <- mv_harmonise_data(exposure_dat, outcome_dat)
  
  #Finally, perform the multivariable MR analysis
  res <- mv_multiple(mvdat)
  
  mv_res<- res$result %>%
    split_outcome() %>% 
    split_exposure() %>% 
    separate(outcome, "outcome", sep="[(]") %>% 
    generate_odds_ratios() %>% 
    select(-id.exposure, -id.outcome) %>% 
    mutate(OR_CI = paste0(round(or,3), " [",round(or_lci95,3) ,":",round(or_uci95,3), "]")) %>% 
    mutate(effect_direction = ifelse(or_lci95 > 1 & or_uci95 >= 1, 'positive',
                                     ifelse(or_lci95 < 1 & or_uci95 <= 1, 'negative', 'overlaps null'))) 
  
  mv_res
}



mvmr_mixed_sources <- function(id1, outcome.id, id2_tophits_file, id2_gwas_file ){
  
  path <- "/Users/ny19205/OneDrive - University of Bristol/Documents - OneDrive/Mini-project2/"
  source(paste0(path, "early-bmi-breast-cancer-mr/functions_mvmr.R"))
  
  # MVMR
  
  instruments1 <- extract_instruments(id1)
  instruments2 <-id2_tophits_file
  
  exposure_list <- list(instruments1, instruments2)
  
  gwas1 <- extract_outcome_data(snps = exposure_list %>% purrr::reduce(bind_rows) %>% pull(SNP), 
                                outcomes = id1)
  
  gwas2 <- id2_gwas_file
  
  
  full_gwas_list <- list(gwas1, gwas2)
  
  
  # create exposure_dat format
  exposure_dat <- get_mv_exposures(exposure_list, full_gwas_list, clump_exposures = F) 
  
  #Next, also extract those SNPs from the outcome.
  outcome_dat <- extract_outcome_data(exposure_dat$SNP, outcome.id)
  
  #Once the data has been obtained, harmonise so that all are on the same reference allele.
  mvdat <- mv_harmonise_data(exposure_dat, outcome_dat)
  
  #Finally, perform the multivariable MR analysis
  res <- mv_multiple(mvdat)
  
  mv_res<- res$result %>%
    split_outcome() %>% 
    split_exposure() %>% 
    separate(outcome, "outcome", sep="[(]") %>% 
    generate_odds_ratios() %>% 
    select(-id.exposure, -id.outcome) %>% 
    mutate(OR_CI = paste0(round(or,2), " [",round(or_lci95,2) ,":",round(or_uci95,2), "]")) %>% 
    mutate(effect_direction = ifelse(or_lci95 > 1 & or_uci95 >= 1, 'positive',
                                     ifelse(or_lci95 < 1 & or_uci95 <= 1, 'negative', 'overlaps null'))) %>% 
    select(exposure, outcome, starts_with("or", ignore.case = T), effect_direction) 
  
  return(mv_res)
  
}



