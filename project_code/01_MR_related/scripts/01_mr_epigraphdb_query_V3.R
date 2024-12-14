library(epigraphdb)
library(dplyr)
library(readr)
source("helper_functions.R")



# outcomes list
bc_gwas <- c( 'ieu-a-1126', 'ieu-a-1127', 'ieu-a-1128') # only using these fro discovery in the new version V3


# query all MR results for the outcomes, not restricting by p-value
query = 
  paste0("
    MATCH (exposure:Gwas)-[mr:MR_EVE_MR]->(outcome:Gwas)
    WHERE outcome.id in ['", paste0(bc_gwas, collapse = "', '"),"'] 
    AND  not exposure.id  in ['", paste0(bc_gwas, collapse = "', '"),"']
    AND (not (toLower(exposure.trait) contains 'breast')) 
    AND mr.pval < 1
    AND exposure.population = 'European'
    with mr, exposure, outcome
    ORDER BY mr.pval 
    RETURN exposure.id, exposure.trait, exposure.sample_size, exposure.sex, exposure.note, exposure.population,
          toInteger(exposure.year) as year, exposure.author as author, exposure.consortium as consortium,
              outcome.id, outcome.sample_size, toInteger(outcome.ncase) as N_case, outcome.year, outcome.nsnp,
              mr.pval, mr.b, mr.se,mr.nsnp,mr.method, mr.moescore
    ") 
full_results<-query_epigraphdb_as_table(query)
length(unique(full_results$exposure.id)) # 2218 -- total number of traits with any result -V3

# calculate CI and get effect direction
full_results<- full_results %>%
  mutate( loci = mr.b - 1.96 * mr.se, 
          upci = mr.b + 1.96 * mr.se,
          or = exp(mr.b), 
          or_loci = exp(loci), 
          or_upci = exp(upci),
          OR_CI = paste0(round(or,2), " [",round(or_loci,2) ,":",round(or_upci,2), "]")) %>% 
  mutate(effect_direction = ifelse(or_loci > 1 & or_upci >= 1, 'positive',
                                   ifelse(or_loci < 1 & or_upci <= 1, 'negative', 'overlaps null'))) %>% 
  mutate(`MR method and score` = paste0(mr.method," / ", mr.moescore)) %>% 
  # add FDR corrected p-val
  arrange(mr.pval) %>% mutate(qval = p.adjust(mr.pval, method = "BH")) 


full_results %>% filter(mr.pval < 0.05) %>% pull(exposure.id) %>% unique() %>% length() # 808
full_results %>% filter(effect_direction != 'overlaps null') %>% pull(exposure.id) %>% unique() %>% length() # 839
full_results %>% filter(qval < 0.05) %>% pull(exposure.id) %>% unique() %>% length() # 378


# save all for supl data
full_results_save<- full_results

write_csv(full_results_save, "01_MR_related/results/mr_evidence_outputs/all_mreve_bc_resultsV3.csv")  # prereq for supl data 1;  V3 - 3 outcomes version 


sub_results <- full_results %>% filter(qval < 0.05) 
length(unique(sub_results$exposure.id)) # 378  in V3

# now re-extract the full MR results (for all outcomes) for those 391 in V3
query = paste0("
      MATCH (exposure:Gwas)-[mr:MR_EVE_MR]->(outcome:Gwas)
      WHERE outcome.id in ['", paste0(bc_gwas, collapse = "', '"),"'] 
      AND  exposure.id  in ['", paste0(sub_results$exposure.id, collapse = "', '"),"'] 
      with mr, exposure, outcome
      ORDER BY exposure.trait
      RETURN exposure.id, exposure.trait, exposure.sample_size, exposure.sex, exposure.note, exposure.population,
      toInteger(exposure.year) as year, exposure.author as author, exposure.consortium as consortium,
              outcome.id, outcome.sample_size, toInteger(outcome.ncase) as N_case, outcome.year, outcome.nsnp,
              mr.pval, mr.b, mr.se, mr.nsnp, mr.method, mr.moescore
      ")

out3<-query_epigraphdb_as_table(query)
dim(out3) 
length(unique(out3$exposure.id)) # 378

write_tsv(out3, "01_MR_related/scripts/app1_MR-EvE_app/data_copy/bc_all_mr_fromCIsV3.tsv")  # main query result (many outcomes)-- saves directly to the app that uses it!
# V3 is based on 3 discovery outcomes
