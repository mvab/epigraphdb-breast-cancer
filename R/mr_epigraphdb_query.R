library(epigraphdb)
source("helper_functions.R")

bc_gwas <- c( 'ieu-a-1126', 'ieu-a-1127', 'ieu-a-1128', 'ieu-a-1129', 'ieu-a-1130', 'ieu-a-1131', 'ieu-a-1132',
              'ieu-a-1133', 'ieu-a-1134', 'ieu-a-1135', 'ieu-a-1136', 'ieu-a-1137', 
              #'ieu-a-1160', 'ieu-a-1161', 'ieu-a-1162', # iCOGS 2015 weird
               'ieu-a-1163', 'ieu-a-1164', 'ieu-a-1165', 'ieu-a-1166', 'ieu-a-1167', 'ieu-a-1168',
              'ukb-a-55', 'ukb-b-16890', 'ukb-d-C3_BREAST_3')


# MR
# get traits that have pval < 1e-05
query = 
  paste0("
    MATCH (exposure:Gwas)-[mr:MR_EVE_MR]->(outcome:Gwas)
    WHERE outcome.id in ['", paste0(bc_gwas, collapse = "', '"),"'] 
    AND  not exposure.id  in ['", paste0(bc_gwas, collapse = "', '"),"']
    AND (not (toLower(exposure.trait) contains 'breast')) 
    AND mr.pval < 1e-05
    with mr, exposure, outcome
    ORDER BY mr.pval 
    RETURN exposure.id, exposure.trait, exposure.sample_size,
            collect(outcome.id) as outcome_ids, 
            collect(mr.pval) as MR_pvals, collect(mr.b) as MR_beta
    ")
    

out<-query_epigraphdb_as_table(query)
dim(out)# 746


# now those traits that appered at least once in something at <1e05, 
# get MR results for those traits with all BC datasets

query = paste0("
      MATCH (exposure:Gwas)-[mr:MR_EVE_MR]->(outcome:Gwas)
      WHERE outcome.id in ['", paste0(bc_gwas, collapse = "', '"),"'] 
      AND  exposure.id  in ['", paste0(out$exposure.id, collapse = "', '"),"'] 
      with mr, exposure, outcome
      ORDER BY exposure.trait
      RETURN exposure.id, exposure.trait, exposure.sample_size, exposure.sex, exposure.note,
      toInteger(exposure.year) as year, exposure.author as author, exposure.consortium as consortium,
              outcome.id, outcome.sample_size, toInteger(outcome.ncase) as N_case, outcome.year, outcome.nsnp,
              mr.pval, mr.b, mr.se, mr.method, mr.moescore
      ")
        
out2<-query_epigraphdb_as_table(query)
dim(out2)#16277

write_tsv(out2, "explore_MR-EvE_app/data_copy/bc_all_mr_madewR.tsv")

     

### ALTRENATIVE

query = 
  paste0("
    MATCH (exposure:Gwas)-[mr:MR_EVE_MR]->(outcome:Gwas)
    WHERE outcome.id in ['", paste0(bc_gwas, collapse = "', '"),"'] 
    AND  not exposure.id  in ['", paste0(bc_gwas, collapse = "', '"),"']
    AND (not (toLower(exposure.trait) contains 'breast')) 
    AND mr.pval < 1
    with mr, exposure, outcome
    ORDER BY mr.pval 
    RETURN exposure.id, exposure.trait, exposure.sample_size, exposure.sex, exposure.note,
          toInteger(exposure.year) as year, exposure.author as author, exposure.consortium as consortium,
              outcome.id, outcome.sample_size, toInteger(outcome.ncase) as N_case, outcome.year, outcome.nsnp,
              mr.pval, mr.b, mr.se, mr.method, mr.moescore
    ") 


full_results<-query_epigraphdb_as_table(query)
dim(full_results)# 49009 # 45702 w/o 2015

full_results<- full_results %>%
  mutate(loci = mr.b - 1.96 * mr.se, 
                upci = mr.b + 1.96 * mr.se,
                or = exp(mr.b), 
                or_loci = exp(loci), 
                or_upci = exp(upci),
                OR_CI = paste0(round(or,2), " [",round(or_loci,2) ,":",round(or_upci,2), "]")) %>% 
  mutate(effect_direction = ifelse(or_loci > 1 & or_upci >= 1, 'positive',
                                   ifelse(or_loci < 1 & or_upci <= 1, 'negative', 'overlaps null'))) %>% 
  mutate(`MR method and score` = paste0(mr.method," / ", mr.moescore)) 

full_results %>% count(effect_direction)   

# negative          4046
# overlaps null    41018
# positive          3945

sub_results <- full_results %>% filter(effect_direction != 'overlaps null') 
length(unique(sub_results$exposure.id)) #2018 # 1970

query = paste0("
      MATCH (exposure:Gwas)-[mr:MR_EVE_MR]->(outcome:Gwas)
      WHERE outcome.id in ['", paste0(bc_gwas, collapse = "', '"),"'] 
      AND  exposure.id  in ['", paste0(sub_results$exposure.id, collapse = "', '"),"'] 
      with mr, exposure, outcome
      ORDER BY exposure.trait
      RETURN exposure.id, exposure.trait, exposure.sample_size, exposure.sex, exposure.note,
      toInteger(exposure.year) as year, exposure.author as author, exposure.consortium as consortium,
              outcome.id, outcome.sample_size, toInteger(outcome.ncase) as N_case, outcome.year, outcome.nsnp,
              mr.pval, mr.b, mr.se, mr.method, mr.moescore
      ")

out3<-query_epigraphdb_as_table(query)
dim(out3) #44631 #40475


write_tsv(out3, "explore_MR-EvE_app/data_copy/bc_all_mr_fromCIs.tsv")


badgwas <- read_tsv("external_files/gwas_that_may_have_wrong_instr.txt", col_names = F)






MATCH (exposure:Gwas)-[mr_eo:MR_EVE_MR]->(outcome:Gwas)
WHERE outcome.id = 'ieu-a-1128'
AND   exposure.id  = 'met-a-362'
with mr_eo, exposure, outcome
match (exposure:Gwas)-[mr_em:MR_EVE_MR]->(mediator:Gwas)
where mr_em.pval < 1e-05
with mr_eo, exposure, outcome, mediator, mr_em
match (mediator:Gwas)-[mr_mo:MR_EVE_MR]->(outcome:Gwas)
where mr_mo.pval < 1e-05
RETURN mediator.trait, mediator.id,
mr_em.b, mr_em.pval, mr_em.method, mr_em.moescore,
mr_mo.b, mr_mo.pval, mr_mo.method, mr_mo.moescore
