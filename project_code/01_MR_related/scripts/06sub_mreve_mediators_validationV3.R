library(tidyverse)
library(TwoSampleMR)
source("helper_functions.R")
source("01_MR_related/scripts/mr_related_functions.R")
source("01_MR_related/scripts/app1_MR-EvE_app/functions.R")


# exposures that passed FDR in BCAC 2017: only really for those I should be looking for mediators
total_effect_fdr <- read_tsv(paste0("01_MR_related/results/mr_evidence_outputs/all_multiple_testingV3.tsv")) %>% 
  filter(!is.na(id.outcome)) %>% 
  filter(method %in% c('Inverse variance weighted', 'Wald ratio')) %>% 
  filter(qval < 0.05)
length(unique(total_effect_fdr$id.exposure)) # 50


results_subset<- read_csv("01_MR_related/results/mr_evidence_outputs/med_extracted_all_r3V3.csv") %>% 
  filter(type %in% c('mediator')) %>% # conf coud be here too if it was extracted
  filter(exposure.id %in% unique(total_effect_fdr$id.exposure)) # keep only those exp (and their meds) that passed FDR

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
exclude_meds<- trait_ss %>% filter(population != "European" | sex == "Males" | grepl("raw", exposure.id))

### select pairs for validation

x<- results_subset %>%
  filter(type == 'mediator') %>% 
  create_exposure_categories() %>% 
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
res_step1_list<-list()
completed = 0

exposures_to_do <- unique(x$exposure.id)[1:103] # 103
#exposures_to_do<-c("prot-a-643",  "prot-a-1369", "prot-a-1317", "prot-a-1540", "prot-a-1541", "prot-a-2824" ,"prot-a-1074", "prot-a-2493", "prot-a-2291")
for (trait_exp in exposures_to_do){
  
  # collected all outcomes to test
  outcomes_to_test <- x %>% filter(exposure.id == !!trait_exp) %>% pull(med.id)
  
  print(paste0(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ",
               trait_exp, " has ", length(outcomes_to_test) , " mediators"))
  
  if (length(outcomes_to_test) > 0) {
    
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
    res_step1_list[[completed]] <- res_step1_oneexposure
  } 
}
#save(res_step1_list,file ="~/Desktop/tmp_res_step1_list.Rdata")

res_step1<- bind_rows(res_step1_list)

res_step1ivw <- res_step1 %>%  
  filter(method %in% c('Inverse variance weighted', 'Wald ratio')) 
# now, put data into a list for each exposure, so that can apply FDR correction for each separately
res_step1ivw_list <- list()
for (i in unique(res_step1ivw$id.exposure)){
  res_step1ivw_list[[i]] <- res_step1ivw %>%
    filter(id.exposure == i) %>% 
    # add FDR corrected p-val
    arrange(pval) %>% mutate(qval = p.adjust(pval, method = "BH")) 
}
res_step1ivw <- bind_rows(res_step1ivw_list) #%>% select(-c(SNP:data_source.outcome))


res_step1other <- res_step1 %>%  
  filter(!method %in% c('Inverse variance weighted', 'Wald ratio')) #%>% select(-c(SNP:data_source.outcome))

res_step1_fdr<- bind_rows(res_step1ivw,res_step1other) %>% arrange(id.outcome, id.exposure, method) %>% select(exposure, everything())
length(unique(res_step1_fdr$id.exposure))

write_tsv(res_step1_fdr,      "01_MR_related/results/mr_evidence_outputs/redone_MRmeds_fulloutput_fulltable_step1_V3.tsv")



# step 2
med_ids <- x %>% pull(med.id) %>% unique()
#med_ids <- c( "prot-a-1074",      "prot-a-1317" ,     "prot-a-1369"  ,    "prot-a-1397" ,     "prot-a-1540"  ,    "prot-a-1541" ,    
#              "prot-a-1898" ,     "prot-a-2291"  ,    "prot-a-2493"   ,   "prot-a-2824"  ,    "prot-a-3123"   ,   "prot-a-643"   ,  "ukb-d-30630_irnt")

res_step2_list <- list()
for (i in 1:length(med_ids)){
  print(i)
  trait = med_ids[i]
  res_step2_list[[i]] <- do_MR_all_BC(trait, exposure_ss_df=trait_ss, protein_regions)
}

res_step2 <- bind_rows( res_step2_list) 

res_step2ivw_backup<- res_step2ivw
res_step2ivw <- res_step2 %>%  
  filter(method %in% c('Inverse variance weighted', 'Wald ratio'))

# now, put data into a list for each outcome, so that can apply FDR correction for each separately
res_step2ivw_list <- list()
for (i in unique(res_step2ivw$id.outcome)){
  res_step2ivw_list[[i]] <- res_step2ivw %>%
    filter(id.outcome == i) %>% 
    # add FDR corrected p-val
    arrange(pval) %>% mutate(qval = p.adjust(pval, method = "BH")) 
}
res_step2ivw <- bind_rows(res_step2ivw_list)


res_step2other <- res_step2 %>%  
  filter(!method %in% c('Inverse variance weighted', 'Wald ratio'))

res_step2_fdr<- bind_rows(res_step2ivw,res_step2other) %>% arrange(id.exposure, method,id.outcome ) %>% select(exposure, everything())

write_tsv(res_step2_fdr,      "01_MR_related/results/mr_evidence_outputs/redone_MRmeds_fulloutput_fulltable_step2_V3.tsv")


####3 supl data 6 ---- no this is not 6; 6 is a at the bottom; ok this is 6b wiht all sens tests
#
    res1<-read_tsv("01_MR_related/results/mr_evidence_outputs/redone_MRmeds_fulloutput_fulltable_step1_V3.tsv")
    res2<-read_tsv( "01_MR_related/results/mr_evidence_outputs/redone_MRmeds_fulloutput_fulltable_step2_V3.tsv")
    twostep<-list()
    twostep[["step1"]] <- res1
    twostep[["step2"]] <- res2
    
    writexl::write_xlsx(twostep, "01_MR_related/results/mr_evidence_outputs/mreve_twostep_all_mr_results_post_validation.xlsx")



## split into MR and sens and save as is

# step 1
res<-read_tsv("01_MR_related/results/mr_evidence_outputs/redone_MRmeds_fulloutput_fulltable_step1_V3.tsv")
redone_MR <- res %>% 
  select(exposure:effect_direction,used_instrument, qval) 
redone_MR_sens <- res %>% 
  select(any_of(c('id.exposure', 'exposure',  'id.outcome','outcome', 'method', 'nsnp',
                  "egger_intercept", "egger_intercept_se","egger_intercept_pval"  , 
                  "Q" , "Q_df" ,"Q_pval",
                  'Fst', 'total_r2', 'correct_causal_direction', 'steiger_pval',  "snp_r2.exposure" ,"snp_r2.outcome" ))) %>% 
  filter(method != 'Weighted median')

write_tsv(redone_MR,      "01_MR_related/results/mr_evidence_outputs/redone_MRmeds_step1_V3.tsv")
write_tsv(redone_MR_sens, "01_MR_related/results/mr_evidence_outputs/redone_MRmeds_sens_step1_V3.tsv") 


# step 2
res<-read_tsv("01_MR_related/results/mr_evidence_outputs/redone_MRmeds_fulloutput_fulltable_step2_V3.tsv")
redone_MR <- res %>% 
  select(exposure:effect_direction, used_instrument, qval) 
redone_MR_sens <- res %>% 
  select(any_of(c('id.exposure', 'exposure','id.outcome', 'outcome', 'method', 'nsnp',
                  "egger_intercept", "egger_intercept_se","egger_intercept_pval"  , 
                  "Q" , "Q_df" ,"Q_pval",
                  'Fst', 'total_r2', 'correct_causal_direction', 'steiger_pval',  "snp_r2.exposure" ,"snp_r2.outcome" ))) %>% 
  filter(method != 'Weighted median')

write_tsv(redone_MR,      "01_MR_related/results/mr_evidence_outputs/redone_MRmeds_step2_V3.tsv")
write_tsv(redone_MR_sens, "01_MR_related/results/mr_evidence_outputs/redone_MRmeds_sens_step2_V3.tsv") 


## now need to find traits that are plausibe mediators when filterign by qvalue in both steps

# Qval
step1_medsQval <-read_tsv("01_MR_related/results/mr_evidence_outputs/redone_MRmeds_fulloutput_fulltable_step1_V3.tsv") %>% 
  filter(method %in% c('Inverse variance weighted', 'Wald ratio')) %>% 
  filter(qval< 0.05) %>% pull(id.outcome) %>% unique() # 330

step2_medsQval <-read_tsv("01_MR_related/results/mr_evidence_outputs/redone_MRmeds_fulloutput_fulltable_step2_V3.tsv") %>% 
  filter(method %in% c('Inverse variance weighted', 'Wald ratio')) %>% 
  filter(qval< 0.05) %>% pull(id.exposure) %>% unique() # 101

plausible_mediators_Qval<- intersect(step1_medsQval, step2_medsQval) # 97

# Pval
step1_medsPval <-read_tsv("01_MR_related/results/mr_evidence_outputs/redone_MRmeds_fulloutput_fulltable_step1_V3.tsv") %>% 
  filter(method %in% c('Inverse variance weighted', 'Wald ratio')) %>% 
  filter(pval< 0.05) %>% pull(id.outcome) %>% unique() # 391

step2_medsPval <-read_tsv("01_MR_related/results/mr_evidence_outputs/redone_MRmeds_fulloutput_fulltable_step2_V3.tsv") %>% 
  filter(method %in% c('Inverse variance weighted', 'Wald ratio')) %>% 
  filter(pval< 0.05) %>% pull(id.exposure) %>% unique() # 192

plausible_mediators_Pval<- intersect(step1_medsPval, step2_medsPval) # 180


# extract validated only

# Qval
step1_validated <- read_tsv("01_MR_related/results/mr_evidence_outputs/redone_MRmeds_fulloutput_fulltable_step1_V3.tsv") %>% 
  select(exposure, id.exposure , id.outcome, beta_CI,b, pval,qval, used_instrument, effect_direction , nsnp, method) %>% 
  filter(method %in% c('Inverse variance weighted', 'Wald ratio')) %>% 
  filter(id.outcome %in% plausible_mediators_Qval)

step2_validated <- read_tsv("01_MR_related/results/mr_evidence_outputs/redone_MRmeds_fulloutput_fulltable_step2_V3.tsv") %>% 
  select(exposure, id.exposure , id.outcome,outcome, OR_CI,or, pval,qval, used_instrument, effect_direction , nsnp, method) %>% 
  filter(method %in% c('Inverse variance weighted', 'Wald ratio')) %>% 
  filter(id.exposure %in% plausible_mediators_Qval)

write_tsv(step1_validated, "01_MR_related/results/mr_evidence_outputs/redone_MRmeds_subsetoutput_ivw_step1_validated_V3_Qval.tsv")
write_tsv(step2_validated, "01_MR_related/results/mr_evidence_outputs/redone_MRmeds_subsetoutput_ivw_step2_validated_V3_Qval.tsv")
step1_validated %>% count(exposure) %>% View() # must be 103-ish: 97

# Pval -
step1_validated <- read_tsv("01_MR_related/results/mr_evidence_outputs/redone_MRmeds_fulloutput_fulltable_step1_V3.tsv") %>% 
  select(exposure, id.exposure , id.outcome,outcome, beta_CI,b, pval,qval, used_instrument, effect_direction , nsnp, method) %>% 
  filter(method %in% c('Inverse variance weighted', 'Wald ratio')) %>% 
  filter(id.outcome %in% plausible_mediators_Pval)

step2_validated <- read_tsv("01_MR_related/results/mr_evidence_outputs/redone_MRmeds_fulloutput_fulltable_step2_V3.tsv") %>% 
  select(exposure, id.exposure , id.outcome,outcome, OR_CI,or, pval,qval, used_instrument, effect_direction , nsnp, method) %>% 
  filter(method %in% c('Inverse variance weighted', 'Wald ratio')) %>% 
  filter(id.exposure %in% plausible_mediators_Pval)

# ---NB this data also have to saved as supl data  - supl 6?? - save unselected cols version 
write_tsv(step1_validated, "01_MR_related/results/mr_evidence_outputs/redone_MRmeds_subsetoutput_ivw_step1_validated_V3_Pval.tsv")
write_tsv(step2_validated, "01_MR_related/results/mr_evidence_outputs/redone_MRmeds_subsetoutput_ivw_step2_validated_V3_Pval.tsv")


#### next: putting everything together


# bring in total effect data
total_effect <- read_tsv(paste0("01_MR_related/results/mr_evidence_outputs/all_multiple_testingV3.tsv")) %>% 
  filter(!is.na(id.outcome)) 
length(unique(total_effect$id.exposure)) # 50

total_effect <- total_effect %>% 
  mutate(total_inst = ifelse(grepl("r2", used_instrument), "cis", "genome-wide")) %>% 
  mutate(total_sign_level = case_when(qval < 0.05 ~ "FDR",
                                      pval < 0.05 ~ "nominal",
                                      TRUE ~ "not")) %>% 
  select( id.exposure, id.outcome,
          "total.OR_CI" = "OR_CI",
          "total.or" = "or",
          "total.pval"  =   "pval"   ,
          "total.qval" =  "qval", 
          "total.effect_direction" = "effect_direction" ,
          "total.nsnp" =  "nsnp"  , 
          "exposure_cat" ,
          'total_inst',
          'total_sign_level')



# re-load and update cols
step1_validated <- read_tsv("01_MR_related/results/mr_evidence_outputs/redone_MRmeds_subsetoutput_ivw_step1_validated_V3_Pval.tsv")
step1_validated <- step1_validated %>% 
  select('id.med' = 'id.outcome',
         'exposure',
         'id.exposure',
         "step1.beta_CI" ="beta_CI" ,
         "step1.beta" =   "b",
         "step1.pval" = "pval",
         "step1.qval"  =  "qval", 
         'step1.used_instrument' ='used_instrument',
         "step1.effect_direction" =  "effect_direction",
         "step1.nsnp" ="nsnp")

step2_validated <- read_tsv("01_MR_related/results/mr_evidence_outputs/redone_MRmeds_subsetoutput_ivw_step2_validated_V3_Pval.tsv")
step2_validated <- step2_validated %>% 
  select('id.med' = 'id.exposure',
         'mediator' = "exposure",
         'id.outcome', 
         'outcome',
         "step2.OR_CI" ="OR_CI" ,
         "step2.or" =   "or",
         "step2.pval" = "pval",
         "step2.qval"  =  "qval", 
         'step2.used_instrument' ='used_instrument',
         "step2.effect_direction" =  "effect_direction",
         "step2.nsnp" ="nsnp") 

steps_joinedQval <- full_join(step1_validated,step2_validated, by = "id.med") %>% 
  select(exposure, id.exposure, mediator, id.med, starts_with("step1"), outcome, id.outcome, starts_with("step2"), everything()) %>% 
  left_join(total_effect, by = c("id.exposure", "id.outcome")) %>% 
  filter(step1.qval < 0.05 & step2.qval < 0.05  & !is.na(total.effect_direction))


steps_joinedPval <- full_join(step1_validated,step2_validated, by = "id.med") %>% 
  select(exposure, id.exposure, mediator, id.med, starts_with("step1"), outcome, id.outcome, starts_with("step2"), everything()) %>% 
  left_join(total_effect, by = c("id.exposure", "id.outcome")) %>% 
  filter(step1.pval < 0.05 & step2.pval < 0.05  & !is.na(total.effect_direction))

med_cats<- steps_joinedPval %>% 
  select(exposure.trait = mediator, exposure.id = id.med, exposure.trait = mediator) %>% 
  create_exposure_categories() %>% 
  select(id.med = exposure.id, med_cat=exposure_cat ) %>% distinct() 
  
steps_joinedPval <- left_join(steps_joinedPval, med_cats)

# no med per trait

comparisonPval<- steps_joinedPval %>% 
  select(exposure_cat,exposure,id.exposure,total.nsnp, total_inst,total_sign_level, id.med) %>% distinct() %>%
  count(exposure_cat,exposure,id.exposure,total.nsnp,  total_inst,  total_sign_level) %>% rename(med_count_P=n) 

comparisonQval<- steps_joinedQval %>% 
  select(exposure_cat,exposure,id.exposure, total.nsnp, total_inst, total_sign_level, id.med) %>% distinct() %>%
  count(exposure_cat,exposure,id.exposure, total.nsnp, total_inst, total_sign_level) %>% rename(med_count_Q=n) 


lit_counts<-read_tsv("02_literature_related//results/literature_outputs/lit_space_statsV3.tsv") %>% 
  select(id.exposure,exposure, unique_triples) %>% 
  mutate(unique_triples =ifelse(is.na(unique_triples), 0, unique_triples))


counts<- full_join(comparisonPval, comparisonQval) %>% left_join(lit_counts) %>% 
         mutate(med_count_P =ifelse(is.na(med_count_P), 0, med_count_P)) %>% 
         mutate(med_count_Q =ifelse(is.na(med_count_Q), 0, med_count_Q))


length(unique(counts$id.exposure)) # 101 yes this is bcac 2017 risk factors only

counts %>% filter(med_count_P > 0) %>% pull(id.exposure) %>% unique() %>% length() # 101 - all exposures have at least 1 putative mediaotr if pval is used
counts %>% filter(med_count_Q > 0) %>% pull(id.exposure) %>% unique() %>% length()  # 64
counts %>% filter(total_sign_level == "FDR") %>% filter(med_count_Q > 0) %>% pull(id.exposure) %>% unique() %>% length()  # 39

counts  %>% write_tsv("01_MR_related/results/mr_evidence_outputs/mediators_counts_per_traitsV3.csv")


counts <- read_tsv("01_MR_related/results/mr_evidence_outputs/mediators_counts_per_traitsV3.csv")
counts_prots <- counts %>% 
  filter(exposure_cat == "Proteins") %>% 
  filter(total_sign_level %in% c("FDR","nominal")) %>%
  #filter( total_inst == "cis") %>% 
  arrange(total_sign_level, total_inst, desc(med_count_P)) 

counts_prots %>% write_tsv("01_MR_related/results/mr_evidence_outputs/mediators_counts_per_traitsV3proteins.tsv")


# saving all mediaotrs per each trait

exposure_to_extract <- unique(steps_joinedPval$id.exposure)

out <- list()
for (i in exposure_to_extract){
  sub <- steps_joinedPval %>%  filter(id.exposure == i)
  out[[i]] <- sub
}

out2 <- out

writexl::write_xlsx(out, "01_MR_related/results/mr_evidence_outputs/med-table-validatedV3.xlsx") #  this Supl table 6a

out[["ukb-d-30760_irnt"]] %>%  filter(step1.qval < 0.05, step2.qval < 0.05) %>% distinct() %>% filter(med_cat != "Lipids") %>% View()

out[["ukb-d-30760_irnt"]] %>% filter(total.pval < 0.05) %>% View()





### old


p <- ggplot(counts, aes(x=med_count)) + 
  geom_histogram()+ facet_wrap(~exposure_cat)


# ad hoc category size checking
counts %>% 
filter(total_sign_level=="FDR") %>%
  filter(!is.na(med_count_Q)) %>% 
  pull(id.exposure) %>% unique() %>% length()# 39/49 have mediators at Q value

counts %>% 
  filter(total_sign_level=="FDR") %>% 
  group_by(exposure_cat) %>%
  mutate(across(everything(), ~replace_na(.x, 0))) %>%  
  summarise(meanP = mean(med_count_P), meanQ = mean(med_count_Q)) %>% View()  # this is Supl Table 7



counts %>% filter(id.exposure %in% c("prot-a-1369", "prot-a-1097",  "prot-b-38", "prot-a-2396", # proteins with lit space
                                     "prot-a-2629", "prot-a-1530",  "prot-a-67", "prot-a-366", 
                                     "prot-a-710", "prot-b-71",  "prot-a-1148", "prot-a-1117", 
                                     "prot-a-2073", "prot-a-3076")) %>% View()

#####
###
###
###


# extracting stuff for each case study
# 
# specify ID:
case_id = "prot-a-1576" # cBMI

# this is already fdr-corrected in both steps; mr-eve
res<- results_subset %>% 
  filter(type %in% c('mediator')) %>% # conf coud be here too if it was extracted
  filter(!med.id %in% exclude_meds$exposure.id) %>% 
  filter(exposure.id == case_id) 

res %>% select(med.trait, med.id, med_cat) %>% distinct() %>% count(med_cat) %>% View() # sum is total from mr-eve


  
mreve_all<-res %>% select(exposure.trait,med.trait, med.id, med_cat, step1.OR_CI=r1.OR_CI, step2.OR_CI=r3.OR_CI, step1.qval=r1.qval, step2.qval=r3.qval, outcome.id) %>% 
  distinct() %>% arrange(med.trait,med_cat)

#mreve_all %>% write_csv("~/Desktop/chronotype_mreve.csv")


res_valid<- steps_joinedPval %>% 
  filter(id.exposure == case_id) 

#res_valid %>% write_csv("~/Desktop/chronotype_mr_mediators_validated.csv")



res_valid %>% select(mediator, id.med, med_cat) %>% distinct() %>% count(med_cat) %>% View() # pval

res_valid %>% filter(step1.qval <0.05 & step2.qval <0.05) %>%  select(mediator, id.med, med_cat) %>% distinct() %>% count(med_cat) %>% View() #qval


res_valid %>% filter(step1.pval <0.05 & step2.pval <0.05) %>% View()

res_valid %>% filter(step1.qval <0.05 & step2.qval <0.05) %>% View()

res_validP<- steps_joinedPval %>% 
  filter(id.exposure == case_id) %>% 
  filter(step1.pval <0.05 & step2.pval <0.05) %>% 
  #filter(!med_cat %in% c("Lipids", "Antrophometric")) %>%  # for CBMI
  distinct() %>% 
  arrange(outcome, med_cat)

res_validP %>% select(id.med, med_cat) %>% distinct() %>% count(med_cat) %>% View()


res_validQ<- steps_joinedPval %>% 
  filter(id.exposure == case_id) %>% 
  filter(step1.qval <0.05 & step2.qval <0.05) %>% 
  #filter(!med_cat %in% c("Lipids")) %>%  # for CBMI
  distinct() %>% 
  arrange(outcome, med_cat)

res_validQ %>% select(id.med, med_cat) %>% distinct() %>% count(med_cat) %>% View()


output_table <- res_validQ %>%
              select("exposure",  "mediator", "outcome",  "step1.beta_CI" , "step1.pval", "step1.qval", "step2.OR_CI"  , "step2.pval", "step2.qval", "mediator_id"='id.med') %>% 
               mutate(outcome = case_when(grepl("ER-", outcome)  ~ "ER- subtype" ,
                                          grepl("ER+" , outcome) ~  "ER+ subtype",
                                                TRUE ~ "Overall BC")) %>% 
              mutate(step1.pval = signif(step1.pval, digits=2)) %>% 
              mutate(step2.pval = signif(step2.pval, digits=2)) %>% 
              mutate(step1.qval = signif(step1.qval, digits=2)) %>% 
              mutate(step2.qval = signif(step2.qval, digits=2)) %>% 
              # pheno names specifc chanages
              mutate(mediator = ifelse(mediator == "Types of physical activity in last 4 weeks: Other exercises (eg: swimming  cycling  keep fit  bowling)",
                                       "Physical activity: exercise to keep fit", mediator)) %>% 
            mutate(exposure = ifelse(exposure == "Immunoglobulin superfamily containing leucine-rich repeat protein 2",
                                     "ISLR2", exposure)) %>% 
              mutate(mediator = paste0(mediator, "\n", mediator_id)) %>% 
              arrange(mediator)

output_table %>%  write_csv(paste0("01_MR_related/results/mr_evidence_outputs/case_studies_V3/", case_id, "_twostepMR_FDR_subset.csv"))
