
##### this makes supl table 1


library(tidyverse)
library(TwoSampleMR)
source("helper_functions.R")
source("01_MR_related/scripts/mr_related_functions.R")
source("01_MR_related/scripts/app1_MR-EvE_app/functions.R")

full_res <-  read_csv("01_MR_related//results/mr_evidence_outputs/all_mreve_bc_resultsV3.csv")
length(unique(full_res$exposure.id)) # 2218


full_res <-  full_res %>%  # from 01_mr_epigraph_query.R script
  # convert MR results to OR
  tidy_display_numbers()%>% 
  mutate(outome = case_when(outcome.id == "ieu-a-1126"~ 'Brest cancer overall',
                            outcome.id == "ieu-a-1127"~ 'Brest cancer ER+',
                            outcome.id == "ieu-a-1128"~ 'Brest cancer ER-')) %>% 
  create_exposure_categories() %>% 
  add_exposure_labels()  %>% 
  select(-c(pval_truncated, log10pval_trunc, empty_col,exposure.ss, tmp,exposure.ss, ukb_tag)) 


dat <-  read_tsv("01_MR_related/scripts/app1_MR-EvE_app/data_copy/bc_all_mr_fromCIsV3.tsv")
dat378 <- dat %>% pull(exposure.id) %>% unique()


x <- read_tsv("01_MR_related//results/mr_evidence_outputs/passed_multiple_testingV3.tsv") %>% # 202
    mutate(bcac = ifelse(grepl("GWAS", outcome), "2017", "2020"))

tidy202 <- x %>% pull(id.exposure) %>% unique()


valid_total_129 <- x %>% filter(pval < 0.05) %>% pull(id.exposure) %>% unique()
valid2017_101 <- x %>% filter(bcac == "2017", pval < 0.05) %>% pull(id.exposure) %>% unique()
valid2020_109 <- x %>% filter(bcac == "2020", pval < 0.05)%>% pull(id.exposure) %>% unique()
fdr_total_63 <- x %>% filter(qval < 0.05) %>% pull(id.exposure) %>% unique()
fdr2017_50 <- x %>% filter(bcac == "2017", qval < 0.05) %>% pull(id.exposure) %>% unique()
fdr2020_52 <- x %>% filter(bcac == "2020", qval < 0.05)%>% pull(id.exposure) %>% unique()




out<-
full_res %>% 
  mutate(traits_378 = ifelse(exposure.id %in% dat378, T,F)) %>% 
  mutate(traits_202 = ifelse(exposure.id %in% tidy202, T,F)) %>% 
  mutate(traits_129 = ifelse(exposure.id %in% valid_total_129, T,F)) %>% 
  mutate(traits_101 = ifelse(exposure.id %in% valid2017_101, T,F)) %>% 
  mutate(traits_109 = ifelse(exposure.id %in% valid2020_109, T,F)) %>% 
  mutate(traits_63 = ifelse(exposure.id %in% fdr_total_63, T,F)) %>% 
  mutate(traits_50 = ifelse(exposure.id %in% fdr2017_50, T,F)) %>% 
  mutate(traits_52 = ifelse(exposure.id %in% fdr2020_52, T,F)) %>% 
  mutate(exposure_cat = ifelse(!exposure.id %in% tidy202, "Other",exposure_cat ))


writexl::write_xlsx(out, "01_MR_related/results/mr_evidence_outputs/mreve_master_table_w_filters_V3.xlsx") # supl data 1


out %>% count(exposure_cat)




#### also metadata: supl data 8

meta <-  read_csv("01_MR_related//results/mr_evidence_outputs/all_mreve_bc_resultsV3.csv") %>% select(exposure = exposure.trait, exposure.id) %>% distinct() %>% 
        left_join(out %>% select(exposure.id, exposure_cat)) %>% 
         mutate(traits_129_validated = ifelse(exposure.id %in% valid_total_129, T,F)) %>% 
        arrange(exposure_cat, desc(traits_129_validated)) %>% select(exposure_cat, everything()) %>% distinct()

ao <- available_outcomes()

meta_upd <- left_join(meta, ao, by = c('exposure.id'='id', 'exposure' = 'trait')) %>% 
   select(-c(study_design, priority,coverage, doi, mr, ontology, group_name)) %>% 
  select(colnames(meta), author, consortium, year,sex,population,sample_size,  ncase, ncontrol, nsnp, unit, sd, everything() )


write_csv(meta_upd, "01_MR_related/results/mr_evidence_outputs/mrever_metadata_for_suppl.csv")



