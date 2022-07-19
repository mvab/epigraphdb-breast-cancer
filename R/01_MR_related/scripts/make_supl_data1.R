

full_res <-  read_csv("01_MR_related/results/mr_evidence_outputs/all_mreve_bc_results.csv")
length(unique(full_res$exposure.id)) # 2332


dat <-  read_tsv("01_MR_related/scripts/app1_MR-EvE_app/data_copy/bc_all_mr_fromCIs.tsv") %>%  # from 01_mr_epigraph_query.R script
  # subset 
  filter(exposure.sex != 'Males') %>% 
  #filter(mr.method != 'Steiger null') %>%  # NA results --- does not make difference to the total number of traits
  # convert MR results to OR
  tidy_display_numbers()%>% 
  # deal with al outcome related changes
  process_bc_outcomes() %>% 
  # create categories of exposure traits
  filter(!grepl("_raw",exposure.id)) %>% 
  create_exposure_categories() %>% 
  add_exposure_labels()  # 1643

list1643<-unique(dat$exposure.id)


selected_cats<- read_tsv("01_MR_related/results/mr_evidence_outputs/tidy_traits_by_cat.tsv") # 905
length(unique(selected_cats$exposure.id)) 
selected_cats %>% count(exposure_cat)
list905<-unique(selected_cats$exposure.id)
cats905<- selected_cats %>% select(exposure.id, trait_category = exposure_cat) %>% distinct()
dim(cats905)

followup <- read_tsv("01_MR_related/results/mr_evidence_outputs/trait_for_followup.tsv") %>% rename(exposure.id=id)
length(unique(followup$exposure.id)) # 309
list309 <- unique(followup$exposure.id)

validated <- read_tsv("01_MR_related/results/mr_evidence_outputs/redone_MR_subsetoutput_ivw.tsv")%>% rename(exposure.id=id.exposure)
length(unique(validated$exposure.id)) # 171
list171<-unique(validated$exposure.id)


antro_blacklist <- c('ieu-a-81','ieu-a-74', "ieu-a-73" ,"ieu-a-79" ,"ieu-a-72" ,"ieu-a-78",
                     'ieu-a-63',  'ieu-a-66', 'ieu-a-60',  'ieu-a-69' , 'ieu-a-61',
                     'ieu-a-54', 'ieu-a-55',  'ieu-a-49' , 'ieu-a-48', 'ieu-a-57' , 'ieu-a-50',
                     'ukb-b-12039', 'ukb-b-2303', 'ieu-a-2', 'ieu-a-835',
                     'ukb-b-18105', 'ieu-a-99', 'ieu-a-105', 'ukb-b-4650', 'ieu-a-95',
                     'ukb-a-248', 'ieu-a-68', 	
                     'ieu-a-101', 'ieu-a-109')



out<-
full_res %>% 
  mutate(effect_in_one_trait_1643_traits = ifelse(exposure.id %in% list1643, T,F)) %>% 
  mutate(in_12_categories_905_traits = ifelse(exposure.id %in% list905, T,F)) %>% 
  left_join(cats905) %>% 
  mutate(consistent_effect_309_traits = ifelse(exposure.id %in% list309, T,F)) %>% 
  mutate(validatated_in_BCAC2017_171_traits = ifelse(exposure.id %in% list171, T,F)) %>% 
  mutate(excluded_in_fig2_as_redundant = ifelse(exposure.id %in% antro_blacklist, T,F))


out %>% select(exposure.id, `validatated_in_BCAC2017_171_traits`) %>% filter(validatated_in_BCAC2017_171_traits==T) %>% distinct() %>% dim()
dim(out)
writexl::write_xlsx(out, "01_MR_related/results/mr_evidence_outputs/mreve_master_table_w_filters.xlsx")



