library(tidyverse)
library(epigraphdb)
source("helper_functions.R")
source("01_MR_related/scripts/mr_related_functions.R")
source("01_MR_related/scripts/app1_MR-EvE_app/functions.R")

# testing example
exposure = "prot-a-670"
outcome = "ieu-a-1127"
pval_threshold = 1e-01
#test<-query_epigraphdb_as_table(mediator_query)
#xx<-tidy_conf_query_output(test, type = "mediator")



#### 

# read effect direction matrix from 03..

traits_df <- read_tsv("01_MR_related/results/mr_evidence_outputs/trait_manual_ivw_subtypes_merged.tsv") %>% 
  select(id.exposure, contains('ieu')) # only for BCAC 2017 cols, because mr-eve is available only for them

# extract pairs of exp-out with effect
traits_all <- bind_rows(
  traits_df %>%filter(`ieu-a-1126` !=0) %>%  select(id.exposure) %>% mutate(id.outcome = "ieu-a-1126"),
  traits_df %>%filter(`ieu-a-1127` !=0) %>%  select(id.exposure) %>% mutate(id.outcome = "ieu-a-1127"),
  traits_df %>%filter(`ieu-a-1128` !=0) %>%  select(id.exposure) %>% mutate(id.outcome = "ieu-a-1128")) %>% distinct()
traits_all %>% pull(id.exposure) %>%  unique() %>% length() # 171

# for each pair collect intermediates (all 4 scenarios)
all_results <- tibble()
for (i in 1:length(traits_all$id.exposure)){
  
  out <- query_and_tidy_conf(exposure = traits_all$id.exposure[i], 
                             outcome = traits_all$id.outcome[i], 
                             pval_threshold =  0.01, ##### keeping almost all ( NB med->BC has no threshold, as we are gonna use prev valiadation as a filter)
                             mediator = T)  # only mediators
  all_results<- bind_rows(all_results, out)
}



# add trait category to intermediate items
all_results_trait_cats <- 
   all_results %>%  
   select(exposure.trait = med.trait, exposure.id = med.id) %>% 
   create_exposure_categories() %>% 
   select(med_cat = exposure_cat, med.trait = exposure.trait, med.id = exposure.id ) %>% distinct()
                  
all_results <- all_results %>% left_join(all_results_trait_cats) 
dim(all_results) # 18693

# filter
results_subset <- all_results %>% 
  filter(!r1.OR_CI %in% c("0 [0:0]", "1 [1:1]")) %>% 
  filter(!r3.OR_CI %in% c("0 [0:0]", "1 [1:1]")) %>% 
  select(exposure.trait, med.trait, outcome.id, r1.OR_CI, r2.OR_CI, r3.OR_CI,  type, med_cat, r1.b, r2.b, r3.b, exposure.id, med.id) %>% 
  distinct() %>% 
  filter(!grepl("arm|leg", med.trait, ignore.case = T)) %>% 
  filter(!med_cat %in% c('Diseases', 'Medical Procedures', 'other', "eye_hearing_teeth", "Education", "Psychology", "CHD")) %>% 
  filter(med.id %in% traits_df$id.exposure)   # this keeps only those mediaotrs that are accepted exposure in main analysis validation
  
dim(results_subset) # 8180

counts_raw <- results_subset %>% select(exposure.trait,exposure.id, med.id) %>% distinct() %>% count(exposure.trait,exposure.id) %>% rename(med_count_preval=n)
mean(counts_raw$med_count_preval)

#write_csv(results_subset, "01_MR_related/mr_evidence_outputs/conf_med_extracted.csv") # p<10e4
write_csv(results_subset, "01_MR_related/results/mr_evidence_outputs/med_extracted_all_r3.csv") # p<0.05 # r3 all is not pval restricted


protein_names <- read_csv("01_MR_related/results/mr_evidence_outputs/protein_names.csv", col_names = F) %>% rename(name = X1, gene = X2)
results_subset<- read_csv("01_MR_related/results/mr_evidence_outputs/med_extracted_all_r3.csv") %>% 
                filter(type %in% c('mediator')) %>%
                create_exposure_categories() %>%
                left_join(protein_names, by = c("med.trait" = 'name')) %>% rename(med.gene = gene) %>% 
                left_join(protein_names, by = c("exposure.trait" = 'name')) %>% rename(exp.gene = gene) %>% 
                select(exposure.trait,exp.gene, med.trait, med.gene, everything()) 
dim(results_subset)# 8180

length(unique(results_subset$exposure.id)) # 162
# next need to validate E->M for mediators (and M->E for confounders)
             
#
#
#
#### -- validation step happens here

#  in 04sub_mreve_mediators_validation.R
#
#
#





## 1) bringing in validated  results for trait-trait
validated <- read_tsv("01_MR_related/results/mr_evidence_outputs/redone_MRmeds_subsetoutput_ivw.tsv") %>% 
              select(id.exp_val = id.exposure , id.out_val = id.outcome, r1.OR_CI_val = OR_CI, r1.nsnp_val=nsnp)

# need to do left_join with conf and med separately
mediators <- results_subset %>% filter(type == 'mediator') %>% left_join(validated, by = c("exposure.id"= 'id.exp_val' , 'med.id'= 'id.out_val'))
confounders <- results_subset %>% filter(type == 'confounder') %>% left_join(validated, by = c("exposure.id"= 'id.out_val', 'med.id'= 'id.exp_val'))

results_validated <-bind_rows(mediators, confounders) # meds only 


## 2) bringing in validated  results for trait-BC
traits_bc_validated <- read_tsv("01_MR_related/results/mr_evidence_outputs/redone_MR_subsetoutput_ivw.tsv") %>%  # from 03sub
                      select(id.trait = id.exposure, id.outcome, OR_CI_val = OR_CI, nsnp_val =nsnp, exposure_cat) 

## 3) joined trait-trait and -bc validation
validated_with_BC <- results_validated %>% 
                    left_join(traits_bc_validated, by = c("exposure.id" = "id.trait", "outcome.id" = "id.outcome" )) %>% rename('r2.OR_CI_val' = 'OR_CI_val', "r2.nsnp_val"="nsnp_val") %>% 
                    left_join(traits_bc_validated, by = c("med.id" = "id.trait", "outcome.id" = "id.outcome" )) %>% rename('r3.OR_CI_val' = 'OR_CI_val', "r3.nsnp_val"="nsnp_val") %>% 
                    filter(!is.na(r1.OR_CI_val)) %>% 
                    filter(!is.na(r2.OR_CI_val)) %>% 
                    filter(!is.na(r3.OR_CI_val)) %>% 
                    select(-r1.OR_CI, -r2.OR_CI,-r3.OR_CI, -r1.b, -r2.b, -r3.b, -exposure_cat.y , -exposure_cat) %>% 
                    rename(exposure_cat = exposure_cat.x)
  
    
dim(validated_with_BC) #4113
 
### add lit size ---- optional / testing
#lit <- read_tsv("02_literature_related/literature_outputs/traits_marked_for_lit_analysis.tsv") %>% select(id, unique_pairs)
#
#validated_with_BC <- validated_with_BC %>%
#  left_join(lit, by = c("exposure.id" = "id")) %>% rename('exposure_lit_pairs' = 'unique_pairs') %>% 
#  left_join(lit, by = c("med.id" = "id")) %>% rename('med_lit_pairs' = 'unique_pairs')
#


# rearrange
validated_with_BC <- validated_with_BC %>%
select(type, outcome.id, 
       exposure.trait,exposure.id,    exposure_cat,
       med.trait,med.id, med_cat,
       #med.gene, #exp.gene,
       #exposure_lit_pairs, med_lit_pairs, 
       r1.OR_CI_val, r1.nsnp_val, r2.OR_CI_val, r2.nsnp_val, r3.OR_CI_val, r3.nsnp_val)
    
       
       
       


counts <- left_join(
  # no conf per trait
  validated_with_BC %>% filter(type == 'confounder')  %>% 
    select(exposure_cat, exposure.trait,exposure.id, med.id) %>% distinct() %>%  count(exposure_cat,exposure.trait,exposure.id) %>% rename(conf_count=n),
  
  # no med per trait
  validated_with_BC %>% filter(type == 'mediator')  %>% 
    select(exposure_cat,exposure.trait,exposure.id, med.id) %>% distinct() %>% count(exposure_cat,exposure.trait,exposure.id) %>% rename(med_count=n)
)

length(unique(validated_with_BC$exposure.id)) # 153 

# no med per trait 153
counts <- validated_with_BC %>% 
  filter(type == 'mediator')  %>% 
  select(exposure_cat,exposure.trait,exposure.id, med.id) %>% distinct() %>% count(exposure_cat,exposure.trait,exposure.id) %>% rename(med_count=n)

mean(counts$med_count)
counts %>% write_tsv("01_MR_related/results/mr_evidence_outputs/mediators_counts_per_traits.csv")

p <- ggplot(counts, aes(x=med_count)) + 
  +     geom_histogram()+ facet_wrap(~exposure_cat)

# ad hoc category checking
x<- counts %>% filter(exposure_cat == 'Lipids')  
mean(x$med_count)# 31


mean(counts$med_count)

# saving all mediaotrs per each trait

exposure_to_extract <- unique(validated_with_BC$exposure.id)

out <- list()

for (i in exposure_to_extract){
  sub <- validated_with_BC %>%  filter(exposure.id == i)
  out[[i]] <- sub
}
writexl::write_xlsx(out, "01_MR_related/results/mr_evidence_outputs/med-table-validated.xlsx")















######




## rules of mediaotrs
results_subset %>%  filter(type == 'mediator') %>% View()

mediators <- results_subset %>% 
  filter(type == 'mediator') %>% 
  # keep only those meds that have proven effect on outcome
 
  mutate(med_version = case_when(
    r1.b < 0 & r2.b < 0 & r3.b > 0 ~ "med1",
    r1.b > 0 & r2.b < 0 & r3.b < 0 ~ "med2",
    r1.b > 0 & r2.b > 0 & r3.b > 0 ~ "med3",
    r1.b < 0 & r2.b > 0 & r3.b < 0 ~ "med4",
    
    r1.b > 0 & r2.b < 0 & r3.b > 0 ~ "med5",
    r1.b < 0 & r2.b > 0 & r3.b > 0 ~ "med6",
    r1.b < 0 & r2.b < 0 & r3.b < 0 ~ "med7",
    r1.b > 0 & r2.b > 0 & r3.b < 0 ~ "med8",
    TRUE ~ "med_other"))
mediators %>% count(med_version)


## rules of conf
confounders <- results_subset %>% 
  filter(type == 'confounder') %>% 
  # keep only those meds that have proven effect on outcome
  filter(med.id %in% traits_df$id.exposure) 






