library(tidyverse)
library(epigraphdb)
source("helper_functions.R")
source("01_MR_related/explore_MR-EvE_app/functions.R")

# testing
exposure = "met-a-362"
outcome = "ieu-a-1128"
pval_threshold = 1e-04
#test<-query_epigraphdb_as_table(mediator_query)
#xx<-tidy_conf_query_output(test, type = "mediator")



tidy_conf_query_output <- function(df, type){
  
  if (dim(df)[1] != 0){
    
    df_OR <- df %>% 
    mutate( r1.loci = r1.b - 1.96 * r1.se, 
            r1.upci = r1.b + 1.96 * r1.se,
            r1.or = exp(r1.b), 
            r1.or_loci = exp(r1.loci), 
            r1.or_upci = exp(r1.upci),
            r1.OR_CI = paste0(round(r1.or,2), " [",round(r1.or_loci,2) ,":",round(r1.or_upci,2), "]")) %>% 
      mutate(r1.effect_direction = ifelse(r1.or_loci > 1 & r1.or_upci >= 1, 'positive',
                                   ifelse(r1.or_loci < 1 & r1.or_upci <= 1, 'negative', 'overlaps null'))) %>% 
      mutate( r2.loci = r2.b - 1.96 * r2.se, 
              r2.upci = r2.b + 1.96 * r2.se,
              r2.or = exp(r2.b), 
              r2.or_loci = exp(r2.loci), 
              r2.or_upci = exp(r2.upci),
              r2.OR_CI = paste0(round(r2.or,2), " [",round(r2.or_loci,2) ,":",round(r2.or_upci,2), "]")) %>% 
      mutate(r2.effect_direction = ifelse(r2.or_loci > 1 & r2.or_upci >= 1, 'positive',
                                   ifelse(r2.or_loci < 1 & r2.or_upci <= 1, 'negative', 'overlaps null'))) %>% 
    
      mutate( r3.loci = r3.b - 1.96 * r3.se, 
              r3.upci = r3.b + 1.96 * r3.se,
              r3.or = exp(r3.b), 
              r3.or_loci = exp(r3.loci), 
              r3.or_upci = exp(r3.upci),
              r3.OR_CI = paste0(round(r3.or,2), " [",round(r3.or_loci,2) ,":",round(r3.or_upci,2), "]")) %>% 
      mutate(r3.effect_direction = ifelse(r3.or_loci > 1 & r3.or_upci >= 1, 'positive',
                                   ifelse(r3.or_loci < 1 & r3.or_upci <= 1, 'negative', 'overlaps null'))) %>% 
      
      select(-contains("loci"), -contains("upci")) %>% 
    
      select("exposure.trait", "exposure.id", "outcome.trait", "outcome.id", "med.trait","med.id" , 
             starts_with("r1"),  starts_with("r2"),starts_with("r3")) %>% 
      filter(r1.effect_direction != 'overlaps null' &r2.effect_direction != 'overlaps null' & r3.effect_direction != 'overlaps null') %>% 
      mutate(type = type)
  } else {
    df_OR <- data.frame()
  }
    
}

query_and_tidy_conf <- function(exposure, outcome,pval_threshold ){
  

  # construct queries
  confounder_query = paste0("
    MATCH (med:Gwas)-[r1:MR_EVE_MR]-> (exposure:Gwas) -[r2:MR_EVE_MR]->(outcome:Gwas) <-[r3:MR_EVE_MR]-(med:Gwas) 
    WHERE exposure.id = '", exposure, "' AND outcome.id = '", outcome,"' AND (not (toLower(med.trait) contains 'breast')) 
    AND r1.pval < ", pval_threshold, " AND r3.pval < ", pval_threshold, " 
    AND med.id <> exposure.id AND med.id <> outcome.id AND exposure.id <> outcome.id AND med.trait <> exposure.trait AND med.trait <> outcome.trait AND exposure.trait <> outcome.trait 
    RETURN exposure {.id, .trait}, outcome {.id, .trait}, med {.id, .trait}, r1 {.b, .se, .pval, .selection, .method, .moescore}, r2 {.b, .se, .pval, .selection, .method, .moescore}, r3 {.b, .se, .pval, .selection, .method, .moescore} ORDER BY r1.p
          ")
  
  collider_query = paste0("
    MATCH (med:Gwas)<-[r1:MR_EVE_MR]- (exposure:Gwas) -[r2:MR_EVE_MR]->(outcome:Gwas) -[r3:MR_EVE_MR]->(med:Gwas) 
    WHERE exposure.id = '", exposure, "' AND outcome.id = '", outcome,"' AND (not (toLower(med.trait) contains 'breast')) 
    AND r1.pval < ", pval_threshold, " AND r3.pval < ", pval_threshold, " 
    AND med.id <> exposure.id AND med.id <> outcome.id AND exposure.id <> outcome.id AND med.trait <> exposure.trait AND med.trait <> outcome.trait AND exposure.trait <> outcome.trait 
    RETURN exposure {.id, .trait}, outcome {.id, .trait}, med {.id, .trait}, r1 {.b, .se, .pval, .selection, .method, .moescore}, r2 {.b, .se, .pval, .selection, .method, .moescore}, r3 {.b, .se, .pval, .selection, .method, .moescore} ORDER BY r1.p
          ")
  
  mediator_query = paste0("
    MATCH (med:Gwas)<-[r1:MR_EVE_MR]- (exposure:Gwas) -[r2:MR_EVE_MR]->(outcome:Gwas) <-[r3:MR_EVE_MR]-(med:Gwas) 
    WHERE exposure.id = '", exposure, "' AND outcome.id = '", outcome,"' AND (not (toLower(med.trait) contains 'breast')) 
    AND r1.pval < ", pval_threshold, " AND r3.pval < ", pval_threshold, " 
    AND med.id <> exposure.id AND med.id <> outcome.id AND exposure.id <> outcome.id AND med.trait <> exposure.trait AND med.trait <> outcome.trait AND exposure.trait <> outcome.trait 
    RETURN exposure {.id, .trait}, outcome {.id, .trait}, med {.id, .trait}, r1 {.b, .se, .pval, .selection, .method, .moescore}, r2 {.b, .se, .pval, .selection, .method, .moescore}, r3 {.b, .se, .pval, .selection, .method, .moescore} ORDER BY r1.p
          ")
  
  reverse_mediator_query = paste0("
    MATCH (med:Gwas)-[r1:MR_EVE_MR]-> (exposure:Gwas) -[r2:MR_EVE_MR]->(outcome:Gwas) -[r3:MR_EVE_MR]->(med:Gwas) 
    WHERE exposure.id = '", exposure, "' AND outcome.id = '", outcome,"' AND (not (toLower(med.trait) contains 'breast')) 
    AND r1.pval < ", pval_threshold, " AND r3.pval < ", pval_threshold, " 
    AND med.id <> exposure.id AND med.id <> outcome.id AND exposure.id <> outcome.id AND med.trait <> exposure.trait AND med.trait <> outcome.trait AND exposure.trait <> outcome.trait 
    RETURN exposure {.id, .trait}, outcome {.id, .trait}, med {.id, .trait}, r1 {.b, .se, .pval, .selection, .method, .moescore}, r2 {.b, .se, .pval, .selection, .method, .moescore}, r3 {.b, .se, .pval, .selection, .method, .moescore} ORDER BY r1.p
          ")
  
  # run queries
  res_list <- list(
    colliders = query_epigraphdb_as_table(collider_query),
    confounders = query_epigraphdb_as_table(confounder_query),
    mediators = query_epigraphdb_as_table(mediator_query),
    rev_mediators = query_epigraphdb_as_table(reverse_mediator_query)
  )
  
  # tidy query outputs
  res_list_tidy <- list(
    colliders = tidy_conf_query_output(res_list$colliders, type = "collider"),
    confounders = tidy_conf_query_output(res_list$confounders, type = "confounder"),
    mediators =  tidy_conf_query_output(res_list$mediators, type = "mediator"),
    rev_mediators = tidy_conf_query_output(res_list$rev_mediators, type = "reverse_mediator")
  )
  
  res_list_tidy_df <- bind_rows(res_list_tidy)
  return(res_list_tidy_df)
  
}


#### 

traits_df <- read_tsv("01_MR_related/mr_evidence_outputs/trait_manual_ivw_subtypes_merged.tsv") 

# wxtract pairs of exp-out with effect
traits_all <- bind_rows(
  traits_df %>%filter(`ieu-a-1126` !=0) %>%  select(id.exposure) %>% mutate(id.outcome = "ieu-a-1126"),
  traits_df %>%filter(`ieu-a-1127` !=0) %>%  select(id.exposure) %>% mutate(id.outcome = "ieu-a-1127"),
  traits_df %>%filter(`ieu-a-1128` !=0) %>%  select(id.exposure) %>% mutate(id.outcome = "ieu-a-1128")) %>% distinct()


# for each pair collect intermediates
all_results <- tibble()
for (i in 1:length(traits_all$id.exposure)){
  
  out <- query_and_tidy_conf(exposure = traits_all$id.exposure[i], 
                             outcome = traits_all$id.outcome[i], 
                             pval_threshold =  1e-04)
  all_results<- bind_rows(all_results, out)
  
}

# add trait category to intermediate items
all_results_trait_cats <- 
  all_results %>%  
 select(exposure.trait = med.trait, exposure.id = med.id) %>% 
 create_exposure_categories() %>% 
 select(med_cat = exposure_cat, med.trait = exposure.trait, med.id = exposure.id ) %>% distinct()
                
all_results <- all_results %>% left_join(all_results_trait_cats) 
dim(all_results)

# filter
results_subset <- all_results %>% 
  filter(!r1.OR_CI %in% c("0 [0:0]", "1 [1:1]")) %>% 
  filter(!r3.OR_CI %in% c("0 [0:0]", "1 [1:1]")) %>% 
  select(exposure.trait, med.trait, outcome.id, r1.OR_CI, r2.OR_CI, r3.OR_CI,  type, med_cat, r1.b, r2.b, r3.b, exposure.id, med.id) %>% 
  distinct() %>% 
  filter(!grepl("arm|leg", med.trait, ignore.case = T)) %>% 
  filter(!med_cat %in% c('Diseases', 'Medical Procedures', 'other', "eye_hearing_teeth"))#, "Education", "Psychology", "CHD")) # other??
dim(results_subset)


write_csv(results_subset, "01_MR_related/mr_evidence_outputs/conf_med_extracted.csv")


## rules of mediaotrs
results_subset %>%  filter(type == 'mediator') %>% View()

mediators <- results_subset %>% 
  filter(type == 'mediator') %>% 
  # keep only those meds that have proven effect on outcome
  filter(med.id %in% traits_df$id.exposure) %>% 
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





#in mediaotrs/conf , we can keep onaly those that do have IVW confirmed effect on BC? 






###


in_final <- read_tsv("01_MR_related/mr_evidence_outputs/trait_manual_ivw_subtypes_merged.tsv") %>% 
  filter(exposure_cat %in% c("Proteins", "Metabolites", "Other biomarkers")) %>% 
  filter(!grepl("LDL|HDL|cholest|trigl", exposure_name, ignore.case = T)) %>% 
  pull(id.exposure)

# check if the third var was also a trai in the final set
tmp %>% mutate(also = ifelse(exposure.id %in% in_final, 1, 0)) %>% View()

query = paste0("
    MATCH (gene1:Gene)-[gene_to_protein:GENE_TO_PROTEIN]->(protein:Protein)-[protein_in_pathway:PROTEIN_IN_PATHWAY]->(pathway:Pathway) 
    WHERE gene1.name in ['MAN1A2', 'FGF7', 'CPXM1', 'IL1RL1', 'TPST2','SULF2','FLNA','DCBLD2','PRDM1','ISLR2']
    RETURN gene1.name, protein.uniprot_id,  pathway.name
    order by pathway.name
      ")


test<-query_epigraphdb_as_table(query)
test %>% count(pathway.name, sort=T)


