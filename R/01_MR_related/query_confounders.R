library(tidyverse)
library(epigraphdb)
source("helper_functions.R")
source("01_MR_related/explore_MR-EvE_app/functions.R")

# testing
exposure = "ukb-b-4650"
outcome = "ieu-a-1126"
pval_threshold = 1e-01
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

traits_df <- read_tsv("01_MR_related/mr_evidence_outputs/trait_manual_ivw_subtypes_merged.tsv") %>%  select(id.exposure, contains('ieu'))

# wxtract pairs of exp-out with effect
traits_all <- bind_rows(
  traits_df %>%filter(`ieu-a-1126` !=0) %>%  select(id.exposure) %>% mutate(id.outcome = "ieu-a-1126"),
  traits_df %>%filter(`ieu-a-1127` !=0) %>%  select(id.exposure) %>% mutate(id.outcome = "ieu-a-1127"),
  traits_df %>%filter(`ieu-a-1128` !=0) %>%  select(id.exposure) %>% mutate(id.outcome = "ieu-a-1128")) %>% distinct()
traits_all %>% pull(id.exposure) %>%  unique() %>% length() # 175

# for each pair collect intermediates
all_results <- tibble()
for (i in 1:length(traits_all$id.exposure)){
  
  out <- query_and_tidy_conf(exposure = traits_all$id.exposure[i], 
                             outcome = traits_all$id.outcome[i], 
                             pval_threshold =  0.01) ##### keeping all 
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
  filter(!med_cat %in% c('Diseases', 'Medical Procedures', 'other', "eye_hearing_teeth", "Education", "Psychology", "CHD")) %>% 
  filter(med.id %in% traits_df$id.exposure)   # this keeps only those mediaotrs that are accepted exposure in main analysis validation
  
dim(results_subset)

# next need to validate E->M for mediators and M->E for confounders

#
#
#




dim(results_subset)

#write_csv(results_subset, "01_MR_related/mr_evidence_outputs/conf_med_extracted.csv") # p<10e4
write_csv(results_subset, "01_MR_related/mr_evidence_outputs/conf_med_extracted_all.csv") # p<0.05



protein_names <- read_csv("01_MR_related/mr_evidence_outputs/protein_names.csv", col_names = F) %>% rename(name = X1, gene = X2)
results_subset<- read_csv("01_MR_related/mr_evidence_outputs/conf_med_extracted_all.csv") %>% 
                filter(type %in% c('confounder', 'mediator')) %>%
                create_exposure_categories() %>%
                left_join(protein_names, by = c("med.trait" = 'name')) %>% rename(med.gene = gene) %>% 
                left_join(protein_names, by = c("exposure.trait" = 'name')) %>% rename(exp.gene = gene) %>% 
                select(exposure.trait,exp.gene, med.trait, med.gene, everything()) 
dim(results_subset)# 10308    17
                

## bringing in validated  results for trait-trait
# id.exposure , id.outcome, OR_CI, effect_direction , nsnp, 
validated <- read_tsv("01_MR_related/mr_evidence_outputs/redone_MRconf_subsetoutput_ivw.tsv") %>% 
              select(id.exp_val = id.exposure , id.out_val = id.outcome, r1.OR_CI_val = OR_CI, r1.nsnp_val=nsnp)

# need to do left_join with conf and med separately

mediators <- results_subset %>% filter(type == 'mediator') %>% left_join(validated, by = c("exposure.id"= 'id.exp_val' , 'med.id'= 'id.out_val'))
confounders <- results_subset %>% filter(type == 'confounder') %>% left_join(validated, by = c("exposure.id"= 'id.out_val', 'med.id'= 'id.exp_val'))

results_validated <-bind_rows(mediators, confounders)
      # remove na here


## bringing in validated  results for trait-BC

traits_bc_validated <- read_tsv("01_MR_related/mr_evidence_outputs/redone_MR_subsetoutput_ivw.tsv") %>%  
                      select(id.trait = id.exposure, id.outcome, OR_CI_val = OR_CI, nsnp_val =nsnp) 


validated_with_BC <- results_validated %>% 
                    left_join(traits_bc_validated, by = c("exposure.id" = "id.trait", "outcome.id" = "id.outcome" )) %>% rename('r2.OR_CI_val' = 'OR_CI_val', "r2.nsnp_val"="nsnp_val") %>% 
                    left_join(traits_bc_validated, by = c("med.id" = "id.trait", "outcome.id" = "id.outcome" )) %>% rename('r3.OR_CI_val' = 'OR_CI_val', "r3.nsnp_val"="nsnp_val") %>% 
                    filter(!is.na(r1.OR_CI_val)) %>% 
                    filter(!is.na(r3.OR_CI_val)) %>% 
                    select(-r1.OR_CI, -r2.OR_CI,-r3.OR_CI, -r1.b, -r2.b, -r3.b)
    
dim(validated_with_BC) #4859   23

## add lit size

lit <- read_tsv("02_literature_related/literature_outputs/traits_marked_for_lit_analysis.tsv") %>% select(id, unique_pairs)

validated_with_BC <- validated_with_BC %>%
  left_join(lit, by = c("exposure.id" = "id")) %>% rename('exposure_lit_pairs' = 'unique_pairs') %>% 
  left_join(lit, by = c("med.id" = "id")) %>% rename('med_lit_pairs' = 'unique_pairs')

# rearrange
validated_with_BC <- validated_with_BC %>%
select(type, outcome.id, 
       exposure.trait, exp.gene, med.trait, med.gene, 
       r1.OR_CI_val, r1.nsnp_val, r2.OR_CI_val, r2.nsnp_val, r3.OR_CI_val, r3.nsnp_val,
       exposure_lit_pairs, med_lit_pairs, 
       exposure_cat, exposure.id,med_cat, med.id)



counts <- left_join(
  # no conf per trait
  validated_with_BC %>% filter(type == 'confounder')  %>% 
    select(exposure_cat, exposure.trait,exposure.id, med.id) %>% count(exposure_cat,exposure.trait,exposure.id) %>% rename(conf_count=n),
  
  # no med per trait
  validated_with_BC %>% filter(type == 'mediator')  %>% 
    select(exposure_cat,exposure.trait,exposure.id, med.id) %>% count(exposure_cat,exposure.trait,exposure.id) %>% rename(med_count=n)
)

# in antro
x<- counts %>% filter(exposure_cat == 'Antrophometric')  
mean(x$conf_count)# 38
mean(x$med_count)# 43




## adhoc checking
test <- validated_with_BC %>%  filter(exposure.trait %in% c("Age at menopause (last menstrual period)")) 
test <- validated_with_BC %>%  filter(exposure.trait %in% c("IGF-1")) 
test <- validated_with_BC %>%  filter(exposure.trait %in% c("Intercellular adhesion molecule 1")) 
test3 <- validated_with_BC %>%  filter(exposure.trait %in% c("Ferritin")) 
test <- validated_with_BC %>%  filter(exposure.trait %in% c("Albumin")) 
test2 <- validated_with_BC %>%  filter(exposure.trait %in% c("C-reactive protein")) 
test4 <- validated_with_BC %>%  filter(exposure.trait %in% c("Fibroblast growth factor 7")) 
test %>% select (med.trait, type) %>% distinct() %>% count(type)

mol_med<- test %>% filter(med_cat %in% c("Proteins", "Metabolites")) %>% select(med.trait, gene) %>% distinct()




# e.g. we'refocussing only on those:

exposure_to_extract<- c("met-c-841",
                        "ukb-d-30770_irnt",
                        "prot-a-670",
                        "prot-a-1397",
                        "prot-a-1148",
                        "prot-b-38",
                        "prot-a-1486",
                        "prot-a-366",
                        "prot-a-1097",
                        "prot-a-710")
out <- list()

for (i in exposure_to_extract){
  sub <- validated_with_BC %>%  filter(exposure.id == i)
  out[[i]] <- sub
}
writexl::write_xlsx(out, "01_MR_related/mr_evidence_outputs/med-conf-table.xlsx")



## collelcting all mr-eve involved things

exp_mol <- validated_with_BC %>% filter(exposure.id %in% exposure_to_extract) %>% 
  filter(exposure_cat %in% c("Proteins", "Metabolites")) %>%
  select(trait = exposure.trait, gene=exp.gene) %>% distinct()
med_mol <- validated_with_BC %>%  filter(exposure.id %in% exposure_to_extract) %>% 
  filter(med_cat %in% c("Proteins", "Metabolites")) %>%
  select(trait = med.trait, gene = med.gene) %>% distinct()

all_mol = bind_rows(exp_mol, med_mol) %>% distinct() %>% 
  filter(!grepl("^X-|bonds|groups", trait))

all_mol %>% write_csv("01_MR_related/mr_evidence_outputs/molecular_traits_involved_in_conf_med_subset_upd.csv")

  











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


