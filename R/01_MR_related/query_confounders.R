library(tidyverse)
library(epigraphdb)
source("helper_functions.R")


exposure = "prot-a-1097"
outcome = "ieu-a-1128"
pval_threshold = 1e-04

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


test<-query_epigraphdb_as_table(mediator_query)
xx<-tidy_conf_query_output(test, type = "mediator")




tidy_conf_query_output <- function(df, type){
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
}





res_list <- list(
  colliders = query_epigraphdb_as_table(collider_query),
  confounders = query_epigraphdb_as_table(confounder_query),
  mediators = query_epigraphdb_as_table(mediator_query),
  rev_mediators = query_epigraphdb_as_table(reverse_mediator_query)
)

res_list_tidy <- list(
  colliders = tidy_conf_query_output(res_list$colliders, type = "collider"),
  confounders = tidy_conf_query_output(res_list$confounders, type = "confounder"),
  mediators =  tidy_conf_query_output(res_list$mediators, type = "mediator"),
  rev_mediators = tidy_conf_query_output(res_list$rev_mediators, type = "reverse_mediator")
)

res_list_tidy_df <- bind_rows(res_list_tidy)
res_list_tidy_df %>% select(exposure.trait, med.trait, outcome.trait, r1.b, r2.b, r3.b, type) %>% View()


# just get the mediaotors:
source("explore_MR-EvE_app/functions.R")
tmp <- res_list_tidy_df %>% #filter(type == 'mediator') %>% 
                  select(exposure.trait = med.trait, 
                         exposure.id = med.id, r3.OR_CI, r3.effect_direction, type) %>% 
                 create_exposure_categories() %>%
                 filter(!exposure_cat %in% c('Diseases', 'Medical Procedures', 'other', "Education", "Psychology", "CHD")) 

tmp %>% select(type, exposure_cat) %>% group_by(type, exposure_cat) %>% count(type, exposure_cat)


in_final <- read_tsv("mr_evidence_outputs/trait_manual_ivw_subtypes_merged.tsv") %>% 
  filter(exposure_cat %in% c("Proteins", "Metabolites", "Other biomarkers")) %>% 
  filter(!grepl("LDL|HDL|cholest|trigl", exposure_name, ignore.case = T)) %>% 
  pull(id.exposure)


tmp %>% mutate(also = ifelse(exposure.id %in% in_final, 1, 0)) %>% View()




query = paste0("
    MATCH (gene1:Gene)-[gene_to_protein:GENE_TO_PROTEIN]->(protein:Protein)-[protein_in_pathway:PROTEIN_IN_PATHWAY]->(pathway:Pathway) 
    WHERE gene1.name in ['MAN1A2', 'FGF7', 'CPXM1', 'IL1RL1', 'TPST2','SULF2','FLNA','DCBLD2','PRDM1','ISLR2']
    RETURN gene1.name, protein.uniprot_id,  pathway.name
    order by pathway.name
      ")


test<-query_epigraphdb_as_table(query)
test %>% count(pathway.name, sort=T)


