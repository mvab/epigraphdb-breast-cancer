library(cowplot)
library(plotly)
library(tidyverse)


#### functions ----
load_and_merge_heatmap_inputs <- function() {
  
  redone_MR <- read_tsv("01_MR_related/results/mr_evidence_outputs/redone_MR_fulloutput.tsv") %>% 
    filter(method %in% c("Inverse variance weighted", "Wald ratio")) %>% 
    select(id.exposure, id.outcome, OR_CI, pval, exposure_name, nsnp) %>% 
    mutate(outcome = case_when(id.outcome =="ieu-a-1126" ~ "BCAC 2017",
                               id.outcome =="ieu-a-1127"    ~ "ER+",   
                               id.outcome =="ieu-a-1128"   ~ "ER-" )) %>% select(-id.outcome)
  
  res_subtypes <- read_tsv("01_MR_related/results/mr_evidence_outputs/mr_subtypes/all_traits_MR_vs_BCAC2020.tsv") %>% 
    filter(method %in% c("Inverse variance weighted", "Wald ratio")) %>% 
    select(id.exposure, outcome, OR_CI, pval,  nsnp) %>% 
    mutate(outcome = case_when(outcome =="Breast cancer BCAC 2020" ~ "BCAC 2020",
                               outcome =="LuminalA ER+PR+HER-"    ~ "Luminal A",   
                               outcome =="LuminalB1 ER+PR+HER-"   ~ "Luminal B1",  
                               outcome =="LuminalB2 ER+PR+HER+" ~ "Luminal B2" ,
                               outcome =="HER2-enriched ER-PR-HER+"  ~ "HER2-enriched" ,
                               outcome =="TNBC ER-PR-HER-"  ~ "TNBC" ,
                               outcome =="TNBC_BRCA1 ER-PR-HER-" ~ "TNBC_BRCA1" ))%>% 
    left_join(redone_MR %>% select(id.exposure, exposure_name )) %>% distinct() 
  
  or_ci_data <- bind_rows(redone_MR,res_subtypes ) %>% rename('exposure.id'= 'id.exposure') 
  
  
  merged <- read_tsv("01_MR_related/results/mr_evidence_outputs/trait_manual_ivw_subtypes_merged.tsv")
  merged <- merged %>% rename("BCAC 2017" = "ieu-a-1126",
                              "ER+" = "ieu-a-1127",     
                              "ER-" = "ieu-a-1128",     
                              "BCAC 2020" = "Breast cancer BCAC 2020",
                              "Luminal A"=  "LuminalA ER+PR+HER-",     
                              "Luminal B1" = "LuminalB1 ER+PR+HER-" ,  
                              "Luminal B2" = "LuminalB2 ER+PR+HER+", 
                              "HER2-enriched" = "HER2-enriched ER-PR-HER+" ,
                              "TNBC" = "TNBC ER-PR-HER-" ,
                              "TNBC_BRCA1"  ="TNBC_BRCA1 ER-PR-HER-") 
  
  
  merged <- merged %>% mutate(exposure_cat = ifelse(exposure == 'Albumin', 'Proteins', exposure_cat))
  
  
  
  
  protein_path_data <- read_tsv("external_files/pathways/proteins_w_pathways.tsv") %>% 
    mutate(name_len = nchar(name)) %>% 
    mutate(name_mix = ifelse(name_len > 36, gene, name)) %>% 
    mutate(name_mix = ifelse(name_mix %in% c("Dual specificity protein kinase CLK2",
                                             "E3 ubiquitin-protein ligase RNF8",
                                             "Extracellular sulfatase Sulf-2"),
                             gene, name_mix)) %>% 
    mutate(name_mix = ifelse(grepl("Vascular endothelial growth factor receptor", name),
                             gsub("Vascular endothelial growth factor receptor", "VEGFR-", name),name_mix))
  
  antro_blacklist <- c('ieu-a-81','ieu-a-74', "ieu-a-73" ,"ieu-a-79" ,"ieu-a-72" ,"ieu-a-78",
                       'ieu-a-63',  'ieu-a-66', 'ieu-a-60',  'ieu-a-69' , 'ieu-a-61',
                       'ieu-a-54', 'ieu-a-55',  'ieu-a-49' , 'ieu-a-48', 'ieu-a-57' , 'ieu-a-50',
                       'ukb-b-12039', 'ukb-b-2303', 'ieu-a-2', 'ieu-a-835',
                       'ukb-b-18105', 'ieu-a-99', 'ieu-a-105', 'ieu-a-109', 'ukb-b-4650', 'ieu-a-95',
                       'ukb-a-248', 'ieu-a-68', 	
                       'ieu-a-101')
  
  # MR pairs that passed MTC
  passed_pairs <- read_tsv("01_MR_related/results/mr_evidence_outputs/passed_multiple_testing.tsv") %>% 
    mutate(outcome = case_when(outcome =="Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis)"  ~ "BCAC'17",
                               outcome =="ER+ Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis)"    ~ "ER+",   
                               outcome =="ER- Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis)"   ~ "ER-",  
                               outcome =="Breast cancer BCAC 2020" ~ "BCAC'20",
                               outcome =="LuminalA ER+PR+HER-"    ~ "Lum A",   
                               outcome =="LuminalB1 ER+PR+HER-"   ~ "Lum B1",  
                               outcome =="LuminalB2 ER+PR+HER+" ~ "Lum B2" ,
                               outcome =="HER2-enriched ER-PR-HER+"  ~ "HER2" ,
                               outcome =="TNBC ER-PR-HER-"  ~ "TNBC" ,
                               outcome =="TNBC_BRCA1 ER-PR-HER-" ~ "TNBC_BRCA1" )) %>% 
    mutate(mtc = ifelse(mtc == "*", "\n*", mtc)) %>% 
    select(exposure.id = id.exposure, outcome, mtc, qval)
  
  return(list(merged = merged,
              or_ci_data = or_ci_data,
              protein_path_data = protein_path_data,
              antro_blacklist = antro_blacklist, 
              passed_pairs = passed_pairs))
}


reshape_heatmap_matrix <- function(data_sub, merged){
  data_sub<- reshape2::melt(data_sub, id.var = 'id.exposure') %>% 
    mutate(value=as.character(value)) %>% 
    arrange( value) %>% 
    left_join(merged %>% select(id.exposure, exposure, exposure_cat, exposure_name)) %>% 
    rename('exposure.trait' = 'exposure','exposure.id'= 'id.exposure') %>% 
    create_exposure_categories() 
}



prepare_data <- function(merged, protein_path_data, antro_blacklist, or_ci_data, passed_pairs){
  
  data_sub <- merged %>% 
    
    filter(!id.exposure %in% antro_blacklist)%>% 
    mutate(exposure = ifelse(grepl("(F)", exposure_name, fixed = T), paste0(exposure, " (F)"), exposure)) %>% 
    mutate(exposure = ifelse(grepl("AdjBMI", exposure_name, fixed = T), paste0(exposure, " - AdjBMI"), exposure)) %>% 
    filter(!id.exposure %in% c("prot-a-2396", 'prot-a-1540')) %>%  # duplicates
    filter(!grepl("Average number|ratio ", exposure, ignore.case = T)) %>% 
    mutate(exposure = gsub('chylomicrons', 'ULDLs', exposure)) %>% 
    select(id.exposure, contains('BCAC'), "ER+", contains("Luminal"), "ER-", "HER2-enriched","TNBC" ) %>% 
    #filter(!(`BCAC 2017` == 0 & `BCAC 2020` == 0 & `ER-` == 0  & `ER+` == 0 & `Luminal A` == 0 & `Luminal B1` == 0 & `Luminal B2` == 0 &`HER2-enriched` == 0 & `TNBC` == 0 ))
    reshape_heatmap_matrix(., merged) %>% 
    mutate(exposure_cat = ifelse(exposure == 'Albumin', 'Proteins', exposure_cat)) %>% 
    mutate(exposure = ifelse(exposure == 'X-11445--5-alpha-pregnan-3beta,20alpha-disulfate', 'X-11445', exposure)) %>% 
    mutate(exposure_cat_sub = exposure_cat) %>% # backup
    mutate(exposure_cat = case_when(exposure_cat == 'Antrophometric' ~ 'Antrophometric traits', 
                                    exposure_cat == 'Metabolites' ~ 'Metabolites', 
                                    exposure_cat == 'Proteins' ~ 'Proteins', 
                                    exposure_cat == 'Lipids' ~ 'Lipids', 
                                    TRUE ~ 'Lifestyle traits')) %>%  # should be after reshaping
    mutate(exposure_cat_sub = replace(exposure_cat_sub, exposure_cat_sub ==  'Proteins', NA)) %>%  # will add subgroups from path data
    left_join(protein_path_data, by = c("exposure.id"="id.exposure")) %>% 
    mutate(gene = factor(gene, levels = unique(protein_path_data$gene))) %>% 
    mutate(exposure_cat_sub = coalesce(exposure_cat_sub, main_path)) 
  
  
  data_tidy <- data_sub %>%
    mutate(exposure = factor(exposure, levels = unique(data_sub$exposure))) %>% 
    rename(outcome = variable) %>%
    mutate(outcome=as.character(outcome)) %>% 
    mutate(value = ifelse(is.na(value), 0 , value)) %>% 
    mutate(effect_direction = ifelse(value == "1" , 'positive',
                                     ifelse(value == "-1", 'negative', 'overlaps null'))) %>% 
    left_join(or_ci_data) %>% distinct() %>% 
    mutate_at(vars('outcome'), funs(case_when(. == "BCAC 2017" ~ "BCAC'17", 
                                              . == "BCAC 2020" ~ "BCAC'20", 
                                              . == "Luminal A" ~ "Lum A", 
                                              . == "Luminal B1" ~ "Lum B1", 
                                              . == "Luminal B2" ~ "Lum B2", 
                                              . == "HER2-enriched" ~ "HER2", TRUE ~ .))) %>% 
    mutate(outcome = factor(outcome, 
                            levels = c("BCAC'17", "BCAC'20","ER+", "Lum A","Lum B1", 
                                       "Lum B2", "ER-", "HER2", "TNBC" ))) %>% 
    left_join(passed_pairs) %>% 
    separate(exposure_name, into = c("tmp", "exposure_details"), sep = "\\(") %>% 
    mutate(exposure_details = ifelse(!grepl("^F|^M", exposure_details), "NA", exposure_details)) %>% 
    mutate(exposure_details = gsub(")", "; ",exposure_details, fixed = T)) %>% 
    mutate(empty_col = "")
  
  print(dim(data_tidy))
  return(data_tidy)
  
  
}


plot_heatmap <- function(data_tidy, font_size = 11, star_size = 7){
  
  font_size_ax = font_size - 2
  
  p<-ggplot(data_tidy, aes( y=exposure, x=outcome, fill = value, 
                            id = exposure.id ,
                            or_ci= OR_CI, pval =  pval, nsnp=nsnp, cat = exposure_cat,
                            text = paste(" ", empty_col,
                                         '</br>Exposure ID: ', exposure.id,
                                         '</br>Exposure: ', exposure,
                                         '</br>Exposure details: ', exposure_details,
                                         '</br>Outcome: ',  outcome,
                                         '</br>P-value: ', format(pval),
                                         '</br>P-value (FDR adjusted): ', format(qval),
                                         '</br>Odds ratio: ', OR_CI,
                                         '</br>nSNPs: ', nsnp)
  )) + 
    geom_tile(colour = "grey") + 
    geom_text(aes(label=mtc,  vjust = 0.27), size = star_size) +
    scale_fill_manual(values=c("#4D9221","white", "#C51B7D"))+ ### normal
    theme_minimal_grid(font_size) +
    panel_border() +
    labs(fill = "Effect direction", x="", y ="")+
    theme(axis.text.y = element_text( size=font_size_ax), #7
          axis.text.x = element_text(angle=40, size=font_size_ax, hjust = 1), legend.position = 'none') #  40 , 6
  
  p<- p + ggforce::facet_col(vars(exposure_cat), scales = "free_y", space = "free") 
  
  return(p)
  
}

plot_heatmap2 <- function(data_tidy, font_size = 11, star_size = 7, col_order = "normal"){
  
  font_size_ax = font_size - 2
  
  if (col_order == 'normal'){
    pal_values <- c( "#4D9221","#88C254", "white","#DE7BB2", "#C51B7D") # g, light g, w, light p, p 
  } else if (col_order == "weird"){
    pal_values <- c( "#88C254", "#4D9221", "white","#DE7BB2", "#C51B7D") #light g, g, w, light p, p 
  }

  
  p<-ggplot(data_tidy, aes( y=exposure, x=outcome, fill = value_mtc, 
                            id = exposure.id ,
                            or_ci= OR_CI, pval =  pval, nsnp=nsnp, cat = exposure_cat,
                            text = paste(" ", empty_col,
                                         '</br>Exposure ID: ', exposure.id,
                                         '</br>Exposure: ', exposure,
                                         '</br>Exposure details: ', exposure_details,
                                         '</br>Outcome: ',  outcome,
                                         '</br>P-value: ', format(pval),
                                         '</br>P-value (FDR adjusted): ', format(qval),
                                         '</br>Odds ratio: ', OR_CI,
                                         '</br>nSNPs: ', nsnp)
  )) + 
    geom_tile(colour = "grey") + 
    geom_text(aes(label=mtc,  vjust = 0.27), size = star_size) +
    scale_fill_manual(values = pal_values)+ 
    
    theme_minimal_grid(font_size) +
    panel_border() +
    labs(fill = "Effect direction", x="", y ="")+
    theme(axis.text.y = element_text( size=font_size_ax), #7
          axis.text.x = element_text(angle=40, size=font_size_ax, hjust = 1), legend.position = 'none') #  40 , 6
  
  p<- p + ggforce::facet_col(vars(exposure_cat), scales = "free_y", space = "free") 
  
  return(p)
  
}



