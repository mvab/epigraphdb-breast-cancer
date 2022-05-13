library(tidyverse)
library(cowplot)
library(plotly)
source("01_MR_related/explore_MR-EvE_app/functions.R")


redone_MR <- read_tsv("01_MR_related/mr_evidence_outputs/redone_MR_fulloutput.tsv") %>% 
              filter(method == "Inverse variance weighted") %>% 
              select(id.exposure, id.outcome, OR_CI, pval, exposure_name, nsnp) %>% 
              mutate(outcome = case_when(id.outcome =="ieu-a-1126" ~ "BCAC 2017",
                                         id.outcome =="ieu-a-1127"    ~ "ER+",   
                                         id.outcome =="ieu-a-1128"   ~ "ER-" )) %>% select(-id.outcome)

res_subtypes <- read_tsv("01_MR_related/mr_evidence_outputs/mr_subtypes/all_traits_MR_vs_BCAC2020.tsv") %>% 
                filter(method == "Inverse variance weighted") %>% 
                select(id.exposure, outcome, OR_CI, pval,  nsnp) %>% 
                mutate(outcome = case_when(outcome =="Breast cancer BCAC 2020" ~ "BCAC 2020",
                                           outcome =="LuminalA ER+PR+HER-"    ~ "Luminal A",   
                                           outcome =="LuminalB1 ER+PR+HER-"   ~ "Luminal B1",  
                                           outcome =="LuminalB2 ER+PR+HER+" ~ "Luminal B2" ,
                                           outcome =="HER2-enriched ER-PR-HER+"  ~ "HER2-enriched" ,
                                           outcome =="TNBC ER-PR-HER-"  ~ "TNBC" ,
                                           outcome =="TNBC_BRCA1 ER-PR-HER-" ~ "TNBC_BRCA1" )) %>% 
              left_join(redone_MR %>% select(id.exposure, exposure_name ))

or_ci_data <- bind_rows(redone_MR,res_subtypes ) %>% rename('exposure.id'= 'id.exposure') 


merged <- read_tsv("01_MR_related/mr_evidence_outputs/trait_manual_ivw_subtypes_merged.tsv")
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
  

protein_path_data <- read_tsv("01_MR_related/pathways/proteins_w_pathways.tsv") %>% 
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
                     'ukb-b-12039', 'ukb-b-2303',
                     'ieu-a-2', 'ieu-a-835',
                     'ukb-b-18105', 'ieu-a-99', 'ieu-a-105', 'ieu-a-109', 'ukb-b-4650', 'ieu-a-95', 'ukb-a-248', 'ieu-a-68', 	
                     'ieu-a-101')





####### SELECTIONS 


## all lifestyle
panel = "lifestyle"
data_sub <- merged %>% 
  filter(!exposure_cat %in% c('Proteins', "Metabolites","Lipids",  "Other biomarkers")) %>%
  select(id.exposure, contains('BCAC'), "ER+", contains("Luminal"), "ER-", "HER2-enriched","TNBC" ) %>% 
  filter(!(`BCAC 2017` == 0 & `BCAC 2020` == 0 & `ER-` == 0  & `ER+` == 0 & `Luminal A` == 0 & `Luminal B1` == 0 & `Luminal B2` == 0 &`HER2-enriched` == 0 & `TNBC` == 0 ))

## all:  proteins ---- extract for pathway analysis
panel = "proteins"
data_sub <- merged %>% 
  filter(exposure_cat %in% c('Proteins', "Other biomarkers")) %>%
  filter(!grepl("LDL|HDL|cholest|trigl", exposure, ignore.case = T)) %>% 
  select(id.exposure, contains('BCAC'), "ER+", contains("Luminal"), "ER-", "HER2-enriched","TNBC" ) %>% 
  filter(!(`BCAC 2017` == 0 & `BCAC 2020` == 0 & `ER-` == 0  & `ER+` == 0 & `Luminal A` == 0 & `Luminal B1` == 0 & `Luminal B2` == 0 &`HER2-enriched` == 0 & `TNBC` == 0 ))
data_sub %>% left_join(merged %>% select(id.exposure, exposure)) %>%  write_tsv("01_MR_related/mr_evidence_outputs/protein_in_final_set.tsv")

## all: metabolites and lipids
panel = 'metabolites'
data_sub <- merged %>% 
  filter(exposure_cat %in% c( "Metabolites","Other biomarkers", "Lipids")) %>%
  filter(!grepl("Average|ratio", exposure, ignore.case = T)) %>% 
  select(id.exposure, contains('BCAC'), "ER+", contains("Luminal"), "ER-", "HER2-enriched","TNBC" ) %>% 
  filter(!(`BCAC 2017` == 0 & `BCAC 2020` == 0 & `ER-` == 0  & `ER+` == 0 & `Luminal A` == 0 & `Luminal B1` == 0 & `Luminal B2` == 0 &`HER2-enriched` == 0 & `TNBC` == 0 ))



##### THE REST

data_sub2<- reshape2::melt(data_sub, id.var = 'id.exposure') %>% 
  mutate(value=as.character(value)) %>% 
  arrange( value) %>% 
  left_join(merged %>% select(id.exposure, exposure, exposure_cat, exposure_name)) %>% 
  rename('exposure.trait' = 'exposure','exposure.id'= 'id.exposure') %>% 
  create_exposure_categories() 

if (panel == "lifestyle"){
  data_sub2<- data_sub2 %>% 
    mutate(exposure_cat = ifelse(exposure_cat == 'Antrophometric', 'Antrophometric traits', 'Lifestyle traits')) %>% 
    filter(!exposure.id %in%   antro_blacklist)%>% 
    mutate(exposure = ifelse(grepl("(F)", exposure_name, fixed = T), paste0(exposure, " (F)"), exposure)) %>% 
    mutate(exposure = ifelse(grepl("AdjBMI", exposure_name, fixed = T), paste0(exposure, " - AdjBMI"), exposure))
  
}  else if (panel == 'proteins'){
  data_sub2<- data_sub2 %>% 
    left_join(protein_path_data, by = c("exposure.id"="id.exposure", "exposure" = "name")) %>% 
    mutate(gene = factor(gene, levels = unique(protein_path_data$gene))) %>% 
    filter(!exposure.id %in% c("prot-a-2396", 'prot-a-1540')) %>%  # duplicates
    mutate(exposure_cat = ifelse(exposure == 'Albumin', 'Proteins', exposure_cat)) %>% 
    select(-exposure_cat) %>% rename(exposure_cat= main_path) %>% 
    mutate(exposure_cat = factor(exposure_cat, levels = c("Immune System", 'Metabolism', 'Signal Transduction','Developmental Biology', "Other", "Not mapped"))) %>% 
    select(-exposure) %>% rename(exposure= name_mix)
  
} else if (panel == 'metabolites'){
  data_sub2<- data_sub2 %>% 
    mutate(exposure = ifelse(exposure == 'X-11445--5-alpha-pregnan-3beta,20alpha-disulfate', 'X-11445', exposure)) %>% 
    mutate(exposure = gsub('chylomicrons', 'ULDLs', exposure))
}

data_sub3 <- data_sub2 %>%
  mutate(exposure = factor(exposure, levels = unique(data_sub2$exposure))) %>% 
  rename(outcome = variable) %>%
  mutate(outcome=as.character(outcome)) %>% 
  mutate(value = ifelse(is.na(value), 0 , value)) %>% 
  mutate(effect_direction = ifelse(value == "1" , 'positive',
                             ifelse(value == "-1", 'negative', 'overlaps null'))) %>% 
  mutate_at(vars('outcome'), funs(case_when(. == "BCAC 2017" ~ "BCAC'17", 
                                           . == "BCAC 2020" ~ "BCAC'20", 
                                           . == "Luminal A" ~ "Lum A", 
                                           . == "Luminal B1" ~ "Lum B1", 
                                           . == "Luminal B2" ~ "Lum B2", 
                                           . == "HER2-enriched" ~ "HER2", TRUE ~ .))) %>% 
  mutate(outcome = factor(outcome, 
                          levels = c("BCAC'17", "BCAC'20","ER+", "Lum A","Lum B1", 
                                     "Lum B2", "ER-", "HER2", "TNBC" ))) %>% 
  left_join(or_ci_data) %>% distinct() 


p<-ggplot(data_sub3, aes( y=exposure, x=outcome, fill = effect_direction, 
                          id = exposure.id ,name = exposure_name,
                          or_ci= OR_CI, pval =  pval, nsnp=nsnp, cat = exposure_cat)) + 
  geom_tile(colour = "grey") + 
  scale_fill_manual(values=c("#4D9221","white", "#C51B7D"))+ ### normal
  theme_minimal_grid(8) +
  panel_border() +
  labs(fill = "Effect direction", x="", y ="")+
  theme(axis.text.y = element_text( size=6), #7
        axis.text.x = element_text(angle=40, size=6, hjust = 1), legend.position = 'none') #  40 , 6

p<- p + ggforce::facet_col(vars(exposure_cat), scales = "free_y", space = "free") 

p

ggplotly(p)


lifestyle <- p
proteins <- p
metabolites <- p


library(cowplot)

leftcol <- plot_grid(lifestyle, metabolites , labels = c('A', 'B'), label_size = 12, ncol=1)
full_plot <- plot_grid(leftcol, proteins , labels = c('', 'C'), label_size = 12, ncol = 2, rel_widths = c(1, 0.9))
full_plot


ggsave(paste0("01_MR_related/heatmap_html/combined_heatmaps.png"),
       plot=full_plot, scale=1, 
       width=17, height=26,
       units=c("cm"), dpi=300, limitsize=F)





##### 
#tmp subsets

to_keep <- c('ukb-d-30770_irnt', 'prot-a-670', 'prot-a-1397', 'prot-a-1097', 'prot-a-831', 'prot-a-1486', 	'prot-a-1148')
to_keep <- c('prot-a-710', 'prot-a-2629', 'prot-a-366', 'prot-a-2395', 'prot-a-1117', 'prot-b-38', 'prot-a-2892', 	'ukb-a-132')
to_keep <- c( "met-c-841", "prot-a-1930", "prot-a-3076", "ieu-a-302", "ukb-d-30770_irnt", "ieu-a-1049", "prot-a-670", "met-a-355", "ukb-a-132",
              "prot-a-3193", "prot-a-1397", "prot-a-655", "prot-a-2720", "ieu-a-1", "prot-a-1148", "met-a-316", "prot-a-67", "prot-b-38", 
              "prot-a-387", "ukb-a-142", "prot-b-71", "prot-a-1486", "prot-a-1542", "prot-a-1166", "prot-a-1313", "prot-a-2395", "prot-a-2007", 
              "prot-b-55", "prot-a-394", "prot-a-2889", "prot-a-366", "prot-a-1097")

to_keep <- c("ieu-a-1096",   "ukb-b-4424")#gr2
to_keep <- c('ieu-a-974', "ukb-b-17422")#g1
to_keep <- c("prot-a-710", "prot-a-366", "prot-a-1486", "prot-b-38")#gr3
to_keep <- c('ukb-d-30770_irnt', 'prot-a-670',  'prot-a-1097', 'prot-a-1397')

## custom
data_sub <- merged %>% 
  filter(id.exposure %in% to_keep) %>%
  filter(!grepl("LDL|HDL|cholest|trigl", exposure, ignore.case = T)) %>% 
  select(id.exposure, contains('BCAC'), "ER+", contains("Luminal"), "ER-", "HER2-enriched","TNBC" ) %>% 
  filter(!(`BCAC 2017` == 0 & `BCAC 2020` == 0 & `ER-` == 0  & `ER+` == 0 & `Luminal A` == 0 & `Luminal B1` == 0 & `Luminal B2` == 0 &`HER2-enriched` == 0 & `TNBC` == 0 ))


## effect in many

lit <-read_tsv("02_literature_related/literature_outputs/traits_marked_for_lit_analysis.tsv") %>% 
  filter(unique_triples > 50) %>%  pull(id)

data_sub <- merged %>%
  filter(!grepl("LDL|HDL|cholest|trigl", exposure, ignore.case = T)) %>% 
  filter(exposure_cat != 'Antrophometric') %>% 
  filter(id.exposure %in% lit) %>% 
  select(id.exposure, `BCAC 2017`:`TNBC`) %>% 
  rowwise() %>%
  mutate(`N_zeros` = sum(c_across(`BCAC 2017`:`TNBC`) == 0)) %>% 
  filter(N_zeros >= 1 & N_zeros <8 )


#####


####

merged %>% select(exposure, id.exposure, exposure_cat) %>% distinct() %>%
  arrange(exposure_cat, exposure) %>% write_tsv("01_MR_related/mr_evidence_outputs/trait_names_for_grouping.tsv")


dat <-  read_tsv("01_MR_related/mr_evidence_outputs/tidy_traits_by_cat.tsv") %>%  # made in explore_mr_results.Rmd
  # update exposure categories 
  create_exposure_categories() %>%  filter(exposure_cat != 'other')

dat %>% filter(exposure_cat == "Alcohol") %>% select(exposure.trait, exposure.id) %>% distinct() %>% View()
