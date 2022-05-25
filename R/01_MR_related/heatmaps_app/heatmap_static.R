library(cowplot)
library(plotly)
library(tidyverse)

source("heatmap_functions.R")


#### load data ----


inputs <- load_and_merge_heatmap_inputs()

saveRDS(inputs, "01_MR_related/heatmaps_app/data/inputs.rds")

readRDS("01_MR_related/heatmaps_app/data/inputs.rds")
merged <- inputs$merged
or_ci_data <- inputs$or_ci_data
protein_path_data <- inputs$protein_path_data
antro_blacklist <- inputs$antro_blacklist
passed_pairs <- inputs$passed_pairs

# prep data ----


data_full<- prepare_data(merged, protein_path_data, antro_blacklist,or_ci_data, passed_pairs)

#### select groups----

if (panel == "lifestyle"){
  data_full<- data_full %>% filter(exposure_cat %in% c("Antrophometric traits", "Lifestyle traits"))
  
}  else if (panel == 'metabolites'){
  
  data_full<- data_full %>% filter(exposure_cat %in% c("Lipids", "Metabolites"))
  
} else if (panel == 'proteins'){
  
  data_full<- data_full %>% filter(exposure_cat %in% c("Proteins"))
}



# plot! ----

# for proteins: with subcats:
data_full <- data_full %>% 
  select(-exposure_cat) %>% rename(exposure_cat=exposure_cat_sub) %>% 
  mutate(exposure_cat_sub = factor(exposure_cat, levels = c("Immune System", 'Metabolism', 'Signal Transduction','Developmental Biology', "Other", "Not mapped"))) %>% 
  select(-exposure) %>% rename(exposure= name_mix)


p <- plot_heatmap(data_full)
p

ggplotly(p)




data_sub %>% left_join(merged %>% select(id.exposure, exposure)) %>%  write_tsv("01_MR_related/mr_evidence_outputs/protein_in_final_set.tsv")






#### for static figure in the paper

lifestyle <- p
proteins <- p
metabolites <- p


leftcol <- plot_grid(lifestyle, metabolites , labels = c('A', 'B'), label_size = 12, ncol=1)
full_plot <- plot_grid(leftcol, proteins , labels = c('', 'C'), label_size = 12, ncol = 2, rel_widths = c(1, 0.9))
full_plot


ggsave(paste0("01_MR_related/heatmap_html/combined_heatmaps.png"),
       plot=full_plot, scale=1, 
       width=17, height=26,
       units=c("cm"), dpi=300, limitsize=F)









# legacy ----
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
