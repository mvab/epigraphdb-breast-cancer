library(cowplot)
library(plotly)
library(tidyverse)

source("01_MR_related/scripts/app2_heatmaps_app/heatmap_functions.R")
source("01_MR_related/scripts/app2_heatmaps_app/functions_copy_from_mreveapp.R")


 #### load data ----


#inputs <- load_and_merge_heatmap_inputs(extCis = "V3cis", ext="V3") # V3 for new data
#saveRDS(inputs, "01_MR_related/scripts/app2_heatmaps_app/data/inputsV3.rds")

inputs<- readRDS("01_MR_related/scripts/app2_heatmaps_app/data/inputsV3.rds")
merged_input <- inputs$merged
or_ci_data <- inputs$or_ci_data
protein_path_data <- inputs$protein_path_data
antro_blacklist <- inputs$antro_blacklist
passed_pairs <- inputs$passed_pairs


# prep data ----

merged_input <-merged_input %>%  mutate(exposure_cat = ifelse(id.exposure %in% c('prot-a-550', 'prot-a-3000'), "Proteins", exposure_cat))

# give better names
#names <- data_full %>% select(exposure_cat_sub,exposure.id, exposure.trait, exposure) %>% distinct()
#write_csv(names, "01_MR_related/results/mr_evidence_outputs/renaming_key_raw.csv")
names_tidy <- read_csv("01_MR_related/scripts/app2_heatmaps_app/data/renaming_key_tidy.csv") %>% select(exposure.id, exposure)# use new exposure column from here
merged<- merged_input %>%
  select(-exposure) %>%
  left_join(names_tidy, by =c("id.exposure" = "exposure.id")) %>% 
  select(exposure, everything())  %>% 
  filter(!is.na(exposure)) %>% 
  mutate(exposure = gsub("*", "", exposure, fixed = T)) 
length(unique(merged$id.exposure)) # 202

data_full<- prepare_data(merged, protein_path_data, antro_blacklist,or_ci_data, passed_pairs) # antro_blacklist is ingnored

# add data hazards column
sensitivity_hazards <- read_tsv("01_MR_related/results/mr_evidence_outputs/all_data_with_sens_filters_hazardsV3.tsv") 
                     # mutate(hazards = case_when(
                     #   hazards == "H" ~  "H  ",
                     #   hazards == "HP" ~  "HP ",
                     #   hazards == "P" ~  " P " ,
                     #   hazards == "X" ~  "__X",
                     #   TRUE ~ hazards
                     # ))
                
data_full <- data_full %>%  left_join(sensitivity_hazards %>% select(-exposure, -method, -exposure_cat)) 


# drop exposures with no effect
exposures_to_drop <- data_full %>% 
                    group_by(exposure.id, value) %>%
                    count(exposure.id, value) %>% 
                    filter(value == 0 & n == 9) %>% # zero in all 9 outcomes
                    pull(exposure.id)

data_full <- data_full %>% filter(!exposure.id %in% exposures_to_drop)

# add column for sharing
data_full <- data_full %>% mutate(value_mtc = ifelse(!is.na(mtc), paste0(value, "*"), value)) %>% 
                        mutate(value_mtc = factor(value_mtc, levels = c("-1*" ,"-1" , "0" ,  "1" ,  "1*" )))




#### select groups and plot----

font_size =8#8
star_size = 4 # 3

# lifestyle
data_sub<- data_full %>% filter(exposure_cat %in% c("Anthropometric traits", "Lifestyle traits"))
data_sub %>% select(exposure.id, exposure_cat) %>% distinct() %>% count(exposure_cat)
lifestyle <- plot_heatmap4(data_sub,font_size = font_size, star_size = star_size/2, col_order = "normal")


# metabolites 
data_sub<- data_full %>% filter(exposure_cat %in% c("Lipids", "Metabolites"))
metabolites <- plot_heatmap4(data_sub,font_size = font_size, star_size = star_size/2)
metabolites

# proteins  
data_sub<- data_full %>% filter(exposure_cat %in% c("Proteins"))
data_sub <- data_sub %>% 
  arrange( value, outcome, main_path) %>%
  # add cis-used label to protein names
  mutate(used_inst_code = ifelse(is.na(used_inst_code), "", used_inst_code)) %>% 
  mutate(name_mix = paste0(gene," ", used_inst_code))  #### final version - using gene names!


data_sub <- data_sub %>%   
  select(-exposure) %>% rename(exposure= name_mix) %>% 
  mutate(exposure = factor(exposure, levels = unique(data_sub$name_mix))) %>% 
  select(-exposure_cat) %>% rename(exposure_cat=exposure_cat_sub) %>% 
  mutate(exposure_cat = case_when(exposure_cat == "Immune System" ~ "Proteins: Immune System",
                                  exposure_cat == "Metabolism" ~ "Proteins: Metabolism",
                                  exposure_cat == "Signal Transduction" ~ "Proteins: Signal Transduction",
                                  #exposure_cat == "Developmental Biology" ~ "Proteins: Developmental Biology",
                                  exposure_cat == "Other" ~ "Proteins: other",
                                  exposure_cat == "Not mapped" ~ "Proteins: not mapped",
                                  TRUE ~ exposure_cat)) %>% 
  mutate(exposure_cat = factor(exposure_cat, levels = c("Proteins: Immune System", 'Proteins: Metabolism',
                                                        'Proteins: Signal Transduction',#'Proteins: Developmental Biology', 
                                                        "Proteins: other", "Proteins: not mapped")))
  
  
#data_sub %>% left_join(merged %>% select(id.exposure, exposure)) %>%  write_tsv("01_MR_related/mr_evidence_outputs/protein_in_final_set.tsv")
proteins <- plot_heatmap4(data_sub,font_size = font_size, star_size = star_size/2)


leftcol <- plot_grid(lifestyle, metabolites , labels = c('A', 'B'), label_size = 12, ncol=1, rel_heights = c(0.7, 0.6))
rightcol <- plot_grid(proteins, NULL , labels = c('', ''),  ncol=1, rel_heights = c(0.9, 0.1))

full_plot <- plot_grid(leftcol, rightcol , labels = c('', 'C'), label_size = 12, ncol = 2, 
                       rel_widths = c(0.98, 0.79)) +panel_border(remove = T)

#full_plotOLD <- plot_grid(leftcol, proteins , labels = c('', 'C'), label_size = 12, ncol = 2, 
#                                            rel_widths = c(1, 1), rel_heights = c(1, 0.6))

ggsave(paste0("01_MR_related/results/heatmap_html/combined_heatmaps2V3cis_genenames.png"),
       plot=full_plot, scale=1, 
       width=18, height=26,  #for paper 
       #width=18, height=32, # for supplement - complete set: to save with white bg,  add theme(panel.background = element_rect(fill='white')) to full plot
       units=c("cm"), dpi=300, limitsize=F)


# calc % not displayed for each group

diplay_summary <- merged %>%  select(exposure, id.exposure, exposure_cat) %>% distinct() %>% 
  mutate(exposure_cat_broad = case_when(exposure_cat %in% c("Antrophometric") ~ "A - Antro",
                                        exposure_cat %in% c("Physical activity" ,
                                                            "Diet and supplements",
                                                            "Reproductive" ,"Alcohol",
                                                            "Smoking" ,"Other biomarkers" , 
                                                            "Sleep",  "Drugs")   ~ "A - Lifestyle",
                                        exposure_cat %in% c("Lipids") ~ "B - Lipids",
                                        exposure_cat %in% c( "Metabolites") ~ "B - Metabolites",
                                        exposure_cat %in% c("Proteins") ~ "C - Proteins",
                                        TRUE ~ NA
                                        )) %>% 
  mutate(included_in_plot = ifelse(id.exposure %in% exposures_to_drop, F, T)) 

diplay_summary %>%  
  group_by(exposure_cat_broad, included_in_plot) %>% 
  summarise(N = n()) %>%
  mutate(freq = round(N / sum(N), 2)) 

#.  exposure_cat_broad included_in_plot     N  freq
#. <chr>              <lgl>            <int> <dbl>
#1 A - Antro          FALSE                1  0.05
#2 A - Antro          TRUE                19  0.95
#3 A - Lifestyle      FALSE               15  0.43
#4 A - Lifestyle      TRUE                20  0.57
#5 B - Lipids         TRUE                18  1   
#6 B - Metabolites    FALSE                6  0.32
#7 B - Metabolites    TRUE                13  0.68
#8 C - Proteins       FALSE               51  0.46
#9 C - Proteins       TRUE                59  0.54





