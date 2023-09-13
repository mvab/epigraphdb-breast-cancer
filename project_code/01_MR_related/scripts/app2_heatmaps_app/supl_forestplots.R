library(tidyverse)
library(cowplot)


inputs<- readRDS("01_MR_related/scripts/app2_heatmaps_app/data/inputsV3.rds")
or_ci_data <- inputs$or_ci_data
passed_pairs <- inputs$passed_pairs
merged <- inputs$merged
antro_blacklist <- inputs$antro_blacklist

merged <- merged %>% 
  select(exposure.id = id.exposure, exposure_cat) %>% 
  mutate(exposure_cat = case_when(exposure_cat == 'Antrophometric' ~ 'Anthropometric traits', 
                                  exposure_cat == 'Metabolites' ~ 'Metabolites', 
                                  exposure_cat == 'Proteins' ~ 'Proteins', 
                                  exposure_cat == 'Lipids' ~ 'Metabolites', 
                                  TRUE ~ 'Lifestyle traits')) 



res<- or_ci_data %>% 
  mutate_at(vars('outcome'), funs(case_when(. == "BCAC 2017" ~ "BCAC'17", 
                                            . == "BCAC 2020" ~ "BCAC'20", 
                                            . == "Luminal A" ~ "Lum A", 
                                            . == "Luminal B1" ~ "Lum B1", 
                                            . == "Luminal B2" ~ "Lum B2", 
                                            . == "HER2-enriched" ~ "HER2", TRUE ~ .))) %>% 
  mutate(OR_CI2 =OR_CI) %>% 
  mutate(OR_CI = gsub(":", " ", OR_CI)) %>% 
  mutate(OR_CI = gsub("[", "", OR_CI, fixed = T)) %>% 
  mutate(OR_CI = gsub("]", "", OR_CI, fixed = T)) %>% 
  separate(OR_CI, into = c("or", "lo_ci", "up_ci"), sep=' ') %>% 
 # rename(OR_CI =OR_CI2) %>% 
  mutate_at(vars(or, lo_ci, up_ci), as.numeric) %>%  
  mutate(OR_CI = paste0(round(or,2), " [", round(lo_ci,2), ":", round(up_ci,2), "]")) %>% 
  left_join(passed_pairs) %>% 
  #filter(!is.na(mtc)) %>% 
  left_join(merged) %>% 
  mutate(outcome = factor(outcome, levels = rev(c("BCAC'17", "BCAC'20" ,
                                                  "ER+"  ,   "ER-" ,  
                                                  "Lum A" ,  "Lum B1" , "Lum B2" , "HER2"  ,  "TNBC"  ))))

# only do this for selected case studies
keep <- c("ukb-b-4650", "ukb-a-11", "ukb-d-30760_irnt", "met-a-500")

res_sub <- res %>% filter(exposure.id %in% keep)



for (i in unique(res_sub$exposure.id)){
  
  mr_results <- res_sub %>% filter(exposure.id == i) 

  pal<- c("#EB5291FF", "#FBBB68FF" ,"#F5BACFFF", "#9DDAF5FF", "#6351A0FF" ,"#FFEA5E", "#1794CEFF",'darkgrey', "#972C8DFF")
  
  p<-ggplot(mr_results, aes(y=outcome, x=or, colour=outcome)) +
    geom_errorbarh(aes(xmin=lo_ci, xmax=up_ci), height=.3) +
    geom_point(size=2)+
    scale_color_manual(values=pal)+
    geom_vline(xintercept=1, linetype='longdash') +
    geom_text(aes(label=OR_CI),hjust=-0.2, vjust=-0.6, size =3, color = '#333232')+
    theme_minimal_vgrid(8, rel_small = 1) +
    scale_y_discrete(position = "right")+
    scale_x_log10()+
    labs(color = "", y = "Outcomes", x = "Odds ratio", title = unique(mr_results$exposure_name) )+
    theme(legend.position = "none",plot.title = element_text(hjust = 0.5))
  
  ggsave(paste0("01_MR_related/results/heatmap_html/forest_plotsV3/forest_plots_FDR_only_", unique(mr_results$exposure.id), ".png"),
         plot=p, scale=1.2, 
         width=10, height=6,
         units=c("cm"), dpi=300, limitsize=F)
  
}




