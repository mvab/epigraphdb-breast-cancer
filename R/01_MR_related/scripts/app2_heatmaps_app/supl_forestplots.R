library(tidyverse)


inputs<- readRDS("01_MR_related/scripts/app2_heatmaps_app/data/inputs.rds")
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
  rename(OR_CI =OR_CI2) %>% 
  mutate_at(vars(or, lo_ci, up_ci), as.numeric) %>%  
  left_join(passed_pairs) %>% 
  filter(!is.na(mtc)) %>% 
  left_join(merged) %>% 
  filter(!exposure.id %in% antro_blacklist)
  



for (i in unique(res$outcome)){
  
  mr_results <- res %>% filter(outcome == i) %>% arrange(exposure_cat, or)

  
  pal<-(c(unname(yarrr::piratepal("pony"))))
  pal[6:7]<-c('darkgrey', "#FFEA5E")
  p<-ggplot(mr_results, aes(y=exposure_name, x=or, colour=exposure_cat)) +
    geom_errorbarh(aes(xmin=lo_ci, xmax=up_ci), height=.2) +
    geom_point(size=1)+
    scale_color_manual(values=pal)+
    geom_vline(xintercept=1, linetype='longdash') +
    geom_text(aes(label=OR_CI),hjust=-0.1, vjust=-0.6, size =3, color = '#333232')+
    theme_minimal_vgrid(8, rel_small = 1) +
    scale_y_discrete(position = "right")+
    facet_wrap( exposure_cat ~ . , ncol=2,  scales="free")+
    labs(color = "", y = unique(mr_results$outcome), x = "Odds ratio" )+
    theme(legend.position = "none", plot.title.position  = "plot")
  p
  
  ggsave(paste0("01_MR_related/results/heatmap_html/forest_plots/forest_plots_FDR_only_", unique(mr_results$outcome), ".png"),
         plot=p, scale=1, 
         width=30, height=30,
         units=c("cm"), dpi=300, limitsize=F)
  
  
  
  
}




