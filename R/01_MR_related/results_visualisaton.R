library(plotly)

redone_MR <- read_tsv("01_MR_related/mr_evidence_outputs/redone_MR_fulloutput.tsv") %>% 
              filter(method == "Inverse variance weighted") %>% 
              select(id.exposure, id.outcome, OR_CI, pval) %>% 
              mutate(outcome = case_when(id.outcome =="ieu-a-1126" ~ "BCAC 2017",
                                         id.outcome =="ieu-a-1127"    ~ "ER+",   
                                         id.outcome =="ieu-a-1128"   ~ "ER-" )) %>% select(-id.outcome)

res_subtypes <- read_tsv("01_MR_related/mr_evidence_outputs/mr_subtypes/all_traits_MR_vs_BCAC2020.tsv") %>% 
                filter(method == "Inverse variance weighted") %>% 
                select(id.exposure, outcome, OR_CI, pval) %>% 
                mutate(outcome = case_when(outcome =="Breast cancer BCAC 2020" ~ "BCAC 2020",
                                           outcome =="LuminalA ER+PR+HER-"    ~ "Luminal A",   
                                           outcome =="LuminalB1 ER+PR+HER-"   ~ "Luminal B1",  
                                           outcome =="LuminalB2 ER+PR+HER+" ~ "Luminal B2" ,
                                           outcome =="HER2-enriched ER-PR-HER+"  ~ "HER2-enriched" ,
                                           outcome =="TNBC ER-PR-HER-"  ~ "TNBC" ,
                                           outcome =="TNBC_BRCA1 ER-PR-HER-" ~ "TNBC_BRCA1" ))

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



## BCAC lifestyle
data_sub <- merged %>% 
  filter(!exposure_cat %in% c('Proteins', "Metabolites", "Other biomarkers")) %>%
  select(id.exposure, contains('BCAC')) %>% 
  filter(!(`BCAC 2017` == 0 & `BCAC 2020` == 0))

## all lifestyle
data_sub <- merged %>% 
  filter(!exposure_cat %in% c('Proteins', "Metabolites", "Other biomarkers")) %>%
  select(id.exposure, contains('BCAC'), "ER+", contains("Luminal"), "ER-", "HER2-enriched","TNBC" ) %>% 
  filter(!(`BCAC 2017` == 0 & `BCAC 2020` == 0 & `ER-` == 0  & `ER+` == 0 & `Luminal A` == 0 & `Luminal B1` == 0 & `Luminal B2` == 0 &`HER2-enriched` == 0 & `TNBC` == 0 ))

## BCAC proteins
data_sub <- merged %>% 
  filter(exposure_cat %in% c('Proteins')) %>%
  select(id.exposure, contains('BCAC')) %>% 
  filter(!(`BCAC 2017` == 0 & `BCAC 2020` == 0))

## BCAC metabolites
data_sub <- merged %>% 
  filter(exposure_cat %in% c('Metabolites')) %>%
  filter(!grepl("LDL|HDL|cholest|trigl", exposure, ignore.case = T)) %>% 
  select(id.exposure, contains('BCAC')) %>% 
  filter(!(`BCAC 2017` == 0 & `BCAC 2020` == 0))

## BCAC lipids
data_sub <- merged %>% 
  filter(exposure_cat %in% c( "Metabolites", "Other biomarkers")) %>%
  filter(grepl("LDL|HDL|cholest|trigl", exposure, ignore.case = T)) %>% 
  select(id.exposure, contains('BCAC')) %>% 
  filter(!(`BCAC 2017` == 0 & `BCAC 2020` == 0))

## BCAC + ER+/ER- lipids
data_sub <- merged %>% 
  filter(exposure_cat %in% c( "Metabolites", "Other biomarkers")) %>%
  filter(grepl("LDL|HDL|cholest|trigl", exposure, ignore.case = T)) %>% 
  select(id.exposure, contains('BCAC'), starts_with("ER")) %>% 
  filter(!(`BCAC 2017` == 0 & `BCAC 2020` == 0 & `ER+` == 0  & `ER-` == 0 ))







data_sub2<- reshape2::melt(data_sub, id.var = 'id.exposure') %>% 
  mutate(value=as.character(value)) %>% 
  arrange(value) %>% 
  left_join(merged %>% select(id.exposure, exposure, exposure_cat)) %>% 
  rename('exposure.trait' = 'exposure','exposure.id'= 'id.exposure') %>% 
  create_exposure_categories()

xx<-c("BCAC 2017" , "BCAC 2020" ,"ER+" ,  "Luminal A", "Luminal B1" ,   "Luminal B2"  ,  "ER-", "HER2-enriched", "TNBC"  )

data_sub3 <- data_sub2 %>%
  mutate(exposure = factor(exposure, levels = unique(data_sub2$exposure))) %>% 
  rename(outcome = variable) %>%
  #mutate(outcome = factor(outcome, levels = xx)) %>% 
  left_join(or_ci_data)
  


p<-ggplot(data_sub3, aes( y=exposure, x=outcome, fill = value, id = exposure.id , or_ci= OR_CI, pval =  pval, cat = exposure_cat)) + 
  geom_tile(colour = "grey") + 
  scale_fill_manual(values=c("#4D9221","white", "#C51B7D"))+
  theme(axis.text.x = element_text(angle=45))+ #, colour = a
  coord_flip()
ggplotly(p)

