
library(enrichR)
library(tidyverse)

alldbs<-listEnrichrDbs()
pathway_dbs <- c("Reactome_2016", "KEGG_2019_Human", "KEGG_2019_Mouse")
gwa_dbs <- c("GWAS_Catalog_2019", "UK_Biobank_GWAS_v1")
tissue_db<- c("GTEx_Tissue_Sample_Gene_Expression_Profiles_up", 
          "GTEx_Tissue_Sample_Gene_Expression_Profiles_down")


#Using the nearest gene mapping, what pathways are enriched for our gene set?



gwas_genes<- 
  c("CBS", #chr21:43067294 
#"DGKG", #3:186288840  #xqtl
#"ARHGEF11", #1:156934840  #xqtl
#"ZNF687", #151281522.  #xqtl
"CPS1",
"SLC6A12", # both
"DMGDH")

## pathway
enrich_pathways<-function(gene_list, pathway_dbs){

  enriched <- enrichr(gene_list, pathway_dbs)

  for (db_name in names(enriched)){
    if( dim(enriched[[db_name]])[1] > 0){
      enriched[[db_name]]$db <- db_name
    } else  {
      # if it's empty, delete it
      enriched[[db_name]] <- NULL
    }
  }
  enriched_df<-bind_rows(enriched)

  if (dim(enriched_df)[1] > 0){
    enriched_df<- enriched_df %>%
      filter(Adjusted.P.value<0.05) %>% 
      separate(Overlap, into = c("found_genes", "total_genes"), sep="/", remove = F)%>% 
      arrange(Odds.Ratio)
  } else{
    enriched_df<-data.frame()
  }
}

enriched_df<-enrich_pathways(gwas_genes, pathway_dbs)
enriched_df_wo_db<-enriched_df %>% select(-db) %>% distinct()
plotEnrich(enriched_df_wo_db, showTerms = 40, numChar = 50, y = "Count", orderBy = "P.value")

#https://www.molbiotools.com/randomgenesetgenerator.html
random_sets<-read_csv("random_gene_sets_4.csv", col_names =F)
random_enrich <- list()
for (i in 501:1000){
  g_list <- random_sets[i,] %>% as_vector() %>% unname()
  random_enrich[[i]] <-enrich_pathways(g_list, pathway_dbs)
  print(paste0("finished ", i, ";  n_path: ", dim(random_enrich[[i]])[1]))
}
save(random_enrich, file="random_enrich.RData")

plotEnrich(random_enrich[[739]], showTerms = 40, numChar = 50, y = "Count", orderBy = "P.value")

  
case_pathways<- unique(enriched_df$Term)
overlap_count<-data_frame()
for (i in 1:length(random_enrich)) {
  path_list <- random_enrich[[i]]$Term
  shared <- intersect(case_pathways, path_list)
  x<-length(shared)
  print(paste0("overlap with ", i, ": ", x))
  overlap_count<-rbind(overlap_count, c(i,x))
}

overlap_count<-overlap_count %>% distinct()
colnames(overlap_count) <- c("gene_set_no", 'pathways_overlap_w_main_geneset')

overlap_count %>% count(pathways_overlap_w_main_geneset, sort=T)

# find which pathways that came up in my gene set too
bind_rows(random_enrich) %>% count(Term, sort=T) %>% mutate(in_my = ifelse(Term %in%case_pathways, 1, 0)) %>% filter(in_my ==1) %>% View()

# tissue

enriched <- enrichr(gwas_genes, tissue_db)

for (i in 1:length(tissue_db)){
  enriched[[i]]$db <-tissue_db[i]
}

enriched_df<-bind_rows(enriched) %>% 
  #filter(Adjusted.P.value<0.05) %>% 
  separate(Overlap, into = c("found_genes", "total_genes"), sep="/", remove = F)%>% 
  arrange(Odds.Ratio) %>% 
  separate(Term, into = c("id", "tissue", "sex", "age"), sep = "\\s",extra = "merge") %>% 
  select(-id)






known_path<- c("BHMT",
               "SLC44A1",
               "CHDH",
               "ALDH7A1",
               "SARDH",
               "DMGDH")


enriched <- enrichr(known_path, dbs)

for (i in 1:length(dbs)){
  enriched[[i]]$db <-dbs[i]
}

enriched_df<-bind_rows(enriched) %>% 
  filter(Adjusted.P.value<0.05) %>% 
  separate(Overlap, into = c("found_genes", "total_genes"), sep="/") %>% 
  arrange(Odds.Ratio)



plotEnrich(enriched[[4]], showTerms = 40, numChar = 50, y = "Count", orderBy = "P.value")







###]



install.packages("pathfindR")
