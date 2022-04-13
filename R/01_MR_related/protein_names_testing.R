library(tidyverse)


proteins_used <- read_tsv("01_MR_related/mr_evidence_outputs/protein_in_final_set.tsv") %>% select(id.exposure, name=exposure)# %>% filter(!id.exposure %in% c("prot-a-2396", 'prot-a-1540'))

protein_names <- read_csv("01_MR_related/mr_evidence_outputs/protein_names.csv", col_names = F) %>% 
                  rename(name = X1, gene = X2) %>% 
                  right_join(proteins_used) %>% drop_na() # 91: w/o 3 drugs
                  
clipr::write_clip(protein_names$gene) ### to reactome

####

#### Enrichment GWAS-eQTL genes in pathways, using g:Profiler

#library(gprofiler2)
#unique(length(protein_names$gene)) # 90 genes
#gostres <- gost(query = protein_names$gene,
#                sources= c('GO:BP', 'GO:MF',  'KEGG', 'REAC'),
#                organism = "hsapiens")
#p <- gostplot(gostres, capped = FALSE, interactive = T)
#p
#
##g:Profiler results as a table:
#
## recall = intersection_size/term_size
#df<- gostres$result %>% arrange(-recall) %>% 
#  select("source", "term_name",p_value,  "term_size", "intersection_size", recall) %>% 
#  mutate(recall=round(recall,2),
#         p_value= scales::scientific(p_value, digits = 2)) 
#kable_it()
#
#gene_to_id <- gostres$meta$genes_metadata$query$query_1$mapping %>% 
#  as_tibble %>% t() %>% as.data.frame() %>% 
#  rownames_to_column("gene")
#
#failed<- gostres$meta$genes_metadata$failed
#
#protein_names %>% filter(gene %in% failed) %>% View()
#
#clipr::write_clip(gene_to_id$V1)
#####
#
#genes_not_in_rectome <- c("CCL14", "CCDC134", "RARRES1", "TXNDC12", "CST8", "IL6RB", "CUZD1", "SULF2", "NELL1", "EVA1C", "MMRN2", "CPXM1", "KIR2DL5A")
#
#ids_not_in_reactome <- c("ENSG00000130656", "ENSG00000196562", "ENSG00000173269", "ENSG00000189171", "ENSG00000223953", "ENSG00000276409", "ENSG00000100147", "ENSG00000059915", "ENSG00000166979", "ENSG00000128536", "ENSG00000057019", "ENSG00000176444", "ENSG00000066056", "ENSG00000173950", "ENSG00000196576", "ENSG00000134247", "ENSG00000117862", "ENSG00000116199", "ENSG00000115661", "ENSG00000138161", "ENSG00000167178", "ENSG00000125815", "ENSG00000118849", "ENSG00000141504", "ENSG00000088882", "ENSG00000165973")
#
#gene_to_id %>% mutate(missing_gene = ifelse(gene %in% genes_not_in_rectome, "missing", "ok"))%>% 
#              mutate(missing_id = ifelse(V1 %in% ids_not_in_reactome, "missing", "ok")) %>% View()



# review reactome res

res <- read_csv("01_MR_related/pathways/result.csv")
res<- res %>% select(`Pathway name`,  `#Entities found`, `Submitted entities found`) %>% arrange(-`#Entities found`)


not_found <- read_csv("01_MR_related/pathways/not_found.csv") %>% rename(gene =1)
dim(not_found) #12

top3_path_df<-tibble()
paths_by_gene<- tibble()
no_data <- c()

for (i in unique(protein_names$gene)){
  print(i)
  sub<- res %>% filter(grepl(i, `Submitted entities found`)) 
  paths_by_gene <- bind_rows(paths_by_gene, sub %>% mutate(gene = i))
  #sub %>% print()
  cat('\n')
  
  if (dim(sub)[1] != 0){
  top3_path <- res %>% 
    filter(grepl(i, `Submitted entities found`)) %>% 
    #filter(!`Pathway name` %in% c("Immune System", "Metabolism", "Signal Transduction")) %>% 
    slice(1:3) %>% 
    #mutate(path_n = paste0(`Pathway name`, " (", `#Entities found` ,")")) %>% 
    pull(`Pathway name`)
  i_paths <- c(i, top3_path)
  top3_path_df <- rbind(top3_path_df, i_paths)
  } else{
    no_data <- c(no_data, i)
  }
}
colnames(top3_path_df) <- c('gene', '1', '2', '3')
length(unique(top3_path_df$gene)) # 64

top3_path_df %>% count(`1`, sort= T)
top3_path_df %>% count(`1`, `2`, `3`, sort= T)

top3_path_df %>% filter(`1` == "Immune System") %>%  count(`2`, sort= T)
top3_path_df %>% filter(`1` == "Metabolism") %>%  count(`2`, sort= T)


# add genes with missing pathwys to official 'not found' list
not_found <- no_data %>% as_tibble() %>% rename(gene = value) %>% 
  bind_rows(not_found)  %>% mutate(main_path = "Not mapped") %>% distinct()
dim(not_found) #25


# tidy up now
path_df <- top3_path_df %>% select(gene, `1`) %>% 
  mutate(main_path = ifelse(!`1` %in% c('Immune System', 'Metabolism', 'Signal Transduction',
                                           'Developmental Biology'), "Other", `1`)) %>% select(gene, main_path)
dim(path_df)# 64


path_df <- bind_rows(path_df, not_found) %>% left_join(protein_names) %>% distinct() %>% drop_na()
dim(path_df)# 89

path_df %>% count(main_path, sort=T)


path_df %>% write_tsv("01_MR_related/pathways/proteins_w_pathways.tsv")

##















BiocManager::install("ReactomeContentService4R")
library(ReactomeContentService4R)


tp53.re <- map2RefEntities("TP53")
str(tp53.re)

# Extract PhysicalEntities of "TP53"
tp53.all.info <- query(tp53.re$dbId)
dim(tp53.all.info)
head(tp53.all.info$physicalEntity, 5)



# Get Pathways associated with "TP53"
tp53.pathways <- map2Events("TP53", resource = "HGNC", species = "human", mapTo = "pathways")
head(tp53.pathways, 5) 



#A recap for all slots of TP53 Reactome object:
  
str(tp53.all.info, max.level = 1)
