
library(gprofiler2)
gostres <- gost(query = eqtl_df$rsid,
                organism = "hsapiens")

# The result is a named list where “result” is a data.frame with the enrichment analysis results
# and “meta” containing a named list with all the metadata for the query.
head(gostres$result)
View(gostres$result)

p <- gostplot(gostres, capped = FALSE, interactive = T)
p

terms<-gostres$result %>% filter(p_value < 0.02) %>% pull(term_id)
pp <- publish_gostplot(p, highlight_terms = terms, 
                       width = NA, height = NA, filename = NULL )
pp



x <- read_tsv("../../../../OneDrive - University of Bristol/Documents - OneDrive/Mini-project2/01_Data/GWAS_tophits/testosterone_bioavailable_tophits.tsv") %>% pull(SNP)
gostres <- gost(query = x,
                #sources = c( "REAC", 'KEGG',"MIRNA", "CORUM", "HP", "HPA", "WP"),
                organism = "hsapiens")

# The result is a named list where “result” is a data.frame with the enrichment analysis results
# and “meta” containing a named list with all the metadata for the query.
head(gostres$result)
View(gostres$result) 


p <- gostplot(gostres, capped = FALSE, interactive = F)
p

terms<-gostres$result %>% filter(recall > 0.20) %>% pull(term_id)
pp <- publish_gostplot(p, highlight_terms = terms, 
                       width = NA, height = NA, filename = NULL )
pp
