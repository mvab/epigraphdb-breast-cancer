library(httr)
library(biomaRt)



variants <- var_to_gene$v.name
instr <- TwoSampleMR::extract_instruments(risk_investigated, p1 = 5e-07, clump=F)
variants <- unique(instr$SNP)

source("helper_functions.R")

eqtl_df <- data.frame()
for (i in 1:length(variants)){
  eqtl_df<-bind_rows(eqtl_df, get_eQTL_for_snp(variants[i]))
  print(paste0("done: ", i))
}
gene_names <- map_ensembl_to_genename(eqtl_df$gene_id)

eqtl_df<-left_join(eqtl_df, gene_names, 
                   by = c('gene_id' = 'ensembl_gene_id')) %>% 
        dplyr::select('hgnc_symbol', everything())





####



extract_mr <- function(outcome_trait, gene_list, qtl_type) {
  endpoint <- "/xqtl/single-snp-mr"
  per_gene <- function(gene_name) {
    params <- list(
      exposure_gene = gene_name,
      outcome_trait = outcome_trait,
      qtl_type = qtl_type,
      pval_threshold = 1e-8
    )
    df <- query_epigraphdb(route = endpoint, params = params, mode = "table")
    df
  }
  res_df <- gene_list %>% map_df(per_gene)
  res_df
}

xqtl_df <- c("pQTL", "eQTL") %>% map_df(function(qtl_type) {
  extract_mr(
    outcome_trait = OUTCOME_TRAIT,
    gene_list = gene_list,
    qtl_type = qtl_type
  ) %>%
    mutate(qtl_type = qtl_type)
})
xqtl_df



extract_instruments('met-a-362')
