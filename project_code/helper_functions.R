# query wrapper function
query_epigraphdb_as_table <- function(query){
  results_subset <- query_epigraphdb(
    route = "/cypher",
    params = list(query = query),
    method = "POST",
    mode = "table")
}


kable_it<-function(df){
  library(kableExtra)
  df %>% 
    kable(.) %>%
    kable_styling()
}

makeforest_plot <- function(res_sub, my_title){
  
  pal<-c(wes_palette("Zissou1"))[c(1,3,5)]
  cols <- c("negative" = pal[1], "overlaps null" = pal[2], "positive" = pal[3])
  
  p<-ggplot(res_sub, 
            aes(y=reorder(outcome.details, -or), x=or, label=outcome.details, colour=effect_direction)) +
    geom_errorbarh(aes(xmin=or_loci, xmax=or_upci), height=.3) +
    geom_point(size=2)+
    geom_text(aes(label=OR_CI),hjust=-0.3, vjust=-0.1, size =3, color = 'darkgrey')+
    theme_light()+
    scale_color_manual(values=cols)+
    scale_y_discrete(position = "left")+
    geom_vline(xintercept=1, linetype='longdash') +
    theme(strip.text = element_text(face = 'bold'))+
    facet_wrap(~outcome, scales = 'free_y', ncol = 1) +
    labs(color = "",y = "", x = "Odds ratio", subtitle="",
         title=paste0("         ", my_title))+
    theme(legend.position = "none")
  return(p)
}



extract_outcome_data_custom <- function(exposure_dat, breast_cancer_id){
  out <- extract_outcome_data(
    snps = exposure_dat$SNP,
    outcome = breast_cancer_id,
    proxies = TRUE,
    rsq = 0.8, maf_threshold = 0.3) 
  return(out)
}

run_local_mr_all_outcomes <- function(trait_investigated, outcomes_to_try){
  
  instruments <- extract_instruments(trait_investigated) 
  
  mr_res_all <-data.frame()
  for (outcome_id in outcomes_to_try){
    outcome_dat <- extract_outcome_data_custom(instruments, outcome_id)
    harmonised_dat <-harmonise_data(exposure_dat = instruments, 
                                    outcome_dat = outcome_dat)
    mr_results <- TwoSampleMR::mr(harmonised_dat, 
                     method_list=c('mr_ivw','mr_egger_regression','mr_weighted_median', 'mr_wald_ratio')) %>% 
      split_outcome() %>% 
      split_exposure() %>% 
      separate(outcome, "outcome", sep="[(]") %>% 
      generate_odds_ratios()
    
    mr_res_all<-bind_rows(mr_res_all, mr_results)
  }
  
  mr_res_all_ivw <- mr_res_all %>%
    filter(method %in% c("Inverse variance weighted", "Wald ratio")) %>% 
    rename(outcome.details=id.outcome,
           or_loci=or_lci95,
           or_upci=or_uci95) %>% 
    mutate(OR_CI = paste0(round(or,2), " [",round(or_loci,2) ,":",round(or_upci,2), "]")) %>% 
    mutate(effect_direction = ifelse(or_loci > 1 & or_upci >= 1, 'positive',
                                     ifelse(or_loci < 1 & or_upci <= 1, 'negative', 'overlaps null'))) 
  return(mr_res_all_ivw)
}


## function to call enrich fro pathways
enrich_dbs<-function(gene_list, dbs, adjpval_filter = 0.05){
  
  enriched <- enrichr(gene_list, dbs)
  # flatten list into a table; handle empty tables
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
      filter(Adjusted.P.value < adjpval_filter) %>% 
      separate(Overlap, into = c("found_genes", "total_genes"), sep="/", remove = F)%>% 
      arrange(Odds.Ratio)
  } else{
    enriched_df<-data.frame()
  }
}



### eQTL

get_eQTL_for_snp <- function(variant){
  
  request = httr::GET(url = "http://www.ebi.ac.uk/eqtl/api/associations", 
                      query = list(
                        variant_id = variant,
                        size = 1000,
                        p_upper = 5e-8)
  )
  stopifnot(request$status_code==200)
  
  response = httr::content(request, as = "text", encoding = "UTF-8")
  variant_assoc = jsonlite::fromJSON(response, flatten = TRUE)$`_embedded`$associations
  
  if (length(variant_assoc)!=0){
    for (i in 1:length(variant_assoc)) {
      variant_assoc[[i]]<-variant_assoc[[i]] %>% purrr::discard(is.null) # drop any null items
    }
    variant_assoc_df <- bind_rows(variant_assoc) %>% dplyr::select(rsid, gene_id, qtl_group, pvalue, everything()) %>% arrange(pvalue)
  }else{
    return(data.frame())
  }
  
  nextq = jsonlite::fromJSON(response, flatten = TRUE)$`_links`
  if( any(names(nextq)=="next")){
    stop("API limits 1000 results, but there might be more -- need to investigate")
  }
  return(variant_assoc_df)
}


map_ensembl_to_genename <- function(values){
  
  #retrieve all genes with their GRCh37 coordinates from biomart
  mart_grch37 = biomaRt::useEnsembl(biomart="ensembl",GRCh=37)
  mart_grch37 = biomaRt::useDataset("hsapiens_gene_ensembl", mart_grch37)
  
  # retrieve gene symbols using biomart (eQTL Catalog returns ensembl)
  mart_query = biomaRt::getBM(mart=mart_grch37,
                              attributes=c("ensembl_gene_id","hgnc_symbol"),
                              filters= c("ensembl_gene_id"), 
                              values=unique(values))
  return(mart_query)
}

