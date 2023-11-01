### functions


# extract literature space R function

extract_literature_space <- function(id){
  
  ## this query is too heavy for R when literature space is huge (like for BC)
  
  # query =  paste0("
  #     MATCH (gwas:Gwas)-[gs1:GWAS_TO_LITERATURE_TRIPLE]->(s1:LiteratureTriple) -[:SEMMEDDB_OBJ]->(st:LiteratureTerm)
  #     WHERE gwas.id = '", id, "'
  #     AND gs1.pval < 0.01
  #     MATCH (s1)-[:SEMMEDDB_SUB]-(st1:LiteratureTerm) 
  #     MATCH (gwas)-[:GWAS_TO_LITERATURE]-(lit:Literature)-[]-(s1)
  #     RETURN lit.id, lit.year,  gwas {.id, .trait}, 
  #     gs1 {.pval, .localCount}, st1 {.name, .type}, s1 {.id, .subject_id, .object_id, .predicate}, st {.name, .type}
  #     ")
  # out<-query_epigraphdb_as_table(query)
  # dim(out)
  
  
  # so we're going to split it in two parts (tested that this approach delivers the same result)
  
  # 1: extract GWAS triples and their terms
  query1 =  paste0("
      MATCH (gwas:Gwas)-[gs1:GWAS_TO_LITERATURE_TRIPLE]->(s1:LiteratureTriple) -[:SEMMEDDB_OBJ]->(st:LiteratureTerm)
      WHERE gwas.id = '", id, "'
      AND gs1.pval < 0.01
      MATCH (s1)-[:SEMMEDDB_SUB]-(st1:LiteratureTerm) 
      RETURN  gwas {.id, .trait}, 
      gs1 {.pval, .localCount}, st1 {.name, .type}, s1 {.id, .subject_id, .object_id, .predicate}, st {.name, .type}
      ")
  # gwas and triples  
  out1<-query_epigraphdb_as_table(query1)
  
  if (dim(out1)[1]==0 ){
    print(paste0("======= ", id, " : empty literature space"))
    return(data.frame())
    
  } else {
    
    dim(out1) # 24239
    triple_id_list <- unique(out1$s1.id) #23278
    length(triple_id_list)
    
    
    # 2: for these triples, extract literature nodes
    query2 =  paste0("
        MATCH (s1:LiteratureTriple)-[]-(lit:Literature)
        WHERE s1.id in ['", paste0(triple_id_list, collapse = "', '"),"'] 
        MATCH (lit)-[:GWAS_TO_LITERATURE]-(gwas:Gwas)
        WHERE gwas.id = '", id, "'
        RETURN gwas.id, lit.id, lit.year, s1 {.id, .subject_id, .object_id, .predicate}
        ")
    
    out2<-query_epigraphdb_as_table(query2) 
    dim(out2)
    
    triple_id_list2 <- unique(out2$s1.id) #23142
    length(triple_id_list2)
    
    litspace <- full_join(out1, out2, 
                          by = c("s1.subject_id"="s1.subject_id",
                                 "s1.predicate"="s1.predicate",
                                 "s1.id" = "s1.id" ,
                                 'gwas.id' = 'gwas.id',
                                 "s1.object_id" ="s1.object_id" )) %>% 
      filter(!is.na(lit.id))
    
    trait <- unique(litspace$gwas.trait)
    
    print(paste0("Lit space of ", trait, " - ", id, ": ", dim(litspace)[1] ))  
    dim(litspace)
    
    litspace<- litspace %>% rowwise() %>% 
      mutate(st1.type = paste0(unlist(st1.type), collapse="/")) %>%  
      mutate(st.type = paste0(unlist(st.type), collapse="/"))
    return(litspace)
  }
  
}



tidy_lit_space <- function(dat){
  
  triples_tidy <- dat %>% 
    #filter(lit.year >=1990) %>% 
    tidy_gwas_to_lit_output() %>%  distinct() 
  
  #length(unique(triples_tidy$lit.id))# unique papers
  
  triples_tidy_count<-triples_tidy %>% 
    select(lit.id, term1, predicate, term2)  %>% distinct() %>% 
    group_by(term1, predicate, term2) %>% 
    count() %>% ungroup() %>% 
    distinct() 
  
  # add type to disease and drug nodes (selective)
  node_types<-triples_tidy %>%  make_node_types()
  
  # find most common things
  #node_type_counts<-node_types %>% count(name, name='size', sort=T)
  
  # add term types to counts
  triples_tidy_count<-triples_tidy_count %>%
    left_join(node_types, by=c('term1'='name')) %>% 
    rename('term1.type'='type', 'term1.type_verbose'='type_verbose') %>%  
    left_join(node_types, by=c('term2'='name')) %>% 
    rename('term2.type'='type', 'term2.type_verbose'='type_verbose') 
  
  # calculate counts by pairs
  pair_counts <- triples_tidy_count %>% 
    group_by(term1, term2) %>% 
    summarise(n_pair = sum(n)) %>%
    ungroup() %>% distinct() %>% 
    mutate(pair = paste0(term1," / ",term2)) # counts by pairs
  
  # extract collapsed pubmed id + year
  pubmed_year <- 
    triples_tidy %>%  select(lit.id, lit.year, term1, predicate, term2)  %>% distinct() %>% 
    mutate(pubmed_year = paste0(lit.id, " (", lit.year, ")")) %>% 
    select(-lit.id, -lit.year) %>% 
    group_by(term1, predicate, term2) %>% 
    summarise(pubmed_year = paste0(pubmed_year, collapse=" / ")) 
  
  triples_tidy_count  <- triples_tidy_count %>% 
    left_join(pair_counts) %>% 
    select(term1, predicate, term2, n_triple = n, n_pair, pair, everything()) %>% 
    left_join(pubmed_year)
  
  
  return(triples_tidy_count)
}




tidy_gwas_to_lit_output <- function(dat){  
  dat %>% 
    select(-one_of(
      'gwas.id', 'gwas.trait','gs1.pval' ,"s1.subject_id", "s1.id", "s1.object_id"  ))%>% 
    rename(term1 = st1.name,
           term2  = st.name,
           predicate = s1.predicate) %>% 
    # tidy up names
    mutate(term1 = gsub(" gene", "", term1), term2 = gsub(" gene", "", term2)) %>% 
    mutate(term1 = gsub(" protein, human", "", term1), term2 = gsub(" protein, human", "", term2)) %>% 
    mutate(term1 = gsub(" protein, mammalian", "", term1), term2 = gsub(" protein, mammalian", "", term2)) %>% 
    mutate(term1 = gsub(", human", "", term1, ignore.case = T), term2 = gsub(", human", "", term2, ignore.case = T)) %>% 
    mutate(term1 = gsub("BRCA1 Protein", "BRCA1", term1, ignore.case = T), term2 = gsub("BRCA1 Protein", "BRCA1", term2, ignore.case = T)) %>% 
    mutate(term1 = gsub("BRCA2 Protein", "BRCA2", term1, ignore.case = T), term2 = gsub("BRCA2 Protein", "BRCA2", term2, ignore.case = T)) %>% 
    mutate(term1 = gsub("BCL-2 Protein", "BCL-2", term1, ignore.case = T), term2 = gsub("BCL-2 Protein", "BCL-2", term2, ignore.case = T)) %>% 
    distinct() %>% 
    filter(!(term2==term1)) 
}  




# add type to disease and drug nodes (selective)
make_node_types <- function(dat){
  
  drug_list <- c('fulvestrant', 'vinorelbine')
  generic_list <-c('Membrane Transport Proteins', 'Small Molecule')
  
  node_types_full<-bind_rows(dat %>% select(name = term1, type =st1.type),
                        dat %>% select(name = term2, type = st.type))
  
  node_types <- node_types_full%>% distinct()
  
  node_type_counts <- node_types %>%  count(name)
  
  terms_w_single_type <- node_type_counts %>% filter(n == 1) 
  
  terms_w_multp_types <- node_type_counts %>% filter(n > 1) 
  
  if (dim(terms_w_multp_types)[1] > 0){
    terms_w_multp_types_to_keep <-  keep_one_type(terms_w_multp_types, node_types_full)
    
    node_types <- bind_rows(
      node_types %>% filter(name %in% terms_w_single_type$name),
      terms_w_multp_types_to_keep %>% select(-n)
    )
  }
  
  
  node_types %>% 
  mutate(type_verbose= case_when(grepl('dsyn', type) ~ 'disease',
                                 #type == "['orch']|['phsu']|['orch', 'phsu']|['hops']|['orch', 'phsu', 'hops']|['orch', 'hops']" ~ 'drug_or_compound',
                                 type == "orch" ~ 'drug_or_compound',
                                 type == "phsu" ~ 'drug_or_compound',
                                 type == "hops" ~ 'drug_or_compound',
                                 type == "orch/phsu" ~ 'drug_or_compound',
                                 type == "orch/phsu/hops" ~ 'drug_or_compound',
                                 type == "orch/hops" ~ 'drug_or_compound',
                                 grepl('mab$',name) ~ 'drug_or_compound',
                                 grepl('tinib',name) ~ 'drug_or_compound',
                                 grepl('cisplatin',name) ~ 'drug_or_compound',
                                 grepl('methotrexate polyglutamate',name) ~ 'drug_or_compound',
                                 grepl('herceptin',name, ignore.case = T) ~ 'drug_or_compound',
                                 grepl('capecitabine',name) ~ 'drug_or_compound',
                                 grepl('fluorouracil',name) ~ 'drug_or_compound',
                                 grepl('cyclophosphamide',name) ~ 'drug_or_compound',
                                 grepl('exemestane',name) ~ 'drug_or_compound',
                                 grepl('ado-trastuzumab emtansine',name) ~ 'drug_or_compound',
                                 grepl('megestrol acetate',name) ~ 'drug_or_compound',
                                 grepl('thiotepa',name) ~ 'drug_or_compound',
                                 grepl('tamoxifen',name) ~ 'drug_or_compound',
                                 grepl('zole$',name) ~ 'drug_or_compound',
                                 grepl('taxel$',name) ~ 'drug_or_compound',
                                 name %in% drug_list ~ 'drug_or_compound',
                                 name %in% generic_list ~ 'generic',
                                 TRUE ~ 'any')) %>% 
  distinct()
}


keep_one_type <- function(terms_w_multp_types, node_types_full){
  
  node_types_count_w_mulp<-node_types_full %>% 
    filter(name %in% terms_w_multp_types$name) %>% 
    group_by(name, type) %>% 
    count() %>% ungroup()
  
  selected_types <- tibble()
  
  for (i in terms_w_multp_types$name){
    
    tmp1 <- node_types_count_w_mulp %>% filter(name == i)
    
    if (length(unique(tmp1$n))!=1){ 
      # types appear with different frequencies:
      
      # sort by freq of type to get the top one
      tmp2<- tmp1 %>% 
        arrange(-n) %>% 
        filter(row_number()==1)
      
    } else {
      # all types appear with the same freq:
      
      # order by type, to keep the longer name
      tmp2<-tmp1 %>%
        mutate(len_type = str_length(type)) %>%
        arrange(-len_type) %>% 
        filter(row_number()==1) %>% 
        select(-len_type)
    }
    selected_types<-bind_rows(selected_types, tmp2)
  }
  
  return(selected_types)
}



tidy_terms_for_viz <- function(df){
  
  df %>% mutate_at(
    vars(one_of('term1', 'term2')),
    funs(case_when(
      . == 'Insulin-Like Growth Factor I' ~ 'IGF1',
      . == 'Insulin-Like Growth Factor II' ~ 'IGF2',
      . == 'Insulin-Like-Growth Factor I Receptor'~ 'IGF1R',
      . == 'Insulin-Like Growth Factor Receptor'~ 'IGFR',
      . == 'Insulin-Like-Growth Factor II Receptor'~ 'IGF2R',
      . == 'Receptor, IGF Type 2'~ 'IGF2R',
      . == 'INS' ~ 'Insulin',
      . == 'INSR' ~ 'Insulin Receptor',
      . == 'Insulin-Like Growth-Factor-Binding Proteins' ~ 'IGFBP',
      . == 'Insulin-Like Growth Factor Binding Protein 3' ~ 'IGFBP3',
      . == 'Insulin-Like Growth Factor Binding Protein 4' ~ 'IGFBP4',
      . == "Insulin-Like Growth-Factor Binding Protein 1" ~ 'IGFBP1',
      . == 'Insulin-Like Growth Factor Binding Protein 5' ~ 'IGFBP5',
      . == "Insulin-Like Growth Factor Binding Protein 6" ~ 'IGFBP6',
      . == 'EGF' ~ 'epidermal growth factor',
      . == 'GHR' ~ 'Growth Hormone Receptor',
      . == 'LEP' ~ 'Leptin',
      . == 'LEPR' ~ 'leptin receptor',
      . == 'leptin' ~ 'Leptin',
      . == 'IRS1' ~  'insulin receptor substrate 1 protein',
      . == 'PRL' ~ 'Prolactin',
      . == 'AR'~ 'Androgen Receptor',
      . == 'PGR'~ 'Progesterone receptor',
      . == 'CRP' ~  'C-reactive protein',
      . == 'Receptors, Steroid' ~ 'Steroid receptor',
      . == 'Receptors, LH' ~ 'luteinizing hormone receptor',
      . == 'Receptors, Progesterone'~ 'Progesterone receptor',
      . == 'MC4R' ~ 'Melanocortin 4 Receptor',
      . == 'MC3R' ~ 'Melanocortin 3 Receptor',
      . == 'RETN' ~ 'resistin',
      
      . == 'Interleukin-1' ~ "IL1", 
      . == 'Interleukin 2 Receptor' ~ "IL2R", 
      . == 'interleukin-1 receptor accessory protein' ~ "IL1RAP", 
      . == 'Interleukin-1 Receptor-Associated Kinase 2' ~ "IL1RAK2", 
      . == 'Interleukin-1 Receptor-Associated Kinases' ~ "IL1RAK", 
      . == 'interleukin-1, beta' ~ "IL1b", 
      . == 'IL1B' ~ "IL1b", 
      . == 'Receptors, Interleukin-1' ~ "IL1", 
      . == 'Interleukin 2 Receptor' ~ "IL2R", 
      . == 'Interleukin-3 Receptor' ~ "IL3R", 
      . == 'Interleukin 2 Receptor, Alpha' ~ "IL2Ra",
      . == 'interleukin-7 receptor, alpha chain' ~ "IL7Ra",
      . == 'Interleukin-1 alpha' ~ "IL1a", 
      . == 'interleukin-10' ~ "IL10", 
      . == 'Interleukin-2' ~ "IL2", 
      . == 'interleukin-4' ~ "IL4", 
      . == 'Interleukin-5' ~ "IL5", 
      . == 'interleukin-6' ~ "IL6", 
      . == 'interleukin-8' ~ "IL8", 
      . == 'Interleukin-12' ~ "IL12", 
      . == 'Interleukin-13' ~ "IL13", 
      . == 'Interleukin-17' ~ "IL17", 
      . == 'interleukin-18' ~ "IL18", 
      . == 'interleukin-22' ~ "IL22", 
      . == 'interleukin-23' ~ "IL23", 
      . == 'Interleukin-18' ~ "IL18", 
      
      . == 'Tumor Necrosis Factors' ~ 'TNF',
      
      
      . == 'glucocorticoid receptor alpha' ~ 'NR3C1',
      . == '1-Phosphatidylinositol 3-Kinase' ~ 'PIK3',
      
      . == 'LIF' ~ 'Leukemia inhibitory factor',
      
      . == 'Tumor Necrosis Factor Receptor Superfamily, Member 10B' ~ 'TNFRSF10B',
      . == 'tumor necrosis factor receptor 1A' ~ 'TNFRSF1A',
      . == 'Tumor necrosis factor receptor 11b' ~ 'TNFRSF11B',

      . == 'Intercellular adhesion molecule 1' ~ "ICAM1", 
      . == 'Intercellular cell adhesion molecule' ~ "ICAM1", 
      . == 'Intercellular Adhesion Molecules' ~ "ICAM1", 
      
      . == 'Toll-Like Receptor 5' ~ 'TLR5',
      . == 'Toll-Like Receptor 2' ~ 'TLR2',
      . == 'toll-like receptor 4' ~ 'TLR4',
      . == 'Toll-Like Receptor 9' ~ 'TLR9',
      . == 'Toll-Like Receptor 6' ~ 'TLR6',
      . == 'Toll-Like Receptor 7' ~ 'TLR7',
      . == 'Toll-like receptors' ~ 'TLR',
      
      
      . == 'Fibroblast Growth Factor 7' ~ 'FGF7',
      . == 'Fibroblast Growth Factor 2' ~ 'FGF2',
      . == 'Fibroblast Growth Factor 1' ~ 'FGF1',
      . == 'Fibroblast Growth Factor Receptor 1' ~ 'FGFR1',
      . == 'Fibroblast Growth Factor Receptor 2' ~ 'FGFR2',
      . == 'Fibroblast Growth Factor Receptors' ~ 'FGFR',
      . == 'Fibroblast Growth Factor' ~ 'FGF',
      
      . == 'Sex Hormone-Binding Globulin' ~ 'SHBG',
      
      . == "Janus kinase 2" ~ "JAK2",
      . == 'ALB' ~ 'Albumin',
      . == 'Albumins' ~ 'Albumin',
      
      . == 'VEGF' ~ 'VEGFA',
      . == 'Vascular Endothelial Growth Factor Receptor-2'~ 'VEGFR2',

      
      . ==  'TFRC' ~ "Transferrin Receptor",
      . ==  'TF' ~ "Transferrin",
      
      . == 'OSM' ~ 'oncostatin M',
      
      . == 'HGF' ~ 'Hepatocyte Growth Factor',
      
      . == 'Proto-Oncogene Proteins c-akt' ~ "AKT1",
      . == 'Epidermal Growth Factor Receptor' ~ "EGFR",
      
      . == 'Urokinase Plasminogen Activator Receptor' ~ 'PLAUR',
      
      . == 'CTF1' ~ "cardiotrophin 1",
      . == 'CD40LG' ~ 'CD40 Ligand',
      
      . == "CCND1" ~ 'Cyclin D1',
      . == "BCL2L11" ~ 'BCL2-Related Protein 11',
      . == 'CLU' ~ 'Clusterin',
      
      . == 'CTSD' ~ "CATHEPSIN D",
      . == 'MSTN' ~ "myostatin",
      . == 'SST' ~ "Somatostatin", 
      . == 'SRC' ~ "src-Family Kinases",
      . == 'TXN' ~ "Thioredoxin",
      . == 'TYR' ~ "tyrosine",
      . == 'VIM' ~ "Vimentin",
      
      . == 'ADIPOQ' ~ "Adiponectin",
      . == 'APOB' ~ "Apolipoproteins B",
      . == 'C3' ~ "Complement 3",

      

      
      TRUE ~ .)))
  
}




get_breast_cancer_triples <- function(bc_triples_tidy_count){
  
  bc_triple1 <- bc_triples_tidy_count %>%  
    filter(term2 %in% c("Breast Diseases","malignant disease")) %>% 
    filter(!predicate %in% c("COEXISTS_WITH", 'NEG_ISA')) %>%
    filter(term1.type_verbose != "drug_or_compound") %>% 
    filter(term1.type != 'inch/phsu') %>% 
    filter(term1 != 'Prostate-Specific Antigen') %>% 
    mutate(term2 = ifelse(term2 == "malignant disease","Breast Diseases", term2 )) %>% 
    select(1:5) %>% 
    group_by(term1,predicate, term2) %>% 
    summarise_all(sum) %>% ungroup()
  
  bc_triple1 %>% select(term1, term2) %>% distinct() %>% dim() # total 46
  
  bc_triple2 <- bc_triples_tidy_count %>%  
    filter(term2 %in% bc_triple1$term1) %>% #& !term1 %in% bc_triple1$term1) %>%
    filter(term1.type_verbose != "drug_or_compound") %>%  
    filter(term1.type != 'inch/phsu') %>%   
    filter(!predicate %in% c("COEXISTS_WITH", 'NEG_COEXISTS_WITH'))
  
  bc_triple2 %>% select(term1, term2) %>% distinct() %>% dim() # total 1291
  
  
  ### TIDYING
  
  bc_triple1_tidy <- bc_triple1 %>% 
    tidy_terms_for_viz() %>% 
    select(term1,term2,n_pair) %>% 
    filter(term1 !=term2) %>% 
    group_by(term1, term2) %>% 
    slice(which.max(n_pair)) %>%
    ungroup()  %>% 
    rename(n=n_pair) 
  
  
  # need to drop reverse connections: keep A-B or B-A depending which one is more common
  bc_triple2_onedir<- bc_triple2 %>% 
    # forward count
    rename(n_pair_f = n_pair) %>% 
    # join col that will show reverse pair count
    left_join(bc_triple2 %>% select(term1,term2,n_pair_b = n_pair), by = c('term1' = 'term2', 'term2' = 'term1')) %>% 
    distinct() %>% 
    # if rel does not exist uin reverse, set it to 0
    mutate(across(n_pair_b, ~replace_na(.x, 0))) %>% 
    # if f more common, keep it, else, keep reverse
    mutate(keep = ifelse(n_pair_f >= n_pair_b ,T,F)) %>% 
    filter(keep==T) %>%  rename(n_pair=n_pair_f)
  
  bc_triple2_tidy <- bc_triple2_onedir %>% 
    tidy_terms_for_viz() %>% 
    # after tidying names you get almost duplicates: keep the one with highest n
    select(term1,term2,n_pair) %>% 
    filter(term1 !=term2) %>% 
    group_by(term1, term2) %>% 
    slice(which.max(n_pair)) %>%
    ungroup()  %>% 
    rename(n=n_pair) 
  
  
  ## single path 
  # pick terms linked to anchor
  bc_x<- bc_triple1_tidy %>% filter(n>1)
  # make sure they are term1 in triple2
  bc_y<- bc_triple2_tidy %>% filter(term2 %in% bc_x$term1 )  %>% filter(n>13) # for viz use 13
  # drop triple 1 term1 that ends up not linked to anything because of previous filteting step
  bc_x2 <- bc_x %>%  filter(term1 %in% bc_y$term2)
  
  bc_twostep_triples<- bind_rows(bc_x2, bc_y)
  bc_s_n<- make_sankey(bc_twostep_triples, fontSize=13)
  
  
  
  bc_triple3 <- bc_triples_tidy_count %>%  
    filter(term2 %in% bc_triple2$term1 & !term1 %in% bc_triple2$term1) %>%
    filter(term1.type_verbose != "drug_or_compound") %>%  
    filter(term1.type != 'inch/phsu') %>%   
    filter(!predicate %in% c("COEXISTS_WITH", 'NEG_COEXISTS_WITH'))
  
  bc_triple3 %>% select(term1, term2) %>% distinct() %>% dim() # total 1767
  
  # need to drop reverse connections: keep A-B or B-A depending which one is more common
  bc_triple3_onedir<- bc_triple3 %>% 
    # forward count
    rename(n_pair_f = n_pair) %>% 
    # join col that will show reverse pair count
    left_join(bc_triple2 %>% select(term1,term2,n_pair_b = n_pair), by = c('term1' = 'term2', 'term2' = 'term1')) %>% 
    distinct() %>% 
    # if rel does not exist uin reverse, set it to 0
    mutate(across(n_pair_b, ~replace_na(.x, 0))) %>% 
    # if f more common, keep it, else, keep reverse
    mutate(keep = ifelse(n_pair_f >= n_pair_b ,T,F)) %>% 
    filter(keep==T) %>%  rename(n_pair=n_pair_f)
  
  bc_triple3_tidy <- bc_triple3_onedir %>% 
    tidy_terms_for_viz() %>% 
    # after tidying names you get almost duplicates: keep the one with highest n
    select(term1,term2,n_pair) %>% 
    filter(term1 !=term2) %>% 
    group_by(term1, term2) %>% 
    slice(which.max(n_pair)) %>%
    ungroup()  %>% 
    rename(n=n_pair) 
  
  bc_triple3_tidy %>% select(term1, term2) %>% distinct() %>% dim() # 1749
  
  
  
  
  bc_triple4 <- bc_triples_tidy_count %>%  
    filter(term2 %in% bc_triple3$term1 & !term1 %in% bc_triple3$term1 & !term1 %in% bc_triple2$term1) %>%
    filter(term1.type_verbose != "drug_or_compound") %>%  
    filter(term1.type != 'inch/phsu') %>%   
    filter(!predicate %in% c("COEXISTS_WITH", 'NEG_COEXISTS_WITH'))
  
  bc_triple4 %>% select(term1, term2) %>% distinct() %>% dim() # total 375
  
  # need to drop reverse connections: keep A-B or B-A depending which one is more common
  bc_triple4_onedir<- bc_triple4 %>% 
    # forward count
    rename(n_pair_f = n_pair) %>% 
    # join col that will show reverse pair count
    left_join(bc_triple2 %>% select(term1,term2,n_pair_b = n_pair), by = c('term1' = 'term2', 'term2' = 'term1')) %>% 
    distinct() %>% 
    # if rel does not exist uin reverse, set it to 0
    mutate(across(n_pair_b, ~replace_na(.x, 0))) %>% 
    # if f more common, keep it, else, keep reverse
    mutate(keep = ifelse(n_pair_f >= n_pair_b ,T,F)) %>% 
    filter(keep==T) %>%  rename(n_pair=n_pair_f)
  
  bc_triple4_tidy <- bc_triple4_onedir %>% 
    tidy_terms_for_viz() %>% 
    # after tidying names you get almost duplicates: keep the one with highest n
    select(term1,term2,n_pair) %>% 
    filter(term1 !=term2) %>% 
    group_by(term1, term2) %>% 
    slice(which.max(n_pair)) %>%
    ungroup()  %>% 
    rename(n=n_pair) 
  
  bc_triple4_tidy %>% select(term1, term2) %>% distinct() %>% dim() # 375
  
  
  bc_triple1_tidy <- bc_triple1_tidy %>% mutate(term2 = ifelse(term2 == 'Breast Diseases', 'Breast cancer', term2))
  
  
  return(list(triple1 = bc_triple1_tidy,
              triple2 = bc_triple2_tidy,
              triple3 = bc_triple3_tidy,
              triple4 = bc_triple4_tidy))
  
}



get_outcome_diesease_triples <- function(bc_triples_tidy_count, outcomes, outcome_name){
  
  bc_triple1 <- bc_triples_tidy_count %>%  
    filter(term2 %in% outcomes) %>% 
    filter(!predicate %in% c("COEXISTS_WITH", 'NEG_ISA')) %>%
    #filter(term1.type_verbose != "drug_or_compound") %>% 
    #filter(term1.type != 'inch/phsu') %>% 
    filter(term1 != 'Prostate-Specific Antigen') %>% 
    mutate(term2 = ifelse(term2 %in% outcomes,outcome_name, term2 )) %>% 
    select(1:5) %>% 
    group_by(term1,predicate, term2) %>% 
    summarise_all(sum) %>% ungroup()
  
  bc_triple1 %>% select(term1, term2) %>% distinct() %>% dim() # total 46
  
  bc_triple2 <- bc_triples_tidy_count %>%  
    filter(term2 %in% bc_triple1$term1) %>% #& !term1 %in% bc_triple1$term1) %>%
    #filter(term1.type_verbose != "drug_or_compound") %>%  
    #filter(term1.type != 'inch/phsu') %>%   
    filter(!predicate %in% c("COEXISTS_WITH", 'NEG_COEXISTS_WITH'))
  
  bc_triple2 %>% select(term1, term2) %>% distinct() %>% dim() # total 1291
  
  
  ### TIDYING
  
  bc_triple1_tidy <- bc_triple1 %>% 
    tidy_terms_for_viz() %>% 
    select(term1,term2,n_pair) %>% 
    filter(term1 !=term2) %>% 
    group_by(term1, term2) %>% 
    slice(which.max(n_pair)) %>%
    ungroup()  %>% 
    rename(n=n_pair) 
  
  
  # need to drop reverse connections: keep A-B or B-A depending which one is more common
  bc_triple2_onedir<- bc_triple2 %>% 
    # forward count
    rename(n_pair_f = n_pair) %>% 
    # join col that will show reverse pair count
    left_join(bc_triple2 %>% select(term1,term2,n_pair_b = n_pair), by = c('term1' = 'term2', 'term2' = 'term1')) %>% 
    distinct() %>% 
    # if rel does not exist uin reverse, set it to 0
    mutate(across(n_pair_b, ~replace_na(.x, 0))) %>% 
    # if f more common, keep it, else, keep reverse
    mutate(keep = ifelse(n_pair_f >= n_pair_b ,T,F)) %>% 
    filter(keep==T) %>%  rename(n_pair=n_pair_f)
  
  bc_triple2_tidy <- bc_triple2_onedir %>% 
    tidy_terms_for_viz() %>% 
    # after tidying names you get almost duplicates: keep the one with highest n
    select(term1,term2,n_pair) %>% 
    filter(term1 !=term2) %>% 
    group_by(term1, term2) %>% 
    slice(which.max(n_pair)) %>%
    ungroup()  %>% 
    rename(n=n_pair) 
  
  
  ## single path 
  # pick terms linked to anchor
  bc_x<- bc_triple1_tidy %>% filter(n>1)
  # make sure they are term1 in triple2
  bc_y<- bc_triple2_tidy %>% filter(term2 %in% bc_x$term1 )  %>% filter(n>13) # for viz use 13
  # drop triple 1 term1 that ends up not linked to anything because of previous filteting step
  bc_x2 <- bc_x %>%  filter(term1 %in% bc_y$term2)
  
  bc_twostep_triples<- bind_rows(bc_x2, bc_y)
  bc_s_n<- make_sankey(bc_twostep_triples, fontSize=13)
  
  
  
  bc_triple3 <- bc_triples_tidy_count %>%  
    filter(term2 %in% bc_triple2$term1 & !term1 %in% bc_triple2$term1) %>%
    #filter(term1.type_verbose != "drug_or_compound") %>%  
    #filter(term1.type != 'inch/phsu') %>%   
    filter(!predicate %in% c("COEXISTS_WITH", 'NEG_COEXISTS_WITH'))
  
  bc_triple3 %>% select(term1, term2) %>% distinct() %>% dim() # total 1767
  
  # need to drop reverse connections: keep A-B or B-A depending which one is more common
  bc_triple3_onedir<- bc_triple3 %>% 
    # forward count
    rename(n_pair_f = n_pair) %>% 
    # join col that will show reverse pair count
    left_join(bc_triple2 %>% select(term1,term2,n_pair_b = n_pair), by = c('term1' = 'term2', 'term2' = 'term1')) %>% 
    distinct() %>% 
    # if rel does not exist uin reverse, set it to 0
    mutate(across(n_pair_b, ~replace_na(.x, 0))) %>% 
    # if f more common, keep it, else, keep reverse
    mutate(keep = ifelse(n_pair_f >= n_pair_b ,T,F)) %>% 
    filter(keep==T) %>%  rename(n_pair=n_pair_f)
  
  bc_triple3_tidy <- bc_triple3_onedir %>% 
    tidy_terms_for_viz() %>% 
    # after tidying names you get almost duplicates: keep the one with highest n
    select(term1,term2,n_pair) %>% 
    filter(term1 !=term2) %>% 
    group_by(term1, term2) %>% 
    slice(which.max(n_pair)) %>%
    ungroup()  %>% 
    rename(n=n_pair) 
  
  bc_triple3_tidy %>% select(term1, term2) %>% distinct() %>% dim() # 1749
  
  
  
  
  bc_triple4 <- bc_triples_tidy_count %>%  
    filter(term2 %in% bc_triple3$term1 & !term1 %in% bc_triple3$term1 & !term1 %in% bc_triple2$term1) %>%
    #filter(term1.type_verbose != "drug_or_compound") %>%  
    #filter(term1.type != 'inch/phsu') %>%   
    filter(!predicate %in% c("COEXISTS_WITH", 'NEG_COEXISTS_WITH'))
  
  bc_triple4 %>% select(term1, term2) %>% distinct() %>% dim() # total 375
  
  # need to drop reverse connections: keep A-B or B-A depending which one is more common
  bc_triple4_onedir<- bc_triple4 %>% 
    # forward count
    rename(n_pair_f = n_pair) %>% 
    # join col that will show reverse pair count
    left_join(bc_triple2 %>% select(term1,term2,n_pair_b = n_pair), by = c('term1' = 'term2', 'term2' = 'term1')) %>% 
    distinct() %>% 
    # if rel does not exist uin reverse, set it to 0
    mutate(across(n_pair_b, ~replace_na(.x, 0))) %>% 
    # if f more common, keep it, else, keep reverse
    mutate(keep = ifelse(n_pair_f >= n_pair_b ,T,F)) %>% 
    filter(keep==T) %>%  rename(n_pair=n_pair_f)
  
  bc_triple4_tidy <- bc_triple4_onedir %>% 
    tidy_terms_for_viz() %>% 
    # after tidying names you get almost duplicates: keep the one with highest n
    select(term1,term2,n_pair) %>% 
    filter(term1 !=term2) %>% 
    group_by(term1, term2) %>% 
    slice(which.max(n_pair)) %>%
    ungroup()  %>% 
    rename(n=n_pair) 
  
  bc_triple4_tidy %>% select(term1, term2) %>% distinct() %>% dim() # 375
  
  
  bc_triple1_tidy <- bc_triple1_tidy %>% mutate(term2 = ifelse(term2 %in% outcomes, outcome_name, term2))
  
  
  return(list(triple1 = bc_triple1_tidy,
              triple2 = bc_triple2_tidy,
              triple3 = bc_triple3_tidy,
              triple4 = bc_triple4_tidy))
  
}




## another method that searches for links in 4 BC triples! from 2 trait triples

overlap_trait_and_bc <- function(trait_twostep_triples, KEY_TERM, n_filter = 1, bc_triples,sankey_font = 13 , outcome = 'Breast cancer'){
  
  bc_triple1_tidy = bc_triples$triple1
  bc_triple2_tidy = bc_triples$triple2
  bc_triple3_tidy = bc_triples$triple3
  bc_triple4_tidy = bc_triples$triple4
  
  triple1 = trait_twostep_triples %>% filter(term1 %in% bc_triple1_tidy$term1 | term2 %in% bc_triple1_tidy$term1) %>% pull(term2) %>% unique()
  triple2 = trait_twostep_triples %>% filter(term1 %in% bc_triple2_tidy$term1 | term2 %in% bc_triple2_tidy$term1) %>% pull(term2) %>% unique()
  triple3 = trait_twostep_triples %>% filter(term1 %in% bc_triple3_tidy$term1 | term2 %in% bc_triple3_tidy$term1) %>% pull(term2) %>% unique()
  triple4 = trait_twostep_triples %>% filter(term1 %in% bc_triple4_tidy$term1 | term2 %in% bc_triple4_tidy$term1) %>% pull(term2) %>% unique()
  
  
  q <- bc_triple4_tidy %>% filter(term1 %in% triple4 ) 
  x <- bc_triple3_tidy %>% filter(term1 %in% q$term2 | term1 %in% triple3) 
  y <- bc_triple2_tidy %>% filter(term1 %in% x$term2 | term1 %in% q$term2 | term1 %in% triple2) 
  z <- bc_triple1_tidy %>% filter(term1 %in% y$term2 |term1 %in% x$term2 | term1 %in% q$term2 | term1 %in% triple1) 
  
  # any trait triples in any BC triples
  shared <- bind_rows(x,y,z, q) %>%  inner_join(trait_twostep_triples, by = c("term1"="term1", "term2"="term2")) 
  
  
  # full 
  a <-bind_rows(x,y,z, q) %>% mutate(group = "BC") %>% 
    bind_rows(., trait_twostep_triples %>% mutate(group = "trait")) %>% distinct() %>% 
    left_join(shared %>% mutate(group2 = "shared"), by = c("term1"="term1", "term2"="term2")) %>% 
    mutate(group = ifelse(!is.na(group2), group2, group)) %>% 
    select(term1, term2, n, group) %>% distinct() %>% 
    mutate(group = as.factor(group)) %>% distinct()
  
  a_sankey<- a %>% 
    filter(!(!term2 %in% term1 & term2 != outcome)) %>% 
    filter(!(!term2 %in% term1 & term2 != outcome))# exclude loose terms
  
  full_sankey<- make_sankey(a_sankey, fontSize=13, colour_links = T)
  
  # exclude  links with specified low n
  
  print(paste0("Using n=", n_filter))
  #to_exl<- a %>%  filter(n <= n_filter & term2 != KEY_TERM) %>%  pull(term2) %>% unique()
  a_sub <- a %>% filter(n > n_filter | term2 == outcome) #%>% filter(!term1 %in% to_exl) 
  
  remove_loose_terms1 =T
  while (remove_loose_terms1) {
    print("term1 while loop filtering")
    a_sub <- a_sub %>% filter(!(!term1 %in% term2 & term1 != KEY_TERM)) 
    tmp <- a_sub %>% filter((!term1 %in% term2 & term1 != KEY_TERM)) 
    if (dim(tmp)[1]==0){
      remove_loose_terms1 =F
    }
  }
  
  a_sub <- a_sub %>% filter(term2 != KEY_TERM)
  
  remove_loose_terms2 =T
  while (remove_loose_terms2) {
    print("term 2 while loop filtering")
    a_sub <- a_sub %>% filter(!(!term2 %in% term1 & term2 != outcome))
    tmp <- a_sub %>% filter((!term2 %in% term1 & term2 != outcome))
    if (dim(tmp)[1]==0){
      remove_loose_terms2 =F
    }
  }
  print(dim(a_sub))
  
  subset_sankey<- make_sankey(a_sub, fontSize=sankey_font, colour_links = T)
  
  
  # get full raw data
  
  terms_list <- unique(c(a_sankey$term1, a_sankey$term2))
  
  return(list(full_sankey = full_sankey,
              subset_sankey = subset_sankey,
              full_sankey_data = a_sankey,
              subset_sankey_data = a_sub,
              terms_list = terms_list))
}




overlap_lifestyle_trait_and_bc <- function(trait_triples,  bc_triples,sankey_font = 13, n_filter = 1 , outcome = "Breast cancer"){
  
  bc_triple1_tidy = bc_triples$triple1
  bc_triple2_tidy = bc_triples$triple2
  bc_triple3_tidy = bc_triples$triple3
  bc_triple4_tidy = bc_triples$triple4
  
  
  bc_triples<-
    bind_rows(
      bc_triple1_tidy %>% mutate(triple = 1),
      bc_triple2_tidy %>% mutate(triple = 2), 
      bc_triple3_tidy %>% mutate(triple = 3),
      bc_triple4_tidy %>% mutate(triple = 4)
    )
  
  
  # any trait triples in any BC triples
  shared <- bc_triples %>%  inner_join(trait_triples, by = c("term1"="term1", "term2"="term2")) %>% select(term1, term2) %>% distinct()
  
  # any trait term2 are term 1 in bc?
  tr_term2_is_bc_term1 <- intersect(trait_triples$term2, bc_triples$term1) 
  
  # get trait triples that can be linked to bc (except via obesity)
  trait_triples_linked_to_bc <- 
    trait_triples %>%
    filter(term2 %in% tr_term2_is_bc_term1) %>% 
    filter(n>0) %>% 
    filter(term2!='Obesity')
  
  # also collect trait triples linked to triples collected to BC space
  terms_A <- trait_triples %>% filter(term2 %in% tr_term2_is_bc_term1) %>% pull(term1)
  triples_X_A <- trait_triples %>% filter(term2 %in% terms_A) %>% 
    filter(!(grepl('Obesity', term2, ignore.case = T) & n <= 3)) 
  
  
  ## get BC chains of triples connected to the extracted trait terms; 
  bc_triples_linked <- bc_triples %>%  filter(term1 %in% tr_term2_is_bc_term1)
  
  bc_triple4_tidy_sub <- bc_triple4_tidy %>% right_join(bc_triples_linked %>% filter(triple == 4), by = c("term1"="term1", "term2"="term2", "n"="n"))
  bc_triple3_tidy_sub <- bc_triple3_tidy %>% right_join(bc_triples_linked %>% filter(triple == 3), by = c("term1"="term1", "term2"="term2", "n"="n"))
  bc_triple2_tidy_sub <- bc_triple2_tidy %>% right_join(bc_triples_linked %>% filter(triple == 2), by = c("term1"="term1", "term2"="term2", "n"="n"))    
  bc_triple1_tidy_sub <- bc_triple1_tidy %>% right_join(bc_triples_linked %>% filter(triple == 1), by = c("term1"="term1", "term2"="term2", "n"="n"))
  
  # make sure all subsequesnt steps are collected 
  bc_triple3_tidy_sub <- bind_rows(bc_triple3_tidy_sub,
                                   bc_triple3_tidy %>% filter(term1 %in% bc_triple4_tidy_sub$term2))
  bc_triple2_tidy_sub <- bind_rows(bc_triple2_tidy_sub,
                                   bc_triple2_tidy %>% filter(term1 %in% bc_triple3_tidy_sub$term2))
  bc_triple1_tidy_sub <- bind_rows(bc_triple1_tidy_sub,
                                   bc_triple1_tidy %>% filter(term1 %in% bc_triple2_tidy_sub$term2))
  
  bc_triples_subset <- bind_rows(bc_triple4_tidy_sub,
                                 bc_triple3_tidy_sub,
                                 bc_triple2_tidy_sub,
                                 bc_triple1_tidy_sub) %>% select(-triple) %>% distinct()
  
  # join trait and BC triples
  trait_tr_and_bc <- bind_rows(trait_triples_linked_to_bc %>% mutate(group = "trait"), 
                               bc_triples_subset %>% mutate(group = outcome),
                               triples_X_A %>% mutate(group = "trait")) %>% 
    left_join(shared %>% mutate(group2 = "shared"), by = c("term1"="term1", "term2"="term2")) %>% 
    mutate(group = ifelse(!is.na(group2), group2, group)) %>% 
    select(term1, term2, n, group) %>% distinct() %>% 
    mutate(group = as.factor(group)) %>% distinct()
  
  # keep the one with hiher n
  trait_tr_and_bc <- 
    trait_tr_and_bc %>% 
    group_by(term1, term2) %>% 
    slice(which.max(n)) %>%
    ungroup() 
  
  
  # need to drop reverse connections: keep A-B or B-A depending which one is more common
  trait_tr_and_bc<- trait_tr_and_bc %>% 
    #filter(term1 %in% c("cytokine", "estrogens") & term2 %in% c("cytokine", "estrogens")) %>% 
    # forward count
    rename(n_pair_f = n) %>% 
    # join col that will show reverse pair count
    left_join(trait_tr_and_bc %>% select(term1,term2,n_pair_b = n),
              by = c('term1' = 'term2', 'term2' = 'term1')) %>% 
    distinct() %>% 
    # if rel does not exist uin reverse, set it to 0
    mutate(across(n_pair_b, ~replace_na(.x, 0))) %>% 
    # if f more common, keep it, else, keep reverse
    mutate(keep = ifelse(n_pair_f >= n_pair_b ,T,F)) %>% 
    filter(keep==T) %>%  rename(n=n_pair_f)
  
  
  remove_loose_terms2 =T
  while (remove_loose_terms2) {
    print("term 2 while loop filtering")
    trait_tr_and_bc <- trait_tr_and_bc %>% filter(!(!term2 %in% term1 & term2 != outcome))
    tmp <- trait_tr_and_bc %>% filter((!term2 %in% term1 & term2 != outcome))
    if (dim(tmp)[1]==0){
      remove_loose_terms2 =F
    }
  }
  
  
  if (n_filter > 1) {
    ## filtering by n
    print(paste0("Using n=", n_filter))
    a_sub <- trait_tr_and_bc %>% filter(n >= n_filter | term2 == outcome) 
    
    remove_loose_terms1 =T
    while (remove_loose_terms1) {
      print("term1 while loop filtering")
      a_sub <- a_sub %>% filter(!(!term1 %in% term2 )) 
      tmp <- a_sub %>% filter((!term1 %in% term2 )) 
      if (dim(tmp)[1]==0){
        remove_loose_terms1 =F
      }
    }
    
    
    remove_loose_terms2 =T
    while (remove_loose_terms2) {
      print("term 2 while loop filtering")
      a_sub <- a_sub %>% filter(!(!term2 %in% term1 & term2 != outcome))
      tmp <- a_sub %>% filter((!term2 %in% term1 & term2 != outcome))
      if (dim(tmp)[1]==0){
        remove_loose_terms2 =F
      }
    }
    
    
    print(dim(a_sub))
  } else{
    a_sub <-tibble()
  }
  
  
  
  
  # get full raw data
  
  terms_list <- unique(c(trait_tr_and_bc$term1, trait_tr_and_bc$term2))
  
  return(list(sankey_data = trait_tr_and_bc,
              sankey_data_filtered = a_sub,
              terms_list = terms_list))
}




extract_two_triples_for_trait <- function(trait_tidy, KEY_TERM,  ignore_terms = c() ){#, keep_only = c()){
  
  trait_tidy <- trait_tidy %>% 
    mutate(term1.type_verbose  = ifelse(grepl("obesity", term1, ignore.case = T), 'antro', term1.type_verbose)) %>%
    mutate(term2.type_verbose  = ifelse(grepl("obesity", term2, ignore.case = T), 'antro', term2.type_verbose)) 
  
  ## triple 1
  
  trait_triple1 <- trait_tidy %>% 
    filter(term1 %in% KEY_TERM ) %>% 
    filter(term2.type_verbose != 'disease') %>% 
    filter(!predicate %in% c("COEXISTS_WITH", 'NEG_ISA', 'NEG_COEXISTS_WITH')) %>%
    filter(term2.type_verbose != "drug_or_compound") %>% 
    select(1:5) %>% 
    group_by(term1, predicate, term2) %>% 
    summarise_all(sum) %>% ungroup()
  
  print("Triple 1")
  print(trait_triple1 %>% select(term1, term2) %>% distinct() %>% dim())
  
  trait_triple1_tidy <- trait_triple1 %>% 
    tidy_terms_for_viz() %>% 
    select(term1,term2,n_pair) %>% 
    filter(term1 !=term2) %>% 
    group_by(term1, term2) %>% 
    slice(which.max(n_pair)) %>%
    ungroup()  %>% 
    rename(n=n_pair)# %>% 
  #filter(!term2 %in% ignore_terms) #%>% 
  #filter(term2 %in% keep_only)
  
  
  
  print("Triple 1 tidy")
  print(trait_triple1_tidy %>% select(term1, term2) %>% distinct() %>% dim())
  
  
  ## triple 2
  trait_triple2 <- trait_tidy %>% 
    filter(term1 %in% trait_triple1$term2) %>% 
    filter(term2.type_verbose != 'disease') %>% 
    filter(!predicate %in% c("COEXISTS_WITH", 'NEG_ISA', 'NEG_COEXISTS_WITH')) %>%
    filter(term2.type_verbose != "drug_or_compound") %>% 
    filter(!term1 %in% KEY_TERM) %>% 
    filter(!term2 %in% KEY_TERM)
  
  print("Triple 2")
  print(trait_triple2 %>% select(term1, term2) %>% distinct() %>% dim() )
  
  # need to drop reverse connections: keep A-B or B-A depending which one is more common
  trait_triple2_onedir<- trait_triple2 %>% 
    # forward count
    rename(n_pair_f = n_pair) %>% 
    # join col that will show reverse pair count
    left_join(trait_triple2 %>% select(term1,term2,n_pair_b = n_pair), by = c('term1' = 'term2', 'term2' = 'term1')) %>% 
    distinct() %>% 
    # if rel does not exist uin reverse, set it to 0
    mutate(across(n_pair_b, ~replace_na(.x, 0))) %>% 
    # if f more common, keep it, else, keep reverse
    mutate(keep = ifelse(n_pair_f >= n_pair_b ,T,F)) %>% 
    filter(keep==T) %>%  rename(n_pair=n_pair_f)
  
  trait_triple2_tidy <- trait_triple2_onedir %>% 
    tidy_terms_for_viz() %>% 
    # after tidying names you get almost duplicates: keep the one with highest n
    select(term1,term2,n_pair) %>% 
    filter(term1 !=term2) %>% 
    group_by(term1, term2) %>% 
    slice(which.max(n_pair)) %>%
    ungroup()  %>% 
    rename(n=n_pair) 
  
  print("Triple 2 tidy")
  print(trait_triple2_tidy %>% select(term1, term2) %>% distinct() %>% dim() )
  
  
  
  ## joining triples
  
  if (dim(trait_triple2_tidy)[1] != 0){
    
    # pick terms linked to anchor
    trait_x<- trait_triple1_tidy %>% filter(n>0)
    # make sure they are term1 in triple2 and not in term2 (so that it does not create loops and levels)
    trait_y<- trait_triple2_tidy %>% filter(term1 %in% trait_x$term2 & !term2 %in% trait_x$term2)  %>% filter(n>0)
    # drop triple 1 term2 that ends up not linked to anything because of previous filteting step
    trait_x2 <- trait_x %>%  filter(term2 %in% trait_y$term1)
    
    trait_twostep_triples<- bind_rows(trait_x2, trait_y)
    
    print("Joined triples")
    print(trait_twostep_triples %>% select(term1, term2) %>% distinct() %>% dim() )
  }else{
    trait_twostep_triples = trait_triple1_tidy
  }
  
  
  return(list(triple1_tidy =trait_triple1_tidy,
              triple2_tidy =trait_triple2_tidy,
              joined_triples = trait_twostep_triples))
  
}


extract_lifestyle_main_triples <- function(trait_tidy ){
  
 
  trait_tidy <- trait_tidy %>% 
    mutate(term1.type_verbose  = ifelse(grepl("obesity", term1, ignore.case = T), 'antro', term1.type_verbose)) %>%
    mutate(term2.type_verbose  = ifelse(grepl("obesity", term2, ignore.case = T), 'antro', term2.type_verbose)) 
  
  
  ## triple 1
  
  trait_triple1 <- trait_tidy %>% 
    filter(term2.type_verbose != 'disease') %>% filter(term1.type_verbose != 'disease') %>% 
    #filter(!predicate %in% c("COEXISTS_WITH", 'NEG_ISA', 'NEG_COEXISTS_WITH')) %>%
    filter(term2.type_verbose != "drug_or_compound") %>% filter(term1.type_verbose != "drug_or_compound") %>% 
    select(1:5) %>% 
    group_by(term1, predicate, term2) %>% 
    summarise_all(sum) %>% ungroup()
  
  print(trait_triple1 %>% select(term1, term2) %>% distinct() %>% dim())
  
  trait_triple1_tidy <- trait_triple1 %>% 
    tidy_terms_for_viz() %>% 
    select(term1,term2,n_pair) %>% 
    filter(term1 !=term2) %>% 
    group_by(term1, term2) %>% 
    slice(which.max(n_pair)) %>%
    ungroup()  %>% 
    rename(n=n_pair)
  
  
  print("Triple 1 tidy")
  print(trait_triple1_tidy %>% select(term1, term2) %>% distinct() %>% dim())
  
  
  
  return(trait_triple1_tidy)
  
  
}







make_sankey <- function(links, fontSize=10, colour_links = F, height=NULL,
                        term_nodes_col = 'grey',
                        shared_col = '#C54F1D', # dark orange
                        trait_col = "#89B6EE", # blue
                        bc_col = "#E6B07F"){  # orange
  
  links <- links %>% rename(source = term1, target = term2, value = n)
  
  nodes <- data.frame(
    name=c(as.character(links$source), as.character(links$target)) %>% 
      unique()
  )
  
  # With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
  links$IDsource <- match(links$source, nodes$name)-1 
  links$IDtarget <- match(links$target, nodes$name)-1
  
  if (colour_links){
    # Add a 'group' column to each node. Here I decide to put all of them in the same group to make them grey
    nodes$group <- as.factor(c("term_node"))
    # Give a color for each group:
    my_color <- paste0('d3.scaleOrdinal() .domain(["trait", "BC", "shared", "term_node"]) .range(["',trait_col ,'", "',bc_col ,'", "',shared_col ,'", "',term_nodes_col ,'"])')
    
    # Make the Network
    p <- sankeyNetwork(Links = links, Nodes = nodes,
                       Source = "IDsource", Target = "IDtarget",
                       Value = "value", NodeID = "name", 
                       fontSize = fontSize,
                       fontFamily = 'sans-serif',
                       colourScale=my_color,
                       LinkGroup="group", 
                       NodeGroup="group",
                       sinksRight=T,
                       height=height)
  } else{ 

  # Make the Network without colouring links
  p <- sankeyNetwork(Links = links, Nodes = nodes,
                     Source = "IDsource", Target = "IDtarget",
                     Value = "value", NodeID = "name", 
                     fontSize = fontSize,
                     fontFamily = 'sans-serif',
                     sinksRight=FALSE,
                     height=height)
  }
  
  p
  
}





select_network <- function(full_data, key_term, network_dir,
                           lowest_count_rel_to_keep = 10,
                           exclude_neutral = F, 
                           include_common_nodes_links = T, 
                           common_nodes = '', 
                           predicate_ref, 
                           rare_node = F){
  
  #    *Network matching methods:*
  #      
  #    * Forwards (2-level, 1 direction): term1 (set) -> term2 (many) -> term1 (many more)
  #    * Backwards (2-level, 1 direction): term2 (many more) <- term 1 (many) <- term2 (set)
  #    
  #    * F+B (1-level, 2 directions) : /[             \[ term1 /]-> term2  \]
  #    /[  term1   <- \[ term2 /]          \]
  #    
  #    
  #    * F+B pair ( two terms match + 1-level, 2 directions) :
  #      termX -> termY -> many 
  #    many <- termX <- termY
  #    
  #    
  #    *Additional network filters:*
  #      
  #    - exclude neutral relationships
  #    - exclude low count relationships (set threshold, default is 10)
  #    - do not included rels from common nodes (e.g. estrogen - get a hairball right away)
  #    - for 'rare' nodes - search the term in both term1 and term2 - use F+B
  
  
  if (!network_dir %in% c('F+B pair', 'F+B multiple') ){
    if (length(key_term) > 1) {
      stop("this network type accepts only 1 key word")
    }
  }
  
  # select the direction of the network   
  if (network_dir == 'forwards'){  
    data_level1<- full_data %>% filter(term1 == key_term) 
    data_level2<- full_data %>% filter(term1 %in% data_level1$term2) 
    
    data_megred <- bind_rows(data_level1, data_level2) %>% distinct()
    if (!include_common_nodes_links){
      print('removing rels of common nodes')
      data_megred<-data_megred %>% filter(!term1 %in% common_nodes)
    }
    
  } else if (network_dir == 'backwards'){  
    data_level1<- full_data %>% filter(term2 == key_term) 
    data_level2<- full_data %>% filter(term2 %in% data_level1$term1) 
    
    data_megred <- bind_rows(data_level1, data_level2) %>% distinct()
    if (include_common_nodes_links == F){
      print('removing rels of common nodes')
      data_megred<-data_megred %>% filter(!term2 %in% common_nodes)
    }
  }  else if (network_dir == 'F+B 1-level'){ 
    data_dir1<- full_data %>% filter(term1 == key_term) 
    data_dir2<- full_data %>% filter(term2 == key_term) 
    
    data_megred <- bind_rows(data_dir1, data_dir2) %>% distinct()
    if (!include_common_nodes_links){
      print('removing rels of common nodes')
      data_megred<-data_megred %>% filter(!term1 %in% common_nodes)
    }
  }  else if (network_dir == 'F+B 2-level'){ 
    data_dir1_level1 <- full_data %>% filter(term1 == key_term) 
    data_dir1_level2 <- full_data %>% filter(term1 %in% data_dir1_level1$term2) 
    data_dir2_level1<- full_data %>% filter(term2 == key_term) 
    data_dir2_level2<- full_data %>% filter(term2 %in% data_dir2_level1$term1) 
    
    data_megred <- bind_rows(data_dir1_level1, data_dir1_level2, data_dir2_level1, data_dir2_level2) %>% distinct()
    
    if (!include_common_nodes_links){
      print('removing rels of common nodes')
      data_megred<-data_megred %>% filter(!term1 %in% common_nodes)
    }
    
  } else if (network_dir == 'F+B pair'){ 
    stopifnot(length(key_term) == 2)
    data_dir1<- full_data %>% filter(term1 %in% key_term) 
    data_dir2<- full_data %>% filter(term2 %in% key_term) 
    
    data_megred <- bind_rows(data_dir1, data_dir2) %>% distinct()
    if (!include_common_nodes_links){
      print('removing rels of common nodes')
      data_megred<-data_megred %>% filter(!term1 %in% common_nodes)
    }
  } else if (network_dir == 'F+B multiple'){ 
    stopifnot(length(key_term) >= 2)
    data_dir1<- full_data %>% filter(term1 %in% key_term) 
    data_dir2<- full_data %>% filter(term2 %in% key_term) 
    
    data_megred <- bind_rows(data_dir1, data_dir2) %>% distinct()
    if (!include_common_nodes_links){
      print('removing rels of common nodes')
      data_megred<-data_megred %>% filter(!term1 %in% common_nodes)
    }
  }  
  
  
  
  # tidy data further
  
  if (exclude_neutral){
    print('removing neutral relationships')
    data_megred <- data_megred %>%  filter(!predicate %in% c('COEXISTS_WITH', 'INTERACTS_WITH', 'ASSOCIATED_WITH')) 
  }
  
  data_megred_tidy<- data_megred %>% 
    select(term1, predicate, term2)  %>% 
    group_by(term1, predicate, term2) %>% 
    count() %>% ungroup() %>% 
    filter_by_count(., lowest_count_rel_to_keep, key_term, rare_node )%>%
    distinct() %>% 
    left_join(predicate_ref %>% select(predicate, direction), by = 'predicate')
  
  print(paste0("Rels in the network subset: ", dim(data_megred_tidy)[1]))
  
  return (data_megred_tidy)
}




filter_by_count <- function(dat, lowest_count_rel_to_keep, key_term='', rare_node=F){
  if (!rare_node){key_term=''}
  # if it's indicate that the key term is a rare node, we will keep all links for it, even below counts threshold
  keep_by_term<-dat %>% filter(term1 %in% key_term | term2 %in% key_term)
  keep_by_count<- dat %>% filter(n >= lowest_count_rel_to_keep)
  out<-bind_rows(keep_by_term, keep_by_count)
  return(out)
}



build_networkD3 <- function(network_subset, node_counts) {
  #http://curleylab.psych.columbia.edu/netviz/netviz2.html
  
  networkData <- data.frame(src = network_subset$term1,
                            target = network_subset$term2, 
                            rel_count = network_subset$n,
                            direction = network_subset$direction, 
                            stringsAsFactors = FALSE)
  direction <- network_subset$direction
  
  # make a nodes data frame out of all unique nodes in networkData
  nodes <- data.frame(name = unique(c(networkData$src, networkData$target)))
  
  # make a group variable where nodes in networkData$src are identified
  nodes$group <- nodes$name %in% networkData$src
  
  # add count how often it appers in network as size
  nodes<- nodes %>% left_join(node_counts, by='name') %>% select(-group) %>% rename('group' =  'type_verbose')
  
  # work out order of occurrence and set colors order
  node_cat_order <- data.frame(cat = nodes$group[!duplicated(nodes$group)])
  category_cols <- data.frame(cat = c("drug_or_compound", "disease", "key_term" , "any" ),
                              col = c("'#F2AD00'", "'#00A08A'", "'green'", "'black'"))
  cols_order <- left_join(node_cat_order, category_cols, by ='cat') %>% pull(col) %>% str_c(.,  collapse = ", ")
  JS_input <- str_c('d3.scaleOrdinal([',cols_order,']);') 
  if(is.na(JS_input)) {JS_input <- "d3.scaleOrdinal(['green','black']);" }
  
  # make a links data frame using the indexes (0-based) of nodes in 'nodes'
  links <- data.frame(source = match(networkData$src, nodes$name) - 1,
                      target = match(networkData$target, nodes$name) - 1)
  links$value <- networkData$rel_count
  
  
  b<-forceNetwork(Links = links, Nodes = nodes, Source = "source",
                  Target = "target", NodeID ="name", Group = "group",
                  opacity = 0.8, opacityNoHover = 0.5, zoom = T,
                  #colourScale = ifelse(nodes[1,"group"] == 'drug_or_compound',
                  #                               JS('d3.scaleOrdinal([ "#F2AD00" ,"black","#00A08A", "green"]);'),
                  #                               JS('d3.scaleOrdinal([ "black","#F2AD00" ,"#00A08A", "green"]);')),
                  colourScale = JS(JS_input),
                  Value = 'value', arrows = T,
                  Nodesize = 'size',
                  charge = -1000, # node repulsion
                  linkDistance = 25,
                  fontSize=24,
                  linkColour = ifelse(direction == 'positive', "#FF0000", 
                                      ifelse(direction == 'negative', "#5BBCD6", 'black'))
  )
  return(b)
}