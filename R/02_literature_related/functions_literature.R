### functions
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


keep_one_type <- function(terms_w_multiple_types, node_types_full){
  
  node_types_count_w_mulp<-node_types_full %>% 
    filter(name %in% terms_w_multiple_types) %>% 
    group_by(name, type) %>% 
    count() %>% ungroup()
  
  selected_types <- tibble()
  for (i in terms_w_multiple_types){
    
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