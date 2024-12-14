library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(TwoSampleMR)

# file from USCS http://genome.ucsc.edu/cgi-bin/hgTables for protein_names.csv as input pasted

locations_all <- read_csv("01_MR_related/results/mr_evidence_outputs/protein_locations_ucsc.csv") %>% 
  filter(!grepl("alt|fix|Y", chrom)) %>% select(-1, -score) %>% distinct()

location_tidy <- tibble()

for (i in unique(locations_all$name2)){
  print(i)
  tmp <- locations_all  %>%  filter(name2 ==  i)
  start<-min(tmp$txStart)
  end<-max(tmp$txEnd)
  chrom<-gsub("chr", "", unique(tmp$chrom))
  out<-tibble(gene = i, chr=chrom, posStart=start, posEnd=end)
  #print(out)
  location_tidy<-bind_rows(location_tidy, out)
}

# get the full list
protein_names <- read_csv("01_MR_related/results/mr_evidence_outputs/protein_names.csv", col_names = c('name', 'gene')) %>% distinct()
names_tidy <- read_tsv("01_MR_related/results/mr_evidence_outputs/protein_in_final_setV3.tsv") %>% distinct()
protein_names_id <- left_join(names_tidy, protein_names , by =c("exposure"="name"))
protein_names_id<- protein_names_id %>% filter(exposure != "Albumin")

location_tidy_names <- right_join(location_tidy, protein_names_id)

# re-extract instruments
inst_list <- list()
for (id in unique(protein_names_id$id.exposure)){
  print(id)
  inst_list[[id]] <- extract_instruments(id)
}

# check total inst, same chr inst
location_tidy_names$total_inst <- NA
location_tidy_names$inst_same_chr <- NA
location_tidy_names$within_1Mb_gene_region <- 0
for (i in 1:length(location_tidy_names$id.exposure)){
  
  id <- location_tidy_names$id.exposure[i]
  print(id)
  total_inst <- nrow(inst_list[[id]])
  
  cis_chr <- location_tidy_names$chr[i]
  
  inst_same_chr <- length(which(inst_list[[id]]$chr.exposure == cis_chr))
  print(paste0("Same chr matches: ",inst_same_chr ))
  
  within_region=0
  if (inst_same_chr > 0){
    pos_snp_same_chr <-  inst_list[[id]] %>% filter(chr.exposure == cis_chr) %>% pull(pos.exposure)
    print(paste("at postitions:", pos_snp_same_chr))
    
    start_1Mb <- location_tidy_names$posStart[i]-1000000
    end_1Mb <- location_tidy_names$posEnd[i]+1000000 
    
    for (j in pos_snp_same_chr){
      if (j >  start_1Mb & j < end_1Mb){ # within 1Mb
        within_region = within_region +1
      }
    }
  }
  
  location_tidy_names$total_inst[i] <- total_inst
  location_tidy_names$inst_same_chr[i] <- inst_same_chr
  location_tidy_names$within_1Mb_gene_region[i] <- within_region
}


write_tsv(location_tidy_names, "01_MR_related/results/mr_evidence_outputs/proteins_gene_region_SNPs.tsv") 


sensitivity_hazards <- read_tsv("01_MR_related/results/mr_evidence_outputs/all_data_with_sens_filters_hazardsV3.tsv") %>% 
    select(id.exposure=exposure.id, outcome, nsnp, egger_intercept:hazards) %>% 
    distinct() %>% filter(outcome == "BCAC'17") %>% 
    mutate(hazards = ifelse(is.na(hazards), "", hazards))

location_tidy_names_w_haz <- left_join(location_tidy_names, sensitivity_hazards, by = "id.exposure")

write_tsv(location_tidy_names_w_haz, "01_MR_related/results/mr_evidence_outputs/proteins_gene_region_SNPs_w_hazards.tsv") 


# also saving all protein instruments

inst_df<- bind_rows(inst_list) %>% 
  select(id.exposure, exposure, SNP, chr.exposure, pos.exposure, everything()) %>% 
  select(-mr_keep.exposure, -pval_origin.exposure, -data_source.exposure)
  
write_tsv(inst_df, "01_MR_related/results/mr_evidence_outputs/proteins_instruments.tsv") 

