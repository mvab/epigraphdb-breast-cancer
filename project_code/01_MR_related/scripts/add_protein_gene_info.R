library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(TwoSampleMR)


##### ideally don't need to do this step if protein_names already contains ids

protein_names <- read_csv("01_MR_related/results/mr_evidence_outputs/protein_names.csv", col_names = c('name', 'gene')) %>% distinct()
# or other file with protein name, gene (and id?)

names_tidy <- read_tsv("01_MR_related/results/mr_evidence_outputs/protein_in_final_setV3.tsv") %>% distinct()
traits <- read_csv("01_MR_related/results/mr_evidence_outputs/all_mreve_bc_resultsV3.csv") %>% 
  select(exposure.id, exposure.trait) %>% distinct() %>% 
  create_exposure_categories() %>% select(-exposure.trait) %>% 
  filter(exposure_cat == "Proteins")

protein_names_id <- full_join(traits, protein_names , by =c("exposure"="name")) %>% distinct()
protein_names_id<- protein_names_id %>% 
  filter(!grepl("albumin",exposure, ignore.case=T)) %>% 
  mutate(exposure.id = ifelse(exposure == "Adiponectin", "ieu-a-1", exposure.id)) %>% 
  filter(!grepl("raw", exposure.id)) %>% 
  arrange(gene)

write_csv(protein_names_id, "01_MR_related/results/mr_evidence_outputs/protein_names_w_ids.csv")

####


#### START HERE

protein_names_id <- read_csv("01_MR_related/results/mr_evidence_outputs/protein_names_w_ids.csv") %>% filter(!is.na(gene))

# file from USCS http://genome.ucsc.edu/cgi-bin/hgTables for protein_names.csv or protein_names_w_ids.csv  as input pasted
locations_all <- read_csv("01_MR_related/results/mr_evidence_outputs/protein_locations_ucsc.csv") %>% 
  filter(!grepl("alt|fix|Y|X", chrom)) %>% select(-1, -score) %>% distinct()

# tidy UCSC output
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


location_tidy_names <- right_join(location_tidy, protein_names_id, by="gene") %>% distinct()

# when redoing: check NA; 
# most will be chr X: FLNA, IL3RA, KLHL13 - so just leave it,
#  but for CCL3L1 and KIR2DL5A manually check and update UCSC file

write_csv(location_tidy_names, "01_MR_related/results/mr_evidence_outputs/protein_gene_regions_ids.csv")



#### when working on literature case studies
#### 

# 1: HDL

# file from USCS http://genome.ucsc.edu/cgi-bin/hgTables 
# gene IDS pasted from  02_literature_related/results/literature_outputs/sankey_terms_storage/lit_terms_HighDensityLipoproteinsNEW_genes.csv

locations_all <- read_csv("01_MR_related/results/case_study_reports_tidy/literature_proteins_ucsc_data/HDL_ucsc_protein_data.csv") %>% 
  filter(!grepl("alt|fix|Y|X", hg38.knownGene.chrom)) %>% select(-1) %>% distinct() %>% arrange(hg38.kgXref.geneSymbol)

exls<- read_csv("01_MR_related/results/case_study_reports_tidy/literature_proteins_ucsc_data/HDL_ucsc_protein_data.csv") %>% 
  filter(grepl("alt|fix|Y|X", hg38.knownGene.chrom)) %>% filter(!hg38.kgXref.geneSymbol %in% locations_all$hg38.kgXref.geneSymbol)

colnames(locations_all) <- c("chrom", 'txStart', 'txEnd', 'name2')


# tidy UCSC output
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

protein_names_id <- read_csv("02_literature_related/results/literature_outputs/sankey_terms_storage/lit_terms_HighDensityLipoproteinsNEW_genes.csv")

location_tidy_names <- left_join(protein_names_id,location_tidy,  by="gene") %>% distinct()


write_csv(location_tidy_names, "02_literature_related/results/literature_outputs/sankey_terms_storage/lit_terms_HighDensityLipoproteinsNEW_genes.csv")



# 1: cBMI

# file from USCS http://genome.ucsc.edu/cgi-bin/hgTables 
# gene IDS pasted from  02_literature_related/results/literature_outputs/sankey_terms_storage/lit_terms_HighDensityLipoproteinsNEW_genes.csv

locations_all <- read_csv("01_MR_related/results/case_study_reports_tidy/literature_proteins_ucsc_data/cBMI_ucsc_protein_data.csv") %>% 
  filter(!grepl("alt|fix|Y|X", hg38.knownGene.chrom)) %>% select(-1) %>% distinct() %>% arrange(hg38.kgXref.geneSymbol)

exls<- read_csv("01_MR_related/results/case_study_reports_tidy/literature_proteins_ucsc_data/cBMI_ucsc_protein_data.csv") %>% 
  filter(grepl("alt|fix|Y|X", hg38.knownGene.chrom)) %>% filter(!hg38.kgXref.geneSymbol %in% locations_all$hg38.kgXref.geneSymbol)

colnames(locations_all) <- c("chrom", 'txStart', 'txEnd', 'name2')


# tidy UCSC output
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

protein_names_id <- read_csv("02_literature_related/results/literature_outputs/sankey_terms_storage/lit_terms_childhood_obesity_NEW_genes.csv")

location_tidy_names <- left_join(protein_names_id,location_tidy,  by="gene") %>% distinct()


write_csv(location_tidy_names, "02_literature_related/results/literature_outputs/sankey_terms_storage/lit_terms_childhood_obesity_NEW_genes.csv")


