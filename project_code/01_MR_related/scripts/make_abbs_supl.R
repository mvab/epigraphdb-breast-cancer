# make abbs supl table


trits_202 <- read_tsv("01_MR_related//results/mr_evidence_outputs/passed_multiple_testingV3.tsv") %>%  select(exposure_cat,id.exposure, exposure) %>% distinct()



key <- read_csv("01_MR_related/scripts/app2_heatmaps_app/data/renaming_key_tidy.csv") %>% 
  rename(exposure_new_name = exposure, exposure=exposure.trait, id.exposure =exposure.id  ) %>% select(-exposure, -exposure_cat_sub)


dat<- left_join(trits_202, key, by="id.exposure") %>% 
      filter(!is.na(exposure_cat)) %>%  # removes some silly proteins
      filter(!(id.exposure =="met-c-841" & exposure == "Albumin"))

dat %>% count(id.exposure) %>% View()


# now add gene names for proteins

proteins <- read_csv("01_MR_related//results/mr_evidence_outputs/protein_names_w_ids.csv") %>% select(gene, id.exposure=exposure.id)

dat <- left_join(dat,proteins) %>% arrange(exposure_cat)

write_csv(dat, "01_MR_related//results/mr_evidence_outputs/abbs_and_genes_suppl.csv")
