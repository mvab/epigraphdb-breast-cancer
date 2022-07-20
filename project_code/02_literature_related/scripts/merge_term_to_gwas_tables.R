
## adhoc script for extracting previously mapped lit terms to opengwas  IDs


prev_terms1<- read_csv("../results/literature_outputs/sankey_terms_storage/terms_to_opengwas_part1.csv")
prev_terms2<- read_csv("../results/literature_outputs/sankey_terms_storage/terms_to_opengwas_part2.csv")
prev_terms3<- read_csv("../results/literature_outputs/sankey_terms_storage/terms_to_opengwas_part3.csv")
prev_terms4<- read_csv("../results/literature_outputs/sankey_terms_storage/terms_to_opengwas_part4.csv")
prev_terms5<- read_csv("../results/literature_outputs/sankey_terms_storage/terms_to_opengwas_part5.csv")
prev_terms6<- read_csv("../results/literature_outputs/sankey_terms_storage/terms_to_opengwas_part6.csv")
prev_terms7<- read_csv("../results/literature_outputs/sankey_terms_storage/terms_to_opengwas_part7.csv")

prev_terms <- bind_rows(prev_terms7,prev_terms6, prev_terms5, prev_terms4, prev_terms3, prev_terms2, prev_terms1) %>% distinct()

new_trait_terms <- read_csv("./results/literature_outputs/sankey_terms_storage/lit_terms_menopause.csv")

merged<- new_trait_terms %>% left_join(prev_terms) %>% arrange(value)

# coalese rows
merged <- merged %>%
  group_by(value) %>%
  fill(everything(), .direction = "downup") %>%
  slice(1)


merged %>% write_csv("./results//literature_outputs/sankey_terms_storage/lit_terms_menopause.csv")




library(TwoSampleMR)
ao <- available_outcomes()
ao <- ao %>% filter(population == 'European')
