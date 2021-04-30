
process_bc_outcomes <- function(dat){
  
  # deal with al outcome related changes
  dat <- dat %>%  
    # create groups
    mutate(outcome = case_when(
    outcome.id %in% c('ieu-a-1127','ieu-a-1132', 'ieu-a-1133', 'ieu-a-1134', 'ieu-a-1167', 'ieu-a-1164', 'ieu-a-1161') ~ 'ER+ postmeno',
    outcome.id %in% c('ieu-a-1128','ieu-a-1135', 'ieu-a-1136', 'ieu-a-1137', 'ieu-a-1166', 'ieu-a-1163', 'ieu-a-1160') ~ 'ER- premeno',
    outcome.id %in% c('ieu-a-1126','ieu-a-1129', 'ieu-a-1130', 'ieu-a-1131', 'ieu-a-1168', 'ieu-a-1165', 'ieu-a-1162', 'ebi-a-GCST004988', 'ebi-a-GCST007236') ~ 'Breast cancer (all)',
    grepl('ukb', outcome.id) ~ 'ER+ postmeno UKB')) %>% 
    # create data version subgroups
    mutate(chip = case_when(
      outcome.id %in% c('ieu-a-1126','ieu-a-1127', 'ieu-a-1128', 'ebi-a-GCST004988') ~ 'Meta',
      outcome.id %in% c('ieu-a-1129', 'ieu-a-1132','ieu-a-1135') ~                     'OncArray',
      outcome.id %in% c('ieu-a-1130', 'ieu-a-1133','ieu-a-1136') ~                     'iCOG2017',
      outcome.id %in% c('ieu-a-1161', 'ieu-a-1160','ieu-a-1162', 'ebi-a-GCST007236') ~ 'iCOG2015',
      outcome.id %in% c('ieu-a-1131', 'ieu-a-1134','ieu-a-1137') ~                     'GWASold1',
      outcome.id %in% c('ieu-a-1166', 'ieu-a-1167','ieu-a-1168') ~                     'GWASold2',
      outcome.id %in% c('ieu-a-1163', 'ieu-a-1164','ieu-a-1165') ~                     'Survival_',
      grepl('ukb', outcome.id) ~                                                       'UKBB')) %>% 
    
    # create outcome label
    mutate(outcome.details =  paste0(outcome.id, " (",format(outcome.sample_size, big.mark = ","))) %>% 
    separate(outcome.details, into=c('outcome.details', 'tmp'), sep=",") %>% 
    #mutate(outcome.details = gsub("ieu-a-|ukb-a-|ukb-b-|ukb-d-|ebi-a-GCST00", '', outcome.details)) %>% 
    mutate(outcome.details = paste0(chip, ": ",outcome.details, "K)")) 
    
    # make outcome.details to be factors, to control in which order they appear
    dat<- dat %>% 
      mutate(outcome.details = factor(outcome.details, levels= c( 
        sort(unique(dat$outcome.details)[grepl('Meta', unique(dat$outcome.details))]),
        sort(unique(dat$outcome.details)[grepl('OncArray', unique(dat$outcome.details))]),
        sort(unique(dat$outcome.details)[grepl('iCOG2017', unique(dat$outcome.details))]),
        sort(unique(dat$outcome.details)[grepl('iCOG2015', unique(dat$outcome.details))]),
        sort(unique(dat$outcome.details)[grepl('GWASold1', unique(dat$outcome.details))]),
        sort(unique(dat$outcome.details)[grepl('GWASold2', unique(dat$outcome.details))]),
        sort(unique(dat$outcome.details)[grepl('Survival_', unique(dat$outcome.details))]),
        sort(unique(dat$outcome.details)[grepl('UKBB', unique(dat$outcome.details))])  ))) 
  
    return(dat)
}
  
create_exposure_categories <- function(dat){
  
  # shorten names for some:
  dat <- dat %>%  
    mutate(exposure = str_replace(exposure.trait, "Treatment/medication code", "Drug"),
         exposure = str_replace(exposure, "Non-cancer illness code  self-reported", "Illness"),
         exposure = str_replace(exposure, "Diagnoses - main ICD10:", "Diagnosis"),
         exposure = str_replace(exposure, "Mineral and other dietary supplements:", "Supplements"),
         exposure = str_replace(exposure, "Vitamin and mineral supplements:", "Vitamins"),
         exposure = str_replace(exposure, "Illness  injury  bereavement  stress in last 2 years", "Stress"),
         exposure = str_replace(exposure, "Types of transport used (excluding work)", "Transport"),
         exposure = str_replace(exposure, "Types of physical activity in last 4 weeks", "Physical activity")) %>% 
    mutate(exposure = gsub("Comparative ", "", exposure)) %>% 
    mutate(exposure = gsub("\\(eg.*)", "", exposure)) %>% 
    mutate(exposure = gsub(", because of other reasons", "", exposure)) %>% 
    mutate(exposure = case_when(grepl("presbyopia", exposure) ~ "Eyesight problems1",
                                grepl("hypermetropia", exposure) ~ "Eyesight problems2",
                                TRUE ~ exposure)) %>% 
    
    ### grouping stuff into categories
    mutate(exposure_cat = case_when(
      grepl("prot", exposure.id) ~ 'protein_measures',
      grepl("met", exposure.id) ~ 'metabolite_measures',
      grepl("Drug|Medication for pain|prescription medications", exposure, ignore.case = T) ~ "drugs",
      grepl("Diagnos|Illness|cancer|neoplasm|glioma|diagnos|death|disorder|eye|carcino|colitis|disease|diabetes|asthma|sclerosis|infarction|neuroblastoma|arthritis|eczema|cholangitis|pain|hernia|Polyarthropathies", exposure, ignore.case = T) ~ "disease",
      grepl("blood|heart|Pulse", exposure, ignore.case = T) ~ "heart_blood_related",
      grepl("operation|operative|Methods of admission|Number of treatments|Spells in hospital|hospital episode", exposure, ignore.case = T) ~ "medical_operations",
      grepl("cylindrical|meridian|asymmetry|glasses|hearing|Corneal|ocular|logMAR|teeth|dental", exposure,  ignore.case = T) ~ "eye_hearing_teeth",
      grepl("vitamin|suppl", exposure, ignore.case = T) ~ "vitamin",
      grepl("waist|hip c|hip r|obesity|trunk|mass|weight|bmi|body size|height|impedance|fat percentage|body fat|Basal metabolic rate", exposure, ignore.case = T) ~ "antrophometric",
      grepl("age at|age started|parous|contraceptive pill|replacement therapy|menopause|menarche|live birth|oophorectomy|hysterectomy|menstrual", exposure, ignore.case = T) ~ "reproductive",
      grepl("alco|wine|spirits|beer", exposure, ignore.case = T) ~ "alcohol",
      grepl("smok|cigar", exposure, ignore.case = T) ~ "smoking",
      grepl("activi|transport |diy|walking|walked|Time spent|Weekly usage of|stair climbing", exposure, ignore.case = T) ~ "activity",
      grepl("sleep|Snoring|chronotype", exposure, ignore.case = T) ~ "sleep_related",
      grepl("intake|diet|food|milk|dairy|coffee|cereal|butter", exposure, ignore.case = T) ~ "diet",
      grepl("cholesterol|glyc", exposure, ignore.case = T) ~ 'metabolite_measures',
      grepl("Albumin|Apoliprotein|Adiponectin|Lipoprotein|reactive protein|Creatinine|Ferritin|Transferrin|transferase|Haemoglobin|Iron|cystatin|Testosterone|Urate|Urea|Glucose|Sodium", exposure, ignore.case = T) ~ 'other_biomarkers',
      grepl("Qualifications|GCSE|Townsend|schooling|College|intelligence|arithmetic", exposure, ignore.case = T) ~ 'education',
      grepl("anxiety|feelings|embarrassment|worr|Bulimia|depressed|guilty|Miserableness|mood|Neuroticism|unenthusiasm|tenseness|Loneliness|self-harm|Risk taking|highly strung", exposure, ignore.case = T) ~ 'psychology',
      TRUE ~ 'other')) 
  
  return(dat)
  
}


add_exposure_labels <-function(dat){
  
  dat<- dat %>% 
  # create tidy exposure sample size
  mutate(exposure.ss = format(exposure.sample_size, big.mark = ",")) %>% 
    separate(exposure.ss, into=c('exposure.ss', 'tmp'), sep=",") %>% 
    mutate(exposure.ss_label = paste0(" [",exposure.ss,"K]")) %>% 
  # add other exposure details
  mutate(exposure = case_when(
    exposure.sex == 'Females' ~ paste0(exposure, " (F) ", coalesce(consortium, author),"/" , year, exposure.ss_label),
    TRUE ~ paste0(exposure, " (M/F) ", coalesce(consortium, author),"/" , year, exposure.ss_label))) %>% 
    # add note
    mutate(exposure = case_when(grepl('Adjusted for BMI', exposure.note) ~ paste0(exposure, " AdjBMI"),
                                TRUE ~ exposure)) 
  
  return(dat)
}

tidy_display_numbers<- function(dat){
  
  dat<- dat %>% 
  mutate(loci = mr.b - 1.96 * mr.se, 
         upci = mr.b + 1.96 * mr.se,
         or = exp(mr.b), 
         or_loci = exp(loci), 
         or_upci = exp(upci),
         OR_CI = paste0(round(or,2), " [",round(or_loci,2) ,":",round(or_upci,2), "]")) %>% 
    mutate(effect_direction = ifelse(or_loci > 1 & or_upci > 1, 'positive',
                                     ifelse(or_loci < 1 & or_upci < 1, 'negative', 'overlaps null'))) %>% 
    # convert all super small pval into 1e-15
    mutate(pval_truncated = ifelse(mr.pval < 1e-15, 1e-15, mr.pval),
           #log transform pvals
           log10pval_trunc = as.integer(-log10(pval_truncated))) %>% 
    # round beta for display
    mutate(mr.b = round(mr.b, 3)) %>% 
    # tidy moe display
    mutate(`MR method and score` = paste0(mr.method," / ", mr.moescore)) 
    
  return(dat)
}




plot_bubble_plot <- function(input, font_size = 10){
  
  input <- input %>%  
   # add effect brackets
    create_beta_ranges()
 
  limit <- max(abs(input$mr.b_col)) * c(-1, 1)
  
  colourCount = length(unique(input$mr.b_col))
  getPalette<-colorRampPalette((brewer.pal(n=7, name = "PiYG")))   
  
  p<-ggplot(input, aes(x= outcome.details, y =reorder(exposure, mr.b), color= mr.b_col_verbose,  size = log10pval_trunc,
                        text = paste('Beta: ', mr.b,
                                    '</br>P-value: ', mr.pval,
                                    '</br>P-value (-log10): ', log10pval_trunc,
                                    '</br>Odds ratio: ', OR_CI,
                                    '</br>MR method and score: ', `MR method and score`,
                                    '</br>Exposure sex: ', exposure.sex)
                        
  )) + 
    geom_point()+
    scale_size(breaks = c(0,3,8,15))+
    #cale_colour_manual(values = getPalette(length(levels(input$mr.b_col_verbose))))+ 
    
    scale_colour_manual(values = getPalette(colourCount))+ 

    scale_y_discrete(position = "left")+
    facet_grid(exposure_cat~outcome, scales = 'free') +
    theme_light()+
    #theme_solarized()+
    labs(size="-log10(pval)", 
         color='beta', x="Breast cancer outcomes", y="Exposures")+
    theme(axis.text.x=element_text(angle=35, hjust=1), 
          axis.text = element_text(size = font_size),
          legend.position = 'right')

  

  return(p)
}

create_beta_ranges <- function(input){
  
  # add effect brackets
  input <- input %>% 
    mutate(mr.b_col = case_when(mr.b > 1 ~ 2,
                                                 mr.b > 0.5 ~ 1,
                                                 mr.b > 0.25 ~ 0.5,
                                                 mr.b > 0.10 ~ 0.25,
                                                 mr.b > 0.01 ~ 0.1,
                                                 mr.b > 0.001 ~ 0.01,
                                                 mr.b < -1 ~ -2,
                                                 mr.b < -0.5 ~ -1,
                                                 mr.b < -0.25 ~ -0.5,
                                                 mr.b < -0.10 ~ -0.25,
                                                 mr.b < -0.01 ~ -0.1,
                                                 mr.b < -0.001 ~ -0.01,
                                                 TRUE ~ 0)) %>% 
    mutate(mr.b_col_verbose = case_when(mr.b > 1 ~    "> 1",
                                        mr.b > 0.5 ~  "0.50:1",
                                        mr.b > 0.25 ~ '0.25:0.50',
                                        mr.b > 0.10 ~ '0.10:0.25',
                                        mr.b > 0.01 ~ '0.001:0.10',
                                        mr.b > 0.001 ~'<  0.001',
                                        mr.b < -1 ~   '< -1',
                                        mr.b < -0.5 ~  "-0.50:-1",
                                        mr.b < -0.25 ~ '-0.25:-0.50',
                                        mr.b < -0.10 ~ '-0.10:-0.25',
                                        mr.b < -0.01 ~ '-0.001:-0.10',
                                        mr.b < -0.001 ~'> -0.001',
                                        TRUE ~ '0')) %>% 
    mutate(mr.b_col_verbose = factor(mr.b_col_verbose, levels=c("> 1",
                                                                "0.50:1",
                                                                '0.25:0.50',
                                                                '0.10:0.25',
                                                                '0.001:0.10',
                                                                '<  0.001',
                                                                '0',
                                                                '> -0.001',
                                                                '-0.001:-0.10',
                                                                '-0.10:-0.25',
                                                                '-0.25:-0.50',
                                                                "-0.50:-1",
                                                                '< -1')))
  return(input)
}
