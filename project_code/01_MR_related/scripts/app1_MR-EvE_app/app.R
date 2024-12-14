library(shiny)
library(shinythemes)
library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(ggplot2)
library(plotly)
library(RColorBrewer)
library(kableExtra)
source("functions.R")


dat <- 
  #read_tsv("data_copy/bc_all_mr_madewR.tsv") %>% 
  read_tsv("data_copy/bc_all_mr_fromCIs.tsv") %>% 
  # subset 
  filter(exposure.sex != 'Males') %>% 
  filter(mr.method != 'Steiger null') %>% 
  # convert MR results to OR
  tidy_display_numbers()%>% 
  # deal ith al outcome related changes
  process_bc_outcomes() %>% 
  # create categories of exposure traits
  filter(!grepl("_raw",exposure.id)) %>% 
  create_exposure_categories() %>%
  # split several hundreds of proteins and metabolites in smaller chunks -- dont need these any more, as I fixed dynamic plot size
  #split_protein_exposures() %>% 
  #split_metabolite_exposures() %>% 
  add_exposure_labels() %>% 
  mutate(exposure_cat = ifelse(exposure_cat == 'Antrophometric', 'Anthropometric', exposure_cat)) %>% 
  mutate(exposure_cat = ifelse(exposure_cat == 'Drugs', 'Medication', exposure_cat))


chip_list <- c("Meta", "OncArray",  "iCOG2017",'GWASv1','GWASv2', 'Survival', 'UKBB')


### APP
ui <- fluidPage(align="center", theme = shinytheme("flatly"),
                
                
                titlePanel("MR-EvE results for breast cancer outcomes"),
                
                htmlOutput("selected_cat"),
                
                br(),
                
                # Output: Tabsets
                tabsetPanel(type = "tabs",
                            tabPanel("Interactive plot", plotlyOutput("bubbleplot2", height = "auto", width = "1000px")),
                            tabPanel("Static plot", plotOutput("bubbleplot1", height = "auto", width = "1000px")),
                            tabPanel("Info", 
                                     
                                     p(style="text-align: left; padding-left: 6cm; padding-right: 6cm; padding-top: 1cm; padding-bottom: 2cm",
                                       strong("About this app / how to read the plots:"),br(),
                                       "This Shiny app presents a visualisation of MR-EvE (Mendelian Randomization “Everything-vs-Everything”) 
                                       results for breast cancer outcomes across 12 exposure trait categories: ",
                                       br(),br(),
                                       "  –  anthropometric traits",br(),
                                       "  –  reproductive traits",br(),
                                       "  –  physical activity",br(),
                                       "  –  sleep traits",br(),
                                       "  –  diet and supplements",br(),
                                       "  –  alcohol",br(),
                                       "  –  smoking",br(),
                                       "  –  medication",br(),
                                       "  –  metabolites",br(),
                                       "  –  lipids",br(),
                                       "  –  proteins",br(),
                                       "  –  other biomarkers",
                                       
                                       br(), br(),
                                       "The plot is split into three outcome types (full sample, ER-, ER+), with specific outcome versions
                                       listed along the bottom. The exposure traits for a selected trait category are listed on the left, 
                                       including the author/consortium and sample size. ",
                                       br(),br(), 
                                       "The circles represent MR ‘best estimate’ effect direction (green – reduces risk, pink - increases 
                                       risk), effect size (the darker the colour, the greater the effect), -log10 p-value (the larger the 
                                       circle, the smaller the p-value). ",
                                       br(), br(),
                                       "The interactive version of the plot (recommended view) allows the user to hover and see the details
                                       of each data point, including the ‘best estimate’ MR method. The widgets below the plot allow the 
                                       user to filter the displayed data, for example, filter by p-value, Mixture of Experts (MoE) score,
                                       exclude all cases where CIs overlap the null, or enter specific trait names to display just them in the plot. ",
                                       br(), br(),
                                       strong("About the source data:"),br(),
                                       
                                       "The MR-EvE estimates for all pairs of exposures and outcomes were extracted from EpiGraphDB 
                                       (epigraphdb.org) (Liu et al 2021). The MR estimates were generated by the MR-MoE method (Hemani et al 2017). ",
                                       br(), br(),
                                       "MR-MoE is a machine learning framework that automates the selection of SNPs and an MR method for 
                                       use in any specific causal analysis. By testing over 20 strategies combining various MR methods and
                                       approaches for instrument filtering and selection, the framework selects the strategy that is most
                                       likely to be correct for a specific MR analysis by predicting the model of pleiotropy. This ‘best 
                                       estimate’ for each pair of traits is then stored as a relationship between GWAS traits in EpiGraphDB. ",
                                       br(),  br(), 
                                       "The table below shows the details about the included breast cancer outcomes. ",
                                       br(), br(),
                                       
                                       strong("Citations:"),br(),
                                       em('Liu et al., "EpiGraphDB: A database and data mining platform for health data science" (2021) Bioinformatics, 37(9), 1304–1311'),br(),
                                       em('Hemani et al., "Automating Mendelian randomization through machine learning to construct a putative causal map of the human phenome" (2017) bioRxiv'),br(),
                                       br(), br(),
                                       ">>> This app is a part of the work presented in ",
                                       em('"Integrating Mendelian randomization and literature-mined evidence for breast cancer risk factors"'), ", Vabistsevits et al 2022 " ,
                                       
                                        br(), br(),
                                       h4("Breast cancer GWAS summary details ", align='center'),
                                       tableOutput("outcome_table"),
                                       
                                      
                                    )
 
                            )
                
                ),
                hr(),
         
                
                fluidRow(
                  helpText("     Set display parameters:"),
                  br(),
                  
                  column(3, align="left",
                         selectInput(inputId ="category",
                                     label = "Trait category", 
                                     choices = list("Anthropometric" = 'Anthropometric',
                                                    "Reproductive" = 'Reproductive',
                                                    "Physical activity" = 'Physical activity',
                                                    "Diet and supplements" = 'Diet and supplements',
                                                    "Alcohol" = 'Alcohol',
                                                    "Smoking" = 'Smoking',
                                                    "Sleep" = 'Sleep',
                                                    "Lipids" = 'Lipids',
                                                    "Metabolites" = 'Metabolites',
                                                    #'Metabolites (pt. 1)'       = 'Metabolites (pt. 1)', 
                                                    #'Metabolites (pt. 2)'       = 'Metabolites (pt. 2)', 
                                                    #'Metabolites (met-) (pt. 3)'       = 'Metabolites (pt. 3)', 
                                                    'Proteins' = "Proteins",
                                                    #'Proteins (prot-) (pt. 1)' = "Proteins (pt. 1)",
                                                    #'Proteins (prot-) (pt. 2)' = "Proteins (pt. 2)",
                                                    #'Proteins (prot-) (pt. 3)' = "Proteins (pt. 3)",
                                                    #'Proteins (prot-) (pt. 4)' = "Proteins (pt. 4)",
                                                    #'Proteins (prot-) (pt. 5)' = "Proteins (pt. 5)",
                                                    #'Proteins (prot-) (pt. 6)' = "Proteins (pt. 6)",
                                                    "Other biomarkers" = 'Other biomarkers',
                                                    "Medication" = 'Medication'), 
                                     
                                     selected = 'Anthropometric'),
                         
                         br(),
                         checkboxGroupInput("ukb_versions", 
                                            p("Display UK Biobank (duplicated) traits from:"), 
                                            choices = c("MRC-IEU", "Neale lab"),
                                            selected = c("MRC-IEU", "Neale lab")),
                         
                         br(),
                         textInput(inputId = 'exposure_contains',
                                   label = "Exposure name contains:", 
                                   value = "")
                  ),
                  column(3,
                         sliderInput(inputId = "moe",
                                     label = "Minimum MoE score:",
                                     min = 0.5,
                                     max = 1,
                                     value = 0.5,
                                     step = 0.05),
                         
                         br(),
                         strong("CIs overlap the null"), 
                         checkboxInput("null_overlap", 
                                       label = "Exclude", 
                                       value = FALSE)
                  ),
                  
                  column(3,
                         sliderInput(inputId = "pval",
                                     label = "Include p-values (-log10 scale):",
                                     min = 0,
                                     max = 15,
                                     value = 0),
                         p("-log10(1e-8) = 8"),
                         
                         br(),
                         textInput(inputId = 'min_beta',
                                   label = "Smallest beta included \n (absolute value)", 
                                   value = 0)
                  ),
                  
                  
                  column(3,align="left",
                         checkboxGroupInput("outcomes", 
                                            p("Include outcomes:"), 
                                            #choices = chip_list,
                                            choices = list("BCAC 2017 meta-analysis" = 'Meta',
                                                           "OncoArray" = "OncArray", 
                                                           "iCOG2017" = "iCOG2017",
                                                           'GWAS v1' =  'GWASv1',
                                                           'GWAS v2' =  'GWASv2',
                                                           'Survival' = 'Survival',
                                                           "UK Biobank" = 'UKBB'),
                                            selected = chip_list[grepl('Meta',chip_list)]),
                         downloadButton('downloadData', 'Download currently displayed data')
                  )
                  
                  
                  
                  
                ),
                hr(),
                
                
                br(), br(),
                helpText("MR-EvE (Mendelian Randomization Everything-vs-Everything) results were extarcted from EpiGraphDB (epigraphdb.org), and were generated using the MR mixture-of-experts model (Hemani et al 2017)"), 
                br(),
                img(src='MRC_IEU_Bristol.png', align='centre', height = '50px')  , 
                br(), br(),
                span(uiOutput("twitter_link"), style="color:grey;font-size:12px;"), br()
                
                
                
)

# Define server logic required to draw a histogram ----
server <- function(input, output) {
  
  dataInput <- reactive({
    dat_sub <- dat %>% 
          filter(exposure_cat == input$category) %>% 
          filter(mr.moescore >= input$moe) %>% 
          filter(log10pval_trunc >= input$pval) %>% 
          filter(chip %in% input$outcomes) %>% 
          filter(mr.b >= as.numeric(input$min_beta) | mr.b <= as.numeric(input$min_beta)*-1 ) %>% 
          filter(grepl(input$exposure_contains, exposure, ignore.case=T)) 
    
    ukb_diff_sources<-dat_sub %>%
      select(exposure.trait, exposure.id, author, consortium, exposure.sample_size, year) %>% 
      filter(author %in% c("Neale", "Ben Elsworth")) %>% distinct() %>% 
      count(exposure.trait) %>% filter(n==2) 
    
    dat_sub <- dat_sub %>% 
      filter(!(exposure.trait %in% ukb_diff_sources$exposure.trait & !ukb_tag %in% input$ukb_versions))
    
    
    ## ad hoc filtering for specific categories
    if (input$category == 'Anthropometric'){
      antro_blacklist <- c('ieu-a-81','ieu-a-74', "ieu-a-73" ,"ieu-a-79" ,"ieu-a-72" ,"ieu-a-78",
                           'ieu-a-63',  'ieu-a-66', 'ieu-a-60',  'ieu-a-69' , 'ieu-a-61',
                           'ieu-a-54', 'ieu-a-55',  'ieu-a-49' , 'ieu-a-48', 'ieu-a-57' , 'ieu-a-50',
                           'ukb-b-12039', 'ukb-b-2303',
                           'ieu-a-2', 'ieu-a-835')
      dat_sub <- dat_sub %>% 
        filter(!grepl("arm|leg|first child", exposure.trait, ignore.case = T)) %>% 
        filter(!exposure.id %in% antro_blacklist)
  
    } else if (input$category %in% c('Physical activity')){
      dat_sub <- dat_sub %>% 
        filter(!grepl("leisure|mental", exposure, ignore.case = T))
      

    } else if (input$category %in% c('Diet and supplements')){  
      dat_sub<- dat_sub %>% filter(exposure_cat %in% c('Diet and supplements')) %>% 
            filter(!grepl("questionnaire", exposure))
      
      
    } else if (input$category %in% c('Other biomarkers')){
    dat_sub <- dat_sub %>% 
      filter(!grepl("FEV",exposure.trait)) 
    
    } else if (grepl('Proteins', input$category)){
      dat_sub <- dat_sub %>% 
        filter(mr.b > 0.01 | mr.b < -0.01)   
    }
    
    # exclude results overlapping null
    if (input$null_overlap){
      dat_sub<-dat_sub %>% filter(effect_direction != 'overlaps null')
    }
    
    # calculate name length in the selected subgroup
    dat_sub <- dat_sub %>% mutate(name_nchar = stringr::str_count(exposure)) # might use this for dynamic width adjustment
      
    return(dat_sub)
  })
  
  
  output$selected_cat <- renderText({ 
    HTML(paste("Showing the results for exposure trait category:", "<b>", input$category, "</b> (use the widget below to change) \n\n" ))
  })

  # static
  output$bubbleplot1 <- renderPlot(
    height = function() max((3 * 10 *length(unique(dataInput()$exposure.id)) + 80), 200),
    {
    plot_bubble_plot(dataInput(), font_size = 12)  
    }
  )
  
  # interactive
  output$bubbleplot2 <- renderPlotly({
    
    plot <- plot_bubble_plot(dataInput(), font_size = 8)  
    plotly::ggplotly(plot , tooltip = c("text"),
                     height = max((3 * 10 *length(unique(dataInput()$exposure.id)) + 100),220)
                     )
    
  })
  
  output$outcome_table <- function(){
    dat %>% create_outcomes_table() %>% 
      knitr::kable("html",
                   format.args = list(big.mark = ",", scientific = FALSE)) %>%
      kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = F) 
      #footnote(general = "BCAC: Breast Cancer Association Consortium") # does not want to display for some reason (worked before!)
  }
  
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste('breast_cancer_MR-EvE_results_for_', tolower(gsub(" ","-",input$category)), "_", Sys.Date(), '.csv', sep='')
    },
    content = function(file) {
      write.csv(dataInput()  %>% select(-pval_truncated, -empty_col, -exposure, -exposure.ss,	-tmp,	-exposure.ss_label,	-ukb_tag,	-name_nchar), 
                file, row.names = FALSE)
    }
  )
  

 
  url <- a("@marina_vab", href="https://twitter.com/marina_vab/")
  output$twitter_link <- renderUI({
    tagList("\nAny problems with the app? Let me know ", url)})
  

}

shinyApp(ui = ui, server = server)
