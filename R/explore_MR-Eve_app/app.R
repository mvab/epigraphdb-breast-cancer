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

####
#current issue with the app:
#- colour scale problem when applying filters -- only become a problem if you filter down to very low number of exposures
#- set plot size dynamically
#- if nrow> x apply filters at display
# include any other trait categories

# reproductive traits??

###

dat <- read_tsv("data_copy/bc_all_mr.tsv") %>% 
  # subset 
  filter(exposure.sex != 'Males') %>% 
  filter(!outcome.id %in% c('ebi-a-GCST007236', 'ebi-a-GCST004988')) %>% 
  filter(mr.method != 'Steiger null') %>% 
  # convert MR results to OR
  tidy_display_numbers()%>% 
  # deal ith al outcome related changes
  process_bc_outcomes() %>% 
  # create categories of exposure traits
  create_exposure_categories() %>%
  # split several hunderds of protins in smaller chunks
  split_protein_exposures() %>% 
  add_exposure_labels()


chip_list <- c("Meta", "OncArray",  "iCOG2017",'iCOG2015','GWASold1','GWASold2', 'Survival')


### APP
ui <- fluidPage(align="center", theme = shinytheme("flatly"),
                
                
                titlePanel("MR-EvE results for breast cancer outcomes"),
                
                htmlOutput("selected_cat"),
                
                br(),
                
                # Output: Tabsets
                tabsetPanel(type = "tabs",
                            tabPanel("Static plot", plotOutput("bubbleplot1", height = "700px", width = "900px")),
                            tabPanel("Interactive plot", plotlyOutput("bubbleplot2", height = "700px", width = "900px"))
                            ),
                
                
                hr(),
                
                fluidRow(
                  helpText("     Set display parameters:"),
                  br(),
                  
                  column(3, align="left",
                         selectInput(inputId ="category",
                                     label = "Trait category", 
                                     choices = list("Antrophometric" = 'Antrophometric',
                                                    "Reproductive" = 'Reproductive',
                                                    "Physical activity" = 'Physical activity',
                                                    "Diet and supplements" = 'Diet and supplements',
                                                    "Alcohol" = 'Alcohol',
                                                    "Smoking" = 'Smoking',
                                                    'Metabolites (met-)'       = 'Metabolites',    
                                                    'Proteins (prot-) (pt. 1)' = "Proteins (pt. 1)",
                                                    'Proteins (prot-) (pt. 2)' = "Proteins (pt. 2)",
                                                    'Proteins (prot-) (pt. 3)' = "Proteins (pt. 3)",
                                                    'Proteins (prot-) (pt. 4)' = "Proteins (pt. 4)",
                                                    'Proteins (prot-) (pt. 5)' = "Proteins (pt. 5)",
                                                    'Proteins (prot-) (pt. 6)' = "Proteins (pt. 6)",
                                                    "Other biomarkers" = 'Other biomarkers',
                                                    "Drugs" = 'Drugs'), 
                                     
                                     selected = 'Antrophometric'),
                         
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
                                            choices = chip_list,
                                            selected = chip_list[!grepl('iCOG2015|GWASold2',chip_list)])
                  )
                  
                ),
                hr(),
                
                br(), br(),
                h4("Breast cancer GWAS summary details ", align='center'),
                tabPanel("Outcomes table", tableOutput("outcome_table")),
                
                hr(),
                br(), br(),
                helpText("MR-EvE (Mendelian Randomization Everything-vs-Everything) results were extarcted from EpiGraphDB (epigraphdb.org), and were generated using the MR mixture-of-experts model (Hemani et al 2017)"), 
                
                img(src='MRC_IEU_Bristol.png', align='centre', height = '60px') 
                
                
                
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
    
    ## ad hoc filtering for specific categories
    if (input$category == 'Antrophometric'){
      dat_sub <- dat_sub %>% 
        filter(!grepl("arm|leg|first child", exposure.trait, ignore.case = T)) 
      
    } else if (input$category %in% c('Other biomarkers')){
    dat_sub <- dat_sub %>% 
      filter(!grepl("_raw",exposure.id))  
    
    } else if (grepl('Proteins', input$category)){
      dat_sub <- dat_sub %>% 
        filter(mr.b > 0.01 | mr.b < -0.01)   
    }
    
    # exclude results overlapping null
    if (input$null_overlap){
      dat_sub<-dat_sub %>% filter(effect_direction != 'overlaps null')
    }

      
    return(dat_sub)
  })
  
  
  output$selected_cat <- renderText({ 
    HTML(paste("Showing the results for exposure trait category:", "<b>", input$category, "</b>"))
  })

  output$bubbleplot1 <- renderPlot({
    
    plot_bubble_plot(dataInput(), font_size = 11)  
    
  })
  
  output$bubbleplot2 <- renderPlotly({
    
    plot <- plot_bubble_plot(dataInput(), font_size = 8)  
    plotly::ggplotly(plot , tooltip = c("text"))
    
  })
  
  output$outcome_table <- function(){
    dat %>% create_outcomes_table() %>% 
    knitr::kable("html",
                 format.args = list(big.mark = ",", scientific = FALSE)) %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = F) %>% 
    footnote(general = "BCAC: Breast Cancer Association Consortium")
  }

 
}

shinyApp(ui = ui, server = server)
