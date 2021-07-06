library(shiny)
library(shinythemes)
library(readr)
library(tidyr)
library(stringr)
library(ggplot2)
library(plotly)
library(RColorBrewer)
#library(kableExtra)
library(dplyr)
source("functions.R")


chip_list <- c("Meta", "OncArray",  "iCOG2017",'iCOG2015','GWASold1','GWASold2', 'Survival')


shinyUI(
    
    fluidPage(align="center", theme = shinytheme("flatly"),
        
        span(strong(uiOutput("tab")), style="color:red"),
              
        titlePanel("MR-EvE results for breast cancer outcomes"),
        
        htmlOutput("selected_cat"),
        
        br(),
        
        # Output: Tabsets
        tabsetPanel(type = "tabs",
                    tabPanel("Static plot", plotOutput("bubbleplot1", height = "700px", width = "900px")),
                    tabPanel("Interactive plot", plotlyOutput("bubbleplot2", height = "700px", width = "900px"))),
        
        
        hr(),
        
        fluidRow(
            helpText("     Set display parameters:"),
            br(),
            
            column(3,align="left",
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
    helpText("BCAC: Breast Cancer Association Consortium"),
    tabPanel("Outcomes table", 
             tableOutput("outcome_table"),
             
     hr(),
     br(), br(),
     helpText("MR-EvE (Mendelian Randomization Everything-vs-Everything) results were extarcted from EpiGraphDB (epigraphdb.org), and were generated using the MR mixture-of-experts model (Hemani et al 2017)"), 
     
    img(src='MRC_IEU_Bristol.png', align='centre', height = '60px')         
    )
    
)
)
