library(shiny)
library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(ggplot2)
library(plotly)
library(RColorBrewer)
source("../functions.R")

####
#current issue with the app:
#- colour scale problem when applying filters
#- survival data not displaying right?
#- add beta widget
#- add 'ovelaps null' widget
#- set plot size dynamically
#- if nrow> x apply filters at display
#
###

dat <- read_tsv("../../query_results/bc_all_mr.tsv") %>% 
  # subset 
  filter(exposure.sex != 'Males') %>% 
  filter(outcome.id != 'ebi-a-GCST007236') %>% 
  # convert MR results to OR
  tidy_display_numbers()%>% 
  # deal ith al outcome related changes
  process_bc_outcomes() %>% 
  # create categories of exposure traits
  create_exposure_categories() %>% 
  add_exposure_labels()


chip_list <- c("Meta", "OncArray",  "iCOG2017",'iCOG2015','GWASold1','GWASold2', 'Survival')


### APP

ui <- fluidPage(
 
  
  titlePanel("MR-EvE results for breast cancer outcomes"),
  
  textOutput("selected_cat"),
  
  br(),
  
  # Output: Tabsets
  tabsetPanel(type = "tabs",
              tabPanel("Static plot", plotOutput("bubbleplot1", height = "700px")),
              tabPanel("Interactive plot", plotlyOutput("bubbleplot2", height = "700px"))),
    
  
  hr(),
    
  fluidRow(
      helpText("     Set display parameters:"),
      br(),
      
    column(3,
      selectInput(inputId ="category",
                  label = "Trait category", 
                  choices = list("Antrophometric" = 'antrophometric',
                                 "Reproductive" = 'reproductive',
                                 "Activity" = 'activity',
                                 "Supplements" = 'vitamin',
                                 "Diet" = 'diet',
                                 "Alcohol" = 'alcohol',
                                 "Smoking" = 'smoking',
                                 'Metabolites' = 'metabolite_measures',
                                 "Proteins" = 'protein_measures',
                                 "Other biomarkers" = 'other_biomarkers'), 
                  
                  selected = 'antrophometric')
      ),
    column(3,
       sliderInput(inputId = "moe",
                   label = "Minimum MoE score:",
                   min = 0.5,
                   max = 1,
                   value = 0.5,
                   step = 0.05)
    ),
       
    column(3,
      sliderInput(inputId = "pval",
                  label = "Include p-values (-log10 scale):",
                  min = 0,
                  max = 15,
                  value = 0),
      p("-log10(1e-8) = 8")
    ),
    
    #column(3,
    #       sliderInput(inputId = "beta",
    #                   label = "Include beta (absolute):",
    #                   min = 0,
    #                   max = 15,
    #                   value = 0)
    #),
   
    column(3,
            checkboxGroupInput("outcomes", 
                               p("Include outcomes:"), 
                               choices = chip_list,
                               selected = chip_list[!grepl('iCOG2015|GWASold2',chip_list)])
    )
    
    
    
    )
    
    
  )



# Define server logic required to draw a histogram ----
server <- function(input, output) {
  
  dataInput <- reactive({
    dat_sub <- dat %>% 
          filter(exposure_cat == input$category) %>% 
          filter(mr.moescore >= input$moe) %>% 
          filter(log10pval_trunc >= input$pval) %>% 
          filter(chip %in% input$outcomes)
    
    ## ad hoc filtering for specific categories
    if (input$category == 'antrophometric'){
      dat_sub <- dat_sub %>% 
        filter(!grepl("arm|leg|first child", exposure.trait, ignore.case = T)) 
    } else if (input$category == 'alcohol'){
      dat_sub <- dat_sub %>% 
         filter(!grepl('dehydrogenase', exposure.trait))
    } else if (input$category == 'protein_measures'){
    dat_sub <- dat_sub %>% 
      filter(!grepl("_raw",exposure.id)) %>% # keep only invr
      filter(mr.b > 0.01 | mr.b < -0.01)   
    }
      
    return(dat_sub)
  })
  
  
  output$selected_cat <- renderText({ 
    paste("Show the results for exposure category:", input$category)
  })

  output$bubbleplot1 <- renderPlot({
    
    plot_bubble_plot(dataInput(), font_size = 11)  
    
  })
  
  output$bubbleplot2 <- renderPlotly({
    
    plot <- plot_bubble_plot(dataInput(), font_size = 8)  
    plotly::ggplotly(plot , tooltip = c("text", "y", "x"))
    
  })
  
}

shinyApp(ui = ui, server = server)