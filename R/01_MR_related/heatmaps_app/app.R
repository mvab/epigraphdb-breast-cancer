library(shiny)
library(shinythemes)
library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(reshape2)
library(ggplot2)
library(plotly)
library(RColorBrewer)
library(cowplot)
library(ggforce)

source("heatmap_functions.R")
source("functions_copy_from_mreveapp.R")
####
#current issue with the app:


###

inputs<- readRDS("data/inputs.rds") # created in heatmap_static.R

merged <- inputs$merged
protein_path_data <- inputs$protein_path_data
antro_blacklist <- inputs$antro_blacklist
passed_pairs <- inputs$passed_pairs
or_ci_data <- inputs$or_ci_data

data_full <- prepare_data(merged, protein_path_data, antro_blacklist,or_ci_data, passed_pairs)

outcome_list <- c("BCAC'17" ,"BCAC'20" ,"ER+" ,"ER-" , "Lum A",  "Lum B1",  "Lum B2" ,  "HER2", "TNBC" )  

### APP
ui <- fluidPage(align="center", theme = shinytheme("flatly"),
                
                
                titlePanel("MR effect direction heatmaps"),
                
                htmlOutput("selected_cat"),
                
                br(),
                
                # Output: Tabsets
                tabsetPanel(type = "tabs",
                            tabPanel("Static plot", plotOutput("heatmap1", height = "auto", width="auto")),
                            tabPanel("Interactive plot", plotlyOutput("heatmap2", height = "auto", width="auto")),
                            #tabPanel("Info", textOutput("info_text")) 
                            tabPanel("Info",

                              p(style="text-align: left; padding-left: 6cm; padding-right: 6cm; padding-top: 1cm; padding-bottom: 2cm",
                                strong("About this app / how to read the plot:"),br(),
                                "Heatmaps of Mendelian randomization (MR) effect direction for lifestyle and anthropometric traits,  metabolites, lipids, proteins 
                                as exposures and breast cancer (BCAC 2017 and 2020, including subtypes) as outcomes.",
                                br(), br(),
                                "The effect direction in MR is represents as colours: pink – positive (causal) effect, green – negative (protective) effect, white – no 
                                evidence of effect, based on 95% confidence intervals. The asterisk indicates results that passed FDR correction.",
                                br(),br(), 
                                strong("About source data:"),br(),
                                
                                "The displayed MR results were estimates using inverse-variance weighted (IVW) method (or Wald ratio if only 1 SNP was available).
                                By hovering over any data point (square) in the interactive version of the plot, you can see each MR analysis details: exposure ID 
                                in OpenGWAS, exposure details (sample size/sex, author/cohort name), OR and CIs, p-value and FDR-adjusted p-value, and number of 
                                SNPs used in each MR analysis.)",
                                br(), br(),
                                "Widget options: select category; split by sub-category (for lifestyle traits and proteains only); don't show exposure that have 
                                no evicdence of effect in any of the outcomes; for proteins only - show full names or abbreveations/gene names; select outcomes to display.",
                                br(), 
                                "All exposure traits are from mixed-sex samples unless otherwise specified (F: female-only) or are female-specific reproductive traits."
                                )
                              )
            
                ),
                
                
                hr(),
                
                
                fluidRow(
                  helpText("     Set display parameters:"),
                  br(),
                  
                  column(3, align="left",
                         selectInput(inputId ="category",
                                     label = "Select trait category", 
                                     choices = list("Antrophometric traits" = 'Antrophometric traits',
                                                    "Lifestyle traits" = 'Lifestyle traits',
                                                    "Metabolites" = 'Metabolites',
                                                    "Lipids" = 'Lipids',
                                                    "Proteins" = 'Proteins'), 
                                     
                                     selected = 'Antrophometric traits'),
                         br(),
                         strong("Sub-categories:"), 
                         checkboxInput("show_subcats", 
                                       label = "Show", 
                                       value = FALSE),
                         br(),
                         strong("Exposures with no effect"), 
                         checkboxInput("show_no_effect_exp", 
                                       label = "Don't show", 
                                       value = TRUE)
                         
                        
                  ),
               
                  column(3, align="left",
                         selectInput(inputId ="exposure_name_type",
                                     label = "Display names (for proteins only)", 
                                     choices = list("Full names" = 'full',
                                                    "Acronyms/genes" = 'gene',
                                                    "Mixed names(capped by length)" = 'mix'), 
                                     selected = 'Full names'),
              
         
                  ),
                  column(3,align="left",
                         checkboxGroupInput("outcomes", 
                                            p("Include outcomes:"), 
                                            choices = outcome_list,
                                            selected = outcome_list)
                  ),
                  

                ),
                hr(),
                
                br(), br(),
                
                hr(),
                br(), br(),
                helpText("help text here"), 
                
                img(src='MRC_IEU_Bristol.png', align='centre', height = '60px')  , 
                br(),
                span(uiOutput("twitter_link"), style="color:grey;font-size:12px;")
                
                
                
)

# Define server logic required to draw a histogram ----
server <- function(input, output) {
  
  dataInput <- reactive({
    dat_sub <- data_full %>% 
      filter(exposure_cat == input$category) %>% 
      filter(outcome %in% input$outcomes) 
    

    if (input$category == 'Proteins'){
      
      if (input$exposure_name_type == 'gene'){
         dat_sub <- dat_sub %>%  select(-exposure) %>% rename(exposure= gene)
      } else if (input$exposure_name_type == 'mix'){
        dat_sub <- dat_sub %>%  select(-exposure) %>% rename(exposure= name_mix)
      } else{
        # leave as as - full
      }
    }
    
    if (input$show_subcats){
      dat_sub <- dat_sub %>%  select(-exposure_cat) %>% rename(exposure_cat= exposure_cat_sub)
      print(unique(dat_sub$exposure_cat))
    }
    
    if(input$show_no_effect_exp){
      
      current_outcomes <- length(unique(dat_sub$outcome))
      
      exposures_to_drop <- dat_sub %>% 
        group_by(exposure.id, value) %>%
        count(exposure.id, value) %>% 
        filter(value ==0 & n ==current_outcomes ) %>% 
        pull(exposure.id)
      
      dat_sub <- dat_sub %>% filter(!exposure.id %in% exposures_to_drop)
      
    }
    
    
    
    return(dat_sub)
  })
  
  
  output$selected_cat <- renderText({ 
    HTML(paste("Showing the results for category:", "<b>", input$category, "</b> (use the widget below to change) \n\n" ))
  })
  
  output$heatmap1 <- renderPlot(
    width = function() 400 + (66 * length(unique(dataInput()$outcome))),
    height = function() 3 * 10 *length(unique(dataInput()$exposure.id)),
    
    #height = function() 3 * nrow(dataInput()),
    res = 96,
    
    {
    plot_heatmap(dataInput(), font_size = 11)  
  })
  
  output$heatmap2 <- renderPlotly(
    
    {
    plot <- plot_heatmap(dataInput(), font_size = 9)  
    plotly::ggplotly(plot , tooltip = c("text"), 
                     width = 400 + (66 * length(unique(dataInput()$outcome))),
                     height = 3 * 10 *length(unique(dataInput()$exposure.id))
                     )
    
  })
  

  
  
  url <- a("@marina_vab", href="https://twitter.com/marina_vab/")
  output$twitter_link <- renderUI({
    tagList("\nAny problems with the app? Let me know ", url)})
  
  
}

shinyApp(ui = ui, server = server)
