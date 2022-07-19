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
#current issues with the app:


###

inputs<- readRDS("data/inputs.rds") # created in heatmap_static.R

merged <- inputs$merged
protein_path_data <- inputs$protein_path_data
antro_blacklist <- inputs$antro_blacklist
passed_pairs <- inputs$passed_pairs
or_ci_data <- inputs$or_ci_data


names_tidy <- read_csv("data/renaming_key_tidy.csv") %>% select(exposure.id, exposure)# use new exposure column from here
merged<- merged %>%
  select(-exposure) %>%
  left_join(names_tidy, by =c("id.exposure" = "exposure.id")) %>% 
  select(exposure, everything())  %>% 
  filter(!is.na(exposure))
data_full <- prepare_data(merged, protein_path_data, antro_blacklist,or_ci_data, passed_pairs) 

# add column for sharing
data_full <- data_full %>% mutate(value_mtc = ifelse(!is.na(mtc), paste0(value, "*"), value)) %>% 
       mutate(value_mtc = factor(value_mtc, levels = c("-1*" ,"-1" , "0" ,  "1" ,  "1*" )))


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
                            tabPanel("Info",

                              p(style="text-align: left; padding-left: 6cm; padding-right: 6cm; padding-top: 1cm; padding-bottom: 2cm",
                                strong("About this app / how to read the plots:"),br(),
                                "This Shiny app presents the heatmap plots of effect direction of Mendelian randomization (MR) estimates 
                                from the analysis of various risk factor traits as exposures (lifestyle traits, anthropometric traits, 
                                metabolites, lipids, proteins) and breast cancer as the outcome (BCAC 2017 and 2020, including subtypes).",
                                br(), br(),
                                "The effect direction in MR is represented as colours: pink – positive (causal) effect, green – negative 
                                (protective) effect, white – no evidence of effect, based on 95% confidence intervals. The asterisk indicates
                                the results that passed the FDR correction.",
                                br(),br(), 
                                strong("About the source data:"),br(),
                                
                                "The MR results are estimated using the inverse-variance weighted (IVW) method (or Wald ratio if only 1 SNP 
                                was available). By hovering over any data point (square) in the interactive version of the plot, you can see
                                the details of each MR analysis: exposure ID in OpenGWAS (gwas.mrcieu.ac.uk), exposure details (sample size/sex, author/cohort name),
                                OR and CIs, p-value and FDR-adjusted p-value, and the number of SNPs).",
                                br(), br(),
                                em("Widget options: "),"select category; split by sub-category (for lifestyle traits and proteins only); don't show exposures 
                                that have no evidence of effect in any of the outcomes; for proteins only - show full names or abbreviations/gene names; 
                                select outcomes to display.",
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
                                     choices = list("Anthropometric traits" = 'Anthropometric traits',
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
                                       value = TRUE),
                  
                         br(),
                         strong("FDR correction"), 
                         checkboxInput("rows_with_fdr", 
                                       label = "Show passed", 
                                       value = F)
                         
   
                  ),
               
                  column(3, align="left",
                         selectInput(inputId ="exposure_name_type",
                                     label = "Display names (for proteins only)", 
                                     choices = list("Full names" = 'full',
                                                    "Abbreviations / gene names" = 'gene',
                                                    "Mixed names (capped by length)" = 'mix'), 
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
                br(),   br(),
                img(src='MRC_IEU_Bristol.png', align='centre', height = '50px')  , 
                br(), br(),
                span(uiOutput("twitter_link"), style="color:grey;font-size:12px;"), br()
                
                
                
)

# Define server logic required to draw a histogram ----
server <- function(input, output) {
  
  dataInput <- reactive({
    dat_sub <- data_full %>% 
      filter(exposure_cat == input$category) %>% 
      filter(outcome %in% input$outcomes) 
    

    if (input$category == 'Proteins'){
      
      if (input$exposure_name_type == 'gene'){
         dat_sub <- dat_sub %>% 
                  arrange( value, outcome, main_path) %>%
                  select(-exposure) %>% rename(exposure= gene)
         dat_sub <- dat_sub %>% 
           mutate(exposure = factor(exposure, levels = unique(dat_sub$exposure)))
         
      } else if (input$exposure_name_type == 'mix'){
        dat_sub <- dat_sub %>%  
                  arrange( value, outcome, main_path) %>%
                  select(-exposure) %>% rename(exposure = name_mix) 
        dat_sub <- dat_sub %>% 
                  mutate(exposure = factor(exposure, levels = unique(dat_sub$exposure)))
      } else{
        # leave as as - full
      }

    }
    
    if (input$show_subcats){
      
      if (input$category == 'Proteins'){
        
        dat_sub <- dat_sub %>% 
          select(-exposure_cat) %>% rename(exposure_cat=exposure_cat_sub) %>% 
          mutate(exposure_cat = factor(exposure_cat, levels = c("Immune System", 'Metabolism', 'Signal Transduction','Developmental Biology', "Other", "Not mapped")))
        
      } else{ # all other cats
        dat_sub <- dat_sub %>%  select(-exposure_cat) %>% rename(exposure_cat= exposure_cat_sub)
      }
      
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
    if(input$rows_with_fdr){
      
      exposures_to_keep <- dat_sub %>% 
        group_by(exposure.id, mtc) %>%
        count(exposure.id, mtc) %>% 
        filter(mtc =="\n*" & n >=1 ) %>% 
        pull(exposure.id)
      
      dat_sub <- dat_sub %>% filter(exposure.id %in% exposures_to_keep)
      
    }
    
    # calculate name length in the selected subgroup
    dat_sub <- dat_sub %>% mutate(name_nchar = stringr::str_count(exposure))
    

    return(dat_sub)
  })
  
  
  output$selected_cat <- renderText({ 
    HTML(paste("Showing the results for category:", "<b>", input$category, "</b> (use the widget below to change) \n\n" ))
  })
  
  output$heatmap1 <- renderPlot(
    width = function() (3 * max(dataInput()$name_nchar)) + max((66 * length(unique(dataInput()$outcome))), 265 ), # set size for trait name based on max length + plot width based on number of selected outcomes or min 265
    height = function() max((3 * 10 * length(unique(dataInput()$exposure.id))), 200), # calculate or use 200
    
    res = 96,
    
    {
    plot_heatmap2(dataInput(), font_size = 11)  
  })
  
  output$heatmap2 <- renderPlotly(
    
    {
    plot <- plot_heatmap2(dataInput(), font_size = 9, star_size = 4)  
    plotly::ggplotly(plot , tooltip = c("text"), 
                     width = (3 * max(dataInput()$name_nchar))  + max((66 * length(unique(dataInput()$outcome))), 265 ), # set size for trait name based on max length + plot width based on number of selected outcomes or min 265
                     height = max((3 * 10 *length(unique(dataInput()$exposure.id))),200) # calculate or use 200
                     )
    
  })
  

  
  
  url <- a("@marina_vab", href="https://twitter.com/marina_vab/")
  output$twitter_link <- renderUI({
    tagList("\nAny problems with the app? Let me know ", url)})
  
  
}

shinyApp(ui = ui, server = server)
