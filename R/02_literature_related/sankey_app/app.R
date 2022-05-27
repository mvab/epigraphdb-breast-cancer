library(shiny)
library(shinythemes)
library(readr)
library(tidyr)
library(stringr)
library(vroom)
library(dplyr)
library(networkD3)



source("functions_literature.R")


# Load data ----

### get breast cancer data ----
dat<- vroom("data/breast_cancer_litspace_prod.csv",show_col_types = FALSE) %>%  mutate(lit.id = as.character(lit.id))

# tidy space 
bc_triples_tidy_count <- tidy_lit_space(dat) %>%    
  mutate(term1.type_verbose  = ifelse(grepl("obesity", term1, ignore.case = T), 'antro', term1.type_verbose)) %>%
  mutate(term2.type_verbose  = ifelse(grepl("obesity", term2, ignore.case = T), 'antro', term2.type_verbose))

bc_triples <- get_breast_cancer_triples(bc_triples_tidy_count)

### get traits data ----
available_traits <- read_csv("data/traits_meta.csv") 
single_ids <- available_traits %>% filter(type == 'single') %>% pull(id)
single_names <- available_traits %>% filter(type == 'single') %>% pull(trait_name)
combined_ids <- available_traits %>% filter(type == 'combined') %>% pull(id)
combined_names <- available_traits %>% filter(type == 'combined') %>% pull(trait_name)


load("data/lit_spaces_finalset_tidy.RData")
trait_tidy_subset1 <- tidy_litspace[single_ids] 
names(trait_tidy_subset1) <- single_names
rm(tidy_litspace)

load("data/lit_spaces_combined_traits_tidy.RData")
trait_tidy_subset2 <- tidy_combined_litspace[combined_ids]
names(trait_tidy_subset2) <-  combined_names
rm(tidy_combined_litspace)

trait_tidy_subset <- append(trait_tidy_subset1, trait_tidy_subset2)


# APP ----

## UI ----

ui <- fluidPage(align="center", theme = shinytheme("flatly"),
                
                
                titlePanel("Sankey plots test"),
                
                htmlOutput("selected_trait"),
                
                br(),
                
                # Output: Tabsets
                tabsetPanel(type = "tabs",
                            tabPanel("Sankey plot", sankeyNetworkOutput("sankey_plot", height = "1200px", width="auto")),
                           
                            #tabPanel("Info", textOutput("info_text")) 
                            tabPanel("Info",
                                     
                                     p(style="text-align: justify;",
                                       strong("How to read the plot:"),br(),
                                       "text to be added",
                                       br(),br(),
                                       strong("About source data:"),br(),
                                       "text to be added")
                            )
                            
                ),
                
                
                hr(),
                
                
                fluidRow(
                  helpText("     Set display parameters:"),
                  br(),
                  
                  column(3, align="left",
                         selectInput(inputId ="trait",
                                     label = "Select trait", 
                                     choices = names(trait_tidy_subset), 
                                     
                                     selected = 'IGF-1'),
                         br(),
                         strong("Display mode (full):"), 
                         checkboxInput("show_subset", 
                                       label = "subset", 
                                       value = FALSE)
              
                         
                         
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

## SERVER ----

server <- function(input, output) {
  
  dataInput <- reactive({
   
    
    
    key_term_expt <- available_traits %>% filter(trait_name == input$trait) %>% pull(key_term_extr) %>% str_split(",")
    trait_tidy <- trait_tidy_subset[[input$trait]]
    
    trait_triples <- extract_two_triples_for_trait(trait_tidy,  
                                                   KEY_TERM = key_term_expt[[1]])
    
    key_term_anchor <- available_traits %>% filter(trait_name == input$trait) %>% pull(key_term_anchor)
    
    
    out <- overlap_trait_and_bc(trait_triples$joined_triples, KEY_TERM =key_term_anchor, n_filter=2, sankey_font = 13, bc_triples)
    
    
    return(out$full_sankey_data)
  })
  
  
  output$selected_trait <- renderText({ 
    HTML(paste("Showing Sankey plot for trait:", "<b>", input$trait, "</b> (use the widget to change) \n\n" ))
  })
  
  
  
  
  
  output$sankey_plot <- renderSankeyNetwork(
   # width = function() 1000,
    #height = function() 1000,
   # res = 96,
    
    {
      make_sankey(dataInput(), fontSize=13, colour_links = T)

    })
  
  
  
 #output$heatmap1 <- renderPlot(
 #  width = function() 400 + (66 * length(unique(dataInput()$outcome))),
 #  height = function() 3 * 10 *length(unique(dataInput()$exposure.id)),
 #  
 #  #height = function() 3 * nrow(dataInput()),
 #  res = 96,
 #  
 #  {
 #    plot_heatmap(dataInput(), font_size = 11)  
 #  })
 #
 #output$heatmap2 <- renderPlotly(
 #  
 #  {
 #    plot <- plot_heatmap(dataInput(), font_size = 9)  
 #    plotly::ggplotly(plot , tooltip = c("text"), 
 #                     width = 400 + (66 * length(unique(dataInput()$outcome))),
 #                     height = 3 * 10 *length(unique(dataInput()$exposure.id))
 #    )
 #    
 #  })
  
  
  
  
  url <- a("@marina_vab", href="https://twitter.com/marina_vab/")
  output$twitter_link <- renderUI({
    tagList("\nAny problems with the app? Let me know ", url)})
  
  
}

shinyApp(ui = ui, server = server)
