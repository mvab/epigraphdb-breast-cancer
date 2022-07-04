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
                
                
                titlePanel("Sankey plot of literature triples overlap between a trait and breast cancer"),
                br(),
                #htmlOutput("selected_trait"),
                
                fluidRow(
                  column(3, align="left",
                         selectInput(inputId ="trait",
                                     label = "Select trait", 
                                     choices = names(trait_tidy_subset), 
                                     
                                     selected = 'IGF-1')),
                  column(3, align="left",
                         radioButtons("display_mode", 
                                      label = "Display mode",
                                      choices = list("Full" = 'full',
                                                     "Subset" = 'subset'), 
                                      selected = "full")),
                  column(3, align="left",
                         selectInput(inputId ="subset_cutoff",
                                     label = "Subset includes links with min size", 
                                     choices = list("2" = 2,
                                                    "3" = 3,
                                                    "4" = 4), 
                                     selected = 2)),
                   column(3, align="left",
                          downloadButton('downloadData', 'Download')),
                                
                         
                  
                ),
                br(),
                
                # Output: Tabsets
                tabsetPanel(type = "tabs",
                            tabPanel("Sankey plot", uiOutput('sankey_plot')),#sankeyNetworkOutput("sankey_plot", width="auto", height = "1000px")),
                           
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
   
    trait_group <- available_traits %>% filter(trait_name == input$trait) %>% pull(group)
    
    if (trait_group == 'molecular'){ 
    
      key_term_expt <- available_traits %>% filter(trait_name == input$trait) %>% pull(key_term_extr) %>% str_split(",")
      key_term_anchor <- available_traits %>% filter(trait_name == input$trait) %>% pull(key_term_anchor)
      
      
      trait_tidy <- trait_tidy_subset[[input$trait]]
      trait_triples <- extract_two_triples_for_trait(trait_tidy,   KEY_TERM = key_term_expt[[1]])
      
      out <- overlap_trait_and_bc(trait_triples$joined_triples, KEY_TERM =key_term_anchor, n_filter=input$subset_cutoff,  bc_triples)
      
      
      if (input$display_mode == 'full'){
          return(out$full_sankey_data)
      } else if (input$display_mode == 'subset'){
          return(out$subset_sankey_data)
      }
    
      
    } else if (trait_group == 'lifestyle'){
      
      trait_tidy <- trait_tidy_subset[[input$trait]]
      trait_triples <- extract_lifestyle_main_triples(trait_tidy)
      
      out <- overlap_lifestyle_trait_and_bc(trait_triples, bc_triples, n_filter = input$subset_cutoff)
      
      
      if (input$display_mode == 'full'){
          return(out$sankey_data)
      } else if (input$display_mode == 'subset'){
          return(out$sankey_data_filtered)
      }

      
    }

  })
  
  
  #output$selected_trait <- renderText({ 
  #  HTML(paste("Showing Sankey plot for trait:", "<b>", input$trait, "</b> \n\n" ))
  #})
  #
  
  output$sankey_plot_prep <- renderSankeyNetwork({
      make_sankey(dataInput(), fontSize=13, colour_links = T)
    })
  #https://stackoverflow.com/questions/58208904/reactive-height-for-sankeynetworkoutput-from-networkd3
  output$sankey_plot <- renderUI({
    h <- 300 * log10(nrow(dataInput()))
    sankeyNetworkOutput("sankey_plot_prep", height = h)
  })
  
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste('sankey-data_', tolower(gsub(" ","-",input$trait)), "_breast-cancer_literature-overlap_", Sys.Date(), '.csv', sep='')
    },
    content = function(file) {
      write.csv(dataInput(), file, row.names = FALSE)
    }
  )
  
  
  
  
  url <- a("@marina_vab", href="https://twitter.com/marina_vab/")
  output$twitter_link <- renderUI({
    tagList("\nAny problems with the app? Let me know ", url)})
  
  
}

shinyApp(ui = ui, server = server)
