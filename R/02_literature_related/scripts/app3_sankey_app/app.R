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
                
                
                titlePanel("Literature spaces overlap"),
                br(),
                htmlOutput("selected_trait"),
                em("See Info tab for details"),
                br(), br(),
                
                fluidRow(
                  column(3, align="left",
                         selectInput(inputId ="trait",
                                     label = "Select trait", 
                                     choices = names(trait_tidy_subset), 
                                     
                                     selected = 'Cardiotrophin-1')),
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

                                     p(style="text-align: left; padding-left: 6cm; padding-right: 6cm; padding-top: 1cm; padding-bottom: 2cm",
                                       
                                       "Sankey plot is a visualisation of literature spaces overlap between a case study trait and breast cancer. Individual
                                       triples are linked via overlapping terms within each space and between two spaces.  ",
                                       br(),br(),
                                       strong("  –  literature space"),  " is a collection of triples related to a specific trait, e.g. breast cancer or a trait,
                                     which is also a GWAS in EpiGraphDB.", 
                                       br(),br(),
                                       strong("  –  literature triple"), " triple is a relationship between two literature terms connected by a verb/predicate that 
                                       describes a relationship between them, e.g. 'subject term 1 – predicate – object term 2'. A subject/object term could be any
                                       biological entity (gene name, drug, phenotype, disease) and a predicate is a verb that represents a relationship between 
                                       the two terms (affects, causes, inhibits, reduces, associated with, etc.), e.g. 'protein A – stimulates – protein B'.",
                                       
                                       "The underlying data comes from SemMedDB, a well-established repository of literature-mined semantic triples mined from
                                       titles and abstracts in PubMed.",
                                       br(),br(),
                                       
                                       strong("  –  triple linkage"),  " is a method of linking triples via overlapping terms into chains (see Figure (a) below). 
                                       This can be done within a literature space or between literature spaces. If a trait is available as a term within its space 
                                       (e.g. IGF1 term in IGF1 literature space), it can be used as an anchor to which sequentially linked triples are attached. 
                                       In the breast cancer space, we directly link triples to the ‘breast cancer’ term or via up to three triples. If no anchoring 
                                       term is available, we can still connect the terms from the breast cancer space to the terms from a trait space, but there 
                                       will not be a clear path from a trait to breast cancer.",
                                       br(),br(),
                                       strong("  –  literature spaces overlap"),  " is a method that was developed for joining two literature spaces together
                                       (see Figures (b) and (c)). It involves triple linkage within one literature space (trait and breast cancer separately)
                                       and then connecting them either via outer triples overlap or any combination of inner triples (see arrows in Figure (b)).  ",
                                       br(),br(),
                                       strong("Why?"), " Overlapping literature spaces allowed us to look for intermediates between a risk factor trait and breast cancer.",
                                       br(),br(),
                                       "See further details on the method in the publication Methods.",
                                       br(),br(),br(),
                                       em("How to read the plot:"),
                                       br(),br(),
                                       "In the plot, trait triples are shown in blue and breast cancer triples in orange. The thickness of the line indicates the 
                                       frequency of the triple pair (i.e. the thicker the line, the more publications mention the specific relationship).
                                       The full plot includes all links; the subset plot shows only those links that have been mentioned in minimum N 
                                       publications (specify in the widget). ",
                                       br(),br(),br(),
                                         
                                         
                                         
                                
                                         
                                       img(src='litmethod_abc_v1.png', align='centre', height = '500px'),br(),
                                       strong("Schemas of literature spaces overlap methods."), "(a) A basic example of two triples linkage via an overlapping term. 
                                       (b) Molecular traits literature space overlap method. Molecular trait literature space (blue box) and breast cancer 
                                       literature space (purple box) both contain triples linked via overlapping terms. The outer terms are restricted 
                                       (‘i.e. ‘anchored’) to the term representing the trait (T: trait) and the term representing breast cancer outcome 
                                       (BC: breast cancer). All intermediate terms are arbitrarily named A-F. Terms T, A, B can be overlapping terms with 
                                       C, D, E, and F, which links the trait and breast cancer literature spaces. The arrows represent all possible paths 
                                       from T to BC alternative to the full path going through all intermediates. (c) Lifestyle traits literature space 
                                       overlap method. Lifestyle trait literature space (green circle) contains triples that cannot be anchored to any specific
                                       term representing the trait. The unlinked triples in the lifestyle trait space are matched to overlap (via B) with any 
                                       triples in the breast cancer literature space with a path to breast cancer. The linked triples are then connected with 
                                       any preceding triples in lifestyle trait space (via A), adding X-A triples into the spaces overlap.")
                              )
                            
                ),
                
                
                hr(),
                
                

                hr(),
                
                br(), br(),
                
                hr(),
                br(), br(),

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
  
  
  output$selected_trait <- renderText({ 
    HTML(paste("Showing Sankey plot for", "<b>", input$trait, "</b> and breast cancer \n\n" ))
  })
  
  
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
