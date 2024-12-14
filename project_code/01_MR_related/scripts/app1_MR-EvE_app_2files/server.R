library(shiny)
library(shinythemes)
library(readr)
library(tidyr)
library(stringr)
library(ggplot2)
library(plotly)
library(RColorBrewer)
library(kableExtra)
library(dplyr)

source("functions.R")


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
    # split several hunderds of proteins in smaller chunks
    split_protein_exposures() %>% 
    add_exposure_labels()

chip_list <- c("Meta", "OncArray",  "iCOG2017",'iCOG2015','GWASold1','GWASold2', 'Survival', 'UKBB')



shinyServer(function(input, output) {
    
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
        paste("Showing the results for exposure trait category:", "<b>", input$category, "</b>")
    })
    
    output$bubbleplot1 <- renderPlot({
        
        plot_bubble_plot(dataInput(), font_size = 11)  
        
    })
    
    output$bubbleplot2 <- renderPlotly({
        
        plot <- plot_bubble_plot(dataInput(), font_size = 8)  
        plotly::ggplotly(plot , tooltip = c("text"))
        
    })
    
    output$outcome_table <- renderTable(
        dat %>% create_outcomes_table() %>% 
            mutate(`Sample size` = format(`Sample size`, big.mark = ",", scientific = FALSE),
                   nSNPs = format(nSNPs, big.mark = ",", scientific = FALSE),
                   `# of cases` = format(`# of cases`, big.mark = ",", scientific = FALSE)),
        striped=T, hover=T)
    
    url <- a("@marina_vab", href="https://twitter.com/marina_vab/")
    output$twitter_link <- renderUI({
        tagList("\nAny problems with the app? Let me know ", url)})
    
    #output$outcome_table <- function(){
     #   dat %>% create_outcomes_table() %>% 
            
            #knitr::kable("html",
            #             format.args = list(big.mark = ",", scientific = FALSE)) %>%
            #kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = F) %>% 
            #footnote(general = "BCAC: Breast Cancer Association Consortium")
    #}
    
    
    
})