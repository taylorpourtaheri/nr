library(shiny)
library(shinythemes)
library(DT)

devtools::load_all()

# load in the deg
de_string <- readRDS('de_string_v11.RDS')
myc_de <- de_string$MYC


# Define UI for application
ui <- fluidPage(
    
    theme = shinytheme("flatly"),
    
    # Application title
    titlePanel(h1("Differential Gene Expression Data Analysis", align = 'center')),
    
    # Sidebar 
    sidebarLayout(
        sidebarPanel(width = 3,
                     textInput(inputId = 'target', label = 'Enter the target gene'),
                     selectInput(inputId = 'method', label = 'Centrality method', choices = c('avg_strength',
                                                                                              'betweenness',
                                                                                              'degree',
                                                                                              'evcent_w',
                                                                                              'strength',
                                                                                              'evcent_uw')),
                     selectInput(inputId = 'logFC', label = 'Log fold change filter', choices = c(0.0, 0.5, 1.0, 1.5, 2.0)),
                     selectInput(inputId = 'pvalue', label = 'P-value filter', choices = c(0.05, 0.10, 0.15, 0.20, 0.25)),
                     selectInput(inputId = 'connected', label = 'Connected subgraph filter', choices = c("TRUE"=1,"FALSE"=0)),
                     selectInput(inputId = 'weighted', label = 'Weighted scoring', choices = c("TRUE"=1,"FALSE"=0)),
                     actionButton(inputId = "build", label = 'Build network')
        )
        ,
        
        mainPanel(
            width = 9,
            navbarPage(title = 'Outputs:',
                       tabPanel('Network Plot',plotOutput('nodePlot')),
                       tabPanel('Target Results',DT::dataTableOutput('targetPerformance')),
                       tabPanel('Table of Genes',DT::dataTableOutput('topGenes')),
                       selectInput("dataset", "Choose a dataset:",
                                   choices = c("Target Results","Table of Genes")),
                       downloadButton("downloadData", "Download")
            )
            
        )
    )
)

# Define server logic 
server <- function(input, output) {
    
    observeEvent(input$build, {
        results <- noderank::centrality_pipeline(deg = myc_de,
                                                 edge_conf_score_min = 950,
                                                 logFC_min = as.numeric(input$logFC),
                                                 pvalue_max = as.numeric(input$pvalue),
                                                 method = input$method,
                                                 causal_gene_symbol = input$target,
                                                 export_network = FALSE,
                                                 sim_method = 'jaccard',
                                                 n_sim = 9999,
                                                 weighted = as.logical(as.numeric(input$weighted)),
                                                 connected_filter = as.logical(as.numeric(input$connected)))
        
        output$nodePlot <- renderPlot({
            plot_graph(results[['network']], method = 'weighted_score', gene_list = c(input$target))
        })
        
        
        performanceNames <- c(names(results[['performance']])[-6])
        
        output$targetPerformance <- DT::renderDataTable({datatable(results[['performance']]) %>%
                formatRound(performanceNames, 3)})
        
        topGnames <- c(names(results[['top_genes']])[-1:-3])
        
        output$topGenes <- DT::renderDataTable({datatable(results[['top_genes']]) %>%
                formatRound(topGnames, 3)}) 
        
    })
}

# Run the application 
shinyApp(ui = ui, server = server)