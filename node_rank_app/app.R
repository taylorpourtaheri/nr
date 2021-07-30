library(shiny)
library(shinythemes)
library(noderank)
library(DT)

devtools::load_all()

# # save the deg results as an excel file
#de_string <- readRDS('data/de_string_v11.RDS')
#myc_de <- de_string$MYC
#openxlsx::write.xlsx(myc_de, 'data/de_string_v11_MYC.xlsx')

# string_db <- readRDS('data/string_db_v11.RDS')
string_ppi <- readRDS('data/string_ppi_v11.RDS')
id_xref <- readRDS('data/biomart_xreference_ppi_genes.RDS')
sim <- readRDS('data/string_ppi_v11_jacc_sim_list_dense.RDS')
lobstr::mem_used()

# sim_mat <- readRDS('data/string_ppi_v11_jacc_sim_mat.RDS')
# lobstr::mem_used()
# sim_mat <- NULL


# Define UI for application
ui <- fluidPage(

    theme = shinytheme("flatly"),

    #titlePanel(h2("noderank\nA node prioritization tool for differential gene expression analysis", align = 'left')),

    titlePanel(div(HTML("<h2><i><em>noderank</em></i> : A node prioritization tool for differential gene expression analysis</h2>"))),

    # Sidebar
    sidebarLayout(
        sidebarPanel(width = 3,

            fileInput(inputId = 'dge_data',
                       label = 'Upload differential gene expression analysis results',
                       multiple = TRUE,
                       accept = c('.xlsx')),
            textInput(inputId = 'target',
                      label = 'Causal gene'),
            selectInput(inputId = 'method',
                        label = 'Centrality method',
                        choices = c('avg_strength', 'betweenness', 'degree',
                                    'evcent_uw', 'evcent_w', 'strength')),
            sliderInput(inputId = 'logFC',
                        label = 'Minimum log2 fold change',
                        min = 0, max = 3,
                        value = 1.5, step = 0.1),
            sliderInput(inputId = 'pvalue', label = 'Maximum p-value',
                        min = 0, max = 0.25,
                        value = 0.05, step = 0.01),
            selectInput(inputId = 'connected',
                        label = 'Return connected component',
                        choices = c("TRUE"=1,"FALSE"=0)),
            selectInput(inputId = 'weighted',
                        label = 'Return weighted score',
                        choices = c("TRUE"=1,"FALSE"=0)),
            actionButton(inputId = "build",
                         label = 'Rank nodes')
            )
        ,

        mainPanel(
            width = 9,
            navbarPage(title = '',
                       tabPanel('Network Plot',plotOutput('nodePlot')),
                       tabPanel('Method Performance',
                                DT::dataTableOutput('targetPerformance'),
                                plotOutput('simulationPlot')),
                       tabPanel('Ranked Genes',DT::dataTableOutput('topGenes')),
                       tabPanel('Authors',fluidPage(
                           tags$h3('This application was authored by:'),
                           tags$a(href = 'https://www.linkedin.com/in/bradhowlett',
                                  'Bradley Howlett'),
                           tags$p('School of Data Science, University of Virginia'),
                           tags$a(href = 'https://www.linkedin.com/in/taylorhderby',
                                  'Taylor Derby Pourtaheri'),
                           tags$p('School of Data Science, University of Virginia'),
                           tags$a(href = 'https://www.linkedin.com/in/patrick-chatfield',
                                  'Patrick Chatfield'),
                           tags$p('School of Data Science, University of Virginia'),
                           tags$a(href = 'www.linkedin.com/in/monish-dadlani-9423ab191',
                                  'Monish Dadlani'),
                           tags$p('School of Data Science, University of Virginia'))
                       )
            ),

            #move this info out of the navbarPage
            #selectInput(inputId = "dataset",
            #            label = "Choose a dataset:",
            #            choices = c("Method Performance", "Ranked Genes")),
            uiOutput("datasetSelector"),
            uiOutput("downloadButton")
            #downloadButton("downloadData", label = "Download")

        )
    )
)


# Define server logic
server <- function(input, output) {

    observeEvent(input$build, {

        output$datasetSelector <- renderUI({
            selectInput(inputId = "dataset",
                        label = "Choose a dataset:",
                        choices = c("Method Performance", "Ranked Genes"))
        })

        output$downloadButton <- renderUI({
        downloadButton('downloadData',label = 'Download')
    })})

    observeEvent(input$build, {

        # read input file
        file <- input$dge_data
        ext <- tools::file_ext(file$datapath)
        req(file)
        validate(need(ext == "xlsx", "Please upload an xlsx file"))
        data <- openxlsx::read.xlsx(file$datapath)

        # generate results
        results <- noderank::centrality_pipeline(deg = data,
                                                 # string_db = string_db,
                                                 id_xref = id_xref,
                                                 ppi = string_ppi,
                                                 sim = sim,
                                                 edge_conf_score_min = 950,
                                                 logFC_min = as.numeric(input$logFC),
                                                 pvalue_max = as.numeric(input$pvalue),
                                                 method = input$method,
                                                 causal_gene_symbol = input$target,
                                                 export_network = FALSE,
                                                 n_sim = 9999,
                                                 weighted = as.logical(as.numeric(input$weighted)),
                                                 connected_filter = as.logical(as.numeric(input$connected)))

        weightedInput <- reactive({
            switch(input$weighted,
                   "0" = input$method,
                   "1" = 'weighted_score')
        })

        output$nodePlot <- renderPlot({
            plot_graph(results[['network']],
                       method = weightedInput(),
                       gene_list = c(input$target))
            })

        output$simulationPlot <- renderPlot({
            ggplot(data = results[['simulation_scores']],
                   aes(x = simulation_mean_score)) +
                geom_histogram() +
                geom_vline(data = results[['performance']],
                           aes(xintercept = mean_score),
                           color = 'red') +
                labs(title = 'Distribution of simulated network scores',
                     subtitle = 'True network score represented by vertical line',
                     x = 'Mean score',
                     y = NULL)
        })




        performance_display <- results[['performance']] %>%
            .[,c('rank', 'mean_score', 'score_pval', 'z_score', 'sample_mean', 'sample_sd')] %>%
            stats::setNames(c('Gene Rank', 'Mean Network Score', 'Network p-value',
                            'Z Score', 'Null Mean Score', 'Null Standard Deviation'))

        performanceNames <- c(names(performance_display))[-1]

        output$targetPerformance <- DT::renderDataTable({datatable(performance_display) %>%
                formatRound(performanceNames, 3)})

        if (as.logical(as.numeric(input$weighted))){

            top_genes_display <- results[['top_genes']] %>%
                .[,c('Symbol', 'logFC', 'AveExpr', 'P.Value', 'adj.P.Val', input$method, 'causal_similarity', 'weighted_score')] %>%
                stats::setNames(c('Symbol', 'Log Fold Change', 'Average Expression',
                           'p-value', 'Adj. p-value', paste(tools::toTitleCase(input$method),' Score'),
                           'Causal Similarity', 'Weighted Score'))
        } else{

            top_genes_display <- results[['top_genes']] %>%
                .[,c('Symbol', 'logFC', 'AveExpr', 'P.Value', 'adj.P.Val', input$method, 'causal_similarity')] %>%
                stats::setNames(c('Symbol', 'Log Fold Change', 'Average Expression',
                                  'p-value', 'Adj. p-value', paste(tools::toTitleCase(input$method),' Score'),
                                  'Causal Similarity'))
        }

        #create a copy without the hyperlink for download
        top_genes_download <- top_genes_display

        #update with hyperlink
        top_genes_display$Symbol <- paste0("<a href = 'https://www.ncbi.nlm.nih.gov/gene/?term=",top_genes_display$Symbol,"'>",top_genes_display$Symbol,"</a>")

        topGnames <- c(names(top_genes_display)[-1])

        output$topGenes <- DT::renderDataTable({datatable(top_genes_display, escape=FALSE) %>%
                formatRound(topGnames, 3)})



        #download the results
        datasetInput <- reactive({
            switch(input$dataset,
                   "Method Performance" = performance_display,
                   "Ranked Genes" = top_genes_download)
        })

        # file_name <- input$dataset
        # file_string <- tolower(gsub(' ', '_', file_name))

        datasetString <- reactive({
            switch(input$dataset,
                   "Method Performance" = 'method_performance',
                   "Ranked Genes" = 'ranked_genes')
        })

        output$downloadData <- downloadHandler(

            filename = function() {
                paste0(datasetString(),'-',input$target,'-',Sys.Date(), ".csv")
            },
            content = function(file) {
                write.csv(datasetInput(), file)
            }
        )

        })

}

lobstr::mem_used()

# Run the application
shinyApp(ui = ui, server = server)
