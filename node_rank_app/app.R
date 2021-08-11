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
#lobstr::mem_used()

# sim_mat <- readRDS('data/string_ppi_v11_jacc_sim_mat.RDS')
# lobstr::mem_used()
# sim_mat <- NULL

# text inputs for documentation
app_description <- HTML(paste('<i>noderank</i> maps differential gene expression results onto',
                              'an established network of protein-protein associations',
                              'and employs network analysis methods to select a ',
                              'prioritized list of important nodes. For data in which',
                              'the mechanism of perturbation is known, each ranked list',
                              'result is scored by the position of the target gene.'))
head1 <- 'Instructions'
step1 <- paste('Upload an Excel file containing differential gene expression results',
               'with gene identifiers represented in the HGNC format. The required columns',
               'include Symbol, the gene symbols in HGNC format, STRING_id, with ',
               'the STRING identifiers for the genes (ENSG)',
               'logFC, containing log2 fold-changes for each gene, and adj.P.Val,',
               'containing the associated p-values.',
               'Specify the logFC and adjusted p-value thresholds for inclusion',
               'in the final network and select the desired connected_component filter',
               'setting. Optionally, a known or hypothesized target gene',
               '(HGNC symbol) may be entered and',
               'weighted=TRUE selected to generate a network score weighted by',
               'structural similarity to the target gene. Select Build.')
method_description <- 'Method'
head2 <- 'Protein Association Network Construction and Pre-processing'
step2 <- paste('First, the protein-protein interaction (ppl) network is downloaded',
               'from the STRING database.',
               'The differential gene expression analysis results are then mapped',
               'onto the ppi network. The ppi network is then filtered to retain',
               'only the nodes that fall within the range of user-defined thresholds.',
               'Optionally, the giant component of the network may be selected to',
               'remove disconnected nodes that may not be relevant to the biological',
               'mechanisms represented in the DEA results. Finally, the network analysis',
               'algorithms can be implemented to rank the nodes.')
head3 <- 'Network and Method Evaluation'
step3 <- paste('If a target gene is know (or suspected), the effectiveness of',
               'the network in identifying nodes that are important to the',
               'biological mechanism can be evaluated by scoring',
               'the structural similarity of the nodes in the final network to',
               'the known target gene. If the weighted argument is TRUE, a',
               'weighted_score, which represents the product of the centrality',
               'score and the similarity score for each node,',
               'will be calculated. If the weighted option is FALSE, only the',
               'centrality score will be reported.',
               'Finally, the overall workflow is evaluated.',
               'T mean performance score (weighted or unweighted) is calculated',
               'across all nodes in the final network and compared to a null',
               'distribution generated from a user-specified number of draws',
               'of the full ppi network. The overall network score is compared',
               'to the scores of the null hypothesis, allowing for the calculation',
               'of a p-value to represent the significance of the score. ')



# Define UI for application
ui <- fluidPage(

    theme = shinytheme("flatly"),

    #titlePanel(h2("noderank\nA node prioritization tool for differential gene expression analysis", align = 'left')),

    titlePanel(div(HTML("<h2><i><em>noderank</em></i> : A node prioritization tool for differential gene expression analysis</h2>")),
               HTML("noderank</title>")),

    # Sidebar
    sidebarLayout(
        sidebarPanel(width = 3,

                     fileInput(inputId = 'dge_data',
                               label = 'Upload differential gene expression analysis results',
                               multiple = TRUE,
                               accept = c('.xlsx')),
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
                                 label = 'Connected component',
                                 choices = c("TRUE"=1,"FALSE"=0)),
                     textInput(inputId = 'target',
                               label = 'Target gene'),
                     selectInput(inputId = 'weighted',
                                 label = 'Weighted score',
                                 choices = c("FALSE"=0,"TRUE"=1)),
                     actionButton(inputId = "build",
                                  label = 'Rank nodes')
        )
        ,

        mainPanel(
            width = 9,
            navbarPage(title = '',
                       tabPanel('Documentation', fluidPage(
                           tags$h4(app_description),
                           tags$h3(head1),
                           tags$p(step1),
                           tags$h3(method_description),
                           tags$h4(head2),
                           tags$p(step2),
                           tags$h4(head3),
                           tags$p(step3)
                       )),
                       tabPanel('Network Plot',plotOutput('nodePlot')),
                       tabPanel(title = 'Method Performance',
                                textOutput('message'),
                                DT::dataTableOutput('targetPerformance'),
                                plotOutput('simulationPlot')),
                       tabPanel('Ranked Genes',DT::dataTableOutput('topGenes')),
                       tabPanel('Authors',
                                fluidPage(
                                    tags$h3('This application was authored by:'),
                                    tags$a(href = 'https://www.linkedin.com/in/bradhowlett',
                                           'Bradley Howlett', target='_blank', rel='noopener noreferrer'),
                                    tags$p('School of Data Science, University of Virginia'),
                                    tags$a(href = 'https://www.linkedin.com/in/taylorhderby',
                                           'Taylor Derby Pourtaheri', target='_blank', rel='noopener noreferrer'),
                                    tags$p('School of Data Science, University of Virginia'),
                                    tags$a(href = 'https://www.linkedin.com/in/patrick-chatfield',
                                           'Patrick Chatfield', target='_blank', rel='noopener noreferrer'),
                                    tags$p('School of Data Science, University of Virginia'),
                                    tags$a(href = 'https://www.linkedin.com/in/monish-dadlani-9423ab191',
                                           'Monish Dadlani', target='_blank', rel='noopener noreferrer'),
                                    tags$p('School of Data Science, University of Virginia'))
                       )
            ),

            #move this info out of the navbarPage
            #selectInput(inputId = "dataset",
            #            label = "Choose a dataset:",
            #            choices = c("Method Performance", "Ranked Genes")),
            uiOutput("h_line"),
            uiOutput("datasetSelector"),
            uiOutput("downloadButton")
            #downloadButton("downloadData", label = "Download")

        )
    )
)


# Define server logic
server <- function(input, output, session) {

    # run
    observeEvent(input$build, {

        output$h_line <- renderUI({
            hr()
        })

        output$datasetSelector <- renderUI({
            selectInput(inputId = "dataset",
                        label = "Choose data to download:",
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

        print(input$target)

        if(input$target == ''){
            x <- NULL
        } else{
            x <- input$target
        }

        print(is.null(x))

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
                                                 causal_gene_symbol = x,
                                                 export_network = FALSE,
                                                 weighted = as.logical(as.numeric(input$weighted)),
                                                 connected_filter = as.logical(as.numeric(input$connected)))

        # set display tab options - target specified
        # output$cond4 = reactive({
        #     length(results) == 4
        # })
        # outputOptions(output, "cond4", suspendWhenHidden = TRUE)

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

        # if target gene specified
        if (length(results) == 4){
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

            # specify columns for top gene table (weighted or not)
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

                output$message <- NULL
            }



            # if no target gene specified
        } else{

            output$simulationPlot <- NULL
            output$targetPerformance <- NULL

            top_genes_display <- results[['top_genes']] %>%
                .[,c('Symbol', 'logFC', 'AveExpr', 'P.Value', 'adj.P.Val', input$method)] %>%
                stats::setNames(c('Symbol', 'Log Fold Change', 'Average Expression',
                                  'p-value', 'Adj. p-value', paste(tools::toTitleCase(input$method),' Score')))

            output$message <- renderText({'A target gene must be specified to evaluate method performance.'})
        }

        #create a copy without the hyperlink for download
        top_genes_download <- top_genes_display

        #update with hyperlink
        top_genes_display$Symbol <- paste0("<a href = 'https://www.ncbi.nlm.nih.gov/gene/?term=",top_genes_display$Symbol,"' target='_blank', rel='noopener noreferrer'>",top_genes_display$Symbol,"</a>")

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

#lobstr::mem_used()

# Run the application
shinyApp(ui = ui, server = server)
