#' Run interactive STsimulator
#'
#' @import dplyr shiny
#' @return A simulated dataset. Either as individual csv files with counts,
#' cell features, or a Giotto object.
#'
#' @export
#'
run_interactive_sCCIgen <- function() {

  ui <- miniUI::miniPage(
    miniUI::gadgetTitleBar("sCCIgen"),
    miniUI::miniContentPanel(
      shiny::tabsetPanel(
        type = "tabs",

        shiny::tabPanel("Download a pre-simulated dataset",
                        shiny::radioButtons(inputId = "dataset",
                                            label = "Select a pre-simulated spatial transcriptomics dataset to download.
                                            These pre-simulated data include 5% spatial differential expressed genes in one cell
                                            type, 5% interaction changed genes in one cell type that are associated with the
                                            proximity of another cell type, and 5% interaction changed genes in cell-type pair that
                                            are associated with each other in the neighboring cells.",
                                            choices = c("Simulated data for 4,751 genes in 4,000 cells of 6 cell types
                                                        in two regions on a unit square, using normal breast data profiled by
                                                        snRNAseq as reference." = "example1",
                                                        "Simulated data for 10,000 genes in 500 spots of 6 cell types, using
                                                        mouse brain data profiled by SeqFISH+ as reference." = "example2",
                                                        "Simulated data for 550 genes in 10,000 cells of 6 cell types, using
                                                        human ovarian cancer data profiled by MERFISH as reference." = "example3"),
                                            width = "100%",
                                            selected = character(0)
                        ), # ends radioButtons for dataset

                        shiny::downloadButton("downloadCount", "Download Counts"),
                        shiny::downloadButton("downloadMeta", "Download Metadata"),
                        shiny::downloadButton("downloadExpr", "Download Expression Pattern"),
                        shiny::downloadButton("downloadParameter", "Download Parameter file")
        ), # ends Download a pre-simulated dataset tab

        shiny::tabPanel("Create the parameter file",
                        shiny::tabsetPanel(
                          type = "tabs",
                          shiny::tabPanel("Step 1: Select your dataset",
                                          shiny::radioButtons(
                                            inputId = "inputdata",
                                            label = "sCCIgen is a real-data based simulator. It accepts a
                                                    variety of input data types including scRNAseq, snRNAseq, or
                                                    single-cell spatially resolved transcriptomics. The simulator
                                                    provides users access to several publicly available datasets
                                                    and also allows them to use their own data, offering flexibility
                                                    in the simulations.

                                                    What data do you want to use as the base for simulation? If the data
                                                    is not yet in your working directory, use the Download buttons
                                                    to get the files in your computer.",
                                            choices = c("Decoy data 1: It includes (1) count matrix for 10 genes by
                                                        1000 cells of 2 cell types, and (2) cell feature matrix for
                                                        annotated cell type." = "fake1",
                                                        "Decoy data 2: It includes (1) count matrix for 10 genes by
                                                        1000 cells of 2 cell types, and (2) cell feature matrix for
                                                        annotated cell type and spatial coordinate." = "fake2",
                                                        "Decoy data 3: It includes (1) count matrix for 10 genes by
                                                        1000 cells of 2 cell types, and (2) cell feature matrix for
                                                        annotated cell type, spatial coordinate, and region." = "fake3",
                                                        "Normal human breast snRNAseq data: It includes (1) count matrix
                                                        for 4751 genes by 5990 cells of 6 cell types (epithelial cell,
                                                        adipocyte, fibroblast, endothelial cell, immune (myeloid) and
                                                        muscle), and (2) cell feature matrix for annotated cell type.
                                                        PMID: 35549429" = "snRNAseq_breast_033023",
                                                        "Normal mouse brain SeqFISH+ data: It includes (1) count matrix
                                                        for 10,000 genes by 511 cells of 6 cell types (excitatory neuron,
                                                        interneuron, astrocyte, microglia, oligodendrocyte and
                                                        endothelial cells), and (2)cell feature matrix including cell
                                                        type annotation and spatial coordinate on 2D (x, y).
                                                        PMID: 35549429" = "SeqFishPlusCortex_033023",
                                                        "Ovarian cancer MERFISH data: It includes (1) count matrix
                                                        for 550 genes by 355,633 cells of 6 cell types (tumor, adipose,
                                                        endo-immune, macro-immune, macrophage, and others), and (2) cell
                                                        feature matrix including cell type annotation and spatial
                                                        coordinate on 2D. Data source: vizgen" = "MERFISH_OV",
                                                        "I want to use my own dataset" = "user_input"),
                                            width = "100%",
                                            selected = character(0)
                                          ), # ends radioButtons for inputdata

                                          shiny::uiOutput("downloadExpression"),

                                          shiny::uiOutput("downloadFeature"),

                                          shiny::uiOutput("userinputexpression_text"),

                                          shiny::HTML(strrep(htmltools::br(), 1)),

                                          shiny::uiOutput("userinputexpression"),

                                          shiny::HTML(strrep(htmltools::br(), 1)),

                                          shiny::uiOutput("user_expression_text"),

                                          shiny::HTML(strrep(htmltools::br(), 1)),

                                          shiny::uiOutput("userinputcellfeature_text"),

                                          shiny::HTML(strrep(htmltools::br(), 1)),

                                          shiny::uiOutput("userinputcellfeature"),

                                          shiny::HTML(strrep(htmltools::br(), 1)),

                                          shiny::uiOutput("user_cellfeature_text"),

                                          shiny::HTML(strrep(htmltools::br(), 1)),

                                          shiny::uiOutput("text_ask_data_directory"),

                                          shiny::HTML(strrep(htmltools::br(), 1)),

                                          shiny::uiOutput("button_ask_data_directory"),

                                          shiny::HTML(strrep(htmltools::br(), 1)),

                                          shiny::uiOutput("path_to_input_dir"),

                                          shiny::HTML(strrep(htmltools::br(), 1)),

                                          shiny::uiOutput("button_read_data"),

                                          shiny::HTML(strrep(htmltools::br(), 2))

                          ), # ends Step 1: Select your dataset

                          shiny::tabPanel("Step 2: Select the simulation parameters",


                                          shiny::HTML(strrep(htmltools::br(), 2)),

                                          shiny::uiOutput("ask_simulate_cells"),

                                          shiny::uiOutput("ask_numbercells"),

                                          shiny::uiOutput("ask_numberregions"),

                                          shiny::uiOutput("ask_custom_cell_type_prop"),

                                          shiny::uiOutput("customcellproportionsSelection"),

                                          shiny::uiOutput("button_customcellproportions"),

                                          shiny::uiOutput("input_file_customcellproportions"),

                                          shiny::uiOutput("text_file_customcellproportions"),

                                          shiny::uiOutput("button_read_file_customcellproportions"),

                                          shiny::uiOutput("ask_custom_interaction"),

                                          shiny::uiOutput("custominteractionSelection"),

                                          shiny::uiOutput("button_custominteraction"),

                                          shiny::uiOutput("button_select_interaction_file"),

                                          shiny::HTML(strrep(htmltools::br(), 1)),

                                          shiny::uiOutput("text_custom_interaction_file"),

                                          shiny::uiOutput("button_read_custom_interaction_file"),

                                          shiny::HTML(strrep(htmltools::br(), 1)),

                                          shiny::uiOutput("ask_even_distribution"),

                                          shiny::uiOutput("ask_windowmethod"),

                                          shiny::uiOutput("ask_overlapcutoff"),

                                          # Parameters for expression profiles

                                          shiny::numericInput(inputId = "depthratio",
                                                              label = "Select the sequencing depth ratio between
                                                              simulated and reference data. By default, simulated data
                                                              have the same depth as the input data; an alteration of this
                                                              parameter is a way to generate batch effects for different simulations.",
                                                              value = 1,
                                                              width = "100%"),

                                          shiny::uiOutput("ask_model_per_region"),

                                          shiny::radioButtons(inputId = "mimiccorrelation",
                                                              label = "Do you want to mimic gene-gene correlation of
                                                              the reference data? Select 'Yes' only if the gene-gene
                                                              correlation is pre-estimated for each cell type.
                                                              STsimulator provides functions/pipelines for estimating
                                                              and saving gene-gene correlation, which can be served
                                                              as input here.",
                                                              choices = c("No, simulate independent genes" = FALSE,
                                                                          "Yes" = TRUE),
                                                              width = "100%"),

                                          shiny::uiOutput("mimiccorrelationSelection"),

                                          shiny::uiOutput("ask_spatialpatterns"),

                                          shiny::uiOutput("spatialpatternsSelection"),

                                          shiny::uiOutput("button_spatialpatterns"),

                                          shiny::uiOutput("button_select_spatialpatterns_file"),

                                          shiny::uiOutput("text_spatialpatterns_file"),

                                          shiny::uiOutput("button_read_spatialpatterns_file"),

                                          shiny::HTML(strrep(htmltools::br(), 1)),

                                          shiny::uiOutput("ask_addcellcellinteraction"),

                                          shiny::uiOutput("addcellcellinteractionSelection"),

                                          shiny::uiOutput("button_addcellcellinteraction"),

                                          shiny::uiOutput("button_select_addcellcellinteractionexpression"),

                                          shiny::uiOutput("text_file_addcellcellinteractionexpression"),

                                          shiny::uiOutput("button_read_file_addcellcellinteractionexpression"),

                                          shiny::HTML(strrep(htmltools::br(), 1)),

                                          shiny::uiOutput("ask_addcellcellinteractionexpression"),

                                          shiny::uiOutput("addcellcellinteractionexpressionSelection1"),

                                          shiny::uiOutput("addcellcellinteractionexpressionSelection2"),

                                          shiny::uiOutput("button_addcellcellinteractionexpression"),

                                          shiny::uiOutput("button_select_cellcellinteractionexpression_neighborhood"),

                                          shiny::uiOutput("text_cellcellinteractionexpression_neighborhood_file"),

                                          shiny::uiOutput("button_read_cellcellinteractionexpression_neighborhood_file"),

                                          shiny::HTML(strrep(htmltools::br(), 1)),

                                          shiny::radioButtons(inputId = "multicell",
                                                              label = "What resolution would you like to use to simulate your data?",
                                                              choices = c("single-cell resolution" = FALSE,
                                                                          "multi-cell resolution" = TRUE),
                                                              width = "100%"),

                                          shiny::uiOutput("multicellSelection"),

                                          shiny::numericInput(inputId = "ndatasets",
                                                              label = "Select the number of simulated data sets",
                                                              value = 1,
                                                              width = "100%"),

                                          shiny::uiOutput("ask_single_seed"),

                                          shiny::textInput(inputId = "outdir",
                                                           label = "Provide the name of the output folder. This folder will be created under your working directory.",
                                                           value = "output_files",
                                                           width = "100%"),

                                          shiny::textInput(inputId = "outname",
                                                           label = "Type the root name of the output files.
                                                           (e.g. example1 will create example1_count_1.tsv, example1_expr_pattern_1.tsv, and
                                                           example1_meta_1.tsv)",
                                                           value = "output",
                                                           width = "100%"),

                                          shiny::textInput(inputId = "parameter_filename",
                                                           label = "Type the name of the parameter file.",
                                                           value = "parameter_file.tsv",
                                                           width = "100%"),


                                          shiny::downloadButton("downloadparams", "Create and Download my parameter file")

                          ) # ends Step 2: Select the simulation parameters



                        ) # ends tabsetPanel inside Run a simulation


        ), # ends Create parameter file
        shiny::tabPanel("Run a simulation",
                        shiny::textInput(inputId = "inputparamfile",
                                         label = "By default, we will look for a parameter_file.tsv file. If your file has a different name, type it here:",
                                         width = "100%",
                                         value = "parameter_file.tsv"),

                        shiny::radioButtons(inputId = "createGiotto",
                                            label = "By the default, the simulation will create the expression, metadata, and pattern output files. Optionally, you can also create a Giotto object.",
                                            choices = c("Run the simulation and create the default output files" = FALSE,
                                                        "Run the simulation and create a Giotto object" = TRUE),
                                            width = "100%"),

                        shiny::uiOutput("createGiottoSelection"),

                        shiny::actionButton(inputId = "runsimulation",
                                            label = "Run simulation",
                                            width = "100%")

        ) # ends Run a simulation


      ) # ends tabsetPanel



    ) # ends miniContentPanel
  ) # ends ui

  server <- function(input, output,session) {

    # Download pre-simulated data
    shiny::observeEvent(input$dataset, {

      output$downloadCount <- shiny::downloadHandler(
        filename = paste0(input$dataset,"_count.tsv"),
        content = function(con) {
          download.file(url = paste0("https://github.com/songxiaoyu/sCCIgen_data/raw/main/example_data/",input$dataset,"_count.tsv"),
                        destfile = con)
        }
      )

      output$downloadMeta <- shiny::downloadHandler(
        filename = paste0(input$dataset,"_meta.tsv"),
        content = function(con) {
          download.file(url = paste0("https://github.com/songxiaoyu/sCCIgen_data/raw/main/example_data/",input$dataset,"_meta.tsv"),
                        destfile = con)
        }
      )

      output$downloadExpr <- shiny::downloadHandler(
        filename = paste0(input$dataset,"_expr_pattern.tsv"),
        content = function(con) {
          download.file(url = paste0("https://github.com/songxiaoyu/sCCIgen_data/raw/main/example_data/",input$dataset,"_expr_pattern.tsv"),
                        destfile = con)
        }
      )

      output$downloadParameter <- shiny::downloadHandler(
        filename = paste0(input$dataset,"_parameter.tsv"),
        content = function(con) {
          download.file(url = paste0("https://github.com/songxiaoyu/sCCIgen_data/raw/main/example_data/",input$dataset,"_parameter.tsv"),
                        destfile = con)
        }
      )

    }) # ends download pre-simulated data


    # Provide or download Rdata

    expression_data_file <- shiny::reactiveVal()
    cell_feature_data_file <- shiny::reactiveVal()
    path_to_input_dir <- shiny::reactiveVal()

    shiny::observeEvent(input$inputdata, {

      if(input$inputdata == "user_input") {

        output$userinputexpression_text <- renderUI({

          shiny::strong("Select the expression file (e.g. expression.Rdata, .tsv, or .csv).
                        Upload a G gene by N cell matrix for expression count data.
                        Instruction: Row names should be unique identifiers of genes.")
        })

        output$userinputexpression <- renderUI({

          shinyFiles::shinyFilesButton(id = "user_expression",
                                       label = "Select the expression file",
                                       title = "Please select a expression file",
                                       multiple = FALSE,
                                       viewtype = "list")
        })

        shinyFiles::shinyFileChoose(input,
                                    id = 'user_expression',
                                    roots = c(wd = getwd()) )


        shiny::observeEvent(input$user_expression, {
          x <- shiny::reactiveVal()
          x_length <- shiny::reactiveVal()

          x(gsub("/wd", "", paste(unlist(input$user_expression), collapse = "/")))
          x_length(length(unlist(input$user_expression)))

          if(x_length() >= 2) {
            expression_data_file(unlist(input$user_expression)[x_length() - 1])
            x_length(length(unlist(input$user_expression)) - 2)
            path_to_input_dir(paste0(getwd(),
                              paste(unlist(input$user_expression)[1:x_length()], collapse = "/")))
          }

        })

        output$user_expression_text <- renderUI({

          shiny::strong({paste("Your expression file is:", expression_data_file())})
        })

        output$userinputcellfeature_text <- renderUI({

          shiny::strong("Select the cellfeature file (e.g. cellfeature.Rdata, .tsv, or .csv).
                Upload a N by K matrix for cell feature data.
                Instruction: Column and row names are not expected.
                Cells should be in the same order as the uploaded expression count data.
                The first column is required, which should be cell type annotation.
                Other columns are optional. If spatial coordinates on x, y axes are provided,
                they should be the 2-3 columns of the input cell feature data.
                If a spatial region variable is provided, it should be the 4th column of the data.",
          )
        })

        output$userinputcellfeature <- renderUI({

          shinyFiles::shinyFilesButton(id = "user_cellfeature",
                                       label = "Select the cell feature file",
                                       title = "Please select a cell feature file",
                                       multiple = FALSE,
                                       viewtype = "list")

        })

        shinyFiles::shinyFileChoose(input,
                                    id = 'user_cellfeature',
                                    roots = c(wd = getwd()) )

        shiny::observeEvent(input$user_cellfeature, {
          x <- shiny::reactiveVal()
          x_length <- shiny::reactiveVal()

          x(gsub("/wd", "", paste(unlist(input$user_cellfeature), collapse = "/")))
          x_length(length(unlist(input$user_cellfeature)))

          if(x_length() >= 2) {
            cell_feature_data_file(unlist(input$user_cellfeature)[x_length() - 1])
          }
        })

        output$user_cellfeature_text <- renderUI({
          shiny::strong({paste("Your cell feature file is:", cell_feature_data_file())})
        })

        output$downloadExpression <- NULL
        output$downloadFeature <- NULL
        output$text_ask_data_directory <- NULL
        output$button_ask_data_directory <- NULL

      } else {

        output$userinputexpression_text <- NULL
        output$userinputexpression <- NULL
        output$user_expression_text <- NULL

        output$userinputcellfeature_text <- NULL
        output$userinputcellfeature <- NULL
        output$user_cellfeature_text <- NULL

        # show download buttons
        output$downloadExpression <- renderUI({
          shiny::downloadButton("downloadselectedexpression", "Download selected expression data")
        })

        output$downloadFeature <- renderUI({
          shiny::downloadButton("downloadselectedcellfeature", "Download selected cellfeature data")
        })

        output$downloadselectedexpression <- shiny::downloadHandler(
          filename = paste0(input$inputdata,"_expr.Rdata"),
          content = function(con) {
            download.file(url = paste0("https://github.com/songxiaoyu/sCCIgen_data/raw/main/InputData/expression_data/",input$inputdata,"_expr.Rdata"),
                          destfile = con)
          })

        output$downloadselectedcellfeature <- shiny::downloadHandler(
          filename = paste0(input$inputdata,"_cellfeature.Rdata"),
          content = function(con) {
            download.file(url = paste0("https://github.com/songxiaoyu/sCCIgen_data/raw/main/InputData/cell_feature_data/",input$inputdata,"_cellfeature.Rdata"),
                          destfile = con)
          })

        expression_data_file(paste0(input$inputdata,"_expr.Rdata"))
        cell_feature_data_file(paste0(input$inputdata,"_cellfeature.Rdata"))

        output$text_ask_data_directory <- renderUI({
          shiny::strong("Where is your input data located? Click the button to select your data directory.")
        })

        output$button_ask_data_directory <- renderUI({
          shinyFiles::shinyDirButton(id = "inputdir",
                                     label = "Select the data folder",
                                     title = "Please select a folder",
                                     multiple = FALSE,
                                     viewtype = "list",
                                     filename = "."
          )
        })

        shinyFiles::shinyDirChoose(input,
                                   id = 'inputdir',
                                   roots = c(wd = getwd()),
                                   allowDirCreate = FALSE)

        shiny::observeEvent(input$inputdir, {
          x <- shiny::reactiveVal()
          x(gsub("/wd", "", paste(unlist(input$inputdir), collapse = "/")))

          path_to_input_dir(paste0(getwd(), x()))
        })

      }

      output$path_to_input_dir <- renderUI({
        shiny::strong({paste("Your selected directory is:", path_to_input_dir() )})
      })

    }) # ends download data

    shiny::observeEvent(path_to_input_dir(), {
      if (!is.null(path_to_input_dir())) {
        output$button_read_data <- shiny::renderUI({
          shiny::actionButton(inputId = "use_inputdata",
                              label = "Use this directory and read my data")
        })
      }
    })

    # Create parameters file

    expression_data_file_type <- shiny::reactiveVal()
    expression_data_cell_types <- shiny::reactiveVal()
    expr_data <- shiny::reactiveVal()
    cellfeature_data <- shiny::reactiveVal()
    ncol_feature_data <- shiny::reactiveVal()

    shiny::observeEvent(input$use_inputdata, {

      if(file.exists(fs::path(path_to_input_dir(), expression_data_file())) &
         file.exists(fs::path(path_to_input_dir(), cell_feature_data_file())) ) {

        expression_data_file_type(tools::file_ext(expression_data_file()))

        if(expression_data_file_type() == "Rdata" ||
           expression_data_file_type() == "RData") {

          load(fs::path(path_to_input_dir(), expression_data_file()), x <- new.env())
          expr_data(get(ls(x), envir = x))
          rm(x)

          load(fs::path(path_to_input_dir(), cell_feature_data_file()), x <- new.env())
          y = get(ls(x), envir = x)
          rm(x)

          colnames(y)[1] = "cell_type"
          cellfeature_data(y)

        } else if(expression_data_file_type() == "csv") {

          expr_data(read.csv(fs::path(path_to_input_dir(), expression_data_file()),
                             row.names = 1))

          cellfeature_data(read.csv(fs::path(path_to_input_dir(), cell_feature_data_file()),
                                    row.names = 1))

        } else if(expression_data_file_type() == "tsv") {

          expr_data(read.delim(fs::path(path_to_input_dir(), expression_data_file()),
                               row.names = 1))

          cellfeature_data(read.delim(fs::path(path_to_input_dir(), cell_feature_data_file()),
                                      row.names = 1))
        }

        expression_data_cell_types(paste(unique(sort(unique(cellfeature_data()[,1]) )),
                                         collapse = ","))

        ncol_feature_data(ncol(cellfeature_data()))

        shiny::showModal(shiny::modalDialog(
          title = "Successful reading of input files",
          "Your input files were successfully loaded.",
          footer = modalButton("OK")
        ))

      } else {

        shiny::showModal(shiny::modalDialog(
          title = "Error reading input files",
          "At least one of your input files was not found in your input directory.
          Verify that the expression and cellfeature files are in the right directory.",
          footer = modalButton("OK")
        ))
      }

    })

    simulate_spatial_data <- shiny::reactiveVal()

    window_method <- shiny::reactiveVal()

    num_regions <- shiny::reactiveVal()
    custom_cell_type_proportions <- shiny::reactiveVal()
    cell_type_proportions <- shiny::reactiveVal()
    custom_props <- shiny::reactiveVal()
    custom_cell_location_interactions <- shiny::reactiveVal()
    cell_location_interactions_df <- shiny::reactiveVal()
    cell_even_distribution <- shiny::reactiveVal()
    region_specific_model <- shiny::reactiveVal()
    region_specific_model("NULL")

    shiny::observeEvent(ncol_feature_data(), {

      if(ncol_feature_data() > 1) {

        output$ask_simulate_cells <- shiny::renderUI({
          shiny::radioButtons(inputId = "simulatecells",
                              label = "Do you want to simulate new cells?",
                              choices = c("Do not simulate new cells, but simulate new expression data
                                          for existing cells" = FALSE,
                                          "I want to simulate new cells" = TRUE),
                              width = "100%")
        })

        shiny::observeEvent(input$simulatecells, {

          simulate_spatial_data(input$simulatecells)

          if(input$simulatecells == TRUE) {

            output$ask_windowmethod <- shiny::renderUI({
              shiny::radioButtons(inputId = "windowmethod",
                                  label = "Select the method for determining the window on existing ST data",
                                  choices = c("network" = "network",
                                              "rectangle" = "rectangle",
                                              "convex" = "convex",
                                              "convex2" = "convex2",
                                              "convex3" = "convex3",
                                              "convex5" = "convex5"),
                                  width = "100%")
            })

            shiny::observeEvent(input$windowmethod, {
              window_method(input$windowmethod)
            })

            output$ask_spatialpatterns <- shiny::renderUI({
              shiny::radioButtons(inputId = "spatialpatterns",
                                  label = "Would you like to add spatial patterns?",
                                  choices = c("No" = FALSE,
                                              "Yes" = TRUE),
                                  width = "100%")
            })

            output$ask_custom_interaction <- shiny::renderUI({
              shiny::radioButtons(inputId = "custominteraction",
                                  label = "Users can select cell type pairs and determine the strength. Strength < 0 indicates cell-cell inhibition; and strength > 0 indicates cell-cell attraction. Would you like to specify the cell-cell location interaction (inhibition/attraction)?",
                                  choices = c("No cell-cell location interaction" = FALSE,
                                              "Specify cell-cell location interaction (inhibition/attraction)" = TRUE),
                                  width = "100%")
            })

            shiny::observeEvent(input$custominteraction, {

              custom_cell_location_interactions(input$custominteraction)

              if(input$custominteraction == TRUE) {

                output$custominteractionSelection <- shiny::renderUI({
                  shiny::textInput(inputId = "locationinteraction",
                                   label = "To specify cell-cell location interaction, enter the cell type pairs
                               followed by the interaction level (suggested value -2 to 2) in a format
                               <cell_type_A>,<cell_type_B>,<value>. Separate the entries by blank space
                               (e.g. 'cell_type_A,cell_type_B,1.2 cell_type_B,cell_type_C,-0.8').

                               Alternatively, select a file separated by commas, where each line is an
                               entry with the format described above.",
                                   width = "100%")
                })

                output$button_custominteraction <- shiny::renderUI({
                  shiny::actionButton("save_custominteraction",
                                      "Use these interactions")
                })

                shiny::observeEvent(input$save_custominteraction, {

                  shiny::observeEvent(input$locationinteraction, {
                    x = unlist(stringr::str_split(input$locationinteraction,
                                                  pattern = " "))

                    for (i in x) {
                      y = unlist(stringr::str_split(i,pattern = ","))
                      if(length(y) != 3) {
                        shiny::showModal(shiny::modalDialog(
                          title = "Error",
                          "Your custom cell-cell location interactions don't have
                      the right format.",
                          footer = modalButton("OK")
                        ))
                      }
                    }

                    cell_location_interactions_df(data.frame(parameters = paste0("cell_interaction_",
                                                                                 1:length(x)),
                                                             value = x))
                  })

                  shiny::showModal(shiny::modalDialog(
                    title = "Success!",
                    "Your custom cell-cell location interactions were saved.",
                    footer = modalButton("OK")
                  ))

                })

                output$button_select_interaction_file <- shiny::renderUI({

                  shinyFiles::shinyFilesButton(id = "file_custom_interaction_file",
                                               label = "Select the cell-cell location interaction file",
                                               title = "Please select a cell-cell location interaction file",
                                               multiple = FALSE,
                                               viewtype = "list")
                })

                shinyFiles::shinyFileChoose(input,
                                            id = "file_custom_interaction_file",
                                            roots = c(wd = getwd()) )

                shiny::observeEvent(input$file_custom_interaction_file, {

                  if(!is.null(input$file_custom_interaction_file)) {
                    x <- paste0(getwd(),
                                gsub("/wd", "",
                                     paste(unlist(input$file_custom_interaction_file),
                                           collapse = "/")))

                    output$text_custom_interaction_file <- shiny::renderUI({
                      shiny::strong({paste("Your selected file is:", x)})
                    })

                    output$button_read_custom_interaction_file <- shiny::renderUI({
                      shiny::actionButton("button_read_custom_interaction_file",
                                          "Read this interaction file")
                    })

                    shiny::observeEvent(input$button_read_custom_interaction_file, {
                      x_df <- read.delim(paste0(getwd(),
                                                gsub("/wd", "",
                                                     paste(unlist(input$file_custom_interaction_file),
                                                           collapse = "/"))),
                                         header = TRUE, sep = ",")

                      if(ncol(x_df) == 3) {
                        x_df$join <- paste(x_df$cell_type_A, x_df$cell_type_B, x_df$value,
                                           sep = ",")

                        cell_location_interactions_df(data.frame(parameters = paste0("cell_interaction_",
                                                                                     1:nrow(x_df)),
                                                                 value = x_df$join))

                        shiny::showModal(shiny::modalDialog(
                          title = "Success!",
                          "Your custom cell location interactions were saved.",
                          footer = modalButton("OK")
                        ))
                      } else {
                        shiny::showModal(shiny::modalDialog(
                          title = "Error",
                          "Your file doesn’t have the right format.",
                          footer = modalButton("OK")
                        ))
                      }

                    })
                  }

                })

              } else {
                output$custominteractionSelection <- NULL
                output$button_custominteraction <- NULL
                output$button_select_interaction_file <- NULL
                output$text_custom_interaction_file <- NULL
                output$button_read_custom_interaction_file <- NULL
              }
            })

          } else {
            output$ask_windowmethod <- NULL
            output$ask_spatialpatterns <- NULL
            output$ask_custom_interaction <- NULL
          }
        })

        num_regions("NULL")

        if(ncol_feature_data() > 3) {
          unique_regions = unique(sort(cellfeature_data()[,4]))

          num_regions(length(unique_regions))

          if(length(unique_regions) > 1) {
            output$ask_model_per_region <- shiny::renderUI({
              shiny::radioButtons(inputId = "model_per_region",
                                  label = "Do you want to model the input expression data
                                  separately for each region?",
                                  choices = c("No" = FALSE,
                                              "Yes" = TRUE),
                                  width = "100%")
            })

            shiny::observeEvent(input$model_per_region, {
              region_specific_model(input$model_per_region)
            })
          } else {
            output$ask_model_per_region <- NULL
          }
        }

        output$ask_numberregions <- NULL
        output$ask_custom_cell_type_prop <- NULL
        output$ask_even_distribution <- NULL
        output$ask_custom_interaction <- NULL
      } else { # when ncol == 1
        output$ask_model_per_region <- NULL
        output$ask_simulate_cells <- NULL
        output$ask_windowmethod <- NULL

        region_specific_model("NULL")
        simulate_spatial_data(TRUE)

        output$ask_numberregions <- shiny::renderUI({
          shiny::numericInput(inputId = "numberregions",
                              label = "Enter the number of regions (suggested: 1-10)",
                              value = 2,
                              width = "100%")
        })

        shiny::observeEvent(input$numberregions, {

          if(!is.null(input$numberregions)) {
            num_regions(input$numberregions)
          } else {num_regions("NULL")}
        })

        output$ask_custom_cell_type_prop <- shiny::renderUI({
          shiny::radioButtons(inputId = "customcellproportions",
                              label = "Choose an option for cell-type proportions in each region.",
                              choices = c("Use all equals to cell-type proportions of the input expression data" = FALSE,
                                          "I want to input custom cell-type proportions" = TRUE),
                              width = "100%"
          )
        })

        shiny::observeEvent(input$customcellproportions, {

          custom_props(input$customcellproportions)

          if(input$customcellproportions == TRUE) { # if TRUE, use custom cell proportions

            output$customcellproportionsSelection <- shiny::renderUI({
              shiny::textInput(inputId = "cellproportions",
                               label = "Enter the custom cell-type proportions in a format
                               <region>,<cell_type>,<proportion>. Separate multiple entries by
                               blank space (e.g '1,cell_type_A,0.25 1,cell_type_B,0.75 2,cell_type_C,1.0').

                               Alternatively, select a file separated by commas, where each line is an
                               entry with the format described above.

                               NOTE: Cell type proportions for a given region must sum up to 1.",
                               width = "100%")
            })

            output$button_customcellproportions <- shiny::renderUI({
              shiny::actionButton("save_customcellproportions", "Use these cell proportions")
            })

            shiny::observeEvent(input$save_customcellproportions, {

              shiny::observeEvent(input$cellproportions, {
                x = unlist(stringr::str_split(string = input$cellproportions,
                                              pattern = " "))

                for(i in x) {
                  y = unlist(stringr::str_split(string = i, pattern = ","))
                  if(length(y) != 3) {
                    shiny::showModal(shiny::modalDialog(
                      title = "Error",
                      "Your custom cell proportions don't have the right format.",
                      footer = modalButton("OK")
                    ))
                  }
                }

                n_total = length(x)
                cell_type_proportions(data.frame(parameters = paste0("cell_type_proportion_", 1:n_total),
                                                 value = x))
              })

              shiny::showModal(shiny::modalDialog(
                title = "Success!",
                "Your custom cell type proportions were saved.",
                footer = modalButton("OK")
              ))
            })

            output$input_file_customcellproportions <- renderUI({

              shinyFiles::shinyFilesButton(id = "file_customcellproportions",
                                           label = "Select the cell proportion file",
                                           title = "Please select a cell proportion file",
                                           multiple = FALSE,
                                           viewtype = "list")

            })

            shinyFiles::shinyFileChoose(input,
                                        id = "file_customcellproportions",
                                        roots = c(wd = getwd()) )

            shiny::observeEvent(input$file_customcellproportions, {

              x <- shiny::reactiveVal()

              if(!is.null(input$file_customcellproportions)) {
                x(paste0(getwd(),
                         gsub("/wd", "", paste(unlist(input$file_customcellproportions), collapse = "/"))))

                output$text_file_customcellproportions <- shiny::renderUI({
                  shiny::strong({paste("Your selected file is:", x())})
                })

                output$button_read_file_customcellproportions <- shiny::renderUI({
                  shiny::actionButton("read_file_customcellproportions", "Read this file")
                })

                shiny::observeEvent(input$read_file_customcellproportions,{
                  x_df <- read.delim(paste0(getwd(),
                                            gsub("/wd", "", paste(unlist(input$file_customcellproportions), collapse = "/"))),
                                     header = TRUE, sep = ",")

                  if(ncol(x_df) == 3) {
                    colnames(x_df) <- c("region", "annotation", "proportion")
                    x_df$join <- paste(x_df$region, x_df$annotation, x_df$proportion,
                                        sep = ",")

                    cell_type_proportions(data.frame(
                      parameters = paste0("cell_type_proportion_", 1:nrow(x_df)),
                      value = x_df$join))

                    shiny::showModal(shiny::modalDialog(
                      title = "Success!",
                      "Your custom cell type proportions were saved.",
                      footer = modalButton("OK")
                    ))

                  } else {
                    shiny::showModal(shiny::modalDialog(
                      title = "Error",
                      "Your file doesn’t have the right format.",
                      footer = modalButton("OK")
                    ))
                  }

                })

              }

            })

          } else {

            output$customcellproportionsSelection <- NULL
            output$button_customcellproportions <- NULL
            output$input_file_customcellproportions <- NULL
            output$text_file_customcellproportions <- NULL
            output$button_read_file_customcellproportions <- NULL

            x = cellfeature_data() %>%
              count(cell_type) %>%
              mutate(proportions = round(n/sum(n),3),
                     proportions2 = paste0(cell_type,",",proportions))

            shiny::observeEvent(input$numberregions, {

              if(!is.null(input$numberregions)) {
                n_total = input$numberregions*length(x$proportions2)

                cell_type_proportions(
                  data.frame(parameters = paste0("cell_type_proportion_", 1:n_total),
                             value = paste0(base::sort(rep(seq_len(input$numberregions),
                                                           length(x$proportions2))),
                                            ",", x$proportions2)) )
              }
            })

          }

        })

        output$ask_custom_interaction <- shiny::renderUI({
          shiny::radioButtons(inputId = "custominteraction",
                              label = "Users can select cell type pairs and determine the strength. Strength < 0 indicates cell-cell inhibition; and strength > 0 indicates cell-cell attraction. Would you like to specify the cell-cell location interaction (inhibition/attraction)?",
                              choices = c("No cell-cell location interaction" = FALSE,
                                          "Specify cell-cell location interaction (inhibition/attraction)" = TRUE),
                              width = "100%")
        })

        shiny::observeEvent(input$custominteraction, {

          custom_cell_location_interactions(input$custominteraction)

          if(input$custominteraction == TRUE) {

            output$custominteractionSelection <- shiny::renderUI({
              shiny::textInput(inputId = "locationinteraction",
                               label = "To specify cell-cell location interaction, enter the cell type pairs
                               followed by the interaction level (suggested value -2 to 2) in a format
                               <cell_type_A>,<cell_type_B>,<value>. Separate the entries by blank space
                               (e.g. 'cell_type_A,cell_type_B,1.2 cell_type_B,cell_type_C,-0.8').

                               Alternatively, select a file separated by commas, where each line is an
                               entry with the format described above.",
                               width = "100%")
            })

            output$button_custominteraction <- shiny::renderUI({
              shiny::actionButton("save_custominteraction",
                                  "Use these interactions")
            })

            shiny::observeEvent(input$save_custominteraction, {

              shiny::observeEvent(input$locationinteraction, {
                x = unlist(stringr::str_split(input$locationinteraction,
                                              pattern = " "))

                for (i in x) {
                  y = unlist(stringr::str_split(i,pattern = ","))
                  if(length(y) != 3) {
                    shiny::showModal(shiny::modalDialog(
                      title = "Error",
                      "Your custom cell-cell location interactions don't have
                      the right format.",
                      footer = modalButton("OK")
                    ))
                  }
                }

                cell_location_interactions_df(data.frame(parameters = paste0("cell_interaction_",
                                                                             1:length(x)),
                                                         value = x))
              })

              shiny::showModal(shiny::modalDialog(
                title = "Success!",
                "Your custom cell-cell location interactions were saved.",
                footer = modalButton("OK")
              ))

            })

            output$button_select_interaction_file <- shiny::renderUI({

              shinyFiles::shinyFilesButton(id = "file_custom_interaction_file",
                                           label = "Select the cell-cell location interaction file",
                                           title = "Please select a cell-cell location interaction file",
                                           multiple = FALSE,
                                           viewtype = "list")
            })

            shinyFiles::shinyFileChoose(input,
                                        id = "file_custom_interaction_file",
                                        roots = c(wd = getwd()) )

            shiny::observeEvent(input$file_custom_interaction_file, {

              if(!is.null(input$file_custom_interaction_file)) {
                x <- paste0(getwd(),
                         gsub("/wd", "",
                              paste(unlist(input$file_custom_interaction_file),
                                    collapse = "/")))

                output$text_custom_interaction_file <- shiny::renderUI({
                  shiny::strong({paste("Your selected file is:", x)})
                })

                output$button_read_custom_interaction_file <- shiny::renderUI({
                  shiny::actionButton("button_read_custom_interaction_file",
                                      "Read this interaction file")
                })

                shiny::observeEvent(input$button_read_custom_interaction_file, {
                  x_df <- read.delim(paste0(getwd(),
                                            gsub("/wd", "",
                                                 paste(unlist(input$file_custom_interaction_file),
                                                       collapse = "/"))),
                                     header = TRUE, sep = ",")

                  if(ncol(x_df) == 3) {
                    colnames(x_df) <- c("cell_type_A", "cell_type_B", "value")

                    x_df$join <- paste(x_df$cell_type_A, x_df$cell_type_B, x_df$value,
                                       sep = ",")

                    cell_location_interactions_df(data.frame(parameters = paste0("cell_interaction_",
                                                                                 1:nrow(x_df)),
                                                             value = x_df$join))

                    shiny::showModal(shiny::modalDialog(
                      title = "Success!",
                      "Your custom cell location interactions were saved.",
                      footer = modalButton("OK")
                    ))
                  } else {
                    shiny::showModal(shiny::modalDialog(
                      title = "Error",
                      "Your file doesn’t have the right format.",
                      footer = modalButton("OK")
                    ))
                  }

                })
              }

            })

          } else {
            output$custominteractionSelection <- NULL
            output$button_custominteraction <- NULL
            output$button_select_interaction_file <- NULL
            output$text_custom_interaction_file <- NULL
            output$button_read_custom_interaction_file <- NULL
          }
        })

        output$ask_even_distribution <- shiny::renderUI({
          shiny::sliderInput(inputId = "evendistribution",
                             label = "Adjusts cells to be evenly distributed on a slide. 0 = uneven distribution, 1 = even distribution.",
                             min = 0,
                             max = 1,
                             value = 0,
                             width = "100%")
        })

        shiny::observeEvent(input$evendistribution, {
          cell_even_distribution(input$evendistribution)
        })

      }
    })

    num_simulated_cells <- shiny::reactiveVal()

    cell_overlap_cutoff <- shiny::reactiveVal()


    shiny::observeEvent(simulate_spatial_data(), {

      if(simulate_spatial_data() == TRUE) {

        output$ask_numbercells <- shiny::renderUI({
          shiny::numericInput(inputId = "numbercells",
                              label = "Enter the number of simulated cells",
                              value = 10000,
                              width = "100%")
        })

        shiny::observeEvent(input$numbercells, {
          num_simulated_cells(input$numbercells)
        })

        output$ask_overlapcutoff <- shiny::renderUI({
          shiny::sliderInput(inputId = "overlapcutoff",
                             label = "Select the overlap cutoff. Cells closer than this cutoff will be considered overlapping to each other and all but one will be removed (range: 0-0.1 of the slide length/width)",
                             min = 0,
                             max = 0.1,
                             value = 0,
                             width = "100%")
        })

        shiny::observeEvent(input$overlapcutoff, {
          cell_overlap_cutoff(input$overlapcutoff)
        })

      } else {
        output$ask_numbercells <- NULL
        output$ask_overlapcutoff <- NULL
      }

    })

    # Parameters for expression profiles

    expr_depth_ratio <- shiny::reactiveVal()

    shiny::observeEvent(input$depthratio, {
      expr_depth_ratio(input$depthratio)
    })

    gene_cor <- shiny::reactiveVal()
    copula_input <- shiny::reactiveVal()

    shiny::observeEvent(input$mimiccorrelation, {

      gene_cor(input$mimiccorrelation)

      if(input$mimiccorrelation == TRUE) {
        output$mimiccorrelationSelection <- shiny::renderUI({
          shiny::textInput(inputId = "genecorfile",
                           label = "Provide a path for uploading a pre-estimated gene-gene correlation file:",
                           width = "100%")
        })

        shiny::observeEvent(input$genecorfile, {
          copula_input(input$genecorfile)
        })

      } else {
        output$mimiccorrelationSelection <- NULL
      }
    })

    output$ask_spatialpatterns <- shiny::renderUI({
      shiny::radioButtons(inputId = "spatialpatterns",
                          label = "Would you like to add spatial patterns?",
                          choices = c("No" = FALSE,
                                      "Yes" = TRUE),
                          width = "100%")
    })

    spatialpatterns_df <- shiny::reactiveVal()

    shiny::observeEvent(input$spatialpatterns, {

      if(input$spatialpatterns == TRUE) {
        output$spatialpatternsSelection <- shiny::renderUI({
          shiny::textInput(inputId = "inputspatialpatterns",
                           label = "Specify the spatial patterns in the format
                           <Region>,<Cell type>,<Gene ID (optional)>,<Gene proportion (optional)>,
                           <Mean effect at log(count + 1) scale>,<SD of effect at log(count + 1) scale>.
                           NOTE: if you don't provide a Gene ID, the gene proportion must be provided.
                           If you're providing a Gene ID, the proportion will be automatically calculated
                           and you can set the gene proportion = NULL.
                           Separate the entries by blank space
                           (e.g. '1,cell_type_A,NULL,0.1,0.5,0 1,cell_type_A,gene_A,NULL,0.5,0')

                           Alternatively, click the button to select a file separated by commma where each
                           line contains an entry with the format mentioned above.",
                           width = "100%")
        })

        output$button_spatialpatterns <- shiny::renderUI({
          shiny::actionButton("save_spatialpatterns",
                              "Use these spatial patterns")
        })

        shiny::observeEvent(input$save_spatialpatterns, {

          shiny::observeEvent(input$inputspatialpatterns, {
            x_vector = unlist(stringr::str_split(input$inputspatialpatterns,
                                                 pattern = " "))

            x_df <- data.frame(parameters = character(),
                               value = character())

            for (i in 1:length(x_vector)) {
              split_elements <- unlist(stringr::str_split(x_vector[i],
                                                          pattern = ","))
              if (length(split_elements) == 6) {
                x = data.frame(parameters = c(paste0("spatial_pattern_",i,"_region"),
                                              paste0("spatial_pattern_",i,"_cell_type"),
                                              paste0("spatial_pattern_",i,"_gene_id"),
                                              paste0("spatial_pattern_",i,"_gene_prop"),
                                              paste0("spatial_pattern_",i,"_mean"),
                                              paste0("spatial_pattern_",i,"_sd")),
                               value =  split_elements)

                x_df = rbind(x_df, x)
              } else {
                shiny::showModal(shiny::modalDialog(
                  title = "Error",
                  "Your custom spatial patterns don't have the right format.",
                  footer = modalButton("OK")
                ))
              }

            }

            spatialpatterns_df(x_df)
          })

          shiny::showModal(shiny::modalDialog(
            title = "Success!",
            "Your custom spatial patterns were saved.",
            footer = modalButton("OK")
          ))
        })

        output$button_select_spatialpatterns_file <- shiny::renderUI({
          shinyFiles::shinyFilesButton(id = "file_spatialpatterns",
                                       label = "Select the spatial patterns file",
                                       title = "Please select a spatial patterns file",
                                       multiple = FALSE,
                                       viewtype = "list")
        })

        shinyFiles::shinyFileChoose(input,
                                    id = "file_spatialpatterns",
                                    roots = c(wd = getwd()) )

        shiny::observeEvent(input$file_spatialpatterns, {

          if(!is.null(input$file_spatialpatterns)) {
            x <- paste0(getwd(),
                        gsub("/wd", "",
                             paste(unlist(input$file_spatialpatterns),
                                   collapse = "/")))

            output$text_spatialpatterns_file <- shiny::renderUI({
              shiny::strong({paste("Your selected file is:", x)})
            })

            output$button_read_spatialpatterns_file <- shiny::renderUI({
              shiny::actionButton("button_read_spatialpatterns_file",
                                  "Read this file")
            })

            shiny::observeEvent(input$button_read_spatialpatterns_file, {
              x_df <- read.delim(file = paste0(getwd(),
                                               gsub("/wd", "",
                                                    paste(unlist(input$file_spatialpatterns),
                                                          collapse = "/"))),
                                 sep = ",",
                                 header = TRUE)

              if(ncol(x_df) == 6) {
                colnames(x_df) <- c("region", "annotation", "gene_id",
                                    "proportion", "mean", "sd")

                x_df2 <- data.frame(parameters = character(),
                                    value = character())

                for (i in seq_len(nrow(x_df))) {
                  x_temp <- data.frame(parameters = c(paste0("spatial_pattern_",i,"_region"),
                                                    paste0("spatial_pattern_",i,"_cell_type"),
                                                    paste0("spatial_pattern_",i,"_gene_id"),
                                                    paste0("spatial_pattern_",i,"_gene_prop"),
                                                    paste0("spatial_pattern_",i,"_mean"),
                                                    paste0("spatial_pattern_",i,"_sd")),
                                      value = c(x_df$region[i],
                                                x_df$annotation[i],
                                                x_df$gene_id[i],
                                                x_df$proportion[i],
                                                x_df$mean[i],
                                                x_df$sd[i]))
                  x_df2 <- rbind(x_df2, x_temp)
                }

                spatialpatterns_df(x_df2)

                shiny::showModal(shiny::modalDialog(
                  title = "Success!",
                  "Your custom spatial patterns were saved.",
                  footer = modalButton("OK")
                ))
              } else {
                shiny::showModal(shiny::modalDialog(
                  title = "Error",
                  "Your file doesn't have the right format.",
                  footer = modalButton("OK")
                ))
              }
            })
          }
        })

      } else {
        output$spatialpatternsSelection <- NULL
        output$button_spatialpatterns <- NULL
        output$button_select_spatialpatterns_file <- NULL
        output$text_spatialpatterns_file <- NULL
        output$button_read_spatialpatterns_file <- NULL
      }
    })


    cellcellinteractions_df <- shiny::reactiveVal()

    output$ask_addcellcellinteraction <- shiny::renderUI({
      shiny::radioButtons(inputId = "addcellcellinteraction",
                          label = "Would you like to add cell-cell interactions - expression
                                       associated with cell-cell distance?",
                          choices = c("No" = FALSE,
                                      "Yes" = TRUE),
                          width = "100%")
    })

    shiny::observeEvent(input$addcellcellinteraction, {

      if(input$addcellcellinteraction == TRUE) {

        if(num_regions() == 1 || num_regions() == "NULL") {

          output$addcellcellinteractionSelection <- shiny::renderUI({
            shiny::textInput(inputId = "cellcellinteractions",
                             label = "Specify the <Perturbed cell type>,<Adjacent cell type>,<Interaction
                              distance threshold (default 0.1)>,<Gene ID (optional)>,<Gene proportion (optional)>,<Mean effect at log(count)
                              scale (default = 0.5)>,<SD of effect at log(count) scale (default = 0)>. NOTE: if you don't provide a Gene ID,
                              the gene proportion must be provided. If you're providing a Gene ID, the proportion will be automatically calculated
                              and you can set the gene proportion = NULL.
                              Separate the entries by blank space (e.g. 'cell_A,cell_B,0.1,NULL,0.1,0.5,0 cell_B,cell_C,0.2,gene_A,NULL,0.5,0').
                              NOTE: Adjacent cell type cannot be the same as peturbed cell type.

                              Alternatively, select a file separated by commas, where each line is an
                              entry with the format described above.",
                           width = "100%")
          })
        } else { # num_regions > 1
          output$addcellcellinteractionSelection <- shiny::renderUI({
            shiny::textInput(inputId = "cellcellinteractions",
                             label = "Specify the <Region>,<Perturbed cell type>,<Adjacent cell type>,<Interaction
                              distance threshold (default 0.1)>,<Gene ID (optional)>,<Gene proportion (optional)>,<Mean effect at log(count)
                              scale (default = 0.5)>,<SD of effect at log(count) scale (default = 0)>. NOTE: if you don't provide a Gene ID,
                              the gene proportion must be provided. If you're providing a Gene ID, the proportion will be automatically calculated
                              and you can set the gene proportion = NULL.
                              Separate the entries by blank space (e.g. '1,cell_A,cell_B,0.1,NULL,0.1,0.5,0 2,cell_B,cell_C,0.2,gene_A,NULL,0.5,0').
                              NOTE: Adjacent cell type cannot be the same as peturbed cell type.

                              Alternatively, select a file separated by commas, where each line is an
                              entry with the format described above.",
                           width = "100%")
          })
        }


        output$button_addcellcellinteraction <- shiny::renderUI({
          shiny::actionButton("save_addcellcellinteraction", "Use these cell-cell interactions")
        })

        shiny::observeEvent(input$save_addcellcellinteraction, {

          shiny::observeEvent(input$cellcellinteractions, {
            x_vector = unlist(stringr::str_split(input$cellcellinteractions,
                                                 pattern = " "))

            x_df <- data.frame(parameters = character(),
                               value = character())

            for (i in 1:length(x_vector)) {

              if(length(x_vector[i] == 7)) {# if the region is missed because the n_regions == 1, add NULL
                x_vector[i] = paste0("NULL,", x_vector[i])
              }

              y = unlist(stringr::str_split(x_vector[i],
                                            pattern = ","))
              if(length(y) == 8) {
                x = data.frame(parameters = c(paste0("spatial_int_dist_",i,"_region"),
                                              paste0("spatial_int_dist_",i,"_cell_type_perturbed"),
                                              paste0("spatial_int_dist_",i,"_cell_type_adj"),
                                              paste0("spatial_int_dist_",i,"_dist_cutoff"),
                                              paste0("spatial_int_dist_",i,"_gene_id1"),
                                              paste0("spatial_int_dist_",i,"_gene_prop"),
                                              paste0("spatial_int_dist_",i,"_mean"),
                                              paste0("spatial_int_dist_",i,"_sd")),
                               value =  y)

                x_df = rbind(x_df, x)

              } else {
                shiny::showModal(shiny::modalDialog(
                  title = "Error",
                  "Your custom cell-cell interactions don't have the right format.",
                  footer = modalButton("OK")
                ))
              }

            }

            cellcellinteractions_df(x_df)
          })

          shiny::showModal(shiny::modalDialog(
            title = "Success!",
            "Your custom cell-cell interactions were saved.",
            footer = modalButton("OK")
          ))
        })

        output$button_select_addcellcellinteractionexpression <- shiny::renderUI({
          shinyFiles::shinyFilesButton(id = "file_cellcellinteractionexpression",
                                       label = "Select the cell-cell interactions - expression
                                       associated file",
                                       title = "Please select a cell-cell interactions file",
                                       multiple = FALSE,
                                       viewtype = "list")
        })

        shinyFiles::shinyFileChoose(input,
                                    id = "file_cellcellinteractionexpression",
                                    roots = c(wd = getwd()) )

        shiny::observeEvent(input$file_cellcellinteractionexpression, {
          if(!is.null(input$file_cellcellinteractionexpression)) {
            x <- paste0(getwd(),
                        gsub("/wd", "",
                             paste(unlist(input$file_cellcellinteractionexpression),
                                   collapse = "/")))

            output$text_file_addcellcellinteractionexpression <- shiny::renderUI({
              shiny::strong({paste("Your selected file is:", x)})
            })

            output$button_read_file_addcellcellinteractionexpression <- shiny::renderUI({
              shiny::actionButton("button_read_cellcellinteractionexpression_file",
                                  "Read this file")
            })

            shiny::observeEvent(input$button_read_cellcellinteractionexpression_file, {
              x_df <- read.delim(file = paste0(getwd(),
                                               gsub("/wd", "",
                                                    paste(unlist(input$file_cellcellinteractionexpression),
                                                          collapse = "/"))),
                                 sep = ",",
                                 header = TRUE)

              if(ncol(x_df) == 7) {
                x_df$region <- rep("NULL", nrow(x_df))
                colnames(x_df) <- c("cell_type_perturbed", "cell_type_adj",
                                    "cutoff", "gene_id", "proportion", "mean", "sd",
                                    "region")
              } else if (ncol(x_df) == 8) {
                colnames(x_df) <- c("region", "cell_type_perturbed", "cell_type_adj",
                                    "cutoff", "gene_id", "proportion", "mean", "sd")
              } else {
                shiny::showModal(shiny::modalDialog(
                  title = "Error",
                  "Your custom file doesn't have the right format.",
                  footer = modalButton("OK")
                ))
              }

              if(ncol(x_df) == 8) {
                x1 <- data.frame(parameters = character(),
                                 value = character())

                for (i in 1:nrow(x_df)) {

                  x2 = data.frame(parameters = c(paste0("spatial_int_dist_",i,"_region"),
                                                 paste0("spatial_int_dist_",i,"_cell_type_perturbed"),
                                                 paste0("spatial_int_dist_",i,"_cell_type_adj"),
                                                 paste0("spatial_int_dist_",i,"_dist_cutoff"),
                                                 paste0("spatial_int_dist_",i,"_gene_id1"),
                                                 paste0("spatial_int_dist_",i,"_gene_prop"),
                                                 paste0("spatial_int_dist_",i,"_mean"),
                                                 paste0("spatial_int_dist_",i,"_sd")),
                                  value = c(x_df$region[i],
                                            x_df$cell_type_perturbed[i],
                                            x_df$cell_type_adj[i],
                                            x_df$cutoff[i],
                                            x_df$gene_id[i],
                                            x_df$proportion[i],
                                            x_df$mean[i],
                                            x_df$sd[i]))

                  x1 = rbind(x1, x2)
                }

                cellcellinteractions_df(x1)

                shiny::showModal(shiny::modalDialog(
                  title = "Success!",
                  "Your interactions were saved.",
                  footer = modalButton("OK")
                ))
              }

            })
          }
        })

      } else {
        output$addcellcellinteractionSelection <- NULL
        output$button_addcellcellinteraction <- NULL
        output$button_select_addcellcellinteractionexpression <- NULL
        output$text_file_addcellcellinteractionexpression <- NULL
        output$button_read_file_addcellcellinteractionexpression <- NULL
      }
    })

    cellcellinteractionsexpression_df <- shiny::reactiveVal()

    output$ask_addcellcellinteractionexpression <- shiny::renderUI({
      shiny::radioButtons(inputId = "addcellcellinteractionexpression",
                          label = "Would you like to add cell-cell interactions - expression
                                      associated with expression of neighboring cells??",
                          choices = c("No" = FALSE,
                                      "Yes" = TRUE),
                          width = "100%")
    })

    shiny::observeEvent(input$addcellcellinteractionexpression, {

      if(input$addcellcellinteractionexpression == TRUE) {

        if(num_regions() == 1 || num_regions() == "NULL") {
          output$addcellcellinteractionexpressionSelection2 <- NULL
          output$addcellcellinteractionexpressionSelection1 <- shiny::renderUI({
            shiny::textInput(inputId = "cellcellinteractionsexpression",
                             label = "Specify the <Perturbed cell type>,<Adjacent cell type>,
                              <Interaction distance threshold (default 0.1)>,<Gene ID 1 (optional)>,
                              <Gene ID 2 (optional)>,<Gene proportion (optional)>,
                              <Bi-directional association (TRUE or FALSE, default = TRUE)>,<Mean effect at log(count)
                              scale (default = 0.5)>,<SD of effect at log(count) scale (default = 0)>.
                              NOTE: if you don't provide the Gene IDs, the proportion of genes with
                              this pattern must be provided, and gene pairs will be randomly generated.
                              If you're providing the Gene IDs, the proportion will be automatically
                              calculated and you can set the gene proportion = NULL.
                              Separate the entries by blank space (e.g. 'cell_A,cell_B,0.1,NULL,NULL,0.1,TRUE,0.5,0
                              cell_B,cell_C,0.2,gene_A,gene_B,NULL,TRUE,0.5,0').
                              NOTE: Adjacent cell type cannot be the same as peturbed cell type.

                              Alternatively, select a file separated by commas, where each line is an
                              entry with the format described above.",
                           width = "100%")
          })
        } else { # num_regions > 1
          output$addcellcellinteractionexpressionSelection1 <- NULL
          output$addcellcellinteractionexpressionSelection2 <- shiny::renderUI({
            shiny::textInput(inputId = "cellcellinteractionsexpression",
                             label = "Specify the <Region><Perturbed cell type>,<Adjacent cell type>,
                           <Interaction distance threshold (default 0.1)>,<Gene ID 1 (optional)>,
                           <Gene ID 2 (optional)>,<Gene proportion (optional)>,
                           <Bi-directional association (TRUE or FALSE, default = TRUE)>,<Mean effect at log(count)
                           scale (default = 0.5)>,<SD of effect at log(count) scale (default = 0)>.
                           NOTE: if you don't provide the Gene IDs, the proportion of genes with
                           this pattern must be provided, and gene pairs will be randomly generated.
                           If you're providing the Gene IDs, the proportion will be automatically
                           calculated and you can set the gene proportion = NULL.
                           Separate the entries by blank space (e.g. '1,cell_A,cell_B,0.1,NULL,NULL,0.1,TRUE,0.5,0
                           2,cell_B,cell_C,0.2,gene_A,gene_B,NULL,TRUE,0.5,0').
                           NOTE: Adjacent cell type cannot be the same as peturbed cell type.

                          Alternatively, select a file separated by commas, where each line is an
                          entry with the format described above.",
                           width = "100%")
          })
        }

        output$button_addcellcellinteractionexpression <- shiny::renderUI({
          shiny::actionButton("save_addcellcellinteractionexpression", "Use these cell-cell interactions")
        })

        shiny::observeEvent(input$save_addcellcellinteractionexpression, {

          shiny::observeEvent(input$cellcellinteractionsexpression, {
            x_vector = unlist(stringr::str_split(input$cellcellinteractionsexpression,
                                                 pattern = " "))

            x_df <- data.frame(parameters = character(),
                               value = character())

            for (i in 1:length(x_vector)) {

              if(num_regions() == 1 || num_regions() == "NULL") { # if the region is missed because the n_regions == 1, add NULL
                x_vector[i] = paste0("NULL,", x_vector[i])
              }

              y = unlist(stringr::str_split(x_vector[i],
                                            pattern = ","))

              if(length(y) == 10) {
                x = data.frame(parameters = c(paste0("spatial_int_expr_",i,"_region"),
                                              paste0("spatial_int_expr_",i,"_cell_type_perturbed"),
                                              paste0("spatial_int_expr_",i,"_cell_type_adj"),
                                              paste0("spatial_int_expr_",i,"_dist_cutoff"),
                                              paste0("spatial_int_expr_",i,"_gene_id1"),
                                              paste0("spatial_int_expr_",i,"_gene_id2"),
                                              paste0("spatial_int_expr_",i,"_gene_prop"),
                                              paste0("spatial_int_expr_",i,"_bidirectional"),
                                              paste0("spatial_int_expr_",i,"_mean"),
                                              paste0("spatial_int_expr_",i,"_sd")),
                               value =  y)

                x_df = rbind(x_df, x)
              } else {
                shiny::showModal(shiny::modalDialog(
                  title = "Error",
                  "Your custom cell-cell interactions don't have the right format.",
                  footer = modalButton("OK")
                ))
              }

            }

            cellcellinteractionsexpression_df(x_df)
          })

          shiny::showModal(shiny::modalDialog(
            title = "Success!",
            "Your custom cell-cell interactions were saved.",
            footer = modalButton("OK")
          ))
        })

        output$button_select_cellcellinteractionexpression_neighborhood <- shiny::renderUI({
          shinyFiles::shinyFilesButton(id = "file_cellcellinteraction_neighborhood",
                                       label = "Select the file",
                                       title = "Please select a file",
                                       multiple = FALSE,
                                       viewtype = "list")
        })

        shinyFiles::shinyFileChoose(input,
                                    id = "file_cellcellinteraction_neighborhood",
                                    roots = c(wd = getwd()) )

        shiny::observeEvent(input$file_cellcellinteraction_neighborhood, {
          if(!is.null(input$file_cellcellinteraction_neighborhood)) {
            x <- paste0(getwd(),
                        gsub("/wd", "",
                             paste(unlist(input$file_cellcellinteraction_neighborhood),
                                   collapse = "/")))

            output$text_cellcellinteractionexpression_neighborhood_file <- shiny::renderUI({
              shiny::strong({paste("Your selected file is:", x)})
            })

            output$button_read_cellcellinteractionexpression_neighborhood_file <- shiny::renderUI({
              shiny::actionButton("button_read_cellcellinteraction_neighborhood_file",
                                  "Read this interaction file")
            })

            shiny::observeEvent(input$button_read_cellcellinteraction_neighborhood_file, {
              x_df <- read.delim(paste0(getwd(),
                                        gsub("/wd", "",
                                             paste(unlist(input$file_cellcellinteraction_neighborhood),
                                                   collapse = "/"))),
                                 header = TRUE, sep = ",")

              if(ncol(x_df) == 9) {
                x_df$region <- rep("NULL", nrow(x_df))
                colnames(x_df) <- c("cell_type_perturbed", "cell_type_adj",
                                    "cutoff", "gene_id1", "gene_id2",
                                    "proportion", "bidirectional", "mean", "sd",
                                    "region")
              } else if(ncol(x_df) == 10) {
                colnames(x_df) <- c("region", "cell_type_perturbed", "cell_type_adj",
                                    "cutoff", "gene_id1", "gene_id2",
                                    "proportion", "bidirectional", "mean", "sd")
              } else {
                shiny::showModal(shiny::modalDialog(
                  title = "Error",
                  "Your file doesn't have the right format.",
                  footer = modalButton("OK")
                ))
              }

              if(ncol(x_df) == 10) {
                x1 <- data.frame(parameters = character(),
                                 value = character())

                for (i in seq_len(nrow(x_df))) {

                  x2 <- data.frame(parameters = c(paste0("spatial_int_expr_",i,"_region"),
                                                  paste0("spatial_int_expr_",i,"_cell_type_perturbed"),
                                                  paste0("spatial_int_expr_",i,"_cell_type_adj"),
                                                  paste0("spatial_int_expr_",i,"_dist_cutoff"),
                                                  paste0("spatial_int_expr_",i,"_gene_id1"),
                                                  paste0("spatial_int_expr_",i,"_gene_id2"),
                                                  paste0("spatial_int_expr_",i,"_gene_prop"),
                                                  paste0("spatial_int_expr_",i,"_bidirectional"),
                                                  paste0("spatial_int_expr_",i,"_mean"),
                                                  paste0("spatial_int_expr_",i,"_sd")),
                                   value = c(x_df$region[i],
                                             x_df$cell_type_perturbed[i],
                                             x_df$cell_type_adj[i],
                                             x_df$cutoff[i],
                                             x_df$gene_id1[i],
                                             x_df$gene_id2[i],
                                             x_df$proportion[i],
                                             x_df$bidirectional[i],
                                             x_df$mean[i],
                                             x_df$sd[i])
                  )

                  x1 <- rbind(x1, x2)
                }

                cellcellinteractionsexpression_df(x1)

                shiny::showModal(shiny::modalDialog(
                  title = "Success!",
                  "Your cell-cell interactions were saved.",
                  footer = modalButton("OK")
                ))
              }
            })
          }
        })

      } else {
        output$addcellcellinteractionexpressionSelection1 <- NULL
        output$addcellcellinteractionexpressionSelection2 <- NULL
        output$button_addcellcellinteractionexpression <- NULL
        output$button_select_cellcellinteractionexpression_neighborhood <- NULL
        output$text_cellcellinteractionexpression_neighborhood_file <- NULL
        output$button_read_cellcellinteractionexpression_neighborhood_file <- NULL
      }
    })


    num_spots <- shiny::reactiveVal()

    shiny::observeEvent(input$multicell, {

      if(input$multicell == TRUE) {

        output$multicellSelection <- renderUI({
          shiny::numericInput(inputId = "nspots",
                              label = "Specify the number of spots",
                              value = 400,
                              width = "100%")
        })

        shiny::observeEvent(input$nspots, {
          num_spots(input$nspots)
        })
      } else {
        output$multicellSelection <- NULL
        num_spots("NULL")}
    })

    num_simulated_datasets <- shiny::reactiveVal()
    parent_simulation_seed <- shiny::reactiveVal()

    shiny::observeEvent(input$ndatasets, {
      num_simulated_datasets(input$ndatasets)
    })

    output$ask_single_seed <- renderUI({
      shiny::numericInput(inputId = "seed",
                          label = "An umbrella seed is required for reproducible simulation. Specify an umbrella seed.",
                          value = 1234,
                          width = "100%")
    })

    shiny::observeEvent(input$seed, {
      parent_simulation_seed(input$seed)
    })

    path_to_output_dir <- shiny::reactiveVal()

    shiny::observeEvent(input$outdir, {
      path_to_output_dir(fs::path(getwd(),input$outdir))
    })

    output_name <- shiny::reactiveVal()

    shiny::observeEvent(input$outname, {
      output_name(input$outname)
    })

    parameter_file_name <- shiny::reactiveVal()

    shiny::observeEvent(input$parameter_filename, {
      parameter_file_name(input$parameter_filename)
    })



    # export parameters file
    output$downloadparams <- shiny::downloadHandler(
      filename = function() {parameter_file_name()},
      content = function(con) {

        param_df <- data.frame(parameters = c("path_to_input_dir",
                                              "expression_data_file",
                                              "cell_feature_data_file",
                                              "expression_data_file_type",
                                              "expression_data_cell_types",
                                              "simulate_spatial_data"
        ),
        value = c(path_to_input_dir(),
                  expression_data_file(),
                  cell_feature_data_file(),
                  expression_data_file_type(),
                  expression_data_cell_types(),
                  simulate_spatial_data()
        )
        )

        if(simulate_spatial_data() == TRUE) {
          df = data.frame(parameters = c("num_simulated_cells",
                                         "cell_overlap_cutoff"),
                          value = c(num_simulated_cells(),
                                    cell_overlap_cutoff()
                          )
          )

          param_df = rbind(param_df, df)

          param_df = rbind(param_df, c("custom_cell_location_interactions",
                                       custom_cell_location_interactions() ))

          if(custom_cell_location_interactions() == TRUE) {
            param_df = rbind(param_df, cell_location_interactions_df())
          }

        }

        if(ncol_feature_data() > 1 & simulate_spatial_data() == TRUE) {
          param_df = rbind(param_df, c("window_method", window_method()))
        }

        if(!is.null(num_regions())) {
          param_df = rbind(param_df, c("num_regions", num_regions()))
        }

        param_df = rbind(param_df, c("region_specific_model",
                                     region_specific_model()))

        if(ncol_feature_data() == 1) {

          param_df = rbind(param_df, c("custom_cell_type_proportions",
                                       custom_props()))

          param_df = rbind(param_df, cell_type_proportions())

          param_df = rbind(param_df,c("cell_even_distribution",
                                      cell_even_distribution() ))
        }

        # simulation parameters for expression profiles
        param_df = rbind(param_df, c("expr_depth_ratio",
                                     expr_depth_ratio() ))

        param_df = rbind(param_df, c("gene_cor",
                                     gene_cor() ))

        if(gene_cor() == TRUE) {
          param_df = rbind(param_df, c("copula_input",
                                       copula_input() ))
        }

        if(!is.null(spatialpatterns_df())) {
          param_df = rbind(param_df, spatialpatterns_df())
        }

        if(!is.null(cellcellinteractions_df())) {
          param_df = rbind(param_df, cellcellinteractions_df())
        }

        if(!is.null(cellcellinteractionsexpression_df())) {
          param_df = rbind(param_df, cellcellinteractionsexpression_df())
        }

        param_df = rbind(param_df, c("num_spots",
                                     num_spots()))

        param_df = rbind(param_df, c("num_simulated_datasets",
                                     num_simulated_datasets()))

        if(!is.null(parent_simulation_seed())) {
          param_df = rbind(param_df, c("parent_simulation_seed",
                                       parent_simulation_seed()))
        }

        param_df = rbind(param_df, c("path_to_output_dir",
                                     path_to_output_dir()))

        param_df = rbind(param_df, c("output_name",
                                     output_name()))

        param_df = rbind(param_df, c("parameter_file",
                                     parameter_file_name()))



        write.table(param_df,
                    file = con,
                    quote = FALSE,
                    row.names = FALSE,
                    sep = "\t")
      }
    ) # ends export parameters file

    # Run simulation

    parameter_file <- shiny::reactiveVal()
    shiny::observeEvent(input$inputparamfile, {
      parameter_file(input$inputparamfile)
    })

    create_Giotto <- shiny::reactiveVal()
    giotto_folder <- shiny::reactiveVal()

    shiny::observeEvent(input$createGiotto, {
      create_Giotto(input$createGiotto)

      if(input$createGiotto) {
        output$createGiottoSelection <- renderUI({
          shiny::textInput(inputId = "giottofolder",
                           label = "Specify the folder name to export the Giotto object",
                           value = "giottobject",
                           width = "100%")
        })

        shiny::observeEvent(input$giottofolder, {
          giotto_folder(input$giottofolder)
        })
      }


      shiny::observeEvent(input$runsimulation, {

        ParaSimulation(input = parameter_file())

        if(create_Giotto()) {

          x_param = ParaDigest(parameter_file())

          x_expression = read.delim(fs::path(x_param$path_to_output_dir,
                                     paste0(x_param$output_name,
                                           "_count_1.tsv")),
                                    row.names = 1)

          x_meta = read.delim(fs::path(x_param$path_to_output_dir,
                              paste0(x_param$output_name,
                                     "_meta_1.tsv"))
                              )

          x_spatlocs = x_meta[,c("Cell", "x.loc", "y.loc")]
          colnames(x_spatlocs) = c("cell_ID", "sdimx", "sdimy")

          x_meta = x_meta[,c("Cell", "annotation", "region")]
          colnames(x_meta)[1] = "cell_ID"

          x = Giotto::createGiottoObject(expression = x_expression,
                                         spatial_locs = x_spatlocs)

          x = Giotto::addCellMetadata(x,
                                      new_metadata = x_meta[-1])

          Giotto::saveGiotto(x, foldername = giotto_folder(), overwrite = TRUE)
        }

        shiny::stopApp()
      })

    }) # ends run simulation



    shiny::observeEvent(input$done, {
      shiny::stopApp()
    })
  }

  shiny::runGadget(ui, server)
}


