# This script contains the server part for the Phosphoproteomics Shiny app. 
# The server script acts based on the inputs from the UI script. 

# Increasing the limit of upload file size from 5 MB to 50 MB
options(shiny.maxRequestSize=500*1024^2)

#temporary directory for unzip raw data folder
outDir <- "rawFolder"

shinyServer(function(input, output, session) {
  
  # a reactive val to store the multiassayexperiment object
  mae <- reactiveVal()
  
  # a reactive variable to stored the processed assay
  processedData <- reactiveVal()

  # a reactive value to store the saved list
  saveList <- reactiveValues(file = list.files("save/"))
  
  # a box showing saved results
  output$saveListBox <- renderUI({
    selectInput("seleFile", label = NULL, saveList$file, size = 5, selectize = FALSE)
  })
  
  # save calculated results
  observeEvent(input$save, {
    if (!is.null(mae())) {
      saveObj <- mae()
      fileName <- paste0(input$text, "_", format(Sys.Date(), "%Y%m%d"), ".Rds")
      saveRDS(saveObj, file = paste0("save/", fileName))
      saveList$file <- unique(c(saveList$file, fileName))
    }
  })
  
  # remove saved file
  observeEvent(input$remove, {
    fileName <- input$seleFile
    file.remove(paste0("save/", fileName))
    saveList$file <- saveList$file[saveList$file != fileName]
  })
  
  # load saved file
  observeEvent(input$load, {
    mae(readRDS(paste0("save/", input$seleFile)))
  })
  
  # zip file
  observeEvent(input$uploadZip, {
    # removing the already existing directory before unzipping
    unlink(outDir, recursive = TRUE)
    dir.create(outDir, showWarnings = FALSE)
    unzip(input$uploadZip$datapath, exdir = outDir, junkpaths = TRUE)
    
    # error check
    fileTable <- NULL
    tryCatch({
      fileTable <- as.data.frame(read.delim(file.path(outDir, "fileTable.txt")))
    },
    error = function(e) {
      showModal(modalDialog(
        title = "Processing zip files failed...",
        "fileTable.txt file not found.",
        easyClose = TRUE,
        footer = NULL
      ))})
    if(!is.null(fileTable))
    {
      tryCatch({
        stopifnot(c("sampleType", "id", "fileName") %in% colnames(fileTable))
        # rendering column annotation option for column annotations
        output$colAnnoBoxPreprocess <- renderUI({
          # excluding type since it's already represented by two assays
          selectInput("colAnnoPreprocess", "Select additional column annotations:",
                      colnames(fileTable)[colnames(fileTable) != "type"],
                      selected = colnames(fileTable)[colnames(fileTable) != "type"], multiple = TRUE)
        })
      },
      error = function(e) {
        showModal(modalDialog(
          title = "Missing columns in fileTable.txt file...",
          "Please make sure fileTable.txt file contains columns with column names: fileName, sampleType and id.",
          easyClose = TRUE,
          footer = NULL
        ))})
    }
  })
  
  # process raw files
  observeEvent(input$processUpload, {
    withProgress(message = 'Processing files', {
      fileTable <- as.data.frame(read.delim(file.path(outDir, "fileTable.txt")))
      #fileTable$fileName <- sub("^", "./rawFolder/", fileTable$fileName)
      fileTable$fileName <- file.path(outDir, fileTable$fileName) #outDir is a variable
      
      tryCatch({
        # for data from Spectronaut
        if (input$tool == "Spectronaut") {
          testData <- SmartPhos::readExperimentDIA(fileTable, annotation_col = input$colAnnoPreprocess)
        }
        # for data from MaxQuant
        else if (input$tool == "MaxQuant") {
          testData <- SmartPhos::readExperiment(fileTable, annotation_col = input$colAnnoPreprocess)
        }
        mae(testData)
      },
      error = function(e) {
        showModal(modalDialog(
          title = "Processing of the uploaded data failed...",
          "Please make sure the correct option for the tool (Spectronaut or MaxQuant) is selected.",
          easyClose = TRUE,
          footer = NULL
        ))})
    })
  })
  
  # read the uploaded object
  observeEvent(input$uploadObject, {
    file <- input$uploadObject
    ext <- tools::file_ext(file$datapath)
    req(file)
    # making sure an rds object is being uploaded
    validate(need(ext %in% c("rds", "RDS", "Rds"), "Please upload an rds object"))
    mae(readRDS(file$datapath))
  })

  # rendering select options in the UI based on the uploaded object
  output$option <- renderUI({
    if (!is.null(mae())) {
      selectInput("assay", "Select one assay", names(experiments(mae())))
    }
  })
  
  # rendering options for the phospho-enriched or non-enriched sample type 
  output$seleAssayBox <- renderUI({
    if (!is.null(mae())) {
      if (input$assay == "Proteome") {
        radioButtons("getPP", "Select the sample type",
                     c("Phospho-enriched" = TRUE, "Non-enriched" = FALSE),
                     selected = FALSE, inline = TRUE)
      }
      else if (input$assay == "Phosphoproteome") {
        radioButtons("getFP", "Select the sample type",
                     c("Phospho-enriched" = FALSE, "Non-enriched" = TRUE),
                     selected = FALSE, inline = TRUE)
      }
    }
  })
  
  # loaded object for boxplot and table
  loadedData <- reactive({
    se <- mae()[[input$assay]]
    colData(se) <- colData(mae()[, colnames(se)])
    if (input$assay == "Phosphoproteome") {
      if (input$getFP) {
        ppe <- se[,se$sampleType == "FullProteome"]
        colData(ppe) <- colData(se)[colnames(ppe),]
      }
      else {
        ppe <- se[,se$sampleType == "Phospho"]
        colData(ppe) <- colData(se)[colnames(ppe),]
      }
      ppe
    }
    else if (input$assay == "Proteome") {
      if (input$getPP) {
        fpe <- se[,se$sampleType == "Phospho"]
        colData(fpe) <- colData(se)[colnames(fpe),]
      }
      else {
        fpe <- se[,se$sampleType == "FullProteome"]
        colData(fpe) <- colData(se)[colnames(fpe),]
      }
      fpe
    }
  })
  
  
  #----------------------------------------------------------------------------------------------------------
  # launching MatrixQCvis from the phosphoproteomics app
  # observeEvent(input$launch_app, {
  #   mae <- readRDS(input$upload$datapath)
  #   se <- mae[[input$assay]]
  #   colData(se) <- colData(mae)
  #   job_env <-  new.env()
  #   job_env$se <- se
  #   rstudioapi::jobRunScript(path = "script.R", importEnv = TRUE)
  # })
  observeEvent(input$launch_app, {
    showModal(modalDialog(
      title = "Launching MatrixQCvis not possible...",
      "This functionality is yet to be implemented in the Shiny App.",
      easyClose = TRUE,
      footer = NULL
    ))
  })
  #-----------------------------------------------------------------------------------------------------------
  
  # rendering options for the phospho-enriched or non-enriched sample type 
  output$seleNormCorrect <- renderUI({
    if (!is.null(mae())) {
      if (input$assay == "Phosphoproteome" & input$getFP == FALSE) {
        checkboxInput("ifNormCorrect","Perform normalization correction (may take long time, no further normalization needed)", value = FALSE)
      }
    }
  })
  
  output$ifAlreadyNormBox <- renderUI({
      if(input$ifNormCorrect) {
          updateRadioButtons(session = getDefaultReactiveDomain(), inputId = "normalize", selected = FALSE)
          radioButtons("ifAlreadyNorm", "Have the data already been normalized by Spectronaut/MaxQuant ?", c("Yes","No"), selected = "No", inline = TRUE)
      }
  })
  
  
  observeEvent(input$processSelection, {
    withProgress(message = 'Processing files', {
      # normalization correction if selected
      # first step, whether to perform phospho normalization correction, if phospho data is selected.    
      if (input$assay == "Phosphoproteome") {
        if (input$ifNormCorrect) {
          maeData <- runPhosphoAdjustment(mae(), 
                                          normalization = ifelse(input$ifAlreadyNorm == "Yes", FALSE, TRUE), #depends on whether the data were already normalized
                                          minOverlap = 3, #at least three overlapped feature between sample pair
                                          completeness = 0.5 #use feature present in at least 50% of the samples
                                          )
          assays(maeData[["Phosphoproteome"]])[["Intensity"]] <- assays(maeData[["Phosphoproteome"]])[["Intensity_adjusted"]]
          assays(maeData[["Phosphoproteome"]])[["Intensity_adjusted"]] <- NULL
          # if (is.na(maeData$adjustFactorPP)) {
          #   showModal(modalDialog(
          #     title = "One sample removed because of low quality...",
          #     easyClose = TRUE,
          #     footer = NULL
          #   ))
          # }
          maeData <- maeData[, !is.na(maeData$adjustFactorPP)]
        }
        else maeData <- mae()
      }
        
      else maeData <- mae()
      # summarizedAssayExperiment object of the selected assay
      se <- maeData[[input$assay]]
      colData(se) <- colData(maeData[, colnames(se)])
      if (input$assay == "Proteome") {
        # function from the utils.R script for the preprocessing of the data
        fp <- preprocessProteome(se, filterList = NULL,
                                 transform = input$transform,
                                 normalize = input$normalize,
                                 getPP = input$getPP,
                                 missCut = input$missFilter,
                                 removeOutlier = strsplit(input$outliers, ",\\s*")[[1]],
                                 impute = input$impute,
                                 scaleFactorTab = NULL)
        processedData(fp)
      }
      else {
        pp <- preprocessPhos(se, filterList = NULL,
                             transform = input$transform,
                             normalize = ifelse(!is.null(input$ifNormCorrect) & input$ifNormCorrect ==FALSE, input$normalize, FALSE), #if normalization correction has been performed, normalization should not be performed again
                             getFP = input$getFP,
                             missCut = input$missFilter,
                             removeOutlier = strsplit(input$outliers, ",\\s*")[[1]],
                             impute = input$impute,
                             scaleFactorTab = NULL)
        processedData(pp)
      }
    })
  })
  
  # text output to show the number of samples and features
  output$dataInfo <- renderUI({
    if (!is.null(processedData())) {
      HTML(sprintf("<b>Number of samples: %s<br/>Number of features: %s<b><br/>",
                   ncol(processedData()), nrow(processedData())))
    }
    else {
      HTML(sprintf("<b>Number of samples: %s<br/>Number of features: %s<b><br/>",
                   ncol(loadedData()), nrow(loadedData())))
    }
  })
  
  # data table
  output$metaData <- DT::renderDataTable({
    if (!is.null(processedData())) {
      colDataTable <- mutate_if(data.frame(colData(processedData())), is.character, as.factor)
      datatable(colDataTable, filter = "top", rownames = FALSE,
                selection = "none", style = "bootstrap")
    }
    else {
      colDataTable <- mutate_if(data.frame(colData(loadedData())), is.character, as.factor)
      datatable(colDataTable, filter = "top", rownames = FALSE,
                selection = "none", style = "bootstrap")
    }
    
  })

  # Plot missing value
  output$missingPlot <- renderPlot({
    countMat <- assay(loadedData())
    plotTab <- tibble(sample = loadedData()$sample, 
                      perNA = colSums(is.na(countMat))/nrow(countMat))
    ggplot(plotTab, aes(x = sample, y = 1-perNA)) +
      geom_bar(stat = "identity") +
      ggtitle("Percentage of sample completeness") +
      ylab("completeness") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0),
            plot.title = element_text(hjust = 0.5, face = "bold"))
  })
  
  # to toggle between show and hide functionality
  toggle("missingPlot")
  observeEvent(input$missingV, {
   toggle("missingPlot") 
  })
  
  output$colorBoxUI <- renderUI({
    selectInput("colorBox", "Select color by:",
                c("none", colnames(colData(mae()))),
                selected = "none")
  })
  
  # plot log ratio
  output$boxPlotLogRatio <- renderPlot({
    plotLogRatio(mae(), normalization = FALSE)
  })
  
  # Plot boxplot
  output$boxPlot <- renderPlot({
    if (!is.null(processedData())) {
      countMat <- assay(processedData())
      countTab <- countMat %>% as_tibble(rownames = "id") %>% 
        pivot_longer(-id) %>%
        filter(!is.na(value))
      meta <- as.data.frame(colData(processedData()))
    }
    else {
      countMat <- assay(loadedData())
      countTab <- countMat %>% as_tibble(rownames = "id") %>% 
        pivot_longer(-id) %>%
        filter(!is.na(value))
      meta <- as.data.frame(colData(loadedData()))
    }
    countTabmeta <- left_join(countTab, meta, by = c('name' = 'sample'))
    g <- ggplot(countTabmeta, aes(x = name, y = value)) +
      ggtitle("Boxplot of intensities") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0),
            plot.title = element_text(hjust = 0.5, face = "bold"))
    if (input$colorBox == "none"){
      g <- g + geom_boxplot()
    }
    else {
      g <- g + geom_boxplot(aes_string(fill = input$colorBox))
    }
    g
  })
  
  ############################################## PCA ########################################################
  
  runPCA <- observeEvent(input$RunPCA, {
    
    if ("imputed" %in% assayNames(processedData())) {
      output$errMsgPCA <- renderText("")
      withProgress(message = "Running principal component analysis, please wait...", {
        pca <- stats::prcomp(t(assays(processedData())[["imputed"]]))
        # variance explained
        varExplained <- pca$sdev^2/sum(pca$sdev^2)
        pcaDf <- as.data.frame(pca[["x"]])
        # download PCA values as tsv file
        output$downloadPCATable <- downloadHandler(
          filename = function() { paste('PCAvalues', '.tsv', sep='') },
          content = function(file) {
            write.table(pcaDf, file = file, quote = FALSE, sep = '\t', col.names = NA)
        })
        # rendering select option for x-axis of PCA plot
        output$xaxisPCAui <- renderUI({
          selectInput("xaxisPCA", "Select the x-axis PC", colnames(pcaDf),
                      selected = "PC1")
        })
        # rendering select option for y-axis of PCA plot
        output$yaxisPCAui <- renderUI({
          selectInput("yaxisPCA", "Select the y-axis PC", colnames(pcaDf),
                      selected = "PC2")
        })
        # rendering color option for PCA plot
        output$colorPCAui <- renderUI({
          selectInput("colorPCA", "Select color by:",
                      c("none",colnames(colData(processedData()))),
                      selected = "none")
        })
        # rendering shape option for PCA plot
        output$shapePCAui <- renderUI({
          selectInput("shapePCA", "Select shape by:",
                      c("none",colnames(colData(processedData()))),
                      selected = "none")
        })

        # adding colData information to pca results
        meta <- as.data.frame(colData(processedData()))
        pcaMeta <- left_join(rownames_to_column(pcaDf),
                             meta, by = c("rowname" = "sample"))
        
        plotpc <- reactive({
          g <- ggplot(pcaMeta, aes(x = !!sym(input$xaxisPCA), y = !!sym(input$yaxisPCA),
                                   text = paste("sample:", meta$sample))) +
            theme_bw() +
            theme(legend.position="top") +
            labs(x=paste0(input$xaxisPCA,": ",
                          round(varExplained[as.numeric(strsplit(input$xaxisPCA, "PC")[[1]][2])]*100, 1), "%"),
                 y=paste0(input$yaxisPCA,": ",
                          round(varExplained[as.numeric(strsplit(input$yaxisPCA, "PC")[[1]][2])]*100, 1), "%")) +
            scale_shape(solid = FALSE)
          
          if (input$colorPCA == "none" & input$shapePCA == "none") {
            g <- g + geom_point(size = 2)
          }
          else if (input$colorPCA == "none") {
            g <- g + geom_point(aes_string(shape = input$shapePCA), size = 2)
          }
          else if (input$shapePCA == "none") {
            g <- g + geom_point(aes_string(color = input$colorPCA), size = 2)
          }
          else {
            g <- g + geom_point(aes_string(color = input$colorPCA,
                                    shape = input$shapePCA),
                                size = 2)
          }
          
        })
        output$pcplot <- renderPlotly({
          ggplotly(plotpc())
        })
      })} 
    
    else {
      showModal(modalDialog(
        title = "PCA not possible...",
        "Please make sure imputation is not selected none in the preprocessing tab.",
        easyClose = TRUE,
        footer = NULL
      ))
    }
    
    # download the PCA plot as PDF file
    output$downPCA <- downloadHandler(
      filename = function() { paste0("pca", '.pdf', sep = '') },
      content = function(file) {
        ggsave(file, plot = plotpc(), 
               device = "pdf",
               width = input$figWidthPCA,
               height = input$figHeightPCA,
               limitsize = FALSE)
      }
    )
  })
  
  ######################################################## heat map ##########################################################################
  
  # rendering column annotation option for heatmap
  output$colAnnoBoxHM <- renderUI({
    if (!is.null(processedData())) {
      selectInput("colAnnoHM", "Select additional column feature:",
                  colnames(colData(processedData())),
                  selected = NULL, multiple = TRUE)
    }
  })
  
  # render options for top variants option
  output$topGenes <- renderUI({
    numericInput("numGenes","Number of genes", value = 100)
  })
  
  output$colCluster <- renderUI({
    numericInput("numClustCol", "Number of column clusters", value = 1)
  })
  
  output$rowCluster <- renderUI({
    numericInput("numClustRow", "Number of row clusters", value = 1)
  })
  
  # for top variants based on standard deviation
  orderID <- reactive({
    exprMat <- assays(processedData())[["imputed"]]
    sds <- apply(exprMat, 1, sd)
    return(names(sort(sds, decreasing = TRUE)))
  })
  
  plotMap <- eventReactive(input$doPlot, {
    
    if (input$chooseType == "Top variant") {
      # plot top variant genes
      setName <- sprintf("Top %s most variant genes", input$numGenes)
      geneIDs <- orderID()[seq(1, as.integer(input$numGenes))]
      exprMat <- assays(processedData())[["imputed"]][geneIDs,]
      geneSymbol <- rowData(processedData()[match(geneIDs, rownames(processedData())),])$Gene
    } 
    else if (input$chooseType == "Differentially expressed") {
      # plot differentially expressed genes
      if(!is.null(filterDE())) {
        setName <- "Differentially expressed genes"
        # prepare the data matrix
        geneIDs <- arrange(filterDE(), stat)$ID
        exprMat <- assays(processedDataSub())[["imputed"]][geneIDs,] 
        geneSymbol <- filterDE()[match(geneIDs, filterDE()$ID),]$Gene
        }
      else {
        showModal(modalDialog(
          title = "Plotting heatmap not possible...",
          "Make sure you have performed Differential expression analysis before.",
          easyClose = TRUE,
          footer = NULL
        ))
      }
    }
    else if (input$chooseType == "Selected time series cluster") {
      # plot genes from the selected cluster
      if(!is.null(selectedCluster())) {
        setName <- input$seleCluster
        # prepare the data matrix
        geneIDs <- unique(selectedCluster()$ID)
        exprMat <- assays(processedDataSub())[["imputed"]][geneIDs,] 
        geneSymbol <- selectedCluster()[match(geneIDs, selectedCluster()$ID),]$Gene
      }
      else {
        showModal(modalDialog(
          title = "Plotting heatmap not possible...",
          "Make sure you have performed Time series clustering before.",
          easyClose = TRUE,
          footer = NULL
        ))
      }
    }
    
    # prepare column annotations
    cd <- as.data.frame(colData(processedData()))
    annCol <- cd[row.names(cd) %in% colnames(exprMat),][c(input$colAnnoHM)]
    #annCol <- cd[colnames(exprMat), input$colAnno, drop=FALSE]
    row.names(annCol) <- colnames(exprMat)
    
    # prepare color scale
    color <- colorRampPalette(c("navy", "white", "firebrick"))(100)
    
    # manual row normalization
    exprMat <- t(scale(t(exprMat)))
    exprMat[exprMat > 4] <- 4
    exprMat[exprMat < -4] <- -4
    
    # plot heatmap
    if (input$chooseType == "Top variant") {
      if (is.null(input$colAnnoHM)) {
        p <- pheatmap(exprMat, color = color,
                      labels_row = geneSymbol,
                      treeheight_row = 0, treeheight_col = 0,
                      main = setName,
                      cutree_cols = input$numClustCol,
                      cutree_rows = input$numClustRow,
                      silient = TRUE)
      }
      else {
        p <- pheatmap(exprMat, color = color,
                      labels_row = geneSymbol,
                      treeheight_row = 0, treeheight_col = 0,
                      main = setName,
                      cutree_cols = input$numClustCol,
                      cutree_rows = input$numClustRow,
                      annotation_col = annCol,
                      silient = TRUE)
      }
    }
    else {
      # first sort the columns by column names
      exprMat <- exprMat[, sort(colnames(exprMat))]
      if (is.null(input$colAnnoHM)) {
        p <- pheatmap(exprMat, color = color,
                      labels_row = geneSymbol,
                      treeheight_row = 0, treeheight_col = 0,
                      main = setName,
                      cluster_rows = FALSE, cluster_cols = FALSE,
                      cutree_cols = input$numClustCol,
                      cutree_rows = input$numClustRow,
                      silient = TRUE)
      }
      else {
        p <- pheatmap(exprMat, color = color,
                      labels_row = geneSymbol,
                      treeheight_row = 0, treeheight_col = 0,
                      main = setName,
                      cluster_rows = FALSE, cluster_cols = FALSE,
                      cutree_cols = input$numClustCol,
                      cutree_rows = input$numClustRow,
                      annotation_col = annCol,
                      silient = TRUE)
      }
    }
    p
  })
  
  # plot the heatmap
  output$plotHM <- renderPlot({
    if (!is.null(plotMap())) {
      output$errMsg1 <- renderText("")
      withProgress(message = "Plotting heatmap, please wait...", value = NULL, {
        grid.draw(plotMap())
      })
    } else {
      output$errMsg1 <- renderText("Please perform differential expression analysis first or load a previous result!")
    }
  })
  
  # download the Heatmap plot as PDF file
  output$downHM <- downloadHandler(
    filename = function() { paste0("heatmap", '.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = plotMap(),
             device = "pdf",
             width = input$figWidthHM,
             height = input$figHeightHM,
             limitsize = FALSE)
  })
  
  ###################################################### differential expression ####################################################
  
  
  # select sample by sample name
  output$ID1Box <- renderUI({
    selectInput("seleID1", "Select sample ID(s) for reference group",
                processedData()$sample, multiple = TRUE)
  })
  
  output$ID2Box <- renderUI({
    selectInput("seleID2", "Select sample ID(s) for target group",
                processedData()$sample, multiple = TRUE)
  })
  
  # the selection box to select treatment for the reference group
  output$treat1Box <- renderUI({
    allTreat <- unique(processedData()$treatment)
    selectInput("seleTreat1", "Select treatment(s) for the reference group", allTreat, multiple = TRUE)
  })
  
  # select time point for reference
  output$time1Box <- renderUI({
    if (!is.null(processedData()$timepoint)) {
      if (!is.null(input$seleTreat1))
        allTime <- unique(processedData()[,processedData()$treatment %in% input$seleTreat1]$timepoint) else
          allTime <- unique(processedData()$timepoint)
        selectInput("seleTime1", "Select time point(s) for the reference group", allTime, multiple = TRUE)
    }
  })
  
  # the selection box to select treatment for the target groups
  output$treat2Box <- renderUI({
    allTreat <- unique(processedData()$treatment)
    selectInput("seleTreat2", "Select treatment(s) for the target group", allTreat, multiple = TRUE)
  })
  
  # the select the comparison time point
  output$time2Box <- renderUI({
    if (!is.null(processedData()$timepoint)) {
      if (!is.null(input$seleTreat2))
        allTime <- unique(processedData()[,processedData()$treatment %in% input$seleTreat2]$timepoint) else
          allTime <- unique(processedData()$timepoint)
        selectInput("seleTime2", "Select time point(s) for the target group", allTime, multiple = TRUE)
    }
  })
  
  # List of sample ID for reference group
  listIDforDE1 <- reactive(
    if (input$seleID) {
      selectedID <- input$seleID1
      selectedID
    } else {
      if (!is.null(processedData()$timepoint))
        selectedID <- processedData()[,processedData()$treatment %in% input$seleTreat1 & processedData()$timepoint %in% input$seleTime1]$sample else
          selectedID <-processedData()[,processedData()$treatment %in% input$seleTreat1]$sample
        selectedID
    }
  )
  
  # List of sample ID for target group
  listIDforDE2 <- reactive(
    if (input$seleID) {
      selectedID <- input$seleID2
      selectedID
    } else {
      if (!is.null(processedData()$timepoint))
        selectedID <- processedData()[,processedData()$treatment %in% input$seleTreat2 & processedData()$timepoint %in% input$seleTime2]$sample else
          selectedID <- processedData()[,processedData()$treatment %in% input$seleTreat2]$sample
        selectedID
    }
  )
  
  output$infoDE <- renderUI({
    if ((!is.null(listIDforDE1())) & (!is.null(listIDforDE2()))) {
      if (length(base::intersect(listIDforDE1(), listIDforDE2())) > 0) {
        HTML(sprintf("<b>Reference group has %s samples<br/>Target group has %s samples<b><br/>WARNING: some sample(s) are present in both groups<b><br/> ",
                   length(listIDforDE1()), length(listIDforDE2())))} else {
        HTML(sprintf("<b>Reference group has %s samples<br/>Target group has %s samples<b><br/> ",
                                  length(listIDforDE1()), length(listIDforDE2())))
        }
    }
  })
  
  output$seleMethodBox <- renderUI({
    allowChoice <- c("limma", "ProDA")
    radioButtons("deMethod", "Select DE method", allowChoice, inline = TRUE)
  })
  
  # a reactive object for subsetting RNAseq dataset
  processedDataSub <- reactive({
    processedData.sub <- processedData()[,processedData()$sample %in% c(listIDforDE1(),listIDforDE2())]
    processedData.sub$comparison <-  ifelse(processedData.sub$sample %in% listIDforDE1(), "reference", "target")
    processedData.sub$comparison <- factor(processedData.sub$comparison, levels = c("reference", "target"))
    processedData.sub
  })
  
  tableDE <- reactiveVal()
  # a reactive value to monitor if plot histogram or boxplot
  ifHistogram <- reactiveValues(value = TRUE)
  # a reactive value to save the clicked row or volcano plot point value
  lastClicked <- reactiveValues()
  # a reactive value to sync the data table with the volcano plot selection
  colorRows <- reactiveValues(
    row_priority = c(),
    row_color = c()
  )
  # a reactive value to highlight the point with row selection
  ptHiglight <- reactiveValues(
    log2FC = NULL,
    pValue = NULL
  )
  
  # reactive event for calculating differential expression
  observeEvent(input$runDE, {
    withProgress(message = "Running DE analysis, please wait...", value = NULL, {
      seqMat <- processedDataSub()
      exprMat <- assays(seqMat)[["Intensity"]]
      colData <- data.frame(colData(seqMat))
      tryCatch({
        # design matrix
        if(is.null(processedDataSub()$subjectID)) {
          design <- model.matrix(~ comparison, data = colData)
        } else {
          design <- model.matrix(~ subjectID + comparison, data = colData)
        }
        resNames <- colnames(design)
        meta <- as.data.frame(elementMetadata(processedDataSub()))
        if (input$deMethod == "limma") {
          fit <- limma::lmFit(exprMat, design = design)
          fit2 <- eBayes(fit)
          resDE <- topTable(fit2, number = Inf, coef=resNames[length(resNames)])
          # get result
          resDE <- merge(resDE, meta, by=0, all=TRUE)
          resDE <- as_tibble(resDE) %>%
            dplyr::rename(log2FC = logFC, stat = t,
                          pvalue = P.Value, padj = adj.P.Val, ID = Row.names) %>%
            select(-c(B, AveExpr)) %>%
            filter(!is.na(padj)) %>%
            arrange(pvalue)
        }
        else if (input$deMethod == "ProDA") {
          fit <- proDA::proDA(exprMat, design = design)
          resDE <- proDA::test_diff(fit, contrast = resNames[length(resNames)])
          # get result
          rownames(resDE) <- resDE[,1]
          resDE[,1] <- NULL
          resDE <- merge(resDE, meta, by=0, all=TRUE)
          resDE <- as_tibble(resDE) %>%
            dplyr::rename(log2FC = diff, stat = t_statistic,
                          pvalue = pval, padj = adj_pval, ID = Row.names) %>%
            select(-c(se, df, avg_abundance, n_approx, n_obs)) %>%
            filter(!is.na(padj)) %>%
            arrange(pvalue)
        }
        tableDE(resDE)
        ifHistogram$value <- TRUE
      },
      error = function(e) {
        showModal(modalDialog(
          title = "Running DE analysis failed...",
          "Make sure you have selected the correct experiment options.",
          easyClose = TRUE,
          footer = NULL
        ))}
      )
    })
  })
  
  # filter options for pvalue, log2FC
  filterDE <- reactive({
    if (!is.null(tableDE())) {
      DEtab <- tableDE()
      if(input$ifAdjusted) {
        DEtab <- filter(DEtab, abs(log2FC) >= input$fcFilter, padj <= as.numeric(input$pFilter))
      } else {
        DEtab <- filter(DEtab, abs(log2FC) >= input$fcFilter, pvalue <= as.numeric(input$pFilter))
      }
      DEtab
    }
  })
  
  # table for differentially expressed genes
  output$DEtab <- DT::renderDataTable({
    if (!is.null(filterDE())) {
      if (!is.null(colorRows$row_color)) {
        datatable(filterDE(), selection = 'single', rownames = FALSE,
                  caption = "Differentially expressed genes") %>%
          formatStyle("ID",
                      target = "row",
                      backgroundColor = styleEqual(colorRows$row_priority, 
                                                   colorRows$row_color, 
                                                   default = 'white')) %>%
          formatRound(c('pvalue',"padj"),digits=3) %>%
          formatRound(c('log2FC','stat'),digits=2)
      }
      else {
        datatable(filterDE(), selection = 'single', rownames = FALSE,
                  caption = "Differentially expressed genes") %>%
          formatRound(c('pvalue',"padj"),digits=3) %>%
          formatRound(c('log2FC','stat'),digits=2)
      }
    }
  })
  
  output$downloadTableUI <- renderUI({
    if (!is.null(filterDE())) {
      downloadLink('downloadTable', 'Download current table')
    }
  })
  
  # a button to download DE gene table as tsv
  output$downloadTable <- downloadHandler(
    filename = function() { paste('DE_Table', '.tsv', sep = '') },
    content = function(file) {
      write.table(filterDE(), file = file, quote=FALSE, sep = '\t', col.names = NA)
    }
  )
  
  # volcano plot
  plotV <- reactive({
    dataVolcano <- data.frame(tableDE())
    dataVolcano$ID <- as.character(dataVolcano$ID)
    plot <- ggplot(dataVolcano, aes(x = log2FC, y = -log10(pvalue), label = Gene, customdata = ID)) +
      geom_vline(xintercept = 0, color = "black", linetype = "solid", size = 0.25) +
      geom_vline(xintercept = as.numeric(input$fcFilter), color = "darkgrey", linetype = "dashed") +
      geom_vline(xintercept = -as.numeric(input$fcFilter), color = "darkgrey", linetype = "dashed") +
      geom_hline(yintercept = -log10(as.numeric(input$pFilter)), color = "darkgrey", linetype = "dashed") +  
      annotate(x = 5.0, y = -log10(as.numeric(input$pFilter))-0.1, label = paste("P-value = ", as.numeric(input$pFilter)),
               geom = "text", size = 3, color = "darkgrey") +
      geom_hline(yintercept = -log10(0.25), color="darkgrey", linetype = "dashed") +  
      annotate(x = 5.0, y = 0.5, label = paste("P-value = ", 0.25),
               geom = "text", size=3, color="darkgrey") +
      geom_point(data = dataVolcano[dataVolcano$log2FC >= as.numeric(input$fcFilter) & dataVolcano$pvalue <= as.numeric(input$pFilter),],
                 color="firebrick3", size = 0.9) +
      geom_point(data = dataVolcano[dataVolcano$log2FC <= -as.numeric(input$fcFilter) & dataVolcano$pvalue <= as.numeric(input$pFilter),],
                 color="navy", size=0.9) +
      geom_point(data = dataVolcano[dataVolcano$pvalue > as.numeric(input$pFilter) | (dataVolcano$log2FC < as.numeric(input$fcFilter) & dataVolcano$log2FC > -as.numeric(input$fcFilter)),], color="darkgrey", size = 0.9) +
      xlab("absolute log2(Quantity) difference") +
      ggtitle("Volcano plot") +
      theme(plot.title = element_text(hjust=0.5))
    plot
  })
  
  output$plotVolcano <- renderPlotly({
    p <- ggplotly(plotV(), source = "volcano") 
    p %>%
      event_register("plotly_click")
    if (!is.null(ptHiglight$log2FC)) {
      p <- add_trace(p, x = ptHiglight$log2FC, y = ptHiglight$pValue,
                     type = "scatter", mode = 'markers',
                     marker = list(size = 10, symbol = "star"))
    }
    p
  })
  
  # observe event when a point in the volcano plot is clicked
  observeEvent(event_data("plotly_click", source = "volcano"),{
    # a point in the volcano plot is clicked, turn the ifHistogram to false
    ifHistogram$value <- FALSE
    d <- event_data("plotly_click", source = "volcano")
    lastInfo <- d$customdata
    lastClicked$geneID <- filterDE()[filterDE()$ID == lastInfo,]$ID
    if (input$assay == "Phosphoproteome") {
      lastClicked$geneSymbol <- filterDE()[filterDE()$ID == lastInfo,]$site
    } else {
      lastClicked$geneSymbol <- filterDE()[filterDE()$ID == lastInfo,]$Gene
    }
    colorRows$row_priority <- filterDE()$ID
    colorRows$row_priority <- c(lastClicked$geneID, colorRows$row_priority[colorRows$row_priority != lastClicked$geneID])
    colorRows$row_color <- lapply(colorRows$row_priority, function(x) ifelse(x %in% lastClicked$geneID, "lightgreen", "white"))
    ptHiglight$log2FC <- filterDE()[filterDE()$ID == lastInfo,]$log2FC
    ptHiglight$pValue <- -log10(filterDE()[filterDE()$ID == lastInfo,]$pvalue)
  })
  # observe event when a row in the DE table is clicked
  observeEvent(input$DEtab_row_last_clicked,{
    # a row is clicked, turn the ifHistogram to false
    ifHistogram$value <- FALSE
    lastInfo <- input$DEtab_row_last_clicked
    lastClicked$geneID <- filterDE()[lastInfo,]$ID
    if (input$assay == "Phosphoproteome") {
      lastClicked$geneSymbol <- filterDE()[lastInfo,]$site
    } else {
      lastClicked$geneSymbol <- filterDE()[lastInfo,]$Gene
    }
    ptHiglight$log2FC <- filterDE()[lastInfo,]$log2FC
    ptHiglight$pValue <- -log10(filterDE()[lastInfo,]$pvalue)
    colorRows$row_priority <- c(lastClicked$geneID, colorRows$row_priority[colorRows$row_priority != lastClicked$geneID])
    colorRows$row_color <- lapply(colorRows$row_priority, function(x) ifelse(x %in% lastClicked$geneID, "lightgreen", "white"))
  })
  
  # a ui to hold the plot on the first panel
  output$ui.plot <- renderUI({
    if (!is.null(tableDE())) {
      #the size of the ui is depend on whether it's a histogram or a boxplot
      if(ifHistogram$value == TRUE) {
        fig.width <- 600
        fig.height <- 400
      } else {
        fig.width <- 600
        fig.height <- 500
      }
      plotlyOutput("plot1", width = paste0(fig.width,"px"),
                   height = paste0(fig.height,"px"))
    }
  })
  
  # Histogram plot or boxplot on the first panel
  output$plot1 <- renderPlotly({
    if (!is.null(tableDE())){
      if (ifHistogram$value == TRUE) {
        # plot the histogram of p values
        p <- ggplot(tableDE(), aes(x = pvalue)) +
          geom_histogram(fill = "grey", col = "blue", alpha=0.7) +
          ggtitle("P value histogram") + theme_bw() +
          theme(plot.title = element_text(hjust = 0.5)) 
        ggplotly(p) %>% config(displayModeBar = F)
      }
      # Box-plot for comparison
      else {
        geneID <- lastClicked$geneID
        geneSymbol <- lastClicked$geneSymbol
       
        seqMat <- processedDataSub()
        exprMat <- assays(seqMat)[["Intensity"]]
        
        if(is.null(processedDataSub()$subjectID)) {
          plotTab <- data.frame(group = seqMat$comparison,
                                value = exprMat[geneID,])
          p <- ggplot(plotTab, aes(x= group, y = value)) 
        } else {
          plotTab <- data.frame(group = seqMat$comparison,
                                value = exprMat[geneID,],
                                subjectID = seqMat$subjectID)
          p <- ggplot(plotTab, aes(x= group, y = value, label = subjectID)) +
            geom_line(aes(group = subjectID), linetype = "dotted", color = "grey50")
        }
        
        p <- p + geom_boxplot(aes(fill = group),
                              width = 0.5, alpha = 0.5,
                              outlier.shape = NA) + 
          geom_point() + 
          ylab("Normalized Intensities") + xlab("") + 
          ggtitle(geneSymbol) + theme_bw() + 
          theme(text=element_text(size=15), 
                plot.title = element_text(hjust = 0.5),
                legend.position = "none",
                axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5, size = 15))
      }
    }
  })
  #################################################### time series clustering ##################################################
  
  ####Widgets
  # the selection box to select the treatment
  output$clusterTreatBox <- renderUI({
    allTreat <- unique(processedData()$treatment)
    selectInput("seleTreat_cluster","Select a treatment", allTreat)
  })
  
  # the selection box to select the reference treatment
  output$clusterTreatRefBox <- renderUI({
    allTreat <- unique(processedData()$treatment)
    allTreat <- allTreat[!allTreat %in% input$seleTreat_cluster]
    selectInput("seleTreat_clusterRef","Select a reference treatment", allTreat)
  })
  
  # link to download the cluster table
  output$downloadClusterTabUI <- renderUI({
    if( !is.null (clusterTabVal())) {
      downloadLink("downloadClusterTab","Download cluster table")
    }
  })
  
  output$downloadClusterTab <- downloadHandler(
    filename = function() { paste('clusterTable', '.tsv', sep='') },
    content = function(file) {
      allClusterFeature <- clusterTabVal() %>% distinct(feature, .keep_all = TRUE) %>% .$feature
      if (!is.null(rowData(processedData())$site)) {
        allClusterSite <- rowData(processedData())[allClusterFeature, "site"]
        allClusterSequence <- rowData(processedData())[allClusterFeature, "Sequence"]
      } else {
        allClusterSite <- rep(NA, length(allClusterFeature))
        allClusterSequence <- rep(NA,length(allClusterFeature))
      }
      clustTab <- clusterTabVal() %>% distinct(feature,.keep_all = TRUE) %>%
        mutate(Gene = rowData(processedData())[allClusterFeature, "Gene"],
               site = allClusterSite,
               Sequence = allClusterSequence) %>%
        select(any_of(c("feature", "Gene", "cluster", "prob", "cNum", "Sequence", "site"))) %>%
        dplyr::rename(id = feature, probability = prob,
                      clusterSize = cNum, PhosphoSite = site)
      write.table(clustTab, file = file, quote = FALSE, sep = '\t', col.names = NA)
    }
  )
  
  # selecting time range
  output$timerangeBox <- renderUI({
    # list the time points available to the selected treatment
    processedDataSub <- processedData()[, processedData()$treatment == input$seleTreat_cluster]
    allTimepoint <- unique(processedDataSub$timepoint)
    # remove time points with only 1 sample
    remove1sampleT <- c()
    for (time in allTimepoint) {
      idPerTime <- processedDataSub[, processedDataSub$timepoint == time]$sample
      if (length(idPerTime) < 2)
        remove1sampleT <- append(remove1sampleT, time)
    }
    allTimepoint <- allTimepoint[!allTimepoint %in% remove1sampleT]
    # if a reference treatment is selected, then only list timepoints shared between the two treatments
    if (input$clusterFor == "logFC" | input$clusterFor == "two-condition expression") {
      processedDataRef <- processedData()[, processedData()$treatment == input$seleTreat_clusterRef]
      timepointRef <- unique(processedDataRef$timepoint)
      # remove time points with only 1 sample
      remove1sampleT <- c()
      for (time in timepointRef) {
        idPerTime <- processedDataRef[, processedDataRef$timepoint == time]$sample
        if (length(idPerTime) < 2)
          remove1sampleT <- append(remove1sampleT, time)
      }
      timepointRef <- timepointRef[!timepointRef %in% remove1sampleT]
      allTimepoint <- intersect(allTimepoint, timepointRef)  # select only the intersection
    }
    checkboxGroupInput("seleTimeRange", "Time points to include",
                       allTimepoint, allTimepoint, inline = TRUE)
  })
  
  exprMatObj <- reactive({
    # Processing with selecting expression (i.e. only 1 condition)
    if (input$clusterFor == "expression") {
      processedDataSub <- processedData()[, processedData()$treatment == input$seleTreat_cluster & 
                                  processedData()$timepoint %in% input$seleTimeRange]
      assayMat <- assay(processedDataSub)
      # Filtering genes with p values from the spline fit test
      if (input$ifFilterFit) {
        if (!is.null(processedDataSub$subjectID)) {
          assayMat <- splineFilter(assayMat, subjectID = processedDataSub$subjectID,
                                   time = processedDataSub$timepoint,
                                   df = length(unique(processedDataSub$timepoint))-1,
                                   pCut = as.numeric(input$pSpline),
                                   ifFDR = input$ifSplineFdr)
        } else {
          assayMat <- splineFilter(assayMat, subjectID = NULL,
                                   time = processedDataSub$timepoint,
                                   df = length(unique(processedDataSub$timepoint))-1,
                                   pCut = as.numeric(input$pSpline),
                                   ifFDR = input$ifSplineFdr)
        }
        processedDataSub <- processedDataSub[rownames(assayMat),]
      }
      exprMat <- lapply(unique(processedDataSub$timepoint), function(tp) {
        rowMedians(assayMat[,processedDataSub$timepoint == tp])
      }) %>% bind_cols() %>% as.matrix()
      
      rownames(exprMat) <- rownames(processedDataSub)
      colnames(exprMat) <- unique(processedDataSub$timepoint)
      
    } else if (input$clusterFor == "logFC") {
      processedDataSub <- processedData()[,processedData()$treatment == input$seleTreat_cluster & 
                                  processedData()$timepoint %in% input$seleTimeRange]
      processedDataRef <- processedData()[,processedData()$treatment == input$seleTreat_clusterRef & 
                                  processedData()$timepoint %in% input$seleTimeRange]
      # Match processedDataSub and processedDataRef by subjectID and timepoint if subjectID is available
      if (!is.null(processedData()$subjectID)) {
        processedDataSub <- processedDataSub[,match(paste0(processedDataRef$subjectID,"_",processedDataRef$timepoint),
                                                    paste0(processedDataSub$subjectID,"_",processedDataSub$timepoint))] # make sure order is the same 
      } else {
        # arrange processedDataSub and processedDataRef by timepoints
        posRef <- posSub <- c()
        allTimePoint <- unique(processedDataRef$timepoint)
        for (time in unique(processedDataRef$timepoint)) {
          pointerSub <- which(processedDataSub$timepoint == time)
          pointerRef <- which(processedDataRef$timepoint == time)
          posSub <- append(posSub, pointerSub)
          posRef <- append(posRef, pointerRef)
        }
        processedDataSub <- processedDataSub[,posSub]
        processedDataRef <- processedDataRef[,posRef]
      } 
      assayMat <- assay(processedDataSub)
      RefMat <- assay(processedDataRef)
      
      #  calculate fold change per matching sample first, apply spline filter, THEN calculate the mean logFC per time point
      
      if (input$ifFilterFit) {
        fcMat <- assayMat - RefMat
        if (!is.null(processedDataSub$subjectID)) {
          fcMat <- splineFilter(fcMat, subjectID = processedDataSub$subjectID,
                                time = processedDataSub$timepoint,
                                df = length(unique(processedDataSub$timepoint))-1,
                                pCut = as.numeric(input$pSpline),
                                ifFDR = input$ifSplineFdr)
        } else { # i.e. if subjectID is not present
          fcMat <- splineFilter(fcMat, subjectID = NULL,
                                time = processedDataSub$timepoint,
                                df = length(unique(processedDataSub$timepoint))-1,
                                pCut = as.numeric(input$pSpline),
                                ifFDR = input$ifSplineFdr)
        }
        processedDataSub <- processedDataSub[rownames(fcMat),]
        exprMat <- lapply(unique(processedDataSub$timepoint), function(tp) {
          rowMeans(fcMat[,processedDataSub$timepoint == tp])
        }) %>% bind_cols() %>% as.matrix()
      } else { # i.e. if choose not to spline filter
        
        # calculate the mean intensities per time point, THEN calculate the logFC
        
        assayMatMean <- lapply(unique(processedDataSub$timepoint), function(tp) {
          rowMeans(assayMat[,processedDataSub$timepoint == tp])
        }) %>% bind_cols() %>% as.matrix()
        RefMatMean <- lapply(unique(processedDataSub$timepoint), function(tp) {
          rowMeans(RefMat[,processedDataRef$timepoint == tp])
        }) %>% bind_cols() %>% as.matrix()
        exprMat <- assayMatMean -  RefMatMean
      }
      rownames(exprMat) <- rownames(processedDataSub)
      colnames(exprMat) <- unique(processedDataSub$timepoint)
      
    } else if (input$clusterFor == "two-condition expression") {
      processedDataSub <- processedData()[,processedData()$treatment %in% c(input$seleTreat_cluster, input$seleTreat_clusterRef) & 
                                  processedData()$timepoint %in% input$seleTimeRange]
      assayMat <- assay(processedDataSub)
      if (input$ifFilterFit) {
        if (!is.null(processedDataSub$subjectID)) {
          assayMat <- splineFilter(assayMat, subjectID = processedDataSub$subjectID,
                                   time = processedDataSub$timepoint,
                                   treatment = processedDataSub$treatment,
                                   refTreatment = input$seleTreat_clusterRef,
                                   df = length(unique(processedDataSub$timepoint))-1,
                                   pCut = as.numeric(input$pSpline),
                                   ifFDR = input$ifSplineFdr)
        } else {
          assayMat <- splineFilter(assayMat, subjectID = NULL,
                                   time = processedDataSub$timepoint,
                                   treatment = processedDataSub$treatment,
                                   refTreatment = input$seleTreat_clusterRef,
                                   df = length(unique(processedDataSub$timepoint))-1,
                                   pCut = as.numeric(input$pSpline),
                                   ifFDR = input$ifSplineFdr)
        }
        processedDataSub <- processedDataSub[rownames(assayMat),]
      }
      processedDataSub$timeTreat <- paste0(processedDataSub$timepoint,"_", processedDataSub$treatment)
      exprMat <- lapply(unique(processedDataSub$timeTreat), function(tt) {
        rowMedians(assayMat[,processedDataSub$timeTreat == tt])
      }) %>% bind_cols() %>% as.matrix()
      
      rownames(exprMat) <- rownames(processedDataSub)
      colnames(exprMat) <- unique(processedDataSub$timeTreat)
    }
    # filter based on variance 
    sds <- apply(exprMat,1,sd)
    varPer <- as.numeric(input$topVarTime)
    exprMat <- exprMat[order(sds, decreasing = TRUE)[seq(1, varPer/100*nrow(exprMat))], ]
    
    # only center when it's for expression
    if (input$clusterFor != "logFC") exprMat <- mscale(exprMat)
    
    # remove NA values
    exprMat <- exprMat[complete.cases(exprMat), ]
    exprMat
  })
  
  clusterPlotVal <- reactiveVal()
  clusterTabVal <- reactiveVal()
  
  # (currently not in use) plot to find out optimal number of clusters 
  observeEvent(input$plotSilhouette, {
    withProgress(message = "Calculating Silhouette and WSS scores, please wait...", {
      d <- exprMatObj()
      p1 <- factoextra::fviz_nbclust(d,  kmeans, c("silhouette"), k.max = 20)
      p2 <- factoextra::fviz_nbclust(d,  kmeans, c("wss"), k.max = 20)
    })
    
    clusterPlotVal(cowplot::plot_grid(p1,p2,NULL,NULL, ncol=2))
    clustNum(6)
  })
  
  # plot time-series clustering result
  observeEvent(input$runCluster, {
    withProgress(message = "Performing cmeans clustering, please wait...", {
      # plot clustering result and error check
      tryCatch({
        if (input$clusterFor != "two-condition expression")
        {clusterRes <- clusterTS(exprMatObj(),
                                 as.numeric(input$seleNumCluster),
                                 pCut = as.numeric(input$seleProbCut))}
        else {
          clusterRes <- clusterTS(exprMatObj(),
                                  as.numeric(input$seleNumCluster),
                                  pCut = as.numeric(input$seleProbCut),
                                  twoCondition = TRUE)
        }
        clusterPlotVal(clusterRes$plot)
        clusterTabVal(clusterRes$cluster)
        clustNum(input$seleNumCluster)
        if (nrow(data.frame(clusterRes$cluster)) == 0) {
          showModal(modalDialog(
            title = "No cluster found...",
            "Try again with a different setting (e.g., time points included, number of clusters, top % variant).",
            easyClose = TRUE,
            footer = NULL
          ))
        }
      },
      error = function(e) {
        showModal(modalDialog(
          title = "Running Time series clustering failed...",
          "Try again with a different setting (e.g., time points included, number of clusters, top % variant).",
          easyClose = TRUE,
          footer = NULL
        ))})
    })
  })
  
  # table showing the proteins or phosphopeptides in a selected cluster
  selectedCluster <- reactive({
    if (is.null(clusterTabVal())) {
      NULL 
    } else {
      selectedTab <- filter(clusterTabVal(), cluster == input$seleCluster) %>% 
        distinct(feature, .keep_all = TRUE) %>%
        mutate(feature = as.character(feature)) 
      clusterData <- rowData(processedData()[selectedTab$feature,])
      # If the data is phosphoproteomic: include the columns indicating phosphorylation site and peptide sequence
      if ((!is.null(clusterData$site)) & (!is.null(clusterData$Sequence))) {
        selectedTab <- selectedTab %>%
          mutate(Sequence = clusterData$Sequence,
                 site = clusterData$site)
      }
      selectedTab <- selectedTab %>% 
        mutate(Gene = clusterData$Gene,
               UniprotID = clusterData$UniprotID,
               prob = formatC(prob, digits=1)) %>%
        select(any_of(c("feature", "Gene", "site", "prob", "cluster", "UniprotID", "Sequence"))) %>%
        arrange(desc(prob)) %>%
        dplyr::rename(ID = feature, probability = prob)
      selectedTab
    }
  })
  
  # selection of a cluster
  output$seleClusterBox <- renderUI({
    if(!is.null(clusterTabVal())) {
      selectInput("seleCluster", "Select a cluster", sort(unique(clusterTabVal()$cluster)))
    }
  })
  
  # report number of genes for clustering
  output$numGeneCluster <- renderText({
    sprintf("Number of genes for clustering: %s", nrow(exprMatObj()))
  })
  
  clustNum <- reactiveVal()
  
  output$clusterPlotUI <- renderUI({
    if (!is.null(clustNum())) {
      plotOutput("clusterPlot",
                 height = 250*ceiling(as.numeric(clustNum())/3),
                 width = 800)
    }
  })
  
  output$clusterPlot <- renderPlot({
    clusterPlotVal()
  })
  
  output$eachClusterTab <- DT::renderDataTable({
    DT::datatable(selectedCluster(), selection = 'single')
  })
  
  # Plot the time-series of individual proteins
  output$clusterTimePlot <- renderPlot({
    lastClicked <- input$eachClusterTab_row_last_clicked
    if (!is.null(lastClicked)) {
      geneID <- selectedCluster()[lastClicked,]$ID
      if (input$assay == "Phosphoproteome")
        geneSymbol <- selectedCluster()[lastClicked,]$site else
          geneSymbol <- selectedCluster()[lastClicked,]$Gene
      
      if (input$clusterFor == "expression") {
        seqMat <- processedData()[,processedData()$treatment == input$seleTreat_cluster & 
                                    processedData()$timepoint %in% input$seleTimeRange]
        yLabText <- "Normalized expression"
      } else if (input$clusterFor == "logFC"){
        seqMat <- processedData()[,processedData()$treatment == input$seleTreat_cluster & 
                                    processedData()$timepoint %in% input$seleTimeRange]
        RefMat <- processedData()[,processedData()$treatment == input$seleTreat_clusterRef & 
                                    processedData()$timepoint %in% input$seleTimeRange]
        # Match seqMat and RefMat by subjectID and timepoint if subjectID is available,
        # otherwise match by only timepoints
        if (!is.null(processedData()$subjectID)) {
          seqMat <- seqMat[,match(paste0(RefMat$subjectID,"_",RefMat$timepoint),
                                  paste0(seqMat$subjectID,"_",seqMat$timepoint))]  # make sure order is the same 
        }
        else {
          # arrange processedDataSub and processedDataRef by timepoints
          posRef <- posSub <- c()
          for (time in unique(RefMat$timepoint)) {
            pointerSub <- which(seqMat$timepoint == time)
            pointerRef <- which(RefMat$timepoint == time)
            posSub <- append(posSub, pointerSub)
            posRef <- append(posRef, pointerRef)
          }
          seqMat <- seqMat[,posSub]
          RefMat <- RefMat[,posRef]
        }
        # compute the fold change
        assay(seqMat) <- assay(seqMat) - assay(RefMat)
        yLabText <- "logFC"
      } else if (input$clusterFor == "two-condition expression") {
        seqMat <- processedData()[,processedData()$treatment %in% c(input$seleTreat_cluster,input$seleTreat_clusterRef) &
                                    processedData()$timepoint %in% input$seleTimeRange]
        yLabText <- "Normalized expression"
      }
      
      # Dotplot of the selected (phospho) protein level
      # time is treated as a numerical variable
      plotTab <- data.frame(time = seqMat$timepoint,
                            value = assay(seqMat)[geneID,], 
                            treatment = as.character(seqMat$treatment)) 
      timeUnit <- str_extract(plotTab$time, "h|min")
      timeUnit <- ifelse(is.na(timeUnit), "", timeUnit)
      if ((any(timeUnit == "h")) & (any(timeUnit == "min"))) {
        plotTab[timeUnit == "min","time"] <- 1/60 * as.numeric(gsub("min", "", plotTab[timeUnit == "min","time"]))
      } 
      plotTab$time <- as.numeric(gsub("h|min", "", plotTab$time))
      p <- ggplot(plotTab, aes(x= time, y = value)) +
        geom_point(aes(color = treatment), size=3) + 
        stat_summary(aes(color=paste("mean",treatment)),fun = mean, geom = "line", linewidth = 2)+
        ylab(yLabText) + xlab("time") + 
        ggtitle(geneSymbol) + theme_bw() + 
        theme(text=element_text(size=15),plot.title = element_text(hjust = 0.5),
              legend.position = "bottom",
              axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5, size=15))
      # connect the dots across time points by subjectID if provided
      if (!is.null(seqMat$subjectID)) {
        plotTab$subjectID <- seqMat$subjectID
        plotTab <- plotTab %>% mutate(patCondi = paste0(subjectID,"_",treatment))
        p <- p + geom_line(aes(group = patCondi), linetype = "dotted", color = "grey50") 
      }
      p
    } else {
      NULL
    }
  })
  
  ################################## Enrichment analysis on DE or time-series clustering result #############################
  
  ####Widgets
  
  # to store the color values used to show gene sets that contain a certain gene
  colorList <- reactiveValues(cols = c())
  
  # a reactive variable to store GSE result
  GSEres <- reactiveValues(resTab = NULL, resObj = NULL)
  
  # a value to check whether enrichment tab is clicked
  clickRecord <- reactiveValues(enrich = FALSE, gene = FALSE, kinase = FALSE)
  
  # function to run GSEA 
  # Note: the analysisMethod is default to Pathway enrichment if the Proteome
  # assay is selected (might change in the future!)
  resGSEA <- observeEvent(input$RunEnrich, {
    if ((input$seleSourceEnrich == "Differential expression") & (input$analysisMethod == "Pathway enrichment" || input$assay == "Proteome")) {
      ################################################################
      ####### Pathway enrichment for differential expression #########
      ################################################################
      if (!is.null(filterDE())) {
        output$errMsg <- renderText("")
        withProgress(message = "Running enrichment analysis, please wait...", {
          # set color list to empty
          colorList$cols <- NULL
          # reading geneset database
          inGMT <- loadGSC(paste0("geneset/",input$sigSet),type="gmt")
          # method
          gseMethod <- input$enrichMethod
          # parameters for GSEA
          if (gseMethod == "GSEA") nPerm <- input$permNum
          # processing differential expression
          corTab <- filterDE() %>%
            arrange(pvalue) %>%
            filter(!duplicated(Gene)) %>%
            arrange(stat)
          # gene level statistics based on user input
          if(input$statType == "stat") {
            myCoef <- data.frame(row.names = corTab$Gene,
                                 stat = corTab$stat,
                                 stringsAsFactors = FALSE)
          } else {
            myCoef <- data.frame(row.names = corTab$Gene,
                                 stat = corTab$log2FC,
                                 stringsAsFactors = FALSE)
          }
          # perform gene set analysis based on selected statistical GSA method
          if (gseMethod == "PAGE") {
            res <- runGSA(geneLevelStats = myCoef,
                          geneSetStat = "page",
                          adjMethod = "fdr",
                          gsc = inGMT,
                          signifMethod = 'nullDist')
          } else if (gseMethod == "GSEA") {
            res <- runGSA(geneLevelStats = myCoef,
                          geneSetStat = "gsea",
                          adjMethod = "fdr",
                          gsc = inGMT,
                          signifMethod = 'geneSampling',
                          nPerm = nPerm)
          }
          
          resTab <- GSAsummaryTable(res)
          colnames(resTab) <- c("Name", "Gene Number", "Stat", "p.up", "p.up.adj",
                                "p.down", "p.down.adj", "Number up", "Number down")
          if(input$ifEnrichFDR) {
            resTab <- filter(resTab,
                             p.up.adj <= input$sigLevel | p.down.adj <= input$sigLevel) %>%
              arrange(desc(Stat))
          } else {
            resTab <- filter(resTab,
                             p.up <= input$sigLevel | p.down <= input$sigLevel) %>%
              arrange(desc(Stat))
          }
          # if all genes are filtered out: show a notification and set
          # GSEres$resTab to NULL to avoid causing error later on
          if (nrow(resTab) > 0) {
            GSEres$resTab <- resTab
          } else {
            GSEres$resTab <- NULL
            showModal(modalDialog(
              title = "No enrichment found...",
              "Try again with a higher p-value cut-off.",
              easyClose = TRUE,
              footer = NULL
            ))}
          GSEres$resTab <- resTab
          GSEres$resObj <- res
          clickRecord$enrich <- FALSE
          clickRecord$gene <- FALSE
        })
      } else {
        output$errMsg <- renderText("Please perform differential expression analysis first or load a previous result!")
      }
    } else if ((input$seleSourceEnrich =="Time series cluster") & (input$analysisMethod =="Pathway enrichment" || input$assay == "Proteome")) {
      ################################################################
      ######### Pathway enrichment for time series cluster ###########
      ################################################################
      if (!is.null(selectedCluster())) {
        output$errMsg <- renderText("")
        withProgress(message = "Running enrichment analysis, please wait...", {
          
          # set color list to empty
          colorList$cols <- NULL
          # reading geneset database
          inGMT <- loadGSC(paste0("geneset/",input$sigSet),type="gmt")
          # Applying the Fisher's exact test with the runFisher function from utils.
          # Other enrichment methods can be added here
          if (input$enrichMethod1 == "Fisher's exact test")
            resTab <- runFisher(unique(selectedCluster()$Gene),
                                reference = unique(rowData(processedData())$Gene),
                                gmtFile = paste0("geneset/", input$sigSet)) %>%
            arrange(pval)
          # Filter by the p-value threshold (input$sigLevel)
          if (input$ifEnrichFDR) {
            resTab <- filter(resTab, padj <= input$sigLevel) %>% arrange(padj)
          } else {
            resTab <- filter(resTab, pval <= input$sigLevel) %>% arrange(pval)
          }
          # if all genes are filtered out: show a notification and set GSEres$resTab
          # to NULL to avoid causing error later on
          if (nrow(resTab) > 0) {
            GSEres$resTab <- resTab
          } else {
            GSEres$resTab <- NULL
            showModal(modalDialog(
              title = "No enrichment found...",
              "Try again with a higher p-value cut-off.",
              easyClose = TRUE,
              footer = NULL
            ))}
          GSEres$resObj <- NULL
          clickRecord$enrich <- FALSE
          clickRecord$gene <- FALSE
        })
      } else {
        output$errMsg <- renderText("Please perform time series clustering first!")
      }
    } else if ((input$seleSourceEnrich == "Differential expression") & (input$analysisMethod == "Phospho-signature enrichment")){
      ################################################################
      ### Phospho-signature enrichment for differential expression ###
      ################################################################
      if ((!is.null(filterDE())) & (input$assay == "Phosphoproteome")) {
        output$errMsg <- renderText("")
        withProgress(message = "Running enrichment analysis, please wait...", {
          # preprocess the input dataframe
          inputTab <- filterDE() %>%
            arrange(pvalue) %>%
            filter(!duplicated(site)) %>%
            arrange(desc(stat))
          # select to use either t-statistic or log2FC
          if (input$statType == "stat") {
            myCoef <- data.frame(row.names = inputTab$site,
                                 stat = inputTab$stat,
                                 stringsAsFactors = FALSE)
          } else if (input$statType == "log2FC") {
            myCoef <- data.frame(row.names = inputTab$site,
                                 stat = inputTab$log2FC,
                                 stringsAsFactors = FALSE)
          }
          # retrieve the phosphodatabase
          ptmSetDb <- read.table(paste0("ptmset/", input$sigSetPTM), header = TRUE, sep = "\t",stringsAsFactors = FALSE)
          
          # perform GSEA
          resTab <- runGSEAforPhospho(geneStat = myCoef, ptmSetDb = ptmSetDb, nPerm =  input$permNum,
                                      weight = 1, correl.type = "rank", statistic = "Kolmogorov-Smirnov", min.overlap = 5) %>%
            as.data.frame()
          colnames(resTab) <- c("Name", "Site.number", "Stat", "Number.pSite.Db", "Number.PTM.site.Db",
                                "pvalue", "Number.up", "Number.down", "padj")
          resTab <- resTab %>%
            select(Name,Site.number,Stat,Number.up,Number.down,Number.pSite.Db,
                   Number.PTM.site.Db, pvalue, padj) # rearrange column order 
          if (input$ifEnrichFDR) {
            resTab <- filter(resTab, padj <= input$sigLevel) %>%
              arrange(desc(Stat))
          } else {
            resTab <- filter(resTab, pvalue <= input$sigLevel) %>%
              arrange(desc(Stat))
          }
          # if all genes are filtered out: show a notification and set
          # GSEres$resTab to NULL to avoid causing error later on
          if (nrow(resTab) > 0) {
            GSEres$resTab <- resTab
          } else {
            GSEres$resTab <- NULL
            showModal(modalDialog(
              title = "No enrichment found...",
              "Try again with a higher p-value cut-off.",
              easyClose = TRUE,
              footer = NULL
            ))}
          GSEres$resTab <- resTab
          #GSEres$resObj <- res
          clickRecord$enrich <- FALSE
          clickRecord$gene <- FALSE
        })
      } else output$errMsg <- renderText("Please perform differential expression analysis first and make sure you have selected the Phosphoproteome assay!")
    } else if ((input$seleSourceEnrich == "Time series cluster") & (input$analysisMethod == "Phospho-signature enrichment")) {
      ################################################################
      ### Phospho-signature enrichment for Time-series clustering ####
      ################################################################
      if ((!is.null(selectedCluster())) & (input$assay == "Phosphoproteome")) {
        output$errMsg <- renderText("")
        withProgress(message = "Running enrichment analysis, please wait...", {
          # set color list to empty
          colorList$cols <- NULL
          # Run the Fisher's exact test for sites in the cluster
          resTab <- runFisher(genes = selectedCluster()$site,
                              reference = rowData(processedData())$site,
                              gmtFile = paste0("ptmset/", input$sigSetPTM),
                              ptm = TRUE) %>%
            rename(Site.number = "Gene.number")
          # Filter by the p-value threshold (input$sigLevel)
          if (input$ifEnrichFDR) {
            resTab <- filter(resTab, padj <= input$sigLevel) %>% arrange(padj)
          } else {
            resTab <- filter(resTab, pval <= input$sigLevel) %>% arrange(pval)
          }
          # if all genes are filtered out: show a notification and set GSEres$resTab
          # to NULL to avoid causing error later on
          if (nrow(resTab) > 0) {
            GSEres$resTab <- resTab
          } else {
            GSEres$resTab <- NULL
            showModal(modalDialog(
              title = "No enrichment found...",
              "Try again with a higher p-value cut-off.",
              easyClose = TRUE,
              footer = NULL
            ))}
          GSEres$resObj <- NULL
          clickRecord$enrich <- FALSE
          clickRecord$gene <- FALSE
        })
      } else output$errMsg <- renderText("Please perform time series clustering first and make sure you have selected the Phopshoproteome assay!")
    }
  })
  
  # show differential expressed genes or enrichment results on table 1 (up table)
  output$enrichTab <- DT::renderDataTable({
    
    if( !is.null (GSEres$resTab)) {
      resTab <- GSEres$resTab 
      
      if (input$seleSourceEnrich == "Differential expression") {
        if (!is.null(colorList$cols)) {  
          # when the bottom table is clicked, color the pathways according to
          # whether it contains the clicked gene 
          datatable(resTab,selection = 'single', caption = "Enriched gene/PTM signature sets") %>%
            formatStyle('Stat',background = styleInterval(c(0), c("lightblue", "pink"))) %>%
            formatStyle('Name', color = styleEqual(resTab$Name, colorList$cols)) %>%
            formatRound(which(colnames(resTab) %in% c("Stat", "pvalue", "padj", "p.up","p.up.adj","p.down","p.down","p.down.adj")), digits=3)
        } else { 
          # when bottom table was not clicked, do not show colors
          datatable(resTab,selection = 'single', caption = "Enriched gene/PTM signature sets") %>% 
            formatStyle('Stat', background = styleInterval(c(0), c("lightblue", "pink"))) %>%
            formatRound(which(colnames(resTab) %in% c("Stat", "pvalue", "padj", "p.up","p.up.adj","p.down","p.down","p.down.adj")), digits = 3)
        }
      } else {
        if (!is.null(colorList$cols)) {  
          # when the bottom table is clicked, color the pathways according to whether it contains the clicked gene 
          datatable(resTab,selection = 'single', caption = "Enriched gene/PTM signature sets") %>%
            formatStyle('Name', color = styleEqual(resTab$Name, colorList$cols)) %>%
            formatRound(c('pval', 'padj'), digits=3)
        } else { 
          # when bottom table was not clicked, do not show colors
          datatable(resTab,selection = 'single', caption = "Enriched gene/PTM signature sets") %>% 
            formatRound(c('pval', 'padj'), digits = 3)
        }
      }
    }
  })
  
  # find gene sets that contain a certain gene
  setGene <- reactive({
    resTab <- GSEres$resTab
    if (input$analysisMethod == "Pathway enrichment" | input$assay == "Proteome")
      setList <- loadGSC(paste0("geneset/",input$sigSet),type="gmt")$gsc[resTab$Name] else {
        setList <- read.table(paste0("ptmset/", input$sigSetPTM), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
        if (input$seleSourceEnrich == "Time series cluster")
          setList <- setList %>% mutate(signature = ifelse(site.direction == "u", paste0(signature,"_upregulated"), paste0(signature, "_downregulated")))
        setList <- setList %>% 
          filter(signature %in% resTab$Name, 
                 site.ptm == "p") %>%
          separate(site.annotation, sep=":", into = c("site", "PubMedID"), extra = "merge", fill="right")
      }
    if (input$seleSourceEnrich == "Differential expression") {
      if (input$analysisMethod == "Pathway enrichment" | input$assay == "Proteome")
        genes <- filterDE()$Gene else
          sites <- filterDE()$site
    } else { # aka if input$seleSourceEnrich == "Time series clustering"
      if (input$analysisMethod == "Pathway enrichment" | input$assay == "Proteome")
        genes <- unique(selectedCluster()$Gene) else
          sites <- selectedCluster()$site
    }
    if (input$analysisMethod == "Pathway enrichment" | input$assay == "Proteome")
      allSets <- sapply(genes, function(geneName) names(setList)[sapply(setList, function(x) geneName %in% x)]) else {
        allSets <- list()
        for (site in sites) {
          if (site %in% setList$site)
            allSets[[site]] <- setList[setList$site == site,"signature"]}
      }
    allSets <- allSets[sapply(allSets, function(x) length(x) != 0)]
    allSets
  })
  
  # list genes that enriched in a certain gene set
  gseaList <- reactive({
    
    setName <- GSEres$resTab[as.integer(input$enrichTab_row_last_clicked), "Name"]
    if (input$assay == "Proteome" | input$analysisMethod == "Pathway enrichment") {
      geneList <- loadGSC(paste0("geneset/", input$sigSet), type = "gmt")$gsc[[setName]]
    } else {
      geneList <- read.table(paste0("ptmset/", input$sigSetPTM), header = T, sep = "\t", stringsAsFactors = F) 
      if (input$seleSourceEnrich == "Time series cluster")
        geneList <- geneList %>% mutate(signature = ifelse(site.direction == "u", paste0(signature,"_upregulated"), paste0(signature, "_downregulated")))
      geneList <- geneList %>%  
        filter(signature == setName, site.ptm == "p") %>%
        separate(site.annotation, sep=":", into = c("site", "PubMedID"), extra = "merge", fill="right")
    }
    
    if (input$seleSourceEnrich == "Differential expression") {
      # Differential expression
      corGene <- filterDE()
      if (input$analysisMethod == "Pathway enrichment" | input$assay == "Proteome") {
        geneTab <- corGene[corGene$Gene %in% geneList,]
        # setNum is the number of gene sets containing the gene in question
        geneTab$setNum <- sapply(geneTab$Gene, function(x) length(setGene()[[x]]))
      } else {
        geneTab <- corGene[corGene$site %in% geneList$site,]
        geneTab <- merge(x = geneTab, y = geneList[,c("site","site.direction","PubMedID")], by = "site", all.x = TRUE)
        # only select sites whose direction of regulation (up- or down-regulated) matches with the database
        geneTab <- geneTab[((geneTab$log2FC >=0) & (geneTab$site.direction == "u")) | ((geneTab$log2FC < 0) & (geneTab$site.direction ==  "d")),]
        geneTab$setNum <- sapply(geneTab$site, function(x) length(setGene()[[x]]))
      }
      geneTab <- select(geneTab, any_of(c("site", "Gene","log2FC", "pvalue", "padj", "setNum", "Sequence", "PubMedID", "ID")))
    } else {
      # time series cluster
      corGene <- selectedCluster()
      if (input$analysisMethod == "Pathway enrichment" | input$assay == "Proteome") {
        geneTab <- corGene[corGene$Gene %in% geneList,]
        # setNum is the number of gene sets containing the gene in question
        geneTab$setNum <- sapply(geneTab$Gene, function(x) length(setGene()[[x]]))
      } else {
        geneTab <- corGene[corGene$site %in% geneList$site,]
        geneTab <- merge(x = geneTab, y = geneList[,c("site","PubMedID")], by = "site", all.x = TRUE)
        geneTab$setNum <- sapply(geneTab$site, function(x) length(setGene()[[x]]))
      }
      geneTab <- select(geneTab, any_of(c("site", "Gene", "cluster", "probability", "setNum", "Sequence", "PubMedID", "ID")))
    }
    geneTab <- geneTab %>% mutate_if(is.numeric, formatC, digits = 2)
    geneTab
  })
  
  # if the enrichtab is clicked, cancel the color
  observeEvent(input$enrichTab_row_last_clicked, {
    colorList$cols <- NULL
    clickRecord$enrich <- TRUE
  })
  
  # get the gene ID as well as set color of the gene sets when click the bottom table
  observeEvent(input$geneTab_row_last_clicked, {
    if (input$analysisMethod == "Pathway enrichment" | input$assay == "Proteome")
      clickSym <- gseaList()[as.integer(input$geneTab_row_last_clicked),"Gene"][[1]] else
        clickSym <- gseaList()[as.integer(input$geneTab_row_last_clicked), "site"][[1]]
    colorList$cols <- sapply(GSEres$resTab$Name, function(x) ifelse(x %in% setGene()[[clickSym]], "red", "black"))
    clickRecord$gene <- TRUE
  })
  
  # using the bottom right table to show genes enriched in the selected set
  output$geneTab <- DT::renderDataTable({
    if (clickRecord$enrich) {
      datatable(select(gseaList(), -ID),
                selection = 'single', rownames = FALSE, 
                caption = "Genes/Phosphosites in the selected set") 
    }
  })
  
  output$plot2 <- renderPlot({
    if (clickRecord$gene) {
      lastClicked <- input$geneTab_row_last_clicked
      geneID <- gseaList()[lastClicked,]$ID
      if (input$assay == "Phosphoproteome")
        geneSymbol <- gseaList()[lastClicked,]$site else
          geneSymbol <- gseaList()[lastClicked,]$Gene
      
      if (!is.null(lastClicked)) {
        
        if (input$seleSourceEnrich == "Differential expression") {
          seqMat <- processedDataSub()
          #assayName <- ifelse(input$deMethod == "limma","voom","vst")
          
          if (input$seleCompare == "treatments per time") {
            # Include patient ID in the plot (subjectID) if provided
            if (is.null(seqMat$subjectID)) {
              plotTab <- data.frame(group = seqMat$comparison,
                                    #  value = assays(processedData()[,colnames(seqMat)])[[assayName]][geneID,])
                                    value = assays(processedData()[,colnames(seqMat)])[["Intensity"]][geneID,])
              p <- ggplot(plotTab, aes(x= group, y = value)) + 
                geom_boxplot(aes(fill = group), 
                             width = 0.5, alpha = 0.5, outlier.shape = NA) + 
                geom_point() 
            } else {
              plotTab <- data.frame(group = seqMat$comparison,
                                    value = assays(processedData()[,colnames(seqMat)])[[assayName]][geneID,], 
                                    subjectID = seqMat$subjectID)
              
              p <- ggplot(plotTab, aes(x= group, y = value, 
                                       label = subjectID)) + 
                geom_boxplot(aes(fill = group), 
                             width = 0.5, alpha = 0.5, outlier.shape = NA) + 
                geom_point() + 
                geom_line(aes(group = subjectID), linetype = "dotted", color = "grey50")
            }
            p <- p + ylab("Normalized expression") + xlab("") + 
              ggtitle(geneSymbol) + theme_bw() + 
              theme(text=element_text(size=15),plot.title = element_text(hjust = 0.5),
                    legend.position = "none",
                    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5, size=15))
            p
            
          } else {
            if (is.null(seqMat$subjectID)) {
              plotTab <- data.frame(group = seqMat$comparison,
                                    value = assays(processedData()[,colnames(seqMat)])[["Intensity"]][geneID,], 
                                    treatment = seqMat$treatment)
              
              p <- ggplot(plotTab, aes(x= group, y = value)) + 
                geom_boxplot(aes(fill = group), 
                             width = 0.5, alpha = 0.5, outlier.shape = NA) + 
                geom_point() 
            } else {
              plotTab <- data.frame(group = seqMat$comparison,
                                    value = assays(processedData()[,colnames(seqMat)])[[assayName]][geneID,], 
                                    subjectID = seqMat$subjectID,
                                    treatment = seqMat$treatment)
              
              p <- ggplot(plotTab, aes(x= group, y = value, 
                                       label = subjectID)) + 
                geom_boxplot(aes(fill = group), 
                             width = 0.5, alpha = 0.5, outlier.shape = NA) + 
                geom_point() + 
                geom_line(aes(group = subjectID), linetype = "dotted", color = "grey50")
            }
            p <- p + ylab("Normalized expression") + xlab("") + 
              ggtitle(geneSymbol) + theme_bw() + facet_wrap(~treatment)+
              theme(text=element_text(size=15),plot.title = element_text(hjust = 0.5),
                    legend.position = "none",
                    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5, size=15))
            p
            
          }
        } else if (input$seleSourceEnrich == "Time series cluster") {
          
          if (input$clusterFor == "expression") {
            seqMat <- processedData()[,processedData()$treatment == input$seleTreat_cluster & processedData()$timepoint %in% input$seleTimeRange]
            yLabText <- "Normalized expression"
          } else if (input$clusterFor == "logFC"){
            seqMat <- processedData()[,processedData()$treatment == input$seleTreat_cluster & processedData()$timepoint %in% input$seleTimeRange]
            RefMat <- processedData()[,processedData()$treatment == input$seleTreat_clusterRef & processedData()$timepoint %in% input$seleTimeRange]
            # Match seqMat and RefMat by subjectID and time point OR by time point and replicate (if subjectID not provided)
            if (!is.null(RefMat$subjectID)) {
              seqMat <- seqMat[,match(paste0(RefMat$subjectID,"_",RefMat$timepoint), paste0(seqMat$subjectID,"_",seqMat$timepoint))] #make sure order is the same 
            } else {
              posRef <- posSub <- c()
              for (time in unique(RefMat$timepoint)) {
                pointerSub <- which(seqMat$timepoint == time)
                pointerRef <- which(RefMat$timepoint == time)
                posSub <- append(posSub, pointerSub)
                posRef <- append(posRef, pointerRef)
              }
              seqMat <- seqMat[,posSub]
              RefMat <- RefMat[,posRef]
            # compute the fold change
            assay(seqMat) <- assay(seqMat) - assay(RefMat)
            }
            yLabText <- "logFC"
          } else if (input$clusterFor == "two-condition expression") {
            seqMat <- processedData()[,processedData()$treatment %in% c(input$seleTreat_cluster,input$seleTreat_clusterRef) & processedData()$timepoint %in% input$seleTimeRange]
            yLabText <- "Normalized expression"
          }
          # Time points are treated as numerical values and not characters (this might change in the future)
          # if both h and min are used as time unit, convert the minute ones to h (this might be removed in the future)
          if ((any(str_ends(seqMat$timepoint, "h"))) & (any(str_ends(seqMat$timepoint, "min")))) {
            ifMinToH <- TRUE
          } else ifMinToH <- FALSE
          if (is.null(seqMat$subjectID)) {
            plotTab <- data.frame(time = seqMat$timepoint,
                                  value = assays(seqMat)[["Intensity"]][geneID,],
                                  treatment = as.character(seqMat$treatment))
            # convert the time column to a numerical variable and convert the minute time points ('min') to hour ('h')
            if (ifMinToH) {
              plotTab$time[str_ends(plotTab$time, "min")] <- 1/60 * as.numeric(gsub("min","",plotTab$time[str_ends(plotTab$time, "min")]))
            }
            plotTab$time <- as.numeric(gsub("h|min","", plotTab$time))
            p <- ggplot(plotTab, aes(x= time, y = value)) +
              geom_point(aes(color = treatment), size=3) 
          } else {
            plotTab <- data.frame(time = seqMat$timepoint,
                                  value = assays(seqMat)[["Intensity"]][geneID,], 
                                  subjectID = seqMat$subjectID,
                                  treatment = as.character(seqMat$treatment)) %>%
              mutate(patCondi = paste0(subjectID,"_",treatment))
            # convert the time column to a numerical variable and convert the minute time points ('min') to hour ('h')
            if (ifMinToH) {
              plotTab$time[str_ends(plotTab$time, "min")] <- 1/60 * as.numeric(gsub("min", "", plotTab$time[str_ends(plotTab$time, "min")]))
            }
            plotTab$time <- as.numeric(gsub("h|min","", plotTab$time))
            p <- ggplot(plotTab, aes(x= time, y = value)) +
              geom_point(aes(color = treatment), size=3) + 
              geom_line(aes(group = patCondi), linetype = "dotted", color = "grey50")
          }
          p <- p + 
            stat_summary(aes(color=paste("mean",treatment)),fun = mean, geom = "line", linewidth = 2) +
            ylab(yLabText) + xlab("time") + 
            ggtitle(geneSymbol) + theme_bw() + 
            theme(text = element_text(size=15), plot.title = element_text(hjust = 0.5),
                  legend.position = "bottom",
                  axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5, size = 15))
          p
        }
      } else { NULL }
    } else { NULL }
  })
  
  # link for download the enriched set list
  output$downloadUI1 <- renderUI({
    if( !is.null (GSEres$resTab)) {
      downloadLink("downloadSet", "Download enriched set list")
    }
  })
  
  output$downloadUI2 <- renderUI({
    if (clickRecord$enrich) {
      downloadLink("downloadGene", "Download gene/phosphosite list")
    }
  })
  
  # a link to download enrichment results as csv
  output$downloadSet <- downloadHandler(
    filename = function() { paste('geneSetList', '.tsv', sep='') },
    content = function(file) {
      write.csv2(GSEres$resTab, file)
    }
  )
  
  # a link to download genes in the clicked set as csv
  output$downloadGene <- downloadHandler(
    filename = function() { paste('geneList', '.tsv', sep='') },
    content = function(file) {
      write.csv2(gseaList(), file)
    }
  )
  
  ####### Not in use: Plot the enrichment map for enrichment analysis################
  
  output$plot4 <- renderPlot({
    if (!is.null(GSEres$resObj)) {
      output$errMsg2 <- renderText("")
      networkPlot(GSEres$resObj, class="distinct", direction = input$setDirection, 
                  significance = input$pCut, adjusted = input$ifAdjCut, overlap = input$setOverlap, 
                  lay = as.numeric(input$layMethod), label = input$labelMethod, cexLabel = input$sizeLable,
                  ncharLabel = 100)
      
    } else output$errMsg2 <- renderText("Please run enrichment analysis first! (Enrichment on time series clusters not supported yet)")
  })
  
  ############################################ Kinase activity inference #############################################
  
  ####Widgets
  
  # reactive variable to store the decoupler network (update itself if change organism)
  decoupler_network <- reactive({
    network <- getDecouplerNetwork(input$speciesRef)
    network
  })
  # reactive variable to store the kinase activity result
  kinaseRes <- reactiveVal()
  
  # Perform the analysis when click the runKinase button
  runKinaseAnalysis <- observeEvent(input$runKinase, {
    if (input$assay == "Phosphoproteome") {
      if (input$seleSourceKinase == "Differential expression") {
        if (!is.null(filterDE())) {
          withProgress(message = "Running kinase activity inference, please wait...", {
            output$errMsgKinase <- renderText("")
            # compute the kinase score or report error (usually because the selected organism was not correct)
            tryCatch({
              scoreTab <- calcKinaseScore(filterDE(), decoupler_network(), statType = input$statTypeKinase, nPerm = input$nPermKinase)
              scoreTab <- scoreTab %>% mutate(padj = p.adjust(p_value, method = "BH")) %>% arrange(p_value)
              kinaseRes(scoreTab)
            }, error = function(e) {
              showModal(modalDialog(
                title = "No kinase found...",
                "Try again with a different organism or a larger phosphosite list.",
                easyClose = TRUE,
                footer = NULL
              ))})
          })
        } else output$errMsgKinase <- renderText("Please do a hypothesis testing first")
      } else {
        
        # Perform the analysis for time-series cluster
        if (!is.null(clusterTabVal())) {
          withProgress(message = "Running kinase activity inference, please wait...",  {
            output$errMsgKinase <- renderText("")
            # Add a column showing phosphosites based on the ID
            clusterData <- clusterTabVal() 
            clusterData <- clusterData[clusterData$cluster == input$seleCluster,]
            allClusterFeature <- clusterData %>% distinct(feature, .keep_all = TRUE) %>% .$feature
            allClusterSite <- data.frame(rowData(processedData())[allClusterFeature, "site"])
            allClusterSite$feature <- allClusterFeature
            clusterData <- clusterData %>% 
              left_join(allClusterSite, by = "feature") %>%
              rename(site = "rowData.processedData....allClusterFeature...site..")
            if (input$seleKinaseTimeMethod == "activity") {
              # compute the kinase ACTIVITY score
              # initiate an empty dataframe to store score result
              scoreTab <- data.frame(source = c(), score = c(), p_value = c(), timepoint = c())
              # get the order of time points
              timeVector <- input$seleTimeRange
              timeUnit <- suppressWarnings(str_extract(timeVector, "h|min"))
              timeUnit <- ifelse(is.na(timeUnit), "", timeUnit)
              # If both h and min are present, divide the min time points by 60
              if ((any(timeUnit == "h")) & (any(timeUnit == "min"))) {
                timeValue <- timeVector
                timeValue[timeUnit == "min"] <- 1/60 * as.numeric(gsub("min", "", timeValue[timeUnit == "min"]))
                timeRank <- rank(as.numeric(gsub("h", "", timeValue)))
              } else {
                timeRank <- rank(as.numeric(gsub("h|min", "", timeVector)))
              }
              timeRank <- as.character(timeRank)
              # try computing the kinase score. If error occurs (usually due to wrong organism selected) then inform the user to try again
              tryCatch({
                if (input$clusterFor == "expression") {
                  # Handling for expression case
                  siteTab <- processedData()[rowData(processedData())$site %in% unique(clusterData$site),
                                             processedData()$treatment == input$seleTreat_cluster & processedData()$timepoint %in% input$seleTimeRange]
                  # an empty vector to store order of timepoints in the heatmap
                  timeOrder = c()
                  # Compute the fold change of a time point with respect to the previous time point
                  for (i in 2:length(timeRank)) {
                    time2 <- timeVector[timeRank == as.character(i)]
                    time1 <- timeVector[timeRank == as.character(i-1)]
                    siteTime1 <- siteTab[,siteTab$timepoint == time1]
                    siteTime2 <- siteTab[,siteTab$timepoint == time2]
                    # match samples by subjectID if provided
                    if (!is.null(siteTab$subjectID)) {
                      siteTime2 <- siteTime2[,match(siteTime1$subjectID,siteTime2$subjectID)]
                      fc <- assay(siteTime2) - assay(siteTime1)
                      clusterDataTime <- data.frame(site = as.character(rowData(siteTab)$site),
                                                    log2FC = rowMeans(fc))
                    } else {
                      fc <- rowMeans(assay(siteTime2)) - rowMeans(assay(siteTime1))
                      clusterDataTime <- data.frame(site = as.character(rowData(siteTab)$site),
                                                    log2FC = fc)
                    }
                    clusterDataTime <- na.omit(clusterDataTime)
                    scoreTabTime <- calcKinaseScore(clusterDataTime, decoupler_network(),statType = "log2FC", nPerm = input$nPermKinase)
                    scoreTabTime$timepoint <- paste0(time2,"_",time1)
                    timeOrder <- append(timeOrder, paste0(time2,"_",time1))
                    scoreTab <- rbind(scoreTab, scoreTabTime)
                  }
                  # make sure the order of timepoints on the heatmap will be correct
                  scoreTab$timepoint <- factor(scoreTab$timepoint, levels = timeOrder)
                } else if (input$clusterFor == "logFC") {
                  # Handling for logFC case
                  for (time in input$seleTimeRange) {
                    clusterDataTime <- clusterData[clusterData$time == time,] %>%
                      select(value, site) %>% rename(log2FC = "value")
                    scoreTabTime <- calcKinaseScore(clusterDataTime, decoupler_network(), statType = "log2FC", nPerm = input$nPermKinase)
                    scoreTabTime$timepoint <- time
                    scoreTab <- rbind(scoreTab, scoreTabTime)
                  }
                  # make sure the time points are in correct order in the heatmap
                  scoreTab$timepoint <- factor(scoreTab$timepoint, levels = timeVector[order(match(timeRank, sort(timeRank)))])
                } else if (input$clusterFor == "two-condition expression") {
                  # Handling for two-condition expression case
                  siteTab <- processedData()[rowData(processedData())$site %in% unique(clusterData$site),
                                             processedData()$treatment == input$seleTreat_cluster & processedData()$timepoint %in% input$seleTimeRange]
                  refTab <- processedData()[rowData(processedData())$site %in% unique(clusterData$site),
                                            processedData()$treatment == input$seleTreat_clusterRef & processedData()$timepoint %in% input$seleTimeRange]
                  # compute the logFC for each time point
                  for (time in input$seleTimeRange) {
                    # if subjectID is present then match samples by subject ID
                    siteTabTime <- siteTab[,siteTab$timepoint == time]
                    refTabTime <- refTab[,refTab$timepoint == time]
                    if ((!is.null(siteTab$subjectID)) & (!is.null(refTab$subjectID))) {
                      siteTabTime <- siteTabTime[,match(refTabTime$subjectID, siteTabTime$subjectID)]
                      fc <- assay(siteTabTime) - assay(refTabTime)
                      clusterDataTime <- data.frame(site = as.character(rowData(siteTabTime)$site),
                                                    log2FC = rowMeans(fc))
                    } else{
                      fc <- rowMeans(assay(siteTabTime)) - rowMeans(assay(refTabTime))
                      clusterDataTime <- data.frame(site = as.character(rowData(siteTabTime)$site),
                                                    log2FC = fc)
                    }
                    clusterDataTime <- na.omit(clusterDataTime)
                    scoreTabTime <- calcKinaseScore(clusterDataTime, decoupler_network(), statType = "log2FC", nPerm = input$nPermKinase)
                    scoreTabTime$timepoint <- time
                    scoreTab <- rbind(scoreTab, scoreTabTime)
                  }
                  # make sure the time points are in correct order in the heatmap
                  scoreTab$timepoint <- factor(scoreTab$timepoint, levels = timeVector[order(match(timeRank, sort(timeRank)))])
                }
                scoreTab <- scoreTab %>% mutate(padj = p.adjust(p_value, method = "BH")) %>%
                  arrange(p_value)
                kinaseRes(scoreTab)
              }, error = function(e) {
                kinaseRes(NULL)
                showModal(modalDialog(
                  title = "No kinase found...",
                  "Try again with a different organism or a larger phosphosite list.",
                  easyClose = TRUE,
                  footer = NULL
                ))})
            } else {
              # compute how likely the kinases are associated with the cluster
              # An error is induced if no kinase is found (wrong organism or list too small)
              tryCatch({
              if (input$seleAssoMethod == "Fisher's exact test") {
                # Fisher's exact test to test kinase association with cluster
                pSiteCluster  <- unique(selectedCluster()$site)
                pSiteRefList <- unique(rowData(processedData())$site)
                pSiteRefList <- pSiteRefList[!pSiteRefList %in% pSiteCluster]
                rtab <- lapply(unique(decoupler_network()$source), function(kinase) {
                  kinaseTarget = as.character(decoupler_network()[decoupler_network()$source == kinase,"target"])
                  RinSet = sum(pSiteRefList %in% kinaseTarget)
                  RninSet = length(pSiteRefList) - RinSet
                  GinSet = sum(pSiteCluster %in% kinaseTarget)
                  GninSet = length(pSiteCluster) - GinSet
                  fmat = matrix(c(GinSet, RinSet, GninSet, RninSet), nrow = 2,
                                ncol = 2, byrow = F)
                  colnames(fmat) = c("inSet", "ninSet")
                  rownames(fmat) = c("genes", "reference")
                  fish = fisher.test(fmat, alternative = "greater")
                  pval = fish$p.value
                  inSet = RinSet + GinSet
                  tibble(source = kinase,
                         `number.pSite.in.cluster`= GinSet, 
                         `number.pSite.by.kinase` = inSet, 
                         p_value = pval)
                }) %>% bind_rows() %>%
                  filter(number.pSite.in.cluster>0) 
              } else { # i.e. if seleAssoMethod == 'FGSEA'
                # GSEA to test kinase association with cluster using FGSEA from decoupleR
                # Data is taken from clusterTabVal() instead of selectedCluster() since
                # the former has unrounded probability values
                selectedTab <- filter(clusterTabVal(), cluster == input$seleCluster) %>% 
                  distinct(feature, .keep_all = TRUE) %>%
                  mutate(feature = as.character(feature)) 
                clusterData <- rowData(processedData()[selectedTab$feature,])
                inputTab <- selectedTab %>%
                  mutate(site = clusterData$site) %>%
                  select(site, prob) %>%
                  column_to_rownames(var = "site") %>%
                  arrange(desc(prob))
                rtab <- decoupleR::run_fgsea(mat = inputTab,
                                               network = decoupler_network(),
                                               minsize = 1,
                                             times = input$nPermKinase) %>%
                  filter(statistic == "fgsea") %>%
                  select(-condition, -statistic) %>%
                  rename(enrich.score = "score")
              }
            # adjusting p-values and filtering based on (adjusted) p-values
            rtab <- rtab %>% 
              mutate(padj = p.adjust(p_value, method = "BH")) %>%
              arrange(p_value)
            if (input$ifKinaseFDR)
              rtab <- rtab %>% filter(padj <= input$pKinase) else
                rtab <- rtab %>% filter(p_value <= input$pKinase)
            # if no kinase: induce an error to show the pop-up window
            if (nrow(rtab) == 0)
              stop("No kinase found... perhaps the wrong organism was chosen or the p-value threshold was too low") else
                kinaseRes(as.data.frame(rtab))
            }, error = function(e) {
                kinaseRes(NULL)
                showModal(modalDialog(
                  title = "No kinase found...",
                  "Try again with a different organism or a larger phosphosite list.",
                  easyClose = TRUE,
                  footer = NULL
                ))})
              }
          })
        } else output$errMsgKinase <- renderText("Please do a time-series clustering (logFC) first")
      }}  else output$errMsgKinase <- renderText("This feature only support phosphoproteome data. Please double check the Preprocessing step!")
  })
  # output to display plot object
  output$plotKinase <- renderPlot({
    if (!is.null(kinaseRes())) {
      scoreTab <- kinaseRes()
      if (input$ifKinaseFDR)
        scoreTab$p_value <- scoreTab$padj
      if (input$seleSourceKinase == "Differential expression") {
        plot <- plotKinaseDE(scoreTab, nTop = input$nTopKinase, pCut = input$pKinase)
        plot
      } else {
        plot <- plotKinaseTimeSeries(scoreTab, pCut = input$pKinase, clusterName = input$seleCluster)
        plot
      }
    }
  })
  
  # output to display table of kinase score
  output$kinaseTab <- DT::renderDataTable({
    if (!is.null(kinaseRes())) {
      if ((input$seleSourceKinase == "Differential expression") | 
          ((input$seleSourceKinase == "Time-series cluster") & 
           (input$seleKinaseTimeMethod == "activity"))){
        # table for kinase activity
        resTab <- kinaseRes() %>% rename(Kinase = "source", Activity = "score")
        datatable(resTab, selection = 'single', rownames = FALSE,
                  caption = "Kinase activity") %>%
          formatStyle('Activity',background=styleInterval(c(0),c("lightblue","pink"))) %>%
          formatRound(c("Activity", "p_value", "padj"), digits = 3)
      } else { # i.e. if doing kinase association for time-series cluster
        # table for kinase association
        resTab <- kinaseRes() %>% rename(Kinase = "source")
        if (input$seleAssoMethod == "Fisher's exact test") {
        datatable(resTab, selection = 'single', rownames = FALSE,
                  caption = paste0("Kinases associated to ", input$seleCluster, ", Fisher's exact test")) %>%
          formatRound(c("p_value", "padj"), digits = 3)
        } else { # i.e. if use FGSEA to estimate kinase association
          datatable(resTab, selection = 'single', rownames = FALSE,
                    caption = paste0("Kinases associated to ", input$seleCluster, ", FGSEA")) %>%
            formatStyle('enrich.score',background=styleInterval(c(0),c("lightblue","pink"))) %>%
            formatRound(c("enrich.score","p_value", "padj"), digits = 3)
        }
      }
    }
  })
  
  # list of phosphosites targeted by the selected kinase
  pSiteList <- reactive({
    kinase <- as.character(kinaseRes()[as.integer(input$kinaseTab_row_last_clicked), "source"])
    pSite <- as.character(decoupler_network()[decoupler_network()$source == kinase, "target"])
    if (input$seleSourceKinase == "Differential expression") 
      pSiteTab <- filterDE()[filterDE()$site %in% pSite,] else 
        pSiteTab <- selectedCluster()[selectedCluster()$site %in% pSite,]
    pSiteTab <- pSiteTab %>% select(any_of(c("Gene","site","log2FC","stat","pvalue","padj","cluster","probability","Sequence")))
    pSiteTab
  })
  # get phosphosites that are targets of the selected kinase in kinaseTab
  observeEvent(input$kinaseTab_row_last_clicked, {
    clickRecord$kinase <- TRUE
  })
  
  # table to show level of phosphosite corresponding to the selected kinase
  output$pSiteTab <- DT::renderDataTable({
    if (clickRecord$kinase) {
      tryCatch({
        datatable(pSiteList(),
                  selection = 'single', rownames = FALSE, 
                  caption = "Phosphosites targeted by the selected kinase") %>% 
          formatRound(c("log2FC", "stat", "pvalue",  "padj"), digits = 3)
      }, error = function(e) {
        datatable(pSiteList(),
                  selection = 'single', rownames = FALSE, 
                  caption = "Phosphosites targeted by the selected kinase") 
      })
    }
  })
  #link for download the kinase score result
  output$downloadUI3 <- renderUI({
    if( !is.null (kinaseRes())) {
      downloadLink("downloadKinase","Download kinase activity inference result")
    }
  })
  #a link to download enrichment results as csv
  output$downloadKinase <- downloadHandler(
    filename = function() { paste('kinaseActivity', '.tsv', sep='') },
    content = function(file) {
      write.csv2(kinaseRes(), file)
    }
  )
  
  #link for download the phosphosite result
  output$downloadUI4 <- renderUI({
    if( !is.null(pSiteList())) {
      downloadLink("downloadPhosphoSite","Download phosphopeptide list")
    }
  })
  #a link to download phosphopeptide result as tsv
  output$downloadPhosphoSite <- downloadHandler(
    filename = function() { paste('phosphopeptide', '.tsv', sep='') },
    content = function(file) {
      write.csv2(pSiteList(), file)
    }
  )
  
})
