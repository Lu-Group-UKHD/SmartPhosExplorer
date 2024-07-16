# helper functions

# get the last symbol of a protein that has multiple gene symbols
getOneSymbol <- function(Gene) {
  outStr <- sapply(Gene, function(x) {
    sp <- str_split(x, ";")[[1]]
    sp[length(sp)]
  })
  names(outStr) <- NULL
  outStr
}


#' @name plotMissing 
#' 
#' @title Plot Missing Data Completeness
#'
#' @description 
#' `plotMissing` generates a bar plot showing the completeness (percentage of non-missing values) for each sample in a SummarizedExperiment object.
#'
#' @param se A SummarizedExperiment object containing the assay data.
#'
#' @return A ggplot object showing the percentage of completeness for each sample.
#'
#' @details
#' This function calculates the percentage of non-missing values for each sample in the provided SummarizedExperiment object. It then generates a bar plot where each bar represents a sample, and the height of the bar corresponds to the completeness (percentage of non-missing values) of that sample.
#'
#' @importFrom SummarizedExperiment assay
#' @importFrom ggplot2 ggplot aes geom_bar ggtitle ylab theme element_text
#' @importFrom tibble tibble
#' @examples
#' # Assuming 'se' is a SummarizedExperiment object:
#' plot <- plotMissing(se)
#' print(plot)
#'
#' @export
plotMissing <- function(se) {
  
  # Extract the assay data from the SummarizedExperiment object
  countMat <- assay(se)
  
  # Create a table with sample names and their corresponding percentage of non-missing values
  plotTab <- tibble(sample = se$sample, 
                    perNA = colSums(is.na(countMat))/nrow(countMat))
  
  # Generate the bar plot using ggplot2
   missPlot <- ggplot(plotTab, aes(x = sample, y = 1-perNA)) +
    geom_bar(stat = "identity") +
    ggtitle("Percentage of sample completeness") +
    ylab("completeness") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0),
          plot.title = element_text(hjust = 0.5, face = "bold"))
   
   return(missPlot)
}


#' @name plotIntensity
#' 
#' @title Plot Intensity Boxplots
#'
#' @description
#' `plotIntensity` generates boxplots of assay intensities for each sample in a SummarizedExperiment object. Optionally, the boxplots can be colored based on a specified metadata column. The function handles missing values by filtering them out before plotting.
#'
#' @param se A SummarizedExperiment object containing the assay data and metadata.
#' @param color A character string specifying the metadata column to use for coloring the boxplots. Default is "none".
#'
#' @return A ggplot object showing boxplots of intensities for each sample.
#'
#' @importFrom SummarizedExperiment assay colData
#' @importFrom ggplot2 ggplot aes geom_boxplot ggtitle theme element_text
#' @importFrom dplyr filter left_join
#' @importFrom tidyr pivot_longer
#' @importFrom tibble as_tibble
#' @examples
#' # Assuming 'se' is a SummarizedExperiment object:
#' plot <- plotIntensity(se, color = "group")
#' print(plot)
#'
#' @export
plotIntensity <- function(se, color = "none") {
  
  # Extract the assay data from the SummarizedExperiment object
  countMat <- assay(se)
  
  # Convert the assay data to a tibble, pivot to long format, and filter out missing values
  countTab <- countMat %>% as_tibble(rownames = "id") %>% 
    pivot_longer(-id) %>%
    filter(!is.na(value))
  
  # Extract metadata from the SummarizedExperiment object
  meta <- as.data.frame(colData(se))
  # Join the count data with metadata
  countTabmeta <- left_join(countTab, meta, by = c('name' = 'sample'))
  
  # Create the ggplot object with boxplots of intensities
  g <- ggplot(countTabmeta, aes(x = name, y = value)) +
    ggtitle("Boxplot of intensities") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0),
          plot.title = element_text(hjust = 0.5, face = "bold"))
  
  # Add color to the boxplots if a valid metadata column is specified
  if (color == "none"){
    g <- g + geom_boxplot()
  }
  else {
    g <- g + geom_boxplot(aes_string(fill = color))
  }
  
  return(g)
}

#' @name plotPCA
#' 
#' @title Plot PCA
#'
#' @description
#' `plotPCA` generates a PCA plot using the results from a PCA analysis and a SummarizedExperiment object. The points on the plot can be colored and shaped based on metadata.
#'
#' @param pca A PCA result object, typically obtained from \code{prcomp}.
#' @param se A SummarizedExperiment object containing the metadata.
#' @param xaxis A character string specifying which principal component to use for the x-axis. Default is `"PC1"`.
#' @param yaxis A character string specifying which principal component to use for the y-axis. Default is `"PC2"`.
#' @param color A character string specifying the metadata column to use for coloring the points. Default is `"none"`.
#' @param shape A character string specifying the metadata column to use for shaping the points. Default is `"none"`.
#'
#' @return A ggplot object showing the PCA plot.
#'
#' @details
#' This function creates a PCA plot using the scores from a PCA result object and metadata from a SummarizedExperiment object. The x-axis and y-axis can be customized to display different principal components, and the points can be optionally colored and shaped based on specified metadata columns.
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom ggplot2 ggplot aes geom_point theme_bw theme labs scale_shape
#' @importFrom dplyr left_join
#' @importFrom tibble rownames_to_column
#' @importFrom rlang sym
#' @examples
#' # Assuming 'pca' is a PCA result object and 'se' is a SummarizedExperiment object:
#' plot <- plotPCA(pca, se, xaxis = "PC1", yaxis = "PC2", color = "group", shape = "type")
#' print(plot)
#'
#' @export
plotPCA <- function(pca, se, xaxis = "PC1", yaxis = "PC2", color = "none", shape = "none") {
  
  # Calculate the proportion of variance explained by each principal component
  varExplained <- pca$sdev^2/sum(pca$sdev^2)
  # Convert the PCA result to a data frame
  pcaDf <- as.data.frame(pca[["x"]])
  # Convert the metadata to a data frame
  meta <- as.data.frame(colData(se))
  # Join the PCA scores with the metadata 
  pcaMeta <- left_join(rownames_to_column(pcaDf),
                       meta, by = c("rowname" = "sample"))
  
  # Create the initial ggplot object with labels for variance explained
  g <- ggplot(pcaMeta, aes(x = !!sym(xaxis), y = !!sym(yaxis),
                           text = paste("sample:", meta$sample))) +
    theme_bw() +
    theme(legend.position="top") +
    labs(x=paste0(xaxis,": ",
                  round(varExplained[as.numeric(strsplit(xaxis, "PC")[[1]][2])]*100, 1), "%"),
         y=paste0(yaxis,": ",
                  round(varExplained[as.numeric(strsplit(yaxis, "PC")[[1]][2])]*100, 1), "%")) +
    scale_shape(solid = FALSE)
  
  # Add points to the plot with optional color and shape aesthetics
  if (color == "none" & shape == "none") {
    g <- g + geom_point(size = 2)
  }
  else if (color == "none") {
    g <- g + geom_point(aes_string(shape = shape), size = 2)
  }
  else if (shape == "none") {
    g <- g + geom_point(aes_string(color = color), size = 2)
  }
  else {
    g <- g + geom_point(aes_string(color = color,
                                   shape = shape),
                        size = 2)
  }
}

#' @name plotHeatmap
#' 
#' @title Plot Heatmap of Intensity assay
#'
#' @description
#' `plotHeatmap` generates a heatmap for intensity assay for different conditions, including top variants, differentially expressed genes, and selected time series clusters.
#'
#' @param type A character string specifying the type of heatmap to plot. Options are `"Top variant"`, `"Differentially expressed"`, and `"Selected time series cluster"`.
#' @param se A SummarizedExperiment object containing the imputed intensity assay.
#' @param data An optional data frame containing additional data for `"Differentially expressed"` and `"Selected time series cluster"` types. Default is `NULL`.
#' @param top An integer specifying the number of top variants to plot. Default is `100`.
#' @param cutCol An integer specifying the number of clusters for columns. Default is `1`.
#' @param cutRow An integer specifying the number of clusters for rows. Default is `1`.
#' @param clustCol A logical value indicating whether to cluster columns. Default is `TRUE`.
#' @param clustRow A logical value indicating whether to cluster rows. Default is `TRUE`.
#' @param annotationCol A character vector specifying the columns in the metadata to use for annotation. Default is `NULL`.
#' @param title A character string specifying the title of the heatmap.
#'
#' @return A pheatmap object showing the heatmap of Intensity data.
#'
#' @details
#' This function creates a heatmap using the Intensity assay from a SummarizedExperiment object. The heatmap can show the top variants based on standard deviation, differentially expressed genes, or selected time series clusters. Row normalization is performed, and the heatmap can include annotations based on specified metadata columns.
#'
#' @importFrom SummarizedExperiment assays colData rowData
#' @importFrom ggplot2 ggplot aes geom_point theme_bw theme labs scale_shape
#' @importFrom dplyr arrange left_join
#' @importFrom pheatmap pheatmap
#' @importFrom tibble rownames_to_column
#' @importFrom rlang sym
#' @examples
#' # Assuming 'se' is a SummarizedExperiment object and 'data' is a data frame:
#' heatmap <- plotHeatmap("Top variant", se, top = 50, annotationCol = "group", title = "Top Variants Heatmap")
#' print(heatmap)
#'
#' @export
plotHeatmap <- function(type, se, data = NULL, top = 100, cutCol = 1, cutRow = 1, clustCol = TRUE, clustRow = TRUE, annotationCol, title) {
  
  # Select the appropriate intensity assay and gene IDs based on the type of heatmap
  if (type == "Top variant") {
    exprMat <- assays(se)[["imputed"]]
    sds <- apply(exprMat, 1, sd)
    orderID <- names(sort(sds, decreasing = TRUE))
    geneIDs <- orderID[seq(1, as.integer(top))]
    exprMat <- exprMat[geneIDs,]
    geneSymbol <- rowData(se[match(geneIDs, rownames(se)),])$Gene
  } 
  else if (type == "Differentially expressed") {
    if(!is.null(data)) {
      geneIDs <- arrange(data, stat)$ID
      exprMat <- assays(se)[["imputed"]][geneIDs,] 
      geneSymbol <- data[match(geneIDs, data$ID),]$Gene
    }
    else {
      print("Please give data argument")
    }
  }
  else if (type == "Selected time series cluster") {
    if(!is.null(data)) {
      geneIDs <- unique(data$ID)
      exprMat <- assays(se)[["imputed"]][geneIDs,] 
      geneSymbol <- data[match(geneIDs, data$ID),]$Gene
    }
    else {
      print("Please give data argument")
    }
  }
  
  # Prepare column annotations from the metadata
  cd <- as.data.frame(colData(se))
  annCol <- cd[row.names(cd) %in% colnames(exprMat),][c(annotationCol)]
  row.names(annCol) <- colnames(exprMat)
  
  # Prepare color scale for the heatmap
  color <- colorRampPalette(c("navy", "white", "firebrick"))(100)
  
  # Perform row normalization and clip extreme values
  exprMat <- t(scale(t(exprMat)))
  exprMat[exprMat > 4] <- 4
  exprMat[exprMat < -4] <- -4
  
  # Plot the heatmap based on the type and whether annotations are provided
  if (type == "Top variant") {
    if (is.null(annotationCol)) {
      p <- pheatmap(exprMat, color = color,
                    labels_row = geneSymbol,
                    treeheight_row = 0, treeheight_col = 0,
                    main = title,
                    cutree_cols = cutCol,
                    cutree_rows = cutRow)
    }
    else {
      p <- pheatmap(exprMat, color = color,
                    labels_row = geneSymbol,
                    treeheight_row = 0, treeheight_col = 0,
                    main = title,
                    cutree_cols = cutCol,
                    cutree_rows = cutRow,
                    annotation_col = annCol)
    }
  }
  else {
    # Sort the columns by their names before plotting
    exprMat <- exprMat[, sort(colnames(exprMat))]
    if (is.null(annotationCol)) {
      p <- pheatmap(exprMat, color = color,
                    labels_row = geneSymbol,
                    treeheight_row = 0, treeheight_col = 0,
                    main = title,
                    cluster_rows = clustRow, cluster_cols = clustCol,
                    cutree_cols = cutCol,
                    cutree_rows = cutRow)
    }
    else {
      p <- pheatmap(exprMat, color = color,
                    labels_row = geneSymbol,
                    treeheight_row = 0, treeheight_col = 0,
                    main = title,
                    cluster_rows = clustRow, cluster_cols = clustCol,
                    cutree_cols = cutCol,
                    cutree_rows = cutRow,
                    annotation_col = annCol)
    }
  }
  return(p)
}


#' @name plotVolcano
#' 
#' @title Plot Volcano Plot for Differential Expression Analysis
#'
#' @description
#' `plotVolcano` generates a volcano plot to visualize differential expression results.
#'
#' @param tableDE A data frame containing differential expression results with columns 'ID', 'log2FC', 'pvalue', and 'Gene'.
#' @param pFilter A numeric value specifying the p-value threshold for significance. Default is `0.05`.
#' @param fcFilter A numeric value specifying the log2 fold-change threshold for significance. Default is `0.5`.
#'
#' @return A ggplot object representing the volcano plot.
#'
#' @details
#' This function creates a volcano plot where differentially expressed genes are categorized as 'Up', 'Down', or 'Not Sig' based on the provided p-value and log2 fold-change thresholds. Points on the plot are color-coded to indicate their expression status.
#'
#' @importFrom ggplot2 ggplot aes geom_vline geom_hline geom_point annotate scale_color_manual xlab ggtitle theme
#' @importFrom dplyr mutate case_when
#' @examples
#' # Assuming 'tableDE' is a data frame containing differential expression results:
#' volcanoPlot <- plotVolcano(tableDE, pFilter = 0.05, fcFilter = 0.5)
#' print(volcanoPlot)
#'
#' @export
plotVolcano <- function(tableDE, pFilter = 0.05, fcFilter = 0.5) {
  
  # Convert the input table to a data frame and ensure the 'ID' column is of type character
  dataVolcano <- data.frame(tableDE)
  dataVolcano$ID <- as.character(dataVolcano$ID)
  # Categorize each gene based on the provided p-value and log2 fold-change thresholds
  dataVolcano <- mutate(dataVolcano, expression = case_when(
    dataVolcano$log2FC >= as.numeric(fcFilter) & dataVolcano$pvalue <= as.numeric(pFilter) ~ "Up",
    dataVolcano$log2FC <= -as.numeric(fcFilter) & dataVolcano$pvalue <= as.numeric(pFilter) ~ "Down",
    dataVolcano$pvalue > as.numeric(pFilter) | (dataVolcano$log2FC < as.numeric(fcFilter) & dataVolcano$log2FC > -as.numeric(fcFilter)) ~ "Not Sig"
  ))
  
  # Create the volcano plot
  v <- ggplot(dataVolcano, aes(x = log2FC, y = -log10(pvalue), label = Gene, customdata = ID)) +
    # Add vertical lines for fold-change thresholds
    geom_vline(xintercept = 0, color = "black", linetype = "solid", size = 0.25) +
    geom_vline(xintercept = as.numeric(fcFilter), color = "darkgrey", linetype = "dashed") +
    geom_vline(xintercept = -as.numeric(fcFilter), color = "darkgrey", linetype = "dashed") +
    # Add horizontal lines for p-value thresholds
    geom_hline(yintercept = -log10(as.numeric(pFilter)), color = "darkgrey", linetype = "dashed") +
    annotate(x = 5.0, y = -log10(as.numeric(pFilter))-0.1, label = paste("P-value = ", as.numeric(pFilter)),
             geom = "text", size = 3, color = "darkgrey") +
    geom_hline(yintercept = -log10(0.25), color="darkgrey", linetype = "dashed") +
    annotate(x = 5.0, y = 0.5, label = paste("P-value = ", 0.25),
             geom = "text", size=3, color="darkgrey") +
    # Plot the points and color them based on their expression status
    geom_point(aes(color = expression), size = 0.9) +
    scale_color_manual(values = c("Up" = "firebrick3", "Down" = "navy", "Not Sig" = "darkgrey")) +
    xlab("absolute log2(Quantity) difference") +
    ggtitle("Volcano plot") +
    theme(plot.title = element_text(hjust=0.5))
  return(v)
}


########################################## Normalization Correction ######################################


#' @name medianNorm
#' 
#' @title Normalize a Matrix Using Median or Mean
#'
#' @description
#' `medianNorm` normalizes the columns of a matrix by either the median or the mean.
#'
#' @param x A numeric matrix to be normalized.
#' @param method A character string specifying the normalization method. Options are `"median"` or `"mean"`. Default is `"median"`.
#'
#' @return A numeric matrix with normalized columns.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item If the `method` is `"median"`, it calculates the median of each column and adjusts by the overall median of these medians.
#'   \item If the `method` is `"mean"`, it calculates the mean of each column and adjusts by the overall mean of these means.
#'   \item It constructs a matrix of these adjusted values and subtracts it from the original matrix to normalize the columns.
#' }
#'
#' @examples
#' # Example usage:
#' x <- matrix(rnorm(20), nrow=5, ncol=4)
#' normalized_x <- medianNorm(x, method = "median")
#' print(normalized_x)
#'
#' @importFrom matrixStats colMedians
#' @export
medianNorm <- function(x, method = "median") {
  if (method == "median") {
    # Calculate the median of each column, ignoring NA values
    mVal <- matrixStats::colMedians(x, na.rm = TRUE)
    # Adjust by the overall median of these medians
    mVal <- mVal - median(mVal, na.rm = TRUE)
  } else if (method == "mean") {
    # Calculate the mean of each column, ignoring NA values
    mVal <- colMeans(x, na.rm = TRUE)
    # Adjust by the overall mean of these means
    mVal <- mVal - mean(mVal, na.rm = TRUE)
  }
  mMat <- matrix(rep(mVal, each = nrow(x)), ncol =ncol(x))
  return(x-mMat)
}

# function for performing normalization of FP and PP samples

#' @name performCombinedNormalization
#' 
#' @title Perform Combined Normalization on MultiAssayExperiment Data
#'
#' @description
#' `performCombinedNormalization` performs combined normalization on proteome and phosphoproteome data from a MultiAssayExperiment object.
#'
#' @param maeData A MultiAssayExperiment object containing proteome and phosphoproteome data.
#'
#' @return A numeric matrix with normalized and log2-transformed data.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Extracts the count matrices for Full Proteome (FP) samples.
#'   \item Combines the proteome and phosphoproteome data into a single matrix.
#'   \item Removes rows with all NA values.
#'   \item Performs median normalization and log2 transformation on the combined matrix.
#' }
#'
#' @examples
#' # Example usage:
#' # Assuming maeData is a MultiAssayExperiment object with appropriate data
#' # normalized_data <- performCombinedNormalization(maeData)
#' # print(normalized_data)
#'
#' @importFrom MultiAssayExperiment assay
#' @export
performCombinedNormalization <- function(maeData) {
  
  # get count matrix from FP (Full Proteome) samples
  setFP <- maeData[,maeData$sampleType %in% c("FullProteome", "FP")]
  protFP <- assay(setFP[["Proteome"]])
  phosFP <- assay(setFP[["Phosphoproteome"]])
  # Combine proteome and phosphoproteome data into a single matrix
  comFP <- rbind(protFP, phosFP)
  comFP <- comFP[rowSums(!is.na(comFP))>0,]
  
  # perform median normalization and log2 transformation
  comFP.norm <- medianNorm(log2(comFP))
  
  return(comFP.norm)
}


#' @name getRatioMatrix
#' 
#' @title Get Ratio Matrix of Phosphoproteome Data
#'
#' @description
#' `getRatioMatrix` calculates the ratio matrix of phosphoproteome data from a MultiAssayExperiment object.
#'
#' @param maeData A MultiAssayExperiment object containing phosphoproteome and full proteome data.
#' @param normalization A logical value indicating whether to perform normalization. Default is `FALSE`.
#' @param getAdjustedPP A logical value indicating whether to use adjusted phosphoproteome data. Default is `FALSE`.
#'
#' @return A numeric matrix representing the ratio of intensity of PP (phosphoproteome) data to FP (full proteome) data.
#'
#' @examples
#' # Example usage:
#' # Assuming maeData is a MultiAssayExperiment object with appropriate data
#' # ratio_matrix <- getRatioMatrix(maeData, normalization = TRUE, getAdjustedPP = TRUE)
#' # print(ratio_matrix)
#'
#' @importFrom MultiAssayExperiment assay assays
#' @export
getRatioMatrix <- function(maeData, normalization = FALSE, getAdjustedPP = FALSE) {
  
  # Ensure the normalization parameter is logical
  stopifnot(is.logical(normalization))
  
  # Extract and log-transform phosphoproteome data based on getAdjustedPP flag
  if (!getAdjustedPP) {
    phosPP <- log2(assay(maeData[,maeData$sampleType %in% c("Phospho", "PP")][["Phosphoproteome"]]))
  } else {
    phosPP <- log2(assays(maeData[,maeData$sampleType %in% c("Phospho", "PP")][["Phosphoproteome"]])[["Intensity_adjusted"]])
  }
  
  # Extract and log-transform full proteome data, with optional normalization
  if (!normalization) {
    phosFP <- log2(assay(maeData[,maeData$sampleType %in% c("FullProteome", "FP")][["Phosphoproteome"]]))
  } 
  else  {
    phosFP <- performCombinedNormalization(maeData)
  }
  
  # Use the sample name without prefix as column name
  colnames(phosFP) <- maeData[,maeData$sampleType %in% c("Phospho", "PP")]$sampleName
  colnames(phosPP) <- maeData[,maeData$sampleType %in% c("FullProteome", "FP")]$sampleName
  
  # Calculate the ratio matrix
  allSmp <- intersect(colnames(phosFP), colnames(phosPP))
  allRow <- intersect(rownames(phosFP), rownames(phosPP))
  ratioMat <- phosPP[allRow, allSmp] - phosFP[allRow, allSmp]
  ratioMat <- ratioMat[rowSums(!is.na(ratioMat)) >0,]
  
  return(ratioMat)
  
}


# plot the log ration of PP/FP intensities


#' @name plotLogRatio
#' 
#' @title Plot Log Ratio of PP/FP (Phosphoproteome to Full Proteome) intensities
#'
#' @description
#' `plotLogRatio` generates a boxplot of the log2 ratio of intensities of phosphoproteome to full proteome data from a MultiAssayExperiment object.
#'
#' @param maeData A MultiAssayExperiment object containing phosphoproteome and full proteome data.
#' @param normalization A logical value indicating whether to perform normalization. Default is `FALSE`.
#'
#' @return A ggplot object representing the boxplot of the log2 ratios.
#'
#' @examples
#' # Example usage:
#' # Assuming maeData is a MultiAssayExperiment object with appropriate data
#' # plot <- plotLogRatio(maeData, normalization = TRUE)
#' # print(plot)
#'
#' @importFrom MultiAssayExperiment assay
#' @importFrom matrixStats colMedians
#' @importFrom tibble as_tibble
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr filter mutate
#' @importFrom ggplot2 ggplot aes geom_boxplot ggtitle xlab ylab geom_hline theme element_text
#' @export 
plotLogRatio <- function(maeData, normalization = FALSE) {
  # Calculate the ratio matrix of phosphoproteome to full proteome data
  ratioMat <- getRatioMatrix(maeData, normalization)
  # Extract and log-transform phosphoproteome data
  phosPP <- log2(assay(maeData[,maeData$sampleType %in% c("Phospho", "PP")][["Phosphoproteome"]]))
  medianPP <- colMedians(phosPP,na.rm = TRUE)
  names(medianPP) <- maeData[,maeData$sampleType %in% c("Phospho", "PP")]$sampleName
  
  # Create a table for plotting
  plotTab <- as_tibble(ratioMat, rownames = "feature") %>%
    pivot_longer(-feature) %>%
    filter(!is.na(value)) %>%
    mutate(medianPP = medianPP[name])
  
  # Generate a ggplot boxplot of the log2 ratios
  ggplot(plotTab, aes(x=name, y=value)) +
    geom_boxplot(aes(fill = medianPP)) +
    ggtitle("Boxplot of Phospho/FullProteome Ratio") +
    xlab("sample") +
    ylab("log2(ratio)") +
    geom_hline(yintercept = median(median(plotTab$value, na.rm=TRUE)), linetype = "dashed", color = "red") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          plot.title = element_text(hjust = 0.5, face = "bold"))
  
}

# check the PP/FP ratio matrix and remove feature that do not meet requirements

#' @name checkRatioMat
#' 
#' @title Check the PP/FP ratio matrix and remove feature that do not meet requirements
#'
#' @description
#' `checkRatioMat` checks the ratio matrix for samples that do not have sufficient overlap of phospho-peptides between enriched (PP) and unenriched (FP) samples.
#'
#' @param ratioMat A numeric matrix representing the ratio of phosphoproteome data to full proteome data.
#' @param minOverlap An integer specifying the minimum number of overlapping peptides required between samples. Default is `3`.
#'
#' @return A character vector of sample names that do not meet the overlap criteria.
#'
#' @examples
#' # Example usage:
#' # Assuming ratioMat is a matrix with appropriate data
#' # excluded_samples <- checkRatioMat(ratioMat, minOverlap = 3)
#' # print(excluded_samples)
#'
#' @export
checkRatioMat <- function(ratioMat, minOverlap = 3) {
  # Initialize a list to keep track of excluded samples
  excludeSampleList <- c()
  
  # Identify samples that don't have any phospho sites detect in both FP and PP samples
  noOverSmp <- colnames(ratioMat)[colSums(!is.na(ratioMat))==0]
  if (length(noOverSmp) >0) {
    warning(paste0("Below samples don't have phopho-peptides detected in both enriched (PP) and unenriched (FP) samples and therefore adjusting factor will set to 0 (no adjustment) for them:\n",
                   paste0(noOverSmp, collapse = ", ")))
    excludeSampleList <- c(excludeSampleList, noOverSmp)
  }
  
  # Remove the identified samples from the ratio matrix
  ratioMat <- ratioMat[, !colnames(ratioMat) %in% noOverSmp]
  
  # Check for samples that do not have enough peptide overlap with other samples
  pairOverlap <- sapply(colnames(ratioMat), function(n) {
    subMat <- ratioMat[!is.na(ratioMat[,n]),]
    minOver <- min(colSums(!is.na(subMat)))
  })
  
  tooFewOverlap <- colnames(ratioMat)[pairOverlap < minOverlap]
  
  if (length(tooFewOverlap) >0) {
    warning(paste0("Below samples don't enough number of overlapped phopho-peptides with other samples and therefore adjusting factor will set to 0 (no adjustment) for them:\n",
                   paste0(tooFewOverlap, collapse = ", ")))
    excludeSampleList <- c(excludeSampleList, tooFewOverlap)
  }
  
  return(excludeSampleList)
}



#' @name runPhosphoAdjustment
#' 
#' @title Run Phospho Adjustment
#'
#' @description
#' `runPhosphoAdjustment` performs phospho adjustment on a MultiAssayExperiment object to normalize the phosphoproteome data.
#'
#' @param maeData A MultiAssayExperiment object containing phosphoproteome and full proteome data.
#' @param normalization A logical value indicating whether to perform normalization. Default is `FALSE`.
#' @param minOverlap An integer specifying the minimum number of overlapping peptides required between samples. Default is `3`.
#' @param completeness A numeric value indicating the required completeness of data for features to be included. Default is `0`.
#' @param ncore An integer specifying the number of cores to use for parallel processing. Default is `1`.
#'
#' @return A MultiAssayExperiment object with adjusted phosphoproteome data.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Defines an optimization function to minimize the sum of squared differences between pairs of samples.
#'   \item Calculates the ratio matrix of phosphoproteome to full proteome data.
#'   \item Subsets features based on completeness criteria.
#'   \item Performs a sanity check to identify and exclude problematic samples.
#'   \item Sets initial values for the adjustment factor based on column medians.
#'   \item Estimates the adjustment factor using parallel optimization.
#'   \item Adjusts the phosphoproteome measurements using the estimated adjustment factor.
#' }
#'
#' @examples
#' # Example usage:
#' # Assuming maeData is a MultiAssayExperiment object with appropriate data
#' # adjusted_maeData <- runPhosphoAdjustment(maeData, normalization = TRUE, minOverlap = 3, completeness = 0.8, ncore = 2)
#' # print(adjusted_maeData)
#'
#' @importFrom MultiAssayExperiment assay assays colData
#' @importFrom matrixStats colMedians
#' @importFrom stats optim
#' @importFrom parallel makeCluster setDefaultCluster stopCluster
#' @importFrom utils combn
#' @export
runPhosphoAdjustment <- function(maeData, normalization = FALSE, minOverlap = 3, completeness = 0, ncore = 1 ) {
  
  # Function to opitmize
  esFun <- function(par, data) {
    comPair <- utils::combn(seq(length(par)), 2)
    sum(((data[comPair[1, ],] + par[comPair[1, ]]) - (data[comPair[2, ],] + par[comPair[2, ]]))^2/rowSums(!is.na(data[comPair[1,],] + data[comPair[2,],])), na.rm = TRUE)
  }
  
  # Get PP/FP ratio matrix
  ratioMat <- getRatioMatrix(maeData, normalization = normalization)
  adjFac <- structure(rep(0, length.out = ncol(ratioMat)), names = colnames(ratioMat))
  
  # Subset features according to completeness in the ratio matrix
  ratioMat <- ratioMat[rowSums(!is.na(ratioMat))/ncol(ratioMat) >= completeness,]
  
  # Sanity check to see if any samples need to be excluded
  excList <- checkRatioMat(ratioMat)
  ratioMat <- ratioMat[, !colnames(ratioMat) %in% excList]
  
  # Set an initial value for B based on col medians of ratioMat, may increase search speed
  colMed <- apply(ratioMat,2, median, na.rm = TRUE)
  iniPar <- median(colMed) - colMed
  
  # Estimating adjusting factor
  cl <- makeCluster(ncore)
  setDefaultCluster(cl = cl)
  optRes <- optimParallel(par=iniPar, fn=esFun, data=t(ratioMat))
  stopCluster(cl)
  
  # Add adjusting factor to sample annotation  
  adjFac[names(optRes$par)] <- optRes$par
  ppName <- colnames(maeData[,maeData$sampleType %in% c("Phospho", "PP")][["Phosphoproteome"]])
  adjFac <- structure(adjFac[maeData[,ppName]$sampleName], names = ppName)
  maeData$adjustFactorPP <- unname(adjFac[match(rownames(colData(maeData)),names(adjFac))])
  # Adjust phospho measurement on PP samples
  phosMat <- assay(maeData[,maeData$sampleType %in% c("Phospho", "PP")][["Phosphoproteome"]])
  phosMat <- t(t(phosMat)*(2^adjFac))
  assays(maeData[["Phosphoproteome"]])[["Intensity_adjusted"]] <- assays(maeData[["Phosphoproteome"]])[["Intensity"]]
  assays(maeData[["Phosphoproteome"]])[["Intensity_adjusted"]][,colnames(phosMat)] <- phosMat
  
  return(maeData)
}


#' @name plotAdjustmentResults
#' 
#' @title Plot Adjustment Results
#'
#' @description
#' `plotAdjustmentResults` generates plots to visualize the results of phosphoproteome adjustment.
#'
#' @param maeData A MultiAssayExperiment object containing phosphoproteome and full proteome data.
#' @param normalization A logical value indicating whether normalization was performed. Default is `FALSE`.
#'
#' @return A list containing:
#' \item{ratioTrendPlot}{A ggplot object showing the line plot of PP/FP ratios for features present in all samples.}
#' \item{ratioBoxplot}{A ggplot object showing the box plot of PP/FP ratios before and after adjustment.}
#' \item{ppBoxplot}{A ggplot object showing the box plot of phosphorylation intensities in PP samples before and after adjustment.}
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Checks if the adjustment factor is present in the sample annotation.
#'   \item Calculates the ratio matrix before and after adjustment.
#'   \item Creates a trend line plot for features present in all samples.
#'   \item Creates box plots of the PP/FP ratios and phosphorylation intensities before and after adjustment.
#' }
#'
#' @examples
#' # Example usage:
#' # Assuming maeData is a MultiAssayExperiment object with appropriate data
#' # plots <- plotAdjustmentResults(maeData, normalization = TRUE)
#' # print(plots$ratioTrendPlot)
#' # print(plots$ratioBoxplot)
#' # print(plots$ppBoxplot)
#'
#' @importFrom MultiAssayExperiment assay assays colData
#' @importFrom dplyr bind_rows filter group_by mutate summarise
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 aes facet_wrap geom_boxplot geom_hline geom_line geom_point ggplot ggtitle theme xlab ylab
#' @export
plotAdjustmentResults <- function(maeData, normalization = FALSE) {
  # Check if the adjustment factor has been applied
  if (!"adjustFactorPP" %in% colnames(colData(maeData))) {
    stop("Phosphorylation measurments have not been adjusted yet. Please perform normalization adjustment using calcAdjustFacotr function first")
  }
  
  # Visualize precursors before and after adjustment
  ratioMat.ori <- getRatioMatrix(maeData, normalization = normalization, getAdjustedPP = FALSE)
  ratioMat.adj <- getRatioMatrix(maeData, normalization = normalization, getAdjustedPP = TRUE)
  ratioPlotTab <- bind_rows(pivot_longer(as_tibble(ratioMat.ori, rownames = "id"), -id, names_to = "sample", values_to = "ratio") %>% mutate(adjustment = "before adjustment"),
                            pivot_longer(as_tibble(ratioMat.adj, rownames = "id"), -id, names_to = "sample", values_to = "ratio") %>% mutate(adjustment = "after adjustment")) %>%
    mutate(adjustment = factor(adjustment, levels = c("before adjustment","after adjustment"))) %>%
    filter(!is.na(ratio))
  
  # For precursors present in all samples
  featureComplete <- rownames(ratioMat.ori)[complete.cases(ratioMat.ori)]
  if (!length(featureComplete) >0) {
    warning("No feature (PP/FP ratio) has been detected in all samples. Raio trend line plot of will not be generated")
    ratioTrendPlot <- NULL 
  } else {
    ratioPlotTab.complete <- filter(ratioPlotTab, id %in% featureComplete)
    ratioTrendPlot <- ggplot(ratioPlotTab.complete, aes(x=sample, y=ratio, color = adjustment)) +
      geom_line(aes(group =id), linetype = "dashed") +
      geom_point() +
      facet_wrap(~adjustment, ncol=1) +
      ggtitle("Line plot of PP/FP ratios available for all samples") + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            legend.position = "none") +
      xlab("") + ylab("log2(PP/FP) ratio") 
  }
  
  # For ratio box plots
  medTab <- group_by(ratioPlotTab, adjustment) %>%
    summarise(medVal = median(ratio, na.rm=TRUE))
  ratioBoxplot <- ggplot(ratioPlotTab, aes(x=sample, y=ratio, fill = adjustment)) +
    geom_boxplot() +
    geom_hline(data=medTab, aes(yintercept = medVal), linetype = "dashed", color = "red") +
    facet_wrap(~adjustment, ncol=1) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          legend.position = "none") +
    xlab("") + ylab("log2(PP/FP) ratio") +
    ggtitle("Box plot of all PP/FP ratios")  
  
  
  # For phosphorylation measurements of PP samples before and after adjustment
  ppMat.adj <- assays(maeData[,maeData$sampleType %in% c("Phospho", "PP")][["Phosphoproteome"]])[["Intensity_adjusted"]] 
  ppMat.ori <- assays(maeData[,maeData$sampleType %in% c("Phospho", "PP")][["Phosphoproteome"]])[["Intensity"]] 
  ppPlotTab <- bind_rows(pivot_longer(as_tibble(ppMat.ori, rownames = "id"), -id, names_to = "sample", values_to = "count") %>% mutate(adjustment = "before adjustment"),
                         pivot_longer(as_tibble(ppMat.adj, rownames = "id"), -id, names_to = "sample", values_to = "count") %>% mutate(adjustment = "after adjustment")) %>%
    mutate(adjustment = factor(adjustment, levels = c("before adjustment","after adjustment")),
           count = log2(count)) %>%
    filter(!is.na(count))
  
  # For phosphorylation intensity box plots
  medTab <- group_by(ppPlotTab, adjustment) %>%
    summarise(medVal = median(count, na.rm=TRUE))
  ppBoxplot <- ggplot(ppPlotTab, aes(x=sample, y=count, fill = adjustment)) +
    geom_boxplot() +
    geom_hline(data = medTab, aes(yintercept = medVal), linetype = "dashed", color = "red") +
    facet_wrap(~adjustment, ncol=1) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          legend.position = "none") +
    xlab("") + ylab("log2(intensity)") +
    ggtitle("Box plot of phosphorylation in PP samples")
  
  return(list(ratioTrendPlot = ratioTrendPlot,
              ratioBoxplot = ratioBoxplot,
              ppBoxplot = ppBoxplot))
}


#' @name preprocessProteome
#' 
#' @title Preprocess Proteome Data
#'
#' @description
#' `preprocessProteome` preprocesses proteome data stored in a SummarizedExperiment object by performing filtering, transformation, normalization, imputation, and batch effect removal.
#'
#' @param seData A SummarizedExperiment object containing proteome data.
#' @param filterList A list of filters to apply on the samples. Default is `NULL`.
#' @param missCut Numeric value specifying the missing value cutoff percentage for filtering features. Default is `50`.
#' @param transform Character string specifying the transformation method ("log2", "vst", "none"). Default is `"log2"`.
#' @param normalize Logical value indicating whether to normalize the data. Default is `FALSE`.
#' @param getPP Logical value indicating whether to retrieve PP samples. Default is `FALSE`.
#' @param removeOutlier Character vector of samples to be removed as outliers. Default is `NULL`.
#' @param impute Character string specifying the imputation method ("QRILC", "MLE", "bpca", "missForest", "MinDet", "none"). Default is `"QRILC"`.
#' @param batch Character vector specifying batch effects to remove. Default is `NULL`.
#' @param verbose Logical value indicating whether to print detailed information. Default is `FALSE`.
#' @param scaleFactorTab Data frame containing scale factors for normalization. Default is `NULL`.
#'
#' @return A SummarizedExperiment object with preprocessed proteome data.
#'
#' @examples
#' # Example usage:
#' # Assuming seData is a SummarizedExperiment object with appropriate data
#' # processedData <- preprocessProteome(seData, filterList = list(sampleType = "FP"), missCut = 30, normalize = TRUE)
#'
#' @importFrom SummarizedExperiment colData rowData assay assays
#' @importFrom dplyr filter mutate
#' @importFrom tidyr pivot_longer
#' @importFrom DEP impute
#' @importFrom limma removeBatchEffect
#' @importFrom missForest missForest
#' @importFrom doParallel registerDoParallel
#' @importFrom doRNG registerDoRNG
#' @export
preprocessProteome <- function(seData, filterList = NULL, missCut = 50,
                               transform = "log2", normalize = FALSE, getPP = FALSE,
                               removeOutlier = NULL, impute = "QRILC", batch = NULL,
                               verbose = FALSE, scaleFactorTab = NULL) {
  
  # Retrieve desired sample type
  if (getPP) {
    # Retrieve PP samples if specified
    fpe <- seData[,seData$sampleType %in% c("Phospho", "PP")]
    colData(fpe) <- colData(seData)[colnames(fpe),]
  } else {
    # Otherwise, retrieve FullProteome samples
    fpe <- seData[,seData$sampleType %in% c("FullProteome", "FP")]
    colData(fpe) <- colData(seData)[colnames(fpe),]
  }
  
  # Remove specified outliers
  if (length(removeOutlier) > 0) {
    if (length(removeOutlier) > 1) {
      for (i in removeOutlier) {
        fpe <- fpe[, !grepl(i, fpe$sample)]
      }
    }
    else {
      fpe <- fpe[, !grepl(removeOutlier, fpe$sample)]
    }
  }
  
  # Apply specified filters
  if (!is.null(filterList)) {
    for (n in names(filterList)) {
      fpe <- fpe[,fpe[[n]] %in% filterList[[n]]]
    }
  }
  
  # Rename columns to sample names
  colnames(fpe) <- fpe$sample
  
  # Process gene names and remove features without symbols
  rowData(fpe)$Gene <- getOneSymbol(rowData(fpe)$Gene)
  fpe <- fpe[!rowData(fpe)$Gene %in% c(NA,""),]
  
  # Filter features based on missing values
  countMat <- assay(fpe)
  missPer <- rowSums(is.na(countMat))/ncol(countMat)*100
  fpeSub <- fpe[missPer < missCut,]
  
  # Apply transformation and normalization
  if (transform=="log2") {
    if (normalize) {
      if (is.null(scaleFactorTab)) {
        assay(fpeSub) <- medianNorm(log2(assay(fpeSub)))
      } else {
        assay(fpeSub) <- log2(t(t(assay(fpeSub))/scaleFactorTab[match(paste0(fpeSub$sample),scaleFactorTab$sample),]$scaleFactor))
      }
    } else {
      assay(fpeSub) <- log2(assay(fpeSub))
    }
  } else if (transform == "vst") {
    if (normalize) {
      if (is.null(scaleFactorTab)) {
        assay(fpeSub) <- vsn::justvsn(assay(fpeSub))
      } else {
        normMat <- t(t(assay(fpeSub))/scaleFactorTab[match(paste0(fpeSub$sample),scaleFactorTab$sample),]$scaleFactor)
        assay(fpeSub) <- vsn::justvsn(normMat, calib="none")
      }
    } else {
      assay(fpeSub) <- vsn::justvsn(assay(fpeSub), calib="none")
    }
  } else if (transform == "none") {
    if (normalize) {
      if (is.null(scaleFactorTab)) {
        assay(fpeSub) <- medianNorm(assay(fpeSub))
      } else {
        assay(fpeSub) <- t(t(assay(fpeSub))/scaleFactorTab[match(paste0(fpeSub$sample),scaleFactorTab$sample),]$scaleFactor)
      }
    } else {
      assay(fpeSub) <- assay(fpeSub)
    }
  }
  
  # Impute missing values
  if (impute != "none") {
    rowData(fpeSub)$name <- rowData(fpeSub)$UniprotID
    rowData(fpeSub)$ID <- rowData(fpeSub)$UniprotID
    if (impute == "QRILC") {
      imp <- DEP::impute(fpeSub, "QRILC")
    }
    else if (impute == "MLE") {
      imp <- DEP::impute(fpeSub, "MLE")
    }
    else if (impute == "bpca") {
      imp <- DEP::impute(fpeSub, "bpca")
    }
    else if (impute == "missForest") {
      doParallel::registerDoParallel(cores = 6)  # set based on number of CPU cores
      doRNG::registerDoRNG(seed = 123)
      mf <- missForest::missForest(t(assay(fpeSub)), parallelize = "forests", maxiter = 2, ntree = 50)
      imp <- t(mf$ximp)
    }
    else if (impute == "MinDet") {
        imp <- DEP::impute(fpeSub, "MinDet")  
    }
    assays(fpeSub)[["imputed"]] <- assay(imp)
    rowData(fpeSub)$name <- NULL
    rowData(fpeSub)$ID <- NULL
    if (verbose) {
      # Show number of samples and features
      print("Number of proteins and samples:")
      print(dim(fpeSub))
      }
  }
  
  # Remove batch effects if specified
  if(!is.null(batch)) {
    if(length(batch) == 1) {
      remBatchImp <- limma::removeBatchEffect(assays(fpeSub)[["imputed"]],
                                              batch = colData(fpeSub)[,batch])
      remBatch <- limma::removeBatchEffect(assay(fpeSub),
                                           batch = colData(fpeSub)[,batch])
    }
    else {
      remBatchImp <- limma::removeBatchEffect(assays(fpeSub)[["imputed"]],
                                              batch = colData(fpeSub)[,batch[1]],
                                              batch2 = colData(fpeSub)[,batch[2]])
      remBatch <- limma::removeBatchEffect(assay(fpeSub),
                                           batch = colData(fpeSub)[,batch[1]],
                                           batch2 = colData(fpeSub)[,batch[2]])
    }
    assays(fpeSub)[["imputed"]] <- assay(remBatchImp)
    assay(fpeSub) <- assay(remBatch)
  }
  
  return(fpeSub)
}

#' @name preprocessPhos 
#' 
#' @title Preprocess Phosphoproteome Data
#'
#' @description
#' `preprocessPhos` preprocesses phosphoproteome data stored in a SummarizedExperiment object by performing filtering, transformation, normalization, imputation, and batch effect removal.
#'
#' @param seData A SummarizedExperiment object containing phosphoproteome data.
#' @param filterList A list of filters to apply on the samples. Default is `NULL`.
#' @param missCut Numeric value specifying the missing value cutoff percentage for filtering features. Default is `50`.
#' @param transform Character string specifying the transformation method ("log2", "vst", "none"). Default is `"log2"`.
#' @param normalize Logical value indicating whether to normalize the data. Default is `FALSE`.
#' @param getFP Logical value indicating whether to retrieve FP samples. Default is `FALSE`.
#' @param removeOutlier Character vector of samples to be removed as outliers. Default is `NULL`.
#' @param assayName Character string specifying the assay name in the SummarizedExperiment object. Default is `NULL`.
#' @param batch Character vector specifying batch effects to remove. Default is `NULL`.
#' @param scaleFactorTab Data frame containing scale factors for normalization. Default is `NULL`.
#' @param impute Character string specifying the imputation method ("QRILC", "MLE", "bpca", "missForest", "MinDet", "none"). Default is `"QRILC"`.
#' @param verbose Logical value indicating whether to print detailed information. Default is `FALSE`.
#'
#' @return A SummarizedExperiment object with preprocessed phosphoproteome data.
#'
#' @examples
#' # Example usage:
#' # Assuming seData is a SummarizedExperiment object with appropriate data
#' # processedData <- preprocessPhos(seData, filterList = list(sampleType = "Phospho"), missCut = 30, normalize = TRUE)
#'
#' @importFrom SummarizedExperiment colData rowData assay assays
#' @importFrom dplyr filter mutate
#' @importFrom tidyr pivot_longer
#' @importFrom DEP impute
#' @importFrom limma removeBatchEffect
#' @importFrom missForest missForest
#' @importFrom doParallel registerDoParallel
#' @importFrom doRNG registerDoRNG
#' @export
preprocessPhos <- function(seData, filterList = NULL, missCut = 50,
                           transform="log2", normalize = FALSE, getFP = FALSE,
                           removeOutlier = NULL, assayName = NULL, batch = NULL,
                           scaleFactorTab = NULL, impute = "QRILC", verbose = FALSE) {
  
  # Retrieve the desired sample type or specified assay
  if (is.null(assayName)) {
    if (getFP) { 
      # Retrieve FullProteome samples if specified
      ppe <- seData[,seData$sampleType %in% c("FullProteome", "FP")]
      colData(ppe) <- colData(seData)[colnames(ppe),]
    } else {
      # Otherwise, retrieve Phospho samples
      ppe <- seData[,seData$sampleType %in% c("Phospho", "PP")]
      colData(ppe) <- colData(seData)[colnames(ppe),]
    }
  } else {
    ppe <- seData[[assayName]]
    colData(ppe) <- colData(seData[,colnames(ppe)])
  }
  
  # Remove specified outliers
  if (length(removeOutlier) > 0) {
    if (length(removeOutlier) > 1) {
      for (i in removeOutlier) {
        ppe <- ppe[, !grepl(i, ppe$sample)]
      }
    }
    else {
      ppe <- ppe[, !grepl(removeOutlier, ppe$sample)]
    }
  }
  
  # Apply specified filters
  if (!is.null(filterList)) {
    for (n in names(filterList)) {
      ppe <- ppe[,ppe[[n]] %in% filterList[[n]]]
    }
  }
  
  # Rename columns to sample names
  colnames(ppe) <- ppe$sample
  # Get last gene name
  rowData(ppe)$Gene <- getOneSymbol(rowData(ppe)$Gene)
  # Get last phosphorylation site
  rowData(ppe)$Residue <- getOneSymbol(rowData(ppe)$Residue)
  rowData(ppe)$Position <- getOneSymbol(rowData(ppe)$Position)
  # Remove features without gene symbols
  ppe <- ppe[!rowData(ppe)$Gene %in% c(NA,""),]
  # Rename phosphorylation sites
  rowData(ppe)$site <- paste0(rowData(ppe)$Gene,"_",rowData(ppe)$Residue,rowData(ppe)$Position)
  # Filter features based on missing values
  countMat <- assay(ppe)
  missPer <- rowSums(is.na(countMat))/ncol(countMat)*100
  ppeSub <- ppe[missPer < missCut,]
  
  # Apply transformation and normalization
  if (transform=="log2") {
    if (normalize) {
      if (is.null(scaleFactorTab)) {
        assay(ppeSub) <- medianNorm(log2(assay(ppeSub)))
      } else {
        assay(ppeSub) <- log2(t(t(assay(ppeSub))/scaleFactorTab[match(paste0(ppeSub$sample),scaleFactorTab$sample),]$scaleFactor))
      }
    } else {
      assay(ppeSub) <- log2(assay(ppeSub))
    }
  } else if (transform == "vst") {
    if (normalize) {
      if (is.null(scaleFactorTab)) {
        assay(ppeSub) <- vsn::justvsn(assay(ppeSub))
      } else {
        normMat <- t(t(assay(ppeSub))/scaleFactorTab[match(paste0(ppeSub$sample),scaleFactorTab$sample),]$scaleFactor)
        assay(ppeSub) <- vsn::justvsn(normMat, calib="none")
      }
    } else {
      assay(ppeSub) <- vsn::justvsn(assay(ppeSub), calib="none")
    }
  } else if (transform == "none") {
    if (normalize) {
      if (is.null(scaleFactorTab)) {
        assay(ppeSub) <- medianNorm(assay(ppeSub))
      } else {
        assay(ppeSub) <- t(t(assay(ppeSub))/scaleFactorTab[match(paste0(ppeSub$sample),scaleFactorTab$sample),]$scaleFactor)
      }
    } else {
      assay(ppeSub) <- assay(ppeSub)
    }
  }
  
  # Impute missing values
  if (impute != "none") {
    rowData(ppeSub)$name <- rowData(ppeSub)$site
    rowData(ppeSub)$ID <- rowData(ppeSub)$site
    if (impute == "QRILC") {
      imp <- DEP::impute(ppeSub, "QRILC")
    }
    else if (impute == "MLE") {
      imp <- DEP::impute(ppeSub, "MLE")
    }
    else if (impute == "bpca") {
      imp <- DEP::impute(ppeSub, "bpca")
    }
    else if (impute == "MinDet") {
        imp <- DEP::impute(ppeSub, "MinDet")  
    }
    else {
      doParallel::registerDoParallel(cores = 6)  # Set number of CPU cores
      doRNG::registerDoRNG(seed = 123)
      mf <- missForest::missForest(t(assay(ppeSub)), parallelize = "forests", maxiter = 2, ntree = 50)
      imp <- t(mf$ximp)
    }
    assays(ppeSub)[["imputed"]] <- assay(imp)
    rowData(ppeSub)$name <- NULL
    rowData(ppeSub)$ID <- NULL
    # Show number of samples and features
    if (verbose) {
      print("Number of proteins and samples:")
      print(dim(ppeSub))
    }
  }
  
  # Remove batch effects if specified
  if(!is.null(batch)) {
    if(length(batch) == 1) {
      remBatchImp <- limma::removeBatchEffect(assays(ppeSub)[["imputed"]],
                                              batch = colData(ppeSub)[,batch])
      remBatch <- limma::removeBatchEffect(assay(ppeSub),
                                           batch = colData(ppeSub)[,batch])
    }
    else {
      remBatchImp <- limma::removeBatchEffect(assays(ppeSub)[["imputed"]],
                                              batch = colData(ppeSub)[,batch[1]],
                                              batch2 = colData(ppeSub)[,batch[2]])
      remBatch <- limma::removeBatchEffect(assay(ppeSub),
                                           batch = colData(ppeSub)[,batch[1]],
                                           batch2 = colData(ppeSub)[,batch[2]])
    }
    assays(ppeSub)[["imputed"]] <- assay(remBatchImp)
    assay(ppeSub) <- assay(remBatch)
  }
  
  return(ppeSub)
}

#' @name addZeroTime
#' 
#' @title Add Zero Timepoint Data to Treatment Subset
#'
#' @description
#' `addZeroTime` adds a zero timepoint to a specific treatment's data subset.
#'
#' @param data A SummarizedExperiment object containing the experimental data.
#' @param treat Character string specifying the treatment to which zero timepoint should be added.
#' @param zeroTreat Character string specifying the treatment representing the zero timepoint.
#' @param timeRange Character vector specifying the timepoints to include for the treatment.
#'
#' @return A SummarizedExperiment object with the zero timepoint added to the specified treatment's data.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Subsets the data for the specified treatment and time range.
#'   \item Subsets the data for the zero timepoint of the specified zero treatment.
#'   \item Combines the assays from the treatment and zero timepoint subsets.
#'   \item Updates the column data to reflect the combined treatment.
#'   \item Returns a SummarizedExperiment object with the combined data.
#' }
#'
#' @examples
#' # Example usage:
#' # Assuming data is a SummarizedExperiment object with appropriate data
#' # result <- addZeroTime(data, treat = "TreatmentA", zeroTreat = "Control", timeRange = c("10min", "20min"))
#'
#' @importFrom SummarizedExperiment colData rowData assay assays elementMetadata SummarizedExperiment
#' @export 
addZeroTime <- function(data, treat, zeroTreat, timeRange) {
  # Subset the data for the specified treatment and time range
  subset1 <- data[, data$treatment == treat & data$timepoint %in% timeRange]
  # Subset the data for the zero timepoint of the specified zero treatment
  subset2 <- data[, data$treatment == zeroTreat & data$timepoint %in% c("0min", "0", "0h")]
  # Combine the assays from the treatment and zero timepoint subsets
  assay <- cbind(assay(subset1), assay(subset2))
  colnames(assay) <- gsub(zeroTreat, treat, colnames(assay))
  
  # Combine the column data from both subsets
  cd1 <- colData(subset1)
  cd2 <- colData(subset2)
  cd <- rbind(cd1, cd2)
  cd$treatment[cd$treatment == zeroTreat] = treat
  cd$sample[cd$sample == zeroTreat] = treat
  rownames(cd) <- gsub(zeroTreat, treat, rownames(cd))
  
  # Retrieve the element metadata from the original data
  emeta <- elementMetadata(data)
  
  return(SummarizedExperiment(assays=SimpleList(intensity=assay), colData = cd, rowData = emeta))
}


#' @name clusterTS
#' 
#' @title Perform Clustering on Time-Series Data
#'
#' @description
#' `clusterTS` performs clustering on time-series data and generates plots for visualization.
#'
#' @param x A numeric matrix with rows as features and columns as time points.
#' @param k An integer specifying the number of clusters.
#' @param pCut A numeric value specifying the probability cutoff for cluster membership. Default is `NULL`.
#' @param twoCondition A logical value indicating if the data contains two conditions. Default is `FALSE`.
#'
#' @return A list containing:
#' \item{cluster}{A tibble with clustering information for each feature.}
#' \item{plot}{A ggplot object for visualizing the clustering results.}
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Sets a seed for reproducibility.
#'   \item Removes rows with missing values.
#'   \item Performs clustering using fuzzy C-means.
#'   \item Filters clusters based on the probability cutoff if provided.
#'   \item Generates plots for visualizing clustering results.
#' }
#'
#' @examples
#' # Example usage:
#' # Assuming x is a numeric matrix with time-series data
#' # result <- clusterTS(x, k = 4, pCut = 0.8, twoCondition = FALSE)
#'
#' @importFrom e1071 cmeans
#' @importFrom dplyr filter left_join mutate arrange group_by ungroup
#' @importFrom tidyr pivot_longer separate
#' @importFrom ggplot2 ggplot aes geom_line scale_color_gradient facet_wrap theme_bw theme element_text unit
#' @importFrom tibble as_tibble tibble
#' @importFrom stringr str_extract str_split
#' @importFrom magrittr %>%
#' @export
clusterTS <- function(x, k, pCut = NULL, twoCondition = FALSE) {
  
  # Set seed for reproducible clustering results
  set.seed(12345)
  
  # Remove rows with NA values
  x.center <- x[complete.cases(x),]
  
  # Perform fuzzy C-means clustering
  res <- e1071::cmeans(x.center, k)

  # Create a tibble with clustering results
  resCluster <- tibble(feature = names(res$cluster),
                       cluster = res$cluster,
                       prob = rowMax(res$membership))
  
  # Filter clusters based on probability cutoff if provided
  if (!is.null(pCut)) resCluster <- filter(resCluster, prob >= pCut)

  if (!twoCondition) {
    # Handle single condition data
    
    # Extract unique time points and determine time unit
    timeVector <- unique(colnames(x.center))
    timeUnit <- str_extract(timeVector, "h|min")
    timeUnit <- ifelse(is.na(timeUnit), "", timeUnit)
    
    # Adjust time values if both hours and minutes are present
    if ((any(timeUnit == "h")) & (any(timeUnit == "min"))) {
      timeValue <- timeVector
      timeValue[timeUnit == "min"] <- 1/60 * as.numeric(gsub("min", "", timeValue[timeUnit == "min"]))
      timeRank <- rank(as.numeric(gsub("h", "", timeValue)))
    } else {
      timeRank <- rank(as.numeric(gsub("h|min", "", timeVector)))
    }
    
    # Determine the order of time points for plotting
    timeOrder <- timeVector[order(match(timeRank, sort(timeRank)))]  
    
    # Prepare data for plotting
    clusterTab <- x.center %>% as_tibble(rownames = "feature") %>%
      pivot_longer(-feature, names_to = "time", values_to = "value") %>%
      left_join(resCluster, by = "feature") %>% filter(!is.na(cluster)) %>%
      arrange(prob) %>%
      mutate(cluster = paste0("cluster",cluster),
             feature = factor(feature, levels = unique(feature))) %>%
      group_by(cluster) %>%
      mutate(cNum = length(unique(feature))) %>%
      mutate(clusterNum = sprintf("%s (%s)",cluster, cNum)) %>%
      ungroup() %>%
      mutate(time = factor(time, levels = timeOrder)) %>%
      mutate(time = droplevels(time))

    # Generate the plot
    p <- ggplot(clusterTab, aes( x = time, y = value, group = feature)) +
      geom_line( aes(col = prob), alpha=0.8) +
      scale_color_gradient(low = "navy", high = "green") +
      facet_wrap(~clusterNum,ncol=3, scales = "free") +
      theme_bw() +
      theme(legend.position = "bottom",
            axis.title = element_text(size = 15),
            axis.text = element_text(size = 12),
            legend.key.width = unit(3,"cm"),
            legend.key.height = unit(0.8, "cm"),
            legend.title = element_text(size = 15, face = "bold"),
            legend.text = element_text(size = 15),
            strip.text = element_text(size=15, face = "bold"))
  } else {
    # Handle data with two conditions
    
    # Extract unique time points and determine time unit
    timeVector <- sapply(colnames(x.center), function(X) unlist(str_split(X, "_"))[1])
    timeVector <- unique(timeVector)
    timeUnit <- str_extract(timeVector, "h|min")
    timeUnit <- ifelse(is.na(timeUnit), "", timeUnit)
    
    # Adjust time values if both hours and minutes are present
    if ((any(timeUnit == "h")) & (any(timeUnit == "min"))) {
      timeValue <- timeVector
      timeValue[timeUnit == "min"] <- 1/60 * as.numeric(gsub("min", "", timeValue[timeUnit == "min"]))
      timeRank <- rank(as.numeric(gsub("h", "", timeValue)))
    } else {
      timeRank <- rank(as.numeric(gsub("h|min", "", timeVector)))
    }
    
    # Determine the order of time points for plotting
    timeOrder <- timeVector[order(match(timeRank, sort(timeRank)))]
    
    # Prepare data for plotting
    clusterTab <- x.center %>% as_tibble(rownames = "feature") %>%
      pivot_longer(-feature, names_to = "timeTreat", values_to = "value") %>%
      left_join(resCluster, by = "feature") %>% filter(!is.na(cluster)) %>%
      arrange(prob) %>%
      mutate(cluster = paste0("cluster",cluster),
             feature = factor(feature, levels = unique(feature))) %>%
      group_by(cluster) %>%
      mutate(cNum = length(unique(feature))) %>%
      mutate(clusterNum = sprintf("%s (%s)",cluster, cNum)) %>%
      ungroup() %>%
      separate(timeTreat,c("time","treatment"),"_", extra = "merge") %>%
      mutate(time = factor(time, levels = timeOrder)) %>%
      mutate(time = droplevels(time)) %>%
      mutate(geneGroup = paste0(feature,treatment))
    
    # Generate the plot
    p <- ggplot(clusterTab, aes( x=time, y= value, group = geneGroup)) +
      geom_line( aes(alpha = prob, color= treatment)) +
      facet_wrap(~clusterNum,ncol=3, scales = "free") +
      theme_bw() +
      theme(legend.position = "bottom",
            axis.title = element_text(size = 15),
            axis.text = element_text(size = 12),
            legend.key.size = unit(1,"cm"),
            legend.title = element_text(size = 15, face = "bold"),
            legend.text = element_text(size = 15),
            strip.text = element_text(size=15, face = "bold"))
  }
  
  return(list(cluster = clusterTab, plot = p))
}


# Function to filter genes based on spline fitting to remove changes that are not
# consistent. If patient IDs are available (recognized by "subjectID" in fileTable),
# then the changes among patients are taken into account, otherwise the
# replicates are considered independent.


#' @name splineFilter
#' 
#' @title Filter Expression Matrix Using Spline Models
#'
#' @description
#' `splineFilter` filters an expression matrix based on spline models fitted to time-series data, optionally considering treatment and subject ID.
#'
#' @param exprMat A numeric matrix of expression data, where rows are features and columns are samples.
#' @param subjectID An optional vector of subject IDs corresponding to columns in `exprMat`. Default is `NULL`.
#' @param time A numeric vector representing the time points corresponding to columns in `exprMat`.
#' @param df An integer specifying the degrees of freedom for the spline basis.
#' @param pCut A numeric value for the p-value cutoff to filter significant features. Default is `0.05`.
#' @param ifFDR A logical value indicating if the false discovery rate (FDR) should be used for filtering. If FALSE, raw p-values are used. Default is `FALSE`.
#' @param treatment An optional vector of treatment labels corresponding to columns in `exprMat`. Default is `NULL`.
#' @param refTreatment An optional reference treatment label for the `treatment` vector. Default is `NULL`.
#'
#' @return A filtered expression matrix containing only the features that meet the significance criteria.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Converts time points from minutes to hours if both units are present.
#'   \item Removes rows with missing values from the expression matrix.
#'   \item Constructs a design matrix for the spline model, optionally including subject IDs and treatments.
#'   \item Fits a linear model using the design matrix and performs empirical Bayes moderation.
#'   \item Extracts significant features based on the specified p-value or FDR cutoff.
#' }
#'
#' @examples
#' # Example usage:
#' # Assuming exprMat is a numeric matrix with expression data
#' # filteredMat <- splineFilter(exprMat, time = timeVec, df = 4, pCut = 0.05, ifFDR = TRUE)
#'
#' @importFrom limma lmFit eBayes topTable
#' @importFrom splines ns
#' @importFrom dplyr filter as_tibble
#' @importFrom stringr str_ends
#' @importFrom tibble rownames_to_column
#' @export
splineFilter <- function(exprMat, subjectID = NULL, time, df, pCut = 0.5, ifFDR = FALSE, treatment = NULL, refTreatment = NULL) {
  # Convert time points from minutes to hours if both units are present
  if ((all(str_ends(time,"h|min"))) & (!all(str_ends(time,"h"))) & (!all(str_ends(time,"min")))) {
    time[str_ends(time, "min")] <- 1/60 * as.numeric(gsub("min","", time[str_ends(time, "min")]))
  }
  time <- as.numeric(gsub("h|min", "", time))
  
  # Remove rows with NA values from the expression matrix
  exprMat <- exprMat[complete.cases(exprMat), ]
  if (is.null(treatment)) {
    # Handle case with one condition or log fold-change
    if (is.null(subjectID)) {
      # Create design matrix without subject IDs
      designTab <- data.frame(row.names = colnames(exprMat))
      designTab$X <- splines::ns(time, df)
      design <- model.matrix(~ 0 + X, data = designTab)
    } else { 
      # Create design matrix with subject IDs
      designTab <- data.frame(row.names = colnames(exprMat), subjectID = subjectID)
      designTab$X <- splines::ns(time, df)
      design <- model.matrix(~ 0 + X + subjectID, data = designTab)
    }
    
    # Fit linear model and perform empirical Bayes moderation
    fit <- lmFit(exprMat, design = design)
    fit2 <- eBayes(fit)
    # Extract results table and filter by p-value or FDR
    resTab <- topTable(fit2, coef = seq(df), number = Inf) %>%
      as_tibble(rownames = "ID") 
    
    if (ifFDR) resTab$p <- resTab$adj.P.Val else resTab$p <- resTab$P.Value
    
    resTab <- filter(resTab, p <= pCut)
    # Return filtered expression matrix
    
    return(exprMat[resTab$ID,])
  } else {
    # Handle case with two conditions
    if (is.null(subjectID)) {
      # Create design matrix without subject IDs
      designTab <- data.frame(row.names = colnames(exprMat), treatment = treatment)
      designTab$treatment <- factor(designTab$treatment, levels = unique(treatment))
      designTab$treatment <- relevel(designTab$treatment, refTreatment)
      designTab$treatment <- droplevels(designTab$treatment)
      designTab$X <- splines::ns(time, df)
      design <- model.matrix(~ 0 + X*treatment, data = designTab)
    } else {
      # Create design matrix with subject IDs
      designTab <- data.frame(row.names = colnames(exprMat), subjectID = subjectID, treatment = treatment)
      designTab$treatment <- factor(designTab$treatment, levels = unique(treatment))
      designTab$treatment <- relevel(designTab$treatment, refTreatment)
      designTab$treatment <- droplevels(designTab$treatment)
      designTab$X <- splines::ns(time, df)
      design <- model.matrix(~ 0 + subjectID + X*treatment, data = designTab)
    }
    
    # Fit linear model and perform empirical Bayes moderation
    fit <- lmFit(exprMat, design = design)
    fit2 <- eBayes(fit)
    
    # Extract results table and filter by p-value or FDR
    resTab <- topTable(fit2, coef = (ncol(design)-df+1):ncol(design), number = Inf) %>%
      as_tibble(rownames = "ID")
    
    if (ifFDR) resTab$p <- resTab$adj.P.Val else resTab$p <- resTab$P.Value
    
    resTab <- filter(resTab, p <= pCut)
    
    # Return filtered expression matrix
    return(exprMat[resTab$ID,])
  }
}

# scaling function for clustering

#' @name mscale
#' 
#' @title Scale and Center a Matrix
#'
#' @description
#' `mscale` scales and centers each row of a matrix, with options for using mean or median, standard deviation or mean absolute deviation, and censoring extreme values.
#'
#' @param x A numeric matrix where rows are features and columns are samples.
#' @param center Logical. If TRUE, the rows are centered by subtracting the mean or median. Default is `TRUE`.
#' @param scale Logical. If TRUE, the rows are scaled by dividing by the standard deviation or mean absolute deviation. Default is `TRUE`.
#' @param censor A numeric vector of length one or two for censoring the scaled values. 
#' If length one, values are censored symmetrically at positive and negative values. 
#' If length two, the first value is the lower limit and the second value is the upper limit. Default is `NULL`.
#' @param useMad Logical. If TRUE, the mean absolute deviation is used for scaling instead of the standard deviation. Default is `FALSE`.
#'
#' @return A scaled and centered numeric matrix with the same dimensions as the input matrix `x`.
#'
#' @details
#' The function allows for flexible scaling and centering of the rows of a matrix:
#' \itemize{
#'   \item If both `center` and `scale` are TRUE, rows are centered and scaled.
#'   \item If only `center` is TRUE, rows are centered but not scaled.
#'   \item If only `scale` is TRUE, rows are scaled but not centered.
#'   \item If neither `center` nor `scale` is TRUE, the original matrix is returned.
#' }
#' The function can also censor extreme values, either symmetrically or asymmetrically, based on the `censor` parameter.
#'
#' @examples
#' # Example usage:
#' # Assuming `dataMatrix` is a numeric matrix with expression data
#' # scaledMatrix <- mscale(dataMatrix, center = TRUE, scale = TRUE, censor = 2)
#'
#' @importFrom stats median sd
#' @importFrom matrixStats rowMads
#' @export
mscale <- function(x, center = TRUE, scale = TRUE, censor = NULL, useMad = FALSE){
  
  # Check if both scaling and centering are requested
  if (scale & center) {
    # Scale using Mean Absolute Deviation (MAD)
    if (useMad) {
      x.scaled <- apply(x, 1, function(y) (y-median(y,na.rm = T))/meanAD(y))
    } else {
      # Scale using Standard Deviation (SD)
      x.scaled <- apply(x, 1, function(y) (y-mean(y,na.rm=T))/sd(y,na.rm = T))
    }
  } else if (center & !scale) {
    # Only center the data
    if (useMad) {
      x.scaled <- apply(x, 1, function(y) (y-median(y,na.rm=T)))
    } else {
      x.scaled <- apply(x, 1, function(y) (y-mean(y,na.rm=T)))
    }
  } else if (!center & scale) {
    # Only scale the data
    if (useMad) {
      x.scaled <- apply(x, 1, function(y) y/meanAD(y))
    } else {
      x.scaled <- apply(x, 1, function(y) y/sd(y,na.rm = T))
    }
  } else {
    # Neither center nor scale
    x.scaled <- t(x)
  }
  
  # Apply censoring if requested
  if (!is.null(censor)) {
    if (length(censor) == 1) {
      # Symmetric censoring
      x.scaled[x.scaled > censor] <- censor
      x.scaled[x.scaled < -censor] <- -censor
    } else {
      # Asymmetric censoring
      x.scaled[x.scaled > censor[2]] <- censor[2]  # Upper limit
      x.scaled[x.scaled < censor[1]] <- censor[1]  # Lower limit
    }
  }
  return(t(as.matrix(x.scaled)))
}

# function to run fisher test for enrichment analysis (time series clustering only)
# Note: inputSet is the gene set database or the PTM set database
# ptm == TRUE Processing the database file will depend on whether the file is
# geneset or ptm set. If ptm == TRUE, the ptmset will be used, and each geneset
# will be split into 2 based on the sites direction of regulation (up or down).


#' @name runFisher
#' 
#' @title Perform Fisher's Exact Test on Gene Sets
#'
#' @description
#' `runFisher` performs Fisher's Exact Test to determine the enrichment of a set of genes within reference gene sets.
#'
#' @param genes A character vector of genes of interest.
#' @param reference A character vector of reference genes.
#' @param inputSet A list containing gene set collections. If `ptm` is TRUE, this should be a data frame with specific columns.
#' @param ptm Logical. If TRUE, perform the test on post-translational modification (PTM) gene sets. Default is `FALSE`.
#'
#' @return A data frame with the results of the Fisher's Exact Test, including the gene set name, the number of genes in the set, set size, p-value, adjusted p-value, and the genes in the set.
#'
#' @details
#' The function can operate in two modes: standard gene sets and PTM-specific gene sets. For PTM-specific gene sets, additional filtering and processing are performed.
#'
#' @examples
#' # Example usage:
#' # Assuming `genesOfInterest` is a vector of genes and `referenceGenes` is a vector of reference genes
#' # and `geneSetCollection` is a list of gene sets
#' # results <- runFisher(genesOfInterest, referenceGenes, geneSetCollection, ptm = FALSE)
#'
#' @importFrom dplyr filter group_by ungroup mutate bind_rows arrange tibble
#' @importFrom tidyr separate
#' @importFrom stats fisher.test p.adjust
#' @export
runFisher <- function (genes, reference, inputSet, ptm = FALSE) {
  
  # Retrieve the gene sets
  if (!ptm) {
    genesets <- inputSet$gsc
    setList <- 1:length(genesets)
  } else {
    # Filter and process the PTM-specific gene sets
    genesets <- inputSet %>%
      filter(!grepl("KINASE", category)) %>%
      dplyr::as_tibble() %>%
      filter(site.ptm == "p") %>%
      group_by(signature) %>%
      filter(n() >= 5) %>%
      ungroup() %>%
      mutate(signature = ifelse(site.direction == "u", paste0(signature,"_upregulated"), paste0(signature, "_downregulated"))) %>%
      separate(site.annotation, sep =  ":", into = c("site", "PubMedID"), extra="merge", fill="right") %>%
      as.data.frame()
    setList <- unique(genesets$signature)
  }
  # Remove genes from the reference set that are in the genes of interest
  reference = reference[!reference %in% genes]
  
  # Apply Fisher's Exact Test to each gene set
  rtab = lapply(setList, function(i) { 
    # Identify the geneset and its name
    if (!ptm) {
      geneset = genesets[[i]]
      nameSet = names(genesets)[i]
    } else {
      geneset = genesets[genesets$signature == i, "site"]
      nameSet = i
    }
    
    # Create the contingency table for Fisher's Exact Test
    RinSet = sum(reference %in% geneset)
    RninSet = length(reference) - RinSet
    GinSet = sum(genes %in% geneset)
    GninSet = length(genes) - GinSet
    fmat = matrix(c(GinSet, RinSet, GninSet, RninSet), nrow = 2,
                  ncol = 2, byrow = F)
    colnames(fmat) = c("inSet", "ninSet")
    rownames(fmat) = c("genes", "reference")
    
    # Perform Fisher's Exact Test
    fish = fisher.test(fmat, alternative = "greater")
    pval = fish$p.value
    inSet = RinSet + GinSet
    
    # Return the result as a tibble
    tibble(Name = nameSet,
           `Gene.number`= GinSet, 
           `Set.size` = inSet, 
           pval = pval,
           Genes = list(intersect(genes, geneset)))
  }) %>% bind_rows() %>%
    filter(Set.size>0) %>%
    mutate(padj = p.adjust(pval, method = "BH")) %>%
    arrange(pval)
  
  return(data.frame(rtab))
}


#' @name clusterEnrich
#' 
#' @title Perform Cluster Enrichment Analysis
#'
#' @description
#' `clusterEnrich` performs enrichment analysis on gene clusters, using Fisher's Exact Test to determine the significance of enrichment for each cluster.
#'
#' @param clusterTab A data frame containing cluster information, where each row corresponds to a gene and its assigned cluster.
#' @param se A SummarizedExperiment object containing gene expression data and metadata.
#' @param inputSet A list or data frame of gene sets to be used for enrichment analysis.
#' @param reference A character vector of reference genes. If NULL, it will be extracted from `se`. Default is `NULL`.
#' @param ptm Logical. If TRUE, the function will perform enrichment analysis on post-translational modification (PTM) gene sets. Default is `FALSE`.
#' @param adj Character. The method for adjusting p-values. Default is `"BH"`.
#' @param filterP Numeric. The p-value threshold for filtering significant results. Default is `0.05`.
#' @param ifFDR Logical. If TRUE, the function will use FDR-adjusted p-values for significance filtering. Default is `FALSE`.
#'
#' @return A list containing two elements:
#' \itemize{
#'   \item `table`: A data frame with enrichment results for each cluster and pathway.
#'   \item `plot`: A ggplot object showing the significance of enrichment for each pathway across clusters.
#' }
#'
#' @details
#' The function first retrieves or computes the reference set of genes or PTM sites. It then performs enrichment analysis for each cluster using the `runFisher` function.
#' The results are filtered based on the p-value threshold and adjusted for multiple testing if `ifFDR` is `TRUE`. The function generates a dot plot where the size and color of the points represent the significance of enrichment.
#'
#' @examples
#' # Example usage:
#' # Assuming `clusterTable` is a data frame with cluster assignments,
#' # `summarizedExp` is a SummarizedExperiment object,
#' # and `geneSetCollection` is a list or data frame of gene sets
#' # results <- clusterEnrich(clusterTable, summarizedExp, geneSetCollection, ptm = FALSE)
#'
#' @importFrom dplyr filter mutate group_by summarise ungroup bind_rows arrange
#' @importFrom ggplot2 ggplot geom_point aes scale_fill_gradient xlab ylab theme element_text
#' @importFrom tidyr separate
#' @importFrom stats p.adjust
#' @export
clusterEnrich <- function(clusterTab, se, inputSet, reference = NULL, ptm = FALSE, adj = "BH", filterP = 0.05, ifFDR = FALSE) {
  
  # If reference is not provided, derive it from the SummarizedExperiment object
  if (is.null(reference)) {
    if (ptm) {
      reference <- rowData(se)$site
    }
    else {
      reference <- unique(rowData(se)$Gene)
    }
  } 
  
  # Perform Fisher's Exact Test for each unique cluster
  resTabFisher <- lapply(unique(clusterTab$cluster), function(cc) {
    # Extract gene IDs for the current cluster
    id <- filter(clusterTab, cluster == cc)$feature
    if (ptm) {
      genes <- unique(elementMetadata(se)[id,]$site)
    }
    else {
      genes <- unique(elementMetadata(se)[id,]$Gene)
    }
    
    # Run Fisher's Exact Test and annotate results with the cluster ID
    eachOut <- runFisher(genes, reference, inputSet, ptm) %>%
      mutate(cluster = cc)
  }) %>% bind_rows()
  
  # Filter results based on significance and prepare for plotting
  if (ifFDR) {
    plotTab <- resTabFisher %>%
      filter(padj <= filterP) %>% arrange(padj) %>%
      mutate(ifSig = padj <= filterP) %>%
      group_by(Name) %>% mutate(atLeast1 = sum(ifSig)>0) %>%
      filter(atLeast1) %>% ungroup()
  }
  else {
    plotTab <- resTabFisher %>%
      filter(pval <= filterP) %>% arrange(pval) %>%
      mutate(ifSig = pval <= filterP) %>%
      group_by(Name) %>% mutate(atLeast1 = sum(ifSig)>0) %>%
      filter(atLeast1) %>% ungroup()
  }
  
  # Create a ggplot object for visualization of enrichment results
  p<- ggplot(plotTab, aes(x=cluster, y=Name, customdata = cluster, key = Name)) +
    geom_point(aes(size =-log10(pval),fill=-log10(pval)), shape = 21, color = "black") +
    scale_fill_gradient(low = "white", high = "red") +
    xlab("Cluster") +
    ylab("Pathway") +
    theme(axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          plot.title = element_text(hjust = 0.5, face = "bold"))
  return(list(table = plotTab, plot = p))
}

###### Helper functions for kinase activity inference ################

#' @name getDecouplerNetwork
#' 
#' @title Load Kinase-Substrate Interaction Network
#'
#' @description
#' `getDecouplerNetwork` loads the kinase-substrate interaction network for a specified species from pre-defined files.
#'
#' @param speciesRef A character string specifying the species. Supported values are "Homo sapiens" and "Mus musculus".
#'
#' @return A data frame containing the kinase-substrate interaction network for the specified species.
#'
#' @details
#' The function reads from tab-separated value (TSV) files located in the `omnipathR_kinase_network` directory.
#' It supports two species: Homo sapiens and Mus musculus.
#' 
#' @examples
#' # Load the human kinase-substrate interaction network
#' human_network <- getDecouplerNetwork("Homo sapiens")
#'
#' # Load the mouse kinase-substrate interaction network
#' mouse_network <- getDecouplerNetwork("Mus musculus")
#'
#' @importFrom utils read.table
#' @export
getDecouplerNetwork <- function(speciesRef) {
   
  # load network of kinase-substrate interaction from omnipathR_kinase_network folder
  if (speciesRef == "Homo sapiens") {
    decoupler_network <- read.table("omnipathR_kinase_network/Homo_sapiens.tsv", sep = "\t", stringsAsFactors = FALSE)
  } else if (speciesRef == "Mus musculus") {
    decoupler_network <- read.table("omnipathR_kinase_network/Mus_musculus.tsv", sep = "\t", stringsAsFactors = FALSE)
  }
  
  # Return the loaded network data
  return(decoupler_network)
}


#' @name calcKinaseScore
#' 
#' @title Calculate Kinase Activity Scores using `decoupleR`
#'
#' @description
#' `calcKinaseScore` calculates kinase activity scores based on input data and a specified network of regulatory relationships (decoupler network).
#'
#' @param resTab A data frame containing the input data with columns `site`, `stat`, and `log2FC`.
#' @param decoupler_network A data frame representing the decoupleR network with columns `source` and `target`.
#' @param corrThreshold A numeric value specifying the correlation threshold for filtering correlated regulons. Default is `0.9.
#' @param statType A character string specifying the type of statistic to use. Options are `"stat"` or `"log2FC"`. Default is `"stat"`.
#' @param nPerm Number of permutations for the null distribution. Default is `100`.
#'
#' @return A data frame with kinase activity scores, including columns for `source`, `score`, and `p_value`.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Removes duplicate rows based on the `site` column.
#'   \item Filters the data to include only those sites present in the `target` column of the `decoupler_network`.
#'   \item Prepares the input table based on the specified `statType`.
#'   \item Intersects the input table with the decoupler network to find common regulons.
#'   \item Checks for correlated regulons and filters out those exceeding the correlation threshold.
#'   \item Calculates kinase activity using a weighted mean approach.
#'   \item Processes the results to handle `NA` values and formats the output.
#' }
#'
#' @examples
#' # Example usage:
#' resTab <- data.frame(
#'   site = c("S1", "S2", "S3", "S4"),
#'   stat = c(1.5, -2.3, 0.7, 1.2),
#'   log2FC = c(0.5, -1.1, 0.3, 0.9)
#' )
#' decoupler_network <- data.frame(
#'   source = c("Kinase1", "Kinase2", "Kinase3"),
#'   target = c("S1", "S2", "S4")
#' )
#' result <- calcKinaseScore(resTab, decoupler_network, corrThreshold = 0.8, statType = "stat")
#' print(result)
#'
#' @importFrom dplyr distinct filter select rename mutate
#' @importFrom decoupleR intersect_regulons check_corr run_wmean
#' @importFrom tibble column_to_rownames
#' 
#' @export
calcKinaseScore <- function(resTab, decoupler_network, corrThreshold = 0.9, statType = "stat", nPerm = 100) {
  # Remove duplicate rows based on the 'site' column and keep all other columns
  resTab <- resTab %>%
    distinct(site, .keep_all = TRUE) %>%
    # Filter the rows where 'site' is present in the 'target' column of decoupler_network
    filter(site %in% decoupler_network$target)
  
  # Prepare the input table based on the specified statType
  if (statType == "stat") {
    inputTab <- resTab %>% select(site, stat) %>% dplyr::rename(t = stat) 
  } else if (statType == "log2FC") {
    inputTab <- resTab %>% select(site, log2FC) 
  }
  rownames(inputTab) <- NULL
  inputTab <- inputTab %>% data.frame() %>% column_to_rownames("site")
  # Intersect the input table with the decoupler network to find common regulons
  decoupler_network <- decoupleR::intersect_regulons(mat = inputTab, 
                                                     network = decoupler_network, 
                                                     .source = source, 
                                                     .target = target, 
                                                     minsize = 5)
  # Check for correlated regulons within the decoupler network and remove interactions with correlation >= threshold
  correlated_regulons <- decoupleR::check_corr(decoupler_network) %>%  #not necessary for now
    dplyr::filter(correlation >= corrThreshold)
  decoupler_network <- decoupler_network %>% 
    dplyr::filter(!source %in% correlated_regulons$source.2)
  # Calculate the kinase score by computing the weighted mean
  kinase_activity <- decoupleR::run_wmean(mat = as.matrix(inputTab), 
                                          network = decoupler_network,
                                          sparse = FALSE,
                                          times = nPerm)
  # Get the wmean statistics, replace NA scores with 0, and replace NA p_value with 1
  kinase_activity <- kinase_activity %>% dplyr::filter(statistic == "wmean") %>%
    select(-statistic, -condition) %>%
    mutate(score = ifelse(is.na(score), 0, score),
           p_value = ifelse(is.na(p_value), 1, p_value))
  return(kinase_activity)
}


# function to plot kinase score result for differential expression
# Parameters ----------------------
## scoreTab: table containing the result of kinase activity inference
## nTop: how many kinase to show for each direction 
## pCut: threshold of p_value to highlight significance in kinase activity
# output: a barplot of kinase score. Those with p_value lower than pCut will be highlighted red, others are grey


#' @name plotKinaseDE
#'
#' @title Plot Kinase score for Differential Expression data
#'
#' @description 
#' `plotKinaseDE` generates a bar plot of the top kinases associated with the differentially expressed genes based on their scores.
#'
#' @param scoreTab A data frame containing kinase scores with columns `source`, `score`, and `p_value`.
#' @param nTop An integer specifying the number of top kinases to plot for each direction. Default is `10`.
#' @param pCut A numeric value specifying the p-value cutoff for significance. Default is `0.05`.
#'
#' @return A ggplot object representing the bar plot of kinase score.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Adds a column for significance based on the p-value cutoff.
#'   \item Adds a column for the sign of the score.
#'   \item Filters out kinases with a score of 0.
#'   \item Selects the top `nTop` kinases by absolute score for each sign of the score.
#'   \item Creates a bar plot with the selected kinases.
#' }
#'
#' @examples
#' # Example usage:
#' scoreTab <- data.frame(
#'   source = c("Kinase1", "Kinase2", "Kinase3", "Kinase4"),
#'   score = c(2.3, -1.5, 0, 3.1),
#'   p_value = c(0.01, 0.2, 0.05, 0.03)
#' )
#' plot <- plotKinaseDE(scoreTab, nTop = 3, pCut = 0.05)
#' print(plot)
#'
#' @importFrom dplyr mutate filter group_by slice_max
#' @importFrom ggplot2 ggplot aes geom_bar scale_fill_manual theme_linedraw theme element_text unit coord_flip ggtitle xlab ylab reorder
#' @export
plotKinaseDE <- function(scoreTab, nTop = 10, pCut = 0.05) {
  plotTab <- scoreTab %>% mutate(significance = ifelse(p_value <= pCut, paste0("p <= ",pCut), paste0("p > ",pCut)),
                                 score_sign = sign(score)) %>%
    # Remove kinases whose scores are 0 in the plot
    filter(score_sign != 0) %>%
    # Group by score sign and select the top nTop kinases by absolute score
    group_by(score_sign) %>% slice_max(abs(score), n = nTop)
  p <- ggplot(plotTab, aes(x = reorder(source, score), y = score)) + 
    geom_bar(aes(fill = significance), stat = "identity") 
  
  # Customize the fill color based on the significance levels
  if (length(unique(plotTab$significance)) == 2) {
    p <- p + scale_fill_manual(values = c("indianred", "lightgrey"), labels = c(paste0("p <= ",pCut), paste0("p > ",pCut)))
  } else if (unique(plotTab$significance) == paste0("p > ",pCut)) {
    p <- p + scale_fill_manual(values = "lightgrey", labels = paste0("p > ",pCut)) 
  } else if (unique(plotTab$significance) == c(paste0("p <= ",pCut))) {
    p <- p + scale_fill_manual(values = "indianred", labels = paste0("p <= ",pCut)) 
  } 
  p <- p + 
    theme_linedraw() +
    theme(axis.title = element_text(face = "bold", size = 15),
          axis.text.y = element_text(size =15),
          axis.text.x = element_text(size = 15),
          plot.title = element_text(size = 17),
          legend.title = element_text(size =15, face= "bold"),
          legend.text = element_text(size =15),
          legend.key.height = unit(0.8, "cm"),
          legend.key.width = unit(0.8, "cm"),
          #panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    xlab("Kinases") + ylab("Kinase score") + 
    ggtitle(paste0("Kinase activity inference, top ", nTop, " kinases")) + 
    coord_flip()
  return(p)
}


#' Plot Kinase Activity Time Series
#'
#' This function creates a heatmap to visualize the result of kinase activity inference for time-series clustering, with significant activity changes marked.
#'
#' @param scoreTab A data frame containing kinase activity scores, p-values, and time points.
#' @param pCut A numeric value specifying the p-value threshold for significance. Default is `0.05`.
#' @param clusterName A character string specifying the name of the cluster for the plot title. Default is `"cluster1"`.
#'
#' @return A ggplot2 object representing the heatmap of kinase activity score.
#'
#' @details
#' The heatmap shows kinase activity scores over different time points. Significant activities (based on the specified p-value threshold) are marked with an asterisk (*). The color gradient represents the activity score, with blue indicating low activity, red indicating high activity, and white as the midpoint.
#' 
#' @examples
#' # Example usage:
#' scoreTab <- data.frame(
#'   timepoint = rep(c("0h", "1h", "2h"), each = 3),
#'   source = rep(c("KinaseA", "KinaseB", "KinaseC"), times = 3),
#'   score = runif(9, -2, 2),
#'   p_value = runif(9, 0, 0.1)
#' )
#' p <- plotKinaseTimeSeries(scoreTab)
#' print(p)
#'
#' @importFrom dplyr mutate rename
#' @importFrom ggplot2 ggplot aes geom_tile geom_text scale_fill_gradient2 scale_x_discrete scale_y_discrete theme_bw ylab xlab ggtitle theme element_text unit
#' @export
plotKinaseTimeSeries <- function(scoreTab, pCut = 0.05, clusterName = "cluster1") {
  
  # Add a significance marker based on the p-value threshold
  plotTab <- dplyr::mutate(scoreTab, sig = ifelse(p_value<=pCut, "*", ""))
  # Rename the score column for better readability in the plot
  plotTab <- plotTab %>% rename(Activity_score = "score")
  
  # Create the heatmap plot
  p <- ggplot(plotTab, aes(x=timepoint, y = source,fill = Activity_score)) +
    geom_tile() +
    geom_text(aes(label = sig), vjust = 0.5, hjust = 0.5, size = 10) +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
    scale_x_discrete(expand = c(0,0)) +scale_y_discrete(expand = c(0,0)) +
    theme_bw() +
    ylab("Kinase") + xlab("Time point") + ggtitle(paste("Kinase activity infererence,", clusterName)) + 
    theme(axis.title = element_text(face = "bold", size = 15),
          axis.text.y = element_text(size =15),
          axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust = 1),
          plot.title = element_text(size = 17),
          legend.title = element_text(size =15, face= "bold"),
          legend.text = element_text(size =15),
          legend.key.height = unit(0.8, "cm"),
          legend.key.width = unit(0.8, "cm"))
  return(p)
}


#### Function to perform PTM-SEA, a modification of the GSEA algorithm to work on databases of site-centric ptm signatures, 
# in which the direction of regulation (up or down) are specified for each site.
# Most of this function was similar to the project.geneset() function from the GitHub page for 
# ssGSEA2.0: https://github.com/broadinstitute/ssGSEA2.0/tree/master. 
#
# The function was modified to be compatible with SmartPhosExplorer and make recognizing phosphosites
# in the signatures stricter. Specifically, a phosphosite is included in a signature if its sign of the 
# test statistics (+ / -) matches the direction of regulation (u / d, abbreviate for `up` and `down`) 
# in the signature. This was not considered in the original ssGSEA2.0 algorithm.
# The function in this script only consider phosphosites (`site.ptm == "p"`). Signatures starting with "KINASE"
# are not considered since they are targets of kinases and hence would be similar to the kinase activity analysis.
# --------------------------------------
# For the publication associated with the algorithm and database please see Krug et al., 2019, at 
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6398202/ 
# For more descriptions of the signature sets please see the PTM signature database website:
# https://proteomics.broadapps.org/ptmsigdb/ 
# --------------------------------------
## Parameters
## geneStat: a dataframe with 1 column listing test statistics (t-stat or logFC) named "stat".
##           row.names are names of the phosphosites.
## ptmSetDb: database of PTM signature set.
## nPerm: Number of permutations
## weight: weighting of the sites based on their test statistic. If weight == 0 then the test statistics do not matter
## correl.type: Can be "rank", "z.score", or "symm.rank"
## statistic: Can be "Kolmogorov-Smirnov" or "area.under.RES"
## min.overlap: Minimum number of sites in the set to be considered for the analysis, default is 5.
# -----------------------------------------
## The resulting enrichment score is normalized to account for differences between signature set sizes


#' @name runGSEAforPhospho
#' 
#' @title Run GSEA for Phosphorylation Data
#'
#' @description 
#' `runGSEAforPhospho` performs Gene Set Enrichment Analysis (GSEA) for phosphorylation data.
#'
#' @param geneStat A data frame containing gene statistics, with gene names as row names and a column named 'stat' for the statistics.
#' @param ptmSetDb A data frame of post-translational modification (PTM) signature sets.
#' @param nPerm An integer specifying the number of permutations for the null distribution.
#' @param weight A numeric value for the weight parameter in the GSEA algorithm. If weight == 0 then the test statistics do not matter. Default is `1`.
#' @param correl.type A character string specifying the correlation type. Options are "rank", "symm.rank", and "z.score". Default is `"rank"`.
#' @param statistic A character string specifying the statistic to be used. Options are "Kolmogorov-Smirnov" and "area.under.RES". Default is `"Kolmogorov-Smirnov"`.
#' @param min.overlap An integer specifying the minimum overlap required between gene sets and the input data. Default is `5`.
#'
#' @return A tibble with enrichment scores and associated statistics for each PTM set.
#'
#' @details
#' This function runs GSEA on phosphorylation data to identify enriched PTM sets. It calculates enrichment scores and p-values for each set, normalizes the scores, and adjusts p-values for multiple testing.
#' 
#' @examples
#' # Example usage:
#' geneStat <- data.frame(stat = runif(100, -2, 2))
#' row.names(geneStat) <- paste0("Gene", 1:100)
#' ptmSetDb <- data.frame(signature = sample(letters, 100, replace = TRUE), category = "example", site.ptm = "p", site.direction = sample(c("u", "d"), 100, replace = TRUE))
#' result <- runGSEAforPhospho(geneStat, ptmSetDb, nPerm = 1000)
#' print(result)
#'
#' @importFrom dplyr mutate rename filter count tibble as_tibble group_by ungroup arrange bind_rows
#' @importFrom tidyr separate
#' @importFrom stats p.adjust
#' @export
runGSEAforPhospho <- function(geneStat, ptmSetDb, nPerm, weight = 1, correl.type = "rank",
                              statistic = "Kolmogorov-Smirnov", min.overlap = 5) {
  
  # Internal function to calculate GSEA enrichment score
  gseaScorePTM <- function (ordered.gene.list, data.expr, gene.set2, 
                            weight = 1, correl.type = "rank", gene.set.direction = NULL,
                            statistic = "Kolmogorov-Smirnov", min.overlap = 5) {
    
    # Function to calculate the enrichment score (ES)
    score <- function(max.ES, min.ES, RES, gaps, valleys, statistic){
      # KM
      if( statistic == "Kolmogorov-Smirnov" ){
        if( max.ES > -min.ES ){
          ES <- signif(max.ES, digits=5)
          arg.ES <- which.max(RES)
        } else{
          ES <- signif(min.ES, digits=5)
          arg.ES <- which.min(RES)
        }
      }
      # AUC
      if( statistic == "area.under.RES"){
        if( max.ES > -min.ES ){
          arg.ES <- which.max(RES)
        } else{
          arg.ES <- which.min(RES)
        }
        gaps = gaps+1
        RES = c(valleys,0) * (gaps) + 0.5*( c(0,RES) - c(valleys,0) ) * (gaps)
        ES = sum(RES)
      }
      return(list(RES=RES, ES=ES, arg.ES=arg.ES))
    } 
    
    n.rows = length(ordered.gene.list)
    
    # Apply weighting to the correlation vector
    if (weight == 0) {
      
      correl.vector <- rep(1, n.rows)
      
    } else if (weight > 0) {
      # If weighting is used (weight > 0), bring 'correl.vector' into the same order as the ordered gene list
      if (correl.type == "rank") {
        correl.vector <- data.expr[ordered.gene.list]
        
      } else if (correl.type == "symm.rank") {
        correl.vector <- data.expr[ordered.gene.list]
        
        correl.vector <- ifelse(correl.vector > correl.vector[ceiling(n.rows/2)],
                                correl.vector,
                                correl.vector + correl.vector - correl.vector[ceiling(n.rows/2)])
      } else if (correl.type == "z.score") {
        x <- data.expr[ordered.gene.list]
        correl.vector <- (x - mean(x))/sd(x)
      }
    }
    
    # Length of gene list is same as the number of rows in input matrix
    N = length(ordered.gene.list)
    
    
    # Sirectionality of the gene set
    if(!is.null(gene.set.direction)){
      
      # Number of 'd' features
      d.idx <- which(gene.set.direction=='d')
      Nh.d <- length(d.idx)
      Nm.d <-  N - Nh.d
      
      # Locations of 'd' features
      tag.d <- sign( match(ordered.gene.list, gene.set2[ d.idx ], nomatch=0) )
      if (weight == 0) {
        ind.d = which(tag.d == 1)} else {
          ind.d = which(tag.d == 1 & correl.vector < 0)}
      number.d = length(ind.d)
      
      # Number of 'u' features
      u.idx <- which(gene.set.direction=='u')
      Nh.u <- length(u.idx)
      Nm.u <-  N - Nh.u
      
      # Locations of 'up' features
      tag.u <- sign( match(ordered.gene.list, gene.set2[ u.idx ], nomatch=0) )
      if (weight == 0) {
        ind.u = which(tag.u == 1)} else {
          ind.u = which(tag.u == 1 & correl.vector >= 0)}
      number.u = length(ind.u)
      
      
      # For up-regulated genes/sites
      if(number.u > 1){
        
        # Extract and apply weighting
        correl.vector.u <- correl.vector[ind.u]
        correl.vector.u <- abs(correl.vector.u)^weight           ## weighting
        
        sum.correl.u <- sum(correl.vector.u)
        
        up.u <- correl.vector.u/sum.correl.u         ## steps up in th random walk
        gaps.u <- (c(ind.u-1, N) - c(0, ind.u))      ## gaps between hits
        down.u <- gaps.u/Nm.u                        ## steps down in the random walk
        
        RES.u <- cumsum(up.u-down.u[1:length(up.u)])  
        
        valleys.u = RES.u-up.u
        
        max.ES.u = suppressWarnings(max(RES.u))
        min.ES.u = suppressWarnings(min(valleys.u))
        
        # Calculate final score
        score.res <- score(max.ES.u, min.ES.u, RES.u, gaps.u, valleys.u, statistic)
        ES.u <- score.res$ES
        arg.ES.u <- score.res$arg.ES
        RES.u <- score.res$RES
        
      } else {
        correl.vector.u <- rep(0, N)
        ES.u=0
        RES.u=0
        arg.ES.u=NA
        up.u=0
        down.u=0
      }
      
      # For down-regulated genes/sites
      if(number.d > 1){  
        # Extract and apply weighting
        correl.vector.d <- correl.vector[ind.d]
        correl.vector.d <- abs(correl.vector.d)^weight           ## weighting
        
        sum.correl.d <- sum(correl.vector.d)
        
        up.d <- correl.vector.d/sum.correl.d
        gaps.d <- (c(ind.d-1, N) - c(0, ind.d))
        down.d <- gaps.d/Nm.d
        
        RES.d <- cumsum(up.d-down.d[1:length(up.d)])               
        valleys.d = RES.d-up.d
        
        max.ES.d = suppressWarnings(max(RES.d))
        min.ES.d = suppressWarnings(min(valleys.d))
        
        # Calculate final score
        score.res <- score(max.ES.d, min.ES.d, RES.d, gaps.d, valleys.d, statistic)
        ES.d <- score.res$ES
        arg.ES.d <- score.res$arg.ES
        RES.d <- score.res$RES
        
      } else {
        correl.vector.d <- rep(0, N)
        ES.d=0
        RES.d=0
        ind.d=NA
        number.d=0
        arg.ES.d=NA
        up.d=0
        down.d=0
      }
  
      # Make sure to meet the min.overlap threshold
      if(Nh.d == 1 & Nh.u < min.overlap | Nh.u == 1 & Nh.d < min.overlap){
        ES.u <- ES.d <- RES.u <- RES.d <- 0
        arg.ES <- arg.ES <- NA
        ind.u <- ind.d <- NULL
      }
      
      # Combine the results
      ES <- ES.u - ES.d
      RES <- list(u=RES.u, d=RES.d)
      arg.ES <- c(arg.ES.u, arg.ES.d)
      correl.vector = list(u=correl.vector.u, d=correl.vector.d)

      ind <- list(u=ind.u, d=ind.d)
      step.up <- list(u=up.u, d=up.d )
      step.down <- list(u=1/Nm.u, d=1/Nm.d)
      gsea.results = list(ES = ES, ES.all = list(u=ES.u, d=ES.d), arg.ES = arg.ES, RES = RES, indicator = ind, correl.vector = correl.vector, step.up=step.up, step.down=step.down,
                          number.u = number.u, number.d = number.d)
      
      
    } else { ## end  if(!is.null(gene.set.direction))
      
      # Without directionality
      Nh <- length(gene.set2)
      Nm <-  N - Nh
      
      
      # Match gene set to data
      tag.indicator <- sign(match(ordered.gene.list, gene.set2, nomatch=0))    # notice that the sign is 0 (no tag) or 1 (tag)
      # Positions of gene set in ordered gene list
      ind = which(tag.indicator==1)
      # 'correl.vector' is now the size of 'gene.set2'
      correl.vector <- abs(correl.vector[ind])^weight
      # Sum of weights
      sum.correl = sum(correl.vector)
      
      # Determine peaks and valleys
      # Divide correl vector by sum of weights
      up = correl.vector/sum.correl     # "up" represents the peaks in the mountain plot
      gaps = (c(ind-1, N) - c(0, ind))  # gaps between ranked pathway genes
      down = gaps/Nm
      
      RES = cumsum(c(up,up[Nh])-down)
      valleys = RES[1:Nh]-up
      
      max.ES = max(RES)
      min.ES = min(valleys)
      
      # Calculate final score
      score.res <- score(max.ES, min.ES, RES[1:Nh], gaps, valleys, statistic)
      
      ES <- score.res$ES
      arg.ES <- score.res$arg.ES
      RES <- score.res$RES
      
      gsea.results = list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = ind, correl.vector = correl.vector, step.up=up, step.down=1/Nm)
    } 
    
    return (gsea.results)
  }
  
  
  # Remove KINASE signature since this is analogous to the kinase activity inference part
  ptmSetDbNoKinase <- ptmSetDb %>%
    filter(!grepl("KINASE", category))
  
  # Get the number of PTM sites for each signature
  ptmSiteCount <- ptmSetDbNoKinase %>%
    count(signature) %>%
    rename(no.PTM.site = "n")
  
  # Preprocessing the geneSetDatabase
  phosphoSetDb <- ptmSetDbNoKinase %>%
    dplyr::as_tibble() %>%
    filter(site.ptm == "p") %>%   
    group_by(signature) %>%
    filter(n() >= 5) %>%
    ungroup() %>%
    separate(site.annotation, sep =  ":", into = c("site", "PubMedID"), extra="merge", fill="right")
  
  # Get the number of phospho sites for each signature
  phosphoSiteCount <- phosphoSetDb %>%
    count(signature) %>%
    rename(no.phospho.site = "n")
  
  # Put input data in a format compatible with gseaScorePTM
  ordered.gene.list <- row.names(geneStat)
  data.expr <- geneStat$stat
  names(data.expr) <- ordered.gene.list
  # Run GSEA for each PTM set
  rtab <- lapply(phosphoSiteCount$signature, function(signature) {
    # Get number of PTM site and phospho site in the database
    nPTMsite = as.numeric(ptmSiteCount[ptmSiteCount$signature == signature, "no.PTM.site"])
    nPpSite = as.numeric(phosphoSiteCount[phosphoSiteCount$signature == signature, "no.phospho.site"])
    signatureSet = phosphoSetDb[phosphoSetDb$signature == signature,]
    gene.set2 = signatureSet$site
    gene.set.direction = signatureSet$site.direction
    gene.set.PMID = signatureSet$PubMedID
    # Calculate the gsea score
    if (sum(row.names(geneStat) %in% gene.set2) < min.overlap) {
      enrichScoreNorm <- enrichScore <- pvalue <- number.u <- number.d <- 0
    } else {
      resGSEA = gseaScorePTM(ordered.gene.list, data.expr =  data.expr, gene.set2 = gene.set2, 
                             weight = weight, correl.type = correl.type,
                             gene.set.direction = gene.set.direction, min.overlap =  min.overlap)
      enrichScore = resGSEA$ES
      if (!is.null(gene.set.direction)) {
        number.u  = resGSEA$number.u
        number.d = resGSEA$number.d
      } else {
        number.u <- number.d <- 0 
      }
      # Calculate the null distribution and pvalue
      if (nPerm == 0) {
        enrichScoreNorm = enrichScore
        pvalue = 1
      } else {
        nullDistES = sapply(1:nPerm,  function(x) gseaScorePTM(sample(ordered.gene.list), data.expr=data.expr, gene.set2=gene.set2, 
                                                               weight, correl.type, gene.set.direction = gene.set.direction, min.overlap = min.overlap)$ES)
        nullDistES = unlist(nullDistES)
        if (enrichScore >= 0) {
          nullDistES.pos = nullDistES[nullDistES >= 0]
          if (length(nullDistES.pos) == 0) nullDistES.pos = 0.5
          posMean = mean(nullDistES.pos)
          enrichScoreNorm = enrichScore/posMean
          s = sum(nullDistES.pos >= enrichScore)/length(nullDistES.pos)
          pvalue = ifelse(s == 0, 1/nPerm, s)
        } else {
          nullDistES.neg = nullDistES[nullDistES < 0]
          if (length(nullDistES.neg) == 0) nullDistES.neg = 0.5
          negMean = mean(nullDistES.neg)
          enrichScoreNorm = enrichScore/negMean
          s = sum(nullDistES.neg <= enrichScore)/length(nullDistES.neg)
          pvalue = ifelse(s == 0, 1/nPerm, s)
        }
      }
    }
    tibble(Name = signature,
           nSite = number.u + number.d,                      # number of phosphosites in the input data
           enrichScore = enrichScoreNorm,                    # normalized enrichment score to correct for differences in signature sizes
           n.P.site.in.Db = nPpSite,                         # number of phosphosites in the database
           n.PTM.site.in.Db = nPTMsite,                      # number of PTM sites in the database
           pvalue = pvalue,
           number.u = number.u,
           number.d = number.d)
  }
  ) %>% bind_rows() %>%
    filter(nSite>= min.overlap, enrichScore!=0) %>%
    mutate(p.adj = p.adjust(pvalue, method = "BH")) %>%
    arrange(pvalue) 
  return(rtab)
}

