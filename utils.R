# Helper functions

# Get the last symbol of a protein that has multiple gene symbols
getOneSymbol <- function(Gene) {
  outStr <- sapply(Gene, function(x) {
    sp <- str_split(x, ";")[[1]]
    sp[length(sp)]
  })
  names(outStr) <- NULL
  outStr
}

########################################## Normalization Correction ######################################

#helper function for median or mean normalization
medianNorm <- function(x, method = "median") {
  if (method == "median") {
    mVal <- matrixStats::colMedians(x, na.rm=TRUE)
    mVal <- mVal - median(mVal, na.rm=TRUE)
  } else if (method == "mean") {
    mVal <- colMeans(x, na.rm=TRUE)
    mVal <- mVal - mean(mVal, na.rm=TRUE)
  }
  mMat <- matrix(rep(mVal, each = nrow(x)), ncol =ncol(x))
  return(x-mMat)
}

#function for performing normalization of FP and PP samples
performCombinedNormalization <- function(maeData) {
  
  # get count matrix from FP samples
  setFP <- maeData[,maeData$sampleType == "FullProteome"]
  protFP <- assay(setFP[["Proteome"]])
  phosFP <- assay(setFP[["Phosphoproteome"]])
  comFP <- rbind(protFP, phosFP)
  comFP <- comFP[rowSums(!is.na(comFP))>0,]
  
  #perform median normalization and log2 transformation
  comFP.norm <- medianNorm(log2(comFP))
  
  return(comFP.norm)
}

removePrefix <- function(name, prefix) {
  name <- gsub(prefix, "", name)
  name <- gsub("^[._ ]","", name) #remove trailing characters after remove "FullProteome","Phospho" suffix or prefix
  name <- gsub("[._ ]$","", name)
  return(name)
}

getRatioMatrix <- function(maeData, normalization = FALSE, getAdjustedPP = FALSE) {
  # Get the ratios of phospho-proteins (or peptides) intensity in PP samples devided by FP samples 
  # maeData: multiAssayExpriment output from SmartPhos package
  # normalization: whether normalization needs to be performed on fullProteome samples. Normalization is necessary normalization has not been performed by Spectronaut.
  stopifnot(is.logical(normalization))
  
  if (!getAdjustedPP) {
    phosPP <- log2(assay(maeData[,maeData$sampleType == "Phospho"][["Phosphoproteome"]]))
  } else {
    phosPP <- log2(assays(maeData[,maeData$sampleType == "Phospho"][["Phosphoproteome"]])[["Intensity_adjusted"]])
  }
  
  if (!normalization) {
    
    phosFP <- log2(assay(maeData[,maeData$sampleType == "FullProteome"][["Phosphoproteome"]]))
    
  } else  {
    
    phosFP <- performCombinedNormalization(maeData)
    
  }
  
  #remove sample type prefix or suffix
  colnames(phosFP) <- removePrefix(colnames(phosFP),"FullProteome")
  colnames(phosPP) <- removePrefix(colnames(phosPP),"Phospho")
  
  #get ratio matrix
  allSmp <- intersect(colnames(phosFP), colnames(phosPP))
  allRow <- intersect(rownames(phosFP), rownames(phosPP))
  ratioMat <- phosPP[allRow, allSmp] - phosFP[allRow, allSmp]
  ratioMat <- ratioMat[rowSums(!is.na(ratioMat)) >0,]
  
  return(ratioMat)
  
}

plotLogRatio <- function(maeData, normalization = FALSE) {
  # Plot the log ration of PP/FP intensities
    
  ratioMat <- getRatioMatrix(maeData, normalization)
  
  phosPP <- log2(assay(maeData[,maeData$sampleType == "Phospho"][["Phosphoproteome"]]))
  medianPP <- colMedians(phosPP,na.rm = TRUE)
  names(medianPP) <- removePrefix(colnames(phosPP), "Phospho")
  
  plotTab <- as_tibble(ratioMat, rownames = "feature") %>%
    pivot_longer(-feature) %>%
    filter(!is.na(value)) %>%
    mutate(medianPP = medianPP[name])
  
  ggplot(plotTab, aes(x=name, y=value)) +
    geom_boxplot(aes(fill = medianPP)) +
    ggtitle("Boxplot of Phospho/FullProteome Ratio") +
    xlab("sample") +
    ylab("log2(ratio)") +
    geom_hline(yintercept = median(median(plotTab$value, na.rm=TRUE)), linetype = "dashed", color = "red") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          plot.title = element_text(hjust = 0.5, face = "bold"))
  
}

checkRatioMat <- function(ratioMat, minOverlap = 3) {
  #check the PP/FP ratio matrix and remove feature that do not meet requirements
    
  excludeSampleList <- c()
  
  #are there any samples that don't have any phospho sites detect in both FP and PP samples?
  noOverSmp <- colnames(ratioMat)[colSums(!is.na(ratioMat))==0]
  if (length(noOverSmp) >0) {
    warning(paste0("Below samples don't have phopho-peptides detected in both enriched (PP) and unenriched (FP) samples and therefore adjusting factor will set to 0 (no adjustment) for them:\n",
                   paste0(noOverSmp, collapse = ", ")))
    excludeSampleList <- c(excludeSampleList, noOverSmp)
  }
  
  ratioMat <- ratioMat[, !colnames(ratioMat) %in% noOverSmp]
  
  #are there any samples don't have enough peptide overlap with other samples?
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

runPhosphoAdjustment <- function(maeData, normalization = FALSE, minOverlap = 3, completeness = 0, ncore = 1 ) {
  
  #function to opitmize
  esFun <- function(par, data) {
    comPair <- utils::combn(seq(length(par)), 2)
    sum(((data[comPair[1, ],] + par[comPair[1, ]]) - (data[comPair[2, ],] + par[comPair[2, ]]))^2/rowSums(!is.na(data[comPair[1,],] + data[comPair[2,],])), na.rm = TRUE)
  }
  
  #get PP/FP ratio matrix
  ratioMat <- getRatioMatrix(maeData, normalization = normalization)
  adjFac <- structure(rep(0, length.out = ncol(ratioMat)), names = colnames(ratioMat))
  
  #subset features according to completeness in the ratio matrix
  ratioMat <- ratioMat[rowSums(!is.na(ratioMat))/ncol(ratioMat) >= completeness,]
  
  #sanity check to see if any samples need to be excluded
  excList <- checkRatioMat(ratioMat)
  ratioMat <- ratioMat[, !colnames(ratioMat) %in% excList]
  
  #set an initial value for B based on col medians of ratioMat, may increase search speed
  colMed <- apply(ratioMat,2, median, na.rm = TRUE)
  iniPar <- median(colMed) - colMed
  
  #estimating adjusting factor
  cl <- makeCluster(ncore)
  setDefaultCluster(cl = cl)
  optRes <- optimParallel(par=iniPar, fn=esFun, data=t(ratioMat))
  stopCluster(cl)
  
  #add adjusting factor to sample annotation  
  adjFac[names(optRes$par)] <- optRes$par
  ppName <- colnames(maeData[,maeData$sampleType == "Phospho"][["Phosphoproteome"]])
  adjFac <- structure(adjFac[removePrefix(ppName,"Phospho")], names = ppName)
  maeData$adjustFactorPP <- unname(adjFac[match(rownames(colData(maeData)),names(adjFac))])
  
  #adjust phospho measurement on PP samples
  phosMat <- assay(maeData[,maeData$sampleType == "Phospho"][["Phosphoproteome"]])
  phosMat <- t(t(phosMat)*(2^adjFac))
  assays(maeData[["Phosphoproteome"]])[["Intensity_adjusted"]] <- assays(maeData[["Phosphoproteome"]])[["Intensity"]]
  assays(maeData[["Phosphoproteome"]])[["Intensity_adjusted"]][,colnames(phosMat)] <- phosMat
  
  return(maeData)
}

plotAdjustmentResults <- function(maeData, normalization = FALSE) {
  
  if (!"adjustFactorPP" %in% colnames(colData(maeData))) {
    stop("Phosphorylation measurments have not been adjusted yet. Please perform normalization adjustment using calcAdjustFacotr function first")
  }
  
  #visualize precursors 
  ratioMat.ori <- getRatioMatrix(maeData, normalization = normalization, getAdjustedPP = FALSE)
  ratioMat.adj <- getRatioMatrix(maeData, normalization = normalization, getAdjustedPP = TRUE)
  ratioPlotTab <- bind_rows(pivot_longer(as_tibble(ratioMat.ori, rownames = "id"), -id, names_to = "sample", values_to = "ratio") %>% mutate(adjustment = "before adjustment"),
                            pivot_longer(as_tibble(ratioMat.adj, rownames = "id"), -id, names_to = "sample", values_to = "ratio") %>% mutate(adjustment = "after adjustment")) %>%
    mutate(adjustment = factor(adjustment, levels = c("before adjustment","after adjustment"))) %>%
    filter(!is.na(ratio))
  
  #for precursors present in all samples
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
  
  #for ratio box plots
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
  
  
  #for phosphomeasurement of PP samples before and after adjustment
  
  ppMat.adj <- assays(maeData[,maeData$sampleType == "Phospho"][["Phosphoproteome"]])[["Intensity_adjusted"]] 
  ppMat.ori <- assays(maeData[,maeData$sampleType == "Phospho"][["Phosphoproteome"]])[["Intensity"]] 
  ppPlotTab <- bind_rows(pivot_longer(as_tibble(ppMat.ori, rownames = "id"), -id, names_to = "sample", values_to = "count") %>% mutate(adjustment = "before adjustment"),
                         pivot_longer(as_tibble(ppMat.adj, rownames = "id"), -id, names_to = "sample", values_to = "count") %>% mutate(adjustment = "after adjustment")) %>%
    mutate(adjustment = factor(adjustment, levels = c("before adjustment","after adjustment")),
           count = log2(count)) %>%
    filter(!is.na(count))
  
  #for ratio box plots
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



###########################################################################################################

# Function to preprocess proteomic data
preprocessProteome <- function(seData, filterList = NULL, missCut = 50,
                               transform = "log2", normalize = FALSE, getPP = FALSE,
                               removeOutlier = NULL, impute = "QRILC",
                               verbose = FALSE, scaleFactorTab = NULL) {
  
  # seData, is the multiAassayExperiment object from SmartPhos
  # filterList, is a list object where users can filter the samples based on sample annotations.
  # missCut, is the maximal missing value percentage allowed
  # getPP, is whether to get PP sample instead of default FP samples
  # removeOutlier, can be a character vector contains the outlier samples to be removed from preprocessing.
  # verbose,  whether to show additional information
  # scaleFactor, for user-specified scale factor. Not useful for the shiny app for now. Remove it or keep it as NULL.
  
  if (getPP) {
    # normally only full proteome samples (FP) are used for protein analysis.
    # But if the user wishes, PP samples can also be retrieved. 
    fpe <- seData[,seData$sampleType == "Phospho"]
    colData(fpe) <- colData(seData)[colnames(fpe),]
  } else {
    fpe <- seData[,seData$sampleType == "FullProteome"]
    colData(fpe) <- colData(seData)[colnames(fpe),]
  }
  
  if (!is.null(removeOutlier)) {
    fpe <- fpe[, !fpe$sample %in% removeOutlier]
  }
  
  # filter for selected criteria
  if (!is.null(filterList)) {
    for (n in names(filterList)) {
      fpe <- fpe[,fpe[[n]] %in% filterList[[n]]]
    }
  }
  # rename column names
  colnames(fpe) <- fpe$sample
  # get last gene name
  rowData(fpe)$Gene <- getOneSymbol(rowData(fpe)$Gene)
  # remove features without symbols
  fpe <- fpe[!rowData(fpe)$Gene %in% c(NA,""),]
  # filter for missing values
  countMat <- assay(fpe)
  missPer <- rowSums(is.na(countMat))/ncol(countMat)*100
  fpeSub <- fpe[missPer < missCut,]
  # nomralization and transformation
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
  
  # imputation 
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
    else {
      doParallel::registerDoParallel(cores = 6)  # set based on number of CPU cores
      doRNG::registerDoRNG(seed = 123)
      mf <- missForest::missForest(t(assay(fpeSub)), parallelize = "forests", maxiter = 2, ntree = 50)
      imp <- t(mf$ximp)
    }
    assays(fpeSub)[["imputed"]] <- assay(imp)
    rowData(fpeSub)$name <- NULL
    rowData(fpeSub)$ID <- NULL
    if (verbose) {
      # show number of samples and features
      print("Number of proteins and samples:")
      print(dim(fpeSub))
      }
  }
  
  return(fpeSub)
}

# Function to preprocess phospho proteomic data
preprocessPhos <- function(seData, filterList = NULL, missCut = 50,
                           transform="log2", normalize = FALSE, getFP = FALSE,
                           removeOutlier = NULL, assayName = NULL,
                           scaleFactorTab = NULL, impute = "QRILC", verbose = FALSE) {
  
  # This function is largely the same as the function above, but it's intended
  # for processing the phosphoproteomic data there is an additional parameter,
  # assayName, which can be used to specify an assay name beyond "Phosphoproteme"
  # to retrieve other assay types, for example, the ratio between phospho and
  # proteome measurement. Currently not useful for the shiny app. 
  
  
  if (is.null(assayName)) {
    if (getFP) { 
      # normally PP samples are used for phosphoproteomic, but if the user wishes, FP sample can be used.
      ppe <- seData[,seData$sampleType == "FullProteome"]
      colData(ppe) <- colData(seData)[colnames(ppe),]
    } else {
      ppe <- seData[,seData$sampleType == "Phospho"]
      colData(ppe) <- colData(seData)[colnames(ppe),]
    }
  } else {
    ppe <- seData[[assayName]]
    colData(ppe) <- colData(seData[,colnames(ppe)])
  }
  
  if (!is.null(removeOutlier)) {
    ppe <- ppe[, !ppe$sample %in% removeOutlier]
  }
  
  # filter for selected criteria
  if (!is.null(filterList)) {
    for (n in names(filterList)) {
      ppe <- ppe[,ppe[[n]] %in% filterList[[n]]]
    }
  }
  
  # rename column names
  colnames(ppe) <- ppe$sample
  # get last gene name
  rowData(ppe)$Gene <- getOneSymbol(rowData(ppe)$Gene)
  # get last phosphorylation site
  rowData(ppe)$Residue <- getOneSymbol(rowData(ppe)$Residue)
  rowData(ppe)$Position <- getOneSymbol(rowData(ppe)$Position)
  # remove features without symbols
  ppe <- ppe[!rowData(ppe)$Gene %in% c(NA,""),]
  # rename site
  rowData(ppe)$site <- paste0(rowData(ppe)$Gene,"_",rowData(ppe)$Residue,rowData(ppe)$Position)
  # filter for missing values
  countMat <- assay(ppe)
  missPer <- rowSums(is.na(countMat))/ncol(countMat)*100
  ppeSub <- ppe[missPer < missCut,]
  
  # nomralization and transformation
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
  
  # imputation
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
    else {
      doParallel::registerDoParallel(cores = 6)  # set based on number of CPU cores
      doRNG::registerDoRNG(seed = 123)
      mf <- missForest::missForest(t(assay(ppeSub)), parallelize = "forests", maxiter = 2, ntree = 50)
      imp <- t(mf$ximp)
    }
    assays(ppeSub)[["imputed"]] <- assay(imp)
    rowData(ppeSub)$name <- NULL
    rowData(ppeSub)$ID <- NULL
    # show number of samples and features
    if (verbose) {
      print("Number of proteins and samples:")
      print(dim(ppeSub))
    }
  }
  return(ppeSub)
}



# Function to plot time-series clustering results
clusterTS <- function(x, k, pCut = NULL, twoCondition = FALSE) {
  
  # x, input matrix, rows as feature, columns as time points
  # k, number of clusters
  
  # remove rows with NA
  x.center <- x[complete.cases(x),]
  res <- e1071::cmeans(x.center, k)

  resCluster <- tibble(feature = names(res$cluster),
                       cluster = res$cluster,
                       prob = rowMax(res$membership))
  
  if (!is.null(pCut)) resCluster <- filter(resCluster, prob >= pCut)

  if (!twoCondition) {
    # get the time unit and the time levels (i.e. order of the timepoints in the plot)
    timeVector <- unique(colnames(x.center))
    timeUnit <- str_extract(timeVector, "h|min")
    timeUnit <- ifelse(is.na(timeUnit), "", timeUnit)
    # If both h and min are present, divide the min time points by 60
    if ((any(timeUnit == "h")) & (any(timeUnit == "min"))) {
      timeValue <- timeVector
      timeValue[timeUnit == "min"] <- 1/60 * as.numeric(gsub("min", "", timeValue[timeUnit == "min"]))
      timeRank <- rank(as.numeric(gsub("h", "", timeValue)))
    } else {
      timeRank <- rank(as.numeric(gsub("h|min", "", timeVector)))
    }
    # the order that the time points will appear in the plot
    timeOrder <- timeVector[order(match(timeRank, sort(timeRank)))]   
    # plot clustering results
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
    # plot clustering for 2 conditions
    # get the time unit and the time levels (i.e. order of the timepoints in the plot)
    timeVector <- sapply(colnames(x.center), function(X) unlist(str_split(X, "_"))[1])
    timeVector <- unique(timeVector)
    timeUnit <- str_extract(timeVector, "h|min")
    timeUnit <- ifelse(is.na(timeUnit), "", timeUnit)
    # If both h and min are present, divide the min time points by 60
    if ((any(timeUnit == "h")) & (any(timeUnit == "min"))) {
      timeValue <- timeVector
      timeValue[timeUnit == "min"] <- 1/60 * as.numeric(gsub("min", "", timeValue[timeUnit == "min"]))
      timeRank <- rank(as.numeric(gsub("h", "", timeValue)))
    } else {
      timeRank <- rank(as.numeric(gsub("h|min", "", timeVector)))
    }
    # the order that the time points will appear in the plot
    timeOrder <- timeVector[order(match(timeRank, sort(timeRank)))]
    # plot clustering result
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
    
    p <- ggplot(clusterTab, aes( x=time, y= value, group = geneGroup)) +
      geom_line( aes(alpha = prob, color= treatment)) +
      # scale_color_gradient(low = "green", high = "red") +
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
splineFilter <- function(exprMat, subjectID = NULL, time, df, pCut, ifFDR, treatment = NULL, refTreatment = NULL) {
  # The time points must have either no unit or in h and/or minute. 
  # If both h and min are present, the minute time points will be converted to h
  if ((all(str_ends(time,"h|min"))) & (!all(str_ends(time,"h"))) & (!all(str_ends(time,"min")))) {
    time[str_ends(time, "min")] <- 1/60 * as.numeric(gsub("min","", time[str_ends(time, "min")]))
  }
  time <- as.numeric(gsub("h|min", "", time))
  
  # remove NA rows
  exprMat <- exprMat[complete.cases(exprMat), ]
  if (is.null(treatment)) {
    # one condition or logFC
    # include subjectID in designTab if provided
    if (is.null(subjectID)) {
      designTab <- data.frame(row.names = colnames(exprMat))
      designTab$X <- splines::ns(time, df)
      design <- model.matrix(~ 0 + X, data = designTab)
    } else { 
      designTab <- data.frame(row.names = colnames(exprMat), subjectID = subjectID)
      designTab$X <- splines::ns(time, df)
      design <- model.matrix(~ 0 + X + subjectID, data = designTab)
    }
    
    fit <- lmFit(exprMat, design = design)
    fit2 <- eBayes(fit)
    resTab <- topTable(fit2, coef = seq(df), number = Inf) %>%
      as_tibble(rownames = "ID") 
    
    if (ifFDR) resTab$p <- resTab$adj.P.Val else resTab$p <- resTab$P.Value
    
    resTab <- filter(resTab, p <= pCut)
    
    return(exprMat[resTab$ID,])
  } else {
    
    # two condition clustering
    if (is.null(subjectID)) {
      designTab <- data.frame(row.names = colnames(exprMat), treatment = treatment)
      designTab$treatment <- factor(designTab$treatment, levels = unique(treatment))
      designTab$treatment <- relevel(designTab$treatment, refTreatment)
      designTab$treatment <- droplevels(designTab$treatment)
      designTab$X <- splines::ns(time, df)
      design <- model.matrix(~ 0 + X*treatment, data = designTab)
    } else {
      designTab <- data.frame(row.names = colnames(exprMat), subjectID = subjectID, treatment = treatment)
      designTab$treatment <- factor(designTab$treatment, levels = unique(treatment))
      designTab$treatment <- relevel(designTab$treatment, refTreatment)
      designTab$treatment <- droplevels(designTab$treatment)
      designTab$X <- splines::ns(time, df)
      design <- model.matrix(~ 0 + subjectID + X*treatment, data = designTab)
    }
    fit <- lmFit(exprMat, design = design)
    fit2 <- eBayes(fit)
    resTab <- topTable(fit2, coef = (ncol(design)-df+1):ncol(design), number = Inf) %>%
      as_tibble(rownames = "ID")
    
    if (ifFDR) resTab$p <- resTab$adj.P.Val else resTab$p <- resTab$P.Value
    
    resTab <- filter(resTab, p <= pCut)
    
    return(exprMat[resTab$ID,])
    
  }
}

# scaling function for clustering
mscale <- function(x, center = TRUE, scale = TRUE, censor = NULL, useMad = FALSE){
  if (scale & center) {
    if (useMad) {
      x.scaled <- apply(x, 1, function(y) (y-median(y,na.rm = T))/meanAD(y))
    } else {
      x.scaled <- apply(x, 1, function(y) (y-mean(y,na.rm=T))/sd(y,na.rm = T))
    }
  } else if (center & !scale) {
    if (useMad) {
      x.scaled <- apply(x, 1, function(y) (y-median(y,na.rm=T)))
    } else {
      x.scaled <- apply(x, 1, function(y) (y-mean(y,na.rm=T)))
    }
  } else if (!center & scale) {
    if (useMad) {
      x.scaled <- apply(x, 1, function(y) y/meanAD(y))
    } else {
      x.scaled <- apply(x, 1, function(y) y/sd(y,na.rm = T))
    }
  } else {
    x.scaled <- t(x)
  }
  
  if (!is.null(censor)) {
    if (length(censor) == 1) {
      x.scaled[x.scaled > censor] <- censor
      x.scaled[x.scaled < -censor] <- -censor
    } else {
      x.scaled[x.scaled > censor[2]] <- censor[2]  # higher limit
      x.scaled[x.scaled < censor[1]] <- censor[1]  # lower limit
    }
  }
  return(t(as.matrix(x.scaled)))
}

# function to run fisher test for enrichment analysis (time series clustering only)
## Note: gmtFile is the directory to the .gmt file for pathway enrichment analysis but is also the directory
##       to the .txt file containing the PTM database, if ptm == TRUE
## Processing the database file will depend on whether the file is geneset or ptm set.
## If ptm == TRUE, the ptmset will be used, and each geneset will be split into 2 based on the sites direction
## of regulation (up or down)
runFisher <- function (genes, reference, gmtFile, ptm = FALSE) {
  # retrieve the database
  if (!ptm) {
    genesets <- piano::loadGSC(gmtFile)$gsc
    setList <- 1:length(genesets)
  } else {
    genesets <- read.table(gmtFile, sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>%
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
  reference = reference[!reference %in% genes]
  
  rtab = lapply(setList, function(i) { # here i is the order number of a set in geneset or the name of the set in  ptm set.
    if (!ptm) {
      geneset = genesets[[i]]
      nameSet = names(genesets)[i]
    } else {
      geneset = genesets[genesets$signature == i, "site"]
      nameSet = i
    }
    RinSet = sum(reference %in% geneset)
    RninSet = length(reference) - RinSet
    GinSet = sum(genes %in% geneset)
    GninSet = length(genes) - GinSet
    fmat = matrix(c(GinSet, RinSet, GninSet, RninSet), nrow = 2,
                  ncol = 2, byrow = F)
    colnames(fmat) = c("inSet", "ninSet")
    rownames(fmat) = c("genes", "reference")
    fish = fisher.test(fmat, alternative = "greater")
    pval = fish$p.value
    inSet = RinSet + GinSet
    tibble(Name = nameSet,
           `Gene.number`= GinSet, 
           `Set.size` = inSet, 
           pval = pval)
  }) %>% bind_rows() %>%
    filter(Set.size>0) %>%
    mutate(padj = p.adjust(pval, method = "BH")) %>%
    arrange(pval)
  
  return(data.frame(rtab))
}

###### Helper functions for kinase activity inference ################

# function to build the prior knowledge phosphorylation network using OmnipathR
# Parameter:
## speciesRef: reference species, currently accepts either "Homo sapiens" or "Mus musculus"
# note: dephosphorylation events are removed, but this could be left as an option?
# Output:
## a network of kinase- phosphorylation site interactions
getDecouplerNetwork <- function(speciesRef) {
   
  # Load network of kinase-substrate interaction from omnipathR_kinase_network folder
  if (speciesRef == "Homo sapiens") {
    decoupler_network <- read.table("omnipathR_kinase_network/Homo_sapiens.tsv", sep = "\t", stringsAsFactors = FALSE)
  } else if (speciesRef == "Mus musculus") {
    decoupler_network <- read.table("omnipathR_kinase_network/Mus_musculus.tsv", sep = "\t", stringsAsFactors = FALSE)
  }
  #print(paste("Loaded", nrow(decoupler_network), "kinase-substrate interactions from OmnipathR"))
}

# function to calculate kinase score using decoupleR
# Parameters --------------------------
## resTab: must have 2 columns, one specifying the statistic to be used (t-statistic or logFC)
## decoupler_network: Network of kinases and their targets derived from Omnipath
## corrThreshold: threshold to check for colinearity in the decoupler_network
## statType: is either "stat" or "log2FC"
## nPerm: number of permutations for the null distribution
# Output ---------------------------
# A dataframe with 4 columns:
## source: name of the kinase
## score: inferred activity score of the kinase
## p_value: p-value from the decoupler::run_wmean function
calcKinaseScore <- function(resTab, decoupler_network, corrThreshold = 0.9, statType = "stat", nPerm = 100) {
  # get differential phosphorylation sites
  resTab <- resTab %>%
    distinct(site, .keep_all = TRUE) %>%
    filter(site %in% decoupler_network$target)
  
  if (statType == "stat") {
    inputTab <- resTab %>% select(site, stat) %>% dplyr::rename(t = stat) 
  } else if (statType == "log2FC") {
    inputTab <- resTab %>% select(site, log2FC) 
  }
  rownames(inputTab) <- NULL
  inputTab <- inputTab %>% data.frame() %>% column_to_rownames("site")
  decoupler_network <- decoupleR::intersect_regulons(mat = inputTab, 
                                                     network = decoupler_network, 
                                                     .source = source, 
                                                     .target = target, 
                                                     minsize = 5)
  # check for colinearity and remove interactions with correlation >= threshold
  correlated_regulons <- decoupleR::check_corr(decoupler_network) %>%  #not necessary for now
    dplyr::filter(correlation >= corrThreshold)
  decoupler_network <- decoupler_network %>% 
    dplyr::filter(!source %in% correlated_regulons$source.2)
  # calculate the kinase score by computing the weighted mean
  kinase_activity <- decoupleR::run_wmean(mat = as.matrix(inputTab), 
                                          network = decoupler_network,
                                          sparse = FALSE,
                                          times = nPerm)
  # get the wmean statistics, replace NA scores with 0, and replace NA p_value with 1
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
plotKinaseDE <- function(scoreTab, nTop = 10, pCut = 0.05) {
  plotTab <- scoreTab %>% mutate(significance = ifelse(p_value <= pCut, paste0("p <= ",pCut), paste0("p > ",pCut)),
                                 score_sign = sign(score)) %>%
    filter(score_sign != 0) %>%  # remove kinases whose scores are 0 in the plot
    group_by(score_sign) %>% slice_max(abs(score), n = nTop)
  p <- ggplot(plotTab, aes(x = reorder(source, score), y = score)) + 
    geom_bar(aes(fill = significance), stat = "identity") 
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

# function to plot the result of kinase activity inference for time-series clustering
# Parameter -------------------------------
## scoreTab: dataframe containing the result of kinase activity inference result for one cluster
## pCut: threshold for highlighting significance
## clusterName: name of the cluster (in the GUI, this is inherited from input$seleCluster)
# Output -----------------------------------
## a heatmap of kinase activity score. Those with p_value < pCut will be highlighted with an asterisk
plotKinaseTimeSeries <- function(scoreTab, pCut = 0.05, clusterName = "cluster1") {
  plotTab <- dplyr::mutate(scoreTab, sig = ifelse(p_value<=pCut, "*", ""))
  plotTab <- plotTab %>% rename(Activity_score = "score")
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
runGSEAforPhospho <- function(geneStat, ptmSetDb, nPerm, weight = 1, correl.type = "rank",
                              statistic = "Kolmogorov-Smirnov", min.overlap = 5) {
  # function from the ssGSEA publication, adapted to work with our format
  ## #############################################################################
  ##
  ##           function to calculate GSEA enrichment score
  ## - apply correlation scheme and weighting
  ## - calculate ES
  ## ############################################################################
  gseaScorePTM <- function (ordered.gene.list, data.expr, gene.set2, 
                            weight = 1, correl.type = "rank", gene.set.direction = NULL,
                            statistic = "Kolmogorov-Smirnov", min.overlap = 5) {
    
    ##################################################################
    ## function to calculate ES score
    score <- function(max.ES, min.ES, RES, gaps, valleys, statistic){
      ## KM
      if( statistic == "Kolmogorov-Smirnov" ){
        if( max.ES > -min.ES ){
          ES <- signif(max.ES, digits=5)
          arg.ES <- which.max(RES)
        } else{
          ES <- signif(min.ES, digits=5)
          arg.ES <- which.min(RES)
        }
      }
      ## AUC
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
    } ## end function score
    
    n.rows = length(ordered.gene.list)
    
    ## #######################################
    ## weighting
    ## #######################################
    if (weight == 0) {
      
      correl.vector <- rep(1, n.rows)
      
    } else if (weight > 0) {
      ## if weighting is used (weight > 0), bring
      ## 'correl.vector' into the same order
      ## as the ordered gene list
      if (correl.type == "rank") {
        ##correl.vector <- data.array[ordered.gene.list, sample.index]
        correl.vector <- data.expr[ordered.gene.list]
        
      } else if (correl.type == "symm.rank") {
        ##correl.vector <- data.array[ordered.gene.list, sample.index]
        correl.vector <- data.expr[ordered.gene.list]
        
        correl.vector <- ifelse(correl.vector > correl.vector[ceiling(n.rows/2)],
                                correl.vector,
                                correl.vector + correl.vector - correl.vector[ceiling(n.rows/2)])
      } else if (correl.type == "z.score") {
        ##x <- data.array[ordered.gene.list, sample.index]
        x <- data.expr[ordered.gene.list]
        correl.vector <- (x - mean(x))/sd(x)
      }
    }
    
    ## length of gene list - equals number of rows in input matrix
    N = length(ordered.gene.list)
    
    ## #####################################
    ## directionality of the gene set
    if(!is.null(gene.set.direction)){
      
      ## number of 'd' features
      d.idx <- which(gene.set.direction=='d')
      Nh.d <- length(d.idx)
      Nm.d <-  N - Nh.d
      ## locations of 'd' features
      tag.d <- sign( match(ordered.gene.list, gene.set2[ d.idx ], nomatch=0) )
      if (weight == 0) {
        ind.d = which(tag.d == 1)} else {
          ind.d = which(tag.d == 1 & correl.vector < 0)}
      number.d = length(ind.d)
      
      ## number of 'u' features
      u.idx <- which(gene.set.direction=='u')
      Nh.u <- length(u.idx)
      Nm.u <-  N - Nh.u
      ## locations of 'up' features
      tag.u <- sign( match(ordered.gene.list, gene.set2[ u.idx ], nomatch=0) )
      if (weight == 0) {
        ind.u = which(tag.u == 1)} else {
          ind.u = which(tag.u == 1 & correl.vector >= 0)}
      number.u = length(ind.u)
      
      
      ########################################
      ## up-regulated part
      ########################################
      if(number.u > 1){
        
        ## extract and apply weighting
        correl.vector.u <- correl.vector[ind.u]
        correl.vector.u <- abs(correl.vector.u)^weight           ## weighting
        
        sum.correl.u <- sum(correl.vector.u)
        
        up.u <- correl.vector.u/sum.correl.u         ## steps up in th random walk
        gaps.u <- (c(ind.u-1, N) - c(0, ind.u))      ## gaps between hits
        down.u <- gaps.u/Nm.u                        ## steps down in the random walk
        
        RES.u <- cumsum(up.u-down.u[1:length(up.u)])  # OLD: RES.u <- cumsum(c(up.u,up.u[Nh.u])-down.u)
        
        valleys.u = RES.u-up.u
        
        max.ES.u = suppressWarnings(max(RES.u))
        min.ES.u = suppressWarnings(min(valleys.u))
        
        ## calculate final score
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
      
      ## ######################################
      ## down-regulated part
      ## ######################################
      if(number.d > 1){  
        ## extract and apply weighting
        correl.vector.d <- correl.vector[ind.d]
        correl.vector.d <- abs(correl.vector.d)^weight           ## weighting
        
        sum.correl.d <- sum(correl.vector.d)
        
        up.d <- correl.vector.d/sum.correl.d
        gaps.d <- (c(ind.d-1, N) - c(0, ind.d))
        down.d <- gaps.d/Nm.d
        
        RES.d <- cumsum(up.d-down.d[1:length(up.d)])               ## RES.d <- cumsum(c(up.d,up.d[Nh.d])-down.d)
        valleys.d = RES.d-up.d
        
        max.ES.d = suppressWarnings(max(RES.d))
        min.ES.d = suppressWarnings(min(valleys.d))
        
        ## calculate final score
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
      ## ############################
      ## make sure to meet the min.overlap
      ## threshold
      if(Nh.d == 1 & Nh.u < min.overlap | Nh.u == 1 & Nh.d < min.overlap){
        ES.u <- ES.d <- RES.u <- RES.d <- 0
        arg.ES <- arg.ES <- NA
        ind.u <- ind.d <- NULL
      }
      
      ## ###########################
      ## combine the results
      ES <- ES.u - ES.d
      RES <- list(u=RES.u, d=RES.d)
      arg.ES <- c(arg.ES.u, arg.ES.d)
      ##tag.indicator <- rep(0, N)
      correl.vector = list(u=correl.vector.u, d=correl.vector.d)
      ##if(!is.null(ind.u))tag.indicator[ind.u] <- 1
      ##if(!is.null(ind.d))tag.indicator[ind.d] <- -1
      ind <- list(u=ind.u, d=ind.d)
      step.up <- list(u=up.u, d=up.d )
      ##                   step.down <- list(u=down.u, d=down.d )
      step.down <- list(u=1/Nm.u, d=1/Nm.d)
      gsea.results = list(ES = ES, ES.all = list(u=ES.u, d=ES.d), arg.ES = arg.ES, RES = RES, indicator = ind, correl.vector = correl.vector, step.up=step.up, step.down=step.down,
                          number.u = number.u, number.d = number.d)
      
      ## ##############################################################
      ##
      ##      original ssGSEA code without directionality
      ##
      ## ##############################################################
      
    } else { ## end  if(!is.null(gene.set.direction))
      
      Nh <- length(gene.set2)
      Nm <-  N - Nh
      
      ## #####################################
      ## match gene set to data
      tag.indicator <- sign(match(ordered.gene.list, gene.set2, nomatch=0))    # notice that the sign is 0 (no tag) or 1 (tag)
      ## positions of gene set in ordered gene list
      ind = which(tag.indicator==1)
      ## 'correl.vector' is now the size of 'gene.set2'
      correl.vector <- abs(correl.vector[ind])^weight
      ## sum of weights
      sum.correl = sum(correl.vector)
      
      #########################################
      ## determine peaks and valleys
      ## divide correl vector by sum of weights
      up = correl.vector/sum.correl     # "up" represents the peaks in the mountain plot
      gaps = (c(ind-1, N) - c(0, ind))  # gaps between ranked pathway genes
      down = gaps/Nm
      
      RES = cumsum(c(up,up[Nh])-down)
      valleys = RES[1:Nh]-up
      
      max.ES = max(RES)
      min.ES = min(valleys)
      
      ## calculate final score
      score.res <- score(max.ES, min.ES, RES[1:Nh], gaps, valleys, statistic)
      
      ES <- score.res$ES
      arg.ES <- score.res$arg.ES
      RES <- score.res$RES
      
      gsea.results = list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = ind, correl.vector = correl.vector, step.up=up, step.down=1/Nm)
    } ## end else
    
    return (gsea.results)
  }
  
  
  # remove KINASE signature since this is analogous to the kinase activity inference part
  ptmSetDbNoKinase <- ptmSetDb %>%
    filter(!grepl("KINASE", category))
  
  # get the number of PTM sites for each signature
  ptmSiteCount <- ptmSetDbNoKinase %>%
    count(signature) %>%
    rename(no.PTM.site = "n")
  
  # preprocessing the geneSetDatabase
  phosphoSetDb <- ptmSetDbNoKinase %>%
    dplyr::as_tibble() %>%
    filter(site.ptm == "p") %>%   
    group_by(signature) %>%
    filter(n() >= 5) %>%
    ungroup() %>%
    separate(site.annotation, sep =  ":", into = c("site", "PubMedID"), extra="merge", fill="right")
  
  # get the number of phospho sites for each signature
  phosphoSiteCount <- phosphoSetDb %>%
    count(signature) %>%
    rename(no.phospho.site = "n")
  
  # put input data in a format compatible with gseaScorePTM
  ordered.gene.list <- row.names(geneStat)
  data.expr <- geneStat$stat
  names(data.expr) <- ordered.gene.list
  # run GSEA for each PTM set
  rtab <- lapply(phosphoSiteCount$signature, function(signature) {
    # get number of PTM site and phospho site in the database
    nPTMsite = as.numeric(ptmSiteCount[ptmSiteCount$signature == signature, "no.PTM.site"])
    nPpSite = as.numeric(phosphoSiteCount[phosphoSiteCount$signature == signature, "no.phospho.site"])
    signatureSet = phosphoSetDb[phosphoSetDb$signature == signature,]
    gene.set2 = signatureSet$site
    gene.set.direction = signatureSet$site.direction
    gene.set.PMID = signatureSet$PubMedID
    # calculate the gsea score
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
      # calculate the null distribution and pvalue
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

