#shiny::runApp("/Users/shubhamagrawal/Documents/work/apps/processSMART/", launch.browser = TRUE, port = 8000)
# library(proDA)
# library(SummarizedExperiment)
# n_samples <- 40
# n_feat <- 1000
# data <- generate_synthetic_data(n_feat, n_conditions = n_samples / 10, 
#                                 n_replicates = n_samples / 4, frac_changed = 0.1)
# a <- data$Y
# colnames(a) <- gsub(colnames(a), pattern = "Condition", replacement = "Sample")
# 
# ## add some treatment-specific effects
# set.seed(1)
# a[, 1:5] <- a[, 1:5] + rnorm(5000, mean = 1.0, sd = 0.5)
# a[, 11:15] <- a[, 11:15] + rnorm(5000, mean = 0.8, sd = 0.5)
# a[, 21:25] <- a[, 21:25] + rnorm(5000, mean = 1.2, sd = 0.5)
# a[, 31:35] <- a[, 31:35] + rnorm(5000, mean = 0.7, sd = 0.5)
# 
# ## create information on the samples
# type_sample <- gsub(data$groups, pattern = "Condition", replacement = "Type")
# trmt_sample <- paste(
#   c(rep("1", 10), rep("2", 10), rep("3", 10), rep("4", 10)),
#   c(rep("A", 5), rep("B", 5)), sep = "_")
# cD <- data.frame(name = colnames(a), type = type_sample, 
#                  treatment = trmt_sample)
# 
# ## create information on the proteins
# rD <- data.frame(spectra = rownames(a))
# 
# ## create se
# se <- SummarizedExperiment(assay = a, rowData = rD, colData = cD)
# 
# ## remove the features that have only NA values
# se <- se[!apply(assay(se), 1, function(row_i) all(is.na(row_i))), ]
MatrixQCvis::shinyQC(se)
