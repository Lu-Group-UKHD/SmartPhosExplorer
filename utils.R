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