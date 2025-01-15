
######################################################
#### PreCSenM (Predictive Cellular Senescence Model)
#### Cellular Senescence score calculation function  
#### Author: Lifei Ma
######################################################

#' @title Cellular Senescence score calculation function 
#' @author Lifei Ma \email{lifei_ma@@126.com}
#'
#' @description
#' \code{PreCSenM} Prediction of Cellular Senescence score from gene expression profiling data.
#' \code{sig.genes} Users must load CS significant genes data
#' 
#' @details
#' multiple samples (Not One) to be predicted as CS score based on their expression profiles of the gene panel.
#' 
#' @param ExprMat Samples to be predicted. (Note: Gene expression values of the samples need to be normalized (scaling),rownames=genes | colnames=samples) )
#' 
#' @return A dataframe with 4 columns:
#' \item{Sample_ID}{Samples to be predicted.}
#' \item{score.01}{Prediction results ("1 or higher" indicates CS-high, and "0 or lower" CS-low).}
#' @examples
#' @export


# load 236 CS significant genes 
# Users must load CS significant genes data
load(file = "sig.genes.Rdata")

# Establish function

PreCSenM <- function(sig.genes, ExprMat) {
  
  ## Check arguments
  if (missing(ExprMat) || !is.numeric(as.matrix(ExprMat)) || !class(ExprMat) %in% c("matrix", "data.frame") || dim(ExprMat)[2] < 2) {
    stop("'ExprMat' is missing, non-numeric, or incorrect format")
  }
  
  if (missing(sig.genes) || !is.data.frame(sig.genes) || !"gene" %in% colnames(sig.genes) || !"Estimate" %in% colnames(sig.genes)) {
    stop("'sig.genes' is missing or incorrectly formatted")
  }
  
  # Expression Matrix 
  data = as.data.frame(ExprMat)
  
  # Check genes
  genes = intersect(rownames(data), sig.genes$gene)
  message("The core CS genes number are: ", length(genes))
  
  if (length(genes) > nrow(sig.genes) * 0.75) {
    sig.genes = sig.genes[sig.genes$gene %in% genes, ]
    data.tmp = data[sig.genes$gene, ]
    data.tmp = data.tmp[order(match(rownames(data.tmp), sig.genes$gene)), ]
    sig.genes = sig.genes[order(match(sig.genes$gene, rownames(data.tmp))), ]
    
    # Calculate score
    score = data.tmp * sig.genes$Estimate
    CS.score = data.frame(score = rowSums(t(score)))
    CS.score$scaled_score = scale(CS.score$score)
    CS.score$score.01 = (CS.score$score - min(CS.score$score)) / (max(CS.score$score) - min(CS.score$score))
    CS.score$sample = rownames(CS.score)
    
    message("----------------------\nCalculation is done!")
    return(CS.score)
  } else {
    stop("Error: Too few genes! Calculation stopped.", call. = FALSE)
  }
}

# Calculate the CS score
#CS.score = PreCSenM(sig.genes = sig.genes, ExprMat = ExprMat)








