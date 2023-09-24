
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
PreCSenM <- function(sig.genes, ExprMat){
  
  ## Check arguments
  if (missing(ExprMat) || !class(ExprMat) %in% c("matrix", "data.frame") || dim(ExprMat)[2] < 3)
    stop("'ExprMat' is missing or incorrect")
  
  # Expression Matrix 
  data = as.data.frame(ExprMat)
  
  # check genes 
  genes = intersect(rownames(data), sig.genes$gene)
  #length(genes)  # <= 236
  message("The core CS genes number are: ",length(genes))
  
  # when CS significant genes > 236*0.75 = 177 
  if (length(genes) > nrow(sig.genes)*0.75) {
    
    # select genes
    sig.genes = sig.genes[sig.genes$gene %in% genes, ]
    
    # select expression data  
    data.tmp = data[genes,]
    table(rownames(data.tmp) == sig.genes$gene)
    
    score = data.tmp * sig.genes$Estimate
    #score = na.omit(score)
    CS.score = t(score)
    
    CS.score=data.frame(score=rowSums(CS.score))
    #CS.score$mean = as.numeric(CS.score$score/nrow(score))
    CS.score$score.scale = scale(CS.score$score)
    CS.score$score.01 = (CS.score$score-min(CS.score$score))/(max(CS.score$score)-min(CS.score$score))
    #CS.score$normalize.oppo = 1-CS.score$normalize 
    CS.score$sample = rownames(CS.score)
    
    message("----------------------\nCalculation is done!")
    
    return(CS.score)
    
  }
  else {
    warning("Error: Too few genes!\n----------------------\nError: Calculation Stopped!!",call. = FALSE)
  }
  
}


# Calculate the CS score
CS.score = PreCSenM(sig.genes = sig.genes, ExprMat = ExprMat)

