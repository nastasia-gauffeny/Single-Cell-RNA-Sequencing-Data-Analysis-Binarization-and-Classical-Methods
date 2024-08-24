#' Get expression matrix from Seurat Object
#'
#' Get expression matrix from Seurat Object
#'
#' @param so Seurat Object
#' @param only.variable Logical indicating whether to include variable features only or not.
#' @param which.assay Seurat assay to get data from. Default is DefaultAssay(so).
#' @param which.data Specify which data to use (refers to slots in Seurat object assay). One of:
#' \itemize{
#' \item "scale" - Default
#' \item "data"
#' }
#' @param use.additional.genes Character vector of additional genes to include (in addition to varibale, if variable flag is specificed). Default is NA.
#' @param as.dense Logical to convert sparse to dense matrix. Only applies if which.data is 'data'. Default is FALSE.
#' @name getExpressionMatrix
#' @author Nicholas Mikolajewicz
#' @return gene x cell expression matrix
#'
getExpressionMatrix_m <- function(so, only.variable = F, which.assay = NULL, which.data = "scale", use.additional.genes = NA, as.dense = F){
  
  # specify assay
  if (is.null(which.assay)) which.assay <- DefaultAssay(so)
  
  # get complete matrix
  if (which.data == "scale"){
    exp.mat.complete <- so@assays[[which.assay]]@scale.data
  } else if (which.data == "data"){
      exp.mat.complete <- (so@assays$RNA$data)
  }
  
  
  if (only.variable){
    var.feat <-  VariableFeatures(so)
    if (length(var.feat) == 0){
      if (which.assay == "SCT"){
        if ("integrated" %in% names(so@assays)){
          var.feat <-  so@assays[["integrated"]]@var.features
        } else {
          var.feat <-  so@assays[["SCT"]]@var.features
        }
      } else {
        stop("Could not find variable features in Seurat Object. ")
      }
    }
    
    exp.mat <- exp.mat.complete[rownames(exp.mat.complete) %in% var.feat, ]
  } else {
    exp.mat <- exp.mat.complete
  }
  
  
  if (!is.na(use.additional.genes)){
    which.missing <- which(!(use.additional.genes %in% rownames(exp.mat)))
    missing.additional.genes <- (use.additional.genes)[which.missing]
    exp.mat.additional <- exp.mat.complete[rownames(exp.mat.complete) %in% missing.additional.genes, ]
    exp.mat <- rbind(exp.mat, exp.mat.additional)
  }
  
  return(exp.mat)
}


