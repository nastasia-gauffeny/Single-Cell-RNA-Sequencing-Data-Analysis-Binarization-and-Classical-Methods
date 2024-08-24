#' Identify expressed genes
#'
#' Identify expressed genes in Seurat object.
#'
#' @param object Seurat Object
#' @param min.pct minimum expressing fraction. Default: 0.1. Ignored if min.cell is specified.
#' @param min.cell minimum number of expressing cells. If specified, min.pct is ignored.
#' @param group Character specifying metadata field to group cells by. If not specified, global expression fraction is evaluated. If specified, group-level gene lists are combined used group.boolean.
#' @param group.boolean Boolean used to combine group genelists. One of "OR" or "AND". Default: "OR". Argument is ignored if 'group' is not specified.
#' @name getExpressedGenes
#' @author Nicholas Mikolajewicz
#' @return vector of gene names
#' @examples
#'
#' split.var <- "seurat_clusters"
#' which.genes <- getExpressedGenes(object = so.query, min.pct = 0.1, group = split.var, group.boolean = "OR")
#'
getExpressedGenes_m <- function(object, min.pct = 0.1, min.cell = NULL, group = NA, group.boolean = "OR"){
  
  # min.groups: specify min.pct fraction grouping
  # NA: min.pct satisfied global
  # character vector speifying metadata feature to group: min.pct satisfied within object meta.data grouping
  
  # min.type.operator
  # AND: must satisfy criteria in ALL subgroup
  # OR: must satisfy critreria in ATLEAST ONE group
  
  emat <- getExpressionMatrix_m(object, which.data = "data")
  
  if (!is.null(min.cell)){
    cfunc <- Matrix::rowSums
    min.pct <- min.cell
  } else {
    cfunc <- Matrix::rowMeans
  }
  
  if (is.na(group)){
    
    # pct.rep <- apply(emat, 1, function(x) mean(x>0))
    pct.rep <- cfunc(emat>0) # faster implementation
    
    expressed.genes <- names(pct.rep)[pct.rep > min.pct]
  } else if (group %in% colnames(object@meta.data)){
    group.var <- as.character(object@meta.data[ ,group])
    u.gv <- unique(group.var)
    expressed.genes.all <- c()
    for (i in 1:length(u.gv)){
      pct.rep <- NULL
      if (!is.null(dim(emat[ ,group.var %in% u.gv[i]]))){
        pct.rep <- cfunc(emat[ ,group.var %in% u.gv[i]]>0)
        expressed.genes.all <- c(expressed.genes.all, names(pct.rep)[pct.rep > min.pct])
      }
    }
    
    df.exp.gene <- data.frame(table(expressed.genes.all))
    
    if (group.boolean == "OR"){
      expressed.genes <- unique(as.character(df.exp.gene$expressed.genes.all))
    } else if (group.boolean == "AND"){
      expressed.genes <- unique(as.character(df.exp.gene$expressed.genes.all[df.exp.gene$Freq == max(df.exp.gene$Freq, na.rm = T)]))
    } else {
      miko_message(paste0("'", group.boolean, "' is not an accepted argument for group.boolean. Must be either 'OR' or 'AND'. 'OR' was used as default argument"))
      expressed.genes <- unique(as.character(df.exp.gene$expressed.genes.all))
    }
    
  } else {
    miko_message(paste0( "'", group,  "' was not found in seurat object. Failed to identify expressed genes. "))
  }
  
  return(expressed.genes)
}

