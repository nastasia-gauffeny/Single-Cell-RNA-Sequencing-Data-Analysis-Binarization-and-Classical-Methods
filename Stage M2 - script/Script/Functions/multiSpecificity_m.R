#' Evaluate specificity of single-cell markers across several cluster resolutions.
#'
#' Evaluate specificity of single-cell markers across several cluster resolutions using co-dependency index-based specificity measure. Consider running multiCluster(...) first.
#'
#' @param object Seurat object with multi-resolution clusters provided in meta data.
#' @param cluster_names Vector specifying names of all cluster configurations found in meta data.
#' @param features If specified, marker specificity analysis is limited to specified features. Otherwise all features are used (more computationally intensive).
#' @param deg_prefilter If TRUE, wilcoxon analysis is performed first to subset DEG features for downstream analysis. Results in faster performance. Default is TRUE.
#' @param geosketch.subset Use GeoSketch method to subsample scRNA-seq data while preserving rare cell states (https://doi.org/10.1016/j.cels.2019.05.003). Logical, T or F (Default F). Recommended if cell type representation is imbalanced.
#' @param cdi_bins Vector specifying binning for CDI-based specificity curve. Must range [0,1]. Default is seq(0, 1, by = 0.01).
#' @param min.pct Minimal expression of features that are considered in specificity analysis. Represents fraction of expression cells and must range [0,1]. Higher values result in faster performance. Default is 0.1.
#' @param n.workers Number of workers used for parallel implementation. Default is 1.
#' @param return_dotplot If TRUE, dot plots visualizing expression of top specific markers are returned. Default is T.
#' @param verbose Print progress. Default is TRUE.
#' @name multiSpecificity
#' @seealso \code{\link{multiCluster}}
#' @author Nicholas Mikolajewicz
#' @return Seurat object
#' @examples
#' # clustering data
#' ms.list <- multiSpecificity(object = so.query, cluster_names = cluster.name, features = NULL, deg_prefilter = T,
#' cdi_bins = seq(0, 1, by = 0.01), min.pct = 0.1,
#' n.workers = 4, return_dotplot = T,  verbose = T)
#' df.auc.spec <- ms.list$specificity_summary
#' qm.res.sum.sum.all <- ms.list$specificity_raw
#' plt.clust.spec <- ms.list$auc_plot
#' plt.auc.spec <- ms.list$resolution_plot
#' plt.auc.dot <- ms.list$dot_plot
multiSpecificity_m <- function(object, cluster_names, features = NULL, deg_prefilter = T, geosketch.subset = F,
                             cdi_bins = seq(0, 1, by = 0.01), min.pct = 0.1, n.workers = 1, return_dotplot = T, verbose = T){
  
  require(parallel);
  require(foreach);
  
  if (verbose){
    mylapply <- pbapply::pblapply
  } else {
    mylapply <- lapply
  }
  if(!("Seurat" %in% class(object))) stop("'object' is not a Seurat object")
  available_clusters <- unique(cluster_names[cluster_names %in% colnames(object@meta.data)])
  if (length(available_clusters) == 0) stop(cluster_names, " not found in 'object' meta data")
  miko_message("Assessing specificity scores for ", length(available_clusters), " unique groupings...", verbose = verbose)
  
  df.group_size <- as.data.frame(apply(object@meta.data[ ,available_clusters], 2, ulength))
  colnames(df.group_size) <- "n"; df.group_size$group = rownames(df.group_size)
  maxClustName <- df.group_size$group[which.max(df.group_size$n)]
  
  # get expressed genes
  expr_genes <- getExpressedGenes_m(object = object, min.pct = min.pct, group = maxClustName)
  if (is.null(features)) features <- rownames(object)
  
  # browser()
  
  # deg prefilter (to improve computational speed)
  if (deg_prefilter){
    
    miko_message("Prefiltering features by presto differential expression analysis...", verbose = verbose)
    
    try({
      
      
      
      if (is.null(n.workers)) n.workers <- 1
      if (n.workers > parallel::detectCores()) n.workers <- parallel::detectCores()
      all.deg.list <- multiDEG(object = object, groups = cluster_names,
                               only_pos = T,
                               nworkers =n.workers,
                               fdr_threshold = 1,
                               logfc_threshold = 0, verbose = T )
      
      deg_top <- unique(unlist(mylapply(all.deg.list, function(x){
        x <- x %>% dplyr::group_by(group) %>% dplyr::top_n(100, auc)
        unique(x$feature)
      })))
      
      features <- features[features %in% deg_top]
      
      
    }, silent = T)
    
  }
  
  features <- features[features %in% expr_genes]
  
  
  # use presto to nominate top AUC and run CDI on subset. faster performance?
  cdi_cluster <- findCDIMarkers(object = object,
                                geosketch.subset = geosketch.subset,
                                features.x = available_clusters,
                                n.workers = n.workers,
                                features.y = features)
  
  cdi_cluster_top <- cdi_cluster %>%
    dplyr::group_by(feature.x) %>%
    dplyr::mutate(cdi_rank = rank(ncdi, ties.method = "random")) %>%
    dplyr::top_n(1, cdi_rank)
  
  df.feature.x <- NULL
  for (i in 1:length(available_clusters)){
    meta.list <- group2list(object = object, group = available_clusters[i])
    group_names <- names(meta.list)
    group_names_full <- paste0(available_clusters[i], "_", names(meta.list))
    df.feature.x <- bind_rows(
      df.feature.x,
      data.frame(
        cluster_name = available_clusters[i],
        entry_name = group_names_full,
        entry_id = group_names,
        resolution = stringr::str_extract(available_clusters[i], "\\d+\\.*\\d*")
      )
    )
  }
  
  n2r <- df.feature.x$resolution
  names(n2r) <- df.feature.x$entry_name
  n2c <- df.feature.x$entry_id
  names(n2c) <- df.feature.x$entry_name
  cdi_cluster_top$resolution <- n2r[cdi_cluster_top$feature.x]
  cdi_cluster_top$cluster <- n2c[cdi_cluster_top$feature.x]
  
  ures <- unique(cdi_cluster_top$resolution )
  ures <- ures[order(ures)]
  
  auc_bins = cdi_bins
  qm.res.sum.all <- NULL
  for (j in 1:(length(auc_bins))){
    qm.res.sum <- cdi_cluster_top %>%
      dplyr::group_by(resolution) %>%
      dplyr::summarize(pdeg = sum((ncdi > auc_bins[j]))/length(ncdi),
                       bin = auc_bins[j], .groups = 'drop')
    qm.res.sum.all <- bind_rows(qm.res.sum.all, qm.res.sum)
    
  }
  
  df.cdi.spec <- NULL
  for (i in 1:length(ures)){
    qm.res.sum.sum.cur <- qm.res.sum.all %>% dplyr::filter(resolution %in% ures[i])
    x = qm.res.sum.sum.cur$bin
    y = qm.res.sum.sum.cur$pdeg
    id <- order(x)
    auc.spec <- sum(diff(x[id])*zoo::rollmean(y[id],2))
    
    df.cdi.spec <- bind_rows(df.cdi.spec,
                             data.frame(
                               res = ures[i],
                               auc = auc.spec
                             ))
  }
  
  df.cdi.spec2 <- df.cdi.spec
  
  plt.cdi.spec <- df.cdi.spec2 %>%
    ggplot(aes(x = as.numeric(res), y = auc)) +
    geom_point(size = 3) + geom_line() +
    theme_miko() +
    labs(x = "Resolution", y = "Specificity Score", title = "CDI Specificity Scores") +
    theme(panel.grid.minor = element_line(colour="grey95", size=0.1),
          panel.grid.major = element_line(colour="grey85", size=0.1),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank()) +
    scale_x_continuous(minor_breaks = seq(0 , max(df.cdi.spec2$res, na.rm = T), 0.1) ,
                       breaks = seq(0, max(df.cdi.spec2$res, na.rm = T), 0.2))
  
  plt.allClustSpec <- qm.res.sum.all %>%
    ggplot(aes(x = bin, y = pdeg, group = resolution , color = resolution)) +
    geom_line() + theme_miko(legend = T) +
    labs(x = "nCDI", y = "Fraction > nCDI Threshold") +
    labs(title = "Cluster Specificity Curves", color = "Resolution")
  
  if (return_dotplot){
    df.feature.x.unique <- unique(df.feature.x[ ,c("cluster_name", "resolution")])
    
    miko_message("Generating dot plots...", verbose = verbose)
    plt_cdi_dot.list <- mylapply(1:nrow(df.feature.x.unique), function(x){
      suppressWarnings({
        suppressMessages({
          ind <- x[[1]]
          
          plt_cdi_dot <- NULL
          try({
            whichres <-df.feature.x.unique$resolution[ind]
            whichclust <- df.feature.x.unique$cluster_name[ind]
            
            cdi_cluster_top2 <- cdi_cluster_top %>% dplyr::filter(resolution %in% whichres)
            cdi_cluster_top2$cluster <- as.numeric(as.character(cdi_cluster_top2$cluster))
            cdi_cluster_top2 <- cdi_cluster_top2 %>% dplyr::arrange(cluster)
            cdi_features <- unique(cdi_cluster_top2$feature.y)
            whichauc <-df.cdi.spec2$auc[df.cdi.spec2$res == whichres]
            plt_cdi_dot <- DotPlot(object = so.query, features = cdi_features, group.by = whichclust) +
              scale_color_gradient2(high = scales::muted("red"), low = scales::muted("blue")) +
              theme(legend.position = "bottom") +
              labs(x = "Genes", y = "Cluster", title = "Top Cluster-Specific Markers",
                   subtitle = paste0("Resolution = ",
                                     whichres, ", Specificity Score = ", signif(whichauc, 3)))
          }, silent= T)
          
          return(plt_cdi_dot)
          
        })
      })
    })
    
    names(plt_cdi_dot.list) <- as.character(df.feature.x.unique$resolution)
    
  } else {
    plt_cdi_dot.list <- NULL
  }
  
  
  return(
    list(
      specificity_summary = df.cdi.spec2,
      specificity_raw = qm.res.sum.all,
      auc_plot = plt.allClustSpec,
      resolution_plot = plt.cdi.spec,
      dot_plot = plt_cdi_dot.list,
      cdi_results = cdi_cluster
    )
  )
  
}

