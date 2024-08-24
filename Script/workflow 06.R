
## Nastasia Gauffeny - nastasia.gauffeny@live.com 

### Pipeline scRNA-seq : binarisation VS méthodes classiques

# 6. Clustering sur données intégrées ####

library(Seurat)
library(ggplot2)
library(ggforce)
library(tidyverse)
library(patchwork)
library(scMiko)
library(viridis)
library(gghalves)
library(radiant.data)
library(scSHC)

# Métriques clustering
library(mclust)
library(clevr)

# Source fonctions scMiko modifiées
for(i in list.files("./Script/Functions")) {
  source(sprintf("./Script/Functions/%s",i))
}
rm(i)

BP_count <- readRDS("./Output/Data/05-BP_count.rds")
BP_bin <- readRDS("./Output/Data/05-BP_bin.rds")
BP_set <- readRDS("./Output/Data/05-BP_set.rds")
custom_colors <- readRDS("./Output/Data/01-custom_colors.rds")

## Clusters de Seurat : count ####

BP_count <- FindNeighbors(BP_count, reduction = "umap", 
                          dims = 1:2, annoy.metric = "manhattan")
BP_count <- FindClusters(BP_count, algorithm = 1, resolution = c(0.01, 0.05, 0.1, 0.2))

## Clusters de Seurat : bin ####

BP_bin <- FindNeighbors(BP_bin, reduction = "umap", 
                        dims = 1:2, annoy.metric = "manhattan")
BP_bin <- FindClusters(BP_bin, algorithm = 1, resolution = c(0.01, 0.05, 0.1, 0.2))

## Choix de la résolution ####

choix_res <- function(seurat_obj, res) {
  
  cat("Tableau des résolutions et nombre de clusters obtenus :\n\n")
  print(sapply(grep("^RNA_snn_", colnames(seurat_obj@meta.data), value = TRUE), 
               function(x) length(unique(seurat_obj@meta.data[, x]))))
  cat("\n")
  
  res_col <- paste0("RNA_snn_res.", res)
  
  Idents(seurat_obj) <- res_col
  seurat_obj@meta.data$seurat_clusters <- seurat_obj@meta.data[[res_col]]
  
  cat("Résolution appliquée dans meta.data$seurat_clusters :", res, "\n")
  
  return(seurat_obj)
}

BP_count <- choix_res(BP_count, 0.01)
BP_bin <- choix_res(BP_bin, 0.01)

## Métriques d'évaluation du clustering ####

# Comptage / binaire
table("clusters count" = BP_count@meta.data$seurat_clusters,
      "clusters bin" = BP_bin@meta.data$seurat_clusters)

adjustedRandIndex(BP_count@meta.data$seurat_clusters,
                  BP_bin@meta.data$seurat_clusters)  
# res 0.01 : 0.9391773
# res 0.05 : 0.9250849
# res 0.1 : 0.9473127
# res 0.2 : 0.5910152

eval_report_clusters(BP_count@meta.data$seurat_clusters,
                     BP_bin@meta.data$seurat_clusters)

# Comptage / réel
table("clusters count" = BP_count@meta.data$seurat_clusters,
      "labels" = BP_set$samplesheet$celltype)

adjustedRandIndex(BP_count@meta.data$seurat_clusters,
                  BP_set$samplesheet$celltype)  
# res 0.01 : 0.9402631
# res 0.05 : 0.9402631
# res 0.1 : 0.9512421
# res 0.2 : 0.6595537

eval_report_clusters(BP_count@meta.data$seurat_clusters,
                     BP_set$samplesheet$celltype)

# Binaire / réel
table("clusters bin" = BP_bin@meta.data$seurat_clusters,
      "labels" = BP_set$samplesheet$celltype)

adjustedRandIndex(BP_bin@meta.data$seurat_clusters,
                  BP_set$samplesheet$celltype)  
# res 0.01 : 0.910417
# res 0.05 : 0.9253429
# res 0.1 : 0.9354802
# res 0.2 : 0.6345641

eval_report_clusters(BP_bin@meta.data$seurat_clusters,
                     BP_set$samplesheet$celltype)

## Alluvial plot ####

plot_Alluvial <- function(seurat_obj) {
  
  # 1. Samples by cluster
  groups <- levels(as.factor(seurat_obj@meta.data$group))
  clusters <- levels(seurat_obj@meta.data$seurat_clusters)
  
  color_assignments <- setNames(
    c(custom_colors$discrete[1:length(groups)], custom_colors$discrete[1:length(clusters)]),
    c(groups,clusters))
  
  data <- seurat_obj@meta.data %>%
    group_by(group, seurat_clusters) %>%
    tally() %>%
    ungroup() %>%
    gather_set_data(1:2) %>%
    dplyr::mutate(
      x = factor(x, levels = unique(x)),
      y = factor(y, levels = unique(y))
    )
  
  # hjust defines whether a label will be aligned to the right (1) or to the left (0); 
  # the nudge_x parameter is used to move the label outside of the boxes
  data_labels <- tibble(
    group = c(
      rep('group', length(groups)),
      rep('seurat_clusters', length(clusters))
    )
  ) %>%
    mutate(
      hjust = ifelse(group == 'group', 1, 0),
      nudge_x = ifelse(group == 'group', -0.1, 0.1)
    )
  
  p1 <- ggplot(data, aes(x, id = id, split = y, value = n)) +
    geom_parallel_sets(aes(fill = seurat_clusters), alpha = 0.75, axis.width = 0.15) +
    geom_parallel_sets_axes(aes(fill = y), color = 'black', axis.width = 0.1) +
    geom_text(
      aes(y = n, split = y), stat = 'parallel_sets_axes', fontface = 'bold',
      hjust = data_labels$hjust, nudge_x = data_labels$nudge_x
    ) +
    scale_x_discrete(labels = c('Group','Cluster'), 
                     expand = expansion(add = c(0.7, 0.2))) +
    scale_fill_manual(values = color_assignments) +
    theme_bw() +
    theme(
      legend.position = 'none',
      axis.title = element_blank(),
      axis.text.x = element_text(face = 'bold', colour = 'black', size = 15),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank()
    )
  
  # 2. Clusters by cell types
  clusters <- levels(seurat_obj@meta.data$seurat_clusters)
  cell_types <- as.vector(sort(unique(seurat_obj@meta.data$label)))
  
  color_assignments <- setNames(
    c(custom_colors$discrete[1:length(clusters)], custom_colors$discrete[1:length(cell_types)]),
    c(clusters,cell_types)
  )
  
  data <- seurat_obj@meta.data %>%
    dplyr::rename(cell_type = label) %>%
    dplyr::mutate(cell_type = factor(cell_type, levels = cell_types)) %>%
    group_by(seurat_clusters, cell_type) %>%
    tally() %>%
    ungroup() %>%
    gather_set_data(1:2) %>%
    dplyr::mutate(
      x = factor(x, levels = unique(x)),
      y = factor(y, levels = c(clusters, cell_types))
    )
  
  data_labels <- tibble(
    group = c(
      rep('seurat_clusters', length(clusters)),
      rep('cell_type', length(cell_types))
    )
  ) %>%
    mutate(
      hjust = ifelse(group == 'seurat_clusters', 1, 0),
      nudge_x = ifelse(group == 'seurat_clusters', -0.1, 0.1)
    )
  
  p2 <- ggplot(data, aes(x, id = id, split = y, value = n)) +
    geom_parallel_sets(aes(fill = seurat_clusters), alpha = 0.75, axis.width = 0.15) +
    geom_parallel_sets_axes(aes(fill = y), color = 'black', axis.width = 0.1) +
    geom_text(
      aes(y = n, split = y), stat = 'parallel_sets_axes', fontface = 'bold',
      hjust = data_labels$hjust, nudge_x = data_labels$nudge_x
    ) +
    scale_x_discrete(labels = c('Cluster','Cell type'), 
                     expand = expansion(add = c(0.2, 0.9))) +
    scale_fill_manual(values = color_assignments) +
    theme_bw() +
    theme(
      legend.position = 'none',
      axis.title = element_blank(),
      axis.text.x = element_text(face = 'bold', colour = 'black', size = 15),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank()
    )
  
  # 3. Combiner
  print(p1 + p2 + plot_layout(ncol = 2))
  return(p1 + p2 + plot_layout(ncol = 2))
  
}

plot1 <- plot_Alluvial(BP_count)
ggsave("Output/Plots/06-AlluvialPlot_count.png", plot1, height = 6, width = 8)

plot2 <- plot_Alluvial(BP_bin)
ggsave("Output/Plots/06-AlluvialPlot_bin.png", plot2, height = 6, width = 8)

## Visualisations sur UMAP des données intégrées ####

## Représentation UMAP des clusters (scMiko)

plot_Label_clusters <- function(seurat_obj, output_file) {
  
  plot1 <- scMiko::cluster.UMAP(so = seurat_obj, group.by = "label", repel = TRUE) + 
    theme(legend.key.size = unit(0.1, 'in'))
  plot2 <- scMiko::cluster.UMAP(so = seurat_obj, group.by = "seurat_clusters") + 
    theme(legend.key.size = unit(0.1, 'in'))
  
  print(plot1 + plot2)
  ggsave(output_file, plot1 + plot2, width = 10, height = 5)
}

plot_Label_clusters(BP_count, "./Output/Plots/06-Label_clusters_count.png")
plot_Label_clusters(BP_bin, "./Output/Plots/06-Label_clusters_bin.png")

## Représentation UMAP complète (Roman Hillje)

plot_UMAP_hillje <- function(seurat_obj) {
  
  plot_umap_by_nCount <- bind_cols(seurat_obj@meta.data, as.data.frame(seurat_obj@reductions[["umap"]]@cell.embeddings)) %>%
    ggplot(aes(umap_1, umap_2, color = nCount_RNA)) +
    geom_point(size = 0.2) +
    theme_bw() +
    scale_color_viridis(
      guide = guide_colorbar(frame.colour = 'black', ticks.colour = 'black'),
      labels = scales::comma,
    ) +
    labs(color = 'Number of\ntranscripts') +
    theme(legend.position = 'left') +
    coord_fixed() +
    annotate(
      geom = 'text', x = Inf, y = -Inf,
      label = paste0('n = ', format(nrow(seurat_obj@meta.data), big.mark = ',', trim = TRUE)),
      vjust = -1.5, hjust = 1.25, color = 'black', size = 2.5
    )
  
  plot_umap_by_sample <- bind_cols(seurat_obj@meta.data, as.data.frame(seurat_obj@reductions[["umap"]]@cell.embeddings)) %>%
    ggplot(aes(umap_1, umap_2, color = group)) +
    geom_point(size = 0.2) +
    theme_bw() +
    scale_color_manual(values = custom_colors$discrete) +
    labs(color = 'Sample') +
    guides(colour = guide_legend(override.aes = list(size = 2))) +
    theme(legend.position = 'right') +
    coord_fixed() +
    annotate(
      geom = 'text', x = Inf, y = -Inf,
      label = paste0('n = ', format(nrow(seurat_obj@meta.data), big.mark = ',', trim = TRUE)),
      vjust = -1.5, hjust = 1.25, color = 'black', size = 2.5
    )
  
  plot_umap_by_cluster <- bind_cols(seurat_obj@meta.data, as.data.frame(seurat_obj@reductions[["umap"]]@cell.embeddings)) %>%
    ggplot(aes(umap_1, umap_2, color = seurat_clusters)) +
    geom_point(size = 0.2) +
    theme_bw() +
    scale_color_manual(
      name = 'Cluster', values = custom_colors$discrete,
      guide = guide_legend(ncol = 2, override.aes = list(size = 2))
    ) +
    theme(legend.position = 'left') +
    coord_fixed() +
    annotate(
      geom = 'text', x = Inf, y = -Inf,
      label = paste0('n = ', format(nrow(seurat_obj@meta.data), big.mark = ',', trim = TRUE)),
      vjust = -1.5, hjust = 1.25, color = 'black', size = 2.5
    )
  
  plot_umap_by_label <- bind_cols(seurat_obj@meta.data, as.data.frame(seurat_obj@reductions[["umap"]]@cell.embeddings)) %>%
    ggplot(aes(umap_1, umap_2, color = label)) +
    geom_point(size = 0.2) +
    theme_bw() +
    scale_color_manual(
      name = 'Label', values = custom_colors$discrete,
      guide = guide_legend(ncol = 2, override.aes = list(size = 2))
    ) +
    theme(legend.position = 'right') +
    coord_fixed() +
    annotate(
      geom = 'text', x = Inf, y = -Inf,
      label = paste0('n = ', format(nrow(seurat_obj@meta.data), big.mark = ',', trim = TRUE)),
      vjust = -1.5, hjust = 1.25, color = 'black', size = 2.5
    )
  
  return(plot_umap_by_nCount + plot_umap_by_sample + 
           plot_umap_by_cluster + plot_umap_by_label +
           plot_layout(ncol = 2))
  
}

plots_umap_count <- plot_UMAP_hillje(BP_count)
ggsave("Output/Plots/06-UMAPs_intég_count.png", plots_umap_count,
       height = 6, width = 10.5)

plots_umap_bin <- plot_UMAP_hillje(BP_bin)
ggsave("Output/Plots/06-UMAPs_intég_bin.png", plots_umap_bin,
       height = 6, width = 10.5)

## Composition des clusters en donneurs

plot_Composition_samples_clusters <- function(seurat_obj, output_file) {
  
  table_samples_by_clusters <- seurat_obj@meta.data %>%
    group_by(group, seurat_clusters) %>%
    summarize(count = n()) %>%
    spread(seurat_clusters, count, fill = 0) %>%
    ungroup() %>%
    mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
    dplyr::select(c('group', 'total_cell_count', everything())) %>%
    arrange(factor(group, levels = levels(seurat_obj@meta.data$group)))
  
  # knitr::kable(table_samples_by_clusters)
  
  temp_labels <- seurat_obj@meta.data %>%
    group_by(group) %>%
    tally()
  
  p1 <- table_samples_by_clusters %>%
    select(-c('total_cell_count')) %>%
    reshape2::melt(id.vars = 'group') %>%
    mutate(group = factor(group, levels = levels(seurat_obj@meta.data$group))) %>%
    ggplot(aes(group, value)) +
    geom_bar(aes(fill = variable), position = 'stack', stat = 'identity') +
    geom_text(
      data = temp_labels,
      aes(x = group, y = Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1),
      color = 'black', size = 2.8
    ) +
    scale_fill_manual(name = 'Cluster', values = custom_colors$discrete) +
    scale_y_continuous(name = 'Number of cells', labels = scales::comma, expand = c(0.01, 0)) +
    coord_cartesian(clip = 'off') +
    theme_bw() +
    theme(
      legend.position = 'left',
      plot.title = element_text(hjust = 0.5),
      text = element_text(size = 16),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      plot.margin = margin(t = 20, r = 0, b = 0, l = 0, unit = 'pt')
    )
  
  table_clusters_by_samples <- seurat_obj@meta.data %>%
    dplyr::rename('cluster' = 'seurat_clusters') %>%
    group_by(cluster, group) %>%
    summarize(count = n()) %>%
    spread(group, count, fill = 0) %>%
    ungroup() %>%
    mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
    select(c('cluster', 'total_cell_count', everything())) %>%
    arrange(factor(cluster, levels = levels(seurat_obj@meta.data$seurat_clusters)))
  
  # knitr::kable(table_clusters_by_samples)
  
  temp_labels <- seurat_obj@meta.data %>%
    group_by(seurat_clusters) %>%
    tally() %>%
    dplyr::rename('cluster' = seurat_clusters)
  
  p2 <- table_clusters_by_samples %>%
    select(-c('total_cell_count')) %>%
    reshape2::melt(id.vars = 'cluster') %>%
    mutate(cluster = factor(cluster, levels = levels(seurat_obj@meta.data$seurat_clusters))) %>%
    ggplot(aes(cluster, value)) +
    geom_bar(aes(fill = variable), position = 'stack', stat = 'identity') +
    geom_text(
      data = temp_labels,
      aes(x = cluster, y = Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1),
      color = 'black', size = 2.8
    ) +
    scale_fill_manual(name = 'Sample', values = custom_colors$discrete) +
    scale_y_continuous(labels = scales::comma, expand = c(0.01, 0)) +
    coord_cartesian(clip = 'off') +
    theme_bw() +
    theme(
      legend.position = 'right',
      plot.title = element_text(hjust = 0.5),
      text = element_text(size = 16),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_blank(),
      plot.margin = margin(t = 20, r = 0, b = 0, l = 10, unit = 'pt')
    )
  
  plot <- p1 + p2 +
    plot_layout(ncol = 2, widths = c(
      seurat_obj@meta.data$group %>% unique() %>% length(),
      seurat_obj@meta.data$seurat_clusters %>% unique() %>% length()
    ))
  
  print(plot)
  
  ggsave(output_file, plot, width = 18, height = 8)
  
}

plot_Composition_samples_clusters(BP_count, "./Output/Plots/06-composition_samples_clusters_count.png")
plot_Composition_samples_clusters(BP_bin, "./Output/Plots/06-composition_samples_clusters_bin.png")

## Nombre de transcrits et gènes exprimés

plot_Nb_transcrits_features <- function(seurat_obj, output_file) {
  
  temp_labels <- seurat_obj@meta.data %>%
    group_by(seurat_clusters) %>%
    tally()
  
  p1 <- ggplot() +
    geom_half_violin(
      data = seurat_obj@meta.data, aes(seurat_clusters, nCount_RNA, fill = seurat_clusters),
      side = 'l', show.legend = FALSE, trim = FALSE
    ) +
    geom_half_boxplot(
      data = seurat_obj@meta.data, aes(seurat_clusters, nCount_RNA, fill = seurat_clusters),
      side = 'r', outlier.color = NA, width = 0.4, show.legend = FALSE
    ) +
    geom_text(
      data = temp_labels,
      aes(x = seurat_clusters, y = -Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1),
      color = 'black', size = 2.8
    ) +
    scale_color_manual(values = custom_colors$discrete) +
    scale_fill_manual(values = custom_colors$discrete) +
    scale_y_continuous(name = 'Number of transcripts', labels = scales::comma, expand = c(0.08, 0)) +
    theme_bw() +
    theme(
      panel.grid.major.x = element_blank(),
      axis.title.x = element_blank()
    )
  
  p2 <- ggplot() +
    geom_half_violin(
      data = seurat_obj@meta.data, aes(seurat_clusters, nFeature_RNA, fill = seurat_clusters),
      side = 'l', show.legend = FALSE, trim = FALSE
    ) +
    geom_half_boxplot(
      data = seurat_obj@meta.data, aes(seurat_clusters, nFeature_RNA, fill = seurat_clusters),
      side = 'r', outlier.color = NA, width = 0.4, show.legend = FALSE
    ) +
    geom_text(
      data = temp_labels,
      aes(x = seurat_clusters, y = -Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1),
      color = 'black', size = 2.8
    ) +
    scale_color_manual(values = custom_colors$discrete) +
    scale_fill_manual(values = custom_colors$discrete) +
    scale_y_continuous(name = 'Number of expressed genes', labels = scales::comma, expand = c(0.08, 0)) +
    theme_bw() +
    theme(
      panel.grid.major.x = element_blank(),
      axis.title.x = element_blank()
    )
  
  plot <- p1 + p2 + plot_layout(ncol = 1)
  print(plot)
  
  ggsave(output_file, plot, height = 7, width = 14)
  
}

plot_Nb_transcrits_features(BP_count, "./Output/Plots/06-Transcripts_features_count.png")
plot_Nb_transcrits_features(BP_bin, "./Output/Plots/06-Transcripts_features_bin.png")

## Similarité entre clusters

plot_clusters_similarity <- function(seurat_obj, output_file) {
  
  require("SingleCellExperiment")
  
  sce <- as.SingleCellExperiment(seurat_obj)
  reducedDim(sce, 'PCA_sub') <- reducedDim(sce, 'PCA')[,1:15, drop = FALSE]
  
  g <- scran::buildSNNGraph(sce, use.dimred = 'PCA_sub')
  ratio <- bluster::pairwiseModularity(g, seurat_obj@meta.data$seurat_clusters, as.ratio = TRUE)
  ratio_to_plot <- log10(ratio+1)
  
  plot <- ratio_to_plot %>%
    as_tibble() %>%
    rownames_to_column(var = 'cluster_1') %>%
    pivot_longer(
      cols = 2:ncol(.),
      names_to = 'cluster_2',
      values_to = 'probability'
    ) %>%
    mutate(
      cluster_1 = as.character(as.numeric(cluster_1) - 1),
      cluster_1 = factor(cluster_1, levels = rev(unique(cluster_1))),
      cluster_2 = factor(cluster_2, levels = unique(cluster_2))
    ) %>%
    ggplot(aes(cluster_2, cluster_1, fill = probability)) +
    geom_tile(color = 'white') +
    geom_text(aes(label = round(probability, digits = 2)), size = 2.5) +
    scale_x_discrete(name = 'Cluster', position = 'top') +
    scale_y_discrete(name = 'Cluster') +
    scale_fill_gradient(
      name = 'log10(ratio)', low = 'white', high = '#c0392b', na.value = '#bdc3c7',
      guide = guide_colorbar(
        frame.colour = 'black', ticks.colour = 'black', title.position = 'left',
        title.theme = element_text(hjust = 1, angle = 90),
        barwidth = 0.75, barheight = 10
      )
    ) +
    coord_fixed() +
    theme_bw() +
    theme(
      legend.position = 'right',
      panel.grid.major = element_blank()
    )
  
  print(plot)
  ggsave(output_file, plot, height = 6, width = 7)
  
}

plot_clusters_similarity(BP_count, "./Output/Plots/06-cluster_sim_count.png")
plot_clusters_similarity(BP_bin, "./Output/Plots/06-cluster_sim_bin.png")

## BuildClusterTree : classification hiérarchique des clusters

plot_cluster_tree <- function(seurat_obj, output_file) {
  
  clusters <- as.data.frame(Idents(seurat_obj))
  
  seurat_obj <- BuildClusterTree(
    seurat_obj,
    dims = 1:15,
    reorder = TRUE,
    reorder.numeric = TRUE
  )
  
  tree <- seurat_obj@tools$BuildClusterTree
  tree$tip.label <- paste0("Cluster ", tree$tip.label)
  
  plot <- ggtree::ggtree(tree, aes(x, y)) +
    scale_y_reverse() +
    ggtree::geom_tree() +
    ggtree::theme_tree() +
    ggtree::geom_tiplab(offset = 1) +
    ggtree::geom_tippoint(color = custom_colors$discrete[1:length(tree$tip.label)], shape = 16, size = 5) +
    ggtree::geom_nodelab(aes(label=node), hjust = - 0.3) + 
    coord_cartesian(clip = 'off') +
    ggtree::geom_nodepoint() +
    theme(plot.margin = unit(c(0,2.5,0,0), 'cm'))
  
  print(plot)
  ggsave(output_file, plot)
  
}

plot_cluster_tree(BP_count, "./Output/Plots/06-clustree_count.png")
plot_cluster_tree(BP_bin, "./Output/Plots/06-clustree_bin.png")

## Silhouette et spécificité CDI (scMiko) ####

run_scMiko <- function(seurat_obj, output_plot_file,
                       list_res = c(0.01, 0.05, 0.1, 0.2)) {
  
  # Clustering multi-résolution
  mc.list <- multiCluster(object = seurat_obj, 
                          resolutions = list_res, 
                          nworkers = 4, 
                          pca_var = 0.9, 
                          group_singletons = T, 
                          algorithm = 1,
                          return_object = F)
  
  plt.umap_by_cluster <- mc.list$plots
  so.query <- mc.list$object
  cr_names <- mc.list$resolution_names
  cluster.name <- mc.list$cluster_names
  assay.pattern <- mc.list$assay_pattern
  
  # Visualisation des configurations de cluster
  plot1 <- cowplot::plot_grid(plotlist = lapply(plt.umap_by_cluster, function(x) {
    x + theme_void() + labs(title = NULL) +
      theme(legend.position = "none", plot.subtitle = element_text(hjust = 0.5))
  }), ncol = 4)
  
  print(plot1)
  ggsave(paste0(output_plot_file, "_cluster_config.png"), plot1,
         bg = "white", height = 6, width = 12)
  
  # Critère de sélection de la résolution basé sur la spécificité
  ms.list <- multiSpecificity_m(object = so.query,
                                cluster_names = cluster.name,
                                features = NULL,
                                deg_prefilter = T,
                                cdi_bins = seq(0, 1, by = 0.01),
                                min.pct = 0.1,
                                n.workers = 4,
                                return_dotplot = T,
                                verbose = T)
  
  df.summary <- ms.list$specificity_summary
  df.raw <- ms.list$specificity_raw
  plt.specificity.umap <- ms.list$umap_plot
  plt.clust.spec <- ms.list$auc_plot
  plt.auc.spec <- ms.list$resolution_plot
  plt.auc.dot <- ms.list$dot_plot
  
  max.auc = max(df.summary$auc)
  speak <- df.summary$res[which.max(df.summary$auc)]
  
  # Visualusation des graphs de spécificité selon la résolution
  plot2 <- cowplot::plot_grid(plt.auc.spec + geom_hline(yintercept = max.auc, linetype = "dashed", color = "tomato") +
                                geom_vline(xintercept = as.numeric(speak), linetype = "dashed", color = "tomato"),
                              plt.clust.spec,
                              nrow = 1, rel_widths = c(1, 1), labels = "AUTO")
  
  print(plot2)
  ggsave(paste0(output_plot_file, "_cluster_specificity.png"), plot2,
         bg = "white", height = 4, width = 6)
  
  # Plot AUC selon la resolution : "Top Cluster-Specific Markers"
  # modification manuelle possible dans la fonction de la resolution 
  # au niveau de plt.auc.dot[["0.05"]] du plot3
  
  plt.auc.dot$'0.01'
  plt.auc.dot$'0.05'
  plt.auc.dot$'0.1'
  plt.auc.dot$'0.2'
  
  # plot3 <- cowplot::plot_grid(plt.auc.dot[["0.05"]] + 
  #                               theme(legend.position = "right", 
  #                                     axis.text.x = element_text(angle = 45, hjust = 1)), 
  #                             nrow = 1, labels = "C")
  
  # print(plot3)
  # ggsave(paste0(output_plot_file, "_clusters_top_markers.png"), plot3,
  #        bg = "white", height = 4, width = 6)
  
  
  # Largeur de silhouette moyenne selon la résolution
  msil_list <- multiSilhouette(object = so.query, groups = cluster.name,
                               assay_pattern = assay.pattern, verbose = T)
  
  # Visualisation du score de silhouette selon la résolution
  plot4 <- cowplot::plot_grid(plotlist = lapply(msil_list$silhouette_plots, function(x) {
    x + theme(plot.subtitle = element_text(hjust = 0.5),
              plot.title = element_text(size=12),
              legend.key.size = unit(0.1, 'in'),
              legend.position="bottom",
              legend.position.inside = c(0, 0.5)) }), nrow = 2)
  # legend.title = element_text(size=12)
  
  print(plot4)
  ggsave(paste0(output_plot_file, "_silhouette.png"), plot4, 
         bg = "white", height = 6, width = 12)
  
}

run_scMiko(BP_count, "./Output/Plots/06-scmiko_count")
run_scMiko(BP_bin, "./Output/Plots/06-scmiko_bin")

## Tableaux récap résolution / label ####

janitor::tabyl(BP_count@meta.data, seurat_clusters, label)
janitor::tabyl(BP_bin@meta.data, seurat_clusters, label)

table(BP_count@meta.data$seurat_clusters, 
      BP_bin@meta.data$seurat_clusters)

# Fin partie 06

saveRDS(BP_count, "./Output/Data/06-BP_count.rds")
rm(BP_count)

saveRDS(BP_bin, "./Output/Data/06-BP_bin.rds")
rm(BP_bin)

saveRDS(BP_set, "./Output/Data/06-BP_set.rds")
rm(BP_set)

rm(custom_colors)
rm(plot1)
rm(plot2)
rm(plots_umap_bin)
rm(plots_umap_count)

gc()

#
