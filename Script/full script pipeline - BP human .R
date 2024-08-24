### Pipeline propre scRNA-seq : binaire VS count

# Changement de résolution : 
# - 04 - Clustering sur données non intégrées : choix_res
# - 06 - Clustering sur données intégrées : choix_res
# - 07 - Expression différentielle : perform_de_analysis
# - 07 - Expression différentielle : visualize_markers

# 1. Préparation jeu de données ####

custom_colors <- list()
custom_colors$discrete <- c(
  '#FFC312','#C4E538','#12CBC4','#FDA7DF','#ED4C67',
  '#F79F1F','#A3CB38','#1289A7','#D980FA','#B53471',
  '#EE5A24','#009432','#0652DD','#9980FA','#833471',
  '#EA2027','#006266','#1B1464','#5758BB','#6F1E51',
  '#40407a','#706fd3','#f7f1e3','#34ace0','#33d9b2',
  '#2c2c54','#474787','#aaa69d','#227093','#218c74',
  '#ff5252','#ff793f','#d1ccc0','#ffb142','#ffda79',
  '#b33939','#cd6133','#84817a','#cc8e35','#ccae62')

## Librairies

library(scRNAseq)
library(Seurat)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(dplyr)
library(forcats)

# Quality Control (QC)
library(BPCells)
library(data.table)
library(purrr)

## Données de comptage : count ####

BP_count <- BaronPancreasData('human')
BP_count <- as.Seurat(BP_count, data = NULL) 
# dim 20125  8569

## Filtre & contrôle qualité (QC) ####

BP_count <- subset(BP_count,
                   subset = nFeature_originalexp > 200
                   & nFeature_originalexp < 2500) 
# dim 20125  7242

BP_count <- AddMetaData(BP_count, BP_count@meta.data$donor, col.name = "group")

plots_QC <- function(seurat_obj, show_feature_scatter = TRUE, y1 = 200, y2 = 2500,
                     feature1 = "nCount_originalexp", feature2 = "nFeature_originalexp", 
                     group.by = "group", show_knee_plot = TRUE, show_vln_plot = TRUE) {
  
  if (show_feature_scatter) {
    plot <- FeatureScatter(seurat_obj, feature1 = feature1, feature2 = feature2, group.by = group.by)
    plot + geom_hline(yintercept = y1, linetype = "dashed") + geom_hline(yintercept = y2, linetype = "dashed")
    print(plot)
    
    ggsave("./Output/Plots/01-QC_FeatureScatter.png", plot, bg = "white")
  }
  
  if (show_knee_plot) {
    plot <- plot_read_count_knee(Matrix::colSums(seurat_obj@assays$originalexp$counts))
    print(plot)
    
    # Comparaison : plot attendu
    # https://bnprks.github.io/BPCells/articles/pbmc3k_files/figure-html/unnamed-chunk-10-1.png
    
    ggsave("./Output/Plots/01-QC_KneePlot.png", plot)
  }
  
  if (show_vln_plot) {
    plot <- VlnPlot(seurat_obj, features = c("nCount_originalexp", "nFeature_originalexp"),
                    ncol = 2) + coord_cartesian(ylim = c(200, 3000)) &
      theme(axis.title.x=element_blank(), 
            axis.text.x=element_blank(), 
            axis.ticks.x=element_blank()) 
    print(plot)
    ggsave("./Output/Plots/01-QC_VlnPlot.png", plot)
  }
}

plots_QC(BP_count)
# FeatureScatter : corrélation de Pearson

## Binarisation de la matrice de comptage : bin ####

BP_bin <- BP_count

binarize <- function(seurat_obj, assay = "originalexp", layer = "counts") {
  extrait_seurat <- LayerData(seurat_obj, assay = assay, layer = layer)
  x <- extrait_seurat@x
  x[x != 1] <- 1
  extrait_seurat@x <- x
  LayerData(seurat_obj, assay = assay, layer = layer) <- extrait_seurat 
  return(seurat_obj)
}

BP_bin <- binarize(BP_bin)

## BP_set ####

BP_set <- NULL
BP_set$binary <- BP_bin@assays$originalexp@counts
BP_set$count <- BP_count@assays$originalexp@counts
BP_samplesheet <- data.frame(celltype = as.factor(BP_count@meta.data$label),
                             colnames = as.factor(BP_count@assays$originalexp@counts@Dimnames[[2]]),
                             group = as.factor(BP_count@meta.data$donor),
                             labelNUM = as.numeric(as.factor(BP_count@meta.data$label)))
BP_set$samplesheet <- BP_samplesheet

## Pieplot de la répartition des vrais labels
# source : https://stackoverflow.com/a/74931586 

plot_Label <- function(df = BP_set$samplesheet, dataset = "", max_x = 40) {
  
  table_freq <- table(df$celltype)
  pie_data <- data.frame(Label = names(table_freq), Count = as.numeric(table_freq))
  pie_data$Percentage <- round((pie_data$Count *100/ sum(pie_data$Count)), 1)
  pie_data$Label <- factor(pie_data$Label, 
                           levels = pie_data$Label[order(pie_data$Percentage, decreasing = TRUE)])
  pie_data$LegendLabels <- paste0(pie_data$Label, " (", pie_data$Count, ")")
  
  plot <- pie_data |>
    mutate(Label = fct_reorder(Label, Percentage, mean)) |>
    ggplot(aes(x = Percentage, y = Label, fill = Label)) +
    geom_bar(aes(x = max_x), stat = "identity", fill = '#e0e0e0') +
    geom_bar(stat = "identity") +
    geom_text(aes(x = Percentage, label = paste(sprintf("%2.1f", round(Percentage, 1)), "%")), 
              adj = 0, nudge_x = 1, size = 3) +
    facet_wrap(~ Percentage < median(Percentage), scales = 'free_y', ncol = 2) +
    scale_x_continuous(breaks = c(0, 10, 20, 30, 40, 50, 60)) +
    ggtitle(paste0("Répartition des vrais types cellulaires\n", dataset, "\n")) +
    labs(x = '', y = '') +
    theme(panel.background = element_blank(),
          strip.text.x = element_blank(), 
          axis.text.y = element_text(face = "bold"),
          panel.grid = element_blank(),
          legend.position = "bottom",
          legend.key.size = unit(0.09, "inch"),
          plot.title = element_text(hjust = 0.5),
          legend.title = element_blank(),
          legend.location = "plot", 
          plot.title.position = "plot"
    ) +
    scale_fill_manual(values = custom_colors$discrete,
                      breaks = levels(pie_data$Label),
                      labels = pie_data$LegendLabels)
  print(plot)
  
  ggsave("./Output/Plots/01-Labels.png", plot, bg = "white",
         height = 4, width = 6)
  
}

plot_Label(dataset = "Baron Pancreas Human")

# Fin partie 01

saveRDS(BP_set, "./Output/Data/01-BP_set.rds")
rm(BP_set)

saveRDS(custom_colors, "./Output/Data/01-custom_colors.rds")
rm(custom_colors)

rm(BP_count)
rm(BP_bin)
rm(BP_samplesheet)

gc()

# 2. Réduction linéaire de dimensions ####

library(Seurat)
library(ggplot2)
library(patchwork)
library(cowplot) # ?
library(ggpubr)

BP_set <- readRDS("./Output/Data/01-BP_set.rds")

## ACP de Seurat : count ####

BP_count <- CreateSeuratObject(counts = BP_set$count, project = "BP")
BP_count <- NormalizeData(BP_count)
BP_count <- FindVariableFeatures(BP_count, selection.method = "vst", nfeatures = 2000)

BP_count <- ScaleData(BP_count, features = VariableFeatures(BP_count))
BP_count <- RunPCA(BP_count,features = VariableFeatures(BP_count))

BP_count <- AddMetaData(BP_count, BP_set$samplesheet$group, col.name = "group")
BP_count <- AddMetaData(BP_count, BP_set$samplesheet$celltype, col.name = "label")

## ACP de Seurat : bin ####

BP_bin <- CreateSeuratObject(counts = BP_set$binary, project = "BP")
BP_bin <- NormalizeData(BP_bin)
BP_bin <- FindVariableFeatures(BP_bin, selection.method = "vst", nfeatures = 2000)

BP_bin <- ScaleData(BP_bin, features = VariableFeatures(BP_count))
BP_bin <- RunPCA(BP_bin,features =  VariableFeatures(BP_count))

BP_bin <- AddMetaData(BP_bin, BP_set$samplesheet$group, col.name = "group")
BP_bin <- AddMetaData(BP_bin, BP_set$samplesheet$celltype, col.name = "label")

## Plots ACP ####

## Variable Feature Plot

plot_VariableFeature <- function(top_label = 10) {
  
  cat("Top", top_label, "variable features:\n",
      "Count: ", head(VariableFeatures(BP_count), 10), "\n")
  # cat("Binaire: ", head(VariableFeatures(BP_bin, 10)), "\n")
  
  plot1 <- VariableFeaturePlot(BP_count) + ggtitle("Variable Feature Plot - count") +
    theme(legend.position = "top")
  
  # plot2 <- VariableFeaturePlot(BP_bin) + ggtitle("Variable Feature Plot - bin")
  # plot <- ggarrange(plot1, plot2, ncol=2, nrow=1, common.legend = TRUE, legend="bottom")
  # print(plot)
  
  print(plot1)
  ggsave("./Output/Plots/02-VariableFeaturePlot.png", plot1, bg = "white",
         height = 4, width = 6)
}

plot_VariableFeature()

## Représentation des observations

BP_set$samplesheet$count_PC1 <- BP_count@reductions$pca@cell.embeddings[,1]
BP_set$samplesheet$count_PC2 <- BP_count@reductions$pca@cell.embeddings[,2]
BP_set$samplesheet$bin_PC1 <- BP_bin@reductions$pca@cell.embeddings[,1]
BP_set$samplesheet$bin_PC2 <- BP_bin@reductions$pca@cell.embeddings[,2]

plot_PCA <- function(df = BP_set$samplesheet){
  
  # Création des plots
  count_PCA <- ggplot(df, aes(count_PC1, count_PC2, col = celltype)) +
    geom_point(size = 0.7, key_glyph = draw_key_pointrange) + 
    theme_minimal() +
    theme(legend.key.size = unit(0.3, 'cm'),
          legend.text = element_text(size=10)) 
  
  bin_PCA <- ggplot(df, aes(bin_PC1, bin_PC2, col = celltype)) +
    geom_point(size = 0.7, show.legend = FALSE) + 
    theme_minimal()
  
  # Légende commune
  PCA_plot <- ggarrange(count_PCA, bin_PCA, ncol=2, nrow=1, common.legend = TRUE, legend="bottom")
  print(PCA_plot)
  
  ggsave("./Output/Plots/02-PCA.png", PCA_plot, bg = "white", 
         height = 4, width = 6)
}

plot_PCA()

## VizDimLoadings

plot_VizDimLoadings <- function() {
  
  plot1 <- VizDimLoadings(BP_count, dims = 1:2, reduction = "pca") + 
    theme(plot.title = element_text(hjust = 1.5)) +
    ggtitle("Visualisation Dim Loadings - count\n")
  
  print(plot1)
  ggsave("./Output/Plots/02-VizDimLoadings_count.png", plot1,
         height = 6, width = 7)
  
  plot2 <- VizDimLoadings(BP_bin, dims = 1:2, reduction = "pca")+ 
    theme(plot.title = element_text(hjust = 1.5)) +
    ggtitle("Visualisation Dim Loadings - bin\n")
  
  print(plot2)
  ggsave("./Output/Plots/02-VizDimLoadings_bin.png", plot2,
         height = 6, width = 7)
}

plot_VizDimLoadings()

# Fin partie 02

saveRDS(BP_count, "./Output/Data/02-BP_count.rds")
rm(BP_count)

saveRDS(BP_bin, "./Output/Data/02-BP_bin.rds")
rm(BP_bin)

saveRDS(BP_set, "./Output/Data/02-BP_set.rds")
rm(BP_set)

gc()

# 3. Réduction de dimensions non linéaire : UMAP ####

library(Seurat)
library(scMiko)
library(ggplot2)
library(ggpubr)
library(viridis)
library(scMiko)
library(reshape2)

BP_count <- readRDS("./Output/Data/02-BP_count.rds")
BP_bin <- readRDS("./Output/Data/02-BP_bin.rds")
BP_set <- readRDS("./Output/Data/02-BP_set.rds")
custom_colors <- readRDS("./Output/Data/01-custom_colors.rds")

## Dimensionalité du dataset ####

plot_Elbow <- function() {
  
  plot1 <- ElbowPlot(BP_count, ndims = 20, reduction = "pca") + 
    geom_vline(xintercept = 15, linetype = "dashed", 
               color = "red", size = 0.8) +
    ggtitle("Elbow plot - count")
  
  plot2 <- ElbowPlot(BP_bin, ndims = 20, reduction = "pca") + 
    geom_vline(xintercept = 15, linetype = "dashed", 
               color = "red", size = 0.8) +
    ggtitle("Elbow plot - bin")
  
  plot <- ggarrange(plot1, plot2)
  print(plot)
  
  ggsave("./Output/Plots/03-ElbowPlot.png", plot, bg = "white",
         height = 4, width = 6)
}

plot_Elbow() 
# 1:15

## UMAP de Seurat sur données non intégrées ####

BP_count <- RunUMAP(BP_count, dims = 1:15, metric = "manhattan")
BP_bin <- RunUMAP(BP_bin, dims = 1:15, metric = "manhattan")

## DimPlot : visualisation des labels et donors sur UMAP

plot_Label_donor <- function(seurat_obj, output_file) {
  plot1 <- DimPlot(seurat_obj, reduction = "umap", group.by = "label") + 
    theme(legend.key.size = unit(0.1, 'in'))
  plot2 <- DimPlot(seurat_obj, reduction = "umap", group.by = "orig.ident") + 
    theme(legend.key.size = unit(0.1, 'in'))
  
  print(plot1 + plot2)
  ggsave(output_file, plot1 + plot2, width = 10, height = 5)
}

plot_Label_donor(BP_count, "./Output/Plots/03-Label_donor_count.png")
plot_Label_donor(BP_bin, "./Output/Plots/03-Label_donor_bin.png")

## Pairwise distances & corrélation ####

BP_set$samplesheet$count_UMAP1 <- BP_count@reductions$umap@cell.embeddings[,1]
BP_set$samplesheet$count_UMAP2 <- BP_count@reductions$umap@cell.embeddings[,2]
BP_set$samplesheet$bin_UMAP1 <- BP_bin@reductions$umap@cell.embeddings[,1]
BP_set$samplesheet$bin_UMAP2 <- BP_bin@reductions$umap@cell.embeddings[,2]

plot_Pairwise_dist <- function(df_set, output_file, cor = "spearman") {
  
  # Echantillon tiré au hasard de 5000 cellules :
  # sel <- sample(df_set$samplesheet$colnames, 5000)
  # df_sel <- df_set$samplesheet[sel,]
  # rownames(df_sel) <- df_sel$colnames
  # df_sel_count <- df_sel[,9:10]
  # df_sel_bin <- df_sel[,11:12]
  
  # Echantillon complet :
  df_sel_count <- df_set$samplesheet[,9:10]
  df_sel_bin <- df_set$samplesheet[,11:12]
  
  # Coordonnées UMAPs
  UMAPs <- list(df_sel_count,
                df_sel_bin)
  
  all_pairwise_dists <- lapply(UMAPs, function(UMAP){
    pairwise_dist <- as.matrix(dist(UMAP))
    pairwise_dist[upper.tri(pairwise_dist)] <- NA
    diag(pairwise_dist) <- NA
    pairwise_dist <- melt(pairwise_dist)
    pairwise_dist <- pairwise_dist[!is.na(pairwise_dist$value), ]
    return(pairwise_dist$value)
  }) |> do.call(what = "cbind")
  
  colnames(all_pairwise_dists) <- c("count","bin")
  all_pairwise_dists <- as.data.frame(all_pairwise_dists)
  
  # Echantillon de 10 000 des distances calculées pour la représentation graphique
  plotData <- all_pairwise_dists[sample(1:nrow(all_pairwise_dists), 10000),]
  plotData <- as.data.frame(plotData)
  plotData$bins <- cut(plotData$count, breaks = seq(0, 24, by = 1))
  
  # Corrélation calculée sur l'intégralité des distances
  cat("Corrélation de", cor, "entre les distances des données de comptage et des données binaires :\n", 
      cor(all_pairwise_dists$count, all_pairwise_dists$bin, method = cor), "\n")
  
  plot <- ggplot(plotData, aes(bins, bin)) + geom_boxplot() + theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    labs(y = "binary-PCA based UMAP pairwise distances", 
         x = "binned count-based UMAP pairwise distances") + 
    theme(plot.title = element_text(size=22))
  
  print(plot)
  ggsave(output_file, plot,
         bg = "white", height = 4, width = 6)
  
}

plot_Pairwise_dist(BP_set, "./Output/Plots/03-Pairwise_dist.png")

# Fin partie 03

saveRDS(BP_count, "./Output/Data/03-BP_count.rds")
rm(BP_count)

saveRDS(BP_bin, "./Output/Data/03-BP_bin.rds")
rm(BP_bin)

saveRDS(BP_set, "./Output/Data/03-BP_set.rds")
rm(BP_set)

rm(custom_colors)
gc()

# 4. Clustering sur données non intégrées ####

library(Seurat)
library(ggplot2)
library(ggforce)
library(tidyverse)
library(patchwork)
library(scMiko)
library(viridis)

# Métriques clustering
library(mclust)
library(clevr)

# Source fonctions scMiko modifiées
for(i in list.files("./Script/Functions")) {
  source(sprintf("./Script/Functions/%s",i))
}
rm(i)

BP_count <- readRDS("./Output/Data/03-BP_count.rds")
BP_bin <- readRDS("./Output/Data/03-BP_bin.rds")
BP_set <- readRDS("./Output/Data/03-BP_set.rds")
custom_colors <- readRDS("./Output/Data/01-custom_colors.rds")

## Clusters de Seurat : count ####

BP_count <- FindNeighbors(BP_count, dims = 1:15, annoy.metric = "manhattan")
BP_count <- FindClusters(BP_count, algorithm = 1, 
                         resolution = c(0.01, 0.05, 0.1, 0.2)) # ajout liste possible

## Clusters de Seurat : bin ####

BP_bin <- FindNeighbors(BP_bin, dims = 1:15, annoy.metric = "manhattan")
BP_bin <- FindClusters(BP_bin, algorithm = 1, 
                       resolution = c(0.01, 0.05, 0.1, 0.2)) # ajout liste possible

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

BP_count <- choix_res(BP_count, 0.05)
BP_bin <- choix_res(BP_bin, 0.05)

## Métriques d'évaluation du clustering ####

# Comptage / binaire
table("clusters count" = BP_count@meta.data$seurat_clusters,
      "clusters bin" = BP_bin@meta.data$seurat_clusters)

adjustedRandIndex(BP_count@meta.data$seurat_clusters,
                  BP_bin@meta.data$seurat_clusters) 
# res 0.01 : 0.7912694
# res 0.05 : 0.9114621
# res 0.1 : 0.9545719
# res 0.2 : 0.7568377

eval_report_clusters(BP_count@meta.data$seurat_clusters,
                     BP_bin@meta.data$seurat_clusters)

# Comptage / réel
table("clusters count" = BP_count@meta.data$seurat_clusters,
      "labels" = BP_set$samplesheet$celltype)

adjustedRandIndex(BP_count@meta.data$seurat_clusters,
                  BP_set$samplesheet$celltype)
# res 0.01 : 0.8732166
# res 0.05 : 0.9342147
# res 0.1 : 0.943493
# res 0.2 : 0.6727237

eval_report_clusters(BP_count@meta.data$seurat_clusters,
                     BP_set$samplesheet$celltype)

# Binaire / réel
table("clusters bin" = BP_bin@meta.data$seurat_clusters,
      "labels" = BP_set$samplesheet$celltype)

adjustedRandIndex(BP_bin@meta.data$seurat_clusters,
                  BP_set$samplesheet$celltype)  
# res 0.01 : 0.7139521
# res 0.05 : 0.8727645
# res 0.1 : 0.9174103
# res 0.2 : 0.7785435

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
ggsave("Output/Plots/04-AlluvialPlot_count.png", plot1, height = 6, width = 8)

plot2 <- plot_Alluvial(BP_bin)
ggsave("Output/Plots/04-AlluvialPlot_bin.png", plot2, height = 6, width = 8)

## Visualisations sur UMAP des données non intégrées ####

## Représentation UMAP des clusters (scMiko)

plot_Label_clusters <- function(seurat_obj, output_file) {
  
  plot1 <- scMiko::cluster.UMAP(so = seurat_obj, group.by = "label", repel = TRUE) + 
    theme(legend.key.size = unit(0.1, 'in'))
  plot2 <- scMiko::cluster.UMAP(so = seurat_obj, group.by = "seurat_clusters") + 
    theme(legend.key.size = unit(0.1, 'in'))
  
  print(plot1 + plot2)
  ggsave(output_file, plot1 + plot2, width = 10, height = 5)
}

plot_Label_clusters(BP_count, "./Output/Plots/04-Label_clusters_count.png")
plot_Label_clusters(BP_bin, "./Output/Plots/04-Label_clusters_bin.png")

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
ggsave("Output/Plots/04-UMAPs_non_intég_count.png", plots_umap_count,
       height = 6, width = 10.5)

plots_umap_bin <- plot_UMAP_hillje(BP_bin)
ggsave("Output/Plots/04-UMAPs_non_intég_bin.png", plots_umap_bin,
       height = 6, width = 10.5)

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
  
  # plot3 <- cowplot::plot_grid(plt.auc.dot[["0.05"]] + theme(legend.position = "right",
  #                               axis.text.x = element_text(angle = 45, hjust = 1)),
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

run_scMiko(BP_count, "./Output/Plots/04-scmiko_count")
run_scMiko(BP_bin, "./Output/Plots/04-scmiko_bin")

# Fin partie 04

saveRDS(BP_count, "./Output/Data/04-BP_count.rds")
rm(BP_count)

saveRDS(BP_bin, "./Output/Data/04-BP_bin.rds")
rm(BP_bin)

saveRDS(BP_set, "./Output/Data/04-BP_set.rds")
rm(BP_set)

rm(custom_colors)
rm(plot1)
rm(plot2)
rm(plots_umap_bin)
rm(plots_umap_count)
gc()

# 5. Intégration : Réduction de dimensions non linéaire UMAP ####

library(harmony)
library(Seurat)
library(scMiko)
library(ggplot2)
library(ggpubr)
library(viridis)
library(scMiko)
library(reshape2)

BP_count <- readRDS("./Output/Data/04-BP_count.rds")
BP_bin <- readRDS("./Output/Data/04-BP_bin.rds")
BP_set <- readRDS("./Output/Data/04-BP_set.rds")
custom_colors <- readRDS("./Output/Data/01-custom_colors.rds")

## Intégration par Harmony selon le donneur (donor ou orig.ident) ####

BP_count <- RunHarmony(BP_count, "group", plot_convergence = "TRUE") # convergé
BP_bin <- RunHarmony(BP_bin, "group", plot_convergence = "TRUE") # convergé

## UMAP de Seurat sur données intégrées ####

BP_count <- RunUMAP(BP_count, reduction = "harmony", dims = 1:15, metric = "manhattan")
BP_bin <- RunUMAP(BP_bin, reduction = "harmony", dims = 1:15, metric = "manhattan")

## DimPlot : visualisation des labels et donors sur UMAP

plot_Label_donor <- function(seurat_obj, output_file) {
  plot1 <- DimPlot(seurat_obj, reduction = "umap", group.by = "label") + 
    theme(legend.key.size = unit(0.1, 'in'))
  plot2 <- DimPlot(seurat_obj, reduction = "umap", group.by = "group") + 
    theme(legend.key.size = unit(0.1, 'in'))
  
  print(plot1 + plot2)
  ggsave(output_file, plot1 + plot2, width = 10, height = 5)
}

plot_Label_donor(BP_count, "./Output/Plots/05-Label_donor_count.png")
plot_Label_donor(BP_bin, "./Output/Plots/05-Label_donor_bin.png")

## Pairwise distances & corrélation ####

BP_set$samplesheet$count_UMAP1 <- BP_count@reductions$umap@cell.embeddings[,1]
BP_set$samplesheet$count_UMAP2 <- BP_count@reductions$umap@cell.embeddings[,2]
BP_set$samplesheet$bin_UMAP1 <- BP_bin@reductions$umap@cell.embeddings[,1]
BP_set$samplesheet$bin_UMAP2 <- BP_bin@reductions$umap@cell.embeddings[,2]

plot_Pairwise_dist <- function(df_set, output_file, cor = "spearman") {
  
  # Echantillon tiré au hasard de 5000 cellules :
  # sel <- sample(df_set$samplesheet$colnames, 5000)
  # df_sel <- df_set$samplesheet[sel,]
  # rownames(df_sel) <- df_sel$colnames
  # df_sel_count <- df_sel[,9:10]
  # df_sel_bin <- df_sel[,11:12]
  
  # Echantillon complet :
  df_sel_count <- df_set$samplesheet[,9:10]
  df_sel_bin <- df_set$samplesheet[,11:12]
  
  # Coordonnées UMAPs
  UMAPs <- list(df_sel_count,
                df_sel_bin)
  
  all_pairwise_dists <- lapply(UMAPs, function(UMAP){
    pairwise_dist <- as.matrix(dist(UMAP))
    pairwise_dist[upper.tri(pairwise_dist)] <- NA
    diag(pairwise_dist) <- NA
    pairwise_dist <- melt(pairwise_dist)
    pairwise_dist <- pairwise_dist[!is.na(pairwise_dist$value), ]
    return(pairwise_dist$value)
  }) |> do.call(what = "cbind")
  
  colnames(all_pairwise_dists) <- c("count","bin")
  all_pairwise_dists <- as.data.frame(all_pairwise_dists)
  
  # Echantillon de 10 000 des distances calculées pour la représentation graphique
  plotData <- all_pairwise_dists[sample(1:nrow(all_pairwise_dists), 10000),]
  plotData <- as.data.frame(plotData)
  plotData$bins <- cut(plotData$count, breaks = seq(0, 24, by = 1))
  
  # Corrélation calculée sur l'intégralité des distances
  cat("Corrélation de", cor, "entre les distances des données de comptage et des données binaires :\n", 
      cor(all_pairwise_dists$count, all_pairwise_dists$bin, method = cor), "\n")
  
  plot <- ggplot(plotData, aes(bins, bin)) + geom_boxplot() + theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    labs(y = "binary-PCA based UMAP pairwise distances", 
         x = "binned count-based UMAP pairwise distances") + 
    theme(plot.title = element_text(size=22))
  
  print(plot)
  ggsave(output_file, plot,
         bg = "white", height = 4, width = 6)
  
}

plot_Pairwise_dist(BP_set, "./Output/Plots/05-Pairwise_dist.png")

# Fin partie 05

saveRDS(BP_count, "./Output/Data/05-BP_count.rds")
rm(BP_count)

saveRDS(BP_bin, "./Output/Data/05-BP_bin.rds")
rm(BP_bin)

saveRDS(BP_set, "./Output/Data/05-BP_set.rds")
rm(BP_set)

rm(custom_colors)
gc()

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

# 7. Expression différentielle ####

library(Seurat)
library(presto)
library(DESeq2)
library(scCustomize)
library(ggplot2)
library(ggplotify)
library(grid)
library(ComplexHeatmap)
library(bseqsc)

BP_count <- readRDS("./Output/Data/06-BP_count.rds")
BP_bin <- readRDS("./Output/Data/06-BP_bin.rds")
BP_set <- readRDS("./Output/Data/06-BP_set.rds")
custom_colors <- readRDS("./Output/Data/01-custom_colors.rds")

## Analyse des gènes différentiellement exprimés - Seurat

perform_de_analysis <- function(seurat_obj, output_file_data, output_file_plot, resolution, 
                                clustered_dotplot = TRUE, filter_pct_diff = TRUE) {
  
  Idents(seurat_obj) <- paste0("RNA_snn_res.", resolution)
  seurat_obj@meta.data$seurat_clusters <- seurat_obj@meta.data[[paste0("RNA_snn_res.", resolution)]]
  
  DEG <- FindAllMarkers(seurat_obj,
                        slot = "counts",
                        only.pos = TRUE,
                        test.use = "wilcox",
                        verbose = TRUE)
  
  if (filter_pct_diff) {
    DEG <- Add_Pct_Diff(DEG)
    DEG <- DEG[DEG$pct_diff > 0.6, ] 
    }
  
  # Clustered dotplot
  if (clustered_dotplot) {
    top_markers <- Extract_Top_Markers(marker_dataframe = DEG, num_genes = 7, 
                                       named_vector = FALSE, make_unique = TRUE,
                                       rank_by = "pct_diff")
    plot <- Clustered_DotPlot(seurat_object = seurat_obj, plot_km_elbow = TRUE,
                              features = top_markers, k = 7)
    
    SSE <- plot[[1]]
    ggsave(paste0(output_file_plot, "_SSE.png"), plot = SSE,
           height = 4, width = 6)
    
    clus_dplot <- plot[[2]]
    plot_clus_dplot <- as.ggplot(grid.grabExpr(draw(clus_dplot)))
    ggsave(paste0(output_file_plot, "_clustered_dotplot.png"), plot_clus_dplot, 
           height = 6, width = 7)
    
  }
  
  cat("\nNombre de gènes différentiellement exprimés:\n")
  print(table(DEG$cluster, DEG$p_val_adj < 0.05))
  cat("\n")
  
  DEG_list <- split(DEG, DEG$cluster)
  
  # Save and return DE results
  saveRDS(DEG_list, output_file_data)
  return(DEG_list)
}

DEG_count <- perform_de_analysis(BP_count, output_file_data = "./Output/Data/07-DEG_count.rds", 
                                 output_file_plot = "./Output/Plots/07-count", 
                                 resolution = 0.01)
DEG_bin <- perform_de_analysis(BP_bin, output_file_data = "./Output/Data/07-DEG_bin.rds", 
                               output_file_plot = "./Output/Plots/07-bin", 
                               resolution = 0.01)

## Gènes marqueurs communs entre count et bin ####

## Liste des gènes marqueurs de référence donnés (library Bseq-sc, Baron et al 2016)
## + ceux obtenus par la littérature :

# ref : A Single-Cell Transcriptome Atlas of the Human Pancreas
# ref : Single-cell transcriptomic and spatial landscapes of the developing human pancreas
# ref : A Single-Cell Transcriptomic Map of the Human and Mouse Pancreas Reveals Inter- and Intra-cell Population Structure

list_genes_ref <- function() {
  
  ref <-  unlist(pancreasMarkers)
  cat("Marqueurs du package Bseq-sc :", length(ref), "\n")
  
  ref_biblio <- c("CPA1", "HES6", "SPINK1", # acinar
                  "IRX1", "IRX2", "GCG", "ARX", "MAFB", # alpha
                  "MAFA", "INS", "NKX6.1", "DLK1", "IAPP", "SIX3", "OLIG1", # beta
                  "LEPR", "SST", "HHEX", "POU3F1", "SIX2", "ESR1", "RXRG",  # delta
                  "PPY", "FEV", # gamma
                  "GHRL", "POU6F2",# epsilon
                  "KRT19", "CFTR", # ductal
                  "VWF", "ADGRL4", # endothelial
                  "SOX10", "CRYAB", "S100B", "NGFR", "PLP1", "PMP22", # schwann
                  "RGS5", "ADIRF", "FABP4", # quiesc
                  "PDGFRA", # activ
                  "PDGFRB", "FGF", "WNT", "TGFB", "PDFGF", # stellate
                  "SDS", # macrophage
                  "TPSAB1", # mast
                  "TRAC", "RAC2", # t cells
                  "PDX1", "ETV1", "MEIS2", # tous sauf alpha
                  "PAX6", # tous mais + faible dans delta & epsilon
                  "NEUROD1", "INSM1", "NKX2.2" # all types
  )
  
  cat("Marqueurs obtenus par bibliographie :", length(ref_biblio), "\n")
  
  ref_tot <- c(ref, ref_biblio)
  ref_tot <- unique(ref_tot)
  
  cat("Marqueurs uniques totaux :", length(ref), "\n")

  return(list(ref_lib = ref, ref_biblio = ref_biblio, ref_tot = ref_tot))
}  # Marqueurs humains

ref <- list_genes_ref()

attr(ref$ref_lib, "names") <- NULL
attr(ref$ref_biblio, "names") <- NULL
attr(ref$ref_tot, "names") <- NULL

## Récupération des gènes les plus différentiellement exprimés

top_n_markers <- function(DEG, n = 7) {
  
  top_n <- lapply(DEG, function(df) {
    if (nrow(df) > n) {
      df <- df[1:n, ] # Que les n premiers genes DE
    }
    return(df) 
  }) 
  
  top_n <- unlist(sapply(top_n, function(x) x$gene)) 
  return(top_n)
}

top_count <- top_n_markers(DEG_count, n = 1)
top_count <- unique(top_count)

top_bin <- top_n_markers(DEG_bin, n = 1)
top_bin <- unique(top_bin)

top_selected_markers <- function(DEG, n = 7) {
  
  all_markers <- unlist(sapply(DEG, function(x) x$gene)) 
  all_markers <- unique(all_markers)
  
  modified_DEG <- lapply(DEG, function(df) {
    if (nrow(df) > n) {
      df <- df[1:n, ] # Que les n premiers genes DE
    }
    return(df) 
  })
  
  lim_markers <- unlist(sapply(modified_DEG, function(x) x$gene)) 
  lim_markers <- unique(lim_markers)
  
  return(list(all_markers = all_markers, lim_markers = lim_markers))
}

# all_markers = DEG totaux par cluster
# lim_markers = top n DEG par cluster

selected_count <- top_selected_markers(DEG_count)
length(selected_count$all_markers) # 181
length(selected_count$lim_markers) # 54

selected_bin <- top_selected_markers(DEG_bin)
length(selected_bin$all_markers) # 136
length(selected_bin$lim_markers) # 35

## Comparaison des marqueurs obtenus et ceux de référence ####

comparer_markers_to_reference <- function(selected, ref) {
  
  cat("Nombre de gènes de référence :", length(ref), "\n")
  cat("Nombre de gènes sélectionnés :", length(selected), "\n\n")
  
  xtab_set <- function(A,B){
    both    <-  union(A,B)
    inA     <-  both %in% A
    inB     <-  both %in% B
    return(table(inA,inB))
  }
  
  xtab_set(selected, ref)
  
}

# Comparaison par rapports aux DEG de référence :
comparer_markers_to_reference(ref = ref$ref_tot, selected = selected_count$all_markers)
comparer_markers_to_reference(ref = ref$ref_tot, selected = selected_bin$all_markers)
comparer_markers_to_reference(ref = ref$ref_tot, selected = selected_count$lim_markers)
comparer_markers_to_reference(ref = ref$ref_tot, selected = selected_bin$lim_markers)

comparer_markers_to_reference(ref = selected_count$all_markers, selected = selected_bin$all_markers)
comparer_markers_to_reference(ref = selected_count$lim_markers, selected = selected_bin$lim_markers)

## Est-ce que ce sont les mêmes genes communs ? ####

get_common_DEG_name <- function(reference = ref$ref_tot,
                                markers_count = selected_count$all_markers,
                                markers_bin = selected_bin$all_markers) {
  sum = 0
  list1 = NULL
  for (i in reference) {
    if (i %in% reference & i %in% markers_count)
    {
      sum <- sum + 1
      list1 <- c(list1, i)
    }}
  cat("Somme des DEG communs à ceux de référence et ceux de comptage :", sum, "\n\n")
  cat("Liste : \n\n", list1, "\n\n")
  
  # return(list1)
  
  sum = 0
  list2 = NULL
  for (i in reference) {
    if (i %in% reference & i %in% markers_bin)
    {
      sum <- sum + 1
      list2 <- c(list2, i)
    }}
  cat("Somme des DEG communs à ceux de référence et ceux de binaire :", sum, "\n\n")
  cat("Liste : \n\n", list2, "\n\n")
  
  # return(list2)
  
  xtab_set <- function(A,B){
    both    <-  union(A,B)
    inA     <-  both %in% A
    inB     <-  both %in% B
    return(table(inA,inB))
  }
  
  cat("Table des DEG communs entre binaire et comptage : \n\n")
  print(xtab_set(list1, list2))
  cat("\n")
  
  return(list(count = list1, bin = list2))
  
}

DEG_to_plot <- get_common_DEG_name()

get_common_DEG_no_ref <- function(markers_count = selected_count$all_markers,
                                  markers_bin = selected_bin$all_markers) {
  
  cat("Longueur de markers_count:", length(markers_count), "\n")
  cat("Longueur de markers_bin:", length(markers_bin), "\n\n")
  
  sum = 0
  list1 = NULL
  for (i in markers_count) {
    if (i %in% markers_count & i %in% markers_bin) {
      sum <- sum + 1
      list1 <- c(list1, i)
    }
  }
  cat("Somme des DEG communs : combien de DEG de binaire sont aussi présents dans comptage :", sum, "\n\n")
  cat("Liste : \n\n", list1, "\n\n")
  
  sum = 0
  list2 = NULL
  for (i in markers_bin) {
    if (i %in% markers_bin & i %in% markers_count) {
      sum <- sum + 1
      list2 <- c(list2, i)
    }
  }
  cat("Somme des DEG communs ; combien de DEG de comptage sont aussi présents dans binaire :", sum, "\n\n")
  cat("Liste : \n\n", list2, "\n\n")
  
  xtab_set <- function(A,B){
    both    <-  union(A,B)
    inA     <-  both %in% A
    inB     <-  both %in% B
    return(table(inA,inB))
  }
  
  cat("Table des DEG communs entre binaire et comptage : \n\n")
  print(xtab_set(list1, list2))
  cat("\n")
  
  return(list(count = list1, bin = list2))
}

# DEG_to_plot <- get_common_DEG_no_ref()

## Visualisation des marqueurs ####

visualize_markers <- function(seurat_obj, res = NULL, features = ref, 
                              output_plot_file = NULL,
                              show_dotplots = FALSE, 
                              show_featureplots = FALSE,
                              show_vlnplots = FALSE, 
                              show_ridgeplots = FALSE,
                              show_UMAP_clusters = FALSE) {
  
  Idents(seurat_obj) <- paste0("RNA_snn_res.", res)
  seurat_obj@meta.data$seurat_clusters <- seurat_obj@meta.data[[paste0("RNA_snn_res.", res)]]
  
  # Rappel du nombre de clusters pour une résolution donnée
  cat("Nombre de clusters pour la résolution", res, ":",
      length(unique(as.factor(seurat_obj@meta.data$seurat_clusters))), "\n")
  
  # UMAP si show_UMAP_clusters = TRUE
  if (show_UMAP_clusters) {
    require(scMiko)
    g1 <- cluster.UMAP(so = seurat_obj, group.by = "seurat_clusters")
    print(g1)
    ggsave(paste0(output_plot_file, "_cluster.UMAP.png"), g1,
           bg = "white", height = 4, width = 6)
  }
  
  # DotPlots si show_dotplots = TRUE
  if (show_dotplots) {
    p1 <- DotPlot(seurat_obj, features = features) + RotatedAxis() + coord_flip()
    print(p1)
    ggsave(paste0(output_plot_file, "_dotplot.png"), p1,
           bg = "white", height = 4, width = 6)
  }
  
  # FeaturePlots si show_featureplots = TRUE
  if (show_featureplots) {
    p1 <- FeaturePlot(seurat_obj, features = features)
    print(p1)
    ggsave(paste0(output_plot_file, "_featureplot.png"), p1,
           bg = "white", height = 6, width = 10)
  }
  
  # ViolinPlot si show_vlnplots = TRUE
  if (show_vlnplots) {
    p1 <- VlnPlot(seurat_obj, features = features)
    print(p1)
    ggsave(paste0(output_plot_file, "_vlnplot.png"), p1,
           bg = "white", height = 6, width = 10)
  }
  
  # RidgePlots si show_ridgeplots = TRUE
  if (show_ridgeplots) {
    p1 <- RidgePlot(seurat_obj, features = features)
    print(p1)
    ggsave(paste0(output_plot_file, "_ridgeplot.png"), p1,
           bg = "white", height = 6, width = 10)
  }
}

## (Rappel) UMAP par clusters
visualize_markers(BP_count, res = 0.01, features = DEG_to_plot$count, 
                  show_UMAP_clusters = TRUE,
                  output_plot_file = "./Output/Plots/07-viz_count")
visualize_markers(BP_bin, res = 0.01, features = DEG_to_plot$bin,
                  show_UMAP_clusters = TRUE,
                  output_plot_file = "./Output/Plots/07-viz_bin")

## DotPlots
visualize_markers(BP_count, res = 0.01, show_dotplots = TRUE,
                  features = top_count, # ou DEG_to_plot$count
                  output_plot_file = "./Output/Plots/07-viz_count")
visualize_markers(BP_bin, res = 0.01, show_dotplots = TRUE,
                  features = top_bin, # ou DEG_to_plot$bin
                  output_plot_file = "./Output/Plots/07-viz_bin")

## Feature Plots
visualize_markers(BP_count, res = 0.01, show_featureplots = TRUE,
                  features = top_count, # ou head(DEG_to_plot$count, 8)
                  output_plot_file = "./Output/Plots/07-viz_count")
visualize_markers(BP_bin, res = 0.01, show_featureplots = TRUE,
                  features = top_bin, # ou head(DEG_to_plot$bin, 8)
                  output_plot_file = "./Output/Plots/07-viz_bin")

## Violin Plots
visualize_markers(BP_count, res = 0.01, show_vlnplots = TRUE,
                  features = top_count, # ou head(DEG_to_plot$count, 8)
                  output_plot_file = "./Output/Plots/07-viz_count")
visualize_markers(BP_bin, res = 0.01, features = top_bin,
                  show_vlnplots = TRUE, # ou head(DEG_to_plot$bin, 8)
                  output_plot_file = "./Output/Plots/07-viz_bin")

## Ridge Plots
visualize_markers(BP_count, res = 0.01, show_ridgeplots = TRUE,
                  features = top_count, # ou head(DEG_to_plot$count, 8)
                  output_plot_file = "./Output/Plots/07-viz_count")
visualize_markers(BP_bin, res = 0.01, show_ridgeplots = TRUE,
                  features = top_bin, # ou head(DEG_to_plot$bin, 8)
                  output_plot_file = "./Output/Plots/07-viz_bin")

## Fin partie 07

saveRDS(BP_count, "./Output/Data/07-BP_count.rds")
rm(BP_count)

saveRDS(BP_bin, "./Output/Data/07-BP_bin.rds")
rm(BP_bin)

saveRDS(BP_set, "./Output/Data/07-BP_set.rds")
rm(BP_set)

rm(custom_colors)
rm(DEG_bin)
rm(DEG_count)
rm(ref)
rm(selected_bin)
rm(selected_count)
rm(DEG_to_plot)
rm(top_bin)
rm(top_count)

gc()

# 8. Fin ####

library(rio)
library(ggplot2)
library(scales)
library(dplyr)
library(ggpubr)

data <- import_list("./Data/BP_human_mm.xlsx")

## Visualisations des critères de choix de résolution et distances #### 

## Nombre de clusters par résolution

df1 <- data[[1]]
df1$Résolution <- as.factor(df1$Résolution)
df1$`Matrice de comptage` <- as.factor(df1$`Matrice de comptage`)
df1$Intégration <- as.factor(df1$Intégration)
df1$Intégration <- factor(df1$Intégration, levels = c("Avant", "Après"))

# Plot

g1 <- ggplot(df1, aes(x = `Matrice de comptage`, y = `Nb clusters`, color = Résolution)) +
  geom_jitter(size = 3, width = 0.05, height = 0.01) +
  scale_y_continuous(name = "Nombre de clusters", breaks = c(seq(4, 17, 1))) +
  scale_x_discrete(name = "Données") +
  labs(title = "Nombre de clusters par résolution\nselon le format des données et l'intégration") +
  facet_wrap(~ Intégration) + 
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 13),
        axis.text.x = element_text(size = 12),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))
plot(g1)

p1 <- ggplot(df1, aes(x = Résolution, y = `Nb clusters`, color = `Matrice de comptage`)) +
  geom_jitter(size = 3, width = 0.05, height = 0.01) +
  scale_y_continuous(name = "Nombre de clusters", breaks = c(seq(4, 17, 1))) +
  scale_x_discrete(name = "Résolution") +
  labs(title = "Nombre de clusters par résolution\nselon le format des données et l'intégration") +
  facet_wrap(~ Intégration) + 
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 13),
        axis.text.x = element_text(size = 12),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) +
  labs(colour = "Données")
plot(p1)

## ARI par résolution

df2 <- data[[2]]
df2$Résolution <- as.factor(df2$Résolution)
df2$`Adjusted Rand Index` <- as.numeric(df2$`Adjusted Rand Index`)
df2$Comparaison <- as.factor(df2$Comparaison)
df2$Comparaison <- factor(df2$Comparaison, levels = c("Comptage-Binaire", "Binaire-Label", "Comptage-Label"))
df2 <- df2 %>%
  mutate(Comparaison = recode(Comparaison, 
                              "Comptage-Binaire" = "Comptage\nBinaire",
                              "Binaire-Label" = "Binaire\nLabel",
                              "Comptage-Label" = "Comptage\nLabel"))
df2$Intégration <- as.factor(df2$Intégration)
df2$Intégration <- factor(df2$Intégration, levels = c("Avant", "Après"))

# Plot

g2 <- ggplot(df2, aes(x = Comparaison, y = `Adjusted Rand Index`, color = Résolution)) +
  geom_jitter(size = 3, width = 0.1, height = 0) +
  scale_y_continuous(name = "ARI", breaks = breaks_pretty()) +
  scale_x_discrete(name ="Clusterings comparés") +
  labs(title = "Index de Rand Ajusté (ARI) par résolution\nselon le format des données et l'intégration") +
  facet_wrap(~ Intégration) + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 12),
        text = element_text(size = 13),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))
plot(g2)

p2 <- ggplot(df2, aes(x = Résolution, y = `Adjusted Rand Index`, color = Comparaison)) +
  geom_jitter(size = 3, width = 0.1, height = 0) +
  scale_y_continuous(name = "ARI", breaks = breaks_pretty()) +
  scale_x_discrete(name = "Résolution") +
  labs(title = "Index de Rand Ajusté (ARI) par résolution\nselon le format des données et l'intégration") +
  facet_wrap(~ Intégration) + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 12),
        text = element_text(size = 13),
        legend.key.spacing.y = unit(2, 'mm'),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) +
  labs(colour = "Clusterings\ncomparés")
plot(p2)

## Largeur moyenne de silhouette par résolution

df3 <- data[[3]]
df3$Résolution <- as.factor(df3$Résolution)
df3$`Largeur moyenne de silhouette` <- as.numeric(df3$`Largeur moyenne de silhouette`)
df3$`Matrice de comptage` <- as.factor(df3$`Matrice de comptage`)
df3$Intégration <- as.factor(df3$Intégration)
df3$Intégration <- factor(df3$Intégration, levels = c("Avant", "Après"))

# Plot

g3 <- ggplot(df3, aes(x = `Matrice de comptage` , y = `Largeur moyenne de silhouette`, color = Résolution)) +
  geom_jitter(size = 3, width = 0.1, height = 0) + 
  ylim(0, 1) +
  scale_x_discrete(name ="Données") +
  labs(title = "Largeur moyenne de silhouette par résolution\nselon le format des données et l'intégration") +
  xlab("Largeur moyenne de silhouette") +
  facet_wrap(~ Intégration) + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 12),
        text = element_text(size = 13),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))
plot(g3)

p3 <- ggplot(df3, aes(x = Résolution , y = `Largeur moyenne de silhouette`, color = `Matrice de comptage`)) +
  geom_jitter(size = 3, width = 0.1, height = 0) + 
  ylim(0, 1) +
  scale_x_discrete(name = "Résolution") +
  labs(title = "Largeur moyenne de silhouette par résolution\nselon le format des données et l'intégration") +
  xlab("Largeur moyenne de silhouette") +
  facet_wrap(~ Intégration) + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 12),
        text = element_text(size = 13),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) +
  labs(colour = "Données")
plot(p3)

## Score de spécificité / AUC (scMiko) par résolution

df4 <- data[[4]]
df4$Résolution <- as.factor(df4$Résolution)
df4$`Score de spécificité` <- as.numeric(df4$`Score de spécificité`)
df4$`Matrice de comptage` <- as.factor(df4$`Matrice de comptage`)
df4$Intégration <- as.factor(df4$Intégration)
df4$Intégration <- factor(df4$Intégration, levels = c("Avant", "Après"))

# Plot

g4 <- ggplot(df4, aes(x = `Matrice de comptage` , y = `Score de spécificité`, color = Résolution)) +
  geom_jitter(size = 3, width = 0.1, height = 0) +
  ylim(0, 1) +
  scale_x_discrete(name ="Données") +
  labs(title = "Score de spécificité par résolution\nselon le format des données et l'intégration") +
  xlab("Score de spécificité") +
  facet_wrap(~ Intégration) + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 12),
        text = element_text(size = 13),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))
plot(g4)

p4 <- ggplot(df4, aes(x = Résolution , y = `Score de spécificité`, color = `Matrice de comptage`)) +
  geom_jitter(size = 3, width = 0.1, height = 0) +
  ylim(0, 1) +
  scale_x_discrete(name = "Résolution") +
  labs(title = "Score de spécificité par résolution\nselon le format des données et l'intégration") +
  xlab("Score de spécificité") +
  facet_wrap(~ Intégration) + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 12),
        text = element_text(size = 13),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) +
  labs(colour = "Données")
plot(p4)

## GGarrange

plot1 <- ggarrange(g1, g2, g3, g4, legend = "right",
          common.legend = TRUE)
print(plot1)
ggsave("./Output/Plots/08-ggplots_critères.png", plot1,
       height = 10, width = 12, bg = "white")

plot2 <- ggarrange(p1, p2, p3, p4, 
          common.legend = FALSE)
print(plot2)
ggsave("./Output/Plots/08-ggplots_critères_2.png", plot2,
       height = 10, width = 12, bg = "white")

rm(data)
rm(df1)
rm(df2)
rm(df3)
rm(df4)
rm(g1)
rm(g2)
rm(g3)
rm(g4)
rm(p1)
rm(p2)
rm(p3)
rm(p4)
rm(plot1)
rm(plot2)
gc()

## Sauvegarder en prenant moins de place ####

BP_bin <- readRDS("./Output/Data/BP_human_0.05/07-BP_bin.rds")
BP_count <- readRDS("./Output/Data/BP_human_0.05/06-BP_count.rds")

### Comment trouver la taille de l'objet : 

print(format(object.size(BP_bin), units = "auto"))    # 1.4 Gb / 1503366400 bytes
print(format(object.size(BP_count), units = "auto"))  # 1.4 Gb / 1503342272 bytes

# Bin

slot_sizes_bin <- lapply(slotNames(BP_bin), function(x) format(object.size(slot(BP_bin, x)), units = "auto"))
# slot_sizes_bin <- lapply(slotNames(BP_bin), function(x) object.size(slot(BP_bin, x)))
names(slot_sizes_bin) <- slotNames(BP_bin)
print(slot_sizes_bin)

# slot_sizes_bin <- lapply(slotNames(BP_bin), function(x) format(object.size(slot(BP_bin, x)), units = "auto"))
slot_sizes_bin <- lapply(slotNames(BP_bin), function(x) object.size(slot(BP_bin, x)))
names(slot_sizes_bin) <- slotNames(BP_bin)
print(slot_sizes_bin)

library(ggplot2)
slot_sizes_df <- data.frame(slot = names(slot_sizes_bin), size = unlist(slot_sizes_bin))
ggplot(slot_sizes_df, aes(x = slot, y = size)) + 
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = "Size (bytes)")

format(object.size(BP_bin@assays$RNA), units = "auto") # 1.4 Gb
format(object.size(BP_bin@assays$RNA$counts), units = "auto") # 142.6 Mb
format(object.size(BP_bin@assays$RNA$data), units = "auto") # 142.6 Mb
format(object.size(BP_bin@assays$RNA$scale.data), units = "auto") # 1.1 Gb
str(BP_bin@assays$RNA$scale.data)

# Count

slot_sizes_count <- lapply(slotNames(BP_count), function(x) format(object.size(slot(BP_count, x)), units = "auto"))
# slot_sizes_count <- lapply(slotNames(BP_count), function(x) object.size(slot(BP_count, x)))
names(slot_sizes_count) <- slotNames(BP_count)
print(slot_sizes_count)

# slot_sizes_count <- lapply(slotNames(BP_count), function(x) format(object.size(slot(BP_count, x)), units = "auto"))
slot_sizes_count <- lapply(slotNames(BP_count), function(x) object.size(slot(BP_count, x)))
names(slot_sizes_count) <- slotNames(BP_count)
print(slot_sizes_count)

library(ggplot2)
slot_sizes_df <- data.frame(slot = names(slot_sizes_count), size = unlist(slot_sizes_count))
ggplot(slot_sizes_df, aes(x = slot, y = size)) + 
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = "Size (bytes)")

## - Option 1 : Bouland et al : Store as bits

# Caution: Directly manipulating Seurat object slots with bit vectors requires 
# careful handling to maintain compatibility with Seurat's functions. 
# You might need to write custom functions for accessing and working with the data

library(bit)
counts <- BP_count@assays$RNA$counts
dim(counts) # 20125 7242
binary <- (counts >= 1)               # logical matrix
List <- vector("list", ncol(binary))  # will store the bit vectors, one for each column (cell)
for(i in 1:length(List)){             # iterates through each column
  message(i)
  newBit <- bit(nrow(binary))         # its length is equal to the number of rows in binary (nb of genes)
  newBit[which(binary[,i])] <- T      # finds the indices (row numbers) where binary[,i] (the gene expression profile for the current cell) is TRUE
  List[[i]] <- newBit                 # if a gene is expressed in a cell, its corresponding bit in the newBit vector will be set to 1 (or TRUE)
}

bitStored <- format(object.size(List), units = "Mb")  # calculates the size of the List object (which now holds all the bit vectors) in megabytes
# 27.9 Mb contre 142.6 Mb

# Sauvegarde un vecteur de bits pour les données binaires
# Ne sauvegarde pas l'objet Seurat :
# - soit pas utilisable autre part qu'au début
# - soit il faudrait recrééer l'objet Seurat + refaire tourner les fonctions 
# des parties précédentes à chaque fois

# - Option 2 : Seurat::DietSeurat

?Seurat::DietSeurat

test_bin <- DietSeurat(BP_bin, dimreducs = c("pca", "umap", "harmony"),
                       layers = c("counts", "data")) # ne retire pas scale.data
View(test_bin)

# Sauvegarde l'objet Seurat
# Ne fait gagner de la place que parce que permet de retirer des slots de l'objet
# Implique, si l'on veut l'appliquer au cours du workflow, 
# de vérifier quels sont les slots absolument nécessaires à l'étape suivante

# - Option 3 : Modifier la sparse matrix

# Convert to a Sparse Logical Matrix (lgCMatrix)
# lgCMatrix uses even less memory than dgCMatrix when storing binary data 
# it only needs to track the positions of TRUE values : logical conversion

BP_bin@assays$RNA$counts <- as(BP_bin@assays$RNA$counts >= 1, "lgCMatrix")
format(object.size(BP_bin@assays$RNA$counts), units = "auto")
# 95.7 Mb contre 142.6 Mb

# - Option 4 : Numerical Precision

# scale.data : 
# Transformation de scale.data en dgcMatrix sparse ? 
# str(BP_bin@assays$RNA$scale.data)
# num [1:20125, 1:7242] -0.0607 -0.4353 -0.1628 0 -0.2111 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:20125] "A1BG" "A1CF" "A2M" "A2ML1" ...
# ..$ : chr [1:7242] "human1_lib1.final_cell_0003" "human1_lib1.final_cell_0006" "human1_lib1.final_cell_0012" "human1_lib1.final_cell_0013" ...

# Consider whether you need the default double-precision floating-point numbers (double) for normalized data. 
# Single-precision (float) can halve memory usage with a small loss of precision.

# Seurat's algorithms (related to dimensionality reduction and clustering), 
# are built on the assumption of double-precision data :
# - Introduce numerical instability, leading to inaccurate results.
# - Cause compatibility issues with certain functions.

# Alternative : 
BP_bin@assays$RNA$scale.data <- as.sparse(BP_bin@assays$RNA$scale.data)
str(BP_bin@assays$RNA$scale.data)
print((object.size(BP_bin)))    # 1.7 Gb / 1838788160 bytes ????????
format(object.size(BP_bin@assays$RNA$scale.data), units = "auto") # 1.4 Gb contre 1.1 Gb

# - Compression : use the compress argument

saveRDS(BP_bin, file = "./Output/Data/test1.rds", compress = "xz") # 129 Mo contre 960 Mo
saveRDS(BP_bin, file = "./Output/Data/test2.rds", compress = "gzip") # 960 Mo contre 960 Mo
saveRDS(BP_bin, file = "./Output/Data/test3.rds", compress = "bzip2") # 386 Mo contre 960 Mo

saveRDS(BP_count, file = "./Output/Data/test1_c.rds", compress = "xz") # 142 Mo contre 973 Mo
saveRDS(BP_count, file = "./Output/Data/test2_c.rds", compress = "gzip") # 973 Mo contre 973 Mo
saveRDS(BP_count, file = "./Output/Data/test3_c.rds", compress = "bzip2") # 400 Mo contre 973

# A FAIRE : 

# - Memory Profiling
# For a more detailed analysis of memory usage during your Seurat workflow, 
# consider using R's memory profiling tools (e.g., Rprof(), memory.profile()).
# https://bookdown.org/rdpeng/rprogdatascience/profiling-r-code.html 


