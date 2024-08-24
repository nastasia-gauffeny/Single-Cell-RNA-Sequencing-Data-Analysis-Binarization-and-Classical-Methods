
## Nastasia Gauffeny - nastasia.gauffeny@live.com 

### Pipeline scRNA-seq : binarisation VS méthodes classiques

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

#
