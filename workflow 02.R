
## Nastasia Gauffeny - nastasia.gauffeny@live.com 

### Pipeline scRNA-seq : binarisation VS méthodes classiques

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

#