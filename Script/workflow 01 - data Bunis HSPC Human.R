
## Nastasia Gauffeny - nastasia.gauffeny@live.com 

### Pipeline scRNA-seq : binarisation VS méthodes classiques

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
library(EnsDb.Hsapiens.v79)
library(singleCellTK)

# Quality Control (QC)
library(BPCells)
library(data.table)
library(purrr)

## Données de comptage : count ####

BP_count <- BunisHSPCData()

# Conversion des annotations des genes : Ensembl -> Gene ID
ensembl.genes <- rownames(BP_count)
geneIDs <- ensembldb::select(EnsDb.Hsapiens.v79, keys= ensembl.genes,
                             keytype = "GENEID", columns = c("SYMBOL","GENEID"))
# count : 33694 ; geneIDs : 33545
length(intersect(ensembl.genes, geneIDs[,2])) # 33545
BP_count <- BP_count[rownames(BP_count) %in% geneIDs[,2], ]
BP_count <- setRowNames(BP_count, geneIDs[,1])

BP_count <- as.Seurat(BP_count, data = NULL)
# dim 33545  5183

## Filtre & contrôle qualité (QC) ####

BP_count <- subset(BP_count,
                   subset = nFeature_originalexp > 200
                   & nFeature_originalexp < 2500) 
# dim 33545  4106

BP_count <- AddMetaData(BP_count, BP_count@meta.data$Sample, col.name = "group")

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
BP_samplesheet <- data.frame(celltype = as.factor(BP_count@meta.data$labels),
                             colnames = as.factor(BP_count@assays$originalexp@counts@Dimnames[[2]]),
                             group = as.factor(BP_count@meta.data$group),
                             labelNUM = as.numeric(as.factor(BP_count@meta.data$labels)),
                             tissue = as.factor(BP_count@meta.data$age)) # Modif : tissue -> age
BP_set$samplesheet <- BP_samplesheet

## Pieplot de la répartition des vrais labels
# source : https://stackoverflow.com/a/74931586 

plot_Label <- function(df = BP_set$samplesheet, dataset = "", max_x = 50) {
  
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
          legend.key.size = unit(0.2, "inch"),
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

plot_Label(dataset = "Bunis HSPC Human")

# Fin partie 01

saveRDS(BP_set, "./Output/Data/01-BP_set.rds")
rm(BP_set)

saveRDS(custom_colors, "./Output/Data/01-custom_colors.rds")
rm(custom_colors)

rm(BP_count)
rm(BP_bin)
rm(BP_samplesheet)
rm(geneIDs)
rm(ensembl.genes)

gc()

#
