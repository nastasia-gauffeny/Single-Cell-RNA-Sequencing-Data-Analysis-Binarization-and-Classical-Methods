
## Nastasia Gauffeny - nastasia.gauffeny@live.com 

### Pipeline scRNA-seq : binarisation VS méthodes classiques

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

#
