library(Seurat)
library(dplyr)
library(readr)
library(tibble)
library(Matrix)
library(qs)
################################################################################

#loading saved object

seurat_obj <- readRDS("seurat_raw.rds")

# Violin plots of features
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA"), 
        ncol = 2, pt.size = 0.1)


FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")


#basic query

summary(seurat_obj$nFeature_RNA)

summary(seurat_obj$nCount_RNA)
#######################################

library(ggplot2)
library(patchwork)

# Extract metadata
meta <- seurat_obj@meta.data

# Plot 1: nFeature_RNA
p1 <- ggplot(meta, aes(x = nFeature_RNA)) +
  geom_density(fill = "#56B4E9", alpha = 0.6) +
  geom_vline(xintercept = c(300, 6000), linetype = "dashed", color = "red", linewidth = 0.8) +
  labs(
    title = "A. Number of Genes per Cell (nFeature_RNA)",
    x = "Genes Detected",
    y = "Density"
  ) +
  theme_bw(base_size = 14)

# Plot 2: nCount_RNA
p2 <- ggplot(meta, aes(x = nCount_RNA)) +
  geom_density(fill = "#009E73", alpha = 0.6) +
  geom_vline(xintercept = c(100000, 3000000), linetype = "dashed", color = "red", linewidth = 0.8) +
  scale_x_log10() +
  labs(
    title = "B. Total RNA Molecules per Cell (nCount_RNA)",
    x = "Total RNA Count (log10)",
    y = "Density"
  ) +
  theme_bw(base_size = 14)

combined_plot <- p1 + p2 + plot_layout(ncol = 2)
combined_plot
ggsave("combined_qc_density_plots.png", plot = combined_plot, width = 12, height = 5, dpi = 300)





