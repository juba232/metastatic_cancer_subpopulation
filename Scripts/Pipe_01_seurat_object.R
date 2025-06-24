library(Seurat)
library(dplyr)
library(readr)
library(tibble)
################################################################################

# Load data

counts <- read.csv("counts_rsem.csv.gz", row.names = 1)
meta <- read.csv("meta_data.csv", row.names = 1)

# Check dimensions match
stopifnot(all(colnames(counts) == rownames(meta)))

# Create Seurat object
seurat_obj <- CreateSeuratObject(counts = counts, meta.data = meta)
# Save object
saveRDS(seurat_obj, "seurat_raw.rds")
