---
title: "TNBC scRNA-seq Metastasis Analysis"
author: "Jubayer Hasan"
output: html_document
---

```
knitr::opts_chunk$set(echo = TRUE)
library(Seurat)
library(dplyr)
library(ggplot2)

counts <- read.csv("data/counts_rsem.csv.gz", row.names = 1)
metadata <- read.csv("data/meta_data.csv", row.names = 1)

```
