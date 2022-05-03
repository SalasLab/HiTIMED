# HiTIMED

HiTIMED: Hierarchical Tumor Immune Microenvironment Epigenetic Deconvolution for accurate cell type resolution in the tumor microenvironment using tumor-type-specific DNA methylation data

The HiTIMED deconvolution estimates proportions up to 17 cell types (tumor, epithelial, endothelial, stromal, basophil, eosinophil, neutrophil, dendritic cell, monocyte, B naïve, B memory, CD4T naïve, CD4T memory, CD8T naïve, CD8T memory, T regulatory, and natural killer cells) in 3 major tumor microenvironment components (tumor, immune, angiogenic).


## Installation

devtools::install_github("SalasLab/HiTIMED")


## Load library 

library(HiTIMED)

load("data/HiTIMED_Library.RDATA")


## Deconvolution function

?HiTIMED_deconvolution
