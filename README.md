# ExTIME

ExTIME: Extended Tumor Immune Micro-Environment cell mixture deconvolution using DNA methylation and a novel tumor-type-specific hierarchical approach. 

The ExTIME deconvolution estimates proportions up to 17 cell types (tumor, epithelial, endothelial, stromal, basophil, eosinophil, neutrophil, dendritic cell, monocyte, B naïve, B memory, CD4T naïve, CD4T memory, CD8T naïve, CD8T memory, T regulatory, and natural killer cells) in 3 major tumor microenvironment components (tumor, immune, angiogenic).


## Installation

devtools::install_github("SalasLab/ExTIME")


## Load library 

library(ExTIME)

load("data/ExTIME_Library.RDATA")


## Deconvolution function

?ExTIME_deconvolution
