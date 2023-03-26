# HiTIMED

HiTIMED: Hierarchical Tumor Immune Microenvironment Epigenetic Deconvolution for accurate cell type resolution in the tumor microenvironment using tumor-type-specific DNA methylation data

The HiTIMED deconvolution estimates proportions up to 17 cell types (tumor, epithelial, endothelial, stromal, basophil, eosinophil, neutrophil, dendritic cell, monocyte, B naïve, B memory, CD4T naïve, CD4T memory, CD8T naïve, CD8T memory, T regulatory, and natural killer cells) in 3 major tumor microenvironment components (tumor, immune, angiogenic). Recommended methylation data preprocess order: pOOBAH --> Background subtraction using oob (noob) --> Dye bias correction.

The manuscript is published at https://translational-medicine.biomedcentral.com/articles/10.1186/s12967-022-03736-6.

This work is recognized by the Neukom Award and we appreciate the support from the Dartmouth Neukom Institute (https://neukom.dartmouth.edu/research/neukom-research-prizes/2022-research-prize-winners).

## Installation
```
devtools::install_github("SalasLab/HiTIMED")
```

## Load library 
```
library(HiTIMED)

data("HiTIMED_Library")
```

## Deconvolution function
```
?HiTIMED_deconvolution
```

## Example
```
data("Example_Beta")
HiTIMED_result<-HiTIMED_deconvolution(Example_Beta,"COAD",6,"tumor")
head(HiTIMED_result)
```


![Figure1-page-001](https://user-images.githubusercontent.com/32206453/169862267-50e498fd-da1c-4625-a424-84de59438446.jpg)
