# HiTIMED

HiTIMED: Hierarchical Tumor Immune Microenvironment Epigenetic Deconvolution for accurate cell type resolution in the tumor microenvironment using tumor-type-specific DNA methylation data

The HiTIMED deconvolution estimates proportions up to 17 cell types (tumor, epithelial, endothelial, stromal, basophil, eosinophil, neutrophil, dendritic cell, monocyte, B naïve, B memory, CD4T naïve, CD4T memory, CD8T naïve, CD8T memory, T regulatory, and natural killer cells) in 3 major tumor microenvironment components (tumor, immune, angiogenic). Recommended methylation data preprocess order: pOOBAH --> Background subtraction using oob (noob) --> Dye bias correction.
The manuscript is published at https://translational-medicine.biomedcentral.com/articles/10.1186/s12967-022-03736-6.
This work is recognized by the Neukom Award and we appreciate the support from the Dartmouth Neukom Institute (https://neukom.dartmouth.edu/research/neukom-research-prizes/2022-research-prize-winners).

## Installation
```
devtools::install_github("SalasLab/HiTIMED")
```

## Load libraries 
```
library(HiTIMED)
library(FlowSorted.Blood.EPIC)
library(dplyr)
```

## Example
```
Example_Beta<-query(ExperimentHub(), "HiTIMED")[["EH8092"]]
HiTIMED_result<-HiTIMED_deconvolution(Example_Beta,"COAD",6,"tumor")
head(HiTIMED_result)
```

## EPICv2 Solution
#### if you're having trouble dowloading v2 annotation files, you can find them at 
####  https://www.dropbox.com/sh/rbxjhq9zalqq58e/AAABR8kKegXKVMeNJV8a2lJRa?dl=0
```
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
library(IlluminaHumanMethylationEPICv2manifest)
```
#### dir should be the folder containing EPICv2 IDATs
```
v2_RGset = read.metharray.exp("dir",recursive = TRUE) 
annotation(v2_RGset)["array"] = "IlluminaHumanMethylationEPICv2" #Update annotation files for v2
annotation(v2_RGset)["annotation"] = "20a1.hg38"
v2_MSet <-preprocessNoob(v2_RGset)
v2_Betas<-getBeta(v2_MSet)
v2_Betas<- sesame::betasCollapseToPfx(v2_Betas)
HiTIMED_deconvolution(v2_Betas,"COAD",6,"tumor")
```

![Figure1-page-001](https://user-images.githubusercontent.com/32206453/169862267-50e498fd-da1c-4625-a424-84de59438446.jpg)
