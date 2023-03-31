### =========================================================================
### HiTIMED metadata
### -------------------------------------------------------------------------
###

meta <- data.frame(
  Title = c("Example_Beta.rda","HiTIMED_Library.rda"),
  Description = c(paste0("The Example_Beta.rda ",
                         "contains a matrix of beta values from three colon ",
                         "adenocarcinoma samples 6285625076_R06C02, ",
                         "6929718079_R01C01, and 6042316015_R02C02. ",
                         "These samples are from publicly available, ",
                         "data source TCGA at https://portal.gdc.cancer.gov/. "),
                  paste0("The HiTIMED_Library.rda ",
                         "contains matrices of the the average ",
                         "DNA methylation values of the probes included ",
                         "in 6 layers of the HiTIMED ",
                         "deconvolution for 20 types of carcinomas. "
                         )),
  BiocVersion = c("3.17"),
  Genome = rep("hg19", 1),
  SourceType = rep("tar.gz", 1),
  SourceUrl = c("https://bit.ly/42HtPK9",paste0("https://bit.ly/42HtPK9, ", "https://bit.ly/3TZpTkl, ",
                     "https://bit.ly/3lBBHMR, ", "https://bit.ly/3z9P8GI, ",
                     "https://bit.ly/40hgaYO"
                     )),
  SourceVersion = "Mar 26 2023",
  Species = "Homo sapiens",
  TaxonomyId = 9606,
  Coordinate_1_based = TRUE,
  DataProvider = paste0("ArrayExpress, ", "TCGA, ","GEO"),
  Maintainer = "Ze Zhang <ze.zhang.gr@dartmouth.edu>",
  RDataClass = c("Beta Matrix") ,
  DispatchClass = c(rep("Rda",1)),
  RDataPath = c(paste0("HiTIMED/",
                       "Example_Beta.rda"),paste0("HiTIMED/",
                                                  "HiTIMED_Library.rda")),
  Tags = "",
  Notes = paste0("Zhang Z et al 2022, ","Farkas SA et al 2013, ", "Zhang W et al 2020, ",
                 "Timp W et al 2014, ", "Lennard K et al 2016")
)

write.csv(meta, file="inst/extdata/metadata.csv", row.names=FALSE)
