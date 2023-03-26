### =========================================================================
### HiTIMED metadata
### -------------------------------------------------------------------------
###

meta <- data.frame(
  Title = c("HiTIMED"),
  Description = c(paste0("The HiTIMED package ",
                         "contains reference libraries derived from Illumina ",
                         "HumanMethylation450K and ",
                         "DNA methylation microarrays (Zhang Z, Salas LA ",
                         "et al. 2022) from public resources.")),
  BiocVersion = c("3.9"),
  Genome = rep("hg19", 1),
  SourceType = rep("tar.gz", 1),
  SourceUrl = paste0("https://bit.ly/42HtPK9, ", "https://bit.ly/3TZpTkl, ",
                     "https://bit.ly/3lBBHMR, ", "https://bit.ly/3z9P8GI, ",
                     "https://bit.ly/40hgaYO"
                     ),
  SourceVersion = "Mar 26 2023",
  Species = "Homo sapiens",
  TaxonomyId = 9606,
  Coordinate_1_based = TRUE,
  DataProvider = paste0("ArrayExpress, ", "TCGA, ","GEO"),
  Maintainer = "Ze Zhang <ze.zhang.gr@dartmouth.edu>",
  RDataClass = c("Beta Matrix") ,
  DispatchClass = c(rep("Rda",1)),
  RDataPath = c(paste0("HiTIMED/",
                       "data")),
  Tags = "",
  Notes = paste0("Zhang Z et al 2022, ","Farkas SA et al 2013, ", "Zhang W et al 2020, ",
                 "Timp W et al 2014, ", "Lennard K et al 2016")
)

write.csv(meta, file="inst/extdata/metadata.csv", row.names=FALSE)
