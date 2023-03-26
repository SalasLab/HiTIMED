#' HiTIMED library CpGs matrix for tumor microenvironment 
#' DNA methylation deconvolution
#'
#' @description
#'     This object contains matrices of the the average
#'     DNA methylation values of the probes included in 6 layers of the HiTIMED
#'     deconvolution for 20 types of carcinomas. These CpGs are used as the 
#'     backbone for deconvolution and were selected because their methylation 
#'     signature differs across the cell subtypes.
#'
#' @format The list contains 240 matrices
#'
#'         The format is:
#'         num [1:332, 1:3] 0.018261532  0.98263707  0.009418320 ...
#'
#' @examples
#' data("HiTIMED_Library")
#' head(HiTIMED_Library)
"HiTIMED_Library"

