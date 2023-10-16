#' @title
#' HiTIMED_deconvolution
#'
#' @name
#' HiTIMED_deconvolution
#'
#' @description
#' The function estimates proportions up to 17 cell types in tumor
#' microenvironment for 20 types of carcinomas.
#'
#' @param
#' tumor_beta Methylation beta matrix or data frame from the bulk tumor samples.
#'
#' @param
#' tumor_type Specify tumor type for microenvironment deconvolution.
#' BLCA: Bladder urothelial carcinoma, BRCA: Breast invasive carcinoma,
#' CESC: Cervical squamous cell carcinoma and endocervical adenocarcinoma,
#' CHOL: Cholangiocarcinoma, COAD: Colon adenocarcinoma, ESCA:
#' Esophageal carcinoma, HNSC: Head and neck squamous cell carcinoma,
#' KIRC: Kidney clear cell renal cell carcinoma, LIHC:
#' Liver hepatocellular carcinoma, LUAD: Lung adenocarcinoma,
#' LUSC: Lung squamous cell carcinoma, OV: Ovarian addenocarcinoma,
#' PAAD: Pancreatic adenocarcinoma,
#' PRAD: Prostate adenocarcinoma, READ: Rectum adenocarcinoma,
#' STAD: Stomach adenocarcinoma, THCA: Thyroid carcinoma,
#' UCEC: Uterine corpus endometrial carcinoma
#'
#' @param
#' h Numeric variable. 
#' Specify the layer of deconvolution in the hierarchical model. Default is 6.
#'
#' @param
#' tissue_type specify whether the tissue is tumor. Default is tumor. If not
#' tumor, the function will preset the tumor purity to 0.
#'
#' @return
#' A matrix with predicted cell proportions in tumor microenvironment.
#'
#' @examples
#' #Step 1: Load example data
#' library(ExperimentHub)
#' Example_Beta<-query(ExperimentHub(), "HiTIMED")[["EH8092"]]
#' #Step 2: Run HiTIMED and show results
#' HiTIMED_result<-HiTIMED_deconvolution(Example_Beta,"COAD",6,"tumor")
#' head(HiTIMED_result)
#'
#' @import  FlowSorted.Blood.EPIC
#'
#' @import  dplyr
#'
#' @import  InfiniumPurify
#'
#' @import  tibble
#' 
#' @import  ExperimentHub
#' 
#' @importFrom minfi preprocessRaw
#'
#' @importFrom minfi getBeta
#'
#' @export

HiTIMED_deconvolution <- function(tumor_beta, tumor_type, h=6, 
                                  tissue_type="tumor"){

  if ((!is(h, "numeric"))) {
    stop(strwrap(sprintf(
      "object is of class '%s', but needs to be of
                                class 'numeric' to use this function",
      class(h)
    ),
    width = 80, prefix = " ",
    initial = ""
    ))
  }
  
  if ((!is(tumor_beta, "matrix")) && (!is(tumor_beta, "data.frame"))) {
    stop(strwrap(sprintf(
      "object is of class '%s', but needs to be of
                                class 'matrix' 'data.frame' 
      to use this function",
      class(tumor_beta)
    ),
    width = 80, prefix = " ",
    initial = ""
    ))
  }
  HiTIMED_Library<-query(ExperimentHub(), "HiTIMED")[["EH8093"]]
  for (i in 1:length(HiTIMED_Library)) {
    HiTIMED_Library[[i]]<-
      HiTIMED_Library[[i]][rownames(HiTIMED_Library[[i]])%in%rownames(tumor_beta),]
  }
  proj2<-as.data.frame(projectCellType_CP(
    tumor_beta[rownames(as.data.frame(
      HiTIMED_Library[[paste0(tumor_type,"_h2")]])),],
      as.matrix(as.data.frame(HiTIMED_Library[[paste0(tumor_type,"_h2")]])), 
      lessThanOne = TRUE))


  proj2[proj2<1e-05]<-0

  for (i in seq_len(nrow(proj2))) {
    z<-1/sum(proj2[i,])
    proj2[i,]<-z*proj2[i,]
  }

  proj3A<-as.data.frame(projectCellType_CP(
    tumor_beta[rownames(as.data.frame(
      HiTIMED_Library[[paste0(tumor_type,"_h3A")]])),],
      as.matrix(as.data.frame(HiTIMED_Library[[paste0(tumor_type,"_h3A")]])), 
      lessThanOne = TRUE))

  proj3A[proj3A<1e-05]<-0

  for (i in seq_len(nrow(proj3A))) {
    z<-1/sum(proj3A[i,])
    proj3A[i,]<-z*proj3A[i,]
  }

  proj3B<-as.data.frame(projectCellType_CP(
    tumor_beta[rownames(as.data.frame(
    HiTIMED_Library[[paste0(tumor_type,"_h3B")]])),],
    as.matrix(as.data.frame(HiTIMED_Library[[paste0(tumor_type,"_h3B")]])), 
    lessThanOne = TRUE))

  proj3B[proj3B<1e-05]<-0

  for (i in seq_len(nrow(proj3B))) {
    z<-1/sum(proj3B[i,])
    proj3B[i,]<-z*proj3B[i,]
  }



  proj4A<-as.data.frame(projectCellType_CP(
    tumor_beta[rownames(as.data.frame(
    HiTIMED_Library[[paste0(tumor_type,"_h4A")]])),],
    as.matrix(as.data.frame(HiTIMED_Library[[paste0(tumor_type,"_h4A")]])), 
    lessThanOne = TRUE))

  proj4A[proj4A<1e-05]<-0

  for (i in seq_len(nrow(proj4A))) {
    z<-1/sum(proj4A[i,])
    proj4A[i,]<-z*proj4A[i,]
  }

  proj4B<-as.data.frame(projectCellType_CP(
    tumor_beta[rownames(as.data.frame(
    HiTIMED_Library[[paste0(tumor_type,"_h4B")]])),],
    as.matrix(as.data.frame(HiTIMED_Library[[paste0(tumor_type,"_h4B")]])), 
    lessThanOne = TRUE))


  proj4B[proj4B<1e-05]<-0

  for (i in seq_len(nrow(proj4B))) {
    z<-1/sum(proj4B[i,])
    proj4B[i,]<-z*proj4B[i,]
  }



  proj5A<-as.data.frame(projectCellType_CP(
    tumor_beta[rownames(as.data.frame(
      HiTIMED_Library[[paste0(tumor_type,"_h5A")]])),],
      as.matrix(as.data.frame(HiTIMED_Library[[paste0(tumor_type,"_h5A")]])), 
      lessThanOne = TRUE))

  proj5A[proj5A<1e-05]<-0

  for (i in seq_len(nrow(proj5A))) {
    z<-1/sum(proj5A[i,])
    proj5A[i,]<-z*proj5A[i,]
  }


  proj5B<-as.data.frame(projectCellType_CP(
    tumor_beta[rownames(as.data.frame(
      HiTIMED_Library[[paste0(tumor_type,"_h5B")]])),],
      as.matrix(as.data.frame(HiTIMED_Library[[paste0(tumor_type,"_h5B")]])), 
      lessThanOne = TRUE))


  proj5B[proj5B<1e-05]<-0

  for (i in seq_len(nrow(proj5B))) {
    z<-1/sum(proj5B[i,])
    proj5B[i,]<-z*proj5B[i,]
  }



  proj5C<-as.data.frame(projectCellType_CP(
    tumor_beta[rownames(as.data.frame(
      HiTIMED_Library[[paste0(tumor_type,"_h5C")]])),],
      as.matrix(as.data.frame(HiTIMED_Library[[paste0(tumor_type,"_h5C")]])), 
      lessThanOne = TRUE))

  proj5C[proj5C<1e-05]<-0

  for (i in seq_len(nrow(proj5C))) {
    z<-1/sum(proj5C[i,])
    proj5C[i,]<-z*proj5C[i,]
  }


  proj5D<-as.data.frame(projectCellType_CP(
    tumor_beta[rownames(as.data.frame(
      HiTIMED_Library[[paste0(tumor_type,"_h5D")]])),],
      as.matrix(as.data.frame(HiTIMED_Library[[paste0(tumor_type,"_h5D")]])), 
      lessThanOne = TRUE))

  proj5D[proj5D<1e-05]<-0

  for (i in seq_len(nrow(proj5D))) {
    z<-1/sum(proj5D[i,])
    proj5D[i,]<-z*proj5D[i,]
  }


  proj6A<-as.data.frame(projectCellType_CP(
    tumor_beta[rownames(as.data.frame(
      HiTIMED_Library[[paste0(tumor_type,"_h6A")]])),],
      as.matrix(as.data.frame(HiTIMED_Library[[paste0(tumor_type,"_h6A")]])), 
      lessThanOne = TRUE))

  proj6A[proj6A<1e-05]<-0

  for (i in seq_len(nrow(proj6A))) {
    z<-1/sum(proj6A[i,])
    proj6A[i,]<-z*proj6A[i,]
  }

  proj6B<-as.data.frame(projectCellType_CP(
    tumor_beta[rownames(as.data.frame(
      HiTIMED_Library[[paste0(tumor_type,"_h6B")]])),],
      as.matrix(as.data.frame(HiTIMED_Library[[paste0(tumor_type,"_h6B")]])), 
      lessThanOne = TRUE))

  proj6B[proj6B<1e-05]<-0

  tumor.sample <- colnames(tumor_beta)
  tumor_beta_iDMC<-tumor_beta[rownames(tumor_beta)%in%rownames(as.data.frame(
    HiTIMED_Library[paste0("iDMC_",tumor_type)])),]
  idmc.dat<-as.data.frame(
    HiTIMED_Library[paste0("iDMC_",tumor_type)])[rownames(tumor_beta_iDMC),]
  purity<-c()

  for (t in tumor.sample) {
    beta.adj <- c(tumor_beta_iDMC[
      idmc.dat[,paste0("iDMC_",tumor_type,".hyper")] == TRUE, t],
                  1 - tumor_beta_iDMC[
                    idmc.dat[,paste0("iDMC_",tumor_type,".hyper")] == FALSE, t])
    pu <- InfiniumPurify:::.get_peak(beta.adj)
    purity[t] <- pu
  }

  purity_iDMC<-as.data.frame(purity)
  proj<-purity_iDMC
  proj$Other<-1-proj$purity
  colnames(proj)[1]<-"Tumor"

  if (tissue_type == "tumor"){
      h1_proj<-proj
  }else{
      proj$Tumor<-0
      proj$Other<-1
      h1_proj<-proj
  }

  identical(rownames(proj),rownames(proj2))
  proj2<-proj2[,-which(colnames(proj2)==paste0(tumor_type," Tumor")),
               drop=FALSE]
  proj2<-proj2/rowSums(proj2)
  proj<-cbind(proj[, -which(colnames(proj) %in% c("Other"))],
              proj[, c("Other")] * proj2)
  colnames(proj)[1]<-"Tumor"
  h2_proj<-proj

  identical(rownames(proj),rownames(proj3A))
  proj3A<-proj3A[,-which(colnames(proj3A) %in% 
                           c(paste0(tumor_type," Tumor"),"Immune")),
                 drop=FALSE]
  proj3A<-proj3A/rowSums(proj3A)
  proj<-cbind(proj[, -which(colnames(proj) %in% c("Angiogenic"))],
              proj[, c("Angiogenic")] * proj3A)
  h3A_proj<-proj

  identical(rownames(proj),rownames(proj3B))
  proj3B<-proj3B[,-which(colnames(proj3B) %in% 
                           c(paste0(tumor_type," Tumor"),"Angiogenic")),
                 drop=FALSE]
  proj3B<-proj3B/rowSums(proj3B)
  proj<-cbind(proj[, -which(colnames(proj) %in% c("Immune"))],
              proj[, c("Immune")] * proj3B)
  h3B_proj<-proj


  identical(rownames(proj),rownames(proj4A))
  proj4A<-proj4A[,-which(colnames(proj4A) %in% c(paste0(tumor_type," Tumor"),
                                                 "Lymphocyte","Angiogenic")),
                 drop=FALSE]
  proj4A<-proj4A/rowSums(proj4A)
  proj<-cbind(proj[, -which(colnames(proj) %in% c("Myeloid"))],
              proj[, c("Myeloid")] * proj4A)
  h4A_proj<-proj


  identical(rownames(proj),rownames(proj4B))
  proj4B<-proj4B[,-which(colnames(proj4B) %in% c(paste0(tumor_type," Tumor"),
                                                 "Myeloid","Angiogenic")),
                 drop=FALSE]
  proj4B<-proj4B/rowSums(proj4B)
  proj<-cbind(proj[, -which(colnames(proj) %in% c("Lymphocyte"))],
              proj[, c("Lymphocyte")] * proj4B)
  h4B_proj<-proj

  identical(rownames(proj),rownames(proj5A))
  proj5A<-proj5A[,-which(colnames(proj5A) %in% c(
    paste0(tumor_type," Tumor"),"Lymphocyte", "Mononuclear","Angiogenic")),
                 drop=FALSE]
  proj5A<-proj5A/rowSums(proj5A)
  proj<-cbind(proj[, -which(colnames(proj) %in% c("Granulocyte"))],
              proj[, c("Granulocyte")] * proj5A)
  h5A_proj<-proj


  identical(rownames(proj),rownames(proj5B))
  proj5B<-proj5B[,-which(colnames(proj5B) %in% 
                           c(paste0(tumor_type," Tumor"),"Lymphocyte", 
                             "Granulocyte","Angiogenic")),
                 drop=FALSE]
  proj5B<-proj5B/rowSums(proj5B)
  proj<-cbind(proj[, -which(colnames(proj) %in% c("Mononuclear"))],
              proj[, c("Mononuclear")] * proj5B)
  h5B_proj<-proj


  identical(rownames(proj),rownames(proj5C))
  proj5C<-proj5C[,which(colnames(proj5C) %in% c("Bnv","Bmem")),
                 drop=FALSE]
  proj5C<-proj5C/rowSums(proj5C)
  proj<-cbind(proj[, -which(colnames(proj) %in% c("Bcell"))],
              proj[, c("Bcell")] * proj5C)
  h5C_proj<-proj

  identical(rownames(proj),rownames(proj5D))
  proj5D<-proj5D[,which(colnames(proj5D) %in% c("CD4T","CD8T")),
                 drop=FALSE]
  proj5D<-proj5D/rowSums(proj5D)
  proj<-cbind(proj[, -which(colnames(proj) %in% c("Tcell"))],
              proj[, c("Tcell")] * proj5D)
  h5D_proj<-proj



  identical(rownames(proj),rownames(proj6A))
  proj6A<-proj6A[,which(colnames(proj6A) %in% c("CD4nv","CD4mem","Treg")),
                 drop=FALSE]
  proj6A<-proj6A/rowSums(proj6A)
  proj<-cbind(proj[, -which(colnames(proj) %in% c("CD4T"))],
              proj[, c("CD4T")] * proj6A)

  h6A_proj<-proj

  identical(rownames(proj),rownames(proj6B))
  proj6B<-proj6B[,which(colnames(proj6B) %in% c("CD8nv","CD8mem")),
                 drop=FALSE]
  proj6B<-proj6B/rowSums(proj6B)
  proj<-cbind(proj[, -which(colnames(proj) %in% c("CD8T"))],
              proj[, c("CD8T")] * proj6B)

  h6B_proj<-proj


  is.nan.data.frame <- function(x)
    do.call(cbind, lapply(x, is.nan))

  proj[is.nan.data.frame(proj)]<-0


  proj$Sum<-round(rowSums(proj),2)

  proj_low<-proj %>% filter(Sum<1)
  ID_low<-rownames(proj_low)
  # empty_list <- vector(mode = "list", length = length(ID_low))
  # names(empty_list)<-ID_low

  proj<-proj[,!colnames(proj)=="Sum"]

  proj[ID_low,]<- h6B_proj[ID_low,]

  proj$Other<-h1_proj$Other
  proj$Immune<-h2_proj$Immune
  proj$Angiogenic<-h2_proj$Angiogenic
  proj$Myeloid<-h3B_proj$Myeloid
  proj$Lymphocyte<-h3B_proj$Lymphocyte
  proj$Granulocyte<-h4B_proj$Granulocyte
  proj$Mononuclear<-h4B_proj$Mononuclear
  proj$Bcell<-h4B_proj$Bcell
  proj$Tcell<-h4B_proj$Tcell
  proj$CD4T<-h5D_proj$CD4T
  proj$CD8T<-h5D_proj$CD8T

  proja<-proj[!rownames(proj)%in%ID_low,]
  proja[is.nan.data.frame(proja)]<-0
  proj[rownames(proja),]<-proja


  if(h==1){
    output<-proj[,c("Tumor","Other")]
  }else{
    if(h==2){
      output<-proj[,c("Tumor","Immune","Angiogenic")]
    }else{
      if(h==3){
        output<-proj[,c("Tumor","Lymphocyte", "Myeloid", "Endothelial", 
                        "Epithelial", "Stromal")]
      }else{
        if(h==4){
          output<-proj[,c("Tumor","Granulocyte", "Mononuclear","Tcell","Bcell",
                          "NK","Endothelial", "Epithelial", "Stromal")]
        }else{
          if(h==5){
            output<-proj[,c("Tumor","Bas","Eos","Neu","Mono","DC","Bnv","Bmem", 
                            "CD4T", "CD8T","NK","Endothelial", "Epithelial", 
                            "Stromal")]
          }else{
            if(h==6){
              output<-proj[,c("Tumor","Endothelial", "Epithelial", "Stromal",
                              "Bnv","Bmem","CD4nv","CD4mem","Treg","CD8nv",
                              "CD8mem","Mono","DC","NK","Bas","Eos","Neu")]
            }}}}}}



  output$Sum<-rowSums(output)

  output_low<-output %>% filter(Sum == "NaN")
  ID_low<-rownames(output_low)
  if(length(ID_low)!=0){
    message(paste0("Recommend lower layer deconvolution for ",
                   toString(ID_low)))}
  return(output[,!colnames(output)=="Sum"]*100)
}
