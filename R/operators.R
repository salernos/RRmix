# Operator Data Frame: Documentation Script
# Stephen Salerno (salerno1212@gmail.com)
# Last Update: May 31, 2016


#---------------------#
# Operator Data Frame #
#---------------------#


#' @docType data
#'
#' @name operators
#'
#' @title Operator Batch Effect Data Frame
#'
#' @description This dataset contains the relative abundances of 265 metabolites across a total of 24 samples.
#' The samples were collected by four different operators, namely A, X, D, and Z.
#' Each operator performed a metabolomics experiment on 6 samples: three untreated (control) replicates, and three treated replicates.
#' Both operators A and X were given samples treated with the drug 2-deoxy-glucose (2DG).
#' Metabolites are in rows (265 metabolites), and samples are in columns (24 samples).
#'
#' @usage operators
#'
#' @format A data frame with abundances of 265 filtered metabolites on 24 samples
#'
#' @source The Locasale Research Group - Department of Pharmacology & Cancer Biology, Duke University - Durham, North Carolina
#'
#' @references NEED
#'
#' @examples
#' dim(operators)      # [1] 265  24
#' summary(operators)  # Quantiles for each metabolite


NULL
