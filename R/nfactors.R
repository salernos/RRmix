# nfactors: Function for Calculating q.in
# Stephen Salerno (salerno1212@gmail.com)
# Last Update: June 2, 2016


#----------#
# nfactors #
#----------#


#' @title Number of factors to be estimated in the RRmix model
#'
#' @description The number of factors to be used in the RRmix model such that
#' \code{var.exp} of the variance in the data is explained.
#'
#' @details This function performs Principal Component Analysis (PCA) on the covariance
#' matrix for the raw data (\code{expr}).  \code{var.exp} sets the criteria for the
#' minimum variance explained by the determined number of factors (\code{q.in}) for the
#' RRmix model. \code{plot = TRUE}  produces a scree plot of the variance explained by
#' each factor.
#'
#' @param expr An \code{n} by \code{G} matrix of \code{G} gene expression variables,
#' or compound abundances (standardized), on \code{n} samples.
#'
#' @param var.exp The minimum cumulative proportion of the variance required to be
#' explained by the determined number of factors (\code{q.in}).  Float on (0,1)
#' (default = 0.8)
#'
#' @param plot boolean.  Produce a scree plot of the variance explained by each
#' factor (default = \code{TRUE}).
#'
#' @return The number of factors that cumulatively explain >= \code{var.exp}
#'
#' @references
#' Gao, C., Tignor, N. L., Salit, J., Strulovici-Barel, Y., Hackett, N. R., Crystal, R. G., and
#' Mezey, J. G. (2014). HEFT: eQTL analysis of many thousands of expressed genes while simultaneously
#' controlling for hidden factors. Bioinformatics, 30(3):369-376.
#'
#' @examples
#' expr <- log(operators[, 1:12])  # Log-Transformed, First Two Operators
#' expr <- t(expr)                 # Transpose Matrix for n x G
#' q    <- nfactors(expr); q       # Number of Factors for 80% Variance
#'
#' @export



nfactors <- function(expr, var.exp = 0.8, plot = TRUE){

  G = ncol(expr)                                              # Number of Compounds

  Y = as.matrix(expr)                                         # Raw Data Matrix

  m = 1/(G-1)*tcrossprod(Y-tcrossprod(rowMeans(Y),rep(1,G)))  # Covariance Matrix

  fit <- princomp(covmat = m, cor=TRUE)                       # PCA on Covariance Matrix

  vars <- cumsum(fit$sdev^2 / sum(fit$sdev^2))                # Cumulative Variance Explained

  if (plot == T){

    plot(fit$sdev^2 / sum(fit$sdev^2), type="b", xlab="Component index",    # Scree Plot
         ylab="% variance explained", main="Scree Plot",
         ylim=c(0,1), pch=16)
    text(fit$sdev^2 / sum(fit$sdev^2),                                      # Variance Explained
         labels=round(as.vector(fit$sdev^2 / sum(fit$sdev^2)), 2), pos=3)
  }

  return(which(vars >= var.exp)[1])

}
