
#-----------------------#
# Package Documentation #
#-----------------------#

#' @name RRmix
#'
#' @title RRmix: Model-Based Classification with Simultaneous Adjustment for Unwanted Variation
#'
#' @docType package
#'
#' @author Muting Wan, Martin T. Wells, James G. Booth, Stephen Salerno
#'
#' @description NEED
#'
#' @details NEED
#'
#' The number of hidden factors \code{q} is determined using the function \code{nfactors},
#' which is based on the method described in Gao et al. (2014). Principal Component Analysis
#' (PCA) is performed on the correlation matrix of the raw data, and the number of factors is
#' chosen based on how many principal components explain \code{var.exp} of the total variation,
#' or through visual examination of the scree plot (\code{plot = TRUE}). Within a reasonable range
#' of values, \code{RRmix} is robust to the value of \code{q} (Wan 2015).
#'
#' @references
#' Gao, C., Tignor, N. L., Salit, J., Strulovici-Barel, Y., Hackett, N. R., Crystal, R. G., and
#' Mezey, J. G. (2014). HEFT: eQTL analysis of many thousands of expressed genes while simultaneously
#' controlling for hidden factors. Bioinformatics, 30(3):369-376.
#'
#' M. Wan, "Model-based classification with applications to high-dimensional data in bioinformatics,"
#' Ph.D. dissertation, Cornell University, 2015.
#'
#' @useDynLib RRmix
#'
#' @importFrom Rcpp evalCpp

NULL

  ##---------------------------------------------------------------------------------------
  ## Initialize
  ##---------------------------------------------------------------------------------------

## function for obtaining shrinkage estimate of gene-specific error variance:
mle <- function(pars, f_g, m_g) {
  alpha <- pars[1]
  betaa <- pars[2]
  fres <- -sum( (f_g/2 -1)*log(m_g) + (f_g/2)*log(f_g/2) + alpha*log(betaa) + lgamma(f_g/2+alpha) - lgamma(f_g/2) - lgamma(alpha) - (f_g/2 +alpha)*log(m_g*f_g/2+betaa))
  return(fres)
}

## compute initial values for HEFTsimp ECM algorithm:
getInitialVals <- function(G.in, n.in, Xc.in, Y.in, SNP.in,  p.0,   q.in){
  X <- as.matrix(cbind(rep(1,n.in), SNP.in))

  ### do no-HiddenEffect lmfit to get SNP-gene-association coefficient for initial values:
  if (ncol(Xc.in)>0){
  noHE.fit <- lm.fit(x = data.matrix(data.frame(intercept=1, SNP=SNP.in,  Xc=Xc.in)), y = Y.in)
  } else {
  noHE.fit <- lm.fit(x = data.matrix(data.frame(intercept=1, SNP=SNP.in)), y = Y.in)
  }
  f_g <- rep(n.in-{2+ncol(Xc.in)}, G.in)
  m_g <- colSums(resid(noHE.fit)^2)/f_g
  assoc.coefs <- coef(noHE.fit)["SNP",]

  ### obtain shrinkage estimate of error variance as initial value:
  alphain <- 2+(mean(m_g))^2/var(m_g)
  betain <- mean(m_g)*(alphain-1)
  est <- nlminb(c(alphain, betain), mle, f_g=f_g, m_g=m_g, lower = c(0.5, 1e-06), upper = c(Inf, Inf))   # same lower, upper limits as LEMMA
  sig2_g.est <- (f_g/2)/(f_g/2 + est[['par']][1] + 1) * m_g + est[['par']][2]/(f_g/2 + est[['par']][1] +1)

  ### obtain initial value for 1)non-null prob and 2)difference between null and non-null mixture components:
  geneid <- match(sort(assoc.coefs),assoc.coefs) #- sort genes by coefficients from noHEfit
  psi.0 <- ifelse( (abs(mean(assoc.coefs[geneid[round(G.in-p.0*G.in+1):G.in]]) - mean(assoc.coefs[geneid[1 : round(G.in-p.0*G.in)]]))
                    > abs(mean(assoc.coefs[geneid[1: round(p.0*G.in)]]) - mean(assoc.coefs[geneid[round(p.0*G.in+1):G.in]]))),
                    mean(assoc.coefs[geneid[round(G.in-p.0*G.in+1):G.in]]) - mean(assoc.coefs[geneid[1 : round(G.in-p.0*G.in)]]),
                    mean(assoc.coefs[geneid[1: round(p.0*G.in)]]) - mean(assoc.coefs[geneid[round(p.0*G.in+1):G.in]])
				 )

  ### obtain initial value for factor loading mx via Spectral decomposition:
  if (ncol(Xc.in)>0){
  Xa <- as.matrix(cbind(X, Xc.in))
  } else {
  Xa <- as.matrix(X)
  }
  H <- Xa %*% tcrossprod(solve(crossprod(Xa)), Xa)
  compS <- function(y, n, H) {
      A <- crossprod(y, diag(1,n) - H)
	  Asqd <- crossprod(A)
	  return(Asqd)
    }
  Ylist <- as.list(data.frame(Y.in))
  S <- Reduce('+', lapply(Ylist, FUN=compS, n=n.in, H=H))/G.in
  r <- eigen(S - diag(mean(sig2_g.est), n.in))
  if (q.in == 1){
    Lam.0 <- matrix(r[['vectors']][,1] * sqrt(r[['values']][1]), nrow=n.in, ncol=1)
  } else {
    Lam.0 <- tcrossprod(r[['vectors']][,1:q.in], diag(sqrt(r[['values']][1:q.in])))
  }

  fres <- list(f_g = f_g, m_g = m_g, est=est, assoc.coefs = assoc.coefs, sig2_g = sig2_g.est,
	psi = psi.0,    Lam = Lam.0)
  return(fres)
}


## function for E-step:
estep.appvec <- function(input, ex, Xc, n, mu.curr, sig20.curr, sig21.curr, Lam.curr, p.curr,  psi.curr, q.in){
  yg <- as.vector(input[1:n])
  sig2g.curr <- input[(n+1)]
  if (length(input)>n+1){
	  betacg.curr <- input[(n+2):length(input)]
		cova1 <- tcrossprod(tcrossprod(ex, diag(c(sig20.curr,sig21.curr))), ex) + tcrossprod(Lam.curr) + diag(sig2g.curr, n)
		cova0 <- tcrossprod(tcrossprod(ex, diag(c(sig20.curr,0))), ex) + tcrossprod(Lam.curr) + diag(sig2g.curr, n)
		num <- p.curr * dmvnrm_arma(matrix(yg,1,n), mu.curr+ex[,2]*psi.curr+Xc%*%betacg.curr, cova1, FALSE )[,1]
		den <- (1-p.curr) * dmvnrm_arma(matrix(yg,1,n), mu.curr+Xc%*%betacg.curr, cova0, FALSE )[,1]
		bg.estep <- num/(num + den)
		bg.estep <- ifelse(is.na(bg.estep)==FALSE, bg.estep, 0)
	  Omega <- matrix(cbind(ex, Lam.curr), n, 2+q.in)
	  Sigma <- solve( diag(c(1/sig20.curr,1/sig21.curr, rep(1, q.in))) + 1/sig2g.curr * crossprod(Omega) )
	  Sigg11.1.estep <- Sigma[1:2, 1:2]                # 2 by 2
	  Sigg21.1.estep <- Sigma[3:{2+q.in}, 1:2]         # q by 2
	  Sigg22.1.estep <- Sigma[3:{2+q.in}, 3:{2+q.in}]  # q by q
	  Gamma <- 1/sig2g.curr * tcrossprod(Sigma, Omega) %*% (yg - mu.curr - Xc%*%betacg.curr - ex[,2]*psi.curr)  #**
	  betag.1.estep <- Gamma[1:2,]
	  Fg.1.estep <- Gamma[3:{2+q.in},]
	  Inv0 <- solve(cova0)	   # for .0 objects, cannot use Woodbury due to degeneracy in 2 by 2 matrix diag(c(sig20, 0))
	  betag.0.estep <-  diag(c(sig20.curr,0))%*% crossprod(ex, Inv0) %*% (yg - mu.curr - Xc%*%betacg.curr)
	  Fg.0.estep <- crossprod(Lam.curr, Inv0) %*% (yg-mu.curr - Xc%*%betacg.curr)
	  Sigg11.0.estep <- diag(c(sig20.curr,0)) - diag(c(sig20.curr,0)) %*% crossprod(ex, Inv0) %*% tcrossprod(ex, diag(c(sig20.curr,0)))
	  Sigg21.0.estep <- (-1) * crossprod(Lam.curr, Inv0) %*% tcrossprod(ex, diag(c(sig20.curr,0)))
	  Sigg22.0.estep <- diag(q.in) - crossprod(Lam.curr, Inv0) %*% Lam.curr
  } else {
		cova1 <- tcrossprod(tcrossprod(ex, diag(c(sig20.curr,sig21.curr))), ex) + tcrossprod(Lam.curr) + diag(sig2g.curr, n)
		cova0 <- tcrossprod(tcrossprod(ex, diag(c(sig20.curr,0))), ex) + tcrossprod(Lam.curr) + diag(sig2g.curr, n)
		num <- p.curr * dmvnrm_arma(matrix(yg,1,n), mu.curr+ex[,2]*psi.curr, cova1, FALSE )[,1]
		den <- (1-p.curr) * dmvnrm_arma(matrix(yg,1,n), mu.curr, cova0, FALSE )[,1]
		bg.estep <- num/(num + den)
		bg.estep <- ifelse(is.na(bg.estep)==FALSE, bg.estep, 0)
	  Omega <- matrix(cbind(ex, Lam.curr), n, 2+q.in)
	  Sigma <- solve( diag(c(1/sig20.curr,1/sig21.curr, rep(1, q.in))) + 1/sig2g.curr * crossprod(Omega) )
	  Sigg11.1.estep <- Sigma[1:2, 1:2]                # 2 by 2
	  Sigg21.1.estep <- Sigma[3:{2+q.in}, 1:2]         # q by 2
	  Sigg22.1.estep <- Sigma[3:{2+q.in}, 3:{2+q.in}]  # q by q
	  Gamma <- 1/sig2g.curr * tcrossprod(Sigma, Omega) %*% (yg - mu.curr - ex[,2]*psi.curr)  #**
	  betag.1.estep <- Gamma[1:2,]
	  Fg.1.estep <- Gamma[3:{2+q.in},]
	  Inv0 <- solve(cova0)	   # for .0 objects, cannot use Woodbury due to degeneracy in 2 by 2 matrix diag(c(sig20, 0))
	  betag.0.estep <-  diag(c(sig20.curr,0))%*% crossprod(ex, Inv0) %*% (yg - mu.curr)
	  Fg.0.estep <- crossprod(Lam.curr, Inv0) %*% (yg-mu.curr)
	  Sigg11.0.estep <- diag(c(sig20.curr,0)) - diag(c(sig20.curr,0)) %*% crossprod(ex, Inv0) %*% tcrossprod(ex, diag(c(sig20.curr,0)))
	  Sigg21.0.estep <- (-1) * crossprod(Lam.curr, Inv0) %*% tcrossprod(ex, diag(c(sig20.curr,0)))
	  Sigg22.0.estep <- diag(q.in) - crossprod(Lam.curr, Inv0) %*% Lam.curr
  }
  fres <- c(bg.estep, betag.1.estep, betag.0.estep, Fg.1.estep, Fg.0.estep,
		c(Sigg11.1.estep), c(Sigg11.0.estep[1]), c(Sigg21.1.estep), c(Sigg21.0.estep[,1]),
		c(Sigg22.1.estep), c(Sigg22.0.estep))
  return(fres)
}


#----------#
# runRRmix #
#----------#

#' @title Emperical Bayes ECM Implementation of RRmix Algorithm
#'
#' @description This function provides a variant on the Expectation-Maximization (EM) Algorithm for the estimation
#' of the \code{RRmix} hierarchical mixture model parameters. Empirical Bayes inference is implemented via an EM
#' algorithm with classification based on the posterior expectation of the latent indicator variables. This model-based
#' classification method with simultaneous adjustment for unwanted variation is adapted to solve association detection
#' problems for high-dimensional data in the presence of unwanted variation.
#'
#' @param Y.in An \code{n} by \code{G} matrix of \code{G} gene expression variables, or compound abundances (standardized), on \code{n} samples.
#'
#' @param SNP.in A \code{G} by \code{1} indicator vector for treatment condition or minor allele presence in a single SNP (standardized).
#'
#' @param Xc.in \code{n} by \code{r} matrix of \code{r} known covariates (standardized, default = \code{matrix(nrow=0,ncol=0)}).
#'
#' @param betac.0 \code{n} by \code{r} matrix of initial values for EM estimated covariate term coefficients (default = \code{matrix(nrow=0,ncol=0)}).
#'
#' @param sig20.0 Initial value for EM estimated first variance component for \eqn{\beta_g | b_g} (default = 1.0).
#'
#' @param sig21.0 Initial value for EM estimated second variance component for \eqn{\beta_g | b_g} (default = 0.1).
#'
#' @param p.0 Initial value for EM estimated proportion of differentially abundant compounds (default = 0.05).
#'
#' @param er_tol.in Convergence criterion for EM algorithm (default = 0.001).
#'
#' @param q.in The number of factors to be estimated (default = 1).
#'
#' @return A list object containing the following named attributes:
#' \item{er.all}{Vector of error/tolerances at each iteration of the EM algorithm until convergence.}
#' \item{lc}{Vector of likelihood values at each iteration of the EM algorithm until convergence.}
#' \item{mu}{Vector of EM estimated sample-level means over genes}
#' \item{sig20}{EM estimated first variance component for \eqn{\beta_g | b_g}.}
#' \item{sig21}{EM estimated second variance component for \eqn{\beta_g | b_g}.}
#' \item{Lam}{Estimated \eqn{n\times q} loading matrix \eqn{(\Lambda)} for the factor analysis component of the model.}
#' \item{sig2_g}{Estimated compound-specific error variance.}
#' \item{b_g}{Vector of the posterior probabilities of differential abundance for each of the \code{G} compounds.}
#' \item{p}{EM estimated proportion of differentially abundant compounds.}
#' \item{psi}{EM estimated main effect of differential abundance.}
#' \item{beta_g.1}{Estimated \eqn{2\times G} coefficient matrix for the primary effect of treatment condition.}
#' \item{beta_g.0}{Initial \eqn{2\times G} coefficient matrix for the primary effect of treatment condition.}
#' \item{F_g.1}{Estimated \eqn{q\times G} factor matrix for the factor analysis component of the model.}
#' \item{F_g.0}{Initial \eqn{q\times G} factor matrix for the factor analysis component of the model.}
#' \item{assoc.coefs}{Vector of \code{G} associated coefficients for each of the \code{G} compounds.}
#' \item{sig20.0}{Initial value for EM estimated first variance component for \eqn{\beta_g | b_g}.}
#' \item{sig21.0}{Initial value for EM estimated second variance component for \eqn{\beta_g | b_g}.}
#' \item{p.0}{Initial value for EM estimated proportion of differentially abundant compounds.}
#' \item{psi.0}{Initial value for EM estimated main effect of differential abundance.}
#'
#' @references
#' M. Wan, "Model-based classification with applications to high-dimensional data in bioinformatics,"
#' Ph.D. dissertation, Cornell University, 2015.
#'
#' @examples
#' expr    <- log(operators[, 1:12])                      # Use log of First Two Operators' Abundances
#' expr    <- t(expr)                                     # Transpose Matrix for n x G
#' trmt    <- c(rep(0,3), rep(1,3), rep(0,3), rep(1,3))   # SNP Vector
#' results <- runRRmix(Y.in=expr, SNP.in=trmt, q.in=1)    # RRmix Results
#'
#' names(results)                                         # Display Result Components
#'
#' @export


runRRmix <- function(Y.in, SNP.in, Xc.in = matrix(nrow=0,ncol=0), betac.0 = matrix(nrow=0,ncol=0),
                     sig20.0 = 1.0, sig21.0 = 0.1, p.0 = 0.05, er_tol.in = 10^(-3), q.in = 1){

  G.in  <- ncol(Y.in)
  n.in  <- nrow(Y.in)
  mu.0 = 1/G.in * as.matrix(Y.in) %*% rep(1,G.in)

  ### obtain output from getInitialVals():

  gIVtemp <- getInitialVals(G.in, n.in, Xc.in, Y.in, SNP.in, p.0, q.in)
  f_g <- gIVtemp[['f_g']]
  m_g <- gIVtemp[['m_g']]
  est <- gIVtemp[['est']]
  assoc.coefs <- gIVtemp[['assoc.coefs']]  #**
  psi <- psi.0 <- gIVtemp[['psi']]         #**
  Lam <- Lam.0 <- gIVtemp[['Lam']]
  sig2_g <- sig2_g.est <- gIVtemp[['sig2_g']]

  mu <- mu.0
  betac <- betac.0
  sig20 <- sig20.0
  sig21 <- sig21.0
  p <- p.0

  X <- as.matrix(cbind(rep(1,n.in), SNP.in))


  lc <- -10^6
  er <- 10^6
  er.all <- er

  k <- 0

  while (is.na(er)==FALSE & er > er_tol.in & k < 10^5) {

     k <- k+1

     if (k >= 2) {
       Y.adj <- Y.in - Lam %*% F_g.1
	   if (ncol(Xc.in)>0) {
			adj.fit <- lm.fit(x = data.matrix(data.frame(intercept=1, SNP=SNP.in, Xc=Xc.in)), y = Y.adj)
	   } else {
			adj.fit <- lm.fit(x = data.matrix(data.frame(intercept=1, SNP=SNP.in)), y = Y.adj)
	   }
       m_g.adj <- colSums(resid(adj.fit)^2)/f_g
       alphain.adj <- 2+(mean(m_g.adj))^2/var(m_g.adj)
       betain.adj <- mean(m_g.adj)*(alphain.adj-1)
       est.adj <- nlminb(c( alphain.adj, betain.adj), mle, f_g=f_g, m_g=m_g.adj, lower = c(0.5, 1e-06), upper = c(Inf, Inf))
       sig2_g <- (f_g/2)/(f_g/2 + est.adj[['par']][1] + 1) * m_g.adj + est.adj[['par']][2]/(f_g/2 + est.adj[['par']][1] +1)
     }

     ###  E step:
	 if (ncol(Xc.in)>0) {
     bigmat <- as.list(data.frame(matrix(rbind(Y.in, matrix(sig2_g, 1, G.in), betac), nrow = n.in+1+ncol(Xc.in), ncol=G.in)))
	 } else {
     bigmat <- as.list(data.frame(matrix(rbind(Y.in, matrix(sig2_g, 1, G.in)), nrow = n.in+1, ncol=G.in)))
	 }
     matout <- lapply(bigmat, FUN=estep.appvec, ex=X, Xc=Xc.in, n=n.in,
						mu.curr=mu, sig20.curr=sig20, sig21.curr=sig21, Lam.curr=Lam, p.curr=p, psi.curr=psi, q.in=q.in)
     b_g <- sapply(matout, function(m) m[1], simplify=TRUE)
     beta_g.1 <- sapply(matout, function(m) m[2:3], simplify=TRUE)
     beta_g.0 <- sapply(matout, function(m) m[4:5], simplify=TRUE)
     F_g.1 <- matrix(sapply(matout, function(m) m[6:{6+q.in-1}], simplify=TRUE), nrow=q.in, ncol=G.in)
     F_g.0 <- matrix(sapply(matout, function(m) m[{6+q.in}:{6+q.in*2-1}], simplify=TRUE), nrow=q.in, ncol=G.in)
     Sig_g11.1 <- sapply(matout, function(x) matrix(x[{6+q.in*2}:{6+q.in*2+2^2-1}], 2, 2), simplify=FALSE)  ; names(Sig_g11.1) <- NULL
     Sig_g11.0 <- sapply(matout, function(x) diag(c(x[{6+q.in*2+2^2}],0)), simplify=FALSE)  ; names(Sig_g11.0) <- NULL
     Sig_g21.1 <- sapply(matout, function(x) matrix(x[{6+q.in*2+2^2+1}:{6+q.in*2+2^2+2*q.in}], q.in, 2), simplify=FALSE)  ; names(Sig_g21.1) <- NULL
     Sig_g21.0 <- sapply(matout, function(x) matrix(c(x[{6+q.in*2+2^2+2*q.in+1}:{6+q.in*2+2^2+2*q.in+q.in}],rep(0,q.in)), q.in, 2), simplify=FALSE)  ; names(Sig_g21.0) <- NULL
     Sig_g22.1 <- sapply(matout, function(x) matrix(x[{6+q.in*2+2^2+2*q.in+q.in+1}:{6+q.in*2+2^2+2*q.in+q.in+q.in^2}], q.in, q.in), simplify=FALSE)  ; names(Sig_g22.1) <- NULL
     Sig_g22.0 <- sapply(matout, function(x) matrix(x[{6+q.in*2+2^2+2*q.in+q.in+q.in^2+1}:{length(x)}], q.in, q.in), simplify=FALSE)  ; names(Sig_g22.0) <- NULL

     ###  M step:
        #------ mu:
		if (ncol(Xc.in) > 0){
        mu <- ( 1/sum(1/sig2_g) *  rowSums( sweep(Y.in - Xc.in%*%betac
				- sweep(X %*% beta_g.1 + Lam %*% F_g.1 + matrix(X[,2]*psi, n.in, G.in), MARGIN=2, b_g, '*')
				- sweep(X %*% beta_g.0 + Lam %*% F_g.0, MARGIN=2, (1-b_g), '*') , MARGIN=2, 1/as.vector(sig2_g), '*') )
			   )
		} else {
        mu <- ( 1/sum(1/sig2_g) *  rowSums( sweep(Y.in
				- sweep(X %*% beta_g.1 + Lam %*% F_g.1 + matrix(X[,2]*psi, n.in, G.in), MARGIN=2, b_g, '*')
				- sweep(X %*% beta_g.0 + Lam %*% F_g.0, MARGIN=2, (1-b_g), '*') , MARGIN=2, 1/as.vector(sig2_g), '*') )
			   )
		}
        #------ sig20, sig21:
        sig20 <- sum( b_g * (beta_g.1[1,]^2 + sapply(Sig_g11.1, function(m) m[1,1]) )
					+ (1-b_g) * (beta_g.0[1,]^2 + sapply(Sig_g11.0, function(m) m[1,1]) ) ) / G.in
        sig21 <- sum( b_g * (beta_g.1[2,]^2 + sapply(Sig_g11.1, function(m) m[2,2])) ) / sum(b_g)
        #------ Lam:
		if (ncol(Xc.in) > 0){
        Lam <- ((tcrossprod({Y.in - Xc.in%*%betac - matrix(mu+X[,2]*psi, n.in, G.in) - X %*% beta_g.1}, sweep(F_g.1, MARGIN=2, b_g/sig2_g, '*'))
                 - Reduce('+', sapply(mapply('*', Sig_g21.1, b_g/sig2_g, SIMPLIFY=FALSE), function(m) tcrossprod(X,m), simplify=FALSE))
                 + tcrossprod({Y.in - Xc.in%*%betac - matrix(mu, n.in, G.in) - X %*% beta_g.0}, sweep(F_g.0, MARGIN=2, (1-b_g)/sig2_g, '*'))
                 - Reduce('+', sapply(mapply('*', Sig_g21.0, (1-b_g)/sig2_g, SIMPLIFY=FALSE), function(m) tcrossprod(X,m), simplify=FALSE))
				 ) %*%  solve(  tcrossprod(F_g.1, sweep(F_g.1, MARGIN=2, b_g/sig2_g, '*'))
                            + Reduce('+', mapply('*', Sig_g22.1, b_g/sig2_g, SIMPLIFY=FALSE))
                            + tcrossprod(F_g.0, sweep(F_g.0, MARGIN=2, (1-b_g)/sig2_g, '*'))
                            + Reduce('+', mapply('*', Sig_g22.0, (1-b_g)/sig2_g, SIMPLIFY=FALSE)) )
			    )
		} else {
        Lam <- ((tcrossprod({Y.in - matrix(mu+X[,2]*psi, n.in, G.in) - X %*% beta_g.1}, sweep(F_g.1, MARGIN=2, b_g/sig2_g, '*'))
                 - Reduce('+', sapply(mapply('*', Sig_g21.1, b_g/sig2_g, SIMPLIFY=FALSE), function(m) tcrossprod(X,m), simplify=FALSE))
                 + tcrossprod({Y.in - matrix(mu, n.in, G.in) - X %*% beta_g.0}, sweep(F_g.0, MARGIN=2, (1-b_g)/sig2_g, '*'))
                 - Reduce('+', sapply(mapply('*', Sig_g21.0, (1-b_g)/sig2_g, SIMPLIFY=FALSE), function(m) tcrossprod(X,m), simplify=FALSE))
				 ) %*%  solve(  tcrossprod(F_g.1, sweep(F_g.1, MARGIN=2, b_g/sig2_g, '*'))
                            + Reduce('+', mapply('*', Sig_g22.1, b_g/sig2_g, SIMPLIFY=FALSE))
                            + tcrossprod(F_g.0, sweep(F_g.0, MARGIN=2, (1-b_g)/sig2_g, '*'))
                            + Reduce('+', mapply('*', Sig_g22.0, (1-b_g)/sig2_g, SIMPLIFY=FALSE)) )
			    )
		}
        #------ p:
        p <- sum(b_g) / G.in
        #------ psi:
		if (ncol(Xc.in) > 0){
        psi <- ( sum(crossprod(Y.in - Xc.in%*%betac - matrix(mu, n.in, G.in) - X %*% beta_g.1 - Lam %*% F_g.1, X[,2]) * (b_g/sig2_g))
				 / (sum(b_g/sig2_g)*as.numeric(crossprod(X[,2])))
				)
		} else {
        psi <- ( sum(crossprod(Y.in - matrix(mu, n.in, G.in) - X %*% beta_g.1 - Lam %*% F_g.1, X[,2]) * (b_g/sig2_g))
				 / (sum(b_g/sig2_g)*as.numeric(crossprod(X[,2])))
				)
		}
        #------ betac:  ncol(Xc) by G
		if (ncol(Xc.in) > 0){
        betac <- tcrossprod(solve(crossprod(Xc.in)), Xc.in) %*% ( Y.in - matrix(mu, n.in, G.in)
		           - sweep(X %*% beta_g.1 + Lam %*% F_g.1 + matrix(X[,2]*psi, n.in, G.in), MARGIN=2, b_g, '*')
				   - sweep(X %*% beta_g.0 + Lam %*% F_g.0, MARGIN=2, (1-b_g), '*')
				 )
		}

     ### Compute complete data log likelihood:
	 if (ncol(Xc.in) > 0){
     lc.temp <- (- n.in/2 * sum(log(sig2_g)) - 1/2 * sum(log(sig20) + b_g * log(sig21))
                 - 1/2 * sum(colSums((Y.in - Xc.in%*%betac - matrix(mu, n.in, G.in) - matrix(X[,2]*psi, n.in, G.in))^2) * b_g/sig2_g)
                   + sum(colSums((Y.in - Xc.in%*%betac - matrix(mu, n.in, G.in) - matrix(X[,2]*psi, n.in, G.in)) * (X %*% beta_g.1 + Lam %*% F_g.1)) * b_g/sig2_g)
                   + sum(b_g) * log(p)
                   - 1/2 * sum(colSums(t(crossprod(beta_g.1, crossprod(X))) * beta_g.1) *  b_g/sig2_g)
                   - 1/2 * sum(colSums(t(crossprod(beta_g.1, diag(c(1/sig20,1/sig21)))) *beta_g.1) *  b_g)
                   - sum(colSums(t(crossprod(beta_g.1, crossprod(X, Lam))) * F_g.1) *  b_g/sig2_g)
                   - 1/2 * sum(colSums(t(crossprod(F_g.1, crossprod(Lam))) * F_g.1) *  b_g/sig2_g)
                   - 1/2 * sum(colSums(F_g.1^2) * b_g)
                   - 1/2 * sum(sapply(Sig_g11.1, function(m) sum(diag(crossprod(X)%*%m)), simplify=TRUE) * b_g/sig2_g)
                   - 1/2 * sum(sapply(Sig_g11.1, function(m) sum(diag(crossprod(diag(c(1/sig20,1/sig21)),m))), simplify=TRUE) * b_g)
                   - sum(sapply(Sig_g21.1, function(m) sum(diag(crossprod(X, Lam)%*%m)), simplify=TRUE) * b_g/sig2_g)
                   - 1/2 * sum(sapply(Sig_g22.1, function(m) sum(diag(crossprod(Lam)%*%m)), simplify=TRUE) * b_g/sig2_g)
                   - 1/2 * sum(sapply(Sig_g22.1, function(m) sum(diag(m)), simplify=TRUE) * b_g)
                 - 1/2 * sum(colSums((Y.in - Xc.in%*%betac - matrix(mu, n.in, G.in))^2) * (1-b_g)/sig2_g)
                   + sum(colSums((Y.in - Xc.in%*%betac - matrix(mu, n.in, G.in)) * (X %*% beta_g.0 + Lam %*% F_g.0)) * (1-b_g)/sig2_g)
                   + sum(1-b_g) * log(1-p)
                   - 1/2 * sum(colSums(t(crossprod(beta_g.0, crossprod(X))) * beta_g.0) *  (1-b_g)/sig2_g)
                   - 1/2 * sum(colSums(t(crossprod(beta_g.0, diag(c(1/sig20,0)))) *beta_g.0) *  (1-b_g))
                   - sum(colSums(t(crossprod(beta_g.0, crossprod(X, Lam))) * F_g.0) *  (1-b_g)/sig2_g)
                   - 1/2 * sum(colSums(t(crossprod(F_g.0, crossprod(Lam))) * F_g.0) *  (1-b_g)/sig2_g)
                   - 1/2 * sum(colSums(F_g.0^2) * (1-b_g))
                   - 1/2 * sum(sapply(Sig_g11.0, function(m) sum(diag(crossprod(X)%*%m)), simplify=TRUE) * (1-b_g)/sig2_g)
                   - 1/2 * sum(sapply(Sig_g11.0, function(m) sum(diag(crossprod(diag(c(1/sig20,0)),m))), simplify=TRUE) * (1-b_g))
                   - sum(sapply(Sig_g21.0, function(m) sum(diag(crossprod(X, Lam)%*%m)), simplify=TRUE) * (1-b_g)/sig2_g)
                   - 1/2 * sum(sapply(Sig_g22.0, function(m) sum(diag(crossprod(Lam)%*%m)), simplify=TRUE) * (1-b_g)/sig2_g)
                   - 1/2 * sum(sapply(Sig_g22.0, function(m) sum(diag(m)), simplify=TRUE) * (1-b_g))
                )
	  } else {
      lc.temp <- (- n.in/2 * sum(log(sig2_g)) - 1/2 * sum(log(sig20) + b_g * log(sig21))
                 - 1/2 * sum(colSums((Y.in - matrix(mu, n.in, G.in) - matrix(X[,2]*psi, n.in, G.in))^2) * b_g/sig2_g)
                   + sum(colSums((Y.in - matrix(mu, n.in, G.in) - matrix(X[,2]*psi, n.in, G.in)) * (X %*% beta_g.1 + Lam %*% F_g.1)) * b_g/sig2_g)
                   + sum(b_g) * log(p)
                   - 1/2 * sum(colSums(t(crossprod(beta_g.1, crossprod(X))) * beta_g.1) *  b_g/sig2_g)
                   - 1/2 * sum(colSums(t(crossprod(beta_g.1, diag(c(1/sig20,1/sig21)))) *beta_g.1) *  b_g)
                   - sum(colSums(t(crossprod(beta_g.1, crossprod(X, Lam))) * F_g.1) *  b_g/sig2_g)
                   - 1/2 * sum(colSums(t(crossprod(F_g.1, crossprod(Lam))) * F_g.1) *  b_g/sig2_g)
                   - 1/2 * sum(colSums(F_g.1^2) * b_g)
                   - 1/2 * sum(sapply(Sig_g11.1, function(m) sum(diag(crossprod(X)%*%m)), simplify=TRUE) * b_g/sig2_g)
                   - 1/2 * sum(sapply(Sig_g11.1, function(m) sum(diag(crossprod(diag(c(1/sig20,1/sig21)),m))), simplify=TRUE) * b_g)
                   - sum(sapply(Sig_g21.1, function(m) sum(diag(crossprod(X, Lam)%*%m)), simplify=TRUE) * b_g/sig2_g)
                   - 1/2 * sum(sapply(Sig_g22.1, function(m) sum(diag(crossprod(Lam)%*%m)), simplify=TRUE) * b_g/sig2_g)
                   - 1/2 * sum(sapply(Sig_g22.1, function(m) sum(diag(m)), simplify=TRUE) * b_g)
                 - 1/2 * sum(colSums((Y.in - matrix(mu, n.in, G.in))^2) * (1-b_g)/sig2_g)
                   + sum(colSums((Y.in - matrix(mu, n.in, G.in)) * (X %*% beta_g.0 + Lam %*% F_g.0)) * (1-b_g)/sig2_g)
                   + sum(1-b_g) * log(1-p)
                   - 1/2 * sum(colSums(t(crossprod(beta_g.0, crossprod(X))) * beta_g.0) *  (1-b_g)/sig2_g)
                   - 1/2 * sum(colSums(t(crossprod(beta_g.0, diag(c(1/sig20,0)))) *beta_g.0) *  (1-b_g))
                   - sum(colSums(t(crossprod(beta_g.0, crossprod(X, Lam))) * F_g.0) *  (1-b_g)/sig2_g)
                   - 1/2 * sum(colSums(t(crossprod(F_g.0, crossprod(Lam))) * F_g.0) *  (1-b_g)/sig2_g)
                   - 1/2 * sum(colSums(F_g.0^2) * (1-b_g))
                   - 1/2 * sum(sapply(Sig_g11.0, function(m) sum(diag(crossprod(X)%*%m)), simplify=TRUE) * (1-b_g)/sig2_g)
                   - 1/2 * sum(sapply(Sig_g11.0, function(m) sum(diag(crossprod(diag(c(1/sig20,0)),m))), simplify=TRUE) * (1-b_g))
                   - sum(sapply(Sig_g21.0, function(m) sum(diag(crossprod(X, Lam)%*%m)), simplify=TRUE) * (1-b_g)/sig2_g)
                   - 1/2 * sum(sapply(Sig_g22.0, function(m) sum(diag(crossprod(Lam)%*%m)), simplify=TRUE) * (1-b_g)/sig2_g)
                   - 1/2 * sum(sapply(Sig_g22.0, function(m) sum(diag(m)), simplify=TRUE) * (1-b_g))
                )
	  }

     ### Check convergence:
     er <- abs(lc.temp-lc[k])/abs(lc[k])
     er.all <- c(er.all, er)
     lc <- c(lc,  lc.temp)

  }

   HEFTmix = list(er.all=er.all, lc=lc,
                  mu = mu, sig20 = sig20, sig21 = sig21, Lam = Lam, sig2_g = sig2_g,
                  b_g = b_g, p = p, psi = psi,
                  beta_g.1 = beta_g.1, beta_g.0 = beta_g.0, F_g.1 = F_g.1, F_g.0 = F_g.0,
                  assoc.coefs = assoc.coefs,
                  sig20.0 = sig20.0, sig21.0 = sig21.0, p.0=p.0, psi.0=psi.0)

  return(HEFTmix)
}
