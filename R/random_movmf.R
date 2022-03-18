#' Provides a function to produce a random movMF model with controlled separation between the clusters.
#' @import wordspace
#' @import movMF
#' @import circular
#' @param n Number of samples.
#' @param d Dimension of samples.
#' @param K Number of components.
#' @param alpha Mixture proportions (by default same proportion).
#' @param kappa Value(s) of kappa.
#' @param separation Separation between components.
#' @param oversample Sample to extract well separated vectors.
#' @param sparsify_mu Parameter to sparsify the directional means.
#' @param verbose Verbose.
#' @md
#' @return returns a list including the following attributes:
#' * X: samples
#' * mu: directional means.
#' * alpha: mixture proportions.
#' * kappa: concentration parameters.
#' * membership: to which cluster the sample belongs to.
#' * separation: the separation between components.
#' * overlap: overlap between components.
#' * sparsity: sparsity of directional means.
#' @export

random_movMF <- function(n, d, K, alpha=rep(1/K,K), kappa=NULL,
                         separation=0.05, oversample=20, sparsify_mu=0,
                         verbose=FALSE) {

  # Check parameters
  if(!is.numeric(n) || n<2){
    stop("n must be numeric and at least equal to 2")
  }
  if(!is.numeric(d) || d<1){
    stop("d must be numeric and at least equal to 1")
  }
  if(K>=n){
    stop("More clusters than distinct samples (rows)")
  }

  centers <- random_hsphere(K*oversample, d)
  prototypes <- spread_subset(centers, K)
  mu <- centers[prototypes$idx,]
  if(sparsify_mu > 0) {
    mu <- sparsify(mu, rate=sparsify_mu, unique=TRUE, verbose=verbose)
  }
  if(is.null(kappa)) {
    kappa <- separating_kappa(mu, prototypes$dIdx, separation)
  } else {
    if(length(kappa)==1) {
      kappa <- rep(kappa, K)
    }
    kappa <- rescale_kappa(mu, prototypes$dIdx, kappa)
    separation <- NA
  }
  the_sample <- sample_movMF(n, list(mu=mu, alpha=alpha, kappa=kappa),
                             compute_overlap=TRUE)
  if(sparsify_mu > 0) {
    sparsity <- mean(abs(the_sample$X)<sqrt(.Machine$double.eps))
  } else {
    sparsity <- 0
  }
  list(X=the_sample$X, mu=mu, alpha=alpha, kappa=kappa, membership=the_sample$membership,
       separation=separation, overlap=the_sample$overlap, sparsity=sparsity)
}

