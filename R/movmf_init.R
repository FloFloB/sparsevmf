#' Provides a function to initialize from prototypes or from clusters a penalized mixtures of von Mises-Fisher.
#' @import wordspace
#' @import skmeans
#' @import aricode
#' @import movMF
#' @import Matrix
#' @param X Data : matrix with column as features.
#' @param proto Directional means.
#' @param cluster Initial clusters.
#' @param shared_kappa If all components share kappa or not.
#' @md
#' @return returns a list including the following attributes:
#' * mu: directional means.
#' * alpha: mixture proportions.
#' * kappa: concentration parameters.
#' @export

movMF_initialisation_internal <- function(X, proto=NULL, cluster=NULL, shared_kappa) {
  if(is.null(proto) & is.null(cluster)) {
    stop("You need to provide prototypes or clusters")
  }
  if(is.null(cluster)) {
    ## assume identical kappas and assign in a binary way all observations to
    ## their "closest" directional mean
    cluster <- apply(tcrossprod(X, proto), 1, which.max)
    K <- nrow(proto)
  } else {
    K <- max(cluster)
  }
  ## the membership matrix
  P <- matrix(0, nrow = nrow(X), ncol= K )
  P[cbind(seq_len(nrow(X)), cluster)] <- 1
  ## alpha
  alpha <- colMeans(P)
  ## mu
  R <- crossprod(P, X)
  R_norm <- wordspace::rowNorms(R)
  if(all(R_norm>.Machine$double.eps)) {
    mu <- sweep(R, 1, R_norm, "/")
    ## kappa
    R_norm_tilde <- R_norm/(alpha*nrow(X))
    if(all(R_norm_tilde<1)) {
      if(shared_kappa) {
        kappa_s <- kappa_Banerjee(sum(R_norm)/nrow(X), ncol(X))
        kappa <- rep(kappa_s, K)
      } else {
        kappa <- kappa_Banerjee(R_norm_tilde, ncol(X))
      }
      return(list(mu = mu, alpha = alpha, kappa = kappa))
    }
  }
  ## impossible configuration
  NULL
}
