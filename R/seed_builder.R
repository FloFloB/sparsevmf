#' Provides a function to create seeds to have a fair comparison between algorithms.
#' @import Matrix
#' @param X Data : matrix with column as features.
#' @param nb_run Number of runs.
#' @param nb_init Number of initialisations.
#' @param ks Number of cluster to be tested.
#' @md
#' @return returns a list of seed.
#' @export

build_seed<-function(X,nb_run=10,nb_init=10,ks=2:6){
  set.seed(0)
  seeds <- vector(mode="list", length=length(ks))
  for(k in ks) {
    print(k)
    seeds[[k]] <- vector(mode="list", length=nb_run)
    for(l in 1:nb_run) {
      PIDs <- matrix(NA,nrow=nb_init, ncol=k)
      for(j in 1:nb_init) {
        for(count in 1:10) {
          PIDs[j,] <- sample.int(nrow(X), k, replace=FALSE)
          pX <- X[PIDs[j,], , drop=FALSE]
          pCluster <- apply(tcrossprod(X, pX), 1, which.max)
          if(length(unique(pCluster))<k) {
            print("Bad sample, generating a new one")
          } else {
            count <- 0
            break
          }
        }
        if(count > 0) {
          stop("Cannot generate sufficiently spread prototypes")
        }
      }
      seeds[[k]][[l]] <- PIDs
    }
  }
  return(seeds)
}
