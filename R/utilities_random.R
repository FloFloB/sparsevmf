## generate n unit vectors on the unit hypersphere in R^d with uniform
## sampling
random_hsphere <- function(n, d) {
  X <- matrix(rnorm(n*d),ncol=d)
  wordspace::normalize.rows(X)
}

## create zeros
## rate: percentage of zeros
## we enforce unitary norm and unique vectors if specified
sparsify <- function(X, rate, unique=FALSE, prec=sqrt(.Machine$double.eps),
                     max_try=10, verbose=FALSE) {
  Xs <- prod(dim(X))
  nb_zeros <- as.integer(round(Xs*rate))
  if(nb_zeros > 0) {
    for(i in 1:max_try) {
      if(verbose) {
        cat("Sparsify try:", i, "\n")
      }
      zeros <- sample(1:Xs, nb_zeros)
      sparse_X <- X
      sparse_X[zeros] <- 0
      sparse_X_norms <- wordspace::rowNorms(sparse_X)
      if(all(sparse_X_norms >= prec)) {
        if(unique) {
          if(nrow(unique(sparse_X))<nrow(sparse_X)) {
            next
          }
        }
        return(sweep(sparse_X,1,sparse_X_norms,"/"))
      }
    }
    warning("Cannot generate a valid sparse matrix")
    NULL
  } else {
    X
  }
}

## take a random subset with maximal separation between the vectors in the
## inner product sense (kmeans++ style)
spread_subset <- function(X, k, random=FALSE) {
  dX <- tcrossprod(X,X)
  edX <- exp(-dX)
  idx <- sample(nrow(X),1)
  remaining <- (1:nrow(X))[-idx]
  for(i in 2:k) {
    weightsM <- edX[idx, remaining,drop=FALSE]
    weights <- apply(weightsM,2,min)
    if(random) {
      new_idx <- sample(remaining,1,prob=weights)
    } else {
      new_idx <- remaining[which.max(weights)]
    }
    idx <- c(idx,new_idx)
    remaining <- setdiff(remaining, new_idx)
  }
  list(idx=idx, dIdx=dIdx <- dX[idx, idx])
}


## compute kappas to enforce a minimal separation between close prototypes
density_separating_kappa <- function(dX, threshold=0.01) {
  diag(dX) <- -1
  closest <- apply(dX,2,max)
  -log(threshold)/(1-closest)
}

overlapping <- function(theta, kappa) {
  2-2*circular::pvonmises(theta, circular::circular(0), kappa=kappa,tol=1e-50)
}

separating_kappa <- function(mu, dMu, threshold=0.05) {
  diag(dMu) <- -Inf
  closest_inner <- pmax(pmin(apply(dMu, 2, max),1),-1)
  ## angles
  closest_angle <- acos(closest_inner)
  kappa <- rep(NA,length(closest_angle))
  for(k in seq_along(closest_angle)) {
    angle <- circular::circular(closest_angle[k]/2)
    kappa[k] <- uniroot(function(x) {overlapping(angle, x) - threshold},
                        c(1e-5, 700))$root
  }
  kappa
}

rescale_kappa <- function(mu, dMu, kappa) {
  diag(dMu) <- -Inf
  closest_inner <- pmax(pmin(apply(dMu, 2, max),1),-1)
  ## angles
  2*kappa/(1-closest_inner)
}

## sample a movMF
sample_movMF <- function(n, movMF_param, compute_overlap=FALSE) {
  K <- length(movMF_param$alpha)
  d <- ncol(movMF_param$mu)
  membership <- sample(1:K,n,replace=TRUE,prob=movMF_param$alpha)
  X <- matrix(NA,nrow=n,ncol=d)
  for(k in 1:K) {
    k.idx <- which(membership==k)
    if(length(k.idx)>0) {
      X[k.idx,] <- movMF::rmovMF(length(k.idx), theta=movMF_param$kappa[k]*movMF_param$mu[k,,drop=FALSE])
    }
  }
  results <- list(X=X, membership=membership)
  if(compute_overlap) {
    best_match <- movMF_hard_assign_internal(X, movMF_param$kappa, movMF_param$mu,
                                             table(membership)/n)
    results$overlap <- 1-sum(diag(table(membership,best_match)))/n
  }
  results
}
