## we use the notation from the paper

## kappa computation
kappa_Banerjee <- function(Rho, d){
  Rho * (d - Rho^2)/(1 - Rho^2)
}

## k-means++ like initialisation
movMF_spread_prototype <- function(X, K) {
  dX <- tcrossprod(X,X)
  edX <- exp(-dX)
  idx <- sample(nrow(X),1)
  remaining <- (1:nrow(X))[-idx]
  for(i in 2:K) {
    weightsM <- edX[idx, remaining,drop=FALSE]
    weights <- apply(weightsM,2,min)
    new_idx <- sample(remaining,1,prob=weights)
    idx <- c(idx,new_idx)
    remaining <- setdiff(remaining, new_idx)
  }
  X[idx, , drop = FALSE]
}

## initialisation from prototypes or from clusters
movMF_initialisation_internal <- function(X, proto=NULL, cluster=NULL, shared_kappa) {
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

## initialisation
movMF_initialisation <- function(X, K, shared_kappa=FALSE, maxtry=10,
                                 mode=c("random", "spread", "skmeans"),
                                 sk_runs=1) {
  mode <- match.arg(mode)
  for(i in 1:maxtry) {
    if(mode=="random") {
      ## k random examples, without replacement
      M <- X[sample.int(nrow(X), K, replace = FALSE), , drop = FALSE]
    } else if (mode=="spread") {
      M <- movMF_spread_prototype(X, K)
    } else {
      sk <- skmeans::skmeans(X, K, method="pclust", control=list(nruns=sk_runs))
      M <- sk$prototypes
    }
    result <- movMF_initialisation_internal(X, M, shared_kappa=shared_kappa)
    if(!is.null(result)) {
      return(result)
    }
  }
  stop("Cannot find a proper initial configuration")
}

## log sum exp trick to reduce overflow
log_row_sums <-  function(L) {
  M <- L[cbind(seq_len(nrow(L)), max.col(L, "first"))]
  M + log(rowSums(exp(L - M)))
}


## soft assign
movMF_soft_assign_internal <- function(X, kappa, mu, alpha) {
  ## log likelihood matrix, ie LM[i,k]<-log alpha[k]+ log C_d(kappa[k]) +
  ## kappa[k]*mu[k,]^Tx[i,]
  ## first we compute the inner product matrix kappa[k]*mu[k,]^Tx[i,]
  inner <- tcrossprod(X, kappa*mu)
  ## log C_d(kappa[k]) using movMF internal functions
  log_cdk <- movMF:::lH(kappa, ncol(X)/2 - 1)
  ## log likelihood matrix
  LM <- sweep(inner, 2, log(alpha) - log_cdk, "+")
  ## marginal log likelihood
  LX <- log_row_sums(LM)
  ## log likelihood and assignment matrix
  list(Tau=exp(LM-LX), L=sum(LX))
}

movMF_soft_assign <- function(movMFResults, X) {
  movMF_soft_assign_internal(X, movMFResults$kappa, movMFResults$mu, movMFResults$alpha)
}

## hardassign
movMF_hard_assign_internal <- function(X, kappa, mu, alpha, Tau=NULL) {
  if(is.null(Tau)) {
    Tau <- movMF_soft_assign_internal(X, kappa, mu, alpha)$Tau
  }
  apply(Tau,1,which.max)
}

movMF_hard_assign <- function(movMFResults, X=NULL) {
  if(!is.null(X)) {
    movMF_hard_assign_internal(X, movMFResults$kappa, movMFResults$mu, movMFResults$alpha)
  } else {
    apply(movMFResults$Tau,1,which.max)
  }
}

movMF_complete_logLikelihood_internal <- function(X, kappa, mu, alpha, Tau) {
  ## hard assignments
  cluster <- apply(Tau, 1, which.max)
  log_cdk <- movMF:::lH(kappa, ncol(mu)/2 - 1)
  inner <- tcrossprod(X, kappa*mu)
  LM <- sweep(inner, 2, log(alpha) - log_cdk, "+")
  sum(LM[cbind(1:length(cluster), cluster)])
}

## EM low level
## Theta is the initial parameter list
## K is induced by this parameter list
## the status is
## 0: converged
## 1: maxiter reached
## 2: not converging to a K cluster solution
## 3: infinite variance
movMF_EM_internal <- function(X, beta, Theta, shared_kappa, maxiter, maxiterfp, prec,
                              kappamax, save_tau, hard_assign, verbose) {
  d <- ncol(X)
  alpha <- Theta$alpha
  mu <- Theta$mu
  kappa <- Theta$kappa
  K <- length(kappa)
  n <- nrow(X)
  logL <- rep(NA,maxiter)
  fpiter <- 0
  status <- 1
  for(i in 1:maxiter) {
    ## E phase
    like_and_Tau <- movMF_soft_assign_internal(X, kappa, mu, alpha)
    ## log likelihood
    L <- like_and_Tau$L
    logL[i] <- L
    if(verbose > 0) {
      cat("Iteration", i, " log likelihood:", L, "\n")
    }
    if(i>1) {
      if(abs(logL[i]-logL[i-1])<prec*(abs(logL[i-1])+prec)) {
        status <- 0
        break
      }
    }
    ## responsabilities
    Tau <- like_and_Tau$Tau

    ## M phase
    ## cluster "centers" R[k,j] <- sum_i Tau[i,k]x[i,j]
    R <- crossprod(Tau, X)
    ## alpha
    cluster_size <- colSums(Tau)
    n_alpha <- cluster_size/n
    if(beta==0) {
      ## simplified updates
      R_norm <-  wordspace::rowNorms(R)
      ## mu
      n_mu <- sweep(R, 1, R_norm, "/")
      R_norm_tilde <- R_norm/cluster_size
      if(any(R_norm_tilde >= 1)) {
        ## infinite variance
        status <- 3
        break
      }
      if(shared_kappa) {
        n_kappa_s <- kappa_Banerjee(sum(R_norm)/n, d)
        n_kappa <- rep(n_kappa_s, K)
      } else {
        n_kappa <- kappa_Banerjee(R_norm_tilde, d)
        ##n_kappa <- movMF:::solve_kappa_Newton_Fourier(R_norm/cluster_size, d)
      }
      n_kappa[n_kappa>kappamax] <- kappamax
      if( verbose > 0) {
        cat("Iteration", i, "kappa =", n_kappa, "\n")
      }
    } else {
      R_abs <- abs(R)
      R_sign <- sign(R)
      ## fixed point stategy
      ## we must start with the mu update for consistency with the
      ## beta=0 case
      n_mu <- mu
      n_kappa <- kappa
      ## we need to make sure the kappas are large enough
      kappa_min <- apply((beta+prec)/R_abs, 1, min)
      if(verbose > 0) {
        cat("Iteration", i, "Kappa bound:",kappa_min ,"\n")
      }
      n_kappa[n_kappa < kappa_min] <- kappa_min[n_kappa < kappa_min]
      n_kappa[n_kappa>kappamax] <- kappamax
      for(j in 1:maxiterfp) {
        fpiter <- fpiter+1
        ## mu
        premu <- n_kappa*R_abs - beta
        premu[premu<0] <- 0
        premu <- premu*R_sign
        premu_norm <- wordspace::rowNorms(premu)
        if(verbose > 1) {
          cat("\tInternal iteration", j, premu_norm, "\n")
          if(anyNA(premu_norm)) {
            print(premu)
            print(R_abs)
          }
        }
        ## when the regularization would set mu_k to zero, this means
        ## that the corresponding component is too wide (or too narrow) and that the
        ## EM is not converging to a solution with K clusters
        ## we stop the iteration and report a convergence error
        if(any(premu_norm<1e-16)) {
          ## we cancel the current update by not updating mu, kappa,
          ## etc.
          status <- 2
          break
        }
        n2_mu <- sweep(premu, 1, premu_norm, "/")
        if(verbose > 1) {
          cat("\tInternal iteration", j, n2_mu,"\n")
        }
        ## kappa
        pre_rho <- rowSums(n2_mu*R)
        rho <- pre_rho/cluster_size
        if(any(rho>=1)) {
          status <- 3
          break
        }
        if(shared_kappa) {
          n2_kappa_s <- kappa_Banerjee(sum(pre_rho)/n, d)
          n2_kappa <- rep(n2_kappa_s, K)
        } else {
          n2_kappa <- kappa_Banerjee(rho, d)
        }
        n2_kappa[n2_kappa < kappa_min] <- kappa_min[n2_kappa < kappa_min]
        n2_kappa[n2_kappa>kappamax] <- kappamax
        if(verbose > 1) {
          cat("\tInternal iteration", j, n2_kappa,"\n")
        }
        delta_kappa <- norm(n2_kappa - n_kappa, type = "2")
        delta_mu <- norm(n2_mu - n_mu, type = "2")
        n_kappa <- n2_kappa
        n_mu <- n2_mu
        if((delta_kappa < prec*(norm(n_kappa, type="2")+prec)) & (delta_mu < prec)) {
          break
        }
      }
      if(status==2) {
        break
      }
    }
    kappa <- n_kappa
    mu <- n_mu
    alpha <- n_alpha
  }
  logL <- logL[!is.na(logL)]
  LogLikelihood <- logL[length(logL)]
  penLogLikelihood <- LogLikelihood
  completeLogLikelihood <- movMF_complete_logLikelihood_internal(X, kappa, mu, alpha, Tau)
  if(beta > 0) {
    penLogLikelihood <- penLogLikelihood - beta*sum(abs(mu))
  }
  if(hard_assign) {
    cluster <- apply(Tau,1,which.max)
  } else {
    cluster <- NA
  }
  if(!save_tau) {
    Tau <- NULL
  }
  list(mu = mu, alpha = alpha, kappa = kappa, Tau = Tau,
       shared_kappa = shared_kappa,
       beta = beta, logL = logL, fpiter = fpiter,
       iter = length(logL),
       logLikelihood = LogLikelihood,
       penLogLikelihood = penLogLikelihood,
       completeLogLikelihood = completeLogLikelihood,
       status = status, n= nrow(X), cluster = cluster)
}

## find next beta
movMF_next_beta <- function(X, movMFResults, prec=sqrt(.Machine$double.eps)) {
  if(is.null(movMFResults$Tau)) {
    Tau <- movMF_soft_assign_internal(X, movMFResults$kappa, movMFResults$mu, movMFResults$alpha)$Tau
  } else {
    Tau <- movMFResults$Tau
  }
  R <- crossprod(Tau, X)
  R_abs <- abs(R)
  premu <- movMFResults$kappa*R_abs - movMFResults$beta
  non_zero <- premu>prec
  if(any(non_zero)) {
    movMFResults$beta+min(premu[non_zero])
  } else {
    Inf
  }
}

## compute the number of zero terms
movMF_nbzero <- function(movMR_res, prec=sqrt(.Machine$double.eps)) {
  sum(abs(movMR_res$mu)<=prec,na.rm=TRUE)
}

## compute ICs and need quantities
## EBIC is implemented with gamma=0.5
movMF_IC <- function(movMR_res, prec=sqrt(.Machine$double.eps)) {
  K <- length(movMR_res$kappa)
  n <- movMR_res$n
  d <- ncol(movMR_res$mu)
  ## kappa and alpha
  if(movMR_res$shared_kappa) {
    ## K-1 alpha and 1 kappa
    kap <- K
  } else {
    ## K-1 alpha and K kappa
    kap <- 2 * K - 1
  }

  mu_zeros <- abs(movMR_res$mu) < prec
  nb_mu_param_theo <- K * (d-1)
  ## we count non zero values per directional mean and keep at least one parameter
  nb_mu_param_non_zero <- sum(apply(mu_zeros, 1, function(x) max(1, d-sum(x)-1)))
  ## we count the variables with at least one non zero coordinate on a
  ## directional mean
  nb_mu_param_non_all_zero <- K * (sum(apply(mu_zeros, 2, sum) < K) - 1)

  list(dBIC=-2*movMR_res$logLikelihood + (nb_mu_param_theo + kap) * log(n),
       zBIC=-2*movMR_res$logLikelihood + (nb_mu_param_non_zero + kap) * log(n),
       vBIC=-2*movMR_res$logLikelihood + (nb_mu_param_non_all_zero + kap) * log(n),
       dAIC=-2*movMR_res$logLikelihood + (nb_mu_param_theo + kap) * 2,
       zAIC=-2*movMR_res$logLikelihood + (nb_mu_param_non_zero + kap) * 2,
       vAIC=-2*movMR_res$logLikelihood + (nb_mu_param_non_all_zero + kap) * 2,
       dRIC=-2*movMR_res$logLikelihood + (nb_mu_param_theo + kap) * 2*log(d),
       zRIC=-2*movMR_res$logLikelihood + (nb_mu_param_non_zero + kap) * 2*log(d),
       vRIC=-2*movMR_res$logLikelihood + (nb_mu_param_non_all_zero + kap) * 2*log(d),
       dRICc=-2*movMR_res$logLikelihood + (nb_mu_param_theo + kap) * 2*(log(d)+log(log(d))),
       zRICc=-2*movMR_res$logLikelihood + (nb_mu_param_non_zero + kap) * 2*(log(d)+log(log(d))),
       vRICc=-2*movMR_res$logLikelihood + (nb_mu_param_non_all_zero + kap) *
         2*(log(d)+log(log(d))),
       dEBIC=-2*movMR_res$logLikelihood + (nb_mu_param_theo + kap) * (log(n)+log(d)),
       zEBIC=-2*movMR_res$logLikelihood + (nb_mu_param_non_zero + kap) * (log(n)+log(d)),
       vEBIC=-2*movMR_res$logLikelihood + (nb_mu_param_non_all_zero + kap) * (log(n)+log(d)),
       counts=list(nb_non_zero=nb_mu_param_non_zero, nb_non_all_zero=nb_mu_param_non_all_zero))
}

## thresholding
movMF_threshold <- function(X, movMFResults, beta,
                            prec=sqrt(.Machine$double.eps)) {
  if(is.null(movMFResults$Tau)) {
    Tau <- movMF_soft_assign_internal(X, movMFResults$kappa, movMFResults$mu, movMFResults$alpha)$Tau
  } else {
    Tau <- movMFResults$Tau
  }
  R <- crossprod(Tau, X)
  R_abs <- abs(R)
  premu <- movMFResults$kappa*R_abs - beta
  premu[premu<=prec] <- 0
  premu <- premu*sign(R)
  wordspace::normalize.rows(premu)
}

movMF_extract_stats <- function(movMFResults) {
  res <- movMFResults
  res[c("mu", "alpha", "Tau", "kappa", "logL")] <- NULL
  res
}

keep_best_IC <- function(model, current, save_tau) {
  if(is.null(current)) {
    the_models <- list()
    criteria <- setdiff(names(model$IC), "counts")
    if(!save_tau) {
      model$Tau <- NULL
    }
    for(cr in criteria) {
      the_models[[cr]] <- model
    }
    list(IC=model$IC[criteria], models=the_models)
  } else {
    for(cr in names(current$IC)) {
      if(model$IC[[cr]] < current$IC[[cr]]) {
        current$models[[cr]] <- model
        current$IC[[cr]] <- model$IC[[cr]]
        if(!save_tau) {
          current$models[[cr]]$Tau <- NULL
        }
      }
    }
    current
  }
}


## extract more statistics from a model
movMF_extract_stats_gt <- function(movMFResults, membership) {
  res <- movMFResults
  res[c("mu", "alpha", "Tau", "kappa", "logL")] <- NULL
  res$ARI <- aricode::ARI(res$cluster, membership)
  res$NMI <- aricode::NMI(res$cluster, membership)
  res
}

movMF_extract_stats_ground_truth <- function(membership) {
  function(movMFResults) {
    movMF_extract_stats_gt(movMFResults, membership)
  }
}


## extract information from a path
movMF_path_summary <- function(movMFPath, with_gt=FALSE) {
  if(with_gt) {
    extractor <- function(x) {
      IC <- x$IC
      c(x$beta, x$logLikelihood, x$penLogLikelihood, x$ARI,
        x$NMI,
        x$sparsity, x$iter, x$fpiter, x$status,
        unlist(IC[1:(length(IC)-1)]), unlist(IC[["counts"]])
      )
    }
    the_names <- c("beta", "logLikelihood", "penLogLikelihood",
                   "ARI", "NMI", "sparsity", "iterations",
                   "internal_iterations", "status")
  } else {
    extractor <- function(x) {
      IC <- x$IC
      c(x$beta, x$logLikelihood, x$penLogLikelihood,
        x$sparsity, x$iter, x$fpiter, x$status,
        unlist(IC[1:(length(IC)-1)]), unlist(IC[["counts"]])
      )
    }
    the_names <- c("beta", "logLikelihood", "penLogLikelihood",
                   "sparsity", "iterations",
                   "internal_iterations", "status")
  }
  pre_res <- as.data.frame(t(sapply(movMFPath, extractor)))
  pre_names <- names(pre_res)
  pre_names[seq_along(the_names)] <- the_names
  names(pre_res) <- pre_names
  for(cl in c("iterations", "internal_iterations", "status",
              "nb_non_zero", "nb_non_all_zero")) {
    pre_res[[cl]] <- as.integer(pre_res[[cl]])
  }
  pre_res
}

