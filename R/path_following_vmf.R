#' Provides a function  the trade-off between the sparsity term and the likelihood one with a simple path following algorithm.
#' @import wordspace
#' @import skmeans
#' @import aricode
#' @import movMF
#' @import Matrix
#' @param X Data : matrix with column as features.
#' @param K The number of mixture components (or clusters).
#' @param beta Penalty parameter.
#' @param Theta Theta is the initial parameter list (alpha, mu, kappa).
#' @param shared_kappa If all components share kappa or not.
#' @param restart Number of restart.
#' @param nb_beta number of steps of the path following approach.
#' @param maxiter The maximum number of iteration.
#' @param maxiterfp The maximum number of iteration for the fix poit strategy.
#' @param prec Precision.
#' @param kappamax Kappa max value.
#' @param mode Method to initialize: "random", "spread", "skmeans".
#' @param sk_runs Number of runs for the initialization for skmeans.
#' @param interpolate interpolate paramters.
#' @param min_rel_inc Minimum increment for beta.
#' @param save_tau Save responsibilities Tau.
#' @param hard_assign Hard assignment.
#' @param save_path save information about the path with 3 options "full", "statistics", "nothing".
#' @param statistics by default internal function movMF_extract_stats to extract ARI and NMI will be used.
#' @param verbose Verbose.
#' @md
#' @return returns two lists depending on the statistic parameter chosen:
#' * First list is composed by the best models along the path selected depending on the criterion.
#' * Second list is a list of models or their principal statistics or only the dense model.
#' @export

movMF_beta_path <- function(X, K, Theta=NULL, shared_kappa=FALSE, restart=NULL,
                            nb_beta=NULL, maxiter=1000, maxiterfp=10,
                            prec=sqrt(.Machine$double.eps), kappamax=1e6,
                            mode="random", sk_runs=1,
                            interpolate = 1, min_rel_inc=0, save_tau=FALSE,
                            hard_assign = TRUE,
                            save_path = c("full", "statistics", "nothing"),
                            statistics = movMF_extract_stats,
                            verbose=0) {
  save_path <- match.arg(save_path)
  verbose <- as.integer(verbose)
  if(!is.null(Theta) & !is.null(restart)) {
    stop("Incompatible options")
  }
  if(K>=nrow(X)){
    stop("More clusters than distinc samples (rows)")
  }
  if(!is.logical(shared_kappa)){
    stop("shared_kappa must be a boolean")
  }
  if(!is.logical(save_tau)){
    stop("save_tau must be a boolean")
  }
  if(!is.logical(hard_assign)){
    stop("hard_assign must be a boolean")
  }
  if(!mode %in% c("random", "spread", "skmeans")){
    stop("mode must be choosed between random, spread or skmeans")
  }
  if(!save_path %in% c("full", "statistics", "nothing")){
    stop("save_path must be choosed between full, statistics or nothing")
  }
  if(!is.numeric(sk_runs) || sk_runs<1){
    stop("sk_runs must be numeric and at least equal to 1")
  }
  if(!is.numeric(maxiter) || maxiter<1){
    stop("maxiter must be numeric and at least equal to 1")
  }
  if(!is.numeric(maxiter) || maxiter<1){
    stop("maxiter must be numeric and at least equal to 1")
  }
  if(!is.numeric(prec) || prec<0){
    stop("prec must be numeric and greater than 0")
  }
  if(!is.numeric(kappamax) || kappamax<0){
    stop("kappamax must be numeric and greater than 0")
  }
  if(!is.numeric(interpolate) || interpolate<1){
    stop("interpolate must be numeric and at least equal to 1")
  }
  if(!is.numeric( min_rel_inc) ||  min_rel_inc<0){
    stop("min_rel_inc must be numeric and greater or equal to 0")
  }


  ## initial configuration
  if(verbose > 0) {
    cat("Initial configuration\n")
  }
  base_model <- movMF_EM(X, K, beta=0, Theta=Theta, shared_kappa=shared_kappa,
                         restart=restart, maxiter=maxiter,
                         maxiterfp=maxiterfp,
                         prec=prec, kappamax=kappamax, mode=mode,
                         sk_runs=sk_runs, save_tau=TRUE, hard_assign = hard_assign,
                         verbose=verbose-1)
  current_model <- base_model
  best_IC <- keep_best_IC(base_model, NULL, save_tau)
  if(save_path=="full") {
    models <- list()
    models[[1]] <- base_model
    if(!save_tau) {
      models[[1]]$Tau <- NULL
    }
  } else if(save_path=="statistics") {
    stats <- list()
    stats[[1]] <- statistics(base_model)
  }
  i <- 1
  need_beta <- TRUE
  last_beta <- 0
  maximum_sparsity <- K*(ncol(X)-1)
  if(is.null(nb_beta)) {
    nb_beta <- 2*maximum_sparsity
  }
  min_rel_inc <- max(min_rel_inc, prec)
  while(i <= nb_beta) {
    ## test for convergence (or lack thereoff)
    if(current_model$converged == 0 ||
       current_model$sparsity==maximum_sparsity) {
      break
    }
    if(need_beta) {
      need_beta <- FALSE
      next_beta <- movMF_next_beta(X, current_model)
      if(is.infinite(next_beta)) {
        break
      }
      if(next_beta - last_beta < min_rel_inc*last_beta) {
        next_beta <- next_beta + min_rel_inc*last_beta
      }
      j <- 1
      if(verbose > 0 && interpolate>1) {
        cat("Computing next beta", next_beta, "\n")
      }
      betas <- seq(last_beta, next_beta, length.out=interpolate+1)
    }
    current_beta <- betas[j+1]
    if(verbose > 0) {
      cat("Beta path step:", i, "internal step:", j, "beta:",
          current_beta, " \n")
    }
    ## threshold small values
    current_model$mu <- movMF_threshold(X, current_model, current_beta, prec)

    ## no more restart!
    new_model <- movMF_EM(X, K, current_beta, current_model, shared_kappa=shared_kappa,
                          maxiter=maxiter, maxiterfp=maxiterfp,
                          prec=prec, kappamax=kappamax,
                          save_tau=TRUE, hard_assign = hard_assign, verbose=verbose-1)
    if(verbose > 0) {
      cat("\tSparsity", new_model$sparsity, "\n")
    }
    if(new_model$status==0) {
      best_IC <- keep_best_IC(new_model, best_IC, save_tau)
    }
    current_model <- new_model
    if(save_path=="full") {
      models[[i+1]] <- new_model
      if(!save_tau) {
        models[[i+1]]$Tau <- NULL
      }
    } else if(save_path=="statistics") {
      stats[[i]] <- statistics(new_model)
    }
    i <- i + 1
    if(j >= interpolate) {
      need_beta <- TRUE
      last_beta <- next_beta
    } else {
      j <- j + 1
    }
  }
  result <- list(best_IC=best_IC)
  result$total_iterations_init <- base_model$total_iterations_init
  if(save_path=="full") {
    result$models <- models
  } else if(save_path=="statistics") {
    result$stats <- stats
  }
  if(save_path == "statistics" || save_path == "nothing") {
    result$dense <- base_model
    if(!save_tau) {
      result$dense$Tau <- NULL
    }
  }
  result
}
