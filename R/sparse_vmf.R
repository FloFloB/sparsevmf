#' Provides a function to build and fit a penalized mixtures of von Mises-Fisher.
#' @import wordspace
#' @import skmeans
#' @import aricode
#' @import movMF
#' @import Matrix
#' @param X Data : matrix with column as features
#' @param K The number of mixture components (or clusters).
#' @param beta Penalty parameter.
#' @param Theta Theta is the initial parameter list.
#' @param shared_kappa If all components share kappa or not.
#' @param restart Number of restart.
#' @param maxiter The maximum number of iteration.
#' @param maxiterfp The maximum number of iteration for the fix poit strategy.
#' @param prec Precision.
#' @param kappamax Kappa max value.
#' @param mode Method to initialize: "random", "spread", "skmeans".
#' @param sk_runs Number of runs for the initialization for skmeans.
#' @param save_tau Save responsibilities Tau.
#' @param hard_assign Hard assignment.
#' @param verbose Verbose.
#' @md
#' @return returns a list including the following attributes:
#' * mu: directional means.
#' * alpha: mixture proportions.
#' * kappa: concentration parameters.
#' * Tau: responsibilities.
#' * shared_kappa: TRUE if the model shares kappa.
#' * beta: penalty parameter.
#' * logL: Log-likelihood.
#' * fpiter: fix point iterations.
#' * iter: iterations.
#' * logLikelihood: best LogLikelihood.
#' * penLogLikelihood: best penalized LogLikelihood.
#' * completeLogLikelihood: best complete LogLikelihood.
#' * status:
#'      + 0: converged,
#'      + 1: maxiter reached,
#'      + 2: not converging to a K cluster solution,
#'      + 3: infinite variance.
#'  * n: number of rows of X.
#'  * cluster: hard assignment.
#'  * converged: nb of status converged.
#'  * infinite: nb of status infinite.
#'  * empty: nb of status not converging to a K cluster solution.
#'  * itermax: iter maximum.
#'  * total_iterations: nb total iteration.
#'  * total_internal_iterations: nb total iteration for fix point strategy.
#'  * all_logLikelihood: all log Likelihood.
#'  * all_penLogLikelihood: all penalized log Likelihood.
#'  * IC: List of Information criteria BIC, AIC, RIC, RICc, EBIC with :
#'      + d: we keep all parameters,
#'      + z: we count non zero values per directional mean and keep at least one parameter,
#'      + v: we count the variables with at least one non zero coordinate on a directional mean.
#'  * sparsity: sparsity of mu.

#' @export

movMF_EM <- function(X, K, beta, Theta=NULL, shared_kappa=FALSE, restart=NULL, maxiter=1000, maxiterfp=10,
                     prec=sqrt(.Machine$double.eps), kappamax=1e6,
                     mode="random", sk_runs=1, save_tau=TRUE,
                     hard_assign=TRUE, verbose=0) {
  if(!is.null(Theta) & !is.null(restart)) {
    stop("Incompatible options")
  }
  verbose <- as.integer(verbose)
  if(is.null(Theta)) {
    if(is.null(restart)) {
      restart <- 1
    }
    irestart <- as.integer(restart)
    if(length(irestart)!=1) {
      stop(paste("Unadapted value for restart:", restart))
    }
    if(irestart<1) {
      stop("restart must be 1 or more")
    }
    bestLike <- -Inf
    bestModel <- NULL
    converged <- 0
    empty <- 0
    itermax <- 0
    infinite <- 0
    all_logLikelihood <- rep(NA, irestart)
    all_penLogLikelihood <- rep(NA, irestart)
    total_iter <- 0
    total_internal_iter <- 0
    for(i in 1:irestart) {
      if(verbose > 0) {
        cat("Initial configuration #", i, "\n")
      }
      Theta <- movMF_initialisation(X, K, shared_kappa=shared_kappa,
                                    mode=mode, sk_runs=sk_runs)
      currentEM <- movMF_EM_internal(X, beta, Theta, shared_kappa=shared_kappa,
                                     maxiter=maxiter,
                                     maxiterfp=maxiterfp,
                                     prec=prec, kappamax=kappamax,
                                     save_tau=save_tau, hard_assign=hard_assign,
                                     verbose=verbose-1)
      all_logLikelihood[i] <- currentEM$logLikelihood
      all_penLogLikelihood[i] <- currentEM$penLogLikelihood
      total_iter <- total_iter + currentEM$iter
      total_internal_iter <- total_internal_iter + currentEM$fpiter
      if(currentEM$status==0 | currentEM$status==1 ) {
        if(currentEM$penLogLikelihood > bestLike) {
          bestLike <- currentEM$penLogLikelihood
          bestModel <- currentEM
        }
        if(currentEM$status==1) {
          if(verbose > 0) {
            print("not converged, increase maxiter")
          }
          itermax <- itermax + 1
        } else {
          converged <- converged + 1
        }
      } else {
        if(currentEM$status==2) {
          empty <- empty + 1
        } else {
          infinite <- infinite + 1
        }
        if(verbose > 0) {
          print("discarding non converging run")
        }
      }
    }
  } else {
    if(verbose > 0) {
      cat("Single initiation configuration\n")
    }
    bestModel <- movMF_EM_internal(X, beta, Theta,
                                   shared_kappa=shared_kappa,
                                   maxiter=maxiter,
                                   maxiterfp=maxiterfp,
                                   prec=prec, kappamax=kappamax,
                                   save_tau=save_tau,
                                   hard_assign=hard_assign, verbose=verbose-1)
    converged <- 0
    infinite <- 0
    empty <- 0
    itermax <- 0
    all_logLikelihood <- bestModel$logLikelihood
    all_penLogLikelihood <- bestModel$penLogLikelihood
    total_iter <- bestModel$iter
    total_internal_iter <- bestModel$fpiter
    if(bestModel$status==0) {
      converged <- 1
    } else if(bestModel$status==1) {
      itermax <- 1
    } else if(bestModel$status==2) {
      empty <- 1
    } else {
      infinite <- 1
    }
  }
  bestModel$converged <- converged
  bestModel$infinite <- infinite
  bestModel$empty <- empty
  bestModel$itermax <- itermax
  bestModel$total_iterations <- total_iter
  bestModel$total_internal_iterations <- total_internal_iter
  bestModel$all_logLikelihood <- all_logLikelihood
  bestModel$all_penLogLikelihood <- all_penLogLikelihood
  ## TODO: this might break when the model did not converge
  bestModel$IC <- movMF_IC(bestModel)
  bestModel$sparsity <- movMF_nbzero(bestModel)/length(bestModel$mu)
  bestModel
}
