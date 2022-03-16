# sparsevmf

<!-- badges: start -->
<!-- badges: end -->

This package implements L1 norm penalization for von Mises-Fisher distribution. Moreover to explore the trade-off between the sparsity term            and the likelihood one, we implement a path following algorithm.

Mathematical explanation for the penalization of von Mises-Fisher distribution will be soon available.

## Installation

You can install the released version of sparsevmf from [GitHub](https://github.com/FloFloB/sparsevmf) with:

``` r
devtools::install_github("FloFloB/sparsevmf")
```

## Example

Functions available in the package are:

* movMF_EM Provides a function to build and fit a penalized mixtures of von Mises-Fisher.
* movMF_beta_path Provides a function  the trade-off between the sparsity term and the likelihood one with a simple path following algorithm.
* movMF_initialisation_internal Provides a function to initialize from prototypes or from clusters a penalized mixtures of von Mises-Fisher.
* build_seed Provides a function to create seeds to have a fair comparison between algorithms.
* plot_data Provides a function to plot the data set according to the result of a mixture of von Mises-Fisher distribution.
* plot_proto Provides a function to plot the directional means of a mixture of von Mises-Fisher distribution.

First, load up your package with:
``` r
library(sparsevmf)
```

Libraries needed for the following example :
``` r
library(tibble)
library(skmeans)
library(Matrix)
```

Example of use of the package functions for comparison with spherical k-means.

``` r
## Data for example
cstr <- R.matlab::readMat("~/R_works/cstr.mat") # available here https://github.com/dbmovMFs/DirecCoclus/tree/master/Data

## Normalize data
cstr_data <- wordspace::normalize.rows(Matrix::Matrix(cstr$fea, sparse=TRUE))

## Improve representation of label
cstr$gnd=cstr$gnd[,1]

## Build seed. For this example, only 1 run with 10 init with ks = 4
ks=4
nb_run=1
nb_init=10
seeds=build_seed(cstr_data,ks=ks,nb_run = nb_run,nb_init=nb_init)

## We keep the best configuration in terms of the quality criterion of each algorithm:

cstr_seeds_current <- seeds[[ks]][[1]]
best_models <- list()
for(h in 1:nrow(cstr_seeds_current)) {
  pId <- cstr_seeds_current[h,]
  pX <- cstr_data[pId, , drop=FALSE]
  pCluster <- apply(tcrossprod(cstr_data, pX), 1, which.max)
  
  # s-kmeans
  cstr_sk <- skmeans(cstr_data, k=k,
                     method="pclust",control=list(start=pCluster))
  if(is.null(best_models[["sk"]])) {
    best_models[["sk"]] <- cstr_sk
  } else if(best_models[["sk"]]$value > cstr_sk$value) {
    best_models[["sk"]] <- cstr_sk
  }
  # movMF with non-shared kappa
   my_movMF_conf <- movMF_initialisation_internal(cstr_data, pX, pCluster, shared_kappa=FALSE)
   cstr_my_movMF <- movMF_EM(cstr_data, K=k, beta=0, Theta=my_movMF_conf)
   if(is.null(best_models[["movMF"]])) {
     best_models[["movMF"]] <- cstr_my_movMF
   } else if(best_models[["movMF"]]$logLikelihood < cstr_my_movMF$logLikelihood) {
     best_models[["movMF"]] <- cstr_my_movMF
   }
  
  # movMF with shared kappa
  my_movMF_conf_shared <- movMF_initialisation_internal(cstr_data, pX, pCluster, shared_kappa=TRUE)
  cstr_my_movMF_shared <- movMF_EM(cstr_data, K=k, beta=0,
                                   Theta=my_movMF_conf_shared,
                                   shared_kappa=TRUE)
  if(is.null(best_models[["movMF_shared"]])) {
    best_models[["movMF_shared"]] <- cstr_my_movMF_shared
  } else if(best_models[["movMF_shared"]]$logLikelihood < cstr_my_movMF_shared$logLikelihood) {
    best_models[["movMF_shared"]] <- cstr_my_movMF_shared
  }
}

# Compare solution in terms of ARI and NMI
one_run <- tibble(method=c("sk", "movMF",
                           "movMF_shared"),
                  criterion=c(best_models[["sk"]]$value,
                              best_models[["movMF"]]$logLikelihood,
                              best_models[["movMF_shared"]]$logLikelihood),
                  ARI=rep(NA,3),
                  NMI=rep(NA,3))
for(j in 1:3) {
  one_run[j,3] <- aricode::ARI(best_models[[j]]$cluster, cstr$gnd)
  one_run[j,4] <- aricode::NMI(best_models[[j]]$cluster, cstr$gnd)
}

```

Example of the path following approach starting from the best model with shared kappa:

``` r
cstr_path <- movMF_beta_path(cstr_data, K=ks,
                            Theta = best_models[["movMF_shared"]],
                            shared_kappa=TRUE,
                            min_rel_inc=0.01,
                            nb_beta=1000, save_path="statistics", verbose=1)
```

Plotting the model's directional means selected by RICc:

``` r
plot_proto(cstr_path$best_IC$models$zRICc,with_color = TRUE)
```

Plotting the data as the model's directional means selected by RICc:

``` r
plot_data(cstr_path$best_IC$models$zRICc,cstr_data,with_color = TRUE)
```

