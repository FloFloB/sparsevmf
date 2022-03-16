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
