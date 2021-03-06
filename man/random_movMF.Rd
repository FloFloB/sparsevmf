% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/random_movmf.R
\name{random_movMF}
\alias{random_movMF}
\title{Provides a function to produce a random movMF model with controlled separation between the clusters.}
\usage{
random_movMF(
  n,
  d,
  K,
  alpha = rep(1/K, K),
  kappa = NULL,
  separation = 0.05,
  oversample = 20,
  sparsify_mu = 0,
  verbose = FALSE
)
}
\arguments{
\item{n}{Number of samples.}

\item{d}{Dimension of samples.}

\item{K}{Number of components.}

\item{alpha}{Mixture proportions (by default same proportion).}

\item{kappa}{Value(s) of kappa.}

\item{separation}{Separation between components.}

\item{oversample}{Sample to extract well separated vectors.}

\item{sparsify_mu}{Parameter to sparsify the directional means.}

\item{verbose}{Verbose.}
}
\value{
returns a list including the following attributes:
\itemize{
\item X: samples
\item mu: directional means.
\item alpha: mixture proportions.
\item kappa: concentration parameters.
\item membership: to which cluster the sample belongs to.
\item separation: the separation between components.
\item overlap: overlap between components.
\item sparsity: sparsity of directional means.
}
}
\description{
Provides a function to produce a random movMF model with controlled separation between the clusters.
}
