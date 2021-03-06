% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distributions.R
\name{greta-distributions}
\alias{greta-distributions}
\alias{normal}
\alias{lognormal}
\alias{bernoulli}
\alias{binomial}
\alias{negative_binomial}
\alias{poisson}
\alias{gamma}
\alias{exponential}
\alias{student}
\alias{beta}
\alias{free}
\alias{flat}
\alias{multivariate_normal}
\alias{wishart}
\title{greta probability distributions}
\usage{
normal(mean, sd, dim = NULL)

lognormal(meanlog, sdlog, dim = NULL)

bernoulli(prob, dim = NULL)

binomial(size, prob, dim = NULL)

negative_binomial(size, prob, dim = NULL)

poisson(lambda, dim = NULL)

gamma(shape, rate, dim = NULL)

exponential(rate, dim = NULL)

student(df, ncp, dim = NULL)

beta(shape1, shape2, dim = NULL)

free(dim = 1)

flat(range, dim = 1)

multivariate_normal(mean, Sigma, dim)

wishart(df, Sigma, dim)
}
\arguments{
\item{mean, meanlog, ncp}{unconstrained parameters}

\item{sd, sdlog, size, lambda, shape, rate, df, shape1, shape2}{positive parameters}

\item{dim}{the dimensions of the variable. For univariate distributions this
can be greater than 1 to represent multiple independent variables. For
multivariate distributions this cannot be smaller than 2.}

\item{prob}{probability parameter (\code{0 < prob < 1})}

\item{range}{a finite, length 2 numeric vector giving the range of values to
which \code{flat} distributions are constrained. The first element must
be lower than the second.}

\item{Sigma}{positive definite variance-covariance matrix parameter}
}
\description{
These probability distributions can be used to define random
  variables in a greta model. They return a greta array object that can be
  combined with other greta arrays to construct a model.
}
\details{
Most of these distributions have non-uniform probability densities,
  however the distributions \code{flat} and \code{free} do not. These can
  therefore be used as parameters in likelihood (rather than Bayesian)
  inference.

  The discrete probability distributions (\code{bernoulli}, \code{binomial},
  \code{negative_binomial}, \code{poisson}) can be used as likelihoods, but
  not as unknown variables.

  Wherever possible, the parameterisation of these distributions matches the
  those in the \code{stats} package. E.g. for the parameterisation of
  \code{negative_binomial()}, see \code{\link{dnbinom}}.
}
\examples{

# a fixed distribution, e.g. for a prior
mu = normal(0, 1)

# an unconstrained, positive parameter sigma
log_sigma = free()
sigma = exp(log_sigma)

# a hierarchical distribution
theta = normal(mu, lognormal(0, 1))

# a vector of 3 variables drawn from the same hierarchical distribution
thetas = normal(mu, sigma, dim = 3)

# a matrix of 12 variables drawn from the same hierarchical distribution
thetas = normal(mu, sigma, dim = c(3, 4))

# a constrained variable with no density (e.g. for a constrained likelihood model)
theta = flat(c(1, 5))

# a multivariate normal variable, with correlation between two elements
Sig <- diag(4)
Sig[3, 4] <- Sig[4, 3] <- 0.6
theta = multivariate_normal(rep(mu, 4), Sig, dim = 4)

# a Wishart variable with the same covariance parameter
theta = wishart(df = 5, Sigma = Sig, dim = 4)
}
