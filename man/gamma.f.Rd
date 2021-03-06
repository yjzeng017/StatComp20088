% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GPbinary.R
\name{gamma.f}
\alias{gamma.f}
\title{Laplace approximation}
\usage{
\method{gamma}{f}(f, hyper, x1, x2, y)
}
\arguments{
\item{f}{the latent process f}

\item{hyper}{the hyper parameters involved in covariance matrix}

\item{x1}{the fitst covariate}

\item{x2}{the second covariate}

\item{y}{the response variable}
}
\value{
return the gradient
}
\description{
A function to approximate the marginal log-likelihood by using Laplace method
}
