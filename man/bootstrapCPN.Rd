% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bootstrapCPN.R
\name{bootstrapCPN}
\alias{bootstrapCPN}
\title{Bootstrap cancer progression networks}
\usage{
bootstrapCPN(
  df,
  model = "CBN",
  algorithm = "DP",
  k = 3,
  epsilon = 0.025,
  N = 100
)
}
\arguments{
\item{df}{A dataframe or matrix of binary mutation data. The columns should be 
mutations and the rows should be patients.}

\item{model}{Whether to fit a network with "conjuctive" interactions (CBN), or
"disjunctive" interactions (DBN).}

\item{algorithm}{The algorithm used to infer the network. Default is
an optimal dynamic programming ("DP") algorithm, which is only feasible when the number of mutations is
smaller than 30. For larger datasets, use the genetic algorithm ("GA"), which 
gives an approximate solution.}

\item{k}{The in-degree bound on the estimated network.}

\item{epsilon}{Penalty term (for mutations not conforming to the estimated network).}

\item{N}{The number of bootstrap iterations.}
}
\value{
A list with components 
\itemize{
\item \code{edgelist} - A vector of length 2*k listing
the vertices in the k inferred edges.
\item \code{score} - The maximum likelihood value.
}
}
\description{
This function performs a nonparameteric
bootstrap to quanity the uncertainty in the 
maximum likelihood network.
}
\examples{
## Load the example data
data("example")

## Fit a conjunctive Bayesian network (CBN)
out <- fitCPN(example, epsilon=0.01)

## Plot the graph
bs <- boots

}
\references{
Nicol, PB, Coombes, KR, Deaver, C, Chkrebtii, O, Paul, S, Toland, AE, Asiaee, A. 
Oncogenetic network estimation with disjunctive Bayesian networks. Comp. Sys. Onco. 2021; 1:e1027.
 https://doi.org/10.1002/cso.21027
}
\author{
Phillip B. Nicol <philnicol740@gmail.com>
}
