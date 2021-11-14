\name{ggm.simulate.pcor}
\alias{ggm.simulate.pcor}

\title{Graphical Gaussian Models: Simulation of Networks}

\usage{
ggm.simulate.pcor(num.nodes, etaA=0.05, stdprec=FALSE)
}

\arguments{
  \item{num.nodes}{number of nodes in the network}
  \item{etaA}{fraction of edges with non-zero partial correlation (default: 0.05)}
  \item{stdprec}{return standardised precision matrix, rather than matrix of partial correlations}
}
\description{
  \code{ggm.simulate.pcor} generates a random matrix of partial correlations that 
  corresponds to a GGM network of a given size (\code{num.nodes})
  with a specified fraction of non-zero edges.  The diagonal entries of the output matrix contain 1.

  If \code{stdprec=TRUE} then the standardised precision matrix is returned instead of the patrix of partial
  correlations.
}

\details{
  The simulation of the partial correlation matrix works by generating a diagonally dominant matrix
  as a positive definite precision matrix (inverse covariance matrix), which is 
  subsequently standardized and transformed into the matrix of partial correlations.
  For the full algorithm see Sch\"afer and Strimmer (2005).
}
\value{
  A positive partial correlation matrix (diagonal 1) with positive definite underlying precision matrix.

  If \code{stdprec=TRUE} then the standardised precision matrix is returned instead.
}
\author{
  Juliane Sch\"afer and
  Korbinian Strimmer (\url{https://strimmerlab.github.io}).
}
\references{
  Sch\"afer, J., and Strimmer, K. (2005).  An empirical Bayes approach to inferring
  large-scale gene association networks. \emph{Bioinformatics} \bold{21}:754-764.
}
\seealso{\code{\link{ggm.simulate.data}},\code{\link{ggm.estimate.pcor}}.}

\examples{
\dontrun{

# load GeneNet library
library("GeneNet")

# generate random network with 40 nodes 
# it contains 780=40*39/2 edges of which 5 percent (=39) are non-zero
true.pcor <- ggm.simulate.pcor(40)
  
# simulate data set with 40 observations
m.sim <- ggm.simulate.data(40, true.pcor)

# simple estimate of partial correlations
estimated.pcor <- cor2pcor( cor(m.sim) )

# comparison of estimated and true values
sum((true.pcor-estimated.pcor)^2)

# a slightly better estimate ...
estimated.pcor.2 <- ggm.estimate.pcor(m.sim)
sum((true.pcor-estimated.pcor.2)^2)

}
}
\keyword{multivariate}
