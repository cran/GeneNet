\name{ggm.simulate.data}
\alias{ggm.simulate.data}

\title{Graphical Gaussian Models: Simulation of Data}

\usage{
ggm.simulate.data(sample.size, pcor)
}

\arguments{
  \item{sample.size}{sample size of simulated data set}
  \item{pcor}{partial correlation matrix}
}
\description{
  \code{ggm.simulate.data} takes a positive definite partial correlation matrix and
  generates an i.i.d. sample from the corresponding standard multinormal distribution 
  (with mean 0 and variance 1).  If the input matrix \code{pcor} is not positive definite
  an error is thrown.
}


\value{
  A multinormal data matrix.
}
\author{
  Juliane Sch\"afer and
  Korbinian Strimmer (\url{https://strimmerlab.github.io}).
}
\references{
  Sch\"afer, J., and Strimmer, K. (2005).  An empirical Bayes approach to inferring
  large-scale gene association networks. \emph{Bioinformatics} \bold{21}:754-764.
}
\seealso{\code{\link{ggm.simulate.pcor}}, \code{\link{ggm.estimate.pcor}}.}

\examples{

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
\keyword{multivariate}
