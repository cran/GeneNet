\name{dcor0}
\alias{dcor0}
\alias{pcor0}
\alias{rcor0}
\alias{ibeta}

\title{Distribution of the Vanishing Correlation Coefficient (rho=0) and Related Functions}

\usage{
dcor0(x, kappa, log=FALSE)
pcor0(q, kappa, lower.tail=TRUE, log.p=FALSE)
rcor0(n, kappa)
ibeta(z, a, b)
}

\arguments{
	\item{x,q}{vector of sample correlations}
	\item{kappa}{the degree of freedom of the distribution (= inverse variance)}
	\item{n}{number of values to generate. If n is a vector, length(n)
                values will be generated}
	\item{log, log.p}{logical vector; if TRUE, probabilities p are given as log(p)}
	\item{lower.tail}{logical vector; if TRUE (default), probabilities are \eqn{P[R <= r]},
	                  otherwise, \eqn{P[R > r]}}
        \item{a,b,z}{numeric vectors}
}

\description{
  Density and distribution function and a random number generator
  of Pearson's correlation coefficient assuming that there is no
  correlation present (\code{rho = 0}).
 
  The distribution has a single parameter, namely the degree of freedom \code{kappa} which
  equals the inverse of the variance of the distribution.
  
  
  The theoretical value of  \code{kappa} depends both on the sample size n and the number p
  of investigated variables.  If a simple correlation coefficient is considered (p=2)
  then the degree of freedom equals \code{kappa = n-1}.
  For p  variables and a partial correlation coefficient conditioned on p-2 variables
  \code{kappa = n-p+1}.
         
  The incomplete beta function \code{ibeta} is needed to compute the distribution function. For
  z=1 the incomplete beta function reduces to the beta function
  (i.e. \code{ibeta(1, a, b) = beta(a, b)}).
  
}
\details{
  For density and distribution functions (as well as a corresponding random number generator)
  of the correlation coefficient for arbitrary values of rho please refer to the
  \code{SuppDists} package by  Bob Wheeler \email{bwheeler@echip.com} (available on CRAN).
  Note that the parameter N in the \code{dPearson} function (and others in the  \code{SuppDists} package)
  corresponds to  \code{N=kappa+1}.
  
}

\value{
  The output values conform to the output from other such functions 
  in \R. \code{dcor0} gives the density and \code{pcor0} 
  the distribution function.
  
  The function \code{ibeta} returns a numeric value.
}

\seealso{\code{\link{cor}}, \code{\link{beta}}.}

\author{
   Korbinian Strimmer (\url{http://strimmerlab.org}).
}

\examples{
# load GeneNet library
library("GeneNet")

# distribution of r for various degrees of freedom
x <- seq(-1,1,0.01)
y1 <- dcor0(x, kappa=7)
y2 <- dcor0(x, kappa=15)
plot(x,y2,type="l", xlab="r", ylab="pdf",
  xlim=c(-1,1), ylim=c(0,2))
lines(x,y1)

# simulated data
r <- rcor0(1000, kappa=7)
hist(r, freq=FALSE, 
  xlim=c(-1,1), ylim=c(0,5))
lines(x,y1,type="l")

# distribution function
pcor0(-0.2, kappa=15)

# incomplete beta function
ibeta(0.4, 1, 3)
ibeta(1, 2, 3)
beta(2, 3)
}
\keyword{distribution}
