\name{cor.fit.mixture}
\alias{cor.fit.mixture}

\title{Fit Mixture Distribution to Estimated Correlation Coefficients}
\usage{
cor.fit.mixture(r, df=7, plot.locfdr=0)
}
\arguments{
  \item{r}{vector of estimated correlations}
  \item{df}{degrees of freedom for the spline fitting the density}
  \item{plot.locfdr}{controls plot option in \code{locfdr} }
}
\description{

  \code{cor.fit.mixture} fits a mixture model 
  
    \code{f(r) = eta0 dcor0(r, kappa) + (1-eta0) etaA fA}
  
  to a vector of empirical partial correlation coefficients.
  This allows to estimate both the degree of freedom \code{kappa} in the
  null-distribution and the proportion eta0 of null r-values.  The alternative distribution
  may be any arbitrary nonparametric distribution which vanishes for values of r near the center r=0. 
  
  \code{cor.fit.mixture} also computes

    \code{etaA fA/f(r)}, 

  i.e. the (empirical) posterior probability that the true correlation is non-zero given the empirical
  correlation r, the degree of freedom of the null-distribution kappa, and the prior eta0
  for the null-distribution. 
}

\details{
  The above functions are useful to determine the null-distribution of edges in a sparse graphical
  Gaussian model, see Schaefer and Strimmer (2005a,b) for more details and an application to infer 
  genetic networks from microarray data.
  
  For details on how to fit the empirical null distribution while at the same time non-parametrically
  estimating the alternative hypothesis see Efron (2004) and the associated R package \code{locfdr}.
}

\value{
  A list object with the following components:
  \item{kappa}{the degree of freedom of the null distribution (see \code{\link{dcor0}})} 
  \item{eta0}{the prior for the null distribution, i.e. the proportion of null r-values}
  \item{prob.nonzero}{empirical posterior probability that the observed correlations
                      are non-zero.}
}


\author{
  Juliane Schaefer (\url{http://www.stat.math.ethz.ch/~schaefer/}) and
  Korbinian Strimmer (\url{http://strimmerlab.org}).
}
\references{
  Efron, B. (2004). Large-scale simulataneous hypothesis testing: the choice
  of a null hypothesis. \emph{JASA} \bold{99}:96-104.

  Schaefer, J., and Strimmer, K. (2005a).  An empirical Bayes approach to inferring
  large-scale gene association networks. \emph{Bioinformatics} \bold{21}:754-764.
  
 Schaefer, J., and Strimmer, K. (2005).  A shrinkage approach to large-scale
 covariance estimation and implications for functional genomics. 
 Statist. Appl. Genet. Mol. Biol. \bold{4}:32.
 (\url{http://www.bepress.com/sagmb/vol4/iss1/art32/})  
  
}

\seealso{\code{\link{dcor0}}, 
\code{\link{kappa2n}}, \code{\link[fdrtool]{fdr.estimate.eta0}}.}

\examples{
# load GeneNet library
library("GeneNet")

# simulate mixture distribution
r <- rcor0(700, kappa=10)
u1 <- runif(200, min=-1, max=0.5)
u2 <- runif(200, min=0.5, max=1)
rc <- c(r, u1, u2)

# estimate kappa and eta0 (=7/9)
c1 <- cor.fit.mixture(r)
c1$eta0
c1$kappa
c2 <- cor.fit.mixture(rc) 
c2$eta0
c2$kappa

}
\keyword{htest}
