\name{ggm.estimate.pcor}
\alias{ggm.estimate.pcor}

\title{Graphical Gaussian Models: Small Sample Estimation of Partial Correlation}

\usage{
ggm.estimate.pcor(x, method = c("static", "dynamic"), \dots)
}

\arguments{
  \item{x}{data matrix (each rows corresponds to one multivariate observation)}
  \item{method}{method used to estimate the partial correlation matrix.
     Available options are "static" (the default) and "dynamic" - both
     are shrinkage methods.}
  \item{\dots}{options passed to \code{\link[corpcor]{pcor.shrink}}
         and to \code{\link[longitudinal]{dyn.pcor}}.}
   
}
\description{
  \code{ggm.estimate.pcor} offers an interface to two related shrinkage estimators
  of partial correlation. Both are fast, statistically efficient, and can be used for
  analyzing small sample data.  
  
  The default method "statics" employs the function 
  \code{\link[corpcor]{pcor.shrink}} whereas the "dynamic" method relies on 
  \code{\link[longitudinal]{dyn.pcor}}.  The difference between the two estimators
  is that the latter takes  the spacings between time points into account
  if the input are multiple time course data (these must be provided as 
  \code{\link[longitudinal]{longitudinal}} object).  
}

\details{
  For details of the shrinkage estimators we refer to Opgen-Rhein and Strimmer (2006a,b)
  and Sch\"afer and Strimmer (2005), as well as to the manual pages of 
  \code{\link[corpcor]{pcor.shrink}} and \code{\link[longitudinal]{dyn.pcor}}.
  
  Previously, this function offered several furthers options.
  The old option called "shrinkage" corresponds to the present "static" option.
  The other old options "observed.pcor", "partial.bagged.cor", and "bagged.pcor"
  are now considered obselete and have been removed.
}

\value{
  An estimated partial correlation matrix.
}
\author{
  Rainer Opgen-Rhein,
  Juliane Sch\"afer, and
  Korbinian Strimmer (\url{http://strimmerlab.org}).
}
\references{
 Opgen-Rhein, R., and K. Strimmer. 2006a. Inferring gene dependency networks from 
 genomic longitudinal data: a functional data approach.
 REVSTAT \bold{4}:53-65. (\url{http://http://www.ine.pt/revstat/}) 

 Opgen-Rhein, R., and K. Strimmer. 2006b. Using regularized dynamic correlation 
 to infer gene dependency networks from time-series microarray data. 
 The 4th International Workshop on Computational Systems Biology,
 WCSB 2006 (June 12-13, 2006, Tampere, Finland). 
 (\url{http://www.cs.tut.fi/wcsb06/})  
  
 Sch\"afer, J., and Strimmer, K. (2005).  A shrinkage approach to large-scale
 covariance estimation and implications for functional genomics. 
 Statist. Appl. Genet. Mol. Biol. \bold{4}:32.
 (\url{http://www.bepress.com/sagmb/vol4/iss1/art32/})  
  
}
\seealso{\code{\link{ggm.simulate.data}}, \code{\link{ggm.estimate.pcor}},
\code{\link[corpcor]{pcor.shrink}}, and \code{\link[longitudinal]{dyn.pcor}}.}

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


\keyword{htest}
