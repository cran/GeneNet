\name{kappa2n}
\alias{kappa2n}
\alias{n2kappa}

\title{Relationship Between Sample Size and the Degree of Freedom of Correlation Distribution}
\usage{
kappa2n(kappa, p=2)
n2kappa(n, p=2)
}
\arguments{
  \item{kappa}{degree of freedom}
  \item{p}{number of variables (p=2 corresponds to simple correlation)}
  \item{n}{sample size}
}
\description{
  The function \code{kappa2n} returns the sample size that
  corresponds to a given degree of freedom kappa, whereas \code{n2kappa}
  converts sample size to the corresponding degree of freedom. 
}

\details{ 
        
  The degree of freedom kappa of the sample distribution of the empirical correlation
  coefficient depends both on the sample size n and the number p of investigated variables, 
  i.e. whether simple or partial correlation coefficients are being considered.
  For p=2 (simple correlation coefficient) the degree of freedom equals \code{kappa = n-1},
  whereas for arbitrary p (with p-2 variables eliminated in the partial correlation coefficient)
  \code{kappa = n-p+1} (see also \code{\link[fdrtool]{dcor0}}).
 

}

\value{
  The sample size n corresponding to a given kappa, or the degree of freedom
  kappa corresponding to a given p.
}


\author{
  Juliane Sch\"afer and
  Korbinian Strimmer (\url{https://strimmerlab.github.io}).
}

\seealso{\code{\link[fdrtool]{dcor0}}.}

\examples{
# load GeneNet library
library("GeneNet")


# sample sizes corresponding to kappa=7
kappa2n(7)     # simple correlation
kappa2n(7, 40) # partial correlation with p=40 variables

# degree of freedom corresponding to n=100
n2kappa(100)
n2kappa(100,40)
}
\keyword{univar}
