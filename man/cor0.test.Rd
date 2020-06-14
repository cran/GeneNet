\name{cor0.test}
\alias{cor0.test}

\title{Test of Vanishing (Partial) Correlation}
\usage{
cor0.test(r, kappa, method=c("student", "dcor0", "ztransform"))
}
\arguments{
  \item{r}{observed correlation}
  \item{kappa}{degree of freedom of the null-distribution}
  \item{method}{method used to compute the p-value}
}
\description{
  \code{cor0.test} computes a p-value for the two-sided test with the null
  hypothesis H0: rho == 0 versus the alternative hypothesis HA: rho != 0.
  
  If \code{method="student"} is selected then the statistic 
  \code{t=r*sqrt((kappa-1)/(1-r*r))} is considered which under H0 is
  student-t distributed with \code{df=kappa-1}.  This method is exact.
  
  If \code{method="dcor0"} is selected then the p-value is computed
  directly from the null distribution of the (partial) correlation
  (see \code{\link[fdrtool]{dcor0}}).
  This method is also exact.
 
  If \code{method="ztransform"} is selected then the p-value is computed
  using the z-transform (see \code{\link{z.transform}}), i.e. using 
  a suitable chosen normal distribution.
  This method returns approximate p-values.
}

\value{
  A p-value.
}


\author{
  Juliane Sch\"afer and
  Korbinian Strimmer (\url{http://strimmerlab.org}).
}

\seealso{ \code{\link[fdrtool]{dcor0}},  \code{\link{kappa2n}}, \code{\link{z.transform}}.}

\examples{
# load GeneNet library
library("GeneNet")

# covariance matrix
m.cov <- rbind(
 c(3,1,1,0),
 c(1,3,0,1),
 c(1,0,2,0),
 c(0,1,0,2)
)

# compute partial correlations
m.pcor <- cor2pcor(m.cov)
m.pcor

# corresponding p-values 
# assuming a sample size of 25, i.e. kappa=22
kappa2n(22, 4)
cor0.test(m.pcor, kappa=22)
cor0.test(m.pcor, kappa=22) < 0.05

# p-values become smaller with larger r 
cor0.test(0.7, 12)
cor0.test(0.8, 12)
cor0.test(0.9, 12)

# comparison of various methods
cor0.test(0.2, 45, method="student")
cor0.test(0.2, 45, method="dcor0")
cor0.test(0.2, 45, method="ztransform")
}
\keyword{htest}
