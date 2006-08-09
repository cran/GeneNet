\name{z.transform}
\alias{z.transform}
\alias{hotelling.transform}

\title{Variance-Stabilizing Transformations of the Correlation Coefficient}

\usage{
z.transform(r)
hotelling.transform(r, kappa)
}

\arguments{
  \item{r}{vector of sample correlations}
  \item{kappa}{degrees of freedom of the distribution of the correlation coefficient}
}
\description{
  \code{z.transform} implements Fisher's (1921)  first-order and Hotelling's (1953) 
  second-order transformations to stabilize the distribution of the correlation coefficient.
  After the transformation the data follows approximately a
  normal distribution  with constant variance (i.e. independent of the mean). 
  
  The Fisher transformation is simply \code{z.transform(r) = atanh(r)}.
  
  
  Hotelling's transformation requires the specification of the degree of freedom \code{kappa} of
  the underlying distribution. This depends on the sample size n used to compute the
  sample correlation and whether simple ot partial correlation coefficients are considered.
  If there are p variables, with p-2 variables eliminated, the degree of freedom is \code{kappa=n-p+1}.
  (cf. also \code{\link{dcor0}}). 
}


\value{
  The vector of transformed sample correlation coefficients.
}
\author{
  Korbinian Strimmer (\url{http://strimmerlab.org}).
}
\references{
  Fisher, R.A. (1921). On the 'probable error' of a coefficient of correlation deduced from
  a small sample. \emph{Metron}, \bold{1}, 1--32.
 
  Hotelling, H. (1953). New light on the correlation coefficient and its transformation.
   \emph{J. Roy. Statist. Soc. B}, \bold{15}, 193--232.

}
\seealso{ \code{\link{dcor0}},  \code{\link{kappa2n}}.}

\examples{
# load GeneNet library
library("GeneNet")

# small example data set 
r <- c(-0.26074194, 0.47251437, 0.23957283,-0.02187209,-0.07699437,
       -0.03809433,-0.06010493, 0.01334491,-0.42383367,-0.25513041)

# transformed data
z1 <- z.transform(r)
z2 <- hotelling.transform(r,7)
z1
z2
}
\keyword{univar}
