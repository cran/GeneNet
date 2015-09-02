\name{arth800}
\alias{arth800}
\alias{arth800.descr}
\alias{arth800.expr}
\alias{arth800.mexpr}
\alias{arth800.name}
\alias{arth800.probe}
\alias{arth800.symbol}

\title{Time Series Expression Data for 800 Arabidopsis Thaliana Genes}

\description{
  This data set describes the temporal expression of 800 genes
  of \emph{A. thaliana} during the diurnal cycle.   The 800 genes are
  a subset of the data presented in Smith et al. (2004) and were
  selected for periodicity according to the method implemented
  in the R package GeneCycle (\url{http://cran.r-project.org/package=GeneCycle}).
}

\usage{
data(arth800)
}

\format{
  \code{arth800.expr} is a \code{\link[longitudinal]{longitudinal}} 
  object with repetitions, and contains the log2 transformed expression data. 

  \code{arth800.mexpr} is a \code{\link[longitudinal]{longitudinal}} 
  object, and contains the mean expression levels of \code{arth800.expr}.

  \code{arth800.descr}, \code{arth800.name}, \code{arth800.probe}, 
  \code{arth800.symbol} are vectors containing additional information about each gene.
}

\source{
  The micoarray experiments were performed in the laboratory of S. Smith (Edinburgh).
  The data are available from the 
  NASCArrays database (\url{http://affymetrix.arabidopsis.info/} under 
  experiment reference number NASCARRAYS-60.
}


\references{
Smith et al. 2004. Diurnal changes in the transcriptom encoding 
enzymes of starch metabolism provide evidence for both transcriptional
and posttranscriptional regulation of starch metabolism in Arabidopsis 
leaves.  Plant Physiol. 136: 2687-2699
}

\examples{
# load GeneNet library
library("GeneNet")

# load data set
data(arth800)

is.longitudinal(arth800.expr)
summary(arth800.expr)

# plot first nine time series
plot(arth800.expr, 1:9)
}
\keyword{datasets}

