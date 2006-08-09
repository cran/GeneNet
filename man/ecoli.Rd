\name{ecoli}
\alias{ecoli}

\title{Microarray Time Series Data for 102 E. Coli Genes 
       Genes}

\description{
  This data set describes the temporal expression of 102 genes
  of \emph{E. Coli} after induction of the expression of SOD
  (recombinant human superoxide dismutase).
}

\usage{
data(ecoli)
}

\format{
  \code{caulobacter} is a \code{\link[longitudinal]{longitudinal}} object
  containing the data from the Schmidt-Heck et al. (2004) experiment.
  Essentially, this is a matrix with with 102 columns (=genes)
  and 9 rows (=time points).  All expression levels are given in
  log2-ratios with respect to the first time point (i.e. the
  induction at time 0).
}

\source{
  The micoarray experiment was performed at the Institute of Applied
  Microbiology, University of Agricultural Sciences of Vienne.
  The data and the experiment is described in Schmidt-Heck et al. (2004).
}


\references{
Schmidt-Heck, W., Guthke, R., Toepfer, S., Reischer, H., Duerrschmid, K.,
and Bayer, K. (2004) Reverse engineering of the stress response during 
expression of a recombinant protein.
In: \emph{Proceedings of the EUNITE 2004 European Symposium on Intelligent
Technologies, Hybrid Systems and their Implementation on Smart Adaptive 
Systems, June 10-12, 2004, Aachen, Germany}, Verlag Mainz, 
Wissenschaftsverlag, Aachen, 2004, 407-441 (ISBN 3-86130-368-X).
}

\examples{
# load GeneNet library
library("GeneNet")

# load data set
data(ecoli)
is.longitudinal(ecoli)

# how many samples and how many genes?
dim(ecoli)
summary(ecoli)
get.time.repeats(ecoli)

# plot first nine time series
plot(ecoli, 1:9)
}
\keyword{datasets}

