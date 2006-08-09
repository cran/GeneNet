\name{ggm.test.edges}
\alias{ggm.test.edges}
\alias{ggm.list.edges}

\title{Graphical Gaussian Models: Testing Edges}

\usage{
ggm.list.edges(r.mat)
ggm.test.edges(r.mat, df=7, plot.locfdr=1)
}

\arguments{
  \item{r.mat}{matrix of partial correlations}
  \item{df}{degrees of freedom for the spline fitting the density (only if fA.type="nonparametric")}
  \item{plot.locfdr}{controls plot option in \code{locfdr} }
}
\description{
 \code{ggm.list.edges} returns a table with all correlations sorted according to magnitude,
  as well as the two node numbers associated with each edge.

  \code{ggm.test.edges} in addition assigns statistical significance to the edges in a GGM network by computing
  p-values, q-values and posterior probabilities for each potential edge.
}

\details{
  A mixture model is fitted to the partial correlations using \code{\link{cor.fit.mixture}}
  Subsequently, two-sided p-values to test non-zero correlation are computed for each edge using 
  \code{\link{cor0.test}}. In addition, corresponding posterior probabilities are
  computed (also using \code{\link{cor.fit.mixture}}.  Finally, to simplify multiple testing q-values
  are computed via  \code{\link[fdrtool]{fdr.control}} whith the specified value of \code{eta0} taken
  into account.

  Theoretical details are explained in Schaefer and Strimmer (2005), along with a simulation
  study and an application to gene expression data.
}
\value{
  A sorted data frame with the following columns:
  
  \item{pcor}{correlation (from r.mat)}
  \item{node1}{first node connected to edge}
  \item{node2}{second node connected to edge}
  \item{pval}{p-value}
  \item{qval}{q-value}
  \item{prob}{probability that edge is nonzero (= 1-local fdr}

  Each row in the data frame corresponds to one edge, and the rows are sorted
  according the absolute strength of the correlation (from strongest to weakest)
  
  \code{ggm.list.edges} returns only the first three columns.
}
\author{
  Juliane Schaefer (\url{http://www.stat.math.ethz.ch/~schaefer/}) and
  Korbinian Strimmer (\url{http://strimmerlab.org}).
}
\references{
  Schaefer, J., and Strimmer, K. (2005).  An empirical Bayes approach to inferring
  large-scale gene association networks. \emph{Bioinformatics} \bold{21}:754-764.
}
\seealso{
\code{\link{cor.fit.mixture}},
\code{\link{cor0.test}},
\code{\link[fdrtool]{fdr.control}},
\code{\link{ggm.estimate.pcor}}.
}

\examples{
# load GeneNet library
library("GeneNet")
 
# ecoli data 
data(ecoli)

# estimate partial correlation matrix 
inferred.pcor <- ggm.estimate.pcor(ecoli)

# p-values, q-values and posterior probabilities for each edge 
#
test.results <- ggm.test.edges(inferred.pcor)

# show best 20 edges (strongest correlation)
test.results[1:20,]

# how many are significant based on FDR cutoff Q=0.05 ?
num.significant.1 <- sum(test.results$qval <= 0.05)
test.results[1:num.significant.1,]

# how many are significant based on "local fdr" cutoff (prob > 0.95) ?
num.significant.2 <- sum(test.results$prob > 0.95)
test.results[1:num.significant.2,]

# parameters of the mixture distribution used to compute p-values etc.
c <- cor.fit.mixture(sm2vec(inferred.pcor))
c$eta0
c$kappa
}
\keyword{htest}
