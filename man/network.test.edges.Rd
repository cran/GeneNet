\name{network.test.edges}
\alias{network.test.edges}
\alias{extract.network}

\title{Graphical Gaussian Models: Assess Significance of Edges (and Directions)}

\usage{
network.test.edges(r.mat, fdr=TRUE, direct=FALSE, plot=TRUE, ...)
extract.network(network.all, method.ggm=c("prob", "qval","number"), 
      cutoff.ggm=0.8, method.dir=c("prob","qval","number", "all"), 
      cutoff.dir=0.8, verbose=TRUE)
}

\arguments{
  \item{r.mat}{matrix of partial correlations}
  \item{fdr}{estimate q-values and local fdr}
  \item{direct}{compute additional statistics for obtaining a partially directed network}
  \item{plot}{plot density and distribution function and (local) fdr values}
  \item{\dots}{parameters passed on to \code{\link[fdrtool]{fdrtool}} }
  \item{network.all}{list with partial correlations and fdr values for all potential edges
                     (i.e. the output of \code{network.test.edges}}
  \item{method.ggm}{determines which criterion is used to select significant partial correlations (default: prob)}
  \item{cutoff.ggm}{default cutoff for significant partial correlations}
  \item{method.dir}{determines which criterion is used to select significant directions (default: prob) }
  \item{cutoff.dir}{default cutoff for significant directions}
  \item{verbose}{print information on the number of significant edges etc.}
}
\description{
 \code{network.test.edges} returns a data frame containing all edges listed
  in order of the magnitude of the partial correlation associated with each edge.
  If \code{fdr=TRUE} then in addition the
  p-values, q-values and posterior probabilities (=1 - local fdr) for each potential 
  edge are computed.

 \code{extract.network} returns a data frame with a subset of significant
 edges. 
}

\details{
  For assessing the significance of edges in the GGM 
  a mixture model is fitted to the partial correlations using \code{\link[fdrtool]{fdrtool}}. 
  This results in (i)  two-sided p-values for the test of non-zero correlation,
  (ii) corresponding posterior probabilities (= 1- local fdr), as well as (iii) tail
  area-based q-values.  See  Sch\"afer and Strimmer (2005) for details.

  For determining putatative directions on this GGM log-ratios of standardized partial 
  variances re estimated, and subsequently the corresponding (local) fdr values are
  computed - see Opgen-Rhein and Strimmer (2007).
}
\value{
  \code{network.test.edges} returns a
    data frame with the following columns:
  
  \item{pcor}{correlation (from r.mat)}
  \item{node1}{first node connected to edge}
  \item{node2}{second node connected to edge}
  \item{pval}{p-value}
  \item{qval}{q-value}
  \item{prob}{probability that edge is nonzero (= 1-local fdr}
  \item{log.spvar}{log ratio of standardized partial variance (determines direction)}
  \item{pval.dir}{p-value (directions)}
  \item{qval.dir}{q-value (directions)}
  \item{prob.dir}{1-local fdr (directions) }


  Each row in the data frame corresponds to one edge, and the rows are sorted
  according the absolute strength of the correlation (from strongest to weakest)
  
  \code{extract.network} processes the above data frame containing all potential edges,
  and returns a dataframe with a subset of edges.  If applicable, an additional
  last column (11) contains additional information on the directionality of an edge.
}
\author{
  Rainer Opgen-Rhein,
  Juliane Sch\"afer,
  Korbinian Strimmer (\url{http://strimmerlab.org}).
}
\references{
  Sch\"afer, J., and Strimmer, K. (2005).  An empirical Bayes approach to inferring
  large-scale gene association networks. \emph{Bioinformatics} \bold{21}:754-764.

  Opgen-Rhein, R., and K. Strimmer. (2007). From correlation to causation networks: 
  a simple approximate learning algorithm and its application to high-dimensional 
  plant gene expression data. \emph{BMC Syst. Biol.} \bold{1}:37.

}
\seealso{
\code{\link{cor0.test}},
\code{\link[fdrtool]{fdrtool}},
\code{\link{ggm.estimate.pcor}}.
}

\examples{
# load GeneNet library
library("GeneNet")
 
# ecoli data 
data(ecoli)

# estimate partial correlation matrix 
inferred.pcor <- ggm.estimate.pcor(ecoli)

# p-values, q-values and posterior probabilities for each potential edge 
#
test.results <- network.test.edges(inferred.pcor)

# show best 20 edges (strongest correlation)
test.results[1:20,]

# extract network containing edges with prob > 0.9 (i.e. local fdr < 0.1)
net <- extract.network(test.results, cutoff.ggm=0.9)
net

# how many are significant based on FDR cutoff Q=0.05 ?
num.significant.1 <- sum(test.results$qval <= 0.05)
test.results[1:num.significant.1,]

# how many are significant based on "local fdr" cutoff (prob > 0.9) ?
num.significant.2 <- sum(test.results$prob > 0.9)
test.results[test.results$prob > 0.9,]

# parameters of the mixture distribution used to compute p-values etc.
c <- fdrtool(sm2vec(inferred.pcor), statistic="correlation")
c$param

}
\keyword{htest}
