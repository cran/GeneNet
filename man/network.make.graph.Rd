\name{network.make.graph}
\alias{network.make.graph}
\alias{network.make.dot}
\alias{node.degree}
\alias{edge.info}
\alias{num.nodes}


\title{Graphical Gaussian Models: Plotting the Network}

\usage{
network.make.dot(filename, edge.list, node.labels, main=NULL, show.edge.labels=FALSE)
network.make.graph(edge.list, node.labels, drop.singles=FALSE)
edge.info(gr)
node.degree(gr)
num.nodes(gr)
}

\arguments{
  \item{filename}{name of file containg the "dot" commands for graphviz}
  
  \item{edge.list}{a data frame, as obtained by \code{\link{network.test.edges}}, listing all
          edges to be included in the graph}
  \item{node.labels}{a vector with labels for each node (will be converted to type character)} 
  \item{main}{title included in plot}
  \item{show.edge.labels}{plot correlation values as edge labels (default: FALSE)}
  \item{drop.singles}{remove unconnected nodes}
  \item{gr}{a graph object}
  \item{\dots}{options passed to \code{plot} functions}
}
\description{
  \code{network.make.dot} converts an edge list as obtained by \code{\link{network.test.edges}}
  into a "dot" file that can directly be used for plotting the network with graphviz.

  \code{network.make.graph} converts an edge list as obtained by \code{\link{network.test.edges}}
  into a graph object. 

  \code{edge.info} shows the edge weights and the edge directions. 

   \code{node.degree} shows number of edges connected to a node (bi-directional/undirected edges are counted only once).

  \code{num.nodes} shows the number of nodes.
}

\details{
  For network plotting the software "graphviz" is employed (\url{http://www.graphviz.org}).

  For the functions \code{network.plot.graph} and \code{network.make.graph}  the "graph" and "Rgraphviz" 
  packages from the  Bioconductor project (\url{http://www.bioconductor.org}) is required. 
 
}
\value{
  \code{network.make.dot} produces a "dot" network description file that
  can directly be fed into graphviz in order to produce a plot of a network.

  \code{network.make.graph} returns a graph object, suitable for plotting with functions from
  the "Rgraphviz" library.
  
  \code{edge.info} returns a list containing vector of weights for all edges contained in a graph, and a vector listing the directions of the edges (using Rgraphviz
   conventions "forward" for directed edge, and "none" for bi-directional/undirected edge). 

  \code{num.nodes} returns the number of nodes.
 
}
\author{
  Juliane Sch\"afer,
  Rainer Opgen-Rhein, and
  Korbinian Strimmer (\url{http://strimmerlab.org}).
}
\seealso{
\code{\link{network.test.edges}}, \code{plot.graph}.
}

\examples{
# load GeneNet library
library("GeneNet")
 
# generate random network with 20 nodes and 10 percent edges (=19 edges)
true.pcor <- ggm.simulate.pcor(20, 0.1)

# convert to edge list 
test.results <- ggm.list.edges(true.pcor)[1:19,]

########  use graphviz directly to produce a plot ##########

# uncomment for actual use!

# nlab <- LETTERS[1:20]
# ggm.make.dot(filename="test.dot", test.results, nlab, main = "A graph") 
# system("fdp -T svg -o test.svg test.dot") # SVG format


########  use Rgraphviz produce a plot ##########

# uncomment for actual use!

# nlab <- LETTERS[1:20]
# gr <- network.make.graph( test.results, nlab) 
# gr 
# num.nodes(gr)
# edge.info(gr)
# gr2 <- network.make.graph( test.results, nlab, drop.singles=TRUE) 
# gr2 
# num.nodes(gr2)
# edge.info(gr2)

# plot network
# NOTE: this requires the installation of the "Rgraphviz" library
# library("Rgraphviz")
# plot(gr, "fdp")
# plot(gr2, "fdp")

## for a full example with beautified Rgraphviz plot see 
## the example scripts provide with GeneNet (e.g. arabidopis-net.R)


}
\keyword{hplot}
