\name{ggm.make.graph}
\alias{ggm.make.graph}
\alias{ggm.make.dot}
\alias{network.make.dot}
\alias{show.edge.weights}

\title{Graphical Gaussian Models: Plotting the Network}

\usage{
ggm.make.dot(filename, edge.list, node.labels, main=NULL, show.edge.labels=FALSE)
network.make.dot(filename, edge.list, node.labels, main=NULL, show.edge.labels=FALSE)
ggm.make.graph(edge.list, node.labels, drop.singles=FALSE)
show.edge.weights(gr)
}

\arguments{
  \item{filename}{name of file containg the "dot" commands for graphviz}
  
  \item{edge.list}{a data frame, as obtained by \code{\link{ggm.test.edges}}, listing all
          edges to be included in the graph}
  \item{node.labels}{a vector with labels for each node (will be converted to type character)} 
  \item{main}{title included in plot}
  \item{show.edge.labels}{plot correlation values as edge labels (default: FALSE)}
  \item{drop.singles}{remove unconnected nodes}
  \item{gr}{a graph object}
  \item{\dots}{options passed to \code{plot} functions}
}
\description{
  \code{ggm.make.dot} converts an edge list as obtained by \code{\link{ggm.test.edges}}
  into a "dot" file that can directly be used for plotting the network with graphviz.

 \code{network.make.dot} is just an alias to \code{ggm.make.dot}.

  \code{ggm.make.graph} converts an edge list as obtained by \code{\link{ggm.test.edges}}
  into a graph object.

  \code{show.edge.weights} summarizes a graph object by prints a vector of weights for all
  edges contained in a graph.
  This function is convenient to gain a first impression of the graph (in particular if
  the "Rgraphviz" library is not installed). 

}

\details{
  For network plotting the software "graphviz" is employed (\url{http://www.graphviz.org}).

  For the functions \code{ggm.plot.graph} and \code{ggm.make.graph}  the "graph" and "Rgraphviz" 
  infrastructure from the  Bioconductor project (\url{http://www.bioconductor.org}) is required. 
 
}
\value{
  \code{ggm.make.dot} produces a "dot" network description file that
  can directly be fed into graphviz in order to produce a plot of a network.

  \code{ggm.make.graph} returns a graph object, suitable for plotting with functions from
  the "Rgraphviz" library.
  
  \code{show.edge.weights} returns a vector of weights for all edges contained in a graph.
  
}
\author{
  Juliane Schaefer (\url{http://www.stat.math.ethz.ch/~schaefer/}),
  Rainer Opgen-Rhein (\url{http://opgen-rhein.de}) and
  Korbinian Strimmer (\url{http://strimmerlab.org}).
}
\seealso{
\code{\link{ggm.test.edges}}, \code{plot.graph}.
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
# gr <- ggm.make.graph( test.results, nlab) 
# gr 
# show.edge.weights(gr) 
# gr2 <- ggm.make.graph( test.results, nlab, drop.singles=TRUE) 
# gr2 

# plot network
# NOTE: this requires the installation of the "Rgraphviz" library
# library("Rgraphviz")
# plot(gr)
}
\keyword{hplot}
