\name{ggm.plot.graph}
\alias{ggm.plot.graph}
\alias{ggm.make.graph}
\alias{show.edge.weights}

\title{Graphical Gaussian Models: Plotting the Network}

\usage{
ggm.make.graph(edge.list, node.labels, drop.singles=FALSE)
show.edge.weights(gr)
ggm.plot.graph(gr, 
    layoutType=c("fdp", "neato", "circo", "dot", "twopi"), 
    show.edge.labels=FALSE, \dots)
}

\arguments{
  \item{edge.list}{a data frame, as obtained by \code{\link{ggm.test.edges}}, listing all
          edges to be included in the graph}
  \item{node.labels}{a vector with labels for each node (will be converted to type character)} 
  \item{drop.singles}{remove unconnected nodes}
  \item{gr}{a graph object}
  \item{layoutType}{type of plot of the graph. Defaults to "fdp".}
  \item{show.edge.labels}{plot correlation values as edge labels (default: FALSE)}
  \item{\dots}{options passed to \code{plot} functions}
}
\description{
  \code{ggm.make.graph} converts an edge list as obtained by \code{\link{ggm.test.edges}}
  into a graph object.

  \code{show.edge.weights} summarizes a graph object by prints a vector of weights for all
  edges contained in a graph.
  This function is convenient to gain a first impression of the graph (in particular if
  the "Rgraphviz" library is not installed). 

  \code{ggm.plot.graph} visualizes the network structure of the graphical Gaussian model 
  using the Rgraphviz network plot package. The correlation coefficients can be 
  printed as edge labels.
}

\details{
  The network plotting functions require the installation of the "graph" and "Rgraphviz" R
  packages. These are available from the Bioconductor website
  (\url{http://www.bioconductor.org}).  Note that it is not necessary to install the
  complete set of Bioconductor packages, only "graph" and "Rgraphviz" are needed by the
  GeneNet package (however, these may in turn require additional packages from Bioconductor).

  \code{ggm.plot.graph} is a utility function to plot a graph in "dfp" layout format.
  The line width and color depends on the relative strength of the partial correlation 
  assigned to an edge: the top 80%-100% edges are drawn in "black" with "lwd=2", 
  the range from 20% to 80% will be drawn in "black" with "lwd=1", and the lower
  20% will be plotted in "grey" with "lwd=1".
  
  The "Rgraphviz" package offers many other options - please consult the manual of
  "Rgraphviz" for the details.

 
}
\value{
  \code{ggm.make.graph} returns a graph object, suitable for plotting with functions from
  the "Rgraphviz" library.
  
  \code{show.edge.weights} returns a vector of weights for all edges contained in a graph.
  
  \code{ggm.plot.graph} plots  the network on the current graphic device.
}
\author{
  Juliane Schaefer (\url{http://www.stat.math.ethz.ch/~schaefer/}) and
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

# generate graph object 
# NOTE: this requires the installation of the "graph" package
# (in the following "try" is used to avoid an error if the library is not installed)
nlab <- LETTERS[1:20]
try( gr <- ggm.make.graph( test.results, nlab) )
try( gr )
try( show.edge.weights(gr) )
try( gr2 <- ggm.make.graph( test.results, nlab, drop.singles=TRUE) )
try( gr2 )

# plot network
# NOTE: this requires the installation of the "Rgraphviz" library
try ( ggm.plot.graph(gr, main = "A graph"))
try ( ggm.plot.graph(gr2, main = "The same graph with singles removed" ))
}
\keyword{hplot}
