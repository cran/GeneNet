\name{ggm.make.igraph}
\alias{ggm.make.igraph}
\alias{network.make.igraph}
\alias{ggm.make.dot}
\alias{network.make.dot}

\title{Graphical Gaussian Models: Plotting the Network}

\usage{
ggm.make.igraph(edge.list, node.labels, show.edge.labels=FALSE)
network.make.igraph(edge.list, node.labels, show.edge.labels=FALSE)
ggm.make.dot(filename, edge.list, node.labels, main=NULL, show.edge.labels=FALSE)
network.make.dot(filename, edge.list, node.labels, main=NULL, show.edge.labels=FALSE)
}

\arguments{
  \item{filename}{name of file containg the "dot" commands for graphviz}
  
  \item{edge.list}{a data frame, as obtained by \code{\link{ggm.test.edges}}, listing all
          edges to be included in the graph}
  \item{node.labels}{a vector with labels for each node (will be converted to type character)} 
  \item{main}{title included in plot}
  \item{show.edge.labels}{show partial correlation as edge labels (default: FALSE)}
  \item{\dots}{options passed to \code{plot} functions}
}
\description{

  \code{ggm.make.igraph} converts an edge list as obtained by \code{\link{ggm.test.edges}}
  into an igraph object.

 \code{network.make.igraph} is just an alias to \code{network.make.igraph}.

  \code{ggm.make.dot} converts an edge list as obtained by \code{\link{ggm.test.edges}}
  into a "dot" file that can directly be used for plotting the network with graphviz.

 \code{network.make.dot} is just an alias to \code{ggm.make.dot}.

}

\details{
  The above functions provide two different ways for plotting networks computed with GeneNet.

\code{ggm.make.igraph} and \code{network.make.igraph}  use the R package \code{igraph} for plotting - see \url{hhttp://igraph.sourceforge.net/}.

\code{ggm.make.dot} and \code{network.make.dot} produce 'dot' files for use with the \code{graphviz} software - see \url{http://www.graphviz.org}.

In the resulting plots, dotted lines indicate negative partial correlation.  The strength of the partial correlation is visualized by the line width and the color of the edge: the strongest 20 percent of all edges are shown with thick black lines, whereas the weakest 20 percent are shown in thin grey lines.
}
\value{
  \code{ggm.make.dot} produces a "dot" network description file that
  can directly be fed into the 'graphviz' in order to produce a plot of a network.

  \code{ggm.make.igraph} returns a graph object, suitable for plotting with functions from
  the \code{igraph} R package.

  
}
\author{
  Korbinian Strimmer (\url{http://strimmerlab.org}).
}
\seealso{
\code{\link{ggm.test.edges}}, \code{plot.igraph}.
}

\examples{
# load GeneNet library
library("GeneNet")
 
# generate random network with 20 nodes and 10 percent edges (=19 edges)
true.pcor <- ggm.simulate.pcor(20, 0.1)

# convert to edge list 
edge.list <- ggm.list.edges(true.pcor)[1:19,]

nlab <- LETTERS[1:20] # node labels

########  use igraph R package produce a plot ##########

igr1 <- ggm.make.igraph(edge.list, nlab) 
plot(igr1, main = "A Random Graph")

igr2 <- ggm.make.igraph(edge.list, nlab, show.edge.labels=TRUE) 
plot(igr2, main = "A Random Graph with Partial Correlations")

# igraph allows to fine-tune the plot
# e.g. smaller edge labels and red nodes:
plot(igr2, main = "A Random Graph with Partial Correlations",
  edge.label.cex=0.7, vertex.color="red")


########  use graphviz to produce a plot ##########

# uncomment for actual use!

# nlab <- LETTERS[1:20]
# ggm.make.dot(filename="randomnet.dot", edge.list, nlab, main = "A graph") 
# system("fdp -T svg -o randomnet.svg randomnet.dot") # SVG format
# system("fdp -T jpg -o randomnet.jpg randomnet.dot") # JPG format

}
\keyword{hplot}
