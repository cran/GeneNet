### ggm.make.igraph  (2013-06-16)
###
###   Construct "igraph" object from given network
###
### Copyright 2013 Korbinian Strimmer
###
###
### This file is part of the `GeneNet' library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 2, or at your option, any later version,
### incorporated herein by reference.
### 
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
### 
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA



# requires the installation of the "igraph" library

network.make.igraph = function(edge.list, node.labels, show.edge.labels = FALSE)
{
  return( ggm.make.igraph(edge.list=edge.list, node.labels=node.labels,
                          show.edge.labels=show.edge.labels) 
  )
}

ggm.make.igraph = function(edge.list, node.labels, show.edge.labels = FALSE)
{
  require("igraph")

  V <- unique(node.labels)
  if (length(V) != length(node.labels)) {
    stop("Duplicate node labels encountered. Node labels must be unique!")
  }
  V <- as.character(V)  

  w = edge.list[,1] # partial correlations
  cutoff = quantile(abs(w), c(0.2, 0.8))
  num.edges = length(w)

  if (ncol(edge.list) == 11) { # with specified direction
    edge.type = edge.list[,11]

    # create new nodes 
    node1 = vector(mode="character", length=num.edges)
    node2 = vector(mode="character", length=num.edges)
    am = integer(num.edges)

    for (i in 1:num.edges) {
      eti = edge.type[i]
      if (eti=="1to2") {
        node1[i] = V[edge.list[i,2]]
        node2[i] = V[edge.list[i,3]]
        am[i] = 2
      }
      if (eti=="2to1") {
        node1[i] = V[edge.list[i,3]]
        node2[i] = V[edge.list[i,2]]
        am[i] = 2
      }
      if (eti=="undirected") { 
        node1[i] = V[edge.list[i,2]]
        node2[i] = V[edge.list[i,3]]
        am[i] = 0 # don't show arrow
      }
    }
  }
  else { # undirected graph
    node1 = V[edge.list[,2]]
    node2 = V[edge.list[,3]]
    am = rep(0, num.edges) # don't show arrows
  }

  ig = graph.edgelist( cbind( node1, node2 ), directed=TRUE)
  E(ig)$arrow.mode=am

  # assign further edge properties
  E(ig)$width = ifelse( abs(w) < cutoff[2], 1.5, 3) 
  E(ig)$color = ifelse( abs(w) < cutoff[1], "grey", "black") 
  E(ig)$lty =   ifelse( w < 0, 3, 1) 
  if(show.edge.labels) {
    E(ig)$label = as.character(round(w, digits = 4))
    E(ig)$label.color = 1
    E(ig)$label.cex = 0.7 
  }

  # b&w vertices
  V(ig)$color=NA
  V(ig)$label.color=1

  return(ig)
}

