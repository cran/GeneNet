### ggm.make.graph  (2007-05-24)
###
###   Construct "graph" object from given network
###
### Copyright 2003-07 Juliane Schaefer and Korbinian Strimmer
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



# requires the installation of the "graph" library

# generate a graph object from and edge list
# (such as obtained from ggm.test.edges) 
ggm.make.graph <- function(edge.list, node.labels, drop.singles=FALSE)
{
    require(graph) # requires the "graph" package from BioC >= 1.5
      
    V <- unique(node.labels)  
    if ( length(V) != length(node.labels) )
    {
       stop("Duplicate node labels encountered. Node labels must be unique!")
    }   
    V <- as.character(V)
     
    
    # create empty graph with no edges
    gR <- new("graphNEL", nodes=V)
   
    # add edges and edge weights (correlations)
    gX <- addEdge(V[edge.list[,2]],
                  V[edge.list[,3]],
                  gR,
                  round(edge.list[,1], digits=2) )
 
 
    if(drop.singles) # remove unconnected nodes
    {
      # nodes with degree > 0
      nd <- nodes(gX)[ degree(gX) > 0 ]
    
      gX <- subGraph(nd, gX)
    }
  
    return(gX)
}



# print vector of edge weights
show.edge.weights <- function(gr)
{
    require(graph) # requires the "graph" package from BioC >= 1.5
  
    em <- graph::edgeMatrix(gr, duplicates=FALSE) # to avoid name clash
      
    return( eWV(gr, em, useNNames = TRUE) )       
   
}


