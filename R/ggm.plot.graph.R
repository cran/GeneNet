### ggm.plot.graph  (2006-08-08)
###
###   Plotting the GGM network
###
### Copyright 2003-06 Juliane Schaefer and Korbinian Strimmer
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
  
    em <- edgeMatrix(gr, duplicates=FALSE)
      
    return( eWV(gr, em, useNNames = TRUE) )       
   
}



# requires installation of the "Rgraphviz" library


ggm.plot.graph <- function(gr, 
    layoutType=c("fdp", "neato", "circo", "dot", "twopi"), 
    show.edge.labels=FALSE,  ...)
{
    layoutType = match.arg(layoutType)
    
    require(Rgraphviz) # requires the "Rgraphviz" package from BioC >= 1.5
  
    # general graph attributes
    gAttrs <- getDefaultAttrs(layoutType="neato")
    gAttrs$edge$color <- "black"    
    gAttrs$node$shape <- "ellipse"
    gAttrs$node$fixedsize <- FALSE
    
    
    #  edge attributes
    em <- edgeMatrix(gr)
    emv <- eWV(gr, em, sep="~", useNNames = TRUE)
    edge.names <- names(emv)
    edge.labels <- as.character(emv)
    eAttrs <- list()
    eAttrs$label <- edge.labels
    names(eAttrs$label) <- edge.names
    
                
    # get Ragraph object (=perform graph layout)
    lg <- agopen(gr, name="test", attrs=gAttrs, edgeAttrs = eAttrs, layoutType=layoutType)
    
    
    # modify Ragraph object to allow different line types
    
    # thresholds for line width and coloring
    cutoff <- quantile(abs(emv), c(0.2, 0.8)) 

    eg <- AgEdge(lg)
    # loop through all edges   
    for (i in 1:length(eg))
    {
      # modify parameters  (it is important to set params for all edges)
      
      w <- as.double( eg[[i]]@txtLabel@labelText )
      
      if (w < 0 )
        eg[[i]]@lty <- 3  # negative values: dotted lines
      else
        eg[[i]]@lty <- 1  # positive values: solid lines

      
      # line thickness and color depends on relative strengh	
      if (abs(w) < cutoff[1]) # lower 20% quantile
      {
         eg[[i]]@lwd <- 1
	 eg[[i]]@color <- "grey"
      }
      else if (abs(w) < cutoff[2]) # from 20%-80%
      {
         eg[[i]]@lwd <- 1
	 eg[[i]]@color <- "black"
      }
      else # top 80%-100%
      {
         eg[[i]]@lwd <- 2
         eg[[i]]@color <- "black"
      }
	 
      # remove edge weight
      if (!show.edge.labels)
        eg[[i]]@txtLabel@labelText <- "" 
	 
	 
    }
    AgEdge(lg) <- eg
   
   
    # finally, plot Ragraph object
   
    plot(lg, ...)
}

