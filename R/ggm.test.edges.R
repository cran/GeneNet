### ggm.test.edges  (2007-02-17)
###
###   Compute p-values, q-values and posterior probabilities for GGM edges
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


##############################

ggm.list.edges <- function(r.mat)
{
   pcor <- sm2vec(r.mat)
   indexes <- sm.index(r.mat)
   colnames(indexes) <- c("node1", "node2")

   result <- cbind(pcor, indexes)

   # sort according to magnitude of correlation 
   sort.idx <- order(-abs(result[,1]))
   result <- result[sort.idx,]
   
   # return as data frame
   return(as.data.frame(result))
}




# assign p-values, q-values and posterior probabilities to each edge
ggm.test.edges <- function(r.mat, ...)
{
   w <- ggm.list.edges(r.mat)
   pc <- w[,1]
    
   # estimate kappa and eta0
   fdr.out <- fdrtool(pc, statistic="correlation", ...)
   pval <- fdr.out$pval
   qval <- fdr.out$qval
   prob <- 1-fdr.out$lfdr
      
   result <- cbind(w, pval, qval, prob)

   
   # return as data frame
   return(as.data.frame(result))
}

