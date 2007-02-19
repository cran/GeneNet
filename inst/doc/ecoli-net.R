#######################################################################
# Note that this note can directly be run in R.
# Version: GeneNet 1.1.0 (February 2007)
#######################################################################


# This reproduces the "Ecoli" network example is from:

# Schaefer, J., and K. Strimmer, K. 2005c.  A shrinkage approach
# to large-scale covariance estimation and implications for
# functional genomics.   Statist. Appl.  Genet. Mol. Biol. 4: 32
# http://www.bepress.com/sagmb/vol4/iss1/art32/



# load GeneNet library
library("GeneNet")

# Example E. Coli data set (102 genes, 9 time points)
data(ecoli)
dim(ecoli)


##### INFER GGM NETWORK #####


###
### Step 1: Estimate partial correlation matrix
###

# there are many options for estimation the partial correlations
# - we recommend to use a shrinkage estimator 

inferred.pcor <- ggm.estimate.pcor(ecoli)
dim(inferred.pcor)


###
### Step 2: Assign p-values, q-values, empirical posterior probabilites to each edge
###

# parameters of the mixture distribution used to compute p-values etc.
c <- fdrtool(sm2vec(inferred.pcor), statistic="correlation")
c$param


test.results <- ggm.test.edges(inferred.pcor)

# show first 20 edges
test.results[1:20,]



###
### Step 3: Decide which edges to include in the network
###

# variant 1:
# how many edges are significant based on FDR cutoff Q=0.05 ?
significant1.idx <- test.results$qval <= 0.05
num.significant.1 <- sum(significant1.idx)
test.results[significant1.idx,]  # list significant edges


# variant 2:
# how many edges are significant based on "local fdr" 0.2 cutoff (prob > 0.80) ?
significant2.idx <- test.results$prob > 0.80
num.significant.2 <- sum(significant2.idx)
test.results[significant2.idx,] # list significant edges


##### PLOT GGM NETWORK #####

#########################################################
## variant 1: produce "dot" file for use with "graphviz"
##            (see http://www.graphviz.org/ )
#########################################################

node.labels <- colnames(ecoli)

# produce dot file (without edge labels)
ggm.make.dot(filename="ecoli.dot", 
             test.results[significant2.idx,], node.labels,
             main="Ecoli Network")

# call graphviz to produce a nice graph
system("fdp -T svg -o ecoli.svg ecoli.dot") # SVG format
system("fdp -T ps2 -o ecoli.eps ecoli.dot") # EPS format



# produce dot file (with edge labels)
ggm.make.dot(filename="ecoli.dot", 
             test.results[significant2.idx,], node.labels,
             show.edge.labels=TRUE, main="Ecoli Network")

# call graphviz to produce a nice graph
system("fdp -T svg -o ecoli.svg ecoli.dot") # SVG format
system("fdp -T ps2 -o ecoli.eps ecoli.dot") # EPS format



#########################################################
## variant 2: use "Rgrapviz" and related packages 
##            from http://www.bioconductor.org
#########################################################



# generate graph object with all significant edges
node.labels <- colnames(ecoli)
gr <- ggm.make.graph( test.results[significant2.idx,], 
         node.labels, drop.singles=TRUE) 
gr 

# print vector of edge weights
show.edge.weights(gr)


# plot network
#pdf(width=12, height=9, file="ecoli1.pdf")
ggm.plot.graph(gr,  main="Ecoli Network", show.edge.labels=FALSE)
#dev.off()


#pdf(width=12, height=9, file="ecoli2.pdf")
ggm.plot.graph(gr,  main="Ecoli Network", show.edge.labels=TRUE)
#dev.off()

