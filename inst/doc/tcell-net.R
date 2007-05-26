######################################################################
# This note can be directly run in R.
# Requires GeneNet 1.2.0 (May 2007)
#######################################################################



# This reproduces the four networks displayed in 
# Opgen-Rhein and Strimmer (2006a,b)

Opgen-Rhein, R., and K. Strimmer. 2006a.  Using regularized dynamic
  correlation to infer gene dependency networks from time-series 
  microarray data.  Proceedings of WCSB 2006 (June 12-13, 2006, 
  Tampere, Finland).
Opgen-Rhein, R., and K. Strimmer. 2006b.  Inferring gene dependency 
  networks from genomic longitudinal data: a functional data approach. 
  REVSTAT 4:53-65.  



#load GeneNet library
library("GeneNet")


# get T cell data 
data(tcell)
tc44 <- combine.longitudinal(tcell.10, tcell.34)


# estimate partial correlations
pc1 <- ggm.estimate.pcor(tc44, lambda=0)                   # static, no shrinkage
pc2 <- ggm.estimate.pcor(tc44, method="dynamic", lambda=0) # dynamic, no shrinkage
pc3 <- ggm.estimate.pcor(tc44)                             # static, with shrinkage         
pc4 <- ggm.estimate.pcor(tc44, method="dynamic")           # dynamic, with shrinkage


# significant edges (local fdr, 0.2 cutoff)

# static, no shrinkage
t1.edges <- ggm.test.edges(pc1)
t1.net <- extract.network(t1.edges) # prob > 0.8
t1.net

# dynamic, no shrinkage
t2.edges <- ggm.test.edges(pc2)
t2.net <- extract.network(t2.edges) # prob > 0.8
t2.net

# static, with shrinkage
t3.edges <- ggm.test.edges(pc3)
t3.net <- extract.network(t3.edges) # prob > 0.8
t3.net

# dynamic, with shrinkage
t4.edges <- ggm.test.edges(pc4)
t4.net <- extract.network(t4.edges) # prob > 0.8
t4.net


######## produce plots using graphviz ###########


# produce dot file (without edge labels)
node.labels <- colnames(tc44)
ggm.make.dot(filename="net1.dot", t1.net, node.labels, main="Static")
ggm.make.dot(filename="net2.dot", t2.net, node.labels, main="Dynamic")
ggm.make.dot(filename="net3.dot", t3.net, node.labels, main="Static + Shrink")
ggm.make.dot(filename="net4.dot", t4.net, node.labels, main="Dynamic + Shrink")

# call graphviz to produce a nice graph
system("fdp -T svg -o net1.svg net1.dot") # SVG format
system("fdp -T svg -o net2.svg net2.dot") # SVG format
system("fdp -T svg -o net3.svg net3.dot") # SVG format
system("fdp -T svg -o net4.svg net4.dot") # SVG format



######## produce plots using Rgraphviz ###########

node.labels <- colnames(tc44)
gr1 <- ggm.make.graph( t1.net, node.labels, drop.singles=TRUE) 
gr2 <- ggm.make.graph( t2.net, node.labels, drop.singles=TRUE)  
gr3 <- ggm.make.graph( t3.net, node.labels, drop.singles=TRUE) 
gr4 <- ggm.make.graph( t4.net, node.labels, drop.singles=TRUE)  

gr1
gr2
gr3
gr4


# plot networks
library("Rgraphviz")
par(mfrow=c(2,2))
plot(gr1, main="Static, no shrinkage")
plot(gr2, main="Dynamic, no shrinkage")
plot(gr3, main="Static, with shrinkage")
plot(gr4, main="Dynamic, with shrinkage")
par(mfrow=c(1,1))

