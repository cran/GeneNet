#######################################################################
# Note that this note can directly be run in R.
# Version: GeneNet 1.1.0 (February 2007)
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
t1 <- ggm.test.edges(pc1)
s1.idx <- t1$prob > 0.80
num.s1 <- sum(s1.idx)
t1[s1.idx,] # 31 significant edges

# dynamic, no shrinkage
t2 <- ggm.test.edges(pc2)
s2.idx <- t2$prob > 0.80
num.s2 <- sum(s2.idx)
t1[s2.idx,] # 55 significant edges

# static, with shrinkage
t3 <- ggm.test.edges(pc3)
s3.idx <- t3$prob > 0.80
num.s3 <- sum(s3.idx)
t1[s3.idx,] # 10 significant edges

# dynamic, with shrinkage
t4 <- ggm.test.edges(pc4)
s4.idx <- t4$prob > 0.80
num.s4 <- sum(s4.idx)
t1[s4.idx,] # 49 significant edges


######## produce plots using graphviz ###########


# produce dot file (without edge labels)
node.labels <- colnames(tc44)
ggm.make.dot(filename="net1.dot", t1[s1.idx,], node.labels, main="Static")
ggm.make.dot(filename="net2.dot", t1[s1.idx,], node.labels, main="Dynamic")
ggm.make.dot(filename="net3.dot", t1[s1.idx,], node.labels, main="Static + Shrink")
ggm.make.dot(filename="net4.dot", t1[s1.idx,], node.labels, main="Dynamic + Shrink")

# call graphviz to produce a nice graph
system("fdp -T svg -o net1.svg net1.dot") # SVG format
system("fdp -T svg -o net2.svg net2.dot") # SVG format
system("fdp -T svg -o net3.svg net3.dot") # SVG format
system("fdp -T svg -o net4.svg net4.dot") # SVG format



######## produce plots using Rgraphviz ###########

node.labels <- colnames(tc44)
gr1 <- ggm.make.graph( t1[s1.idx,], node.labels, drop.singles=TRUE) 
gr2 <- ggm.make.graph( t2[s2.idx,], node.labels, drop.singles=TRUE)  
gr3 <- ggm.make.graph( t3[s3.idx,], node.labels, drop.singles=TRUE) 
gr4 <- ggm.make.graph( t4[s4.idx,], node.labels, drop.singles=TRUE)  

gr1
gr2
gr3
gr4


# plot networks
ggm.plot.graph(gr1, main="Static, no shrinkage", layoutType="neato")
ggm.plot.graph(gr2, main="Dynamic, no shrinkage", layoutType="neato")
ggm.plot.graph(gr3, main="Static, with shrinkage", layoutType="neato")
ggm.plot.graph(gr4, main="Dynamic, with shrinkage", layoutType="neato")


