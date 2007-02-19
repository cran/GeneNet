######################################################################
# Note that this note can directly be run in R.
# Version: GeneNet 1.1.0 (February 2007)
#######################################################################


# Reconstruction of gene association network for 800 periodic
# genes from Arabidopsis thaliana

# This 800 gene data was used also in:
# Opgen-Rhein and Strimmer. 2007.  Learning causal networks from 
# systems biology time course data: an effective model selection 
# procedure for the vector autoregressive process. 
# BMC Bioinformatics 8 (suppl.) in press 


# Source: 
# Smith et al. 2004. Diurnal changes in the transcriptom encoding 
# enzymes of starch metabolism provide evidence for both transcriptional
# and posttranscriptional regulation of starch metabolism in Arabidopsis 
# leaves.  Plant Physiol. 136: 2687-2699


# This example was suggested by:
# Papapit Ingkasuwan, Division of Biotechnology,
# School of Bioresources and Technology,
# King Mongkut's University of Technology Thonburi, Bangkok, Thailand




########################
# step 1: inspect data #
########################



library("GeneNet")
data("arth800")
summary(arth800.expr)

# plot time series
plot(arth800.expr, 1:9)

# inspect pairwise scatter plots
panel.cor <- function(x, y, digits=2, prefix="", cex.cor)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y))
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex * r)
}
pairs(arth800.expr[,1:9], lower.panel=panel.smooth, upper.panel=panel.cor)



########################################
# step 2: compute partial correlations #
########################################

# the difference between "static" and "dynamic" is that the latter
# takes into account the spacings between the time points
 
pcor.stat <- ggm.estimate.pcor(arth800.expr, method = "static")
pcor.dyn <- ggm.estimate.pcor(arth800.expr, method = "dynamic")

# put partial correlation in a table (in order but without testing)
results.stat <- ggm.list.edges(pcor.stat)
results.dyn <- ggm.list.edges(pcor.dyn)

# put partial correlation in a table (with testing)
#results.stat <- ggm.test.edges(pcor.stat)
#sum(results.stat$prob > 0.80) # 7163 significant
#sum(results.stat$qval < 0.05) # 5559 significant

#results.dyn <- ggm.test.edges(pcor.dyn)
#sum(results.dyn$prob > 0.80) # 6102 significant
#sum(results.dyn$qval < 0.05) # 4583 significant


########################################################
# step 3: construct graph containing the 150 top edges #
########################################################

### plot network using graphviz

node.labels <- as.character(1:ncol(arth800.expr))
ggm.make.dot(filename="arthstat.dot", results.stat[1:150,], node.labels, main="Arabidopsis Network (Static)")
ggm.make.dot(filename="arthdyn.dot", results.dyn[1:150,], node.labels, main="Arabdiopsis Network (Dynamic)")

system("fdp -T svg -o arthstat.svg arthstat.dot") # SVG format
system("fdp -T svg -o arthdyn.svg arthdyn.dot") # SVG format



### plot network using Rgraphviz

node.labels <- as.character(1:ncol(arth800.expr))
gr.stat <- ggm.make.graph( results.stat[1:150,], node.labels, drop.singles=TRUE) 
gr.stat 
#A graphNEL graph with undirected edges
#Number of Nodes = 109
#Number of Edges = 150

gr.dyn <- ggm.make.graph( results.dyn[1:150,], node.labels, drop.singles=TRUE) 
gr.dyn 
#Number of Nodes = 107
#Number of Edges = 150


pdf(file="static.pdf", height=30, width=30)
ggm.plot.graph(gr.stat,  main="Arabidopis Network (static)", show.edge.labels=FALSE)
dev.off()

pdf(file="dynamic.pdf", height=30, width=30)
ggm.plot.graph(gr.dyn,  main="Arabidopsis Network (dynamic)", show.edge.labels=FALSE)
dev.off()



##########################
# step 4: interpretation #
##########################

# some of the discovered hubs
sort(degree(gr.stat), decreasing=TRUE)[1:10]
#570  81 558 539 726  47 198 422 783 793
# 22  13  12   9   9   8   8   8   8   8
sort(degree(gr.dyn), decreasing=TRUE)[1:10]
#570  81 783  47 422 558 452 539 738 272
# 20  17  10   9   9   9   8   8   8   7

arth800.descr[570]
#[1] "AP2 transcription factor - like protein"

arth800.descr[81]
#[1] "ATRPAC43; DNA binding / DNA-directed RNA polymerase; 

arth800.descr[558]
#[1]"structural constituent of ribosome;

arth800.descr[539]
#[1] "DNA binding / transcription factor; 

arth800.descr[783]
#[1] "RNA binding / RNA methyltransferase;


# look at the complete graph (requires testing in step 2 above)
gr.all.dyn <- ggm.make.graph( results.dyn[results.dyn$prob > 0.80,], 
  node.labels, drop.singles=TRUE)
gr.all.dyn
#A graphNEL graph with undirected edges
#Number of Nodes = 669 
#Number of Edges = 6102 

sort(degree(gr.all.dyn), decreasing=TRUE)[1:10]
#570  81 558 539 627 783 198 738  20 111 
#178 172 138 134 115 114 112  94  92  92 

