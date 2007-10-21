######################################################################
# This note can be directly run in R.
# Requires GeneNet 1.2.1 (October 2007)
#######################################################################


# Reconstruction of partially directed gene association network for 
# 800 periodic genes from Arabidopsis thaliana from:

# Opgen-Rhein, R., and K. Strimmer. 2007. From correlation to causation 
# networks: a simple approximate learning algorithm and its application 
# to high-dimensional plant gene expression data.
# BMC Syst. Biol. 1: 37.
 

# Original source of the data: 
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

pcor.dyn <- ggm.estimate.pcor(arth800.expr, method = "dynamic")


###########################################################
# step 3: assign (local) fdr values to all possible edges #
###########################################################

arth.edges <- network.test.edges(pcor.dyn,direct=TRUE)

dim(arth.edges)


########################################################
# step 4: construct graph containing the 150 top edges #
########################################################

arth.net <- extract.network(arth.edges, method.ggm="number", cutoff.ggm=150)


#####################################
# step 5: plot graph using graphviz #
#####################################

### plot network using graphviz (recommended way)

node.labels <- as.character(1:ncol(arth800.expr))
network.make.dot(filename="arthdyn.dot", arth.net, node.labels, main="Arabdiopsis Network")
system("fdp -T svg -o arthdyn.svg arthdyn.dot") # SVG format


### plot network using Rgraphviz

library("graph")

gr <- ggm.make.graph( arth.net, node.labels, drop.singles=TRUE) 
gr 
#Number of Nodes = 107
#Number of Edges = 150

library("Rgraphviz")
plot(gr, "neato")



############################
# step 6: further analysis #
############################

# some of the discovered hubs
sort(degree(gr), decreasing=TRUE)[1:10]
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

