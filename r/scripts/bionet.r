library("Rgraphviz")
library("BioNet")

##load DIP ppi
ig <- loadNetwork.sif(network.file.name)
data <- read.delim(deg.file.name, sep="\t", stringsAsFactors=FALSE, row.names = 1)
pval <- data[["pval"]]
names(pval) <- rownames(data)
combined.dataset <- data
subnet <- subNetwork(rownames(combined.dataset), ig)
subnet <- rmSelfLoops(subnet)

pval <-na.omit(pval)
# pval <- pval[!pval==1]
fb <- fitBumModel(pval, plot = TRUE, starts = 10)
scores <- scoreNodes(subnet, fb, fdr = .05) #0.9999999999
# scores = -log10(pval)
module <- runFastHeinz(subnet, scores)
# logFC <- dataLym$diff
# names(logFC) <- dataLym$label
# plotModule(module, scores = scores) # , diff.expr = logFC
module.genes <- nodes(module)
bg.genes <- nodes(subnet)