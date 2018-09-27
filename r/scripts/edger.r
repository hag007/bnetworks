library(edgeR)

y <- DGEList(counts=data,group=group)
y <- calcNormFactors(y)
design <- model.matrix(~group)
y <- estimateDisp(y,design)
et <- exactTest(y)
edgeR_results = topTags(et, n = nrow(data))$table
# rownames(edgeR_results) <- genes
result <- data.frame(edgeR_results[order(edgeR_results$FDR),])
