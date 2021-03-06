library(scater)
library(scran)
library(edgeR)

options(digits=3)
dge <- read10X(
  mtx="GSM2510617_P7-matrix.mtx.gz",
  genes="GSM2510617_P7-genes.tsv.gz",
  barcodes="GSM2510617_P7-barcodes.tsv.gz",
  DGEList=TRUE)

dim(dge)
head(dge$genes)

ann <- alias2SymbolUsingNCBI(dge$genes$Symbol, required.columns = c("GeneID", "Symbol"), gene.info.file = "Mus_musculus.gene_info.gz") 
dge$genes <- cbind(dge$genes, Official=ann$Symbol, GeneID=ann$GeneID)
head(dge$genes)
ann

dimnames(dge$counts)<-list(dge$genes$Official, colnames(dge$counts))

mito<- grep("^mt-", dge$genes$Symbol)

percent.mito <- colSums(dge$counts[mito, ])/dge$samples$lib.size
nGenes <- colSums(dge$counts != 0)

dge$samples <- cbind(dge$samples, percent.mito=percent.mito, nGenes=nGenes)
head(dge$samples)

par(mfrow=c(1,2))
plot(dge$samples[,c("lib.size","nGenes")], pch=16, cex=0.7)
plot(dge$samples[,c("lib.size","percent.mito")], pch=16, cex=0.7)
o <- order(rowSums(dge$counts), decreasing = TRUE)

dge <- dge[o,]
f1 <- rowSums(dge$counts > 0) >= ncol(dge)*0.01
f2 <- !is.na(dge$genes$Official)
f3 <- !duplicated(dge$genes$Official)
filter <- f1 & f2 & f3
table(filter)

dge.filtered <- dge[filter,]

rownames(dge.filtered) <- dge.filtered$genes$Official


sce <- SingleCellExperiment(list(counts=dge.filtered$counts))
sce


clst <- quickCluster(sce, method="igraph", min.mean=0.5)
table(clst)

sce <- computeSumFactors(sce, cluster=clst)
summary(sizeFactors(sce))

libSize <- dge$samples$lib.size
plot(libSize/1000, sizeFactors(sce), log="xy", pch=16, cex=0.7, xlab="library size (x1000)", ylab="Size fact")




sce <-logNormCounts(sce)



var <- modelGeneVar(sce)
top.var<-var[order(var$bio, decreasing = TRUE),]
top.var<-as.data.frame(top.var)
head(top.var, n=10L)

top <-1000
hvg<-rownames(top.var)[1:top]

means <- rowMeans(logcounts(sce))
vars<- rowVars(logcounts(sce))
fit <- fitTrendVar(means, vars)
names(vars)<-names(means)


plot(means, vars,pch=16,cex=0.7)
curve(fit$trend(x), col="blue", add=TRUE, lwd=2)
points(means[hvg], vars[hvg], col="red", pch=16, cex=0.7)

plotExpression(sce, features=rownames(top.var)[1:10])

sce <- denoisePCA(sce, technical=var$tech, subset.row=hvg, min.rank=10, max.rank=30)

ncol(reducedDim(sce, "PCA"))


colData(sce) <- cbind(colData(sce), log10LibSize=log10(libSize))
plotPCASCE(sce, colour_by = "log10LibSize", ncomponents = 3)

set.seed(100)
sce <- runTSNE(sce, use_dimred="PCA", perplexity=100, theta=0)

plotTSNE(sce, colour_by="log10LibSize")


p1 <- plotTSNE(sce, colour_by="Krt14") + ggtitle("Krt14")
p2 <- plotTSNE(sce, colour_by="Krt18") + ggtitle("Krt18")
multiplot(p1, p2, cols=2)


snn_k60 <- buildSNNGraph(sce, k=60, use.dimred="PCA")
wt_k60 <- igraph::cluster_walktrap(snn_k60)
clst_k60 <- factor(wt_k60$membership)
table(clst_k60)

snn_k70 <- buildSNNGraph(sce, k=70, use.dimred="PCA")
wt_k70 <- igraph::cluster_walktrap(snn_k70)
clst_k70 <- factor(wt_k70$membership)
table(clst_k70)

sce$Cluster <- clst_k60
k60 <- plotTSNE(sce, colour_by="Cluster") + ggtitle("SNN clutering with k=60")

sce$Cluster <- clst_k70
k70 <- plotTSNE(sce, colour_by="Cluster") + ggtitle("SNN clutering with k=70")
multiplot(k60, k70, cols=2)

hclst <- quickCluster(sce, subset.row=hvg, assay.type="logcounts", method="hclust", min.size=50, min.mean=0.5)
table(hclst)

sce$Cluster <- factor(hclst)
plotTSNE(sce, colour_by="Cluster") + ggtitle("Hierarchical clustering")

markers <- findMarkers(sce, group=sce$Cluster, direction="up")

marker.set <- as.data.frame(markers[["1"]])
head(marker.set, n=20L)


top10 <- lapply(markers, "[", 1:10, )
GSet <- as.vector(sapply(top10, "rownames"))
GSet <- GSet[!duplicated(GSet)]
GSet

library(pheatmap)
plotHeatmap(sce, features=GSet, exprs_values="logcounts", color=colorRampPalette(c("blue","white","red"))(100), columns=order(sce$Cluster), labels_col="", colour_columns_by="Cluster", cluster_cols=FALSE, center=TRUE, symmetric=TRUE, zlim=c(-3, 3))
