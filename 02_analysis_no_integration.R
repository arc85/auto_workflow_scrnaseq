##Standard Seurat Workflow - No Integration
##TC Feb 6 2020

## @knitr pcaPlot

dat <- NormalizeData(dat)
dat <- FindVariableFeatures(dat)
dat <- ScaleData(dat)
dat <- RunPCA(dat)

pc.res <- dat@reductions$pca@stdev
pc.indices <- which(diff(pc.res)/pc.res[2:50]>-0.01)
pcs.to.include <- round(median(pc.indices[1:3]))

plot(pc.res)
abline(v=pcs.to.include)


## @knitr cluster

#Clustering and UMAP projection
dat <- FindNeighbors(dat,dims=1:pcs.to.include)
dat <- FindClusters(dat,resolution=c(0.3,0.5,0.7,1.0))
dat <- RunUMAP(dat,dims=1:pcs.to.include)