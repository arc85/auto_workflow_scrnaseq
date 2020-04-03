##Standard Seurat Workflow - Integration
##TC Feb 6 2020

## @knitr pcaPlot

dat.list <- SplitObject(dat,split.by="chemistry")
dat.list <- lapply(dat.list,function(x) {
	
	x <- NormalizeData(x,verbose=F)
	x <- FindVariableFeatures(x,verbose=F)
	
})

features.use <- SelectIntegrationFeatures(dat.list)
dat.list <- lapply(dat.list,function(x) {
	
	x <- ScaleData(x,features=features.use,verbose=F)
	x <- RunPCA(x,features=features.use,verbose=F)
	
})

anchors <- FindIntegrationAnchors(dat.list,reference=1,reduction="rpca",dims=1:30)
dat.integrated <- IntegrateData(anchorset=anchors,dims=1:30)

dat.integrated <- ScaleData(dat.integrated,verbose=F)
dat.integrated <- RunPCA(dat.integrated,verbose=F)

pc.res <- dat.integrated@reductions$pca@stdev
pc.indices <- which(diff(pc.res)/pc.res[2:50]>-0.01)
pcs.to.include <- round(median(pc.indices[1:3]))

plot(pc.res)
abline(v=pcs.to.include)


## @knitr cluster

#Clustering and UMAP projection
dat.integrated <- FindNeighbors(dat.integrated,dims=1:pcs.to.include)
dat.integrated <- FindClusters(dat.integrated,resolution=c(0.3,0.5,0.7,1.0))
dat.integrated <- RunUMAP(dat.integrated,dims=1:pcs.to.include)
