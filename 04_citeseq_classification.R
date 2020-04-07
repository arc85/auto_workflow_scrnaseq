
## @knitr cell_unhashing

#Condition for no cell unhashing
if (unhashing=="NO") {

	print("No cell hashing present")

} else {

#Condition for cell unhashing

#Read in data and set up unhashing lists
cell_unhash <- read.csv("cell_unhashing_identities.csv")

samples <- levels(as.factor(cell_unhash$sample))
samples_dir <- paste(samples,"_Cell_Hash",sep="")

plot_results <- unhash_results <- vector("list",length=length(samples))

for (i in 1:length(plot_results)) {
	plot_results[[i]] <- vector("list",length=2)
	names(plot_results[[i]]) <- c("pre_doublet_removal","post_doublet_removal")
}

names(plot_results) <- names(unhash_results) <- samples

#For each sample, read in cells and identify cutpoints for each hash

for (a in 1:length(samples)) {

	cell_unhash_sub <- cell_unhash %>% filter(sample==samples[a])

	citeseq <- Read10X(paste("./citeseq_matrices/",samples_dir[a],"/read_count/",sep=""),gene.column=1)
	citeseq <- citeseq[!rownames(citeseq)=="unmapped",]

	cite.log <- log1p(citeseq)
	kmeans.res <- cutpoints <- vector("list",length=nrow(cell_unhash_sub))

	for (i in 1:nrow(cell_unhash_sub)) {

		kmeans.res[[i]] <- kmeans(cite.log[i,],centers=2)
		cutpoints[[i]] <- mean(kmeans.res[[i]]$center)

	}

	#Plot to cutpoints

	cite.log.trans <- data.frame(t(data.frame(cite.log)))
	combs.for.plot <- combn(nrow(cell_unhash_sub),2)

	plot_list <- vector("list",length=ncol(combs.for.plot))

	for (i in 1:length(plot_list)) {

		plot.index1 <- combs.for.plot[1,i]
		plot.index2 <- combs.for.plot[2,i]

		plot_list[[i]] <-ggplot(cite.log.trans,aes_string(colnames(cite.log.trans)[plot.index1],colnames(cite.log.trans)[plot.index2])) +
		geom_point(size=0.1) +
		geom_hline(yintercept=cutpoints[[plot.index2]]) +
		geom_vline(xintercept=cutpoints[[plot.index1]])

	}

	plot_results[[a]][[1]] <- plot_list

	#Identify and remove doublets

	citeseq_res <- matrix(data=0,nrow=nrow(cite.log),ncol=ncol(cite.log))

	for (i in 1:nrow(cite.log)) {

		index <- cite.log[i,]>cutpoints[[i]]
		citeseq_res[i,index] <- 1

	}

	citeseq_sums <- colSums(citeseq_res)
	doublet.index <- citeseq_sums>1
	citeseq_res[,doublet.index] <- 0

	citeseq_use <- colSums(citeseq_res)
	use.index <- citeseq_use>0
	cite.log.use <- cite.log[,use.index]

	#Plot hashes and cutpoints without doublets

	cite.log.trans <- data.frame(t(data.frame(cite.log.use)))

	plot_list <- vector("list",length=ncol(combs.for.plot))

	for (i in 1:length(plot_list)) {

		plot.index1 <- combs.for.plot[1,i]
		plot.index2 <- combs.for.plot[2,i]

		plot_list[[i]] <- ggplot(cite.log.trans,aes_string(colnames(cite.log.trans)[plot.index1],colnames(cite.log.trans)[plot.index2])) +
		geom_point(size=0.1) +
		geom_hline(yintercept=cutpoints[[plot.index2]]) +
		geom_vline(xintercept=cutpoints[[plot.index1]])

	}

	plot_results[[a]][[2]] <- plot_list

	#Identify cell unhashed sample

	sample.ident <- apply(citeseq_res,2,function(x) if (any(x>0)) {
	which(x>0) } else {
		0
	})

	level_key <- as.character(paste(samples[a],cell_unhash_sub$individual_sample,sep="_"))
	names(level_key) <- seq_along(1:nrow(cell_unhash_sub))
	sample.ident <- as.factor(recode(sample.ident,!!!level_key,.default="NA"))
	names(sample.ident) <- colnames(cite.log)

	unhash_results[[a]] <- sample.ident

}

#Plot cutpoint results from unhashing

for (i in 1:length(plot_results)) {
	for (a in 1:length(plot_results[[i]])) {
		print(
			wrap_plots(plot_results[[i]][[a]]) +
				plot_annotation(title=paste(names(plot_results)[i],if_else(a==1," with doublets"," final ids"),sep=""))
		)
	}
}

#Incorporate unhashing into metadata for each sample
dat.split <- SplitObject(dat,split.by="sample")

for (i in 1:length(dat.split)) {
	if (i==1) {

	dat.in.unhash <- colnames(dat.split[[i]])[colnames(dat.split[[i]]) %in% names(unhash_results[[i]])]
	unhash.in.dat <- names(unhash_results[[i]])[names(unhash_results[[i]]) %in% colnames(dat.split[[i]])]
	cells.to.use <- intersect(dat.in.unhash,unhash.in.dat)
	dat.split[[i]] <- dat.split[[i]][,cells.to.use]
	unhash_results[[i]] <- unhash_results[[i]][cells.to.use]
	dat.split[[i]]$unhashed.samples <- unhash_results[[i]]

	} else {

		names(unhash_results[[i]]) <- paste(i,names(unhash_results[[i]]),sep="_")
		dat.in.unhash <- colnames(dat.split[[i]])[colnames(dat.split[[i]]) %in% names(unhash_results[[i]])]
		unhash.in.dat <- names(unhash_results[[i]])[names(unhash_results[[i]]) %in% colnames(dat.split[[i]])]
		cells.to.use <- intersect(dat.in.unhash,unhash.in.dat)
		dat.split[[i]] <- dat.split[[i]][,cells.to.use]
		unhash_results[[i]] <- unhash_results[[i]][cells.to.use]
		dat.split[[i]]$unhashed.samples <- unhash_results[[i]]

	}
}

#Recombine Seurat objects

dat1 <- dat.split[[1]]
dat.other <- dat.split[2:length(dat.split)]
dat <- merge(dat1,dat.other)

} #ends toplevel else
