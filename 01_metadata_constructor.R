

metadata.constructor <- function(raw.data,metadata) {

	if (is.null(dim.raw.dat)) {

		data.to.parse <- raw.data[[1]]

	} else {

		data.to.parse <- raw.data

	}

	num.cells.sample <- vector("logical",length=num.samples)

	for (i in 1:num.samples) {

	if (i==1) {
		length1 <- length(grep("^[A-z]",colnames(data.to.parse)))
		num.cells.sample[i] <- length1
	} else {
		length1 <- length(grep(paste("^",i,sep=""),colnames(data.to.parse)))
		num.cells.sample[i] <- length1
	}

	}

	meta.frame <- data.frame(matrix(data=NA,nrow=sum(num.cells.sample),ncol=ncol(metadata)))
	colnames(meta.frame) <- colnames(metadata)
	rownames(meta.frame) <- colnames(data.to.parse)

	for (i in 1:ncol(meta.frame)) {
		for (z in 1:length(num.cells.sample)) {

			if (z==1) {

				data.to.add <- as.character(metadata[z,i])
				meta.frame[1:num.cells.sample[z],i] <- data.to.add

			} else {

				start.count <- num.cells.sample[z-1]+1
				end.count <- cumsum(num.cells.sample[1:z])[z]

				data.to.add <- as.character(metadata[z,i])

				meta.frame[start.count:end.count,i] <- data.to.add

			}


		}

	}

	meta.frame <- data.frame(lapply(meta.frame,as.factor))
	rownames(meta.frame) <- colnames(data.to.parse)
	return(meta.frame)

}
