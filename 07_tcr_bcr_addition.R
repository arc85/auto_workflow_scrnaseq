
## @knitr tcr_bcr_addition

#Read in TCR and BCR for each sample

tcr.bcr.files <- list.files("./tcr_bcr_files",recursive=T)

samples <- unique(sapply(strsplit(tcr.bcr.files,split="/"),function(x) x[1]))

tcr.files <- vector("list",length=length(samples))
names(tcr.files) <- samples

if (any(grepl("BCR",tcr.bcr.files))) {
	bcr.files <- vector("list",length=length(samples))
	names(bcr.files) <- samples
}

for (i in 1:length(tcr.files)) {
	tcr.to.read <- grep("TCR",tcr.bcr.files)[i]
	tcr.files[[i]] <- read.csv(paste("./tcr_bcr_files/",tcr.bcr.files[tcr.to.read],sep=""))
	tcr.files[[i]] <- tcr.files[[i]] %>% filter(is_cell=="True",full_length=="True",productive=="True",chain=="TRB") %>% mutate(cell.barcode=sapply(strsplit(as.character(barcode),split="-",),function(x) x[1]))
}

if (any(grepl("BCR",tcr.bcr.files))) {

for (i in 1:length(bcr.files)) {
	bcr.to.read <- grep("BCR",tcr.bcr.files)[i]
	bcr.files[[i]] <- read.csv(paste("./tcr_bcr_files/",tcr.bcr.files[bcr.to.read],sep=""))
	bcr.files[[i]] <- bcr.files[[i]] %>% filter(is_cell=="True",full_length=="True",productive=="True",chain=="IGH") %>% mutate(cell.barcode=sapply(strsplit(as.character(barcode),split="-",),function(x) x[1]))
}

}


#Split data to add per sample TCR and BCR data - remove duplicated barcodes for now

dat.split <- SplitObject(dat,split.by="sample")

#Add TCR data

for (i in 1:length(samples)) {

	if (i==1) {

		sample.cell.names <- colnames(dat.split[[i]])
		tcr.metadata <- data.frame(cell.barcode=sample.cell.names)
		tcr.metadata <- left_join(tcr.metadata,tcr.files[[i]] %>% filter(!duplicated(cell.barcode)) %>%  select(cell.barcode,cdr3),by="cell.barcode")
		rownames(tcr.metadata) <- tcr.metadata$cell.barcode
		tcr.metadata <- tcr.metadata %>% select(-cell.barcode)
		dat.split[[i]]$tcrb.cdr3 <- tcr.metadata

	} else {

		sample.cell.names <- colnames(dat.split[[i]])
		sample.cell.names <- sapply(strsplit(sample.cell.names,split="_"),function(x) x[2])
		tcr.metadata <- data.frame(cell.barcode=sample.cell.names)
		tcr.metadata <- left_join(tcr.metadata,tcr.files[[i]] %>% filter(!duplicated(cell.barcode)) %>%  select(cell.barcode,cdr3),by="cell.barcode")
		rownames(tcr.metadata) <- paste(i,"_",tcr.metadata$cell.barcode,sep="")
		tcr.metadata <- tcr.metadata %>% select(-cell.barcode)
		dat.split[[i]]$tcrb.cdr3 <- tcr.metadata

	}

}


##Add BCR data, if present

if (any(grepl("BCR",tcr.bcr.files))) {

for (i in 1:length(samples)) {

	if (i==1) {

		sample.cell.names <- colnames(dat.split[[i]])
		bcr.metadata <- data.frame(cell.barcode=sample.cell.names)
		bcr.metadata <- left_join(bcr.metadata,bcr.files[[i]] %>% filter(!duplicated(cell.barcode)) %>%  select(cell.barcode,cdr3),by="cell.barcode")
		rownames(bcr.metadata) <- bcr.metadata$cell.barcode
		bcr.metadata <- bcr.metadata %>% select(-cell.barcode)
		dat.split[[i]]$bcr.igh.cdr3 <- tcr.metadata

	} else {

		sample.cell.names <- colnames(dat.split[[i]])
		sample.cell.names <- sapply(strsplit(sample.cell.names,split="_"),function(x) x[2])
		bcr.metadata <- data.frame(cell.barcode=sample.cell.names)
		bcr.metadata <- left_join(bcr.metadata,bcr.files[[i]] %>% filter(!duplicated(cell.barcode)) %>%  select(cell.barcode,cdr3),by="cell.barcode")
		rownames(bcr.metadata) <- bcr.metadata$cell.barcode
		bcr.metadata <- bcr.metadata %>% select(-cell.barcode)
		dat.split[[i]]$bcr.igh.cdr3 <- tcr.metadata

	}

}

}

##Recombine objects

if (length(dat.split)==1) {

	dat <- dat.split

} else {

	dat1 <- dat.split[[1]]
	dat.other <- dat.split[2:length(dat.split)]
	dat <- merge(dat1,dat.other)

}
