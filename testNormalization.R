#Test to find out how the samples are normalized
#Want to know if we can directly compare between samples and/or between studies

#Load common functions
source("Rtcga.R")

###############################################
## Settings
###############################################
nCores <- 15

###############################################
## Packages
###############################################
library(ggplot2)
library(gridExtra)
library(doMC)
registerDoMC(nCores)

studies <- findStudies()
samples <- findSamples()

ggplot.clean <- function()
{
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), legend.key.size = unit(0.8, "lines"), axis.line = element_line(colour = "grey50"))
}
read.table.big <- function(path)
{
	#detects colClasses on first 5 rows to speed up a big read in
	tab5rows <- read.table(file=path, comment.char="", header=TRUE, stringsAsFactors=FALSE, sep="\t", quote="", nrows = 5)
	classes <- sapply(tab5rows, class)
	tabAll <- read.table(file=path, comment.char="", header=TRUE, stringsAsFactors=FALSE, sep="\t", quote="", colClasses = classes)
	tabAll
}

calcSampleMeans <- function(study,samples)
{
	#all means for all genes to each sample table

	#get all tumor and normal samples for the given study
	samples.study <- samples[which((samples$Disease==study)&(samples$Sample.Type=="Solid Tissue Normal" | samples$Sample.Type=="Primary solid Tumor")),]
	
	#join in file names for the raw data as well
	

	##########################################
	#Calculate means for every sample
	##########################################

	#go file by file and collect the sample mean
	samples.study.means <- foreach(i=1:nrow(samples.study), .combine=rbind, .inorder=TRUE, .verbose=FALSE) %dopar%
	{
		sample.data <- read.table.big(as.character(samples.study[i,]$file))		
		mean.norm <- mean(sample.data$normalized_count)

		sample.data.raw <- read.table.big(as.character(samples.study[i,]$file.raw))	
		mean.raw <- mean(sample.data.raw$raw_count)

		data.frame(mean.norm=mean.norm,mean.raw=mean.raw)
	}

	#quantile normalize all data and find means for these
	library(preprocessCore)

	#create matrix of all raw data for all samples
	samples.data.raw <- foreach(i=1:nrow(samples.study), .combine=cbind, .inorder=TRUE, .verbose=FALSE) %dopar%
	{
		sample.data.raw <- read.table.big(as.character(samples.study[i,]$file.raw))
		as.vector(sample.data.raw$raw_count)
	}
	samples.data.raw.normal <- samples.data.raw[,samples.study$Sample.Type=="Solid Tissue Normal"]
	samples.data.raw.tumor <- samples.data.raw[,samples.study$Sample.Type=="Primary solid Tumor"]

	#quant norm for all samples
	quant.all <- normalize.quantiles(samples.data.raw,copy=TRUE)
	mean.quant.all <- colMeans(quant.all)

	#combine means all into one dataframe
	samples.study.out <- data.frame(samples.study,mean.norm=samples.study.means$mean.norm,mean.raw=samples.study.means$mean.raw,mean.quant=mean.quant.all,stringsAsFactors=FALSE)

	#quant norm for just tumor
	quant.tumor <- normalize.quantiles(samples.data.raw.tumor,copy=TRUE)
	mean.quant.tumor <- colMeans(quant.tumor)
	samples.study.out$mean.quant.tumor <- NA
	samples.study.out[samples.study.out$Sample.Type=="Primary solid Tumor",]$mean.quant.tumor <- mean.quant.tumor

	#quant norm for just normals
	quant.norm <- normalize.quantiles(samples.data.raw.normal,copy=TRUE)
	mean.quant.normal <- colMeans(quant.norm)
	samples.study.out$mean.quant.normal <- NA
	samples.study.out[samples.study.out$Sample.Type=="Solid Tissue Normal",]$mean.quant.normal <- mean.quant.normal

	samples.study.out
}

plotSampleMeans <- function(study,samples)
{
	samples.study <- samples[which((samples$Disease==study)&(samples$Sample.Type=="Solid Tissue Normal" | samples$Sample.Type=="Primary solid Tumor")),]

	##########################################
	#plot tumors - tcga normalized
	##########################################

	#read in matrix, with each col a sample and each row a gene
	samples.data.norm <- foreach(i=1:nrow(samples.study), .combine=cbind, .inorder=TRUE, .verbose=TRUE) %dopar%
	{
		sample.data <- read.table.big(as.character(samples.study[i,]$file))
		m <- as.matrix(sample.data$normalized_count)
		rownames(m) <- sample.data$gene_id
		colnames(m) <- as.character(samples.study[i,]$UUID)
		m
	}
	#convert it to a df for easy plotting/subsetting
	data.plot <- melt(samples.data.norm)
	names(data.plot) <- c("gene","UUID","level")
	data.plot$sample <- factor(data.plot$UUID,labels=1:ncol(samples.data.norm))

	#add in tissue types

	types <- data.frame(UUID=samples.study$UUID,tissue=samples.study$Sample.Type,stringsAsFactors=FALSE)
	types[types$tissue=="Solid Tissue Normal",]$tissue <- "Normal"
	types[types$tissue=="Primary solid Tumor",]$tissue <- "Tumor"

	data.plot$tissue <- types[match(data.plot$UUID,types$UUID),]$tissue

	#try plotting just a subset first
	data.plot.sub <- data.plot[1:(nrow(samples.data.norm)*75),]

	png(filename=paste(study,".samples.allgenes.tcga_norm.png",sep=""),width=1500,height=800,res=120)
	p0 <- ggplot(data.plot.sub, aes(sample,level,fill=tissue)) + geom_boxplot(outlier.shape = NA)
	ylim1 = boxplot.stats(data.plot.sub$level)$stats[c(1, 5)]
	p1 = p0 + coord_cartesian(ylim = c(0,ylim1*1.05))
	p1 + ggplot.clean()
	dev.off()

	##########################################
	#plot tumors - raw
	##########################################
	samples.data.raw <- foreach(i=1:nrow(samples.study), .combine=cbind, .inorder=TRUE, .verbose=TRUE) %dopar%
	{
		sample.data <- read.table.big(as.character(samples.study[i,]$file.raw))
		m <- as.matrix(sample.data$raw_count)
		rownames(m) <- sample.data$gene_id
		colnames(m) <- as.character(samples.study[i,]$UUID)
		m
	}
	#convert it to a df for easy plotting/subsetting
	data.plot <- melt(samples.data.raw)
	names(data.plot) <- c("gene","UUID","level")
	data.plot$sample <- factor(data.plot$UUID,labels=1:ncol(samples.data.raw))

	#add in tissue types

	types <- data.frame(UUID=samples.study$UUID,tissue=samples.study$Sample.Type,stringsAsFactors=FALSE)
	types[types$tissue=="Solid Tissue Normal",]$tissue <- "Normal"
	types[types$tissue=="Primary solid Tumor",]$tissue <- "Tumor"

	data.plot$tissue <- types[match(data.plot$UUID,types$UUID),]$tissue

	#try plotting just a subset first
	data.plot.sub <- data.plot[1:(nrow(samples.data.raw)*75),]

	png(filename=paste(study,".samples.allgenes.raw.png",sep=""),width=1500,height=800,res=120)
	p0 <- ggplot(data.plot.sub, aes(sample,level,fill=tissue)) + geom_boxplot(outlier.shape = NA)
	ylim1 = boxplot.stats(data.plot.sub$level)$stats[c(1, 5)]
	p1 = p0 + coord_cartesian(ylim = c(0,ylim1*1.05))
	p1 + ggplot.clean()
	dev.off()

	##########################################
	#plot tumors - quantile norm
	##########################################
	quant.all <- normalize.quantiles(samples.data.raw,copy=TRUE)
	colnames(quant.all) <- colnames(samples.data.raw)
	rownames(quant.all) <- rownames(samples.data.raw)

	#convert it to a df for easy plotting/subsetting
	data.plot <- melt(quant.all)
	names(data.plot) <- c("gene","UUID","level")
	data.plot$sample <- factor(data.plot$UUID,labels=1:ncol(samples.data.raw))

	#add in tissue types

	types <- data.frame(UUID=samples.study$UUID,tissue=samples.study$Sample.Type,stringsAsFactors=FALSE)
	types[types$tissue=="Solid Tissue Normal",]$tissue <- "Normal"
	types[types$tissue=="Primary solid Tumor",]$tissue <- "Tumor"

	data.plot$tissue <- types[match(data.plot$UUID,types$UUID),]$tissue

	#try plotting just a subset first
	data.plot.sub <- data.plot[1:(nrow(samples.data.raw)*75),]

	png(filename=paste(study,".samples.allgenes.quantile.png",sep=""),width=1500,height=800,res=120)
	p0 <- ggplot(data.plot.sub, aes(sample,level,fill=tissue)) + geom_boxplot(outlier.shape = NA)
	ylim1 = boxplot.stats(data.plot.sub$level)$stats[c(1, 5)]
	p1 = p0 + coord_cartesian(ylim = c(0,ylim1*1.05))
	p1 + ggplot.clean()
	#+ facet_grid(. ~ tissue)
	dev.off()

	#plot normals
	#plot pairs

}
plotSampleMeans("COAD",samples)
plotSampleMeans("PRAD",samples)
