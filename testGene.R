#Plots for any given gene
#Want to know if we can directly compare between samples and/or between studies

#Load common functions
source("rtcga.R")
source("GeneInfo.R")

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

ann <- readUCSCAnnotation(genome="hg19",data.path="GeneInfo/")

studies <- findStudies()
samples <- findSamples()

testGene <- function(samples,study,gene)
{
	mysamples <- mysamples[which(mysamples$Disease==mystudy),]

	normal.files <- as.character(mysamples[mysamples$Sample.Type=="Solid Tissue Normal",]$file)
	normal <- foreach(i=1:length(normal.files),.combine=c,.verbose=TRUE) %dopar%
	{
		data <- read.table.big(normal.files[i])

		strsplit(data$gene_id,"\\|")

		s <- function(x)
		{
			strsplit(x,"\\|")[[1]][1]
		}
		genes <- unlist(lapply(data$gene_id,FUN=s))

		data$name <- genes

		data[data$name==mygene,]$normalized_count
	}

	tumor.files <- as.character(mysamples[mysamples$Sample.Type=="Primary solid Tumor",]$file)
	tumor <- foreach(i=1:length(tumor.files),.combine=c,.verbose=TRUE) %dopar%
	{
		data <- read.table.big(tumor.files[i])

		strsplit(data$gene_id,"\\|")

		s <- function(x)
		{
			strsplit(x,"\\|")[[1]][1]
		}
		genes <- unlist(lapply(data$gene_id,FUN=s))

		data$name <- genes

		data[data$name==mygene,]$normalized_count
	}

	#t-test
	myt <- t.test(tumor,normal)

	#boxplot
	normal <- data.frame(group="normal",level=normal)
	tumor <- data.frame(group="tumor",level=tumor)
	plot.data <- rbind(tumor,normal)

	png(filename=paste("New.Boxplot.",mygene,".",mystudy,".png",sep=""),res=150,width=1000,height=1000)
	p0 <- ggplot(plot.data, aes(group,level,fill=group)) + geom_boxplot(outlier.shape = NA)
	ylim1 = boxplot.stats(plot.data$level)$stats[c(1, 5)]
	p1 = p0 + coord_cartesian(ylim = c(0,ylim1*1.05))
	p1 + ggplot.clean()
	dev.off()
}

getGeneData <- function(ann,samples,study,gene)
{
	#return matrix of data points for gene and each isoform as list with one for tumor one for normal
	#rows are samples and cols are first the gene and then each isoform

	samples.tumor <- samples[which((samples$Disease==study)&(samples$Sample.Type=="Primary solid Tumor")),]
	samples.normal <- samples[which((samples$Disease==study)&(samples$Sample.Type=="Solid Tissue Normal")),]

	##get combined gene data

	#tumor
	tumor.files <- as.character(samples.tumor$file)
	tumor.data <- foreach(i=1:length(tumor.files),.combine=c,.verbose=TRUE) %dopar%
	{
		data <- read.table.big(tumor.files[i])

		s <- function(x)
		{
			strsplit(x,"\\|")[[1]][1]
		}
		genes <- unlist(lapply(data$gene_id,FUN=s))

		data$name <- genes

		data[data$name==gene,]$normalized_count
	}
	matrix.tumor <- matrix(tumor.data,dimnames=list(samples.tumor$Barcode,gene))

	#normal
	normal.files <- as.character(samples.normal$file)
	normal.data <- foreach(i=1:length(normal.files),.combine=c,.verbose=TRUE) %dopar%
	{
		data <- read.table.big(normal.files[i])

		s <- function(x)
		{
			strsplit(x,"\\|")[[1]][1]
		}
		genes <- unlist(lapply(data$gene_id,FUN=s))

		data$name <- genes

		data[data$name==gene,]$normalized_count
	}
	matrix.normal <- matrix(normal.data,dimnames=list(samples.normal$Barcode,gene))

	##get isoform data
	
	#query annotation for number/list of isoforms
	gene.isoforms <- getAll(ann,gene)
	#strip off version numbers
	s <- function(x){strsplit(x,"\\.")[[1]][1]}
	gene.isoforms.names <- unlist(lapply(gene.isoforms$name,FUN=s))

	#tumor - isoforms
	tumor.files.iso <- as.character(samples.tumor$file.iso)
	tumor.data.iso <- foreach(i=1:length(tumor.files.iso),.combine=rbind,.verbose=TRUE) %dopar%
	{
		data <- read.table.big(tumor.files.iso[i])

		s <- function(x){strsplit(x,"\\.")[[1]][1]}
		isos <- unlist(lapply(data$isoform_id,FUN=s))

		data$name <- isos

		#should have same order as gene.isoforms.names
		data.sub <- data[match(gene.isoforms.names,data$name),]$normalized_count
	}
	colnames(tumor.data.iso) <- gene.isoforms.names
	rownames(tumor.data.iso) <- samples.tumor$Barcode

	#normal - isoforms
	normal.files.iso <- as.character(samples.normal$file.iso)
	normal.data.iso <- foreach(i=1:length(normal.files.iso),.combine=rbind,.verbose=TRUE) %dopar%
	{
		data <- read.table.big(normal.files.iso[i])

		s <- function(x){strsplit(x,"\\.")[[1]][1]}
		isos <- unlist(lapply(data$isoform_id,FUN=s))

		data$name <- isos

		#should have same order as gene.isoforms.names
		data.sub <- data[match(gene.isoforms.names,data$name),]$normalized_count
	}
	colnames(normal.data.iso) <- gene.isoforms.names
	rownames(normal.data.iso) <- samples.normal$Barcode

	#combine gene with isoforms to make one matrix for each type
	matrix.tumor.out <- cbind(matrix.tumor,tumor.data.iso)
	matrix.normal.out <- cbind(matrix.normal,normal.data.iso)

	#join both into named list
	list("tumor"=matrix.tumor.out,"normal"=matrix.normal.out)
}

plotGene <- function(samples.data.matrix)
{
	#Produce a plot after we have extracted the data for the given gene and its isoforms
}


samples.data.matrix <- getGeneData(ann,samples,"BRCA","SMCHD1")
plotGene(samples.data.matrix)


samples.data.matrix <- getGeneData(ann,samples,"PRAD","ABCA1")
plotGene(samples.data.matrix)
