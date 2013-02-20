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
library(reshape)
library(preprocessCore)
library(doMC)
registerDoMC(nCores)

###############################################
## Local Functions
###############################################
testGenes <- function()
{
	#run t-test on all genes available - want ranked list of genes for a given study
}
testGeneIsoforms <- function()
{
	#run t-test on all gene isoforms available, using annotation to connect IDs back to gene names
}
testGenesPaired <- function()
{
	#run t-test on all genes available - want ranked list of genes for a given study
}
testGeneIsoformsPaired <- function()
{
	#run t-test on all gene isoforms available, using annotation to connect IDs back to gene names
}

plotGene <- function(samples,samples.data,study,gene)
{
	#subset the correct study
	samples.study <- samples[samples$Disease==study,]

	#subset the tumor and normal samples
	samples.tumor.rows <- samples.study$Sample.Type=="Primary solid Tumor"
	samples.tumor <- samples.study[samples.tumor.rows,]

	samples.normal.rows <- samples.study$Sample.Type=="Solid Tissue Normal"
	samples.normal <- samples.study[samples.normal.rows,]

	#subset the desired gene
	split.pipe <- function(x){strsplit(x,"\\|")[[1]][1]}
	genes.names <- unlist(lapply(colnames(samples.data$genes.norm),FUN=split.pipe))
	genes.col <- genes.names==gene

	#subset the desired isoforms
	split.dot <- function(x){strsplit(x,"\\.")[[1]][1]}
	isoforms.names <- unlist(lapply(colnames(samples.data$isoforms.norm),FUN=split.dot))
	
	#get isoforms for the gene symbol from annotation
	gene.isoforms <- getAll(ann,gene)
	isoforms.this.gene <- unlist(lapply(gene.isoforms$name,FUN=split.dot))

	isoforms.cols <- isoforms.names %in% isoforms.this.gene

	#create subsetted matrices
	genesub <- matrix(samples.data$genes.norm[samples.tumor.rows,genes.col])
	colnames(genesub) <- gene
	samples.data.tumor <- cbind(genesub,samples.data$isoforms.norm[samples.tumor.rows,isoforms.cols])

	genesub <- matrix(samples.data$genes.norm[samples.normal.rows,genes.col])
	colnames(genesub) <- gene
	samples.data.normal <- cbind(genesub,samples.data$isoforms.norm[samples.normal.rows,isoforms.cols])

	#run t-tests
	ts <- foreach(i=1:ncol(samples.data.tumor),.combine=c) %do%
	{
		myt <- t.test(samples.data.tumor[,i],samples.data.normal[,i])
		myt$p.value
	}

	isos <- colnames(samples.data.normal)
	isos <- isos[order(isos)]
	ts <- round(ts[order(ts)],digits=8)

	table.data <- data.frame()
	table.data <- rbind(table.data,ts)
	names(table.data) <- isos
	rownames(table.data) <- "pvalue"

	#recast data into a table
	plot.data1 <- melt(samples.data.tumor)
	plot.data1$tissue <- "tumor"

	plot.data2 <- melt(samples.data.normal)
	plot.data2$tissue <- "normal"

	plot.data <- rbind(plot.data1,plot.data2)
	names(plot.data) <- c("sample","isoform","level","tissue")

	#produce plot
	p1 <- ggplot(plot.data, aes(isoform,level,fill=tissue)) + geom_boxplot() + ggplot.clean() + scale_fill_manual(values=c("#009E73","#0072B2")) + labs(title=paste(gene," in ",study,sep=""))

	table.grob <- arrangeGrob(tableGrob(table.data, gpar.coretext=gpar(fontsize=10), gpar.coltext=gpar(fontsize=10), gpar.rowtext=gpar(fontsize=10)))
	grid.arrange(p1, table.grob, ncol=1,heights=c(9/10,1/10))

}



doQuantileNorm <- function(samples,study)
{
	#perform quantile normalization on the raw data and then use this data for getGeneData
	samples.study <- samples[which((samples$Disease==study)&(samples$Sample.Type=="Solid Tissue Normal" | samples$Sample.Type=="Primary solid Tumor")),]

	samples.data.raw <- foreach(i=1:nrow(samples.study), .combine=cbind, .inorder=TRUE, .verbose=TRUE) %dopar%
	{
		sample.data <- read.table.big(as.character(samples.study[i,]$file.iso.raw))
		m <- as.matrix(sample.data$raw_count)
		rownames(m) <- sample.data$isoform_id
		colnames(m) <- as.character(samples.study[i,]$UUID)
		m
	}

	quant.all <- normalize.quantiles(samples.data.raw,copy=TRUE)

}

getGeneDataQuantile <- function(samples,study,gene)
{

}

################################
# Initial Data Loading
################################
ann <- readUCSCAnnotation(genome="hg19",data.path="GeneInfo/")

studies <- findStudies()
samples <- findSamples()

################################
# Test Select Genes in BRCA
################################

#samples.brca <- samples[which(samples$Disease=="BRCA"),]
#samples.data.brca <- readSampleData(samples.brca)
#save(samples.brca,samples.data.brca,file="output/brca.RData")
load("output/brca.RData")

mystudy <- "BRCA"

mygene <- "SMCHD1"
png(filename=paste("output/",mystudy,".",mygene,".png",sep=""),res=120,width=1000,height=800)
plotGene(samples.brca,samples.data.brca,study=mystudy,gene=mygene)
dev.off()

mygene <- "ABCA1"
png(filename=paste("output/",mystudy,".",mygene,".png",sep=""),res=120,width=1000,height=800)
plotGene(samples.brca,samples.data.brca,study=mystudy,gene=mygene)
dev.off()

mygene <- "BRCA1"
png(filename=paste("output/",mystudy,".",mygene,".png",sep=""),res=120,width=2200,height=800)
plotGene(samples.brca,samples.data.brca,study=mystudy,gene=mygene)
dev.off()

mygene <- "APP"
png(filename=paste("output/",mystudy,".",mygene,".png",sep=""),res=120,width=1800,height=800)
plotGene(samples.brca,samples.data.brca,study=mystudy,gene=mygene)
dev.off()

mygene <- "ATM"
png(filename=paste("output/",mystudy,".",mygene,".png",sep=""),res=120,width=1800,height=800)
plotGene(samples.brca,samples.data.brca,study=mystudy,gene=mygene)
dev.off()

mygene <- "GAPDH"
png(filename=paste("output/",mystudy,".",mygene,".png",sep=""),res=120,width=1800,height=800)
plotGene(samples.brca,samples.data.brca,study=mystudy,gene=mygene)
dev.off()

mygene <- "ACTB"
png(filename=paste("output/",mystudy,".",mygene,".png",sep=""),res=120,width=1800,height=800)
plotGene(samples.brca,samples.data.brca,study=mystudy,gene=mygene)
dev.off()

################################
# Test SMCHD1 in PRAD
################################

#samples.prad <- samples[which(samples$Disease=="PRAD"),]
#samples.prad.data <- readSampleData(samples.prad)
#save(samples.prad,samples.data.prad,file="output/prad.RData")
load("output/prad.RData")

mystudy <- "PRAD"

mygene <- "SMCHD1"
png(filename=paste("output/",mystudy,".",mygene,".png",sep=""),res=120,width=1000,height=800)
plotGene(samples.prad,samples.data.prad,study=mystudy,gene=mygene)
dev.off()

mygene <- "ABCA1"
png(filename=paste("output/",mystudy,".",mygene,".png",sep=""),res=120,width=1000,height=800)
plotGene(samples.prad,samples.data.prad,study=mystudy,gene=mygene)
dev.off()

mygene <- "BRCA1"
png(filename=paste("output/",mystudy,".",mygene,".png",sep=""),res=120,width=2200,height=800)
plotGene(samples.prad,samples.data.prad,study=mystudy,gene=mygene)
dev.off()

mygene <- "APP"
png(filename=paste("output/",mystudy,".",mygene,".png",sep=""),res=120,width=1800,height=800)
plotGene(samples.prad,samples.data.prad,study=mystudy,gene=mygene)
dev.off()

mygene <- "ATM"
png(filename=paste("output/",mystudy,".",mygene,".png",sep=""),res=120,width=1800,height=800)
plotGene(samples.prad,samples.data.prad,study=mystudy,gene=mygene)
dev.off()

mygene <- "GAPDH"
png(filename=paste("output/",mystudy,".",mygene,".png",sep=""),res=120,width=1800,height=800)
plotGene(samples.prad,samples.data.prad,study=mystudy,gene=mygene)
dev.off()

mygene <- "ACTB"
png(filename=paste("output/",mystudy,".",mygene,".png",sep=""),res=120,width=1800,height=800)
plotGene(samples.prad,samples.data.prad,study=mystudy,gene=mygene)
dev.off()

