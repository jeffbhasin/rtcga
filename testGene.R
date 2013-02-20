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
plotGene <- function(study,gene,samples.data.matrix)
{
	#Produce a plot after we have extracted the data for the given gene and its isoforms

	#run t-tests
	ts <- foreach(i=1:ncol(samples.data.matrix$tumor),.combine=c) %do%
	{
		myt <- t.test(samples.data.matrix$tumor[,i],samples.data.matrix$normal[,i])
		myt$p.value
	}

	isos <- colnames(samples.data.matrix$normal)
	isos <- isos[order(isos)]
	ts <- round(ts[order(ts)],digits=8)

	table.data <- data.frame()
	table.data <- rbind(table.data,ts)
	names(table.data) <- isos
	rownames(table.data) <- "pvalue"

	#recast data into a table
	plot.data1 <- melt(samples.data.matrix$tumor)
	plot.data1$tissue <- "tumor"

	plot.data2 <- melt(samples.data.matrix$normal)
	plot.data2$tissue <- "normal"

	plot.data <- rbind(plot.data1,plot.data2)
	names(plot.data) <- c("sample","isoform","level","tissue")

	#produce plot
	#p0 <- ggplot(plot.data, aes(isoform,level,fill=tissue)) + geom_boxplot(outlier.shape = NA)
	#ylim1 = boxplot.stats(plot.data$level)$stats[c(1, 5)]
	#p1 = p0 + coord_cartesian(ylim = c(0,ylim1*1.05))
	#p1 + ggplot.clean() + scale_fill_manual(values=c("#009E73","#0072B2")) + labs(title=paste(gene," in ",study,sep=""))

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
# Test SMCHD1 in BRCA
################################
samples.brca <- samples[which(samples$Disease=="BRCA"),]
samples.brca.data <- readSampleData(samples.brca)
save(samples.brca,samples.brca.data,file="output/brca.RData")
load("output/brca.RData")


################################
# Test SMCHD1 in PRAD
################################
samples.prad <- samples[which(samples$Disease=="PRAD"),]
samples.prad.data <- readSampleData(samples.prad)
save(samples.prad,samples.prad.data,file="output/prad.RData")
load("output/prad.RData")


#regular plots
mystudy <- "BRCA"
mygene <- "SMCHD1"
samples.data.matrix <- getGeneData(ann,samples,mystudy,mygene)
png(filename=paste("output/",mystudy,".",mygene,".png",sep=""),res=120,width=1000,height=800)
plotGene(mystudy,mygene,samples.data.matrix)
dev.off()

mystudy <- "BRCA"
mygene <- "ABCA1"
samples.data.matrix <- getGeneData(ann,samples,mystudy,mygene)
png(filename=paste("output/",mystudy,".",mygene,".png",sep=""),res=120,width=1000,height=800)
plotGene(mystudy,mygene,samples.data.matrix)
dev.off()

mystudy <- "BRCA"
mygene <- "BRCA1"
samples.data.matrix <- getGeneData(ann,samples,mystudy,mygene)
png(filename=paste("output/",mystudy,".",mygene,".png",sep=""),res=120,width=1000,height=800)
plotGene(mystudy,mygene,samples.data.matrix)
dev.off()

mystudy <- "PRAD"
mygene <- "SMCHD1"
samples.data.matrix <- getGeneData(ann,samples,mystudy,mygene)
png(filename=paste("output/",mystudy,".",mygene,".png",sep=""),res=120,width=1000,height=800)
plotGene(mystudy,mygene,samples.data.matrix)
dev.off()

mystudy <- "PRAD"
mygene <- "ABCA1"
samples.data.matrix <- getGeneData(ann,samples,mystudy,mygene)
png(filename=paste("output/",mystudy,".",mygene,".png",sep=""),res=120,width=1000,height=800)
plotGene(mystudy,mygene,samples.data.matrix)
dev.off()

mystudy <- "PRAD"
mygene <- "BRCA1"
samples.data.matrix <- getGeneData(ann,samples,mystudy,mygene)
png(filename=paste("output/",mystudy,".",mygene,".png",sep=""),res=120,width=1600,height=800)
plotGene(mystudy,mygene,samples.data.matrix)
dev.off()
