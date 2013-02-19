#Plots for any given gene
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

testGene <- function(mysamples,mystudy,mygene)
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
testGene(samples,"BRCA1","SMCHD1")
testGene(samples,"BRCA1","ABCA1")
testGene(samples,"PRAD","ABCA1")

