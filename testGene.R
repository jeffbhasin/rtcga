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
## Local Functions - move to rtcga.R when done
###############################################
testGenes <- function(samples,samples.data,study,ann)
{
	#run t-test on all genes available - want ranked list of genes for a given study

	#subset the correct study
	samples.study <- samples[samples$Disease==study,]

	#subset the tumor and normal samples
	samples.tumor.rows <- (samples$Sample.Type=="Primary solid Tumor")&(samples$Disease==study)
	samples.tumor <- samples[samples.tumor.rows,]

	samples.normal.rows <- (samples.study$Sample.Type=="Solid Tissue Normal")&(samples$Disease==study)
	samples.normal <- samples[samples.normal.rows,]
	
	#create subsetted matrices
	samples.data.tumor <- list("genes.norm"=samples.data$genes.norm[samples.tumor.rows,], "isoforms.norm"=samples.data$isoforms.norm[samples.tumor.rows,])
	samples.data.normal <- list("genes.norm"=samples.data$genes.norm[samples.normal.rows,], "isoforms.norm"=samples.data$isoforms.norm[samples.normal.rows,])

	#run t-tests
	tt.gene <- foreach(i=1:ncol(samples.data.tumor$genes.norm),.combine=c,.verbose=TRUE) %dopar%
	{
		myt <- t.test(samples.data.tumor$genes.norm[,i],samples.data.normal$genes.norm[,i])
		myt$p.value
	}
	split.pipe <- function(x){strsplit(x,"\\|")[[1]][1]}
	genes.names <- unlist(lapply(colnames(samples.data.tumor$genes.norm),FUN=split.pipe))
	split.pipe2 <- function(x){strsplit(x,"\\|")[[1]][2]}
	entid <- unlist(lapply(colnames(samples.data.tumor$genes.norm),FUN=split.pipe2))
	tt.gene <- data.frame(gene.id=genes.names,entrez.id=entid,p.value=tt.gene)

	tt.isoform <- foreach(i=1:ncol(samples.data.tumor$isoforms.norm),.combine=c,.verbose=TRUE) %dopar%
	{
		myt <- t.test(samples.data.tumor$isoforms.norm[,i],samples.data.normal$isoforms.norm[,i])
		myt$p.value
	}
	split.dot <- function(x){strsplit(x,"\\.")[[1]][1]}
	tt.iso.base <- unlist(lapply(colnames(samples.data.tumor$isoforms.norm),FUN=split.dot))
	tt.isoform <- data.frame(isoform.base.id=tt.iso.base,isoform.id=colnames(samples.data.tumor$isoforms.norm),p.value=tt.isoform)

	#get gene names for isoform IDs from the annotation
	iso.base <- unlist(lapply(ann$name,FUN=split.dot))
	ann.sub <- data.frame(isoform.base.id=iso.base, isoform.id=ann$name,gene.id=ann$geneSymbol,description=ann$description)

	tt.isoform.join <- join(tt.isoform,ann.sub,by="isoform.base.id")
	
	test.out = list("genes.t.test"=tt.gene,"isoforms.t.test"=tt.isoform.join)
	test.out
}

testGenesPaired <- function(samples,samples.data,study)
{
	#TODO
	#run t-test on all genes available - want ranked list of genes for a given study
}

plotGene <- function(samples,samples.data,study,gene)
{
	#subset the correct study
	samples.study <- samples[samples$Disease==study,]

	#subset the tumor and normal samples
	samples.tumor.rows <- (samples$Sample.Type=="Primary solid Tumor")&(samples$Disease==study)
	samples.tumor <- samples[samples.tumor.rows,]

	samples.normal.rows <- (samples.study$Sample.Type=="Solid Tissue Normal")&(samples$Disease==study)
	samples.normal <- samples[samples.normal.rows,]

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

plotGenePaired <- function()
{
	#TODO
}


doQuantileNorm <- function(samples, samples.data, study)
{
	#perform quantile normalization on the raw data and then use this data for getGeneData

	#subset the correct study
	samples.study <- samples[samples$Disease==study,]
	samples.study.rows <- samples$Disease==study

	#create subsetted matrix
	samples.data.study <- list("genes.raw"=samples.data$genes.raw[samples.study.rows,], "isoforms.raw"=samples.data$isoforms.norm[samples.study.rows,])

	quant.all.genes <- normalize.quantiles(samples.data.study$genes.raw,copy=TRUE)
	colnames(quant.all.genes) <- colnames(samples.data.study$genes.raw)
	rownames(quant.all.genes) <- rownames(samples.data.study$genes.raw)

	quant.all.isoforms <- normalize.quantiles(samples.data.study$isoforms.raw,copy=TRUE)
	colnames(quant.all.isoforms) <- colnames(samples.data.study$isoforms.raw)
	rownames(quant.all.isoforms) <- rownames(samples.data.study$isoforms.raw)

	mine.out <- list("genes.norm"=quant.all.genes,"isoforms.norm"=quant.all.isoforms)
	mine.out
}

doLogTransform <- function(samples, samples.data,study)
{
	#subset the correct study
	samples.study <- samples[samples$Disease==study,]
	samples.study.rows <- samples$Disease==study

	#create subsetted matrix
	samples.data.study <- list("genes.norm"=samples.data$genes.norm[samples.study.rows,], "isoforms.norm"=samples.data$isoforms.norm[samples.study.rows,])

	log.all.genes <- log10(samples.data.study$genes.norm)
	log.all.isoforms <- log10(samples.data.study$isoforms.norm)

	mine.out <- list("genes.norm"=log.all.genes,"isoforms.norm"=log.all.isoforms)
	mine.out

}

################################
# Initial Data Loading
################################
ann <- readUCSCAnnotation(genome="hg19",data.path="GeneInfo/")

studies <- findStudies()
samples <- findSamples()

#samples.brca <- samples[which(samples$Disease=="BRCA"),]
#samples.data.brca <- readSampleData(samples.brca)
#save(samples.brca,samples.data.brca,file="output/brca.RData")
load("output/brca.RData")
samples.data.brca.quant <- doQuantileNorm(samples.brca,samples.data.brca,"BRCA")
samples.data.brca.log <- doLogTransform(samples.brca,samples.data.brca,"BRCA")
samples.data.brca.quant.log <- doLogTransform(samples.brca,samples.data.brca.quant,"BRCA")

################################
# Test All Genes in BRCA
################################
test.brca <- testGenes(samples.brca,samples.data.brca,study="BRCA",ann=ann)
write.csv(test.brca$genes.t.test,file="output/BRCA.t-test.genes.csv", row.names=FALSE)
write.csv(test.brca$isoforms.t.test,file="output/BRCA.t-test.isoforms.csv", row.names=FALSE)

################################
# Test All Genes in BRCA - Quantile Normalized
################################
test.brca <- testGenes(samples.brca,samples.data.brca.quant,study="BRCA",ann=ann)
write.csv(test.brca$genes.t.test,file="output/BRCA.t-test.genes-quantile_norm.csv", row.names=FALSE)
write.csv(test.brca$isoforms.t.test,file="output/BRCA.t-test.isoforms-quantile_norm.csv", row.names=FALSE)

################################
# Test All Genes in BRCA - Log Scale
################################
test.brca <- testGenes(samples.brca,samples.data.brca.log,study="BRCA",ann=ann)
write.csv(test.brca$genes.t.test,file="output/BRCA.t-test.genes-log10.csv", row.names=FALSE)
write.csv(test.brca$isoforms.t.test,file="output/BRCA.t-test.isoforms-log10.csv", row.names=FALSE)

################################
# Test All Genes in BRCA - Quantile Normalized Log Scale
################################
test.brca <- testGenes(samples.brca,samples.data.brca.quant.log,study="BRCA",ann=ann)
write.csv(test.brca$genes.t.test,file="output/BRCA.t-test.genes-quantile_norm_log10.csv", row.names=FALSE)
write.csv(test.brca$isoforms.t.test,file="output/BRCA.t-test.isoforms-quantile_norm_log10.csv", row.names=FALSE)

################################
# Test All Genes in PRAD
################################
test.prad <- testGenes(samples.prad,samples.data.prad,study="PRAD",ann=ann)
write.csv(test.prad$genes.t.test,file="output/PRAD.t-test.genes.csv", row.names=FALSE)
write.csv(test.prad$isoforms.t.test,file="output/PRAD.t-test.isoforms.csv", row.names=FALSE)

################################
# Test+Plot Select Genes in BRCA
################################
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
# Test+Plot Select Genes in BRCA - log scale
################################
mystudy <- "BRCA"

mygene <- "SMCHD1"
png(filename=paste("output/",mystudy,".",mygene,"-log.png",sep=""),res=120,width=1000,height=800)
plotGene(samples.brca,samples.data.brca.log,study=mystudy,gene=mygene)
dev.off()

mygene <- "ABCA1"
png(filename=paste("output/",mystudy,".",mygene,"-log.png",sep=""),res=120,width=1000,height=800)
plotGene(samples.brca,samples.data.brca.log,study=mystudy,gene=mygene)
dev.off()

mygene <- "BRCA1"
png(filename=paste("output/",mystudy,".",mygene,"-log.png",sep=""),res=120,width=2200,height=800)
plotGene(samples.brca,samples.data.brca.log,study=mystudy,gene=mygene)
dev.off()

mygene <- "APP"
png(filename=paste("output/",mystudy,".",mygene,"-log.png",sep=""),res=120,width=1800,height=800)
plotGene(samples.brca,samples.data.brca.log,study=mystudy,gene=mygene)
dev.off()

mygene <- "ATM"
png(filename=paste("output/",mystudy,".",mygene,"-log.png",sep=""),res=120,width=1800,height=800)
plotGene(samples.brca,samples.data.brca.log,study=mystudy,gene=mygene)
dev.off()

mygene <- "GAPDH"
png(filename=paste("output/",mystudy,".",mygene,"-log.png",sep=""),res=120,width=1800,height=800)
plotGene(samples.brca,samples.data.brca.log,study=mystudy,gene=mygene)
dev.off()

mygene <- "ACTB"
png(filename=paste("output/",mystudy,".",mygene,"-log.png",sep=""),res=120,width=1800,height=800)
plotGene(samples.brca,samples.data.brca.log,study=mystudy,gene=mygene)
dev.off()

################################
# Test+Plot Select Genes in BRCA - quantile normalized
################################
mystudy <- "BRCA"

mygene <- "SMCHD1"
png(filename=paste("output/",mystudy,".",mygene,"-quantile_norm.png",sep=""),res=120,width=1000,height=800)
plotGene(samples.brca,samples.data.brca.quant,study=mystudy,gene=mygene)
dev.off()

mygene <- "ABCA1"
png(filename=paste("output/",mystudy,".",mygene,"-quantile_norm.png",sep=""),res=120,width=1000,height=800)
plotGene(samples.brca,samples.data.brca.quant,study=mystudy,gene=mygene)
dev.off()

mygene <- "BRCA1"
png(filename=paste("output/",mystudy,".",mygene,"-quantile_norm.png",sep=""),res=120,width=2200,height=800)
plotGene(samples.brca,samples.data.brca.quant,study=mystudy,gene=mygene)
dev.off()

mygene <- "APP"
png(filename=paste("output/",mystudy,".",mygene,"-quantile_norm.png",sep=""),res=120,width=1800,height=800)
plotGene(samples.brca,samples.data.brca.quant,study=mystudy,gene=mygene)
dev.off()

mygene <- "ATM"
png(filename=paste("output/",mystudy,".",mygene,"-quantile_norm.png",sep=""),res=120,width=1800,height=800)
plotGene(samples.brca,samples.data.brca.quant,study=mystudy,gene=mygene)
dev.off()

mygene <- "GAPDH"
png(filename=paste("output/",mystudy,".",mygene,"-quantile_norm.png",sep=""),res=120,width=1800,height=800)
plotGene(samples.brca,samples.data.brca.quant,study=mystudy,gene=mygene)
dev.off()

mygene <- "ACTB"
png(filename=paste("output/",mystudy,".",mygene,"-quantile_norm.png",sep=""),res=120,width=1800,height=800)
plotGene(samples.brca,samples.data.brca.quant,study=mystudy,gene=mygene)
dev.off()

################################
# Test+Plot Select Genes in BRCA - quantile normalized log scale
################################
mystudy <- "BRCA"

mygene <- "SMCHD1"
png(filename=paste("output/",mystudy,".",mygene,"-quantile_norm_log10.png",sep=""),res=120,width=1000,height=800)
plotGene(samples.brca,samples.data.brca.quant.log,study=mystudy,gene=mygene)
dev.off()

mygene <- "ABCA1"
png(filename=paste("output/",mystudy,".",mygene,"-quantile_norm_log10.png",sep=""),res=120,width=1000,height=800)
plotGene(samples.brca,samples.data.brca.quant.log,study=mystudy,gene=mygene)
dev.off()

mygene <- "BRCA1"
png(filename=paste("output/",mystudy,".",mygene,"-quantile_norm_log10.png",sep=""),res=120,width=2200,height=800)
plotGene(samples.brca,samples.data.brca.quant.log,study=mystudy,gene=mygene)
dev.off()

mygene <- "APP"
png(filename=paste("output/",mystudy,".",mygene,"-quantile_norm_log10.png",sep=""),res=120,width=1800,height=800)
plotGene(samples.brca,samples.data.brca.quant.log,study=mystudy,gene=mygene)
dev.off()

mygene <- "ATM"
png(filename=paste("output/",mystudy,".",mygene,"-quantile_norm_log10.png",sep=""),res=120,width=1800,height=800)
plotGene(samples.brca,samples.data.brca.quant.log,study=mystudy,gene=mygene)
dev.off()

mygene <- "GAPDH"
png(filename=paste("output/",mystudy,".",mygene,"-quantile_norm_log10.png",sep=""),res=120,width=1800,height=800)
plotGene(samples.brca,samples.data.brca.quant.log,study=mystudy,gene=mygene)
dev.off()

mygene <- "ACTB"
png(filename=paste("output/",mystudy,".",mygene,"-quantile_norm_log10.png",sep=""),res=120,width=1800,height=800)
plotGene(samples.brca,samples.data.brca.quant.log,study=mystudy,gene=mygene)
dev.off()

################################
# Test+Plot Select Genes in PRAD
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

