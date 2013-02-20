#Common Functions for Working with TCGA Data

library(plyr)
library(stringr)
library(reshape)
library(foreach)

ggplot.clean <- function()
{
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), legend.key.size = unit(0.8, "lines"), axis.line = element_line(colour = "grey50"))
}

head.matrix <- function(mat,n=5)
{
	mat[1:n,1:n]
}

read.table.big <- function(path)
{
	#detects colClasses on first 5 rows to speed up a big read in
	tab5rows <- read.table(file=path, comment.char="", header=TRUE, stringsAsFactors=FALSE, sep="\t", quote="", nrows = 5)
	classes <- sapply(tab5rows, class)
	tabAll <- read.table(file=path, comment.char="", header=TRUE, stringsAsFactors=FALSE, sep="\t", quote="", colClasses = classes)
	tabAll
}

findStudies <- function()
{
	#############################################################
	#Find Downloaded RNAseq Data
	#############################################################
	#read in list of all disease types
	studies <- read.csv(file="input/diseaseStudy.csv",header=TRUE)
	names(studies) <- c("key","disease")

	#how many do we have data for?
	downloaded <- data.frame(path=list.dirs(path="input",recursive=FALSE))
	downloaded$key <- downloaded$path
	downloaded$key <- str_replace(downloaded$key,".*edu_","")
	downloaded$key <- str_replace(downloaded$key,"\\.Illum.*","")

	data <- join(studies,downloaded,type="inner")
	data
}

findSamples <- function()
{
	#############################################################
	#Index Sample Flat Files and Convert Their IDs in File Names to Metadata
	#############################################################
	#read file that maps UUIDs to Barcodes
	meta <- read.csv(file="input/uuidBrowser.csv",header=TRUE,comment.char="",colClasses=rep("character",7),stringsAsFactors=FALSE)

	##do recursive search for all RNA-seq data files - normalized results by gene
	files <- list.files(path="input",recursive=TRUE,pattern="\\.genes.normalized_results$",full.names=TRUE)
	length(files)

	#lookup each UUID
	uuids <- str_extract(files,"[a-f0-9]{8}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{12}")

	#join them in
	matched <- join(data.frame(file=files,UUID=uuids),meta)


	##do recursive search for all RNA-seq data files - raw results by gene
	files <- list.files(path="input",recursive=TRUE,pattern="\\.genes.results$",full.names=TRUE)
	length(files)

	#lookup each UUID
	uuids <- str_extract(files,"[a-f0-9]{8}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{12}")

	#join them in
	matched2 <- join(matched,data.frame(file.raw=files,UUID=uuids))

	##do recursive search for all isoform data files
	files <- list.files(path="input",recursive=TRUE,pattern="\\.isoforms.normalized_results$",full.names=TRUE)
	length(files)

	#lookup each UUID
	uuids <- str_extract(files,"[a-f0-9]{8}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{12}")

	#join them in
	matched3 <- join(matched2,data.frame(file.iso=files,UUID=uuids))

	##do recursive search for all raw isoform data files
	files <- list.files(path="input",recursive=TRUE,pattern="\\.isoforms.results$",full.names=TRUE)
	length(files)

	#lookup each UUID
	uuids <- str_extract(files,"[a-f0-9]{8}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{12}")

	#join them in
	matched4 <- join(matched3,data.frame(file.iso.raw=files,UUID=uuids))

	#return final data.frame	
	matched4
}

readSampleData <- function(samples,study)
{
	#single unified reader function to read in all data from disk once and only once so it can be resued from R in memory
	#if user has memory limits, they can simply send a subset of samples to this reader function based on their needs
	#returns a list of matricies broken down by data types
	#each matrix row corresponds to a row in the input samples data.frame - so this is the key to getting back metadata for future subsetting

	genes.norm.matrix <- foreach(i=1:nrow(samples), .combine=rbind, .inorder=TRUE, .verbose=TRUE) %dopar%
	{
		sample.data <- read.table.big(as.character(samples[i,]$file))
		m <- matrix(sample.data$normalized_count,byrow=TRUE,nrow=1)
		colnames(m) <- sample.data$gene_id
		rownames(m) <- as.character(samples[i,]$UUID)
		m
	}
	isoforms.norm.matrix <- foreach(i=1:nrow(samples), .combine=rbind, .inorder=TRUE, .verbose=TRUE) %dopar%
	{
		sample.data <- read.table.big(as.character(samples[i,]$file.iso))
		m <- matrix(sample.data$normalized_count,byrow=TRUE,nrow=1)
		colnames(m) <- sample.data$isoform_id
		rownames(m) <- as.character(samples[i,]$UUID)
		m
	}
	genes.raw.matrix <- foreach(i=1:nrow(samples), .combine=rbind, .inorder=TRUE, .verbose=TRUE) %dopar%
	{
		sample.data <- read.table.big(as.character(samples[i,]$file.raw))
		m <- matrix(sample.data$raw_count,byrow=TRUE,nrow=1)
		colnames(m) <- sample.data$gene_id
		rownames(m) <- as.character(samples[i,]$UUID)
		m
	}
	isoforms.raw.matrix <- foreach(i=1:nrow(samples), .combine=rbind, .inorder=TRUE, .verbose=TRUE) %dopar%
	{
		sample.data <- read.table.big(as.character(samples[i,]$file.iso.raw))
		m <- matrix(sample.data$raw_count,byrow=TRUE,nrow=1)
		colnames(m) <- sample.data$isoform_id
		rownames(m) <- as.character(samples[i,]$UUID)
		m
	}

	data.out <- list("genes.norm"=genes.norm.matrix, "isoforms.norm"=isoforms.norm.matrix, "genes.raw"=genes.raw.matrix, "isoforms.raw"=isoforms.raw.matrix)
	data.out
}
