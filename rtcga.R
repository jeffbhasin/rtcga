#Common Functions for Working with TCGA Data

library(plyr)
library(stringr)
library(reshape)
library(foreach)

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

findStudies <- function()
{
	#############################################################
	#Find Downloaded RNAseq Data
	#############################################################
	#read in list of all disease types
	studies <- read.csv(file="diseaseStudy.csv",header=TRUE)
	names(studies) <- c("key","disease")

	#how many do we have data for?
	downloaded <- data.frame(path=list.dirs(path="BulkDownload",recursive=FALSE))
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
	meta <- read.csv(file="uuidBrowser.csv",header=TRUE,comment.char="",colClasses=rep("character",7),stringsAsFactors=FALSE)

	##do recursive search for all RNA-seq data files - normalized results by gene
	files <- list.files(path="BulkDownload",recursive=TRUE,pattern="\\.genes.normalized_results$",full.names=TRUE)
	length(files)

	#lookup each UUID
	uuids <- str_extract(files,"[a-f0-9]{8}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{12}")

	#join them in
	matched <- join(data.frame(file=files,UUID=uuids),meta)


	##do recursive search for all RNA-seq data files - raw results by gene
	files <- list.files(path="BulkDownload",recursive=TRUE,pattern="\\.genes.results$",full.names=TRUE)
	length(files)

	#lookup each UUID
	uuids <- str_extract(files,"[a-f0-9]{8}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{12}")

	#join them in
	matched2 <- join(matched,data.frame(file.raw=files,UUID=uuids))

	#return final data.frame	
	matched2
}
