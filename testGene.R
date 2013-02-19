#For any gene symbol, tell if overexpressed in cancer vs normal both accross types and for each type, if data is available

library(plyr)
library(stringr)
library(reshape)
library(foreach)

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

#############################################################
#Count Sample Types (tumor vs normal) for each study
#############################################################
#read file that maps UUIDs to Barcodes
meta <- read.csv(file="uuidBrowser.csv",header=TRUE,comment.char="",colClasses=rep("character",7),stringsAsFactors=FALSE)

#do recursive search for all RNA-seq data files
files <- list.files(path="BulkDownload",recursive=TRUE,pattern="\\.genes.normalized_results$",full.names=TRUE)
length(files)

#lookup each UUID
uuids <- str_extract(files,"[a-f0-9]{8}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{12}")
matched <- join(data.frame(file=files,UUID=uuids),meta)
nrow(matched)

#count how many UUIDs weren't in the metadata table
matched[is.na(matched$Disease),]

#filter out NAs
matched <- matched[!is.na(matched$Disease),]

#aggregate counts
#total
overall <- table(matched$Sample.Type)
#by disease
rnaseq<-ddply(matched, c("Disease", "Sample.Type"),nrow,.drop=FALSE)
rnaseq <- cast(rnaseq,Disease~Sample.Type)
rnaseq$key <- rnaseq$Disease
rnaseq <- join(studies,rnaseq,type="inner")
rnaseq$key <- NULL
rnaseq$Total <- rowSums(rnaseq[,3:ncol(rnaseq)])

#now can you get me counts of matchings? we need to good at sample ID overlap for each sample type from meta
#look for pairs
#extract the participiant ID from the barcode
#matched$patientid <- substr(matched$Barcode,1,12)

#perform a join of subset for Normal code and subset for Tumor code
normals <- matched[matched$Sample.Type=="Solid Tissue Normal",]
tumors <-matched[matched$Sample.Type=="Primary solid Tumor",]

#check for duplicate patientids
nrow(normals)
nrow(tumors)
length(levels(factor(normals$patientid)))
length(levels(factor(tumors$patientid)))
tumors[duplicated(tumors$patientid),]

#prefix names before doing the join
names(normals) <- paste("normal_",names(normals),sep="")
names(tumors) <- paste("tumor_",names(tumors),sep="")

normals$patientid <- substr(normals$normal_Barcode,1,12)
tumors$patientid <- substr(tumors$tumor_Barcode,1,12)

#join to make data frame with one entry for each combination of shared sampleids
pairs <- join(tumors,normals,by="patientid",type="inner")

#create summary table of number of pairs by cancer type
rnaseq_pairs <- ddply(pairs, c("tumor_Disease"),nrow,.drop=FALSE)
names(rnaseq_pairs) <- c("Disease","Pairs")
rnaseq <- join(rnaseq,rnaseq_pairs)
rnaseq$Pairs[is.na(rnaseq$Pairs)] <- 0

#do final row of sums
sumrow <- data.frame("Total","",t(colSums(rnaseq[,3:ncol(rnaseq)])))
names(sumrow) <- names(rnaseq)
rnaseq <- rbind(rnaseq,sumrow)

#also want to know how much data is available for each one: tumor, normal, matching
#make a summary table of everything that's available
write.csv(rnaseq,file="TCGA-RNAseq-Avail.csv",row.names=FALSE)


#############################################################
#T-Test/Boxplot of Overall Average Expression in Tumor vs Normal
#############################################################

#mygene = "SMCHD1"

#bump out list of all files for normal
#write.table(matched[matched$Sample.Type=="Solid Tissue Normal",]$file, file="normal.tmp.txt", append=FALSE, quote=FALSE, row.names=FALSE, col.names=FALSE)

#bump out list of all files for tumor
#write.table(matched[matched$Sample.Type=="Primary solid Tumor",]$file,file="tumor.tmp.txt", append=FALSE, quote=FALSE, row.names=FALSE, col.names=FALSE)

#use a bash grep to extract all data for a given gene quickly
#system(paste("cat normal.tmp.txt | xargs cat | grep '",mygene,"' > normal.data.tmp.txt",sep=""))
#system(paste("cat tumor.tmp.txt | xargs cat | grep '",mygene,"' > tumor.data.tmp.txt",sep=""))

#read in the data we extracted
#normal <- read.table(file="normal.data.tmp.txt",header=FALSE)$V2
#tumor <- read.table(file="tumor.data.tmp.txt",header=FALSE)$V2

#t-test
#t.test(tumor,normal)

#boxplot
#normal <- data.frame(group="normal",level=normal)
#tumor <- data.frame(group="tumor",level=tumor)
#all <- rbind(tumor,normal)

#png(filename=paste("Boxplot.",mygene,".png",sep=""),res=150,width=1000,height=1000)
#boxplot(level~group,data=all, main=paste(mygene," RNAseq Expression",sep=""), xlab="Group", ylab="Exp Level")
#dev.off()

testGene <- function(mygene)
{
	#bump out list of all files for normal
	write.table(matched[matched$Sample.Type=="Solid Tissue Normal",]$file, file="normal.tmp.txt", append=FALSE, quote=FALSE, row.names=FALSE, col.names=FALSE)

	#bump out list of all files for tumor
	write.table(matched[matched$Sample.Type=="Primary solid Tumor",]$file,file="tumor.tmp.txt", append=FALSE, quote=FALSE, row.names=FALSE, col.names=FALSE)
	#use a bash grep to extract all data for a given gene quickly
	system(paste("cat normal.tmp.txt | xargs cat | grep '",mygene,"|' > normal.data.tmp.txt",sep=""))
	system(paste("cat tumor.tmp.txt | xargs cat | grep '",mygene,"|' > tumor.data.tmp.txt",sep=""))

	#read in the data we extracted
	normal <- read.table(file="normal.data.tmp.txt",header=FALSE)$V2
	tumor <- read.table(file="tumor.data.tmp.txt",header=FALSE)$V2

	#t-test
	myt <- t.test(tumor,normal)

	#boxplot
	normal <- data.frame(group="normal",level=normal)
	tumor <- data.frame(group="tumor",level=tumor)
	all <- rbind(tumor,normal)

	png(filename=paste("Boxplot.",mygene,".png",sep=""),res=150,width=1000,height=1000)
	boxplot(level~group,data=all, main=paste(mygene," RNAseq Exp (p=",myt$p.value,")",sep=""), xlab="Group", ylab="Exp Level", outline=FALSE)
	dev.off()

		system("rm normal.tmp.txt tumor.tmp.txt normal.data.tmp.txt tumor.data.tmp.txt")
}


testGene("SMCHD1")
testGene("BRCA1")
testGene("MYC")
testGene("PAX9")
testGene("SMURF1")
testGene("GAPDH")
testGene("ABCA1")


testGeneInStudy <- function(mygene,mystudy)
{
	#bump out list of all files for normal
	write.table(matched[matched$Sample.Type=="Solid Tissue Normal",]$file, file="normal.tmp.txt", append=FALSE, quote=FALSE, row.names=FALSE, col.names=FALSE)

	#bump out list of all files for tumor
	write.table(matched[matched$Sample.Type=="Primary solid Tumor",]$file,file="tumor.tmp.txt", append=FALSE, quote=FALSE, row.names=FALSE, col.names=FALSE)
	#use a bash grep to extract all data for a given gene quickly
	system(paste("cat normal.tmp.txt | grep '",mystudy,"' | xargs cat | grep '",mygene,"|' > normal.data.tmp.txt",sep=""))
	system(paste("cat tumor.tmp.txt | grep '",mystudy,"' | xargs cat | grep '",mygene,"|' > tumor.data.tmp.txt",sep=""))

	#read in the data we extracted
	normal <- read.table(file="normal.data.tmp.txt",header=FALSE)$V2
	tumor <- read.table(file="tumor.data.tmp.txt",header=FALSE)$V2

	#t-test
	myt <- t.test(tumor,normal)

	#boxplot
	normal <- data.frame(group="normal",level=normal)
	tumor <- data.frame(group="tumor",level=tumor)
	all <- rbind(tumor,normal)

	png(filename=paste("Boxplot.",mygene,".",mystudy,".png",sep=""),res=150,width=1000,height=1000)
	boxplot(level~group,data=all, main=paste(mygene," in ",mystudy," RNAseq Exp (p=",myt$p.value,")",sep=""), xlab="Group", ylab="Exp Level", outline=FALSE)
	dev.off()

		system("rm normal.tmp.txt tumor.tmp.txt normal.data.tmp.txt tumor.data.tmp.txt")
}
testGeneInStudy("ABCA1","PRAD")
testGeneInStudy("ABCA1","BRCA")


#test in all studies and make a table of t-tests and one big boxplot
testGeneInAllStudies <- function(mygene)
{
	mydata <- data.frame()
	mysummary <- data.frame()
	#do for each study we know about
	foreach(mystudy=as.character(data$key),.combine="c") %do%
	{
		#if either group has no data, skip this sample
		if((rnaseq[rnaseq$Disease==mystudy,]$"Primary solid Tumor" == 0) | (rnaseq[rnaseq$Disease==mystudy,]$"Solid Tissue Normal") == 0)
		{
			print(paste(mystudy,"has no data"))
		}		
		else 
		{
		#bump out list of all files for normal
		write.table(matched[matched$Sample.Type=="Solid Tissue Normal",]$file, file="normal.tmp.txt", append=FALSE, quote=FALSE, row.names=FALSE, col.names=FALSE)
		#bump out list of all files for tumor
		write.table(matched[matched$Sample.Type=="Primary solid Tumor",]$file,file="tumor.tmp.txt", append=FALSE, quote=FALSE, row.names=FALSE, col.names=FALSE)
		#use a bash grep to extract all data for a given gene quickly
		system(paste("cat normal.tmp.txt | grep '",mystudy,"' | xargs cat | grep '",mygene,"|' > normal.data.tmp.txt",sep=""))
		system(paste("cat tumor.tmp.txt | grep '",mystudy,"' | xargs cat | grep '",mygene,"|' > tumor.data.tmp.txt",sep=""))

		#read in the data we extracted
		normal <- data.frame()
		normal <- read.table(file="normal.data.tmp.txt",header=FALSE)$V2
		tumor <- data.frame()
		tumor <- read.table(file="tumor.data.tmp.txt",header=FALSE)$V2

		#clean up
		system("rm normal.tmp.txt tumor.tmp.txt normal.data.tmp.txt tumor.data.tmp.txt")


		#t-test
		myt <- "NA"
		if((length(normal) > 0) & (length(tumor)>0))
		{
			myt <- t.test(tumor,normal)$p.value
		}


		direction <- "-"
		if(mean(tumor)>mean(normal))
		{
			direction <- "+"
		}

		normal2 <- data.frame(group=paste(mystudy,"normal",sep="_"),level=normal)
		tumor2 <- data.frame(group=paste(mystudy,"tumor",sep="_"),level=tumor)
		all <- rbind(tumor2,normal2)

		mydata<-rbind(mydata,all)
		mysummary<-rbind(mysummary, data.frame(mystudy,length(tumor),mean(tumor),length(normal),mean(normal),direction,myt))
		}
	}

	#make a big boxplot for all studies
	png(filename=paste("Boxplot.",mygene,".AllStudies.png",sep=""),width=1500,height=1000)
	boxplot(level~group,data=mydata, main=paste(mygene," RNAseq Exp All Studies",sep=""), xlab="Group", ylab="Exp Level", outline=FALSE)
	dev.off()

	#write CSV summary of all studies
	names(mysummary) <- c("key","nTumor","meanTumorExp","nNormal","meanNormalExp","direction","p")
	mysummary <- join(mysummary,data)
	write.csv(mysummary,file=paste("Summary.",mygene,".AllStudies.csv",sep=""),row.names=FALSE)

}
testGeneInAllStudies("SMCHD1")
testGeneInAllStudies("ABCA1")
testGeneInAllStudies("GAPDH")
testGeneInAllStudies("BRCA1")
testGeneInAllStudies("BRCA2")

testGenePairs <- function(mygene)
{
	#make boxplot of ratios (fold changes) for all paired samples looking at tumor/normal

	#bump out list of all files for normal
	write.table(pairs$normal_file, file="normal.tmp.txt", append=FALSE, quote=FALSE, row.names=FALSE, col.names=FALSE)

	#bump out list of all files for tumor
	write.table(pairs$tumor_file,file="tumor.tmp.txt", append=FALSE, quote=FALSE, row.names=FALSE, col.names=FALSE)

	#use a bash grep to extract all data for a given gene quickly
	system(paste("cat normal.tmp.txt | xargs cat | grep '",mygene,"|' > normal.data.tmp.txt",sep=""))
	system(paste("cat tumor.tmp.txt | xargs cat | grep '",mygene,"|' > tumor.data.tmp.txt",sep=""))

	#read in the data we extracted
	normal <- read.table(file="normal.data.tmp.txt",header=FALSE)$V2
	tumor <- read.table(file="tumor.data.tmp.txt",header=FALSE)$V2

	log2ratios <- log2(tumor/normal)
	ratios <- data.frame(group="log2(tumor/normal)",ratio=log2ratios)

	png(filename=paste("Boxplot.",mygene,".pairs.png",sep=""),res=150,width=1000,height=1000)
	boxplot(ratio~group,data=ratios, main=paste(mygene," RNAseq Exp Pairs (n=",nrow(ratios),")",sep=""), xlab="Group", ylab="log2(tumor exp/normal exp)", outline=FALSE)
	dev.off()
}
testGenePairs("SMCHD1")
testGenePairs("ABCA1")
testGenePairs("BRCA2")

#############################################################
#Run Global Tumor vs Normal T-Tests for EVERY GENE just for Kicks
#############################################################

#want a CSV output I can sort by p-value
#we may have to import it all to SQL in order to do this easily

#############################################################
#T-Test/Boxplot of Per-Disease Expression in Tumor vs Normal
#############################################################

#############################################################
#Some way to handle matched samples?
#############################################################

#how to ID them? need to find patient ids that have more than one sample type
