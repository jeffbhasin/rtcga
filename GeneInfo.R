#Given a Gene Symbol, Find Basic Facts About the Gene From the Annotation

#load libraries
library(plyr)
library(GenomicRanges)

#########################################
# Functions
#########################################
read.table.ucsc <- function(path)
{
	#one function to read in all tables exported from UCSC the same way
	read.table(file=path, comment.char="", header=TRUE, stringsAsFactors=FALSE, sep="\t", quote="")
}
read.table.ucsc.big <- function(path)
{
	#detects colClasses on first 5 rows to speed up a big read in
	tab5rows <- read.table(file=path, comment.char="", header=TRUE, stringsAsFactors=FALSE, sep="\t", quote="", nrows = 5)
	classes <- sapply(tab5rows, class)
	tabAll <- read.table(file=path, comment.char="", header=TRUE, stringsAsFactors=FALSE, sep="\t", quote="", colClasses = classes)
	tabAll
}
readEnsemblAnnotation <- function()
{
	#downloaded from UCSC genome browser "tables" page
	ensGene.path <- "input/ensGene.txt"
	ensemblToGeneName.path <- "input/ensemblToGeneName.txt"
	ccdsInfo.path <- "input/ccdsInfo.txt"

	#read in annotation tables
	ensGene <- read.table(file=ensGene.path, comment.char="", header=TRUE, stringsAsFactors=FALSE)
	ensemblToGeneName <- read.table(file=ensemblToGeneName.path, comment.char="", header=TRUE, stringsAsFactors=FALSE)
	names(ensemblToGeneName) <- c("name","symbol")
	ccdsInfo <- read.table(file=ccdsInfo.path, comment.char="", header=TRUE, stringsAsFactors=FALSE)

	#how many pairs of transcript IDs and names?
	nNamePairs <- nrow(ensemblToGeneName)
	nNamePairs

	#how many unique named genes?
	nUniqueGeneSymbols <- length(levels(factor(ensemblToGeneName$value)))
	nUniqueGeneSymbols

	#how many unique transcript IDs?
	nUniqueIDs <- length(levels(factor(ensemblToGeneName$name)))
	nUniqueIDs

	#how many unique trancript IDs in ensGene?
	ensGene.nUniqueIDs <- length(levels(factor(ensGene$name)))
	ensGene.nUniqueIDs

	#do all ensGene transcript IDs match those in the name map table?
	table(ensGene$name %in% ensemblToGeneName$name)

	#combine names into ensGene
	ens.ann <- join(ensGene,ensemblToGeneName,by="name",type="left")

	#how many consensus coding sequences are there?
	nccds <- length(levels(factor(ccdsInfo$X.ccds)))
	nccds

	ccdsInfo <- data.frame(name=ccdsInfo$mrnaAcc,ccds=ccdsInfo$X.ccds,stringsAsFactors=FALSE)

	#combine consensus coding sequence into ensGene - labels which variant is consensus for gene symbols that appear twice
	ens.ann <- join(ens.ann,ccdsInfo,by="name",type="left")	

	#convert coordinates to 1-based from UCSC's 0-based backend
	ens.ann$txStart <- ens.ann$txStart + 1

	ens.ann$cdsStart <- ens.ann$cdsStart + 1

	#not going to add the offet to exons at this point, better to do so when we actually parse the exons to save time

	ens.ann
}

downloadUCSCAnnotation <- function()
{
	#http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/knownGene.txt.gz
	#download latest versions, gunzip, then diff with latest and tell me if there has been an update or not, if there has been, save these new versions with the data appended
}

readUCSCAnnotation <- function(genome="hg19",data.path="")
{
	#load downloaded tables
	knownGene.path <- paste(data.path,genome,".knownGene.txt",sep="")
	knownCanonical.path <- paste(data.path,genome,".knownCanonical.txt",sep="")
	kgXref.path <- paste(data.path,genome,".kgXref.txt",sep="")
	cytoBand.path <- paste(data.path,genome,".cytoBand.txt",sep="")

	knownGene <- read.table.ucsc(knownGene.path)
	knownCanonical <- read.table.ucsc(knownCanonical.path)
	kgXref <- read.table.ucsc(kgXref.path)
	cytoBand <- read.table.ucsc(cytoBand.path)

	#join in gene symbols
	kg.sub <- knownGene
	names(kg.sub)[1] <- "name" 
	kgXref.sub <- kgXref
	names(kgXref.sub)[1] <- "name"

	kg.ann.1 <- join(kg.sub,kgXref.sub,by="name",type="left")

	#join in canonical annotation (largest coding seq from nearest cluster)
	kgCanon.sub <- data.frame(name=knownCanonical$transcript,canonical="1",stringsAsFactors=FALSE)
	
	kg.ann.2 <- join(kg.ann.1,kgCanon.sub,by="name",type="left")

	kg.ann <- kg.ann.2

	#join in cytological bands based on gene position
	kg.ranges <- GRanges(seqnames=kg.ann$chrom,ranges=IRanges(start=kg.ann$txStart,end=kg.ann$txEnd),strand=kg.ann$strand,id=kg.ann$name)
	cb.ranges <- with(cytoBand, GRanges(seqnames=X.chrom,ranges=IRanges(start=chromStart,end=chromEnd),band=name,gstain=gieStain))
	
	#band genes start in
	kg.ranges$CytoBandStartIndex <- findOverlaps(kg.ranges,cb.ranges,select="first")
	#band genes end in (in case they span more than one)
	kg.ranges$CytoBandEndIndex <- findOverlaps(kg.ranges,cb.ranges,select="last")
	#total bands spanned
	kg.ranges$CytoBandsSpan <- countOverlaps(kg.ranges,cb.ranges)
	#create DF from GR metadata
	cyto.map <- data.frame(name=kg.ranges$id,count=kg.ranges$CytoBandsSpan,start=kg.ranges$CytoBandStartIndex,end=kg.ranges$CytoBandEndIndex,stringsAsFactors=FALSE)

	cyto.map$startname <- "NA"
	cyto.map[is.na(cyto.map$start)==FALSE,]$startname <- cb.ranges[cyto.map[is.na(cyto.map$start)==FALSE,]$start]$band
	cyto.map$endname <- "NA"
	cyto.map[is.na(cyto.map$start)==FALSE,]$endname <- cb.ranges[cyto.map[is.na(cyto.map$start)==FALSE,]$end]$band
	cyto.map <- data.frame(name=cyto.map$name,cytoBandSpan=cyto.map$count,cytoBandStart=cyto.map$startname,cytoBandEnd=cyto.map$endname,stringsAsFactors=FALSE)

	#join cyto map data back into main annotation by id
	kg.ann <- join(kg.ann,cyto.map,by="name",type="left")


	#convert 0-based start coordinates to 1-based
	kg.ann$txStart <- kg.ann$txStart + 1

	kg.ann$cdsStart <- kg.ann$cdsStart + 1

	kg.ann
}

getAll <- function(ann,gene.symbol)
{
	ann[ann$geneSymbol==gene.symbol,]
}

getGenePosition <- function(ann,gene.symbol)
{
	#returns coordinates for canonical splice variant
	this.ann <- ann[which((ann$canonical=="1")&(ann$geneSymbol==gene.symbol)),]
	data.frame(this.ann$chrom,this.ann$txStart,this.ann$txEnd)
}

getGeneCytoBands <- function(ann,gene.symbol)
{
	this.ann <- ann[which((ann$canonical=="1")&(ann$geneSymbol==gene.symbol)),]
	data.frame(this.ann$cytoBandSpan,this.ann$cytoBandStart,this.ann$cytoBandEnd)

}

getGeneSize <- function(ann,gene.symbol)
{
	#size of canonical splice variant
	this.ann <- ann[which((ann$canonical=="1")&(ann$geneSymbol==gene.symbol)),]
	this.ann$txEnd-this.ann$txStart	
}

getExonNumber <- function(ann,gene.symbol)
{
	this.ann <- ann[which((ann$canonical=="1")&(ann$geneSymbol==gene.symbol)),]
	this.ann$exonCount
}

getIntrons <- function(ann,gene.symbol)
{
	#list of intron coordinates
	this.ann <- ann[which((ann$canonical=="1")&(ann$geneSymbol==gene.symbol)),]
	
	#introns should just be the negative ranges of the exons
	exons <- getExons(ann,gene.symbol)

	exons.ranges <- IRanges(start=exons$exonStart,end=exons$exonEnd)
	exons.gaps <- gaps(exons.ranges)

	data.frame(intronStart=start(exons.gaps),intronEnd=end(exons.gaps),intronSize=width(exons.gaps)+1)
}

getExons <- function(ann,gene.symbol)
{
	#list of exon coordinates
	this.ann <- ann[which((ann$canonical=="1")&(ann$geneSymbol==gene.symbol)),]
	starts <- unlist(strsplit(this.ann$exonStarts,","))
	starts <- as.numeric(starts)+1
	ends <- as.numeric(unlist(strsplit(this.ann$exonEnds,",")))
	sizes <- ends-starts	

	data.frame(exonStart=starts,exonEnd=ends,exonSize=sizes)
}

getFlankingGenes <- function(ann,gene.symbol,flankSize=1000)
{
	this.ann <- ann[which((ann$canonical=="1")&(ann$geneSymbol==gene.symbol)),]
	chr <- this.ann$chrom
	start <- as.numeric(this.ann$txStart) - as.numeric(flankSize)
	end <- as.numeric(this.ann$txEnd) + as.numeric(flankSize)

	#make Granges for all other canonical genes on this chr
	ann.canon <- ann[which(ann$canonical=="1"),]
	ann.ranges <- with(ann.canon, GRanges(seqnames=chrom,ranges=IRanges(start=txStart, end=txEnd), name=name, geneSymbol=geneSymbol, strand=strand))

	test.range <- GRanges(seqnames=chr,ranges=IRanges(start=start,end=end))

	subsetByOverlaps(ann.ranges,test.range)
	

}

getSpliceVariants <- function(ann,gene.symbol)
{
	this.ann <- ann[ann$geneSymbol==gene.symbol,]
	data.frame(name=this.ann$name,txStart=this.ann$txStart,txEnd=this.ann$txEnd,cdsStart=this.ann$cdsStart,cdsEnd=this.ann$cdsEnd,exonCount=this.ann$exonCount,refmRNA=this.ann$mRNA,canonical=this.ann$canonical)
}

readRepeatMasker <- function(genome)
{
	rmsk.path <- paste(genome,".rmsk.txt",sep="")
	rmsk <- read.table.ucsc.big(rmsk.path)
	rmsk
}

getRepeats <- function(rmsk,chr,start,end)
{
	#microsatellites, dups, repeatmask percent? what to look at?
	#using data from RepeatMasker

	#subset to the right chr
	rmsk.sub <- rmsk[rmsk$genoName==chr,]

	#subset to all regions within the range
	rmsk.region <- rmsk.sub[(rmsk.sub$genoStart>=start) & (rmsk.sub$genoEnd<=end),]

	rmsk.out <- rmsk.region

	rmsk.out
}

getRepeatsIntronic <- function()
{
	#use Granges to intersect intron ranges with the repeat ranges
}

getRepeatsExonic <- function()
{

}

getRepeatsSummary <- function(repeats)
{
	ddply(repeats,.(repClass),nrow)
}

getConservation <- function()
{
	#one number to summarize between species? percent id? alignment scores?
}

getExonConservation <- function()
{
}

getIntronConservation <- function()
{
}

getMouseOrtholog <- function(ann,gene.symbol)
{
	mmBlastTab.path <- "input/hg19.mmBlastTab.txt"
	mmBlastTab <- read.table.ucsc(mmBlastTab.path)

	#take best scoring match for our geneid
	this.ann <- ann[which((ann$canonical=="1")&(ann$geneSymbol==gene.symbol)),]
	matches <- mmBlastTab[mmBlastTab$X.query==this.ann$name,]
	matches[order(matches$identity),]
	mouse_id <- matches[1,]$target

	mouse_id	
}

getGeneSymbol <- function(ann,gene.ucsc.id)
{
	this.ann <- ann[ann$name==gene.ucsc.id,]
	this.ann$geneSymbol
}

#########################################
# Usage
#########################################

#human stuff for APP
#ann <- readUCSCAnnotation("hg19")
#gene.symbol <- "APP"

#getAll(ann,gene.symbol)
#getGenePosition(ann,gene.symbol)
#getGeneCytoBands(ann,gene.symbol)
#getGeneSize(ann,gene.symbol)
#getExonNumber(ann,gene.symbol)

#introns <- getIntrons(ann,gene.symbol)
#introns[order(introns$intronSize),]

#getFlankingGenes(ann,gene.symbol,flankSize=1e6)

#rmsk <- readRepeatMasker("hg19")
#repeats <- getRepeats(rmsk,getGenePosition(ann,gene.symbol)[,1],getGenePosition(ann,gene.symbol)[,2],getGenePosition(ann,gene.symbol)[,3])
#repeats.summary <- getRepeatsSummary(repeats)

#mouse stuff for APP ortholog
#mouse.ann <- readUCSCAnnotation("mm10")
#mouse.id <- getMouseOrtholog(ann,gene.symbol)
#mouse.gene <- getGeneSymbol(mouse.ann,mouse.id)

#getAll(mouse.ann,mouse.gene)
#getGenePosition(mouse.ann,mouse.gene)
#getGeneCytoBands(mouse.ann,mouse.gene)
#getGeneSize(mouse.ann,mouse.gene)
#getExonNumber(mouse.ann,mouse.gene)

#introns <- getIntrons(mouse.ann,mouse.gene)
#introns[order(introns$intronSize),]

#getFlankingGenes(mouse.ann,mouse.gene,flankSize=1e6)


#mouse.rmsk <- readRepeatMasker("mm10")
#repeats <- getRepeats(mouse.rmsk,getGenePosition(mouse.ann,mouse.gene)[,1],getGenePosition(mouse.ann,mouse.gene)[,2],getGenePosition(mouse.ann,mouse.gene)[,3])
#repeats.summary <- getRepeatsSummary(repeats)

