options(echo=T)
args <- commandArgs(trailingOnly = TRUE)

### Load necessary packages ###
if(!("csaw" %in% installed.packages())) {source("http://bioconductor.org/biocLite.R"); biocLite("csaw")}
if(!("Gviz" %in% installed.packages())) {source("http://bioconductor.org/biocLite.R"); biocLite("Gviz")}
if(!("GenomicRanges" %in% installed.packages())) {source("http://bioconductor.org/biocLite.R"); biocLite("GenomicRanges")}
if(!("GenomicAlignments" %in% installed.packages())) {source("http://bioconductor.org/biocLite.R"); biocLite("GenomicAlignments")}


genes <- args[1] # The name of the file containing the genes to be plotted
genes <- unlist(read.table(genes, stringsAsFactors=F))
window <- as.numeric(args[2]) # Smoothing window size - the larger the window, the smoother the tracks
upperlim <- as.numeric(args[3]) # upper limit for y axis
upstream <- as.numeric(args[4])  # length of the upstream region to plot
downstream <- as.numeric(args[5]) # length of the downstream region to plot
pathToBam <- args[6] # path to the directory containing bam files and their indices
useNormalisedLibrarySize <- args[7] # Logical statement on whether to use normalised library sizes

bams <- list.files(pathToBam, ".bam$", full.names = T)
dir.create(paste0(pathToBam, "/TrackFiles"), showWarnings = F)

## Calculate normalisation factors for plotting - This is only done once.
if(length(list.files(pathToBam, "ChIP.library.size.csv")) == 0){
	require(csaw)
	load("repeat.masker.rda")
	param <- readParam(minq=30, pe = "none", discard=repeat.masker, dedup = T)
	binned.10kb <- windowCounts(bam.files=bams, bin=TRUE, width=10000L, param=param)
	binned.10kb
	normfacs <- normOffsets(binned.10kb)
	normfacs
	ChIP.library.size <- data.frame(LibrarySize = binned.10kb$totals, NormFactors=normfacs, EffectiveLibSize=binned.10kb$totals*normfacs)
	rownames(ChIP.library.size) <- bams
	write.csv(ChIP.library.size, file=paste0(pathToBam, "/", "ChIP.library.size.csv"))
}

require(Gviz)
require(GenomicRanges)
require(GenomicAlignments)

## Set the genome to mouse mm10 annotation
gen <- "mm10"
## Load annotation track and annotation data
load("trTrack.mm10.rda")
load("SubreadAnno.mm10.rda")
anno <- anno[ !is.na(anno$Symbol),]

# Read in the csaw-normalised library sizes and allocate them to appropriate bam files
lib.size <- read.csv(paste0(pathToBam, "/", "ChIP.library.size.csv"))
if (useNormalisedLibrarySize) {
	lib.size <- lib.size$EffectiveLibSize[match(bams, lib.size$X)] # get library size of experimental samples
} else {
	lib.size <- lib.size$LibrarySize[match(bams, lib.size$X)] # get library size of experimental samples
}
names(lib.size) <- bams
cols <- c("darkred", "darkblue", "darkgreen", "darkorange", "#ff00ff", "maroon", "forestgreen", "dodgerblue")

plot.chip.tracks <- function(target, window=window, flank.left=upstream, flank.right=downstream, upperlim=upperlim, anno.track="squish"){
  stacking(trTrack) <- anno.track
  # Get the coordinates base on the target gene
  target.start <- anno$Start[ anno$Symbol == target]
  target.end <- anno$End[ anno$Symbol == target]
  target.range <- c((min(c(target.start, target.end))), (max(c(target.start, target.end))))
  ## Set the flanking regions to desired fraction of the gene length
  target.range[1] <- target.range[1] - round(abs(target.range[1] - target.range[2])*flank.left)
  target.range[2] <- target.range[2] + round(abs(target.range[1] - target.range[2])*flank.right)
  ## Set the chromosome
  chr <- anno$Chr[ anno$Symbol == target][1]
  gtrack <- GenomeAxisTrack(fontcolor="black", col="black")
  chromosome(trTrack) <- chr
  ## Make a GRanges object to cover the region of interest
  z <- GRanges(chr,IRanges(target.range[1], target.range[2]))
  ## Read in the region of interest from all the bam files
  param <- ScanBamParam(which=z, flag=scanBamFlag(isDuplicate=F))
  x <- lapply(bams, function(x) readGAlignments(x, param=param))
  names(x) <- bams
  ## Calculate the read coverage over the region of interest
  xcov <- lapply(x, coverage)
  
  # Adjusting for library size, either raw or normalised
  xcov <- lapply(bams, function(x) xcov[[x]]/lib.size[x]*10^7)
  names(xcov) <- bams
 #  ## Calculate average coverage for each group
 #  groups <- as.factor(sapply(bams, function(x) {unlist(strsplit(split = "_", x = x))[3]}))
 #  xcov.avg <- list()
 #  for (gr in levels(groups)){
 #    print(paste("Averaging coverage for", gr))
 #    # If there are no replicates, skip averaging for that group
 #    if (sum(groups %in% gr) == 1) {xcov.avg[[gr]] <- xcov[[which(groups %in% gr)[1]]]; break}
 #    # Add up coverages for all the replicates from a given group
 #    xcov.sum <- xcov[[which(groups %in% gr)[1]]]
 #    for (i in 2:sum(groups %in% gr)){
 #      xcov.sum <- xcov.sum + xcov[[which(groups %in% gr)[i]]]
 #    }
 #    # Divide the sum coverage by the number of replicates
 #    xcov.avg[[gr]] <- xcov.sum / sum(groups %in% gr)
 #  }
 # # xcov <- lapply(bams, function(x) xcov[[x]]/xcov[[reference]])
 #  #(xcov[[1]] + xcov[[2]])
 #  #xcov[[1]] <- (xcov[[1]] + xcov[[2]])/2
 #  #xcov[[2]] <- (xcov[[3]] + xcov[[4]])/2
 #  #xcov <- xcov[1:2]
 #  grps <- names(xcov.avg) <- gsub(".bam", "", names(xcov.avg))
 #  xcov <- xcov.avg
 # 
 ## Convert the coverage object into GRanges object
 xgr <- lapply(xcov, function(x) as(x, "GRanges"))
 ## Set the colours for Control and KO panels respectively
 clr <- cols[1:length(bams)]
 names(clr) <- bams
 
 ## Generate DataTrack objects for reads
 # tracks <- lapply(bams, function(x) DataTrack(xgr[[x]][xgr[[x]] %over% z], background.title=clr[x], genome=gen, name=gsub(".bam", "", gsub(".*/", "", x)), chromosome = chr))
 tracks <- lapply(bams, function(x) DataTrack(xgr[[x]][xgr[[x]] %over% z], background.title=clr[x], genome=gen, name=gsub(".bam*", "", gsub(".*/", "", x)), chromosome = chr))

	## Plot the tracks along with the genome axis and annotation track - ideally would want to collapse the annotation track
 plotTracks(c(gtrack, tracks, trTrack), 
 		   from = target.range[1], 
 		   to = target.range[2],
 		   background.panel = "#fffff6",
 		   lwd=2, col.line="darkblue",
 		   ylim=c(0, upperlim),
 		   showId=T,
 		   geneSymbols=T,
 		   type=c("h"),
 		   col.histogram = "#f5a122",
 		   fill.histogram = "#f5a122", showTitle = T, window = -1, windowSize = window)
}

for (g in genes){
	print(paste0("Plotting ", g, " ..."))
	pdf(file=paste0(pathToBam, "/TrackFiles/", g, ".pdf"), width=8, height=6)
	plot.chip.tracks(target = g, window = window, flank.left = upstream, flank.right = downstream, upperlim = upperlim)
	dev.off()
}
