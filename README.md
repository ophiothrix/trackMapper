# trackMapper
***An executable wrapper for Bioconductor's Gviz package to plot coverage tracks for ChIP-seq experiments***

## Input variables
The script requires two input variables: the list of **bam** files and the list of **genes** to be plotted.

**bam** files are specified as the path to the directory containing the bam files. This path is specified in `MapNormTracks.param.txt` file, along with other parameters. There is currently no way to specify exact bam files - all files in the directory will be plotted. The bam files have to be position sorted and indexed. Use `samtools index` to generate index files.
**genes** are specified in `genes.txt` file, one gene symbol per line.

## Output
The script generates a pdf file for each of the genes specified in `TrackFiles` directory in the directory containing bam files. It also generates a csv file containing raw and normalised library sizes (see below) for each of the samples.

## Functionality
When the script is run for the first time, it counts the reads across the whole genome into 10kb bins. These bin counts, along with the total number of reads (library size) are used later to adjust for differences in library sizes between the samples. Additionally, the script calculates the normalisation factors based on the above bin counts. This normalisation accounts for "composition bias" as defined in the User Guide to `csaw` package. The assumption is made that on average the number of reads in these large bins is unchanged between the samples. This is suitable for experiments where overall changes in mark abundance are expected between the conditions. This "composition bias" normalisation is turned OFF by default, but can be enabled in the parameters file. By default, the coverage will be adjusted for raw library sizes.
The total library sizes, normalisation factors and normalised library sizes are saved in `ChIP.library.size.csv` in the directory containing bam files. This is only done during the first run. If you need to change anything, e.g. you add a new bam file, delete this file so it can be generated anew.

For each gene to be plotted, the script extracts all the reads overlapping the genes from the provided bam files and calculates coverage at each position over the target region. The coverage values are then normalised for library size and/or composition bias. Normalised coverages are plotted using Gviz package over the target region.

## Limitations and notes
1. The annotation files are only provided for mm10 reference genome. Therefore the script is currently only available to be used with bam files aligned to mm10 annotation. I might add hg19 annotation in the future. But for now feel free to branch the repository or get in touch if you have any specific requirements (i.e. specific reference genome).
2. There is limited space to plot the file name in the title panel. So try to make your file names as short as practical.
3. The current script plots each of the bam files separately. It is possible to average the coverages between the replicates. But since it's impossible to predict how the replicated are denoted in file names, this functionality is currently unavailable.

## Usage
1. Specify genes over which the coverage is to be plotted in `genes.txt` file.
2. Specify correct `pathToBam` argument in `MapNormTracks.param.txt` file.
3. If necessary, set other parameters as necessary, see below.

## Parameter options
The following parameters are specified inside `MapNormTracks.param.txt` file.

*genes* - string specifying the path to the file with gene symbols to be plotted. Defaults to `genes.txt` file in the script directory.

*pathToBam* - string specifying path to the directory containing bam files and their indices

*window* - integer specifying the size of smoothing window - the larger the number, the smoother the tracks.

*upperlim* - integer specifying the upper limit for y axis of the coverage tracks.

*upstream* - numeric specifying the size of the flanking region upstream of the gene, expressed as the fraction of the gene length. E.g. 0.3 - will plot the region equal to 30% of the gene length upstream of the gene. 

*downstream* - same as above but for the region downstream of the gene.

*useNormalisedLibrarySize* - logical specifying whether library size should be normalised or not.
