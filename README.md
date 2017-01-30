# trackMapper
An executable wrapper for Bioconductor's Gviz package to plot coverage tracks for ChIP-seq experiments

The script takes aligned reads in the form of bam files and their indices and produces coverage tracks.

Normalisation options.
By default, the coverage is normalised by raw library sizes, this can be changed in the parameters file to use compositional bias normalisation as defined in the User Guide to csaw package - i.e. the reads are counted into 10kb bins and the assumption is made that on average the number of reads in these large bins is unchanged between the samples. 
