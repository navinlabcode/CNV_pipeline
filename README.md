CNV Pipeline Algorithm:
========================================================================================

* Copy number was detected from sequence read density using the variable binning method. Briefly, copy number is calculated from read density by dividing the genome into ‘bins’ and counting the number of unique reads in each interval. To determine interval sizes we simulated sequence reads by sampling 200 million sequences of length 48 from the human reference genome (HG19/NCBI37) and introduced single nucleotide errors with a frequency encountered during Illumina sequencing. These sequences were mapped back to the human reference genome using BWA and filtered for unique mappings. We assigned a number of bins to each chromosome based on the proportion of simulated reads mapped. We then divided each chromosome into bins with an equal number of simulated reads. This resulted in 12,508 genomic bins with no bins crossing chromosome boundaries. The median genomic length spanned by each bin is 220 kb. This variable binning efficiently reduces false deletion events when compared to uniform length-fixed bins. Large bins were filtered to remove false-positive amplifications in the centromeric and telomeric regions. We then applied Loess normalization to correct for GC bias. The copy number profiles were segmented using the Kolmogorov–Smirnov (KS) statistical test.
-------------

***Requirements: Bowtie, Samtools, GNU Parallel ***
-------------
 Features:

        1. it can support both hiseq and nextseq data.

        2. The segmentations for CBS can run in paralell.  

        3. it can analysis the CNV based on different bin size resolution from 100k to 10M.

        4. it generates the QC results by dashboard

-------------
Usage:

step1: specify the path to all fastq files, for example:

  	find /path/ "*.gz" -type f >fastq.input

step2: run as follows


	sh run_CNA.sh [--fastq] <fastq_path_file> [--output] <outdir_path> [--sample] <tissue sample name> [--res] <resolution of bin size> [--cpu] < Cpu number> [--undo_prune] <the proportional increase number for the CBS segment function> [alpha] <significance levels for the CBS test [filter_CellWithEmptyBin] < filter the cells with empty reads in bins > [--makeFig] <plot the CNV profile per cell>  1>process.info 2>process.log

	example:
	sh run_CNA.sh -f fastq_input -o output
	or
	sh run_CNA.sh --fastq fastq_nextseq --output output_nextseq --sample SPA --filter_CellWithEmptyBin 0.05
