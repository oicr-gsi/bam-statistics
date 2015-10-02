# bam-statistics
Scripts that generate statistics for BAM files


## bamSequenceInfo
Reads information from a BAM file and generates a TSV with summary statistics about the following:

* Read lengths, calculated with samtools
* Normalized coverage, calculated by Picard's GcBiasMetrics
* Base quality distributions using FastQC

The full output from GcBiasMetrics and FastQC is generated and preserved for convenience.

*Dependencies*

* samtools
* R
* Picard tools 1.130+
* fastqc
* Reference fasta

*Usage:*

    perl bamSequenceInfo.pl --bam <path to bam file> --ref <path to reference file> [--help for more options] 
    
    Defaults are in [square brackets]
        ref     path to genome reference fasta
        bam     path to BAM file
        samtools        path to samtools executable. [samtools]
        picard  path to picard 1.30+ JAR file. [picard.jar]
        fastqc  path to fastqc executable [fastqc]
        r-path  path to R executable [R]
        help    print this message

*Output:*

* (BAM)-GcBiasMetrics-chart.pdf - chart from Picard's CollectGcBiasMetrics
* (BAM)-GcBiasMetrics-output.txt - full output from Picard's CollectGcBiasMetrics
* (BAM).tsv - tsv with summary statistics of read lengths, normalized coverage, and base quality distribution
* (BAM)_fastq.zip - results from FastQC

Example output file:

    	Min     Q1      Median  Mean    Q3      Max     StdDev
    Read Length     	101     101     101     101     101     101     0
    Normalized Coverage     17.00   21.00   21.00   20.96   21.00   23.00   0.6468507
    Base Quality Distribution       25.75   25.89   26.71   26.91   27.71   28.94   1.096683


