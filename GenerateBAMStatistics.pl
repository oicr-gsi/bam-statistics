#!/usr/bin/perl

# This is a script to pull some stats from a BAM file in hopes to find which sequencing technology created it
use strict;
use warnings;
use Cwd 'abs_path';
use File::Basename;
use File::Path qw(make_path remove_tree);

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Set up the script
#
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Set up location of BAM file and grab the filename
my $BAMFile = $ARGV[0]; # File path of BAM file
chomp($BAMFile);
my $BAMFileName = fileparse($BAMFile, ".bam");

# Set up reference fasta
my $RefFasta = $ARGV[1];
chomp($RefFasta);

# Set up output directory
my $OutputDir = "bam-stats";

if (-e $OutputDir) {
	remove_tree($OutputDir);
}

make_path($OutputDir) or die "$@";

# Set up output file
my $OutputFile = $OutputDir . "/" . $BAMFileName . ".csv";

# Print basic script information
print "*****GenerateBAMStatistics.pl*****\n";
print "#####Parameters#####\n";
print "\t$BAMFile\n";
print "\t$RefFasta\n";
print "\n";

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Run initial tools to ensure that all data is generated before analysis of data
#
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Run Picard CollectGcBiasMetrics ----------------------------------------------------------------
print "Running Picard CollectGcBiasMetrics\n";

# Set up variables
my $ChartOutput = $OutputDir . "/" . $BAMFileName . "-GC-chart.pdf";
my $SummaryOutput = $OutputDir . "/" . $BAMFileName . "-GC-summary.txt";
my $GCOutput = $OutputDir . "/" . $BAMFileName . "-GC-output.txt";

# Run tool
`picard-tools CollectGcBiasMetrics CHART_OUTPUT=$ChartOutput SUMMARY_OUTPUT=$SummaryOutput WINDOW_SIZE=100 INPUT=$BAMFile OUTPUT=$GCOutput REFERENCE_SEQUENCE=$RefFasta`;

# Run Picard CollectInsertSizeMetrics -----------------------------------------------------------
print "Running Picard CollectInsertSizeMetrics\n";

# Set up variables
my $InsertHistogramOutput = $OutputDir . "/" . $BAMFileName . "-Insert-histogram";
my $InsertOutput = $OutputDir . "/" . $BAMFileName . "-Insert-output.txt";

# Run tool
`picard-tools CollectInsertSizeMetrics HISTOGRAM_FILE=$InsertHistogramOutput INPUT=$BAMFile OUTPUT=$InsertOutput`;

# Run FastQC ------------------------------------------------------------------------------------
print "Running FastQC\n";

# Set up variables
my $FastQCOutputZip = $OutputDir . "/" . $BAMFileName . "_fastqc.zip";
my $FastQCFile = $OutputDir . "/" . $BAMFileName . "_fastqc/fastqc_data.txt";

# Run tool
#`fastqc $BAMFile`;
#`mv $FastQCOutputZip $OutputDir`;
#`unzip $FastQCOutputZip -d $OutputDir`; # Extract data

# Tools completed
print "All tools have completed!\n";

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Set up for general variables used by all stats
#
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

my $RComplete = $OutputDir . "/R-Complete.txt"; # Created once R is done running
my $LineCount = 1;

# Open output file for writing (Concatenation)
open my $OUTPUT_FH, ">>", $OutputFile or die "Cannot open '$OutputFile'\n";

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Collect Read Length Stats
#
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

print "Collecting statistics on Read Length\n";

# Set up variables
my $ReadLengthOutput = $OutputDir . "/" . $BAMFileName . "-read-length.txt"; # Contains the length of all reads
my $ReadLengthStatistics = $OutputDir . "/" . $BAMFileName . "-read-length-statistics.txt"; # Contains that statistics of Read Length
my $ReadLengthGraph = $OutputDir . "/" . $BAMFileName . "-read-length-graph.png";

# Generate statistics
`samtools view $BAMFile | cut -f10 | awk '{print length}' | sort > $ReadLengthOutput`;

`R -q -e "x <- read.csv('$ReadLengthOutput', header = F); sink('$ReadLengthStatistics', append=TRUE); summary(x); sd(x[ , 1]);d<-scan('$ReadLengthOutput'); png('$ReadLengthGraph'); hist(d, col='#2196f3',xlab='Read Length (bp)', main='Distribution of Read Lengths'); sink('$RComplete');"`;

sleep 2 while not -e $RComplete;

`rm $RComplete`;

# Parse output from R and store to output CSV
open my $READ_LENGTH_FH, "<", $ReadLengthStatistics or die "Cannot open '$ReadLengthStatistics'\n";

while (<$READ_LENGTH_FH>) {
	chomp($_);

	$_ =~ s/\s//g;

        if ($LineCount == 1) {
                $LineCount += 1;
                next;
        } elsif ($LineCount == 8) {
                $_ =~ s/^\[.*\]//g;
        } else {
                $_ =~ s/^.*://g;
        }

        print $OUTPUT_FH "$_,";
        $LineCount += 1;

}

close $READ_LENGTH_FH;

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Collect Insert Size stats
#
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

print "Collecting statistics on Insert Size\n";

# Set up variables
$LineCount = 1;
my $InsertSizeData = $OutputDir . "/" . $BAMFileName . "-Insert-size.tsv";
my $InsertSizeStats = $OutputDir . "/" . $BAMFileName . "-Insert-size-stats.txt";

# Generate statistics
`tail -n +12 $InsertOutput > $InsertSizeData`;

`R -q -e "x <- read.csv('$InsertSizeData', header = F); sink('$InsertSizeStats', append=TRUE); summary(x);sd(x[ , 1]);sink('$RComplete');"`;

sleep 2 while not -e $RComplete;

`rm $RComplete`;
exit;
# Parse output from R and store to output CSV
open my $INSERT_SIZE_FH, "<", "$InsertSizeStats" or die "Cannot open '$InsertSizeStats'\n";

while (<$INSERT_SIZE_FH>) {
        chomp($_);
        $_ =~ s/\s//g;

        if ($LineCount == 1) {
                $LineCount += 1;
                next;
        } elsif ($LineCount == 8) {
                $_ =~ s/^\[.*\]//g;
        } else {
                $_ =~ s/^.*://g;
        }

        print $OUTPUT_FH "$_,";
        $LineCount += 1;
}

close $INSERT_SIZE_FH;


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Collect Coverage bias stats
#
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

print "Collecting statistics on Coverage Bias\n";

# Set up variables
my $NormalizedCoverageFile = $OutputDir . "/" . $BAMFileName . "-normalized-coverage.txt";
my $NormalizedCoverageStats = $OutputDir . "/" . $BAMFileName . "-normalized-coverage-stats.txt"; 

# Generate statistics
`tail -n +8 $GCOutput | cut -f5 > $NormalizedCoverageFile`;

`R -q -e "x <- read.csv('$NormalizedCoverageFile', header = F); sink('$NormalizedCoverageStats', append=TRUE); summary(x);sd(x[ , 1]);sink('$RComplete');"`;

sleep 2 while not -e $RComplete;

`rm $RComplete`;

# Parse output from R and store to output CSV
open my $COVERAGE_BIAS_FH, "<", "$NormalizedCoverageStats" or die "Cannot open '$NormalizedCoverageStats'\n";

while (<$COVERAGE_BIAS_FH>) {
        chomp($_);
        $_ =~ s/\s//g;

        if ($LineCount == 1) {
                $LineCount += 1;
                next;
        } elsif ($LineCount == 8) {
                $_ =~ s/^\[.*\]//g;
        } else {
                $_ =~ s/^.*://g;
        }

        print $OUTPUT_FH "$_,";
        $LineCount += 1;
}

close $COVERAGE_BIAS_FH;

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# General closing and cleanup
#
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close $OUTPUT_FH;

