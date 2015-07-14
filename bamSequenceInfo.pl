#!/usr/bin/perl

# This is a script to pull some stats from a BAM file in hopes to find which sequencing technology created it
use strict;
use warnings;
use Cwd 'abs_path';
use File::Basename;

# Set up location of BAM file
my $BAMFile = $ARGV[0]; # File path of BAM file
chomp($BAMFile);
my $BAMFileName = fileparse($BAMFile, ".bam");

my $Technology = $ARGV[1]; # Type of technology
my $ReadLengthFile = "$BAMFileName-Read-Length.txt";
my $ReadLengthStats = "$BAMFileName-Read-Length-Stats.txt";

# Set up output file
my $OutputFile = $BAMFileName . ".csv";

# Print basic script information
print "*****bamSequenceInfo.pl*****\n";
print "#####Parameters:\n";
print "\t$BAMFile\n";
print "\n";

# Generate statistics for read length of the given BAM file.  Uses Samtools and R

print "#####Generating statistics for read length\n";

system(`samtools view $BAMFile | cut -f10 | awk '{print length}' | sort > $ReadLengthFile`);

#`R -q -e "x <- read.csv('$ReadLengthFile', header = F); sink('$ReadLengthStats', append=TRUE); summary(x); sd(x[ , 1]);d<-scan('$ReadLengthFile'); png('$ReadLengthGraph'); hist(d, col='#2196f3',xlab='Read Length (bp)', main='Distribution of Read Lengths'); sink('RComplete.file');"`;
`R -q -e "x <- read.csv('$ReadLengthFile', header = F); sink('$ReadLengthStats', append=TRUE); summary(x);sd(x[ , 1]);sink('RComplete.file');"`;

sleep 2 while not -e "RComplete.file";

my $LineCount = 1;

open my $READ_LENGTH_FH, "<", $ReadLengthStats or die "Cannot open '$ReadLengthStats'\n";
open my $OUTPUT_FH, ">>", $OutputFile or die "Cannot open '$OutputFile'\n";

while (my $line = <$READ_LENGTH_FH>) {
	chomp($line);
	$line =~ s/\s//g;

	if ($LineCount == 1) {
		$LineCount += 1;
		next;
	} elsif ($LineCount == 8) {
		$line =~ s/^\[.*\]//g;
	} else {
		$line =~ s/^.*://g;
	}
	
	print $OUTPUT_FH "$line,";
	$LineCount += 1;
}
close $READ_LENGTH_FH;

print "#####Cleaning up files\n";
`rm $ReadLengthFile`;
`rm RComplete.file`;
`rm $ReadLengthStats`;

# Generate statistics for coverage bias
print "#####Generating statistics for coverage bias\n";
`picard-tools CollectGcBiasMetrics CHART_OUTPUT=$BAMFileName-chart.pdf SUMMARY_OUTPUT=$BAMFileName-summary.txt WINDOW_SIZE=100 INPUT=$BAMFile OUTPUT=$BAMFileName-GC-Bias.txt REFERENCE_SEQUENCE=hg19_random.fa`;

`tail -n +8 $BAMFileName-GC-Bias.txt | cut -f5 > $BAMFileName-Normalized-Coverage.txt`;

`R -q -e "x <- read.csv('$BAMFileName-Normalized-Coverage.txt', header = F); sink('$BAMFileName-Normalized-Coverage-stats.txt', append=TRUE); summary(x);sd(x[ , 1]);sink('RComplete2.file');"`;

sleep 2 while not -e "RComplete2.file";

$LineCount = 1;
open my $COVERAGE_BIAS_FH, "<", "$BAMFileName-Normalized-Coverage-stats.txt" or die "Cannot open '$BAMFileName-Normalized-Coverage.txt'\n";

while (my $line = <$COVERAGE_BIAS_FH>) {
	chomp($line);
        $line =~ s/\s//g;

        if ($LineCount == 1) {
                $LineCount += 1;
                next;
        } elsif ($LineCount == 8) {
                $line =~ s/^\[.*\]//g;
        } else {
                $line =~ s/^.*://g;
        }

        print $OUTPUT_FH "$line,";
        $LineCount += 1;
}

close $COVERAGE_BIAS_FH;

print "#####Cleaning up files\n";
`rm $BAMFileName-summary.txt`;
`rm $BAMFileName-Normalized-Coverage.txt`;
`rm $BAMFileName-Normalized-Coverage-stats.txt`;
`rm RComplete2.file`;
`rm $BAMFileName-GC-Bias.txt`;

# Generate Statistics for base quality score distribution
print "#####Generating statistics for base quality score distribution...\n";
`fastqc $BAMFile`;
`unzip $BAMFileName"_fastqc.zip"`;

my $FastQCOutputFile = $BAMFileName . "_fastqc/fastqc_data.txt";
my $BaseQualityDist = "$BAMFileName-Base-Quality-Dist.txt";

open my $FASTQC_FH, '<', $FastQCOutputFile or die "Can't open '$FastQCOutputFile'\n";
open my $BASE_QUALITY_FH, '>', $BaseQualityDist or die "Can't write to '$BaseQualityDist'\n";

while (<$FASTQC_FH>) {
	chomp($_);
        if (/>>Per base sequence quality/../>>END_MODULE/) {
                next if />>Per base sequence quality/ || />>END_MODULE/;
		my @Columns = split (/\t/, $_);
		if ($Columns[1] ne "Mean"){
                	print $BASE_QUALITY_FH "$Columns[1]\n";
		}
        }
}
close $BASE_QUALITY_FH;
close $FASTQC_FH;

`rm $BAMFileName-BaseQualityDist-stats.txt`;

`R -q -e "x <- read.csv('$BaseQualityDist', header = F); sink('$BAMFileName-BaseQualityDist-stats.txt', append=TRUE); summary(x);sd(x[ , 1]);sink('RComplete3.file');"`;

sleep 2 while not -e "RComplete3.file";

$LineCount = 1;
open my $BASE_QUALITY_STATS_FH, "<", "$BAMFileName-BaseQualityDist-stats.txt" or die "Cannot open '$BAMFileName-BaseQualityDist-stats.txt'\n";

while (my $line = <$BASE_QUALITY_STATS_FH>) {
       chomp($line);
        $line =~ s/\s//g;

        if ($LineCount == 1) {
                $LineCount += 1;
                next;
        } elsif ($LineCount == 8) {
                $line =~ s/^\[.*\]//g;
        } else {
                $line =~ s/^.*://g;
        }

        print $OUTPUT_FH "$line,";
        $LineCount += 1;
}

close $BASE_QUALITY_STATS_FH;

print "#####Cleaning up files\n";
`rm RComplete3.file`;
`rm $BAMFileName-BaseQualityDist-stats.txt`;
`rm $BAMFileName"_fastqc.zip"`;
`rm -r $BAMFileBam"_fastqc"`;
`rm $BAMFileName"_fastqc.html"`;
`rm $BaseQualityDist`;

print $OUTPUT_FH "$Technology\n";
close $OUTPUT_FH;
