#!/usr/bin/perl

# This is a script to pull some stats from a BAM file in hopes to find which sequencing technology created it
use strict;
use warnings;
use Getopt::Long;
use Cwd 'abs_path';
use File::Basename;

# Set up location of BAM file
my $BAMFile;    # File path of BAM file
my $ref;        # File path of reference
my $help;
my $samtools = "samtools";
my $picard   = "picard.jar";
my $fastqc   = "fastqc";
my $rpath    = "R";
my $OutputFile;

sub usage {
    return
"Usage: perl bamSequenceInfo.pl --bam <path to bam file> --ref <path to reference file> [--help for more options] \n\n";
}

sub options {
    print "Defaults are in [square brackets]\n";
    print "\tref	path to genome reference fasta\n";
    print "\tbam	path to BAM file\n";
    print "\tsamtools	path to samtools executable. [samtools]\n";
    print "\tpicard	path to picard 1.30+ JAR file. [picard.jar]\n";
    print "\tfastqc	path to fastqc executable [fastqc]\n";
    print "\tr-path	path to R executable [R]\n";
    print "\toutput	name or path of output file\n";
    print "\thelp	print this message\n";
}

my $result = GetOptions(
    "bam=s"      => \$BAMFile,
    "ref=s"      => \$ref,
    "samtools=s" => \$samtools,
    "picard=s"   => \$picard,
    "fastqc=s"   => \$fastqc,
    "r-path=s"   => \$rpath,
    "output=s"   => \$OutputFile,
    "help"       => \$help
) or die( "Error in command line arguments\n", usage() );

if ( defined $help ) {
    print usage();
    options();
    exit;
}

if ( not defined $BAMFile ) {
    print "BAM file must be specified\n";
    die(usage);
}
if ( not defined $ref ) {
    print "Reference file must be specified\n";
    die usage();
}

my $exitcode = system("$samtools -h >/dev/null 2>&1");
if ( $exitcode == 127 ) {
    die("samtools is not found at '$samtools'");
}
$exitcode = system("$fastqc -h >/dev/null 2>&1");
if ( $exitcode == 127 ) {
    die("fastqc is not found at '$fastqc'");
}
$exitcode = system("$rpath -h >/dev/null 2>&1");
if ( $exitcode == 127 ) {
    die("R is not found at '$rpath'");
}

###########################################################

chomp($BAMFile);
my $BAMFileName = fileparse( $BAMFile, ".bam" );

#my $Technology      = $ARGV[1];                             # Type of technology
my $ReadLengthFile  = "$BAMFileName-Read-Length.txt";
my $ReadLengthStats = "$BAMFileName-Read-Length-Stats.txt";

# Set up output file
if ( not defined $OutputFile ) {
    $OutputFile = $BAMFileName . ".csv";
}

# Print basic script information
print "*****bamSequenceInfo.pl*****\n";
print "#####Parameters:\n";
print "\t$BAMFile\n";
print "\n";

# Generate statistics for read length of the given BAM file.  Uses Samtools and R

print "#####Generating statistics for read length\n";

system(
`$samtools view $BAMFile | cut -f10 | awk '{print length}' | sort > $ReadLengthFile`
);

`R -q -e "x <- read.csv('$ReadLengthFile', header = F); sink('$ReadLengthStats', append=TRUE); summary(x);sd(x[ , 1]);sink('RComplete.file');"`;

sleep 2 while not -e "RComplete.file";

my $LineCount = 1;

open my $READ_LENGTH_FH, "<", $ReadLengthStats
  or die "Cannot open '$ReadLengthStats'\n";
open my $OUTPUT_FH, ">", $OutputFile or die "Cannot open '$OutputFile'\n";
print $OUTPUT_FH "\tMin\tQ1\tMedian\tMean\tQ3\tMax\tStdDev\n";
print $OUTPUT_FH "Read Length";
while ( my $line = <$READ_LENGTH_FH> ) {
    chomp($line);
    $line =~ s/\s//g;

    if ( $LineCount == 1 ) {
        $LineCount += 1;
        next;
    }
    elsif ( $LineCount == 8 ) {
        $line =~ s/^\[.*\]//g;
    }
    else {
        $line =~ s/^.*://g;
    }

    print $OUTPUT_FH "\t$line";
    $LineCount += 1;
}
close $READ_LENGTH_FH;

print "#####Cleaning up files\n";
`rm $ReadLengthFile`;
`rm RComplete.file`;
`rm $ReadLengthStats`;

# Generate statistics for coverage bias
print "#####Generating statistics for coverage bias\n";
`java -jar $picard CollectGcBiasMetrics CHART_OUTPUT=$BAMFileName-GcBiasMetrics-chart.pdf SUMMARY_OUTPUT=$BAMFileName-GcBiasMetrics-summary.txt WINDOW_SIZE=100 INPUT=$BAMFile OUTPUT=$BAMFileName-GcBiasMetrics-output.txt REFERENCE_SEQUENCE=$ref`;

`tail -n +8 $BAMFileName-GcBiasMetrics-output.txt | cut -f5 > $BAMFileName-Normalized-Coverage.txt`;

`R -q -e "x <- read.csv('$BAMFileName-Normalized-Coverage.txt', header = F); sink('$BAMFileName-Normalized-Coverage-stats.txt', append=TRUE); summary(x);sd(x[ , 1]);sink('RComplete2.file');"`;

sleep 2 while not -e "RComplete2.file";

$LineCount = 1;
open my $COVERAGE_BIAS_FH, "<", "$BAMFileName-Normalized-Coverage-stats.txt"
  or die "Cannot open '$BAMFileName-Normalized-Coverage.txt'\n";
print $OUTPUT_FH "\nNormalized Coverage";
while ( my $line = <$COVERAGE_BIAS_FH> ) {
    chomp($line);
    $line =~ s/\s//g;

    if ( $LineCount == 1 ) {
        $LineCount += 1;
        next;
    }
    elsif ( $LineCount == 8 ) {
        $line =~ s/^\[.*\]//g;
    }
    else {
        $line =~ s/^.*://g;
    }

    print $OUTPUT_FH "\t$line";
    $LineCount += 1;
}

close $COVERAGE_BIAS_FH;

print "#####Cleaning up files\n";
`rm $BAMFileName-GcBiasMetrics-summary.txt`;
`rm $BAMFileName-Normalized-Coverage.txt`;
`rm $BAMFileName-Normalized-Coverage-stats.txt`;
`rm RComplete2.file`;

#`rm $BAMFileName-GcBiasMetrics-output.txt`;

# Generate Statistics for base quality score distribution
print "#####Generating statistics for base quality score distribution...\n";
`$fastqc $BAMFile`;
`unzip $BAMFileName"_fastqc.zip"`;

my $FastQCOutputFile = $BAMFileName . "_fastqc/fastqc_data.txt";
my $BaseQualityDist  = "$BAMFileName-Base-Quality-Dist.txt";

open my $FASTQC_FH, '<', $FastQCOutputFile
  or die "Can't open '$FastQCOutputFile'\n";
open my $BASE_QUALITY_FH, '>', $BaseQualityDist
  or die "Can't write to '$BaseQualityDist'\n";

while (<$FASTQC_FH>) {
    chomp($_);
    if ( />>Per base sequence quality/ .. />>END_MODULE/ ) {
        next if />>Per base sequence quality/ || />>END_MODULE/;
        my @Columns = split( /\t/, $_ );
        if ( $Columns[1] ne "Mean" ) {
            print $BASE_QUALITY_FH "$Columns[1]\n";
        }
    }
}
close $BASE_QUALITY_FH;
close $FASTQC_FH;

#`rm $BAMFileName-BaseQualityDist-stats.txt`;

`R -q -e "x <- read.csv('$BaseQualityDist', header = F); sink('$BAMFileName-BaseQualityDist-stats.txt', append=TRUE); summary(x);sd(x[ , 1]);sink('RComplete3.file');"`;

sleep 2 while not -e "RComplete3.file";

$LineCount = 1;
open my $BASE_QUALITY_STATS_FH, "<", "$BAMFileName-BaseQualityDist-stats.txt"
  or die "Cannot open '$BAMFileName-BaseQualityDist-stats.txt'\n";
print $OUTPUT_FH "\nBase Quality Distribution";
while ( my $line = <$BASE_QUALITY_STATS_FH> ) {
    chomp($line);
    $line =~ s/\s//g;

    if ( $LineCount == 1 ) {
        $LineCount += 1;
        next;
    }
    elsif ( $LineCount == 8 ) {
        $line =~ s/^\[.*\]//g;
    }
    else {
        $line =~ s/^.*://g;
    }

    print $OUTPUT_FH "\t$line";
    $LineCount += 1;
}

close $BASE_QUALITY_STATS_FH;

print "#####Cleaning up files\n";
`rm RComplete3.file`;
`rm $BAMFileName-BaseQualityDist-stats.txt`;

#`rm $BAMFileName"_fastqc.zip"`;
`rm -r $BAMFileName"_fastqc"`;

#`rm $BAMFileName"_fastqc.html"`;
`rm $BaseQualityDist`;

#print $OUTPUT_FH "$Technology\n";
print $OUTPUT_FH "\n";
close $OUTPUT_FH;
