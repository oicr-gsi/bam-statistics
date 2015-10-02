#!/usr/bin/perl

use strict;
use warnings;
use Cwd 'abs_path';
use File::Basename;
use File::Path qw(make_path remove_tree);
use Getopt::Long;

###################################################################################################
# CoverageAndBaseQualityStats.pl
# This script pulls out coverage and quality statistics.
# Unfortunately it only works with WG as it looks at every position in the reference genome
#
# Uses:
# SAMTools
# BEDTools

########## Set up Paramaters ######################################################################

# File path of BAM file
my $BAMFile;

# File path of Fasta file
my $FAFile;

my $samtools = "samtools";
my $bedtools = "bedtools";
my $OutputDir;
my $help;

sub usage {
    return
"Usage: perl CoverageAndBaseQualityStats.pl --bam <path to bam file> --ref <path to reference file> [--help for more options] \n\n";
}

sub options {
    print "Defaults are in [square brackets]\n";
    print "\tref            path to genome reference fasta\n";
    print "\tbam            path to BAM file\n";
    print "\tsamtools       path to samtools executable. [samtools]\n";
    print "\tbedtools       path to bedtools executable. [bedtools]\n";
    print "\toutput-dir     name or path of output directory\n";
    print "\thelp           print this message\n";
}

my $result = GetOptions(
    "bam=s"        => \$BAMFile,
    "ref=s"        => \$FAFile,
    "samtools=s"   => \$samtools,
    "bedtools=s"   => \$bedtools,
    "output-dir=s" => \$OutputDir,
    "help"         => \$help
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
if ( not defined $FAFile ) {
    print "Reference file must be specified\n";
    die usage();
}

$BAMFile = abs_path($BAMFile);
$FAFile  = abs_path($FAFile);

my $exitcode = system("$samtools -h >/dev/null 2>&1");
if ( $exitcode == 127 ) {
    die("samtools is not found at '$samtools'");
}

$exitcode = system("$bedtools -h >/dev/null 2>&1");
if ( $exitcode == 127 ) {
    die("bedtools is not found at '$bedtools'");
}

####################################################################################################

######### Set up output directory and files #######################################################
my $BAMFileName = fileparse( $BAMFile, ".bam" );
if ( not defined $OutputDir ) {
    $OutputDir = $BAMFileName . "-coverage_and_base_quality";
}

# Remove old directory if it exists
if ( -e $OutputDir ) {
    remove_tree($OutputDir);
}

# Make output directory
make_path($OutputDir) or die "$@";

# BEDTools genomecov output file
my $BEDToolsOutput = "$OutputDir/BEDToolsGenomeCov.tsv";
my $GenomeCov      = "$OutputDir/genomeCov.tsv";

# Base Quality Output file
my $BaseQualityReport = "$OutputDir/BaseQualityReport.tsv";

my $AvgBQ  = 0;
my $AvgCov = 0;

####### Functions ###############################################################################

# Calls BEDTools GenomeCov Tool
# Pre: BAM file and Fasta reference file
# Post: Genome Coverage statistics
sub bedTools {
    `$bedtools genomecov -ibam $BAMFile -g $FAFile > $BEDToolsOutput`;
    `awk '{if (\$1 == "genome"){print \$0}}' $BEDToolsOutput > $GenomeCov`;
    `cat $GenomeCov`;
    `rm $BEDToolsOutput`;
}

# Determines the percent of bases at each possible quality value and stores to file
# Pre: BAM file
# Post: Base Quality Distribution file
sub baseQuality {
    my $Quality    = $ARGV[0];
    my @Qualities  = (0) x 94;
    my $TotalBases = 0;
    my $Count      = 93;      # Max Value of BQ is 94 (or 93, double check this)
    my $Perc       = 0;
    $AvgBQ = 0;

    open my $BASE_QUALITY_FH, '>', $BaseQualityReport
      or die "Can't write to '$BaseQualityReport'\n";

    foreach my $m (`$samtools view $BAMFile | cut -f11`) {
        chomp($m);

        foreach ( split //, $m ) {
            my $Qual = ord($_) - 33;
            $Qualities[$Qual] += 1;
            $TotalBases       += 1;
            $AvgBQ            += $Qual;
        }
    }

    my @ReverseQualities = reverse(@Qualities);
    foreach my $n (@ReverseQualities) {
        $Perc += $n;
        my $val = $Perc / $TotalBases * 100;
        print $BASE_QUALITY_FH "$Count\t$val\n";
        $Count--;
    }

    for ( my $i = 0 ; $i < 94 ; $i++ ) {
        $AvgBQ += $i * $Qualities[$i];
    }

    print "$AvgBQ	$TotalBases\n";
    $AvgBQ = $AvgBQ / $TotalBases;
    close $BASE_QUALITY_FH;

}

# Determines percent of bases larger than or equal to the given coverage
# Pre: Genome Coverage file from BEDTools GenomeCov, Coverage Value
# Post: Percent of Bases >= Given Coverage
sub percentBasesCovGTE {
    my $GTE        = $_[0];
    my $TotalBases = `head -1 $GenomeCov | cut -f4`;

    return
`awk 'BEGIN { sum=0 } {if(\$2 >= $GTE) {sum+=\$3/\$4;}} END {print sum*100}' $GenomeCov`;
}

# Determines percent of bases with Coverage larger than a preset set of values
# Pre: None
# Post: Prints percents to screen
sub analyzeBedTools {
    print "%Cov>0 : " . percentBasesCovGTE(1);
    print "%Cov>=8 : " . percentBasesCovGTE(8);
    print "%Cov>=30 : " . percentBasesCovGTE(30);
    print "%Cov>=50 : " . percentBasesCovGTE(50);
}

# Determines percent of bases larger than or equal to the given base quality
# Pre: Base Quality Report from baseQuality(), Base Quality Value
# Post: Percent of Bases >= Given Base Quality
sub percentBaseQualityGTE {
    my $GTE = $_[0];
    return `awk '{if(\$1 == $GTE) {print \$2}}' $BaseQualityReport`;
}

# Determines percent of bases with Base Quality larger than a preset set of values
# Pre: None
# Post: Prints percents to screen
sub analyzeBaseQuality {
    print "%BC>=0 : " . percentBaseQualityGTE(0);
    print "%BC>=20 : " . percentBaseQualityGTE(20);
    print "%BC>=30 : " . percentBaseQualityGTE(30);
}

######## Call BEDTools genomecov ################################################################
print "Collecting Coverage Statistics...\n";
bedTools();

####### Call Base Quality tool #################################################################
print "Collecting Base Quality Statistics...\n";
baseQuality();

####### Analyze BEDTools genomecov data ########################################################
print "Analyzing Coverage Statistics...\n";
analyzeBedTools();

####### Analyze Base Quality data #################################################################
print "Analyzing Base Quality Statistics...\n";
analyzeBaseQuality();

print "Script completed in ";
print time - $^T;
print " seconds\n";
