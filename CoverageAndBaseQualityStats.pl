#!/usr/bin/perl

use strict;
use warnings;
use Cwd 'abs_path';
use File::Basename;
use File::Path qw(make_path remove_tree);

###################################################################################################
# CoverageAndBaseQualityStats.pl
# This script pulls out coverage and quality statistics.
# 
# Uses:
# BioAlcidae
# BEDTools


########## Set up Paramaters ######################################################################
# File path of BAM file
my $BAMFile = $ARGV[0];
$BAMFile = abs_path( $BAMFile );

# File path of Fasta file
my $FAFile = $ARGV[1];
$FAFile = abs_path( $FAFile );

######### Set up output directory and files #######################################################
my $BAMFileName = fileparse($BAMFile, ".bam");
my $OutputDir = $BAMFileName . "-coverage_and_base_quality";

# Remove old directory if it exists
if (-e $OutputDir) {
        remove_tree($OutputDir);
}

# Make output directory
make_path($OutputDir) or die "$@";

# BEDTools genomecov output file
my $BEDToolsOutput = "$OutputDir/BEDToolsGenomeCov.tsv";
my $GenomeCov = "$OutputDir/genomeCov.tsv";

# Base Quality Output file
my $BaseQualityReport = "$OutputDir/BaseQualityReport.tsv";

######## Call BEDTools genomecov ################################################################
#bedTools();

####### Analyze BEDTools genomecov data ########################################################
analyzeBedTools();

####### Call Base Quality tool #################################################################
baseQuality();

####### Analyze Base Quality data #################################################################
analyzeBaseQuality();

####### Functions ###############################################################################
bedTools {
	`bedtools genomecov -ibam $BAMFile -g $FAFile > $BEDToolsOutput`;
	`awk '{if (\$1 == "genome"){print \$0}}' $BEDToolsOutput > $GenomeCov`;
}

baseQuality {
	# call base quality and store to $BaseQualityReport
}

percentBasesCovGTE {
	my $GTE = $_[0];
	my $TotalBases = `head -1 $GenomeCov | cut -f4`;

	return `awk 'BEGIN { sum=0 } {if(\$2 >= $GTE) {sum+=\$3;}} END {print sum/$TotalBases}' $GenomeCov`;
}

percentBasesCovLTE {
	my $LTE = $_[0];
	my $TotalBases = `head -1 $GenomeCov | cut -f4`;

        return `awk 'BEGIN { sum=0 } {if(\$2 <= $LTE) {sum+=\$3;}} END {print sum/$TotalBases}' $GenomeCov`;
}

analyzeBedTools {
	print percentBasesCovGTE(3);
	print percentBasesCovGTE(0);	
}

percentBasesQualityGTE {
	my $GTE = $_[0];
	return `awk '{if(\$1 == $GTE) {print \$2}}'`;
}

percentBasesQualityGTE {
	my $LTE = $_[0];
	return `awk '{if(\$1 == $LTE) {print \$2}}'`;
}

analyzeBaseQuality {

}
