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
# SAMTools
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
#if (-e $OutputDir) {
#        remove_tree($OutputDir);
#}

# Make output directory
#make_path($OutputDir) or die "$@";

# BEDTools genomecov output file
my $BEDToolsOutput = "$OutputDir/BEDToolsGenomeCov.tsv";
my $GenomeCov = "$OutputDir/genomeCov.tsv";

# Base Quality Output file
my $BaseQualityReport = "$OutputDir/BaseQualityReport.tsv";

my $AvgBQ = 0;
my $AvgCov = 0;

####### Functions ###############################################################################
sub bedTools {
	`bedtools genomecov -ibam $BAMFile -g $FAFile > $BEDToolsOutput`;
	`awk '{if (\$1 == "genome"){print \$0}}' $BEDToolsOutput > $GenomeCov`;
	`cat $GenomeCov`;
	`rm $BEDToolsOutput`;
}

sub baseQuality {
	my $Quality = $ARGV[0];
	my @Qualities = (0)x94;
	my $TotalBases = 0;
	my $Count = 93;
	my $Perc = 0;
	$AvgBQ = 0;	

	open my $BASE_QUALITY_FH, '>', $BaseQualityReport or die "Can't write to '$BaseQualityReport'\n";

	foreach my $m (`samtools view $BAMFile | cut -f11`) {
		chomp($m);
		
		foreach (split //, $m) {
			my $Qual = ord($_) - 33;
			$Qualities[$Qual] += 1;
			$TotalBases += 1;
			$AvgBQ += $Qual;
		}
	}
	
	my @ReverseQualities = reverse(@Qualities);
	foreach my $n (@ReverseQualities) {
        	$Perc += $n;
        	my $val = $Perc / $TotalBases * 100;
        	print $BASE_QUALITY_FH "$Count\t$val\n";
        	$Count--;
	}
	
	for (my $i = 0; $i < 94; $i++) {
		$AvgBQ += $i * $Qualities[$i];
	}

	print "$AvgBQ	$TotalBases\n";
	$AvgBQ = $AvgBQ / $TotalBases;
	close $BASE_QUALITY_FH;

}

sub percentBasesCovGTE {
	my $GTE = $_[0];
	my $TotalBases = `head -1 $GenomeCov | cut -f4`;

	return `awk 'BEGIN { sum=0 } {if(\$2 >= $GTE) {sum+=\$3/\$4;}} END {print sum*100}' $GenomeCov`;
}

sub percentBasesCovLTE {
	my $LTE = $_[0];
	my $TotalBases = `head -1 $GenomeCov | cut -f4`;

        return `awk 'BEGIN { sum=0 } {if(\$2 <= $LTE) {sum+=\$3/\$4;}} END {print sum*100}' $GenomeCov`;
}

sub analyzeBedTools {
	print "%Cov>0 : " . percentBasesCovGTE(1);	
	print "%Cov>=8 : " . percentBasesCovGTE(8);
	print "%Cov>=30 : " . percentBasesCovGTE(30);
	print "%Cov>=50 : " . percentBasesCovGTE(50);
}

sub percentBaseQualityGTE {
	my $GTE = $_[0];
	return `awk '{if(\$1 == $GTE) {print \$2}}' $BaseQualityReport`;
}

sub percentBaseQualityLTE {
	my $LTE = $_[0];
	return `awk '{if(\$1 == $LTE) {print \$2}}' $BaseQualityReport`;
}

sub analyzeBaseQuality {
	print "%BC>=0 : " . percentBaseQualityGTE(0);
	print "%BC>=20 : " . percentBaseQualityGTE(20);
	print "%BC>=30 : " . percentBaseQualityGTE(30);
	print "%BC>=" . int($AvgBQ) . " : " . percentBaseQualityGTE(int($AvgBQ));
}

######## Call BEDTools genomecov ################################################################
print "Collecting Coverage Statistics...\n";
bedTools();

####### Analyze BEDTools genomecov data ########################################################
print "Analyzing Coverage Statistics...\n";
analyzeBedTools();

####### Call Base Quality tool #################################################################
print "Collecting Base Quality Statistics...\n";
baseQuality();

####### Analyze Base Quality data #################################################################
print "Analyzing Base Quality Statistics...\n";
analyzeBaseQuality();

print "Script completed in ";
print time - $^T;
print " seconds\n";
