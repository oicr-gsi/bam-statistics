#!/usr/bin/perl

# This is a script to calculate read depth of every base over target boundaries 
# provided by two different BED files, one upstream and one downstream
# Note: This script uses R! 
#
# USAGE:
#
# perl readDepthOverTargetBoundary.pl /PATH/TO/YOUR/BAMFILE.bam  \
#	/PATH/TO/YOUR/UPSTREAM_BEDFILE.bed \
# 	/PATH/TO/YOUR/DOWNSTREAM_BEDFILE.bed
#
#

use strict;
use warnings;
use Cwd 'abs_path';
use File::Basename;
use List::MoreUtils qw(pairwise minmax);
#use Math::Round qw(nearest);

# Add your Samtools path here!!! Must be version 1.3 or higher
my $SAMTOOLS = "/u/sbundhoo/tools/samtools/samtools";

# Set up location of BAM file
my $BAMFile = $ARGV[0]; # File path of BAM file
chomp($BAMFile);
my $BAMFileName = fileparse($BAMFile, ".bam");

# Read BED files  
my $UPSTREAMBEDFile = $ARGV[1];
my $DOWNSTREAMBEDFile = $ARGV[2];


my $ReadDepthFileUpstream = "$BAMFileName-Read-Depth-Upstream.txt";
my $ReadDepthFileDownstream = "$BAMFileName-Read-Depth-Downstream.txt";
#my $UpstreamLocation = "$BAMFileName-Upstream-Location.txt";
my $UpstreamIndices = "$BAMFileName-Upstream-Indices.txt";
my $DownstreamIndices = "$BAMFileName-Downstream-Indices.txt";

# Set up output file
#my $OutputFile = $BAMFileName . ".csv";

# Print basic script information
print "*****readDepth.pl*****\n";
print "#####Parameters:\n";
print "\t$BAMFile\n";
print "\t$UPSTREAMBEDFile\n";
print "\t$DOWNSTREAMBEDFile\n";
print "\n";


#### Getting Read Depth of Bases over Upstream Region Intervals ####
## Generating Read Depth 
print "Calculating read depth for upstream regions\n";
system(`$SAMTOOLS depth -a -b $UPSTREAMBEDFile $BAMFile | cut -f2,3 > $ReadDepthFileUpstream`);
print "Samtools Depth - Done\n";
print "Getting Locations\n";

## Get upstream location
`cat $ReadDepthFileUpstream | cut -f1 > $BAMFileName-UpstreamLocation.txt`;
#`R -q -e "x <- read.csv('$BAMFileName-Upstream-Location.txt', header = F); sink('$BAMFileName-Upstream-Location-stats.txt', append=TRUE); summary(x);sd(x[ , 1]);sink('RComplete2.file');"`;
# sleep 2 while not -e "RComplete2.file";
print "Changing locations to indices for upstream boundary regions\n";

## Change locations into indices 
open my $UPSTREAM_LOCATION_FH, '<', "$BAMFileName-UpstreamLocation.txt" or die "Cannot open '$BAMFileName-UpstreamLocation.txt'\n";
open my $OUTPUT_FH, '>>', "$BAMFileName-UpstreamIndices.txt" or die "Cannot write to '$BAMFileName-UpstreamIndices.txt'\n";

my $LineCount = 0;
my $count = 999; #length of target intervals 

while (my $line = <$UPSTREAM_LOCATION_FH>) {
	if ($LineCount == 0) {
		print $OUTPUT_FH "$count\n";
		$LineCount = int($line);
		$count -= 1;
	}
	else {
		my $val = int($line);
		my $new_val = $val - 1;
		if ($LineCount == $new_val) {
			print $OUTPUT_FH "$count\n";
			$count -= 1;
			$LineCount = int($line);
		}
		else {
			print $OUTPUT_FH "999\n";
			$count = 998;
			$LineCount = int($line);
		}
	}
}

close $OUTPUT_FH;
close $UPSTREAM_LOCATION_FH;

`cat $ReadDepthFileUpstream | cut -f2 > $BAMFileName-UpstreamReadDepth.txt`;
`paste $BAMFileName-UpstreamIndices.txt $BAMFileName-UpstreamReadDepth.txt > $BAMFileName-Upstream.txt`;
print "Read depth for upstream regions is complete\n";

##CLEAN UP files
print "Cleaning up\n";
`rm $BAMFileName-UpstreamIndices.txt $BAMFileName-UpstreamLocation.txt`;
`rm $BAMFileName-Read-Depth-Upstream.txt $BAMFileName-UpstreamReadDepth.txt`;


################################################################################

#### Getting Read Depth of Bases over Downstream Region Intervals ####
## Generating Read Depth 

print "Calculating read depth for downstream regions\n";
system(`$SAMTOOLS depth -a -b $DOWNSTREAMBEDFile $BAMFile | cut -f2,3 > $ReadDepthFileDownstream`);
print "Samtools Depth - Done\n";
print "Getting locations\n";

## Get upstream location
`cat $ReadDepthFileDownstream | cut -f1 > $BAMFileName-DownstreamLocation.txt`;

print "Changing locations to indices for downstream boundary regions\n";

## Change locations into indices 
open my $DOWNSTREAM_LOCATION_FH, '<', "$BAMFileName-DownstreamLocation.txt" or die "Cannot open '$BAMFileName-DownstreamLocation.txt'\n";
open my $D_OUTPUT_FH, '>>', "$BAMFileName-DownstreamIndices.txt" or die "Cannot write to '$BAMFileName-DownstreamIndices.txt'\n";

my $D_LineCount = 0;
my $D_count = 0; #length of target intervals 

while (my $line = <$DOWNSTREAM_LOCATION_FH>) {
	if ($D_LineCount == 0) {
		print $D_OUTPUT_FH "$D_count\n";
		$D_LineCount = int($line);
		$D_count += 1;
	}
	else {
		my $d_val = int($line);
		my $d_new_val = $d_val - 1;
		if ($D_LineCount == $d_new_val) {
			print $D_OUTPUT_FH "$D_count\n";
			$D_count += 1;
			$D_LineCount = int($line);
		}
		else {
			print $D_OUTPUT_FH "0\n";
			$D_count = 1;
			$D_LineCount = int($line);
		}
	}
}

close $D_OUTPUT_FH;
close $DOWNSTREAM_LOCATION_FH;

`cat $ReadDepthFileDownstream | cut -f2 > $BAMFileName-DownstreamReadDepth.txt`;
`paste $BAMFileName-DownstreamIndices.txt $BAMFileName-DownstreamReadDepth.txt > $BAMFileName-Downstream.txt`;
print "Read depth for downstream regions is complete\n";

##CLEAN UP files
print "Cleaning up\n";
`rm $BAMFileName-DownstreamIndices.txt $BAMFileName-DownstreamLocation.txt`;
`rm $BAMFileName-Read-Depth-Downstream.txt $BAMFileName-DownstreamReadDepth.txt`;

#### Merging the two files $BAMFileName-Upstream.txt and $BAMFileName-Downstream.txt ####
print "Merging started\n";
`cat $BAMFileName-Upstream.txt $BAMFileName-Downstream.txt > $BAMFileName-Complete.txt`;
#`sort -nk1 $BAMFileName-Complete.txt > $BAMFileName-Complete-Sorted.txt`;

my @hdr;
my %sums = ();
my %count = ();
my $key;

open my $MERGED_FH, '<', "$BAMFileName-Complete.txt" or die "Cannot open '$BAMFileName-Complete.txt'\n";
open my $F_OUTPUT_FH, '>>', "$BAMFileName-Merged.csv" or die "Cannot write to '$BAMFileName-Merged.csv'\n";


while (defined($_ = <$MERGED_FH>)) {
  chomp $_;
  my @F = split(' ', $_, 0);

  # UGLY: hardcoded to expect exactly 1 header row
  if ($. == 1) {
    @hdr = @F;
    next;
  }
  # sum column-wise, grouped by first column
  $key = shift @F;
  if ( exists $sums{$key} ) {
    $sums{$key} = [ pairwise { $a + $b } @{ $sums{$key} }, @F];
  }
  else {
    $sums{$key} = \@F;
  }

  $count{$key}++;
}

my %avgs = ();
# find the column averages, and the global max of those averages
for $key ( keys %sums ) {
  $avgs{$key} = [ map { $_ / $count{$key} } @{ $sums{$key} } ];
}
# print the results
print $F_OUTPUT_FH "indices,avg\n";
for $key ( sort keys %avgs ) {
  print $F_OUTPUT_FH join ",", $key, @{ $avgs{$key} };
  print $F_OUTPUT_FH "\n";
}

close $F_OUTPUT_FH;
close $MERGED_FH;

#`sed -i 's/\t/,/g' $BAMFileName-Merged.txt`;
#`echo 'distance,Dataset_1,Dataset_2,Dataset_3,Dataset_4,Dataset_5,Dataset_6' > header.csv`;
#`cat header.csv $BAMFileName-Merged.txt > $BAMFileName.csv`;

print "Merging completed\n";

## CLEAN UP
print "Cleaning up\n";
`rm $BAMFileName-Upstream.txt $BAMFileName-Downstream.txt`;
`rm $BAMFileName-Complete.txt`;

## $BAMFileName.csv can be used to generate scatter plot in R
#### Plotting the curve in R workspace ####

#`R "png('Horizon_GM24143_Ly_C_PE_301_EX_chr21.png'); 
#dat <- read.csv('Horizon_GM24143_Ly_C_PE_301_EX_chr21-Merged.csv', header = T); 
#ind <- dat$indices; 
#avg <- dat$avg; 
#g <- plot(ind, avg, main="Read Depth Across Target Boundary Regions", ylab="Average Read Depth", xlab="distance (bp) from target boundary");
#dev.off();"`;
