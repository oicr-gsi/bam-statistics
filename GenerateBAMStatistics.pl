#!/usr/bin/perl
# This is a script to pull some stats from a BAM file in hopes to find which sequencing technology created it
use strict;
use warnings;
use Cwd 'abs_path';
use File::Basename;
use File::Path qw(make_path remove_tree);
use POSIX qw/ceil/;

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
my $OutputDir = "$BAMFileName-stats";

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

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Useful functions
#
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Gets the location of a median of a set of numbers of size N
sub getMedianLoc {
	my $Size = $_[0];
	if ($Size % 2 == 0) {
		return ($Size + 1) / 2;
	} else {
		return ceil($Size / 2);
	}
}

# Generates a data summary file based on an input file
sub DataSummary {
	my $File = $_[0];
	my $OutputFileName = $_[1];
	my $Mean = 0;
	my $Min = 0;
	my $Max = 0;
	my $StdDev = 0;
	my $FirstQuartile = 0;
	my $FirstQuartileLower = 0;
	my $FirstQuartileUpper = 0;
	my $Median = 0;
	my $MedianLower = 0;
	my $MedianUpper = 0;
	my $ThirdQuartile = 0;
	my $ThirdQuartileLower = 0;
	my $ThirdQuartileUpper = 0;
	my $Count = 1;
	my $LineCount = `wc -l $File | cut -d " " -f1`;
	my $EvenFirstQuartile = 1;
	my $EvenThirdQuartile = 1;
	my $FirstQuartileVal = 0;
	my $FirstQuartileLowerVal = 0;
	my $FirstQuartileUpperVal = 0;
	my $MedianVal = 0;
	my $MedianLowerVal = 0;
	my $MedianUpperVal = 0;
	my $ThirdQuartileVal = 0;
	my $ThirdQuartileLowerVal = 0;
	my $ThirdQuartileUpperVal = 0;
 

	open my $INPUT_FH, '<', $File or die "Can't open '$File'\n";

	$Median = $LineCount / 2;

	if ($LineCount % 2 == 0) {
		$MedianUpper = $Median + 1;
		$MedianLower = $Median;
	
		$FirstQuartile = getMedianLoc($MedianLower - 1);
		$ThirdQuartile = getMedianLoc($LineCount - $MedianUpper) + $MedianUpper;

	} else {
		$Median = ceil($Median);
		$FirstQuartile = getMedianLoc($Median - 1);
		$ThirdQuartile = getMedianLoc($LineCount - $Median) + $Median;
	}

	if (index($FirstQuartile, ".5") != -1) {
		$FirstQuartileLower = $FirstQuartile - 0.5;
		$FirstQuartileUpper = $FirstQuartile + 0.5;
		$EvenFirstQuartile = 0;
	}

	if (index($ThirdQuartile, ".5") != -1) {
	        $ThirdQuartileLower = $ThirdQuartile - 0.5;
	        $ThirdQuartileUpper = $ThirdQuartile + 0.5;
		$EvenThirdQuartile = 0;
	}

	while (<$INPUT_FH>) {
		chomp($_);
		if ( $_ eq "" ) {
			next;
		}
		if ($Count == 1) {
			$Min = $_;
			$Max = $_;
		} elsif ($Count == $LineCount) {
			$Max = $_;
		}
	
		if ($LineCount % 2 == 0) {
			if ($Count == $MedianLower) {
				$MedianLowerVal = $_;
			} elsif ($Count == $MedianUpper) {
				$MedianUpperVal = $_;
			}
		} else {
			if ($Count == $Median) {
				$MedianVal = $_;
			}
		}
	
		if (index($FirstQuartile, ".5") != -1) { 
			if ($Count == $FirstQuartileLower)  {
				$FirstQuartileLowerVal = $_;
			} elsif ($Count == $FirstQuartileUpper) {
				$FirstQuartileUpperVal = $_;
			}
		} else {
			if ($Count == $FirstQuartile) {
				$FirstQuartileVal = $_;
			}
		}
	
		if (index($ThirdQuartile, ".5") != -1) {
			if ($Count == $ThirdQuartileLower)  {
	                        $ThirdQuartileLowerVal = $_;
	                } elsif ($Count == $ThirdQuartileUpper) {
	                        $ThirdQuartileUpperVal = $_;
	                }
		} else {
			if ($Count == $ThirdQuartile) {
				$ThirdQuartileVal = $_;
			}
		}
		
		$Mean += $_/$LineCount;
		
		$Count += 1;
	}
	
	seek $INPUT_FH, 0, 0;
	
	while (<$INPUT_FH>) {
		chomp($_);
		$StdDev += ($_ - $Mean) ** 2;
	}
	
	$StdDev = sqrt($StdDev / $LineCount);
	
	if ($LineCount % 2 == 0) {
		$MedianVal = ($MedianLowerVal + $MedianUpperVal) / 2;
	}
	
	if ($EvenFirstQuartile == 0) {
		$FirstQuartileVal = ($FirstQuartileUpperVal + $FirstQuartileLowerVal) / 2;
	}
	
	if ($EvenThirdQuartile == 0) {
	        $ThirdQuartileVal = ($ThirdQuartileUpperVal + $ThirdQuartileLowerVal) / 2;
	}
	
	close $INPUT_FH;
	
	open my $STAT_FH, '>', $OutputFileName or die "Can't write to '$OutputFileName'\n";

	print $STAT_FH "Min:$Min\n";
	print $STAT_FH "FirstQuartile:$FirstQuartileVal\n";
	print $STAT_FH "Median:$MedianVal\n";
	print $STAT_FH "Mean:$Mean\n";
	print $STAT_FH "ThirdQuartile:$ThirdQuartileVal\n";
	print $STAT_FH "Max:$Max\n";
	print $STAT_FH "StdDev:$StdDev\n";
	
	close $STAT_FH;
}


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
#`picard-tools CollectGcBiasMetrics CHART_OUTPUT=$ChartOutput SUMMARY_OUTPUT=$SummaryOutput WINDOW_SIZE=100 INPUT=$BAMFile OUTPUT=$GCOutput REFERENCE_SEQUENCE=$RefFasta`;

#`java -Xmx2g -jar /oicr/local/analysis/sw/picard/picard-tools-1.90/CollectGcBiasMetrics.jar CHART_OUTPUT=$ChartOutput SUMMARY_OUTPUT=$SummaryOutput WINDOW_SIZE=100 INPUT=$BAMFile OUTPUT=$GCOutput REFERENCE_SEQUENCE=$RefFasta`;

# Run Picard CollectInsertSizeMetrics -----------------------------------------------------------
print "Running Picard CollectInsertSizeMetrics\n";

# Set up variables
my $InsertHistogramOutput = $OutputDir . "/" . $BAMFileName . "-Insert-histogram";
my $InsertOutput = $OutputDir . "/" . $BAMFileName . "-Insert-output.txt";

# Run tool
#`picard-tools CollectInsertSizeMetrics HISTOGRAM_FILE=$InsertHistogramOutput INPUT=$BAMFile OUTPUT=$InsertOutput`;

#`java -Xmx2g -jar /oicr/local/analysis/sw/picard/picard-tools-1.90/CollectInsertSizeMetrics.jar HISTOGRAM_FILE=$InsertHistogramOutput INPUT=$BAMFile OUTPUT=$InsertOutput`;

# Run FastQC ------------------------------------------------------------------------------------
print "Running FastQC\n";

# Set up variables
my $FastQCOutputZip = $OutputDir . "/" . $BAMFileName . "_fastqc.zip";
my $FastQCFile = $OutputDir . "/" . $BAMFileName . "_fastqc/fastqc_data.txt";

# Run tool
#`/oicr/local/analysis/sw/fastqc/0.9.1/fastqc -o $OutputDir $BAMFile`;
#`rm $FastQCOutputZip`;

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

print $OUTPUT_FH "RL Min,RL 1st Quartile,RL Median,RL Mean,RL 3rd Quartile,RL Max,RL Std. Dev.,IS Min,IS 1st Quartile,IS Median,IS Mean,IS 3rd Quartile,IS Max,IS Std. Dev.,Cov Min,Cov Median,Cov Mean,Cov Max,Cov Std. Dev.,Avg Cov. below GC 10,Avg Cov above GC 75,BQ Min,BQ 1st Quartile,BQ Median,BQ Mean,BQ 3rd Quartile,BQ Max,BQ Std. Dev.,Project\n";

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Collect Read Length Stats
#
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

print "Collecting statistics on Read Length\n";

# Set up variables
my $ReadLengthOutput = $OutputDir . "/" . $BAMFileName . "-read-length.txt"; # Contains the length of all reads
my $NormReadLengthOutput = $OutputDir . "/" . $BAMFileName . "-read-length-norm.txt";
my $ReadLengthStatistics = $OutputDir . "/" . $BAMFileName . "-read-length-statistics.txt"; # Contains that statistics of Read Length
my $ReadLengthGraph = $OutputDir . "/" . $BAMFileName . "-read-length-graph.png";
my $ReadLengthAvg = 0;
my $NumOfReads = 0;
my $FirstQuartilePos = 0;
my $MedianPos = 0;
my $ThirdQuartilePos = 0; 

# Generate statistics
my $ReadLengthString = `samtools view $BAMFile | awk 'BEGIN{sum=0;count=0;} {if (\$10 != "") {len=length(\$10);sum+=len;count+=1;}} END{print sum/count " " count}'`;
($ReadLengthAvg,$NumOfReads)= split(/ /, $ReadLengthString);
$MedianPos = getMedianLoc($NumOfReads);

if ($MedianPos%2==0) {
	$FirstQuartilePos = getMedianLoc($MedianPos - 1);
	$ThirdQuartilePos = $MedianPos + $FirstQuartilePos;
} else {
	$FirstQuartilePos = getMedianLoc($MedianPos - 0.5 - 1);
	$ThirdQuartilePos = $MedianPos + $FirstQuartilePos + 0.5;
}


`samtools view $BAMFile | awk '{if (\$10 != "") print length(\$10);}' | sort -n | awk 'BEGIN{mean=$ReadLengthAvg;Lines=$NumOfReads;stddev=0;line=0;min=0;max=0;median=$MedianPos;fq=$FirstQuartilePos;tq=$ThirdQuartilePos;medianVal=0;fqVal=0;tqVal=0;} {line+=1; stddev+=(\$0-mean)**2;if (line==1) {min=\$0;max=\$0}; if (line==Lines) {max=\$0}; if (line==median) {medianVal=\$0}; if (median-0.5==line) {medianVal+=\$0} else if (median+0.5==line) {medianVal+=\$0; medianVal/=2} if (line==fq) {fqVal=\$0}; if (fq-0.5==line) {fqVal+=\$0} else if (fq+0.5==line) {fqVal+=\$0; fqVal/=2} if (line==tq) {tqVal=\$0}; if (tq-0.5==line) {tqVal+=\$0} else if (tq+0.5==line) {tqVal+=\$0; tqVal/=2}} END {print "Min:"min; print "FirstQuartile:"fqVal; print "Median:"medianVal; print "Mean:"mean; print "ThirdQuartile:"tqVal; print "Max:"max; print "StdDev:"(stddev/Lines)**0.5 }' > $ReadLengthStatistics`;

# Parse output and store to output CSV
open my $READ_LENGTH_STATS_FH, "<", $ReadLengthStatistics or die "Couldn't open '$ReadLengthStatistics'\n";

while (<$READ_LENGTH_STATS_FH>) {
	chomp($_);
	my $NormalizedValue;
	$_ =~ s/\s//g;

	$_ =~ s/^.*://g;
	$NormalizedValue = $_/$ReadLengthAvg;
        print $OUTPUT_FH "$NormalizedValue,";
        $LineCount += 1;

}
close $READ_LENGTH_STATS_FH;

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Collect Insert Size stats
#
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

print "Collecting statistics on Insert Size\n";

# Set up variables
$LineCount = 1;
my $InsertSizeData = $OutputDir . "/" . $BAMFileName . "-Insert-size.tsv";
my $TempInsertSizeData = $OutputDir . "/" . $BAMFileName . "-Insert-size-temp.tsv";
my $InsertSizeStats = $OutputDir . "/" . $BAMFileName . "-Insert-size-stats.txt";
my $InsertSizeDataExpanded = $OutputDir . "/" . $BAMFileName . "-Insert-size-exp.tsv";
my $InsertSizeAvg = 0;
$NumOfReads = 0;

# Generate statistics
`tail -n +12 $InsertOutput > $InsertSizeData`;

open my $INSERT_VALUES_FH, '<', $InsertSizeData or die "Can't read '$InsertSizeData'\n";
open my $INSERT_VALUES_EXP_FH, '>', $InsertSizeDataExpanded or die "Can't write to '$InsertSizeDataExpanded'\n";

while (<$INSERT_VALUES_FH>) {
	chomp($_);
	if (length($_) > 0) {
		my @Line = split(/\t/, $_);
		print $INSERT_VALUES_EXP_FH "$Line[0]\n" x $Line[1];
	}
}

close $INSERT_VALUES_FH;
close $INSERT_VALUES_EXP_FH;

open my $INSERT_VALUES_EXP_R_FH, '<', $InsertSizeDataExpanded or die "Could not open '$InsertSizeDataExpanded'\n";

while (<$INSERT_VALUES_EXP_R_FH>) {
	chomp($_);
	if (length($_) > 0) {
		$NumOfReads += 1;
		$InsertSizeAvg += $_;
	}
}


$InsertSizeAvg = $InsertSizeAvg / $NumOfReads;

seek $INSERT_VALUES_EXP_R_FH, 0, 0;

open my $INSERT_VALUES_NORM_FH, ">", $TempInsertSizeData or die "Can't open '$TempInsertSizeData'\n";

while (<$INSERT_VALUES_EXP_R_FH>) {
	chomp($_);
	my $Avg;
	
	if ($_ ne "") {
		$Avg = $_ / $InsertSizeAvg;
		print $INSERT_VALUES_NORM_FH "$Avg\n";
	}
}

close $INSERT_VALUES_NORM_FH;
close $INSERT_VALUES_EXP_R_FH;

`mv $TempInsertSizeData $InsertSizeDataExpanded`;
`rm $InsertSizeData`;

DataSummary($InsertSizeDataExpanded, $InsertSizeStats);

# Parse output from R and store to output CSV
open my $INSERT_SIZE_FH, "<", "$InsertSizeStats" or die "Cannot open '$InsertSizeStats'\n";

while (<$INSERT_SIZE_FH>) {
        chomp($_);
        $_ =~ s/\s//g;

        $_ =~ s/^.*://g;

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
$LineCount = 1;
my $NormalizedCoverageFile = $OutputDir . "/" . $BAMFileName . "-normalized-coverage.txt";
my $NormalizedCoverageStats = $OutputDir . "/" . $BAMFileName . "-normalized-coverage-stats.txt"; 
my $NormalizedCoverageAndGCFile = $OutputDir . "/" . $BAMFileName . "-normalized-coverage-gc.txt";
my $AvgCovBelowGC10 = 0;
my $AvgCovAboveGC75 = 0;

# Generate statistics
`tail -n +8 $GCOutput | cut -f5 | sort -n | awk '{if (length(\$0) > 0) print \$0;}'> $NormalizedCoverageFile`;
`tail -n +8 $GCOutput | cut -f1,5 | awk '{if (length(\$0) > 0) print \$0;}'> $NormalizedCoverageAndGCFile`;

DataSummary($NormalizedCoverageFile, $NormalizedCoverageStats);

# Parse output from R and store to output CSV
open my $COVERAGE_BIAS_FH, "<", "$NormalizedCoverageStats" or die "Cannot open '$NormalizedCoverageStats'\n";

while (<$COVERAGE_BIAS_FH>) {
        chomp($_);
        $_ =~ s/\s//g;

        if ($LineCount == 2 or $LineCount == 5) {
                $LineCount += 1;
                next;
        } else {
		$_ =~ s/^.*://g;
        }

        print $OUTPUT_FH "$_,";
        $LineCount += 1;
}

close $COVERAGE_BIAS_FH;

# Calculate extra stats
open my $COVERAGE_AND_GC_FH, "<", $NormalizedCoverageAndGCFile or die "Can't open '$NormalizedCoverageAndGCFile'\n";

my $Sum = 0;
$LineCount = 1;
while (<$COVERAGE_AND_GC_FH>) {
	chomp($_);
	if ($_ ne "") { 
		my @Line = split(/\t/, $_);
		if ($Line[0] <= 10) {
			$Sum += $Line[1];
		} else {
			last;
		}
	}
	$LineCount += 1;
}

$Sum = $Sum / $LineCount;

print $OUTPUT_FH  "$Sum,";

seek $COVERAGE_AND_GC_FH, 0, 0;

$Sum = 0;
$LineCount = 0;
while (<$COVERAGE_AND_GC_FH>) {
        chomp($_);
	if ($_ ne "") {
        	my @Line = split(/\t/, $_);
        	if ($Line[0] >= 75) {
        	        $Sum += $Line[1];
        	} else {
			next;
        	}
        }
	$LineCount += 1;
}

$Sum = $Sum / $LineCount;

print $OUTPUT_FH  "$Sum,";

close $COVERAGE_AND_GC_FH;

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Collect Base Quality stats
#
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
print "Collecting statistics on base qualities\n";
# Set up variables
$LineCount = 1;
my $BaseQualityDist = $OutputDir . "/" . $BAMFileName . "-Base-quality-dist.txt";
my $TempBaseQualityDist = $OutputDir . "/" . $BAMFileName . "-Base-quality-dist-temp.txt";
my $BaseQualityDistStats = $OutputDir . "/" . $BAMFileName . "-Base-quality-dist-stats.txt";

# Generate statistics
open my $FASTQC_FH, '<', $FastQCFile or die "Can't open '$FastQCFile'\n";
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

`sort -n $BaseQualityDist | awk '{if (length(\$0) > 0) print \$0;}' > $TempBaseQualityDist`;
`mv $TempBaseQualityDist $BaseQualityDist`;

# Parse output from R and store to output CSV
DataSummary($BaseQualityDist, $BaseQualityDistStats);

open my $BASE_QUALITY_STATS_FH, "<", "$BaseQualityDistStats" or die "Cannot open '$BaseQualityDistStats'\n";

while (<$BASE_QUALITY_STATS_FH>) {
       chomp($_);
        $_ =~ s/\s//g;

        $_ =~ s/^.*://g;

        print $OUTPUT_FH "$_,";
        $LineCount += 1;
}

close $BASE_QUALITY_STATS_FH;

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# General closing and cleanup
#
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
print $OUTPUT_FH "\n";
close $OUTPUT_FH;

