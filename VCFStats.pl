#!/usr/bin/perl

use strict;
use warnings;
use File::Basename;
use File::Path qw(make_path remove_tree);

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Set up the script
#
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Set up location of VCF file and grab the filename
my $VCFFile = $ARGV[0];    # File path of VCF file
chomp($VCFFile);
my $VCFFileName = fileparse( $VCFFile, ".vcf" );

# Set up location of Fasta file
my $FastaFile = $ARGV[1];
chomp($FastaFile);

my $TotalBases =
  `grep -v ">" $FastaFile | sed 's/N//g' | sed '/^\\s*\$/d' | wc -c`;
chomp($TotalBases);
my $TotalVariants;
my $AltAlleleCount        = 0;
my $CombinedAlleleCount   = 0;
my $AlleleFreqMean        = 0;
my $AlleleFreqStdDev      = 0;
my $ALLELE_COUNT          = "AC";
my $ALLELE_NUMBER         = "AN";
my $ALLELE_FREQ           = "AF";
my $VariantOccurrenceRate = 0;
my $HeterozygousVariants  = 0;
my $Count                 = 0;
my $SampleCount           = 0;

foreach (`bcftools view -H $VCFFile`) {
    my @Line     = split( /\t/, $_ );
    my @InfoLine = split( /;/,  $Line[7] );
    $Count += 1;
    $SampleCount = 0;

    foreach my $i (@InfoLine) {

        my $NumericalValue = $i;
        $NumericalValue =~ s/[^.\d]//g;

        if ( $i =~ /^$ALLELE_FREQ/ ) {
            $AlleleFreqMean += $NumericalValue;
        }
        elsif ( $i =~ /^$ALLELE_NUMBER/ ) {
            $CombinedAlleleCount += $NumericalValue;
        }
        elsif ( $i =~ /^$ALLELE_COUNT/ ) {
            $AltAlleleCount += $NumericalValue;
        }
    }

    my $ColumnCount = scalar @Line;

    for ( my $e = 9 ; $e < $ColumnCount ; $e++ ) {
        my @Sample = split( /:/, $Line[$e] );
        my $FirstAllele  = substr( $Sample[0], 0, 1 );
        my $SecondAllele = substr( $Sample[0], 2, 1 );
        if ( $FirstAllele != $SecondAllele ) {
            $HeterozygousVariants += 1;
            print "$Sample[0]\n";
        }
    }
}

$AlleleFreqMean = 1 - ( $AlleleFreqMean / $Count );
$VariantOccurrenceRate = $AltAlleleCount / $TotalBases;
$HeterozygousVariants /= $CombinedAlleleCount;

foreach (`bcftools view -H $VCFFile`) {
    my @Line     = split( /\t/, $_ );
    my @InfoLine = split( /;/,  $Line[7] );

    foreach my $i (@InfoLine) {
        my $NumericalValue = $i;
        $NumericalValue =~ s/[^.\d]//g;

        if ( $i =~ /^$ALLELE_FREQ/ ) {
            $AlleleFreqStdDev +=
              ( ( 1 - $NumericalValue ) - $AlleleFreqMean )**2;
        }
    }

}

$AlleleFreqStdDev = sqrt( $AlleleFreqStdDev / $Count );

print "Allele Balance Mean = $AlleleFreqMean\n";
print "Allele Balance Std. Deviation = $AlleleFreqStdDev\n";
print "Average Variant Occurence Rate = $VariantOccurrenceRate\n";
print "Ratio of Heterozygous Variants = $HeterozygousVariants\n";
