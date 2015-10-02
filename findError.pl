#!/usr/bin/perl

use strict;
use warnings;
use File::Basename;
use Term::ANSIColor qw(:constants);
use Getopt::Long;

# Get BAM file name
my $BAMFile;
my $FastaFile;
my $OutputFile;
my $help;
my $samtools="samtools";

sub usage {
    return
"Usage: perl findErrors.pl --bam <path to bam file> --ref <path to reference file> [--help for more options] \n\n";
}

sub options {
    print "Defaults are in [square brackets]\n";
    print "\tref            path to genome reference fasta\n";
    print "\tbam            path to BAM file\n";
    print "\tsamtools       path to samtools executable. [samtools]\n";
    print "\toutput     name or path of output file\n";
    print "\thelp           print this message\n";
}

my $result = GetOptions(
    "bam=s"        => \$BAMFile,
    "ref=s"        => \$FastaFile,
    "samtools=s"   => \$samtools,
    "output=s" => \$OutputFile,
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
if ( not defined $FastaFile ) {
    print "Reference file must be specified\n";
    die usage();
}

my $exitcode = system("$samtools -h >/dev/null 2>&1");
if ( $exitcode == 127 ) {
    die("samtools is not found at '$samtools'");
}

chomp($FastaFile);
chomp($BAMFile);
my $BAMFileName = fileparse( $BAMFile, ".bam" );

if (not defined $OutputFile) {
$OutputFile = "$BAMFile.tsv";

}



# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ExpandCIGAR function
# Pre: Cigar string (ex. 90M5I5M)
# Post: Expands CIGAR string
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sub ExpandCIGAR {
    my $Input       = $_[0];
    my $CigarValues = "";
    my $CigarString;
    my $CigarNum;
    my $PendingString;

    for ( $Input =~ /\d+[a-zA-Z]/g ) {
        $PendingString = "";
        $CigarNum      = substr( $_, 0, -1 );
        $CigarString   = substr( $_, -1 );
        for ( my $i = 0 ; $i < $CigarNum ; $i++ ) {
            $PendingString .= $CigarString;
        }
        $CigarValues .= $PendingString;
    }

    return $CigarValues;
}

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ASCIIToQual function
# Pre: ASCII value
# Post: Base Quality
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sub ASCIIToQual {
    my $Input   = $_[0];
    my $Quality = ord($Input) - 33;
    my $Value   = 10**( -$Quality / 10 );
    return $Value;
}

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Main Code Body
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Grab columns 5 and 10 from BAM file (cigar and sequence)
my $CIGAR_Sequence = $BAMFileName . "_CIGAR_Seq.tsv";

#`$samtools view $BAMFile | cut -f3,4,6,10 > $CIGAR_Sequence`;
`cut -f3,4,6,10,11 $BAMFile > $CIGAR_Sequence`;

open my $C_S_FH, "<", $CIGAR_Sequence
  or die "Can't read from '$CIGAR_Sequence'\n";

my $Global_D_Count        = 0;
my $Global_I_Count        = 0;
my $Global_Mismatch_Count = 0;
my $D_Count               = 0;
my $I_Count               = 0;
my $Mismatch_Count        = 0;
my $SeqOffset;
my $RefPos;
my $BaseSum            = 0;
my $count              = 0;
my $ExpectedErrorRate  = 0;
my $ReadErrorRate      = 0;
my $ChrPosFile         = "chr_pos.tsv";
my $CommonVariantCount = 0;

while (<$C_S_FH>) {
    chomp($_);
    $count += 1;
    $D_Count        = 0;
    $I_Count        = 0;
    $Mismatch_Count = 0;

    my ( $Chr, $Pos, $CIGAR, $Seq, $Qual ) = split( /\t/, $_ );
    my $CIGARChar;
    my $SeqChar;
    my $RefChar;

    for ( $CIGAR =~ /\d+D/g ) {
        $D_Count += 1;
    }

    $CIGAR = ExpandCIGAR($CIGAR);

    my $Length = length($CIGAR);

    $BaseSum += length($Seq);

    # Look at each position
    $SeqOffset = 0;
    $RefPos    = $Pos;

    #	print "$count $Seq\n";
    #	print "$count ";
    for ( my $i = 0 ; $i < $Length ; $i++ ) {
        $CIGARChar = substr( $CIGAR, $i, 1 );
        $SeqChar = substr( $Seq, $i + $SeqOffset, 1 );

        if ( $CIGARChar eq "I" ) {
            $I_Count += 1;

            #			print GREEN, "I";
            #			print RESET;
        }
        elsif ( $CIGARChar eq "M" or $CIGARChar eq "=" or $CIGARChar eq "X" ) {

            #			$Chr =~ s/chr//;
            if ( index( lc $Chr, "chr" ) != -1 ) {
                $RefChar =
                  `$samtools faidx $FastaFile "$Chr:$RefPos-$RefPos" | tail -1`;
            }
            else {
                $RefChar =
`$samtools faidx $FastaFile "chr$Chr:$RefPos-$RefPos" | tail -1`;
            }
            chomp($RefChar);
            if ( lc $SeqChar ne lc $RefChar ) {    # Different from reference
                if ( `grep "$Chr	$RefPos" $Chr"_"$ChrPosFile | wc -l` == 0 ) {
                    $Mismatch_Count += 1;
                }
                else {
                    $CommonVariantCount += 1;
                }

                #				print RED, (uc $RefChar);
                #				print RESET;
            }
            else {
                #				print uc $RefChar;
            }
            $RefPos += 1;
        }
        elsif ( $CIGARChar eq "D" ) {
            $SeqOffset -= 1;
            $RefPos += 1;
        }
        elsif ( $CIGARChar eq "S" ) {

            #			print "S";
        }
        elsif ( $CIGARChar eq "H" ) {
            $SeqOffset -= 1;
        }
        elsif ( $CIGARChar eq "N" ) {
            $RefPos += 1;
        }
        elsif ( $CIGARChar eq "P" ) {
            $RefPos += 1;
        }
        else {
            $RefPos += 1;
        }
    }

    #	print "\n\n";

    # Find expected error rate
    $ReadErrorRate = 0;
    for ( my $i = 0 ; $i < length($Qual) ; $i++ ) {
        my $QualChar = substr( $Qual, $i, 1 );
        $ReadErrorRate += ASCIIToQual($QualChar);

    }
    $ExpectedErrorRate += $ReadErrorRate;

    #	print "There are $D_Count D blocks.\n";
    #	print "There are $I_Count I's.\n";
    #	print "There are $Mismatch_Count mismatches.\n\n";
    $Global_Mismatch_Count += $Mismatch_Count;
    $Global_I_Count        += $I_Count;
    $Global_D_Count        += $D_Count;
}
close $C_S_FH;

$ExpectedErrorRate = $ExpectedErrorRate / $BaseSum;

print "$Global_Mismatch_Count\t$Global_I_Count\t$Global_D_Count\t$BaseSum\n";
print
`echo "$Global_Mismatch_Count\t$Global_I_Count\t$Global_D_Count\t$BaseSum" | awk '{print \$1/\$4 " " \$2/\$4 " " \$3/\$4}'`;
print "Expected error rate: $ExpectedErrorRate (Expected # of errors : "
  . ( $ExpectedErrorRate * $BaseSum ) . ")\n";
print "Common variants $CommonVariantCount\n";
