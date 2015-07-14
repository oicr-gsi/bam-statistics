#!/usr/bin/perl

use strict;
use warnings;
use File::Basename;
use Term::ANSIColor qw(:constants);

# Get BAM file name
my $BAMFile = $ARGV[0];
chomp($BAMFile);
my $BAMFileName = fileparse($BAMFile, ".bam");

my $FastaFile = $ARGV[1];
chomp($FastaFile);

my $OutputFile = "$BAMFile.tsv";

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ExpandCIGAR function
# Pre: Cigar string (ex. 90M5I5M)
# Post: Expands CIGAR string 
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sub ExpandCIGAR{
	my $Input = $_[0];
	my $CigarValues = "";
	my $CigarString;
	my $CigarNum;
	my $PendingString;

	for ($Input =~ /\d+[a-zA-Z]/g) {
       		$PendingString = "";
        	$CigarNum = substr($_, 0, -1);
        	$CigarString = substr($_, -1);
        	for (my $i=0; $i < $CigarNum; $i++) {
        	        $PendingString .= $CigarString;
        	}
        	$CigarValues .= $PendingString;
	}

	return $CigarValues;
}

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Main Code Body
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Grab columns 5 and 10 from BAM file (cigar and sequence)
my $CIGAR_Sequence = $BAMFileName . "_CIGAR_Seq.tsv";

#`samtools view $BAMFile | cut -f3,4,6,10 > $CIGAR_Sequence`;
`cut -f3,4,6,10 $BAMFile > $CIGAR_Sequence`;

open my $C_S_FH, "<", $CIGAR_Sequence or die "Can't read from '$CIGAR_Sequence'\n";

my $Global_D_Count = 0;
my $Global_I_Count = 0;
my $Global_Mismatch_Count = 0;
my $D_Count = 0;
my $I_Count = 0;
my $Mismatch_Count = 0;
my $SeqOffset;
my $RefPos;
my $BaseSum = 0;
my $count = 0;

while (<$C_S_FH>) {
	chomp($_);
	$count += 1;
	$D_Count = 0;
	$I_Count = 0;
	$Mismatch_Count = 0;

	my ($Chr, $Pos, $CIGAR, $Seq) = split (/\t/, $_);
	my $CIGARChar;
	my $SeqChar;
	my $RefChar;

	for ($CIGAR =~ /\d+D/g) {
		$D_Count += 1;
	}

	$CIGAR = ExpandCIGAR($CIGAR);

	my $Length = length($CIGAR);

	$BaseSum += length($Seq);

	# Look at each position
	$SeqOffset = 0;
	$RefPos = $Pos;
	print "$count $Seq\n";
	for (my $i = 0; $i < $Length; $i++) {
		$CIGARChar = substr($CIGAR, $i, 1);
		$SeqChar = substr($Seq, $i + $SeqOffset, 1);
		
		if ($CIGARChar eq "I") {
			$I_Count += 1;
			print GREEN, "I";
			print RESET;
		} elsif ($CIGARChar eq "M" or $CIGARChar eq "=" or $CIGARChar eq "X") {
#			$Chr =~ s/chr//;
			if (index(lc $Chr, "chr") != -1) {
				$RefChar = `samtools faidx $FastaFile "$Chr:$RefPos-$RefPos" | tail -1`;
			} else {
				$RefChar = `samtools faidx $FastaFile "chr$Chr:$RefPos-$RefPos" | tail -1`;
			}
			chomp($RefChar);
			if (lc $SeqChar ne lc $RefChar) { # Different from reference
				$Mismatch_Count += 1;
				print RED, (uc $RefChar);
				print RESET;
			} else {
				print uc $RefChar;
			}
			$RefPos += 1;
		} elsif ($CIGARChar eq "D") {
			$SeqOffset -= 1;
			$RefPos += 1;
		} elsif ($CIGARChar eq "S") {
			print "S";
		} elsif ($CIGARChar eq "H") {
			$SeqOffset -= 1;
		} elsif ($CIGARChar eq "N") {
			$RefPos += 1;
		} elsif ($CIGARChar eq "P") {
			$RefPos += 1;			
		} else {
			$RefPos += 1;
		}
	} 

	print "\n\n";

#	print "There are $D_Count D blocks.\n";
#	print "There are $I_Count I's.\n";
#	print "There are $Mismatch_Count mismatches.\n\n";
	$Global_Mismatch_Count += $Mismatch_Count;
	$Global_I_Count += $I_Count;
	$Global_D_Count += $D_Count;
}
close $C_S_FH;

print "$Global_Mismatch_Count\t$Global_I_Count\t$Global_D_Count\t$BaseSum\n";
print `echo "$Global_Mismatch_Count\t$Global_I_Count\t$Global_D_Count\t$BaseSum" | awk '{print \$1/\$4*100 " " \$2/\$4*100 " " \$3/\$4*100}'`;
