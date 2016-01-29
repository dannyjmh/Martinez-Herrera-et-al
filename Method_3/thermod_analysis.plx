#!/usr/bin/perl
use strict;
use warnings;
use Tie::File::AsHash;
use List::MoreUtils qw(uniq);
require 'library.pl';
$| = 1;

# Launches RNAcofold with the seedmatches and constraints annotated in the 'position' file

unless($#ARGV == 4){		# Print help
	print <DATA>;
	exit 0
}

print "Do the thermodynamic analysis using the whole cDNA sequences (press)"


foreach my $i(0..2){												
	transc_one_line($ARGV[$i]);	
}

tie my %mirnas, 'Tie::File::AsHash', "temp_$ARGV[0]", split => '\t' or die "$!: $ARGV[0]\n";
shift(@ARGV);

tie my %cdnas, 'Tie::File::AsHash', "temp_$ARGV[0]", split => '\t'  or die "$!: $ARGV[0]\n";
shift(@ARGV);

tie my %threes, 'Tie::File::AsHash', "temp_$ARGV[0]", split => '\t'  or die "$!: $ARGV[0]\n";
shift(@ARGV);

my $coords = get_coords($ARGV[0]);										# Find the 3'UTR location in the cDNA
shift(@ARGV);

my ($seed,$target,$pos,$mirna,$start_seed,$struc,%calculations);
open my $posit,'<',$ARGV[0] or die "$!: $ARGV[0]";						# Read the 'position'-like file
while(<$posit>){
	chomp;
	if(/^>/){
		(undef,$seed,$target,$pos) = split /\t/;
		$pos += $coords->{$target} if $coords->{$target};
	}
	else{
		my @mirdata = split(/;/);
		@mirdata = uniq(@mirdata);
		foreach my $mirdata(@mirdata){
			($mirna,$start_seed,$struc) = split(/,/,$mirdata);
		}
		
		$calculations{"${seed};${mirna};${target};${start_seed};${pos};$struc"} = $.;
	}	
	
}
close $posit or die "$!: $ARGV[0]";


send_2ry(\%calculations,\%cdnas,\%threes,\%mirnas);						# Launch RNAcofold using cDNA sequences
send_2ry_70(\%calculations,\%cdnas,\%threes,\%mirnas);					# Launch RNAcofold using seedmatch and 70nt-long flanking sequences


__DATA__

	Launches RNAcofold with the seedmatches and constraints annotated in the 'position' file
	
Usage: ./thermod_analysis.plx <miRNAs> <cDNAs secuences> <3'UTR sequences> <3'UTRs start> <'position'-like file>


Intructions of use
------------------
1- We strongly recommend using running the script input_nspks.plx with the data prior to using thermod_analysis.plx

2- The subroutines file 'library.pl' must be located either on this folder or any of those contained in @INC.
   Alternatively, its path can be specified in line number 4,

	 e.g. "require /<my_path>/library.pl"

3- The input sequence files (microRNAs, cDNAs and 3'UTRs) must be in FASTA format. There must not be "LRG" nor "Sequences unavailable" entries.

4- This script makes use of all available CPU cores. The user can specify this number by adjusting the value of the variable '$cores' of the subroutines file 'library.pl' (lines number 1079 and 1144).
   Default value: 5.
   
5- The <3'UTRs start> file must be a tab-delimited, two-column file containing the start position of the 3'UTRs on the cDNA sequence of the transcript,
   like the output file of the script 'get_coordinates.plx', e.g.
		Gene|Transcript	3UTR start
		ENSDARG00000000001|ENSDART00000000004	1380
		ENSDARG00000000018|ENSDART00000000019	1643
		ENSDARG00000000018|ENSDART00000138183	1587
		ENSDARG00000000019|ENSDART00000124452	671
		ENSDARG00000000068|ENSDART00000000069	1293

6- Last argument to the script must be a file like the output of the subroutine &enrich (file 'position'):

	><\tab><seedmatch><\tab><target transcript><\tab><position of the site in the 3'UTR>
	<miRNA n>,<first nt of the seed,0-based count>;....;<miRNA n>,<first nt of the seed,0-based count>

	e.g.
	
	>	CTCCATTAC	ENSMUSG00000006378_ENSMUST00000006544	833
	mmu-miR-1951,0;mmu-miR-6988-5p,0
	
7- The output of this script are two types of files, named:
	a) 'constraint<number>', e.g. 'constraint2'. These files contain the sequences of the miRNA and the target sequence as well as the hybrid structure in dot-bracket notation,
		e.g. 

			GCCCCUGGGCCUAUCCUAGAA&CCCTAGCCAGTGCCGTCCCCACGGTGAGCTGGCAGAACAGGGGAACCACCGAACTAAAGCCCAATTTGCTTGCATTTTCTGATGA
			.((((((..............&.....................................))))))..........................................
	b) 'partners<number>', e.g. 'partners2'. These are 5-column, tab-delimited files that contain for each interaction in that order: 
		seedmatch, miRNA, target, position of the seedmatch in the cDNA and the seedmatch type.
		The seedmatch type can be of four groups: 'w'(seedmatch including wobbles), 'm'(seedmatch including mismatches), 'wm'(seedmatch including both wobbles and mismatches)
		and 'nowm' (seedmatch that does not include either wobbles nor mismatches).

8- By default, this script calculates constrained binding energies with RNAcofold, first using the whole target cDNA (line number ) and secondly
 by creating a target sequence composed of the seedmatch and flanking sequences of 70 nucleotides (line number ). Comment any of those lines to supress th calculation.

Best regards,
Danny(dmherrera@unav.es)


