#!/usr/bin/perl
use strict;
use warnings;
require 'library.pl';
use Math::Combinatorics;
use Math::Cartesian::Product;
use Tie::File;
use Tie::File::AsHash;
use POSIX;
use Array::Utils qw(array_minus);
use List::Util qw(sum);
use List::MoreUtils qw(uniq);
use POSIX qw(isdigit);
$| = 1;

# Script to calculate the enrichment of miRNAs and their seeds in a group of sequences as a whole, e.g.: all 3'UTRs of a species


unless($#ARGV == 2){		# Print help
		print <DATA>;
		exit 0;
}


my $mir_ver = $1 if $ARGV[0] =~ /(miRB\d+)/;
my $organism = (split /_/,$ARGV[1])[0];
my $ens_ver = $1 if $ARGV[1] =~ /(Ens\d+)/;

transc_one_line($ARGV[0]);
background_one_line($ARGV[1]);

tie my @background, 'Tie::File', "temp_$ARGV[1]" or die "$!: $ARGV[1]\n";
pop(@ARGV);

tie my %mirnas, 'Tie::File::AsHash', "temp_$ARGV[-1]", split => '\t' or die "$!: $ARGV[-1]\n";
pop(@ARGV);

foreach my $file(<temp_*>){unlink($file)};

print "Extracting the seedmatches....";
my ($seeds,$seed_type,$compiled_seeds) = get_seedmatches(\%mirnas);			# Create the seedmatches of every miRNA
print "DONE!\n";

my $tot_length = sum( map { length($_) } @background);								# Total length of input sequences
print "Calculating miRNAs' and seeds' enrichment....";

my $background_nspks = enrich_background(\%mirnas,\@background,$seeds,$compiled_seeds,$tot_length);	# Calculate miRNAs' and seedmatches' nspk values
print "DONE!\n";

untie @background or die "$!: background\n";
untie %mirnas or die "$!: miRNAs\n";

my $outfile = ($organism && $mir_ver && $ens_ver) ? "enrichment_${organism}_${mir_ver}_$ens_ver" : "enrichment_${ARGV[0]}-$ARGV[1]";
	
open my $res, ">", "enrichment_${organism}_${mir_ver}_${ens_ver}_6-8" or die "$!: enrichment_${organism}_${mir_ver}_$ens_ver";
print $res "miRNA|seed\tNSPKg\n";
foreach my $mir_seeds(sort keys %$background_nspks){
	print $res "$mir_seeds\t$background_nspks->{$mir_seeds}\n"	

}	
close $res or die "$!: enrichment_${organism}_${mir_ver}_$ens_ver";
print "Assessed miRNAs' and seedmatches' enrichment in ",$organism," usando ",$mir_ver," y ",$ens_ver,"\n";


__DATA__

	Calculates the enrichment of miRNAs and their seeds in a group of sequences as a whole, e.g.: all 3'UTRs of a species
	
Usage: ./background_nspk.plx <miRNAs> <background sequence file>

Intructions of use
------------------
1- The subroutines file 'library.pl' must be located either on this folder or any of those contained in @INC.
   Alternatively, its path can be specified in the line number 4,

	 e.g. "require /<my_path>/library.pl"

2- The target sequences file must preferibly be named as follows: <three-letters species code>_Ens<Ensembl version number>_<rest of the name>. 
   The miRNA sequences filename must contain the miRBase version as follows: "miRB<version>".
	 
	 e.g.: mmu_Ens74_all_3UTRs.fa (target sequences filename), mmu_miRB20.fa (miRNA sequences filename)
	
   Otherwise, the results file will be named "enrichment_<input file 1>-<input file 2>".

3- The input sequence files (microRNAs as well as target sequences) must be in FASTA format. There must not be "LRG" nor "Sequences unavailable" entries.

4- This script makes use of all available CPU cores. The user can specify this number by adjusting the value of the variable '$cores' (line number 621).
   Default value: 5.
 
5- The output is a two-column, tab-delimited file containing all miRNAs and seedmatches, as well as their NSPK values in the background sequences.

Best regards,
Danny(dmherrera@unav.es)


