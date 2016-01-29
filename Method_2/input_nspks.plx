#!/usr/bin/perl
use strict;
use warnings;
require 'library.pl';
require Tie::File;
require Tie::File::AsHash;
require Math::Combinatorics;
use Math::Cartesian::Product qw(cartesian);
use Array::Utils qw(array_minus);
use Statistics::Basic qw(stddev mean);
use BSD::Resource qw(setrlimit);
use List::Util qw(sum);
use List::MoreUtils qw(uniq minmax);
use Scalar::Util qw(looks_like_number);
use POSIX qw(isdigit);
$| = 1;


my $test_sets = 10000;												# Random, input-sized test datasets used to resample the background
setrlimit("RLIMIT_NOFILE",$test_sets+100,$test_sets+100);			#Increase maximum file open limit 
																	#(the user must be allowed to do this in /etc/security/limits.conf)

#Script to calculate the enrichment of miRNAs and their seeds in a list of mRNAs

unless ($#ARGV == 2){			# Print help
	print <DATA>;
	exit 0;
}


print "STARTING THE ENRICHMENT IN $ARGV[1]\n";


my $list = $ARGV[1];

my($sp,$mir_ver) = (split /_/,(split /\./,$ARGV[0])[0])[0,1];
my $ens_ver = $1 if $ARGV[3] =~ /(Ens\d+)/;


print "Getting the background nspk....";
my $background_nspk = background_values($ARGV[-1]);	
print "DONE!\n";
pop(@ARGV);

clean_gnome($ARGV[-1]);
tie my @background, 'Tie::File', "temp_$ARGV[-1]" or die "$!: $ARGV[-1]\n";
pop(@ARGV);
																		
foreach my $i(0..$#ARGV){												
	transc_one_line($ARGV[$i]);						# Put the sequences one per line
}
srand;

tie my %input, 'Tie::File::AsHash', "temp_$ARGV[-1]", split => '\t'  or die "$!: $ARGV[-1]\n";
pop(@ARGV);
my $tot_length = length(join('',values %input));					# Total length of input sequences

tie my %mirnas, 'Tie::File::AsHash', "temp_$ARGV[-1]", split => '\t' or die "$!: $ARGV[-1]\n";
pop(@ARGV);

foreach my $file(<temp_*>){unlink($file)};

print "Extracting the seedmatches....";
my ($seeds,$seed_type,$compiled_seeds) = get_seedmatches(\%mirnas);

open my $file,'>',"seedms_miRNAs_structure_${sp}_$mir_ver" or die "$!: seedms_miRNAs_structure_${sp}_$mir_ver";
print $file "seedmatch\tmiRNA\tw-m-wm-no\n";
foreach my $type(keys %$seed_type){
	for(my $i=0;$i<$#{$seed_type->{$type}}-1;$i+=3){
		print $file $seed_type->{$type}[$i],"\t",$seed_type->{$type}[$i+1],"\t";
		if($seed_type->{$type}[$i+2]){
			my $str = join('x',map{$_ ? $_ : ()}@{$seed_type->{$type}});
			if($str =~ /:/ && $str =~ /\*/){
				print $file "wm\n"
			}elsif($str !~ /:/){
				print $file "m\n"
			}elsif($str =~ /\*/){
				print $file "w\n"
			}
		}else{
			print $file "no\n"
		}
	}
}
close $file or die "!: seedms_miRNAs_structure_${sp}_$mir_ver";

print "DONE!\n";

print "Calculating miRNAs' and seeds' enrichment in the input....";
enrich(\%mirnas,\%input,$seeds,$seed_type);
print "DONE!\n";

my ($input_fe,$input_nspk) = get_count($seeds,$tot_length,$background_nspk,\%mirnas);	# Calculate miRNAs' and seedmatches' FE

print "Calculating p-values and z-scores....";
get_nspks($seeds,\@background,scalar keys %input,$compiled_seeds,$test_sets);
my ($pvalues,$zscores) = p_and_z($seeds,$input_fe,$background_nspk,$test_sets);
print "DONE!\n";

untie @background or die "$!: 3UTRome\n";
untie %input or die "$!: input\n";
untie %mirnas or die "$!: miRNAs\n";

open $file,">","enrichment_${list}_${ens_ver}_$mir_ver" or die "$!: enrichment_${list}_${ens_ver}_$mir_ver\n";
print $file "miRNA|Seedm\tZ-score\tp-value\n";
foreach my $seed_mir(sort keys %$zscores){
	print $file "$seed_mir\t$zscores->{$seed_mir}\t$pvalues->{$seed_mir}\n";	
}
close $file or die "$!: enrichment_${list}_${ens_ver}_$mir_ver\n";

open $file, ">", "seedtypes_enrichment_${list}" or die "$!: seedtypes_enrichment_${list}";
print $file "Seedtype\tAvg Z-Score\n";
foreach my $type(keys %$seed_type){
	my @real;
	foreach my $elem(@{$seed_type->{$type}}){
		if ($zscores->{$elem}){
			$zscores->{$elem} = 0 unless looks_like_number($zscores->{$elem});
			push(@real,$zscores->{$elem})
		}
	}
	print $file $type,"\t",mean(@real),"\n"
}
close $file or die "$!: seedtypes_enrichment_${list}";

__DATA__

	Calculates the enrichment (Fold-Enrichment, Z-score and p-value) of miRNAs and their seeds in a group of sequences, e.g.: set of differentially expressed genes of an experiment
	
Usage: ./input_nspks.pl <miRNAs> <target sequences> <background miRNA enrichment file>

Intructions of use
------------------
1- The subroutines file 'library.pl' must be located either on this folder or any of those contained in @INC.
   Alternatively, its path can be specified in line number 4,

	 e.g. "require /<my_path>/library.pl"

2- The target sequences file must preferibly be named as follows: <three-letters species code>_Ens<Ensembl version number>_<rest of the name>. 
   The miRNA sequences filename must contain the miRBase version as follows: "miRB<version>".
	 
	 e.g.: mmu_Ens74_all_3UTRs.fa (target sequences filename), mmu_miRB20.fa (miRNA sequences filename)
	
   Otherwise, the results file will be named "enrichment_<input file 1>-<input file 2>".

3- The input sequence files (microRNAs as well as target sequences) must be in FASTA format. There must not be "LRG" nor "Sequences unavailable" entries.

4- The 'background miRNA enrichment' file is the output of background_nspk.plx. It is a tab-delimited, two-column file containing all miRNAs and seedmatches
   with their respective value of enrichment in the background sequences.

5- This script makes use of all available CPU cores. The user can specify this number by adjusting the value of the variable '$cores' in the subroutines file 'library.pl' (lines number 561 and 703).
   Default value: 5.
   
6- The number of input-sized test dataset files used to resample the background and calculate the p-values can be modified
    by changing the value of the variable $test_sets in line number 19.
    
7- The outputs of this script are two files, named:
	a) 'enrichment_<target file>_<Ensembl version>_<miRBase version>', e.g. 'enrichment_rattus-genes-microarray23.fa_Ens72_miRB20'.
		This is a three-column, tab-delimited file containing all miRNAs and seedmatches, as well as their fold enrichment Z-score in the input target sequences and the p-values of these determinations.
	b)	'seedtypes_enrichment_<target file>', e.g. 'seedtypes_enrichment_rattus-genes-microarray23.fa'. 
		This is a two-column, tab-delimited file containing the average fold enrichment Z-score for each seedmatch type, in the input target sequences.

Best regards,
Danny(dmherrera@unav.es)


