#!/usr/bin/perl
#Script para sacar la posición de inicio de la 3'UTR en el cDNA
use strict;
use warnings;
require 'library.pl';
use Tie::File::AsHash;
$| = 1;


unless($#ARGV == 1){			# Print help
	print <DATA>;
	exit 0
}
die "Usage: ./get_3utr_coords.plx <3'UTRs> <cDNAs>" unless $#ARGV == 1;

`rm 3utr_start*`;
`rm 3utr_NotIncDNA_*`;

my($organism,$ensver) = (split /_/,$ARGV[1])[0,1];

foreach my $file(@ARGV){transc_one_line($file)};
tie my %three,'Tie::File::AsHash',"temp_$ARGV[0]",split => "\t" or die "$!: temp_ARGV[0]";

my @threes = keys %three;
my $cores = 5;
my $group = 0;
my $bin = int(scalar @threes/$cores);
my @packs;
for my $i(0..$#threes){
	if ($i == $group * $bin){
		$group++;
		push(@{$packs[$group-1]},"3utr_starts$group")
	}
	push(@{$packs[$group-1]},[$threes[$i],$three{$threes[$i]}]);

}

my $sons = 1;
for my $procs(0..$#packs){
	if (fork() == 0){
		tie my %cdnas,'Tie::File::AsHash',"temp_$ARGV[1]",split => "\t" or die "$!: temp_ARGV[1]";
		my $filename = shift(@{$packs[$procs]});
		my $notname = 'not_'.$filename;
		open my $file,'>',$filename or die "$!: $filename";
		open my $not,'>',$notname or die "$!: $filename";
		foreach my $k(0..$#{$packs[$procs]}){
			my $name = $packs[$procs][$k][0];
			my $seq = $packs[$procs][$k][1];
			if($cdnas{$name}){
				if($cdnas{$name} =~ /$seq/){
					print $file $name,"\t",$-[0],"\n"
				}else{
					print $file $name,"sequence\n"
				}		
			}else{
				print $not $name,"\t","existence\n"
			}
		}
		close $file or die "$!: $filename";
		$notname or die "$!: $filename";
		untie %cdnas;
		exit()
	}else{
		$sons++
	}
}
wait() for 1..$sons;
untie %three;

`cat 3utr_starts* >> 3utr_start_${organism}_$ensver`;
`rm 3utr_starts*`;

`cat not_* >> 3utr_NotIncDNA_${organism}_$ensver`;
`rm not_*`;


__DATA__

	Launches RNAcofold with the seedmatches and constraints annotated in the 'position' file
	
Usage: ./get_coordinates.plx <3'UTRs> <cDNAs>


Intructions of use
------------------

1- The subroutines file 'library.pl' must be located either on this folder or any of those contained in @INC.
   Alternatively, its path can be specified in line number 4,

	 e.g. "require /<my_path>/library.pl"

2- The sequences files must preferibly be named as follows: <three-letters species code>_Ens<Ensembl version number>_<rest of the name>. 
    
	 e.g.: mmu_Ens74_all_3UTRs.fa (target sequences filename), mmu_miRB20.fa (miRNA sequences filename)

   Otherwise, the results file will be named "enrichment_<input file 1>-<input file 2>"

3- The input sequence files must be in FASTA format. There must not be "LRG" nor "Sequences unavailable" entries.

4- This script makes use of all available CPU cores. The user can specify this number in the variable '$cores' (line number 25).
   Default value: 5.

5- The outputs of this script    

Best regards,
Danny(dmherrera@unav.es)

