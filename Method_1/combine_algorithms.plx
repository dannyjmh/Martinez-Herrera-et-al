#!/usr/bin/perl
use strict;
use warnings FATAL => 'all';
require 'library.pl';
use IO::Handle;
use POSIX;
use List::MoreUtils qw(uniq);
$| = 1;

# Parses the outputs of miRNA prediction algorithms and combine them

print "Deleting the following old combination files (rm *S*) in 10 seconds:\n";
system "ls *S*";
sleep(10);
`rm *S*`;

unless(@ARGV){				# Print help
	print <DATA>;
	exit 0
}

unless(-d $ARGV[0]){
	print "The selected path does not exist\n";
	exit 0
}

my $path= pop(@ARGV);
$path.= "/" unless $path=~ /\/$/;
$path= "/$path" unless $path=~ /^\//;

my @ARGV= reverse glob "${path}*pita_results.tab ${path}*findtar* ${path}*miranda* ${path}*rnahybrid* ${path}*targetscan*";
my($sp,$ens_ver,$mir_ver);

foreach my $file(0..$#ARGV){
	if ($ARGV[$file] =~ /pita/){
		parse_pita($ARGV[$file])
	}elsif ($ARGV[$file] =~ /findtar/){
		my $program = $&;
		parse_findtar($program,$ARGV[$file])
	}elsif($ARGV[$file] =~ /miranda/){
		parse_miranda($&,$ARGV[$file])
	}elsif($ARGV[$file] =~ /rnahybrid/){
		my $a= 2;
		parse_rnahybrid($&,$ARGV[$file],$a)
	}elsif($ARGV[$file] =~ /targetscan/){
		parse_targetscan($ARGV[$file])
	}else{
		die "$file does not contain any algorithm name\n";
	}
	
	unless($ARGV[$file] =~ /pita/){				# Versions of Ensembl and miRBase
		if (!$ens_ver or !$mir_ver){
			($sp,$ens_ver,$mir_ver) = (split /_/,((split /\//,$ARGV[$file])[-1]))[3,4,5];
		}
	}	
	if ($file == $#ARGV){
		unless ($mir_ver && $ens_ver && $sp){
			print "Type the three-letter species ID and versions of Ensembl and miRBase, in that order:\n";
			$sp = <STDIN>;chomp($sp);
			$ens_ver = <STDIN>;chomp($ens_ver);
			$mir_ver = <STDIN>;chomp($mir_ver);	
		}
	}	
}

	# Do the algorithms output combination
	
@ARGV = reverse <*_SITES_3utr>;
summary(\@ARGV,$sp,$ens_ver,$mir_ver) if $#ARGV > 0;


__DATA__

	Parses the outputs of miRNA prediction algorithms and combines them
	
Usage: ./combine_algorithms.plx <path to the files>

Intructions of use
------------------
1- The subroutines file 'library.pl' must be located either on this folder or any of those contained in @INC.
   Alternatively, its path can be specified in the line number 4,

	 e.g. "require /<my_path>/library.pl"

2- This script parses the output of the prediction algorithms: miRanda, PITA, FindTar, TargetScan and RNAhybrid. The output filenames must contain the name of the algorithm,

	e.g. results_findtar_zebrafish
	
3- The outputs of this script are files named as follows:
	<algorithm>_SITES_3utr:	canonical sites predicted by the algorithm
	<algorithm>_NON-CANONICAL_SITES: non-canonical sites predicted by the algorithm
	<algorithm>_NON-CANONICAL_STRUC: non-canonical sites miRNA:seedmatch structure

Best regards,
Danny(dmherrera@unav.es)


