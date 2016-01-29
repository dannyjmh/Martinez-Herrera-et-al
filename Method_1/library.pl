
sub background_values{														#Get the genome nspks
	my(%genome_nspk);
	open my $file,"<",$_[0] or die "$!: $_[0]";
	while(<$file>){
		next if $. == 1;
		chomp;
		my @line = split /\t/;
		$line[0] = $line[0] =~ /_/ ? (split /_/,$line[0])[0] : $line[0];
		$genome_nspk{$line[0]} = $line[1]			
	}
	close $file or die "$!: $_[0]";
	
	return \%genome_nspk;
}

sub background_one_line{
																		#Put genome sequences one per line
																		#(No need for sequence IDs. Genome sequences are only
																		#needed to create random test_sets)
	my $data;	
	open my $file,"<", $_[0] or die "$!: $_[0]";
	open my $dest, ">", "temp_$_[0]" or die "$!: temp_$_[0]";
	while(<$file>){
		die "Sequences must be in FASTA format\n" if ($.==1 && !/^>/);
		die "LRG in line $. of $_[0]. Please remove them before doing the analysis\n" if /LRG/;
		die "Invalid entry in line $.. Please remove them before doing the analysis\n" if /Sequence\s+unavailable/;
		chomp;
		if (/^>/){
			print $dest "\n" unless ($. == 1);
			next
		}else{
			print $dest "$_";
		}
	}
	close $file or die "$!: $_[0]";
	close $dest or die "$!: temp_$_[0]";
}

sub transc_one_line{
																#Put input sequences one per line
	my $data;	
	open my $file,"<", $_[0] or die "$!: $_[0]";
	open my $dest, ">", "temp_$_[0]" or die "$!: temp_$_[0]";
	while(<$file>){
		die "Sequences must be in FASTA format\n" if ($.==1 && !/^>/);
		die "LRG in line $. of $_[0]. Please remove them before doing the analysis\n" if /LRG/;
		die "Invalid entry in line $.. Please remove them before doing the analysis\n" if /Sequence\s+unavailable/;
		chomp;
		if (/^>/){
			print $dest "\n" unless ($. == 1);
			if(/\|/){
				my @data = (split /\|/,(split />/)[1])[0,1];
				$data = join('_', @data);
			}else{
				$data = (split />/)[1];
			}
			print $dest "$data\t";
		}elsif(/^\w+/){
			print $dest "$_";
		}else{
			next
		}
	}
	close $file or die "$!: $_[0]";
	close $dest or die "$!: temp_$_[0]";
}

sub get_seedmatches{							# Create the seedmatches of every miRNA
	my (%seeds,%wobbles,%feats,%misma,%types);	
	my @sizes = (6..10);						#Seedmatch sizes
	my @limits = (								#Seed start,end
			"1,6",
			"2,7",
			"3,8",
			"1,7",
			"2,8",
			"1,8",
			"2,9",
			"1,9",
			"2,10",
			"1,10",
			"2,11",
			"3,12",
	);

											#Assign seed start and end positions to each seedmatch size
	@feats{@sizes} = ([@limits[0,1,2]],[@limits[3,4]],[@limits[5,6]],[@limits[7,8]],[@limits[9,10,11]]);
	
	@wobbles{@sizes} = (1,1,2,2,3);			#wobbles allowed for each seed size (Didiano and Hobert,2006)
	@misma{@sizes} = (0,0,1,1,2);			#mismatches allowed for each seed size (~ Grosswendt 2014, Helwak 2013,Friedman 2009, Vella 2004 y Lal 2009)

	%types = (
		$limits[0] => "6nt(1)",
		$limits[1] => "6nt(2)",
		$limits[2] => "6nt(3)",
		$limits[3] => "7nt(1)",
		$limits[4] => "7nt(2)",
		$limits[5] => "8nt(1)",
		$limits[6] => "8nt(2)",
		$limits[7] => "9nt(1)",
		$limits[8] => "9nt(2)",
		$limits[9] => "10nt(1)",
		$limits[10] => "10nt(2)",
		$limits[11] => "10nt(3)",
	);

	
	my (%seed_type,$seedm,$mirnas);
	$mirnas = pop(@_);

	foreach my $mir(keys %$mirnas){
		my $seq = $mirnas->{$mir};
		my @sps = split(//,$seq);
		local ($a,$b);
		foreach my $size(sort {$a<=>$b} keys %feats){
			UPPER: foreach my $i(0..$#{$feats{$size}}){
				my @boundaries = split(/,/,$feats{$size}[$i]);
														
										#Get the miRNA nt's between the seed start and end
				$seedm = substr $seq,$boundaries[0]-1,$boundaries[1]-$boundaries[0]+1;
				$seedm =~ tr/AUCGT/TAGCA/;
				$seedm = my $seedm_copy = scalar reverse $seedm;
				
				my @variants;			#Cartesian product to create all variants of wobbles and mismatches
				cartesian {push(@variants,[@_])} [0..$wobbles{$size}], [0..$misma{$size}];
				
											
				my $done = 0;				#Check if 7mer-1A has been created
				my $mer7;					
				CORE:foreach my $combo(0..$#variants){	#Send to &get_wobbles and &get_mismatches according to the @variants array values

					next if $seedm_copy !~ /A|C/ && $misma{$size} == 0 && $combo < $#variants;
					
					my ($wobbs,$positions) = get_wobbles($seedm_copy,$variants[$combo]->[0],$mer7);

					my($mis,$struc) = get_mismatches($wobbs,$variants[$combo]->[1],$positions);
					
					foreach my $k(0..$#$mis){
						my $seed = $mis->[$k];
						push(@{$seeds{$seed}},$mir);
						
						if ($done == 1){	#If the 7mer-1A has been created
							push(@{$seed_type{"7mer-1A"}},$seed,$mir,$struc->{$seed});
							push(@{$seeds{$seed}},0);
							next UPPER if ($k == $#$mis && $combo == $#variants)
						}else{
							push(@{$seed_type{$types{$feats{$size}[$i]}}},$seed,$mir,$struc->{$seed});
							push(@{$seeds{$seed}},$boundaries[0]-1);								
						}
					}
				}

				if ($feats{$size}[$i] eq '2,7'){#Create the 7mer-1A(6nt(2)+A)
					$seedm .= 'A';
					$seedm_copy = $seedm;
					undef(@variants);
					cartesian {push(@variants,[@_])} [0..1], [0];
					$done++;
					$mer7++;
					goto CORE
				}
			}	
		}						
		
	}
	
	my (@temp,@compiled_seeds);
																		#Compile the seedmatches for faster lookup
	my @seeds = keys %seeds;
	@temp = map {
	    my $code = "sub { return scalar (() = shift =~ /\Q$_\E/g) }";
	    my $sub = eval $code;
	    die "Unable to cache sub $code: $@" if $@;
	    $sub
	}@seeds;
	
						#@compiled seeds = [compiled 1][normal1]......[compiled n][normal n]
	foreach my $i(0..$#temp){$compiled_seeds[$i] = [$temp[$i],$seeds[$i]]}
	return \%seeds,\%seed_type,\@compiled_seeds
}

sub get_seedmatches_genome{							# Create the seedmatches of every miRNA
	my (%seeds,%wobbles,%feats,%misma,%types);	
	my @sizes = (6..10);						#Seedmatch sizes
	my @limits = (								#Seed start and end
			"1,6",
			"2,7",
			"3,8",
			"1,7",
			"2,8",
			"1,8",
			"2,9",
			"1,9",
			"2,10",
			"1,10",
			"2,11",
			"3,12",
	);

											#Assign seed start and end positions to each seedmatch size
	@feats{@sizes} = ([@limits[0,1,2]],[@limits[3,4]],[@limits[5,6]],[@limits[7,8]],[@limits[9,10,11]]);
	
	@wobbles{@sizes} = (1,1,2,2,3);			#wobbles allowed for each seed size (Didiano and Hobert,2006)
	@misma{@sizes} = (0,0,1,1,1);			#mismatches allowed for each seed size (~ Grosswendt 2014, Helwak 2013,Friedman 2009, Vella 2004 y Lal 2009)

	%types = (
		$limits[0] => "6nt(1)",
		$limits[1] => "6nt(2)",
		$limits[2] => "6nt(3)",
		$limits[3] => "7nt(1)",
		$limits[4] => "7nt(2)",
		$limits[5] => "8nt(1)",
		$limits[6] => "8nt(2)",
		$limits[7] => "9nt(1)",
		$limits[8] => "9nt(2)",
		$limits[9] => "10nt(1)",
		$limits[10] => "10nt(2)",
		$limits[11] => "10nt(3)",
	);

	
	my (%seed_type,$seedm,$mirnas);
	$mirnas = pop(@_);
	
	my %structure;						#Store the seed's pairing scheme for the constrained thermodynamic calculations

	foreach my $mir(keys %$mirnas){
		my $seq = $mirnas->{$mir};
		my @sps = split(//,$seq);
		local ($a,$b);
		foreach my $size(sort {$a<=>$b} keys %feats){
			UPPER: foreach my $i(0..$#{$feats{$size}}){
				my @boundaries = split(/,/,$feats{$size}[$i]);
														
								#Get the miRNA nt's between the seed start and end
				$seedm = substr $seq,$boundaries[0]-1,$boundaries[1]-$boundaries[0]+1;
				$seedm =~ tr/AUCGT/TAGCA/;
				$seedm = my $seedm_copy = scalar reverse $seedm;
				
				my @variants;			#Cartesian product to create all variants of wobbles and mismatches
				cartesian {push(@variants,[@_])} [0..$wobbles{$size}], [0..$misma{$size}];
				
											
				my $done = 0;				#Check if 7mer-1A has been created
				my $mer7 = 0;
				CORE:foreach my $combo(0..$#variants){	#Send to &get_wobbles and &get_mismatches according to the @variants array values
					my ($wobbs,$positions) = get_wobbles($seedm_copy,$variants[$combo]->[0],$mer7);

					my($mis,$struc) = get_mismatches($wobbs,$variants[$combo]->[1],$positions);
					
					
					foreach my $seed(@$mis){
						push(@{$seeds{$seed}},$mir);
						push(@{$seed_type{$types{$feats{$size}[$i]}}},$seed);
					}

					if ($combo == $#variants && $done == 1){
						foreach my $seed(@$mis){
							push(@{$seed_type{"7mer-1A"}},$seed);
						}
						next UPPER;								
					}
					
				}
				if ($feats{$size}[$i] eq '2,7'){				#Create the 7mer-1A(6nt(2)+A)
					$seedm .= 'A';
					$seedm_copy = $seedm;
					undef(@variants);
					cartesian {push(@variants,[@_])} [0..1], [0];
					$done++;
					$mer7++;
					goto CORE
				}
			}	
		}
	}
	my (@temp,@compiled_seeds);
																	#Compile seedmatches for faster lookup
	my @seeds = keys %seeds;
	@temp = map {
	    my $code = "sub { return scalar (() = shift =~ /\Q$_\E/g) }";
	    my $sub = eval $code;
	    die "Unable to cache sub $code: $@" if $@;
	    $sub
	}@seeds;
	
						#@compiled seeds = [compiled 1][normal1]......[compiled n][normal n]
	foreach my $i(0..$#temp){$compiled_seeds[$i] = [$temp[$i],$seeds[$i]]}

	foreach my $keys(keys %seeds){
		@{$seeds{$keys}} = uniq(@{$seeds{$keys}});	
	}

	return \%seeds,\%seed_type,\@compiled_seeds
}


sub get_seedmatches_one_seedtype{							# Create seedmatches of one type
	my (%seeds,%wobbles,%feats,%misma,%types,$seedt,$wm);
	
	my (%seed_type,$seedm,$mirnas);
	($mirnas,$seedt,$wm,@_) = @_;	
		
	my @sizes = (6..10);						#Seedmatch sizes
	my @limits = (								#Seed start and end
			"1,6",
			"2,7",
			"3,8",
			"1,7",
			"2,8",
			"1,8",
			"2,9",
			"1,9",
			"2,10",
			"1,10",
			"2,11",
			"3,12",
	);
											#Assign seed start and end positions to each seedmatch size
	@feats{@sizes} = ([@limits[0,1,2]],[@limits[3,4]],[@limits[5,6]],[@limits[7,8]],[@limits[9,10,11]]);
		
	if ($wm eq 'wm'){
		@wobbles{@sizes} = (1,1,2,2,3);			#wobbles allowed for each seed size (Didiano and Hobert,2006)
		@misma{@sizes} = (0,0,1,1,1);			#mismatches allowed for each seed size (~ Grosswendt 2014, Helwak 2013,Friedman 2009, Vella 2004 y Lal 2009)
	}elsif($wm eq 'nowm'){
		@wobbles{@sizes} = (0,0,0,0,0);
		@misma{@sizes} = (0,0,0,0,0);
	}elsif($wm eq 'w'){
		@wobbles{@sizes} = (1,1,2,2,3);
		@misma{@sizes} = (0,0,0,0,0);
	}elsif($wm eq 'm'){
		@wobbles{@sizes} = (0,0,0,0,0);
		@misma{@sizes} = (0,0,1,1,1);
	}

	%types = (
		$limits[0] => "6nt_1",
		$limits[1] => "6nt_2",
		$limits[2] => "6nt_3",
		$limits[3] => "7nt_1",
		$limits[4] => "7nt_2",
		$limits[5] => "8nt_1",
		$limits[6] => "8nt_2",
		$limits[7] => "9nt_1",
		$limits[8] => "9nt_2",
		$limits[9] => "10nt_1",
		$limits[10] => "10nt_2",
		$limits[11] => "10nt_3",
	);

	my $check_seedtype = join(',',sort values %types). ',7mer-1A';
	die "The selected seedtype is invalid. It must be one of: ",$check_seedtype,"\n" unless $check_seedtype =~ /$seedt/;

	foreach my $mir(keys %$mirnas){
		my $seq = $mirnas->{$mir};
		my @sps = split(//,$seq);
		local ($a,$b);
		foreach my $size(sort {$a<=>$b} keys %feats){
			UPPER: foreach my $i(0..$#{$feats{$size}}){
				next unless $seedt eq $types{$feats{$size}[$i]};
				my @boundaries = split(/,/,$feats{$size}[$i]);
														
								#Get the miRNA nt's between the seed start and end
				$seedm = substr $seq,$boundaries[0]-1,$boundaries[1]-$boundaries[0]+1;
				$seedm =~ tr/AUCGT/TAGCA/;
				$seedm = my $seedm_copy = scalar reverse $seedm;
				
				my @variants;			#Cartesian product to create all variants of wobbles and mismatches
				cartesian {push(@variants,[@_])} [0..$wobbles{$size}], [0..$misma{$size}];
											
				my $done = 0;				#Check if 7mer-1A has been already created
				my $mer7;					
				CORE:foreach my $combo(0..$#variants){	# Create the seedmatch variants containing wobbles and mismatches
					my ($wobbs,$positions) = get_wobbles($seedm_copy,$variants[$combo]->[0],$mer7);

					my($mis,$struc) = get_mismatches($wobbs,$variants[$combo]->[1],$positions);
					
					foreach my $k(0..$#$mis){
						my $seed = $mis->[$k];
						push(@{$seeds{$seed}},$mir);
						if ($done == 1){								
							push(@{$seed_type{"7mer-1A"}},$seed,$mir,$struc->{$seed});
							push(@{$seeds{$seed}},0);
							next UPPER if ($k == $#$mis && $combo == $#variants)
						}else{
							push(@{$seed_type{$types{$feats{$size}[$i]}}},$seed,$mir,$struc->{$seed});
							push(@{$seeds{$seed}},$boundaries[0]-1);
						}
					}						
				}

				if ($feats{$size}[$i] eq '2,7' && $seedt eq '7mer-1A'){				#Create the 7mer-1A(6nt(2)+A)
					$seedm .= 'A';
					$seedm_copy = $seedm;
					undef(@variants);
					cartesian {push(@variants,[@_])} [0..1], [0];
					$done++;
					$mer7++;
					goto CORE
				}
			}	
		}			
	}
	
	return \%seeds,\%seed_type
}


sub get_wobbles{														#Insert the wobbles
	my ($seedm,$trigger,$mer7);											# trigger: number of As to replace every time
	($seedm,$trigger,$mer7,@_) = @_;
	
	my @exchange = split(//,$seedm);
	my (@wobbles,%positions,@numbers);

	my $uracils = () = $seedm =~ /A|C/g;
	
	$uracils-- if $mer7;

	if($trigger && $uracils){
		
		my %repl = (
			'A' => 'G',
			'C' => 'T'
			);
			
		for my $change(qw(A C)){
			while($seedm =~ /$change/g){
				push(@numbers,$-[0])
			}
		}

		if ($mer7){															#Do not replace the 7mer-1A last 'A'
	 		@numbers = sort @numbers;
			pop(@numbers);
		}

												#Replace the triggers
		my $combination = Math::Combinatorics->new(count => $trigger,
													data => [@numbers],
												);
		while(my @combo = $combination->next_combination){
			my @chain = @exchange;
			for my $pos(@combo){
				$chain[$pos] = $repl{$chain[$pos]}
			}
			push(@wobbles,join('',@chain));

									#Find the positions of replacements. Store them in $positions{seed}= [pos 1]....[pos n]									
			my $diff = $seedm ^ $wobbles[-1];
			while ($diff =~ /[^\0]/g){
				push(@{$positions{$wobbles[-1]}},$-[0]);
			}
		}		
	}else{
		push(@wobbles, $seedm);
		$positions{$wobbles[-1]}[0] = 20;	#If no uracils or no substitution programmed for this position
											#assign an unreachable number. '0' not allowed because of the use
											#of the operator array_minus
	}
	
	@wobbles = uniq(@wobbles);

	return \@wobbles,\%positions
}

sub get_mismatches{									#Insert the mismatches
	my ($exchange,$mism_number,$positions,$struc);
	($exchange,$mism_number,$positions,@_) = @_;

	my %repl = (
		"A" => ["C","G"],
		"U" => ["C","T","G"],
		"T" => ["C","G"],
		"C" => ["A","T"],
		"G" => ["A","T"]
	);
	
	my @mismatches;

	foreach my $seed(@$exchange){
		unless($positions->{$seed}[0] == 20){
			foreach my $pos(@{$positions->{$seed}}){							#Annotate the mismatch in the seedmatch structure
				$struc->{$seed}[$pos] = ':';
			}
		}
		if($mism_number){	
			my @seeds = my @seeds_orig = split(//,$seed);
			my @numbers = (0..$#seeds-1);
			@numbers = array_minus(@numbers, @{$positions->{$seed}});
			
			my $combination = Math::Combinatorics->new(count => $mism_number,
														data => [@numbers],
													);												
			while(my @combo = $combination->next_combination){
				if($mism_number == 1){								#If only 1 mismatch to insert
					foreach my $alt(@{$repl{$seeds[$combo[0]]}}){	#For each of the mismatch variants
						$seeds[$combo[0]] = $alt;
						push(@mismatches,join('',@seeds));
						@{$struc->{$mismatches[-1]}} = @{$struc->{$seed}} if $struc->{$seed};	#Annotate the wobbles of the seedmatch in its structure variable
						$struc->{$mismatches[-1]}[$combo[0]] = '*';								#Annotate the mismatch
						@seeds = @seeds_orig;													#Back to the original seedmatch before next insertions
					}
				}else{												#If inserting 2 mismatches
					foreach my $alt(@{$repl{$seeds[$combo[0]]}}){
						$seeds[$combo[0]] = $alt;
						my @seeds_sub1 = @seeds;						#Keep a copy of the single-mismatched seedmatch
						foreach my $alt1(@{$repl{$seeds[$combo[1]]}}){	#For each variant of the second mismatch
							$seeds[$combo[1]] = $alt1;
							push(@mismatches,join('',@seeds));
							@{$struc->{$mismatches[-1]}} = @{$struc->{$seed}} if $struc->{$seed};   #Annotate the wobbles of the seedmatch in its structure variable
							foreach my $pos(@combo){
								$struc->{$mismatches[-1]}[$pos] = '*';
							}
							@seeds = @seeds_sub1;
						}
						@seeds = @seeds_orig;													#Back to the original seedmatch before next insertions
					}
				}				
			}
		}else{
			push(@mismatches,$seed);
			$struc->{$seed} = 1 if $positions->{$seed}[0] == 20;
		}
	}
	@mismatches = uniq(@mismatches);
	
	return \@mismatches,$struc
}

sub get_struc{															#annotate structure of the seedmatch: struc{seed}{mirna} = [position_feature]
																		#eg.: struc{ATGCTA}{mmu-miR-201} = [2_:][4_*]
	my $seed_type = pop(@_);
	my %struc;
	
	foreach my $type(keys %$seed_type){
		for(my $i=0;$i<=$#{$seed_type->{$type}}-2;$i+=3){
			my @conv;
			for my $j(0..$#{$seed_type->{$type}[$i+2]}){
				if($seed_type->{$type}[$i+2][$j]){
					push(@conv,"${j}_$seed_type->{$type}[$i+2][$j]")
				}
			}
			$struc{$seed_type->{$type}[$i]}{$seed_type->{$type}[$i+1]} = join('-',@conv);
		}
	}
	
	return \%struc
}

sub enrich{ 															#Calculate miRNAs' and seedmatches' nspk in the input
	my($mirnas,$input,$seeds,$seed_type);
	($mirnas,$input,$seeds,$seed_type,@_) = @_;
	
	`rm position` if -e 'position';
	
	my $struc = get_struc($seed_type);									
	
	my @seeds = keys %$seeds;
	my $cores = 5;
	my $group = 0;
	my $bin = int((scalar keys %$seeds)/$cores);
	my @packs;
	for my $i(0..$#seeds){
		if ($i == $group * $bin){
			$group++;
			push(@{$packs[$group-1]},"positions$group")
		}
		push(@{$packs[$group-1]},$seeds[$i]);

	}
		
	my @keys = keys %$input;
	
	my $sons = 1;
	for my $procs(0..$cores-1){		
		if (fork() == 0){
			my $filename = shift(@{$packs[$procs]});
			open my $file,'>',$filename or die "$!: $filename";
			foreach my $i(0..$#keys){
				my $seq = $keys[$i];
				my $str = $input->{$seq};
				foreach my $seed(@{$packs[$procs]}){
					while($str =~ /$seed/g){
						print $file ">\t",$seed,"\t",$seq,"\t",$-[0],"\n";
						foreach my $i(0..$#{$seeds->{$seed}}){
							print $file $seeds->{$seed}[$i];							
							if(exists $seeds->{$seed}[$i+1]){
								if(isdigit($seeds->{$seed}[$i+1])){
									print $file ','
								}else{
									print $file ',',$struc->{$seed}{$seeds->{$seed}[$i-1]}
								}
							}else{
								print $file ',',$struc->{$seed}{$seeds->{$seed}[$i-1]}
							}
						}print $file "\n";
					}
				}					
			}
			close $file or die "$!: $filename";
			exit();
		}else{
			$sons++		
		}
	}
				
	wait() for 1..$sons;
	
	`cat positions* >> position`;
	`rm positions*`;

}

sub enrich_background{ 														#Calculate miRNAs' and seedmatches' nspk in the genome
	my($mirnas,$genome,$seeds,$compiled_seeds,$tot_length);
	($mirnas,$genome,$seeds,$compiled_seeds,$tot_length,@_) = @_;
	
	my $number = scalar @$compiled_seeds;
	my $cores = 5;
	my $group = 0;
	my $bin = int($number/$cores);
	my @packs;
	my $seqs = join('x',@$genome);
	study($seqs);
	
	for my $i(0..$number-1){
		if ($i == $group * $bin){
			$group++;
			push(@{$packs[$group-1]},"nspks$group")
		}
		push(@{$packs[$group-1]},[@{$compiled_seeds->[$i]}]);
	}

	my $sons = 1;
	for my $procs(0..$#packs){
		if (fork() == 0){
			my $filename = shift(@{$packs[$procs]});
			open my $file,'>',$filename or die "$!: $filename";				
			foreach my $k(0..$#{$packs[$procs]}){
				my $comp = $packs[$procs][$k][0];
				my $seed = $packs[$procs][$k][1];
				my $nspk = $comp->($seqs);					#Count the total number of sites per seedmatch and 3'UTR
				print $file ">\t",$seed,"\t",$nspk,"\n";
				@{$seeds->{$seed}} = uniq(@{$seeds->{$seed}});
				print $file join(',',@{$seeds->{$seed}}),"\n";
			}
			close $file or die "$!: $filename";
			exit()
		}else{
			$sons++
		}
	}
	wait() for 1..$sons;

	`cat nspks* >> nspk`;
	`rm nspks*`;
	
	
	my %count;
	my $temp;
	open my $file,'<','nspk' or die "$!: nspk";
	while(<$file>){
		chomp;
		if(/^>/){			
			my @data = split /\t/;
			$temp = $data[2];
			$count{$data[1]} += $data[2];
		}else{
			my @mirs = split(/,/);
			foreach my $mirna(@mirs){
				$count{$mirna} += $temp;
			}undef($temp);
		}
	}
	close $file  or die "$!: nspk";

	foreach my $elem(keys %count){
		$count{$elem} *= 1000/$tot_length;
	}
	
	return \%count

}

sub get_nspks{															#Create the 10000 random, input-sized test sets from the genome
	my($seeds,$genome,$lines,$compiled_seeds,$test_sets);
	($seeds,$genome,$lines,$compiled_seeds,$test_sets,@_) = @_;

	my @test_sets;
	for my $j(0..$test_sets-1){
		open my $file,">","temp_test_set$j" or die "$!: temp_test_set$j\n";
		for my $k(1..$lines){
			my $number = scalar(@$genome);
			print $file "$genome->[int(rand($number))]\n";
		}
		close $file or die "$!: temp_test_set$j\n";
	}
	
	#------------------------------
	my $number = scalar @$compiled_seeds;
	my $cores = 5;
	my $group = 0;
	my $bin = int($number/$cores);
	my @packs;
	for my $i(0..$number-1){
		if ($i == $group * $bin){
			$group++;
			push(@{$packs[$group-1]},"nspks$group")
		}
		push(@{$packs[$group-1]},[@{$compiled_seeds->[$i]}]);

	}

	my $sons = 1;
	for my $procs(0..$#packs){
		if (fork() == 0){			
			for my $j(0..$test_sets-1){		
				tie @{$test_sets[$j]}, 'Tie::File', "temp_test_set$j" or die "Can't tie: temp_test_set$j";
				unlink("temp_test_set$j");
			}			
			my $filename = shift(@{$packs[$procs]});
			open my $file,'>',$filename or die "$!: $filename";				
			foreach my $k(0..$#{$packs[$procs]}){
				my $comp = $packs[$procs][$k][0];
				my $seed = $packs[$procs][$k][1];
				my @nspks;
				foreach my $j(0..$#test_sets){
					my $seqs = join('X',@{$test_sets[$j]});
					study($seqs);
					my $tot_length = length($seqs) - $lines + 1;
					my $nspk = $comp->($seqs) * 1000/$tot_length;
					push(@nspks,$nspk)									#nspks of every seedmatch in the test sets
				}
				print $file $seed,"\t",join("\t",@nspks),"\n";
			}
			close $file or die "$!: $filename";
			exit()
		}else{
			$sons++
		}
	}
	wait() for 1..$sons;

	`cat nspks* >> nspk`;
	`rm nspks*`;
	
}

sub p_and_z{															#Calculate p-values and z-scores
	my ($seeds,$input_fe,$genome_nspks,$test_sets);
	($seeds,$input_fe,$genome_nspks,$test_sets,@_) = @_;
	
	
	tie my @nspks, 'Tie::File','nspk' or die "$!: nspk";
	unlink("nspk");
	
	my (%zscores,%pvalues,%mir_nspks);

	local($a,$b);

	#-------------------seeds
	
	for my $values(@nspks){
		
		my @data = split(/\t/,$values);
		my $seed = shift(@data);

		my $mir;
		for my $i(0..$#data){
			for (my $j=0;$j<=$#{$seeds->{$seed}};$j+=2){
				$mir = $seeds->{$seed}[$j];
				$mir_nspks{$mir}[$i] += $data[$i]						#Add up miRNAs nspks
			}
					
			if ($data[$i] != 0 && $genome_nspks->{$seed}){
				$data[$i] = log($data[$i]/$genome_nspks->{$seed})		#nspk -> FE
			}else{
				$data[$i] = 0
			}
		}		
		
		@data = sort {$b <=> $a} @data;
		my $stdev = stddev(@data);
		
		unless($stdev == 0){
			my $mean = mean(@data);

			$zscores{$seed} = ($input_fe->{$seed} - $mean)/$stdev;
		
			$pvalues{$seed} = 0;
			foreach my $elem(@data){
				if($elem >= $input_fe->{$seed}){
					$pvalues{$seed}++
				}else{					
					last
				}
			}
			
			if ($pvalues{$seed} == 0){
				$pvalues{$seed} = '<'. 1/$test_sets;
			}elsif($pvalues{$seed} == $test_sets){
				$pvalues{$seed} = "~1"
			}else{
				$pvalues{$seed} = sprintf("%.5f",$pvalues{$seed}/$test_sets);				
			}		
		}else{
			$zscores{$seed} = $pvalues{$seed} = 'undef';			 
		}
	}
	
	#-------------------seeds
	
	#-------------------miRNAs
	
	foreach my $mir(keys %mir_nspks){
		foreach my $value(@{$mir_nspks{$mir}}){
			if ($value != 0 && $genome_nspks->{$mir}){
				$value = log($value/$genome_nspks->{$mir})			#nspk -> FE
			}else{
				$value = 0
			}
		}		
		
		my $stdev = stddev(@{$mir_nspks{$mir}});
		
		unless($stdev == 0){
			my $mean = mean(@{$mir_nspks{$mir}});		
			$zscores{$mir} = sprintf("%.5f",($input_fe->{$mir} - $mean)/$stdev);	
		
			@{$mir_nspks{$mir}} = sort {$b <=> $a} @{$mir_nspks{$mir}};
			
			$pvalues{$mir} = 0;
			foreach my $value(@{$mir_nspks{$mir}}){
				$value = 0 unless $value;				
				if($value >= $input_fe->{$mir}){
					$pvalues{$mir}++
				}else{
					last
				}
			}	
			if ($pvalues{$mir} == 0){
				$pvalues{$mir} = '<'. 1/$test_sets;
			}elsif($pvalues{$mir} == $test_sets){
				$pvalues{$mir} = "~1"
			}else{
				$pvalues{$mir} = sprintf("%.5f",$pvalues{$mir}/$test_sets);				
			}		
		}else{			
			$zscores{$mir} = $pvalues{$mir} = 'undef'
		}
	}
	#-----------------miRNAs

	untie @nspks;
	return \%pvalues,\%zscores
}

sub get_coords{													#Get start position of the 3'UTRs in their cDNAs
	my %coords;

	open my $file, "<", $_[0] or die "$!: $_[0]\n";
	while(<$file>){
		next if $.== 1;
		chomp;
		s/\|/_/g;		
		my @data = />/ ? split(/\t/, (split />/)[1]):split(/\t/);
		$coords{$data[0]} = $data[1];
	}
	close $file or die "$!: $_[0]\n"; 
	
	return \%coords
}

sub make_array_struc{											#Create structure array
	my @struc;
	my @array = split(/-/,$_[0]);
	foreach my $elem(@array){
		my @elems = split(/_/,$elem);
		$struc[$elems[0]] = $elems[1]
	}
	
	return \@struc
}

sub dot_bracket{												#Write the dot-bracket notation of the miRNA-target pair
	my($seed,$mirna,$start_seed,$pos,$cdna,$str,$struc);
	($seed,$mirna,$start_seed,$pos,$cdna,$str,@_) = @_;
	
	if(ref($str) eq 'ARRAY'){
		$struc = $str 
	}else{
		$struc = make_array_struc($str)
	}	
	
	my $separator = '&';
	my $pair;
	my $end_seed = $start_seed+length($seed)-1;
	my @base = split(//,'.' x length($mirna).$separator.'.' x length($cdna));

									#Distance between last nt of the seed and the first nt of the seedmatch
	my $in_between = length($mirna) - 1 - $end_seed + length($separator) + $pos;	
	for(my $i=$start_seed;$i<=$end_seed;$i++){
		$pair = $i+($end_seed-$i)*2+$in_between+1;
		if ($struc && $struc->[$i] && $struc->[$i] eq '*'){
			$base[$i] = $base[$pair] = 'x';
		}else{
			$base[$i] = '(';
								#closing parenthesis at a distance twice that between $i and the end of the seed + $in_between
			$base[$pair] = ')';
		}
	}	
	my $constr = join('',@base)."\n";
	return $constr;
}

sub dot_bracket_70{		#Dot-bracket notation of -70nts........miRNA..........70nts
	my($seed,$mirna,$start_seed,$pos,$cdna,$struc);
	($seed,$mirna,$start_seed,$pos,$cdna,$struc,@_) = @_;
		
	my $separator = '&';
	my $pair;
	my $end_seed = $start_seed+length($seed)-1;
	
#Check if there are at least 70nt upstream the seedmatch. If not, extend downstream or add polyA							
	my $get_left = 70 - $pos;
	my $add_to_right = 0;
	
	if($get_left > 0){
		$add_to_right = $get_left + $start_seed + 1;
		$get_left = 70 - $get_left - $start_seed - 1;
	}else{
		$get_left = 70;
	}
										
	my $get_right += $add_to_right;		#nt's upstream
	my $polyA;							#polyA
	
	my $plus70 = $pos + length($mirna) - $start_seed + 70 + $get_right - length($cdna);
	if($plus70 > 0){
		$polyA = 'A' x $plus70;
		$get_right = 70 + $get_right - $plus70;
	}else{
		$get_right += 70;
	}

	my $seedm_region = substr($cdna,$pos-$get_left-$start_seed,$get_left);			#nt's upstream the seedmatch
	$seedm_region .= substr($cdna,$pos-$start_seed,length($mirna));					#miRNA
	$seedm_region .= substr($cdna,$pos-$start_seed+length($mirna),$get_right);		#nt's downstream the miRNA
	$seedm_region .= $polyA if $polyA;												#la polyA
	
	my @base = split(//,'.' x length($mirna).$separator.'.' x length($seedm_region));

									#Distance between last nt of the seed and the first nt of the seedmatch
	my $pos_in_70 = $-[0] if($seedm_region =~ /$seed/);
	my $in_between = length($mirna) - 1 - $end_seed + length($separator) + $pos_in_70;
	for(my $i=$start_seed;$i<=$end_seed;$i++){
		$pair = $i+($end_seed-$i)*2+$in_between+1;
		if ($struc && $struc->[$i] && $struc->[$i] eq '*'){
			$base[$i] = $base[$pair] = 'x';
		}else{
			$base[$i] = '(';
								#closing parenthesis at a distance twice that between $i and the end of the seed + $in_between
			$base[$pair] = ')';
		}
	}	
	my $constr = join('',@base)."\n";
	return $constr,$seedm_region;
}

sub get_count{											#Get the total number of sites per miRNA and seed
	my($seeds,$count,$tot_length,$genome_nspk,$mirnas);
	($seeds,$tot_length,$genome_nspk,$mirnas,@_) = @_;

	my %count;
	my @keys = keys %$seeds;
	@count{@keys} = (0) x scalar @keys;
	open my $file,'<','position' or die "$!: position";
	while(<$file>){
		if(/^>/){
			my $seed = (split /\t/)[1];
			$count{$seed}++
		}
		else{
			chomp;
			my @mirdata = split(/;/);
			@mirdata = uniq(@mirdata);
			foreach my $mirdata(@mirdata){
				my @mirs = split(/,/,$mirdata);
				my $mirna = shift(@mirs);
				$count{$mirna}++;
			}
		}
	}
	close $file  or die "$!: position";

	my %input_fe;			
	foreach my $seed_mir(keys %count){							#Calculate FE of miRNAs and seeds in the input
		$count{$seed_mir} *= 1000/$tot_length;					#nspks yet
		if ($count{$seed_mir} != 0 && $genome_nspk->{$seed_mir} != 0){
			$input_fe{$seed_mir} = log($count{$seed_mir}/$genome_nspk->{$seed_mir});	#FEs already
		}else{
			$input_fe{$seed_mir} = 0;
		}
	}

	foreach my $mir(keys %$mirnas){
		$input_fe{$mir} = 0 unless $input_fe{$mir}
	}

	foreach my $sem(keys %$seeds){
		$input_fe{$sem} = 0 unless $input_fe{$sem}
	}

	return \%input_fe,\%count
}

sub get_count_one_seedtype{						#Get the total number of sites per miRNA and seed of one seedtype only
	my($seeds,$cdnas,$three_utrs,$mirnas,$seed_type,$coords,$sites,$fixed_start);
	($seeds,$cdnas,$three_utrs,$mirnas,$seed_type,$coords,$sites,$fixed_start,@_) = @_;
	
	my %all_data;
	foreach my $typ(keys %$seed_type){						
		for (my $i=0;$i<$#{$seed_type->{$typ}}-2;$i+=3){
			$all_data{"${seed_type->{$typ}[0]}_${seed_type->{$typ}[1]}"} = $seed_type->{$typ}[2]
												#$all_data{"seed_miRNA"} = structure
		}						
	}	
	
	my %count;
	my %calculations;
	my ($seed,$target,$pos,$mir);
	open my $file,'<',$sites or die "$!: $sites";
	while(<$file>){
		chomp;
		if(/^>/){
			($seed,$target,$pos) = (split /\t/)[1..3];

			unless($seeds->{$seed}){
				undef($seed);
				next
			}
			$pos += $coords->{$target} if($coords->{$target});
		}else{
			next unless($seed);
			my @line = split(/;/);
											#check if if the seedmatch in file 'position' is reported for that miRNA
											#in %seeds, with that exact start position
			my $check_mir = join(',',@{$seeds->{$seed}});
			my @mirs;
			foreach my $elem(@line){
				push(@mirs,$elem) if ($check_mir =~ /$elem/)
			}
			next unless @mirs;
											#-------------------------------------------------------------------------
			for my $i(0..$#mirs){
				my @pairs = split(/,/,$mirs[$i]);
				$mir = $pairs[0];
				my $start_seed = $pairs[1];
				next unless $start_seed == $fixed_start;				
				
				my $struc = $all_data{"${seed}_$mir"};

				$calculations{"${seed};${mir};${target};${start_seed};${pos};$struc"} = $.;									
			}
		}
	}
	close $file or die "$!: $sites";
		
	return \%calculations
}

sub send_2ry{								# Send the thermodynamic calculations (whole cDNA)
	my ($calculations,$cdnas,$three_utrs,$mirnas);
	   ($calculations,$cdnas,$three_utrs,$mirnas,@_) = @_;
	
	my @calcs = keys %$calculations;
	my $cores = 5;
	my $group = 0;
	my $bin = int((scalar @calcs)/$cores);
	my @packs;
	my ($file,$part);
	
	for my $i(0..$#calcs){

		if ($i == $group * $bin){
			close $file if $file;
			$group++;
			open $file, '>',"constraint$group" or die "$!: constraint$group"; 	#constraints file for RNAcofold
			open $part, '>',"partners$group" or die "$!: partners$group";		#miRNAs and cDNAs namefile
			push(@{$packs[$group-1]},"constraint$group","res_constraint$group")
		}
		
		my ($seed,$mir,$target,$start_seed,$pos,$struc) = split(/;/,$calcs[$i]);
		my $seq;

		if($cdnas->{$target}){
			$seq = $cdnas->{$target} 
		}else{
			$seq = $three_utrs->{$target}
		}

		my $constraint = dot_bracket($seed,$mirnas->{$mir},$start_seed,$pos,$seq,$struc);
		
		print $file $mirnas->{$mir}.'&'.$seq,"\n",$constraint;
		print $part $seed,"\t",$mir,"\t",$target,"\t",$pos,"\t";
				
		my $struc_scalar = join('x',map{$_ ? $_ : ()}@$struc);
		unless($struc_scalar){
			print $part "nowm\n"
		}elsif($struc_scalar =~ /:/ && $struc_scalar =~ /\*/){
			print $part "wm\n"
		}elsif($struc_scalar !~ /:/){
			print $part "m\n"
		}elsif($struc_scalar =~ /\*/){
			print $part "w\n"
		}
		
	}
	close $file or die "$!: constraint$group";
	close $part or die "$!: partners$group";
	
	my $sons = 1;
	for my $procs(0..$#packs){		
		if (fork() == 0){
			my $container = shift(@{$packs[$procs]});
			my $results = shift(@{$packs[$procs]});
			system "RNAcofold --noPS -p0 --noconv -C < $container > $results";
			exit();
		}else{
			$sons++
		}
	}
	wait for 1..$sons;
}


sub send_2ry_70{								# Send the thermodynamic calculations (w/ 70nt-long flanking sequences)
	my ($calculations,$cdnas,$three_utrs,$mirnas);
	   ($calculations,$cdnas,$three_utrs,$mirnas,@_) = @_;
	
	my @calcs = keys %$calculations;
	my $cores = 5;
	my $group = 0;
	my $bin = int((scalar @calcs)/$cores);
	my @packs;
	my ($file,$part);
	
	for my $i(0..$#calcs){

		if ($i == $group * $bin){
			close $file if $file;
			$group++;
			open $file, '>',"constraint$group" or die "$!: constraint$group"; 	#Constraints file for RNAcofold
			open $part, '>',"partners$group" or die "$!: partners$group";		#miRNAs and cDNAs namefile
			push(@{$packs[$group-1]},"constraint$group","res_constraint$group")
		}
		
		my ($seed,$mir,$target,$start_seed,$pos,$struc) = split(/;/,$calcs[$i]);
		my $seq;

		if($cdnas->{$target}){
			$seq = $cdnas->{$target} 
		}else{
			$seq = $three_utrs->{$target}
		}
		my ($constraint,$site70) = dot_bracket_70($seed,$mirnas->{$mir},$start_seed,$pos,$seq,$struc);
		
		print $file $mirnas->{$mir}.'&'.$site70,"\n",$constraint;
		print $part $seed,"\t",$mir,"\t",$target,"\t",$pos,"\t";
				
		my $struc_scalar = join('x',map{$_ ? $_ : ()}@$struc);
		unless($struc_scalar){
			print $part "nowm\n"
		}elsif($struc_scalar =~ /:/ && $struc_scalar =~ /\*/){
			print $part "wm\n"
		}elsif($struc_scalar !~ /:/){
			print $part "m\n"
		}elsif($struc_scalar =~ /\*/){
			print $part "w\n"
		}
		
	}
	close $file or die "$!: constraint$group";
	close $part or die "$!: partners$group";
	
	my $sons = 1;
	for my $procs(0..$#packs){		
		if (fork() == 0){
			my $container = shift(@{$packs[$procs]});
			my $results = shift(@{$packs[$procs]});
			system "RNAcofold --noPS -p0 --noconv -C < $container > $results";
			exit();
		}else{
			$sons++
		}
	}
	wait for 1..$sons;
}


sub do_distance{												#Filtering sites by their relative position in the 3'UTR and expression pattern
	my($dgs,$trigger,$pos3utr,$long_3utr,$label);
	($dgs,$trigger,$pos3utr,$long_3utr,@_) = @_;			
	
	$label = 'miRexpr_';
	my(@keys,@vals);
	if($trigger eq 'si'){
		$label .= 'dist';
		my @pos_filtered;
		foreach my $name(keys %$dgs){
			if ($pos3utr->{$name} > 15){						#15 nt's filter (Grimson et al., 2007)
				push(@pos_filtered,$name);
				pop (@pos_filtered) if ($long_3utr->{$name} >= 1300 && $pos3utr->{$name} >= $long_3utr->{$name}/10 * 4 && $pos3utr->{$name} <= $long_3utr->{$name}/10 * 6);	
			}
		}		
		@keys = sort { $dgs->{$b} <=> $dgs->{$a} } @pos_filtered;
		@vals = @$dgs{@keys};

	}else{
		$label .= 'nodist';
		@keys = sort { $dgs->{$b} <=> $dgs->{$a} } keys(%$dgs);
		@vals = @$dgs{@keys};
	}
	
	my @data = (
		[@keys],
		[@vals]
	);	
	return \@data,$label
	
}

sub summary{							#Do the combination of the outputs of prediction algorithms
	
	my $mir_ver = pop(@_);
	my $ens_ver = pop(@_);
	my $sp = pop(@_);
	my $files = pop(@_);
	
	my($program,%sites,%progs);
		
	open my $coinc,">","Summary_3utr_${sp}_${ens_ver}_${mir_ver}" or die "$!: Summary_3utr_${sp}_${ens_ver}_${mir_ver}";
	print $coinc "miRNA\t3'UTR\tStart\tEnd\tSeedmatch\tSeedmatch_type\t3'Pairing\t#site(transc)\t#sites(gene)\tPrograms\n";
		
	foreach my $file(@{$files}){
		print "Analyzing $file\n";sleep(1);
		$program = $1 if $file =~ /^(\w+)_SITES/;
		open my $handle, "<", $file or die "$!: $file";
		while(<$handle>){
			next if $.== 1;
			my @line= split(/\t/);
			my $mirna = $line[1];			
			(my $gene, my $transc) = (split /\|/,$line[2])[0,1];
			$sites{$mirna}{$gene}{$transc}{"${line[3]}\t$line[4]"} = "$line[5]\t$line[6]\t$line[7]";;	#%sites{miRNA}{target_start_end}=Seed-match\tSeed-match_type\t3'Pairing
			
			push(@{$progs{$mirna}{$gene}{$transc}{"${line[3]}\t$line[4]"}}, $program);
		}
		close $handle or die "$!: $file";
	}
	
	my %sitenumber; 
	foreach my $mirna(sort keys %progs){
		foreach my $gene(sort keys %{$progs{$mirna}}){
			foreach my $transc(keys %{$progs{$mirna}{$gene}}){
				foreach my $coords(keys %{$progs{$mirna}{$gene}{$transc}}){
					if ($#{$progs{$mirna}{$gene}{$transc}{$coords}} > 0){
						$sitenumber{$mirna}{$gene}++;
						$sitenumber{$mirna}{$transc}++;					
					}
				}
			}
		}
	}
	
	foreach my $mirna(sort keys %progs){
		foreach my $gene(sort keys %{$progs{$mirna}}){
			foreach my $transc(keys %{$progs{$mirna}{$gene}}){
				foreach my $coords(keys %{$progs{$mirna}{$gene}{$transc}}){
					
					@{$progs{$mirna}{$gene}{$transc}{$coords}} = uniq(@{$progs{$mirna}{$gene}{$transc}{$coords}});
					
					if ($#{$progs{$mirna}{$gene}{$transc}{$coords}} > 0){
						chomp(@{$progs{$mirna}{$gene}{$transc}{$coords}});
						
						my $a = $sites{$mirna}{$gene}{$transc}{$coords};
						print $coinc "$mirna\t${gene}|$transc\t$coords\t$a\t$sitenumber{$mirna}{$transc}\t$sitenumber{$mirna}{$gene}\t";

						
						foreach my $i (0..$#{$progs{$mirna}{$gene}{$transc}{$coords}}){
							print $coinc "$progs{$mirna}{$gene}{$transc}{$coords}[$i]";
							if ($i == $#{$progs{$mirna}{$gene}{$transc}{$coords}}){
								print $coinc "\n"
							}else{
								print $coinc ","
							}
						}
					}
				}
			}
		}
	}
	close $coinc or die "$!: Summary_3utr_${sp}_${ens_ver}_${mir_ver}";	
}

sub create_match{					#Draw the structure of the seed:seedmatch hybrid
	my @seq;
	$seq[0] = shift(@_);
	$seq[2] = shift(@_);
	
	my @temp_utr= split(//,$seq[2]);
	my @temp_mir= split(//,$seq[0]);
  
	my @temp_match;
	for (my $j=0;$j<=$#temp_mir;$j++){
		if (($temp_utr[$j] eq "A" && $temp_mir[$j]=~ /U|T/) || ($temp_mir[$j] eq "A" && $temp_utr[$j]=~ /U|T/) || ($temp_utr[$j] eq "C" && $temp_mir[$j] eq "G") || ($temp_utr[$j] eq "G" && $temp_mir[$j] eq "C")){		#Si hay nt's complementarios escribe un "|"
			$temp_match[$j]= "|";
		}elsif(($temp_utr[$j] eq "G" && $temp_mir[$j]=~ /U|T/) || ($temp_mir[$j] eq "G" && $temp_utr[$j]=~ /U|T/)){#If G facing U
			$temp_match[$j]= ":";
		}else{									#If blank space
			$temp_match[$j]= " ";
		}
	}
	$seq[1]= join('',@temp_match);
	
	return $seq[1];	
}

sub smatch_beg_end_findtar{							#FindTar. Determine start and end positions of the seedmatch
	my $positions = shift(@_);
	my $to_get_pos2 = shift(@_);
	my $store_mir = shift(@_);
	
	for my $q(0..$#{$to_get_pos2}){                #Seedmatch start
		if ($to_get_pos2->[$q] =~ /[A-U]/ && $store_mir->[$q] eq '|'){
			$positions->[0] += $q;last;
		}
	}
	my $end_cds;
	for (my $q=$#{$to_get_pos2};$q>=0;$q--){                #Seedmatch end
		if ($to_get_pos2->[$q] =~ /[A-U]/ && $store_mir->[$q] eq '|'){
			$end_cds = $q+$positions->[0];last;
		}
	}	
	
	return ($positions->[0],$end_cds);
}



sub smatch_beg_end_rnah{							#RNAHybrid. Determine start and end positions of the seedmatch

	my $to_get_pos2 = pop(@_);
	my $to_get_pos1 = pop(@_);
	my $positions = pop(@_);	
	my $n=0;
	my $blank = 1;
	
	for my $q(0..$#{$to_get_pos1}){		#Seedmatch start
		if ($to_get_pos1->[$q] eq " " && $to_get_pos2->[$q] =~ /[A-U]/){
			$positions->[0] += $q-$n;last;
		}elsif($to_get_pos1->[$q] eq " " && $to_get_pos2->[$q] eq " "){
			$n++ 
		}
	}
	my $end_cds;
	for (my $q=$#{$to_get_pos1};$q>=0;$q--){ 		#Seedmatch end
		if ($to_get_pos1->[$q] eq " " && $to_get_pos2->[$q] =~ /[A-U]/){
			$end_cds = $q+$positions->[0];last;
		}
	}
	for my $q(0..$#{$to_get_pos1}){
		if($to_get_pos1->[$q] eq " " && $to_get_pos2->[$q] eq " "){
			$blank++;
		}
	}$end_cds -= $blank;

	return ($positions->[0],$end_cds);

}


sub smatch_beg_end_miranda{								#miRanda. Determine start and end positions of the seedmatch

	my $to_get_pos2 = pop(@_);
	my $positions = pop(@_);

	for my $q(0..$#{$to_get_pos2}){                #Seedmatch start
	  if ($to_get_pos2->[$q] =~  /[A-U]/ ){
		$positions->[0] += $q;last;
	  }   
	}
	
	my $end_cds;
	for (my $q=$#{$to_get_pos2};$q>=0;$q--){                #Seedmatch end
	  if ($to_get_pos2->[$q] =~ /[A-U]/){
		$end_cds = $q+$positions->[0];last
	  }   
	}

	return ($positions->[0],$end_cds);
}

sub parse_pita {										#Parse the output of PITA algorithm
	
	my($sm_type,@line,$target,$mirna,$start,$end,$dg,$site);

	print "Reading $_[0]\n";
	open my $fh_res, ">","pita_SITES_3utr" or die "$!: pita_SITES_3utr";
	print $fh_res "Line#\tmiRNA\t3utr\tStart\tEnd\tSeedmatch\tSeedmatch_type\t3'Pairing\tdG\tTotal_Score\tProgram\n";		
		
	open my $file,"<",$_[0] or die "$!:$_[0]";
	while(<$file>){
		next if $.==1;		
		@line = split /\t/;
		$target=$line[0];
		$mirna=$line[1];
		$site = "3utr_${mirna}_$target";		
								
		next if $line[4] =~ /:1/;		
		
		$start=$line[3];
		$end= $line[2]-1;
		if ($end-$start == 5){
			$sm_type = "6mer"
		}elsif ($end-$start == 6){
			$sm_type = "7mer"
		}elsif ($end-$start == 7){
			$sm_type = "8mer"
		}		
		
		$dg=$line[$#line];
		chomp($dg);
		print $fh_res "$.\t$mirna\t$target\t$start\t$end\t\t$sm_type\t?\t$dg\t\tpita\n";
	}
	close $file or die "$!: $_[0]";
	close $fh_res or die "$!: pita_SITES_3utr"
}

sub parse_miranda{										#Parse the output of miRanda algorithm
	
	my $program = shift(@_);
	my (@seq,@positions,$mirna,$target,$tot_score,$dg,$site,@transc);
	
	print "Reading $_[0]\n";	
	open my $fh_res, ">","${program}_SITES_3utr" or die "$!: ${program}_SITES_3utr";
	print $fh_res "Line#\tmiRNA\tTarget\tStart\tEnd\tSeedmatch\tSeedmatch_type\t3'Pairing\tdG\tTotal_Score\tProgram\n";

    open my $fh_struc, ">", "${program}_3utr_NON-CANONICAL_STRUCT" or die "$!: ${program}_NON-CANONICAL_STRUCT";

    open my $fh_rare, ">", "${program}_3utr_NON-CANONICAL_SITES" or die "$!: ${program}_NON-CANONICAL_SITES";
    print $fh_rare "Line#\tmiRNA\t3utr\tMatch_pos\tMismatch_pos\tWobble_pos\tdG\tTotal_Score\tProgram\n";
	
	my @temp;
	open my $file, "<", $_[0] or die "$!:$_[0]";
	while(<$file>){
		if (/^\s+Query:\s+3\'\s(\S+)\s5\'/){
			$seq[0] = (split, /\s+/)[2];
			$seq[0] = uc $seq[0];
		}elsif(/^\s+Ref:\s+5\'\s(\S+)\s3\'/){
			$seq[2] = (split, /\s+/)[2];
			$seq[2] = uc $seq[2];
			
			$seq[1] = create_match($seq[0],$seq[2]);		#Crear pairing scheme line
		}elsif(/^>\w+/){			
			($mirna,$target,$tot_score,$dg,$temp[0],$positions[0]) = split(/\t/);
			$mirna =~ s/>//;
			$positions[0] = (split /\s/,$positions[0])[0];
			
			$site= "${mirna}_$target";
			my @to_get_pos2= split(//,$seq[2]);
			
			
			for(my $j=0;$j<=$#to_get_pos2;$j++){         #Get the seeedmatch
				if ($to_get_pos2[$j]=~ /\w/){
				 push(@transc,$to_get_pos2[$j]);
				}
			}
			for (my $j=1;$j<=$#transc;$j++){              #Fill in @positions with the seed positions
				$positions[$j]= $positions[$j-1]+1;
			}
			undef(@transc);
			my $wobbs = $seq[1] =~ tr/://;
			my $mirseq = $seq[1] =~ tr/ //;
			$mirseq -= $seq[0] =~ tr/-//;
			
			for my $b(0..$#to_get_pos2){
				$to_get_pos2[$b] = "\u$to_get_pos2[$b]"
			}
			$seq[2] = "\U$seq[2]";
			
			($positions[0],my $end_cds) =  smatch_beg_end_miranda(\@positions,\@to_get_pos2);
			my @data= site_type_from_struct($program,\@seq,\@positions,\@to_get_pos2);
				
			
			if (@data){
				if ($data[0] && !$data[1]){
					my @print = split(/\t/,$data[0]);
					print $fh_res "$.\t$mirna\t$target\t";
					for my $l(0..4){print $fh_res "$print[$l]\t"}
					print $fh_res "$dg\t$tot_score\t$program\n";		
				}
				elsif ($data[1]){
					for my $i(0..1){	
						my @print = split(/\t/,$data[$i]);
						if($i == 0){	
							print $fh_rare "$.\t$mirna\t$target\t";
							for my $l(0..2){print $fh_rare "$print[$l]\t"}			
							print $fh_rare "$program\n";
						}else{
							print $fh_struc "Line#: $.\tmiRNA:${mirna}\tTarget:${target}\t";
							for my $l(0..2){print $fh_struc "$print[$l]\t"}
							print $fh_struc "dG:$dg,Program: $program\n$seq[0]\n$seq[1]\n$seq[2]\n\n"
						}	
					}
				}
				undef(@data);
			}			
		}
	}
	close $file or die "$!:$_[0]";
	close $fh_res or die "$: ${program}_SITES_3utr";
	close $fh_rare or die "$!: ${program}_NON-CANONICAL_SITES";
	close $fh_struc or die "$!: ${program}_NON-CANONICAL_STRUC";
}


sub parse_targetscan { 								#Parse the output of TargetScan algorithm
	my (%site_nu);
	
	print "Reading $_[0]\n";
	open my $fh_res, ">","targetscan_SITES_3utr" or die "$!: targetscan_SITES_3utr";
	print $fh_res "Line#\tmiRNA\t3utr\tStart\tEnd\tSeedmatch\tSeedmatch_type\t3'Pairing\tdG\tTotal_Score\tProgram\n";		
	
	open my $file,"<", $_[0] or die "$!: $_[0]";
	while(<$file>){
		next if $.== 1;
		s/>//;
		my @line = split /\t/;
		my $target= $line[0]; my $mirna= $line[1]; my $start= $line[3]; my $end= $line[4]; my $type= $line[8];
		my $site= "${mirna}_$target";			
		print $fh_res "$.\t$mirna\t$target\t$start\t$end\t\t$type\t?\t\t\ttargetscan\n";
	}
	close $file or die "$!:$_[0]";
	close $fh_res or die "$!:targetscan_SITES_3utr";
   
}

sub site_type_from_struct{	#Analize the miRNA-site structure and determine seedmatch type
	my(@mir,@match,@mrna,$i,@mismatch_loc,@wobble_loc,@match_loc,%match_pos,$end_match_row,$salida,$un_site,$un_str);
	my $match_nu = my $mismatch_nu = my $wobble_nu = my $match_pos = my $wobble_pos = my $mismatch_pos = 0;
	my($mismatch_loc,$wobble_loc,$match_loc,@in_a_row);
	my $program = shift(@_);
	my $seq = shift(@_);
	my $positions = shift(@_);
	my $to_get_pos2= pop(@_);
	my $to_get_pos1= shift(@_);	

																	# Longest stretch of matches(*in_a_row)
	@mir= split(//,$seq->[0]);
	@match= split(//,$seq->[1]);
	@mrna= split(//,$seq->[2]);
	
	my $in_a_row = 0;
	for ($i=0;$i<=$#match;$i++){          # Find the position of the first match, mismatch or wobble
		if ($match[$i] eq "|"){
			$match_nu++;
			$match_pos= $#match-$i+1;
			$match_loc.= "$match_pos-";
			$in_a_row++;
			$end_match_row=$i;
			if($i== $#match){
				push(@in_a_row,$in_a_row);undef($in_a_row)
			}
		}elsif ($match[$i] eq ":"){
			$wobble_nu++;
			$wobble_pos= $#match-$i+1;
			$wobble_loc.= "$wobble_pos-";
			if($in_a_row){
				push(@in_a_row,$in_a_row);undef($in_a_row)
			}
		}elsif ($mir[$i] eq "-" || $mrna[$i] eq "-" || $match[$i] eq " "){
			$mismatch_nu++;
			$mismatch_pos= $#match-$i+1;
			$mismatch_loc.= "$mismatch_pos-";
			if($in_a_row){
				push(@in_a_row,$in_a_row);undef($in_a_row)
			}
		}
	}

	my $beg_match_row= $end_match_row-$in_a_row[-1]+1;

	if ($wobble_loc){chop($wobble_loc)}else{$wobble_loc=" "};if ($match_loc){chop($match_loc)}else{$match_loc=" "};
	if ($mismatch_loc){chop($mismatch_loc)}else{$mismatch_loc=" "};
	my @in_a_row_ok= @in_a_row;
	$in_a_row = $in_a_row_ok[-1];

	my $compl_pair = "-";

	if ($in_a_row >= 6){
		my $j=0;						
		if($program =~ /rnahybrid/){
			for($i=0;$i<=$#{$to_get_pos2};$i++){
				if ($to_get_pos2->[$i] =~ /[A-U]/ || $to_get_pos1->[$i] =~ /[A-U]/){
					$match_pos{$i}= $positions->[$j] if ($match[$i] eq "|");          
					$j++
				}
			}
		}else{
			for($i=0;$i<=$#{$to_get_pos2};$i++){
				if ($to_get_pos2->[$i] =~ /[A-U]/){
					$match_pos{$i}= $positions->[$j] if ($match[$i] eq "|");
					$j++
				}   
			}   
		}			
		if(($mismatch_pos >= 7 || $mismatch_pos <= 1) && ($wobble_pos >= 7 || $wobble_pos<= 1) && $match_pos<= 2){
		
			# 3' Complementary pairing
			
			if($match[$#match-13] eq "|" && $match[$#match-14] eq "|"){	#Complementarity with miRNA nts 14 and 15
				if($match[$#match-11] eq "|" && $match[$#match-12] eq "|" && $match[$#match-15] eq "|" && $match[$#match-16] eq "|"){
					$compl_pair = "12-17";
				}elsif($match[$#match-11] eq "|" && $match[$#match-12] eq "|" && $match[$#match-15] eq "|"){
					$compl_pair = "12-16";
				}elsif($match[$#match-12] eq "|" && $match[$#match-15] eq "|" && $match[$#match-16] eq "|"){
					$compl_pair = "13-17";
				}elsif($match[$#match-11] eq "|" && $match[$#match-12] eq "|"){
					$compl_pair = "12-15";
				}elsif($match[$#match-12] eq "|" && $match[$#match-15] eq "|"){
					$compl_pair = "13-16";
				}elsif($match[$#match-15] eq "|" && $match[$#match-16] eq "|"){
					$compl_pair = "14-17";
				}
			}
			if($in_a_row== 6 && $match[$#match-1] eq "|" && $match[$#match-2] eq "|" && $match[$#match-3] eq "|" && $match[$#match-4] eq "|" && $match[$#match-5] eq "|"){
				if ($match[$#match] eq "|"){
					my $smseq;
					for (my $r=5;$r>=0;$r--){
						$smseq .= $mrna[$#mrna-$r];				
					}
					$salida= "$match_pos{$beg_match_row}\t$match_pos{$end_match_row}\t$smseq\t6mer\t$compl_pair";
					return $salida;
				}elsif ($match[$#match-6] eq "|" && $mrna[$#mrna] ne "A"){
					my $smseq;
					for (my $r=6;$r>=1;$r--){
						$smseq .= $mrna[$#mrna-$r];				
					}
					$salida= "$match_pos{$beg_match_row}\t$match_pos{$end_match_row}\t$smseq\t6mer\t$compl_pair";
					return $salida;
				}elsif ($match[$#match-6] eq "|" && $mrna[$#mrna] eq "A"){
					$match_pos{$end_match_row}++;
					my $smseq;
					for (my $r=6;$r>=1;$r--){
						$smseq .= $mrna[$#mrna-$r];				
					}
					$salida= "$match_pos{$beg_match_row}\t$match_pos{$end_match_row}\t$smseq\t7mer-A1\t$compl_pair";
					return $salida;
				}
			}elsif($in_a_row== 7 && ($match[$#match-1] eq "|" && $match[$#match-2] eq "|" && $match[$#match-3] eq "|" && $match[$#match-4] eq "|" && $match[$#match-5] eq "|" && $match[$#match-6] eq "|")){
				if ($match[$#match] eq "|"){
					my $smseq;
					for (my $r=6;$r>=0;$r--){
						$smseq .= $mrna[$#mrna-$r];				
					}
					$salida= "$match_pos{$beg_match_row}\t$match_pos{$end_match_row}\t$smseq\t7mer-m8\t$compl_pair";
					return $salida;
				}elsif ($match[$#match-7] eq "|"){
					my $smseq;
					for (my $r=7;$r>=1;$r--){
						$smseq .= $mrna[$#mrna-$r];				
					}
					$salida= "$match_pos{$beg_match_row}\t$match_pos{$end_match_row}\t$smseq\t7mer-m8\t$compl_pair";
					return $salida;
				}
			}elsif ($in_a_row==8 && $match[$#match-1] eq "|" && $match[$#match-2] eq "|" && $match[$#match-3] eq "|" && $match[$#match-4] eq "|" && $match[$#match-5] eq "|" && $match[$#match-6] eq "|" && $match[$match_pos-7] eq "|"){
				if ($match[$#match] eq "|"){
					my $smseq;
					for (my $r=7;$r>=0;$r--){
						$smseq .= $mrna[$#mrna-$r];				
					}
					$salida= "$match_pos{$beg_match_row}\t$match_pos{$end_match_row}\t$smseq\t8mer\t$compl_pair";
					return $salida;
				}elsif ($match[$#match-8] eq "|"){
					my $smseq;
					for (my $r=8;$r>=1;$r--){
						$smseq .= $mrna[$#mrna-$r];				
					}
					$salida= "$match_pos{$beg_match_row}\t$match_pos{$end_match_row}\t$smseq\t8mer\t$compl_pair";
					return $salida;
				}
			}elsif($in_a_row>8){
				$salida= "$match_pos{$beg_match_row}\t$match_pos{$end_match_row}\t$smseq\tNON-C\t$compl_pair";
				$un_site= "$match_loc\t$mismatch_loc\t$wobble_loc";
				$un_str= "No. matches:${match_nu}\tNo. mismatches:${mismatch_nu}\tNo. wobbles:${wobble_nu}\t";
				return ($un_site,$un_str);
			}
		}	# General conditions of mismatches, wobbles and match_position
	}	# If more than 5 matches in a row
}

sub parse_rnahybrid {     	#Parse the output of RNAHybrid algorithm
	my $program = shift(@_);
	my $a = pop(@_);
	my $tot_score= " ";
	my ($mirseq,@seq,@positions,$mirna,$target);
	
	print "Reading $_[0]\n";	
	open my $fh_res, ">","${program}_SITES_3utr" or die "$!: ${program}_SITES_3utr";
	print $fh_res "Line#\tmiRNA\tTarget\tStart\tEnd\tSeedmatch\tSeedmatch_type\t3'Pairing\tdG[kCal/mol]\tTotal_Score\tProgram\n";

    open my $fh_struc, ">", "${program}_3utr_NON-CANONICAL_STRUCT" or die "$!: ${program}_NON-CANONICAL_STRUCT";

    open my $fh_rare, ">", "${program}_3utr_NON-CANONICAL_SITES" or die "$!: ${program}_NON-CANONICAL_SITES";
    print $fh_rare "Line#\tmiRNA\t3utr\tMatch_pos\tMismatch_pos\tWobble_pos\tdG\tTotal_Score\tProgram\n";

	my (@to_get_pos1,$dg,@mirseq);
	open my $file, "<", $_[0] or die "$!: $_[0]";
	while(<$file>){
		next if /^$/;
		my @line = split /\s/;
		if (/^target:/){
			$target = $line[1];
		}elsif(/^miRNA :/){
			$mirna = $line[2];
		}elsif(/^mfe/){
			$dg = $line[1];
		}elsif(/^position/){
			undef(@positions);
			$positions[0] = $line[2];
		}elsif(/^target\s5\'\s(.*)\s3\'/){
			@to_get_pos1= split(//,$1);
		}elsif (/^\s{10}(.*?)\s{3}$/){			                          
			$seq[$a]= $1;
			$a-= 2;
			next if $a == 0;
		}elsif(/^miRNA\s+3\'\s(.*)\s5\'/){
			@mirseq = split(/\s+/,$1);
			$mirseq = join("",@mirseq);$a-= 2;
		}
		if ($a == -4){		
			$a = 2;
			my $site= "${mirna}_$target";

			$seq[1] = create_match($seq[0],$seq[2]);	#Create the pairing scheme line
			
			my @to_get_pos2= split(//,$seq[2]);
			
			my @transc;
			for(my $j=0;$j<=$#to_get_pos2;$j++){		#Get the seedmatch sequence
				if ($to_get_pos1[$j]=~ /[A-U]/){
					push(@transc,$to_get_pos1[$j]);
				}elsif($to_get_pos2[$j]=~ /[A-U]/){
					push(@transc,$to_get_pos2[$j]);
				}
			}
			for (my $j=1;$j<=$#transc;$j++){		#Fill in @positions with the positions of the site
				$positions[$j]= $positions[$j-1]+1
			}
				
			my($length_site,@length_site);
			for (my $j=0;$j<=$#to_get_pos2;$j++){
				if ($to_get_pos2[$j]=~/[A-U]/){
					$length_site.= $to_get_pos2[$j]
				}else{
					push(@length_site,$length_site) if ($length_site);undef($length_site)
				}
				push(@length_site,$length_site) if ($length_site);
			}
			if (@length_site && length($length_site[$#length_site]) > 5) { #If seedmatch is longer than 5 nt's

				my $wobbs = $seq[1] =~ tr/://;
				($positions[0],my $end_cds) =  smatch_beg_end_rnah(\@positions,\@to_get_pos1,\@to_get_pos2);
				my @data = site_type_from_struct($program,\@seq,\@positions,\@to_get_pos1,\@to_get_pos2);
					
				if (@data){
					if ($data[0] && !$data[1]){
						my @print = split(/\t/,$data[0]);
						print $fh_res "$.\t$mirna\t$target\t";
						for my $l(0..4){print $fh_res "$print[$l]\t"}
						print $fh_res "$dg\t$tot_score\t$program\n";					
					}
					elsif ($data[1]){
						for my $i(0..1){	
							my @print = split(/\t/,$data[$i]);
							if($i == 0){	
								print $fh_rare "$.\t$mirna\t$target\t";
								for my $l(0..2){print $fh_rare "$print[$l]\t"}			
								print $fh_rare "$program\n";
							}else{
								print $fh_struc "line#: $.\tmiRNA:${mirna}\tTarget:${target}\t";
								for my $l(0..2){print $fh_struc "$print[$l]\t"}
								print $fh_struc "dG:$dg,Program: $program\n$seq[0]\n$seq[1]\n$seq[2]\n\n"
							}	
						}
					}
					undef(@data);
				}
			}
		}
	}
	close $file or die "$!: $_[0]";
	close $fh_res or die "$!: rnahybrid_SITES_3utr";
    close $fh_rare or die "$!: rnahybrid_3utr_NON-CANONICAL_SITES";
    close $fh_struc or die "$!: rnahybrid_3utr_NON-CANONICAL_STRUCT"
}

sub parse_findtar{							#Parse the output of FindTar algorithm

	my $program = shift(@_);

	my (@seq,@positions,$mirna,$target,$tot_score,$dg,$site,@transc,@line);
	
	print "Reading $_[0]\n";	
	open my $fh_res, ">","${program}_SITES_3utr" or die "$!: ${program}_SITES_3utr";
	print $fh_res "Line#\tmiRNA\tTarget\tStart\tEnd\tSeedmatch\tSeedmatch_type\t3'Pairing\tdG\tTotal_Score\tProgram\n";

    open my $fh_struc, ">", "${program}_3utr_NON-CANONICAL_STRUCT" or die "$!: ${program}_NON-CANONICAL_STRUCT";

    open my $fh_rare, ">", "${program}_3utr_NON-CANONICAL_SITES" or die "$!: ${program}_NON-CANONICAL_SITES";
    print $fh_rare "Line#\tmiRNA\tTarget\tMatch_pos\tMismatch_pos\tWobble_pos\tdG\tTotal_Score\tProgram\n";

	
	open my $file, $_[0] or die "$!: $_[0]";
	while(<$file>){
		next if /^(#|======)/;	    
		if(/findtar\sfinal\sstructure:/){
			for(1..4){<$file>};
		}else{
			@line = split /\s+/;
			if (/^\S+\s\S+$/ && !/:/){
				$mirna = $line[0];
				$target = $line[1];				
				$site = "${mirna}_$target";
			}elsif(/^Position/){
				$positions[0]= $line[1];chop($positions[0]);
				$tot_score = $line[4];
				$dg = $line[6];
			}elsif(/^miRNA/){
				$seq[0]= $line[1];
			}elsif(/^\s{7}/){
				$seq[1]= $line[1];
				$seq[1]=~ s/\*/ /g;
			}elsif(/3\'UTR:/){
				$seq[2]= $line[1];
				my @to_get_pos2= split(//,$seq[2]);
				undef(@transc);   
			
				for(my $j=0;$j<=$#to_get_pos2;$j++){         #Get the seedmatch
				if ($to_get_pos2[$j]=~ /[A-U]/){
					push(@transc,$to_get_pos2[$j]);
				}
				}
				for (my $j=1;$j<=$#transc;$j++){              #Fill in @positions with the seed positions
					$positions[$j]= $positions[$j-1]+1;
				}
				my $wobbs = $seq[1] =~ tr/://;
				my $mirseq = $seq[1] =~ tr/ //;
				$mirseq -= $seq[0] =~ tr/-//;
				
				my @data = site_type_from_struct($program,\@seq,\@positions,\@to_get_pos2);
				
				if (@data){
					if ($data[0] && !$data[1]){
						my @print = split(/\t/,$data[0]);
						print $fh_res "$.\t$mirna\t$target\t";
						for my $l(0..4){print $fh_res "$print[$l]\t"}
						print $fh_res "$dg\t$tot_score\t$program\n";			
					}
					elsif ($data[1]){
						for my $i(0..1){	
							my @print = split(/\t/,$data[$i]);
							if($i == 0){	
								print $fh_rare "$.\t$mirna\t$target\t";
								for my $l(0..2){print $fh_rare "$print[$l]\t"}			
								print $fh_rare "$program\n";
							}else{
								print $fh_struc "Line#: $.\tmiRNA:${mirna}\tTarget:${target}\t";
								for my $l(0..2){print $fh_struc "$print[$l]\t"}
								print $fh_struc "dG:$dg,Program: $program\n$seq[0]\n$seq[1]\n$seq[2]\n\n"
							}	
						}
					}
					undef(@data);
				}
			}
		} 
	}
	close $file or die "$!: $_[0]";
	close $fh_res or die "$: ${program}_SITES_3utr";
	close $fh_rare or die "$!: ${program}_NON-CANONICAL_SITES";
	close $fh_struc or die "$!: ${program}_NON-CANONICAL_STRUC";	
}

1
