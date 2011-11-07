package exonerate;
use strict;
use Data::Dumper;
my $exhaustive = 1;
my $debug = 1;

sub new {#{{{
	print join(" ", (caller(0))[0..3]), "\n" if $debug;
	my $self_tmp = [];
	my $self;
	my ($class, $dna, $prot, $path) = @_;
	$path = '/tmp' unless defined $path;
	print "Creating new exonerate object\n";
	my $pid = $$;
	my $protname = $pid . '_protein';
	my $dnaname = $pid . '_dna';

	# some exonerate options, to be configured later (TODO)
	my $exonerate_model = $exhaustive ? 'protein2genome:bestfit' : 'protein2genome';
	my $exonerate_exhaustive = $exhaustive ? '--exhaustive 2> /dev/null' : '';
	# roll your own output for exonerate
	my $exonerate_ryo = "Score: %s\n%V\n>%qi_%ti_[%tcb-%tce]_cdna\n%tcs//\n>%qi[%qab:%qae]_query\n%qas//\n>%ti[%tab:%tae]_target\n%tas//\n";
	my $exonerate_cmd = "exonerate --ryo '$exonerate_ryo' $exonerate_exhaustive --model $exonerate_model --verbose 0 --showalignment no --showvulgar no";
	my $translate_cmd = "fastatranslate -F 1";

	# print the two seqs to files in path
	open(my $dnafh, ">$path/$dnaname") or die "Fatal: Could not open $path/$dnaname for writing: $!\n";
	print $dnafh ">$dnaname\n$dna";
	close $dnafh;
	open(my $protfh, ">$path/$protname") or die "Fatal: Could not open $path/$protname for writing: $!\n";
	print $protfh ">$protname\n$prot";
	close $protfh;	#and that's that

	# now run exonerate on the two sequences!
	# saving the whole output in an array
	print 'running: $exonerate_cmd ' . "$path/$protname" . ' ' . "$path/$dnaname" . " 2> /dev/null\n";
	$self_tmp = [`$exonerate_cmd $path/$protname $path/$dnaname 2> /dev/null`];
	if (scalar @{$self_tmp} == 0) {
		warn "Warning: Alignment program about to crash later: exonerate returned nothing. Exit code: $?\n" 
	}
	foreach (@$self_tmp) {	# remove all newlines from eol
		chomp;
	}
	$self->{gw} = $self_tmp;
	$self->{gw_count} = scalar @$self_tmp;
	$self->{nt_seq} = $dna;
	$self->{prot_seq} = $prot;
	$self->{protname} = $protname;
	$self->{dnaname} = $dnaname;
	$self->{get_indel} = 1;
	$self->{indels} = _GetIndels($self_tmp);
	$self->{tmpdir} = $path;
	print '$self:', Dumper($self) if $debug;

	bless($self, $class);

	print join(" ", (caller(0))[0..3]), ", leaving\n" if $debug;
	return $self;
}#}}}

#-------------------------------------------------- 
sub score {#{{{
	print join(" ", (caller(0))[0..3]), "\n" if $debug;
	my $self = shift;
	my $score;
	foreach my $line ($self->{gw}) {
		if ($line =~ /^Score: (\d+\.?\d+)/) {
			$score = $1;
			last;
		}
	}
	print join(" ", (caller(0))[0..3]), ", leaving\n" if $debug;
	return $score;
}#}}}

##################
#mp copypasta from hamstr
#mp this function is never called?!
#--------------------------------------------------
# sub protein {#{{{
# 	print join(" ", (caller(0))[0..3]), "\n" if $debug;
# 	my $self = shift;
# 	my $gw = $self->{gw};
# 	my $prot = '';
# 	for (my $i = 0; $i < @$gw; $i++) {
# 		if ($gw->[$i] =~ /^>.*_query$/) { #the protein seq starts
# 			my $count = 1;
# 			while ($gw->[$i+$count] ne '//') {
# 				my $protpart = $gw->[$i+$count];
# 				chomp $protpart;
# 				$prot .= $protpart;
# 				$count++;
# 			}
# 		}
# 		elsif (length $prot > 0) {
# 			last;
# 		}
# 	}
# 	print join(" ", (caller(0))[0..3]), ", leaving\n" if $debug;
# 	return($prot);
# }#}}}
#-------------------------------------------------- 

# Sub: translate_cdna
# Input: scalar string cDNA sequence
# Returns: scalar string translated cDNA sequence
sub translate_cdna {#{{{
	print join(" ", (caller(0))[0..3]), "\n" if $debug;
	my ($cdna, $self) = @_;
	my $tmpdir = $self->{tmpdir};
	my $cdnatmpfile = "$tmpdir/$$" . '_cdna.tmp';
	open( my $tmpfh, ">$cdnatmpfile" ) or die "Fatal: Could not open tmpfile $cdnatmpfile for writing: $!";
	print $tmpfh '>' . $self->{protname}, "\n", $cdna;	# cdna isn't enough, fastatranslate requires a valid fasta file including header
	close( $tmpfh ) or die "Fatal: Could not close tmpfile $cdnatmpfile: $!";
	print "running: fastatranslate -F 1 $cdnatmpfile\n" 
		if $debug;
	my $translate_result = [`fastatranslate -F 1 $cdnatmpfile`];
	die "Fatal: cDNA translation failed: $!" 
		unless (scalar @$translate_result > 0);
	warn "Warning: Alignmentprog about to crash later: fastatranslate died: $!"
		if scalar @{$translate_result} == 0;
	#--------------------------------------------------
	# unlink $cdnatmpfile	# delete temp file
	# 	or warn "Warning: Could not delete cDNA tmp file $cdnatmpfile: $!";
	#-------------------------------------------------- 

	shift @$translate_result;	# remove the first line (header)
	my $translated_cdna;
	# concatenate lines because some people (including me) like one-liner seqs better 
	foreach (@$translate_result) {
		chomp;
		last if /^>/;
		$translated_cdna .= $_;
	}
	print join(" ", (caller(0))[0..3]), ", leaving\n" if $debug;
	return $translated_cdna;
}#}}}

#mp copypasta from hamstr
sub translation {#{{{
	print join(" ", (caller(0))[0..3]), "\n" if $debug;
	my $self = shift;
	my $finish = 0;
	my $cdna_seq = '';
	my @transtmp;

	## step 1: extract the relevant info from the genewise output
	for (my $i = 0; $i < $self->{gw_count}; $i++) {
		#if ($self->{gw}->[$i] =~ />.*.tr/) {# a translated bit starts #mp !!regex will also match ">424324_bartralle.sp" because the . is not escaped :P 
		if ($self->{gw}->[$i] =~ />.*_cdna$/) {	#mp corrected regex
			print '$self->{gw}->[$i] is: ' , $self->{gw}->[$i], "\n" if $debug;
			while ($self->{gw}->[$i] !~ '//') {
				push @transtmp, $self->{gw}->[$i];
				$i++;
			}
			last; # end the for loop since nothing left to be done
		}
	}
	
	## step two: get the sequences
	my $count = -1;
	my $trans = [ ];
	for (my $i = 0; $i < scalar @transtmp; $i++) {
		if ($transtmp[$i] =~ /^>/) {
			$count++;
			$trans->[$count]->{seq} = ''; # initialize
			if ($transtmp[$i] =~ /.*\[([0-9]+):([0-9]+)\].*/) {
				$trans->[$count]->{start} = $1;
				$trans->[$count]->{end} = $2;
				}
		}
		else {
			$trans->[$count]->{seq} .= $transtmp[$i];
		}
	}

	## step 3: connect the fragments
	if (@$trans == 1) {
		$cdna_seq = $trans->[0]->{seq};
	}
	else {
		for (my $i = 0; $i < @$trans; $i++) {
			$cdna_seq .= $trans->[$i]->{seq};
			if ($i < (@$trans - 1)) {
				my $missing = $trans->[$i+1]->{start} - $trans->[$i]->{end} -1;	#mp never used 
				$cdna_seq .= 'X';
			}
		}
	}
	my $translated_seq = &translate_cdna($cdna_seq, $self);
	print join(" ", (caller(0))[0..3]), ", leaving\n" if $debug;
	return($translated_seq);
}#}}}

#mp copypasta from hamstr
sub codons {#{{{
	print join(" ", (caller(0))[0..3]), "\n" if $debug;
	my $self = shift;
	my $finish = 0;
	my $codon_seq = '';
	my @transtmp;

	## step 1: extract the relevant info from the exonerate output
	for (my $i = 0; $i < $self->{gw_count}; $i++) {
		if ($self->{gw}->[$i] =~ /^>.*_target$/) {# the codons set starts
			while ($self->{gw}->[$i] !~ '//') {
				push @transtmp, $self->{gw}->[$i];
				$i++;
			}
			last; # end the for loop since nothing left to be done
		}
	}
	
	## step two: get the sequences
	my $count = -1;
	my $trans = [ ];
	for (my $i = 0; $i < @transtmp; $i++) {
		if ($transtmp[$i] =~ />/) {
			$count++;
			$trans->[$count]->{seq} = ''; # initialize
			if ($transtmp[$i] =~ /.*\[([0-9]+):([0-9]+)\].*/) {	# this info needs to be in the exonerate output as well
				$trans->[$count]->{start} = $1;
				$trans->[$count]->{end} = $2;
			}
		}
		else {
			$transtmp[$i] =~ tr/a-z/A-Z/;	# upper-case all sequence lines
			$trans->[$count]->{seq} .= $transtmp[$i];
		}
	}	# $count never exceeds 1 because there is never more than one target seq in the exonerate output

	## step 3: connect the fragments
	if (@$trans == 1) {
		print '@$trans is 1, only one seq', "\n" if $debug;
		$codon_seq = $trans->[0]->{seq};
	}
	else {	# under what circumstances is there more than one target seq in the genewise/exonerate output?
		print '@$trans is not 1, more than one seq', "\n" if $debug;
		for (my $i = 0; $i < @$trans; $i++) {
			$codon_seq .= $trans->[$i]->{seq};
			if ($i < (@$trans - 1)) {
				my $indel = '';
				my $missing = $trans->[$i+1]->{start} - $trans->[$i]->{end} -1;

				## now decide whether the nts that did not got translated are masked by
				## 'N' or whether they will be represented as lower case letters
				if ($self->{get_indel}) {
					$indel = substr($self->{nt_seq}, $trans->[$i]->{end}, $missing);
					$indel =~ tr/A-Z/a-z/;
				}
				else {
					$indel = 'N' x $missing;
				}
				## now append gap characters until the frame is recovered. Note that the gap
				## characters are added to the end of the indel-part. Thus, the codons are
				## not considered.
				while (length($indel)%3 != 0) {
					$indel .= '-';
				}

				$codon_seq .= $indel;
			}
		}
	}
	print join(" ", (caller(0))[0..3]), ", leaving\n" if $debug;
	return ($codon_seq);
}#}}}

# is this ever called? apparently not.
#--------------------------------------------------
# sub protein_borders {#{{{
# 	print join(" ", (caller(0))[0..3]), "\n" if $debug;
# 	my $self = shift;
# 	my $gw = $self->{gw};
# 	for (my $i = 0; $i < @$gw; ++$i) {
# 		if ($gw->[$i] =~ /Score/) {
# 			my ($start, $end) = $gw->[$i+1] =~ /.*$self->{protname}\s+([0-9]+)\s+([0-9]+).*/;
# 			print join(" ", (caller(0))[0..3]), ", leaving\n" if $debug;
# 			return ($start, $end);
# 		}
# 		else {
# 			die "Fatal: Unable to determine protein boundaries. Check exonerate command.\n";
# 		}
# 	}
# }#}}}
#-------------------------------------------------- 

# and this? apparently not, either
#--------------------------------------------------
# sub cdna_borders {#{{{
# 	print join(" ", (caller(0))[0..3]), "\n" if $debug;
# 	my $self = shift;
# 	my $gw = $self->{gw};
# 	for (my $i = 0; $i < @$gw; ++$i) {
# 		if ($gw->[$i] =~ /Score/) {
# 			my ($start, $end) = $gw->[$i+1] =~ /.*$self->{dnaname}\s+([0-9]+)\s+([0-9]+).*/;
# 			print join(" ", (caller(0))[0..3]), ", leaving\n" if $debug;
# 			return ($start, $end);
# 		}
# 		else {
# 			die "Fatal: Unable to determine cDNA boundaries. Check exonerate command.\n";
# 		}
# 	}
# }#}}}
#-------------------------------------------------- 

# really return indels, not this f'ed up regex shit from Hamstr
# TODO fix this up, get exonerate to output indel info
sub _GetIndels {#{{{
	print join(" ", (caller(0))[0..3]), "\n" if $debug;
	my $gw = shift;
	my $indel;
	for (my $i = 0; $i < @$gw; ++$i) {
		if ($gw->[$i] =~ /Indels: /) {
			print join(" ", (caller(0))[0..3]), ", leaving\n" if $debug;
			return 1;
		}
		else {
			print join(" ", (caller(0))[0..3]), ", leaving\n" if $debug;
			return 1;
		}
	}
}#}}}

# return fuckin' true!
1;
