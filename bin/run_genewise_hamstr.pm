package run_genewise_hamstr;
use strict;
my $debug = 0;
$ENV{'WISECONFIGDIR'} =  "../wisecfg";
# this module runs genewise on a DNA sequence and a protein sequence
# and then allows to parse this result.
# the constructor creates an object containing a reference to an array
# containing the file content

# LAST Modified: 11.01.2010 renamed the file names for the genewise run to avoid overwriting of files when multipe runs are performed in parallel on the same sequence file

1;
sub new {	#mp creates a new instance of a genewise result object#{{{
	print join(" ", (caller(0))[0..3]), "\n" if $debug;
	my $self_tmp = [];
	my $self;
	my ($class, $dna, $prot, $path) = @_;	#mp $prot ^= $refseq
	if (!defined $path) {
		$path = '/tmp';
	}
	my $pid=$$;
	# the file names
	my $protname = $pid.'_protein';
	my $dnaname = $pid . '_dna';
	## print the two sequences to default path /tmp/
	open (DNA, ">$path/$dnaname") or die "could not open $path/$dnaname for writing: $!\n";
	print DNA ">$dnaname\n$dna";
	close DNA;
	open (PROTEIN, ">$path/$protname") or die "could not open $path/$protname for writing: $!\n";
	print PROTEIN ">$protname\n$prot";
	close PROTEIN;

	## run genewise on the two sequences
  `echo \$WISECONFIGDIR`;
    
	my $genewise_cmd = "genewise -trans -cdna -pep -sum $path/$protname $path/$dnaname";
	print "Running $genewise_cmd\n";
	$self_tmp = [`$genewise_cmd`];
	die "genewise failed with error code: $?\n" if ($? > 0);	#mp die if genewise crashes 
	for (my $i = 0; $i < @$self_tmp; $i++) {
		chomp $self_tmp->[$i];
		$self_tmp->[$i] =~ s/\s{1,}$//;	#mp this also removes the newlines, you dork :P
	}

	$self->{gw} = $self_tmp;
	$self->{nt_seq} = $dna;
	$self->{prot_seq} = $prot;
	$self->{protname} = $protname;
	$self->{dnaname} = $dnaname;
	$self->{gw_count} = @$self_tmp;
	$self->{get_indel} = 1; ## per default the indel-part is recovererd, rather than masked by 'N'. See code for details
	$self->{indels} = _GetIndels($self_tmp);
	bless ($self, $class);
	return $self;
}#}}}
#################
## sub score extract the score for the alignment
sub score {#{{{
	print join(" ", (caller(0))[0..3]), "\n" if $debug;
	my $self = shift;
	my $score;
	for (my $i = 0; $i < $self->{gw_count}; $i ++) {
		if ($self->{gw}->[$i] =~ /^(\d{1,}\.{0,1}\d{0,}).*/) {
			$score = $1;
			last;
		}
	}
	return ($score);
}#}}}
##################
#mp this function is never called?!
sub protein {#{{{
	print join(" ", (caller(0))[0..3]), "\n" if $debug;
	my $self = shift;
	my $gw = $self->{gw};
	my $prot = '';
	for (my $i = 0; $i < @$gw; $i++) {
		if ($gw->[$i] =~ />.*\.pep/) { #the protein seq starts
			my $count = 1;
			while ($gw->[$i+$count] ne '//') {
				my $protpart = $gw->[$i+$count];
				chomp $protpart;
				$prot .= $protpart;
				$count ++;
			}
		}
		elsif (length $prot > 0) {
			last;
		}
	}
	return($prot);
}#}}}
##################
sub translation {#{{{
	print join(" ", (caller(0))[0..3]), "\n" if $debug;
	my $self = shift;
	my $finish = 0;
	my $translated_seq = '';
	my @transtmp;

	## step 1: extract the relevant info from the genewise output
	for (my $i = 0; $i < $self->{gw_count}; $i++) {
		#if ($self->{gw}->[$i] =~ />.*.tr/) {# a translated bit starts #mp !!regex will also match ">424324_bartralle.sp" because the . is not escaped :P 
		if ($self->{gw}->[$i] =~ />.*\.tr/) {	#mp corrected regex
			print $self->{gw}->[$i], "\n" if $debug;
			while ($self->{gw}->[$i] !~ '//') {
				push @transtmp, $self->{gw}->[$i];
				$i++;
			}
			last; # end the for loop since nothing left to be done
		}
	}
	
	## step two: get the sequences
	my $count = -1;
	my $trans;
	for (my $i = 0; $i < @transtmp; $i++) {
		if ($transtmp[$i] =~ />/) {
			$count++;
			$trans->[$count]->{seq} = ''; # initialize
			if ($transtmp[$i] =~ /.*\[([0-9]{1,}):([0-9]{1,})\].*/) {
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
		$translated_seq = $trans->[0]->{seq};
	}
	else {
		for (my $i = 0; $i < @$trans; $i++) {
			$translated_seq .= $trans->[$i]->{seq};
			if ($i < (@$trans - 1)) {
				my $missing = $trans->[$i+1]->{start} - $trans->[$i]->{end} -1;	#mp never used 
				$translated_seq .= 'X';
			}
		}
	}
	return($translated_seq);
}#}}}

##################
sub codons {#{{{
	print join(" ", (caller(0))[0..3]), "\n" if $debug;
	my $self = shift;
	my $finish = 0;
	my $codon_seq = '';
	my @transtmp;

	## step 1: extract the relevant info from the genewise output
	for (my $i = 0; $i < $self->{gw_count}; $i++) {
		if ($self->{gw}->[$i] =~ />.*sp$/) {# the codons set starts
			while ($self->{gw}->[$i] !~ '//') {
				push @transtmp, $self->{gw}->[$i];
				$i++;
			}
			last; # end the for loop since nothing left to be done
		}
	}
	
	## step two: get the sequences
	my $count = -1;
	my $trans;
	for (my $i = 0; $i < @transtmp; $i++) {
		if ($transtmp[$i] =~ />/) {
			$count++;
			$trans->[$count]->{seq} = ''; # initialize
			if ($transtmp[$i] =~ /.*\[([0-9]{1,}):([0-9]{1,})\].*/) {
				$trans->[$count]->{start} = $1;
				$trans->[$count]->{end} = $2;
			}
		}
		else {
			$transtmp[$i] =~ tr/a-z/A-Z/;
			$trans->[$count]->{seq} .= $transtmp[$i];
		}
	}	#mp $count never exceeds 1 because there is never more than one .sp seq in the genewise output

	## step 3: connect the fragments
	if (@$trans == 1) {
		print '@$trans is 1, only one seq', "\n" if $debug;
		$codon_seq = $trans->[0]->{seq};
	}
	else {	#mp under what circumstances is there more than one target seq in the genewise/exonerate output?
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
				## now append gap characters until the frame is recovered. Not that the gap
				## characters are added to the end of the indel-part. Thus, the codons are
				## not considered.
				while (length($indel)%3 != 0) {
					$indel .= '-';
				}

				$codon_seq .= $indel;
			}
		}
	}
	return ($codon_seq);
}#}}}
###########################
sub protein_borders {#{{{
	print join(" ", (caller(0))[0..3]), "\n" if $debug;
  my $self = shift;
  my $gw = $self->{gw};
  for (my $i = 0; $i < @$gw; $i++) {
    if ($gw->[$i] =~ /Bits.*introns$/) {
      #my ($start, $end) = $gw->[$i+1] =~ /.*$self->{protname}\s{1,}([0-9]{1,})\s{1,}([0-9]{1,}).*/;	#mp useless use of quantifying braces in regex
      my ($start, $end) = $gw->[$i+1] =~ /.*$self->{protname}\s+([0-9]+)\s+([0-9]+).*/;	#mp made regex more clear
			print join(" ", (caller(0))[0..3]), ", leaving\n" if $debug;
      return($start, $end);
    }
    else {
      die "no protein-start and end could not be determnined. Check genewise command\n";
    }
  }
}#}}}
##########################
sub cdna_borders {#{{{
	print join(" ", (caller(0))[0..3]), "\n" if $debug;
  my $self = shift;
  my $gw = $self->{gw};
  for (my $i = 0; $i < @$gw; $i++) {
    if ($gw->[$i] =~ /Bits.*introns$/) {
      #my ($start, $end) = $gw->[$i+1] =~ /.*$self->{dnaname}\s{1,}([0-9]{1,})\s{1,}([0-9]{1,}).*/;	#mp useless use of quantifying braces in regex
      my ($start, $end) = $gw->[$i+1] =~ /.*$self->{dnaname}\s+([0-9]+)\s+([0-9]+).*/;	#mp made regex more clear
			print join(" ", (caller(0))[0..3]), ", leaving\n" if $debug;
      return($start, $end);
    }
    else {
      die "no cdna-start and end could not be determnined. Check genewise command\n";
    }
  }
}#}}}
##########################
sub _GetIndels {#{{{
	print join(" ", (caller(0))[0..3]), "\n" if $debug;
  my $gw = shift;
  my $indel;
  for (my $i = 0; $i < @$gw; $i++) {
    if ($gw->[$i] =~ /Bits/) {
      #$indel = $gw->[$i+1] =~ /.*([0-9]{1,})/;	#mp WTF is this regex? This will match the first number in the gw output :P
			#mp also, uncomfortable use of quantifiers
      $indel = $gw->[$i+1] =~ /^.*([0-9]+)\s[0-9]+$/;	#mp corrected regex
			print "Subroutine _GetIndels matched: $1\n";
      return($indel);
    }
  }
}#}}}
