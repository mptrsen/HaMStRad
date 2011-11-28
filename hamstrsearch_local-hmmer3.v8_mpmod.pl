#!/usr//bin/env perl
# PROGRAMNAME: hamstrsearch_local.pl

# Copyright (C) 2009 INGO EBERSBERGER, ingo.ebersberger@univie.ac.at
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published
# by the Free Software Foundation; either version 3 of the License
# or any later version.

# This program is distributed in the hope that it will be useful
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
# You should have received a copy of the GNU General Public License
# along with this program; If not, see http://www.gnu.org/licenses
 
use strict;
use warnings;	#mp
use lib '../bin/';
use lib '../bin/Bio';
use lib '/dmp-work/dmp/tools/hamstrsearch_local/bin';
use Getopt::Long;
use Bio::SearchIO;
use Bio::Search::Hit::BlastHit;
use run_genewise_hamstr;
use exonerate;	#mp
use Data::Dumper;	#mp

# PROGRAM DESCRIPTION: See bottom of this file.
######################## start main #############################
my $version = "hamstrsearch_local-hmmer3.v8.pl\n";
### EDIT THE FOLLOWING LINES TO CUSTOMIZE YOUR SCRIPT
## note: all the variables can also be set via the command line
my $hmmsearchprog = 'hmmsearch'; #program for the hmm search
my $wiseprog = 'genewise';

########## EDIT THE FOLLOWING TWO LINES TO CHOOSE YOUR BLAST PROGRAM ##########
#my $blast_prog = 'blastall';
my $blast_prog = 'blastp';
###############################################################################

my $alignmentprog = 'clustalw2';
my $hmmpath = '../core_orthologs'; #path where the hmms are located
my $blastpath = '../blast_dir'; #path to the blast-dbs
my $tmpdir = 'tmp';
my $wiseconfigdir = '../wisecfg';
my $eval = 1; # eval cutoff for the hmm search
my $idsep = '__'; #character used to replace whitespaces in the sequence header with (flag -longhead)

my $hmm_dir = 'hmm_dir';
my $fa_dir  = 'fa_dir';
##############################
my $pid = $$;
my $help;
my $seq2store_file='';
my $cds2store_file='';
my $hmm;
my @hmms;
my $fa;
my $fafile;
my @seqs2store;
my @cds2store;
my $ep2eg;
my $estfile;
my $aln;
my $idfile;
my $taxon_check = 0;
my $hmmset;
my $hmmsearch_dir;
my $dbfile = ''; # the file hmmsearch is run against #mp (the EST file)
my $dbfile_short;
my $taxon_file;
my $refspec_string;
my @refspec = qw();
my @primer_taxa;
my $refspec_name = '';
my $taxon_global;
my $fileobj;
my $fa_dir_neu = '';
my $cdna_dir = '';
my $cdnaobj;
my $gwrefprot;
my $seqtype;
my $align;
my $rep;
my $estflag;
my $proteinflag;
my $refseq;
my $strict;
my $relaxed;
my $refspec_final = '';
my $concat;
my $seqs2store_file;
my $append;
my $longhead;
my $check = 1;
my @log = qw();
my $filter = 'T';
my $bhh;
my $exhaustive;
my $debug;
my $use_exonerate;
my $skipcount = 0;


if (@ARGV==0) {
	$help = 1;
}
## help message
my $helpmessage = "USAGE: $0 -sequence_file=<> -hmmset=<> -taxon=<>  -refspec=<> [-est|-protein] [-hmm=<>] [-representative] [-h]\nFor more information, see perldoc $0.\n";

#mp get cmd line args
GetOptions ("h"        => \$help,
            "hmm=s"    => \$hmm,
            "est"    => \$estflag,
            "protein"=> \$proteinflag,
            "sequence_file=s" => \$dbfile,
            "fasta_file=s" => \$fafile,
            "hmmset=s" => \$hmmset,
            "hmmpath=s" => \$hmmpath,
            "taxon_file=s" => \$taxon_file,
            "taxon=s"  => \$taxon_global,
            "eval_limit=s" => \$eval,
            "refspec=s" => \$refspec_string,
            "estfile=s" => \$estfile,
            "representative" => \$rep,
      	    "strict" => \$strict,
	          "relaxed" => \$relaxed,
	          "append" => \$append,
      	    "filter=s" => \$filter,
	          "longhead" => \$longhead,
	          "rbh" => \$bhh,
            "blastpath=s" => \$blastpath,
						"d" => \$debug,
						"use_exonerate" => \$use_exonerate);


if ($help) {
  print $helpmessage;
  exit;
}
#####determine the hostname#######
#mp who the fuck cares for the hostname?
push @log, "VERSION\n\t$version\n";
my $hostname = `hostname`;
chomp $hostname;
print "$version\n";
print "hostname is $hostname\n";
push @log, "HOSTNAME\n\t$hostname\n";
#################################
$wiseprog = $use_exonerate ? 'exonerate' : 'genewise';

## 1) check if all information is available to run HaMStR
($check, @log) = &checkInput();
if ($check == 0) {
  print "$helpmessage";
  print "#######################################\n";
  print "There was an error running $version\n";
  print join "\n", @log;
  exit;
}
else {
  open (OUT, ">hamstrsearch.log") or die "could not open logfile\n";
  print OUT join "\n", @log;
  close OUT;
}
### read in of the core-ortholog sequences
#mp literally, read in the entire core-orthologs fasta file.
#mp I am so not going to mimic this.
my $co_seqs = parseSeqfile("$fafile");

## 2) loop through the hmms
## process each hmm file separately
for (my $i = 0; $i < @hmms; $i++) {
  $fileobj = undef;
  my @seqs = qw();
  my @newseqs = qw();## var to contain the sequences to be added to the orthologous cluster
  my @newcds = qw();
	my @newcdna = qw();	#mp added newcdna array: stores cdna for output
  my $hmm = $hmms[$i];
  my $hmmout = $hmm;
  $hmmout =~ s/\.hmm/\.out/;

  ## 3) run the hmm search
  if (!(-e "$hmmsearch_dir/$hmmout")) {
    print "now running $hmmsearchprog using $hmm\n";
#	print "$hmmsearchprog $hmm_dir/$hmm $dbfile >$hmmsearch_dir/$hmmout";
    !`$hmmsearchprog $hmm_dir/$hmm $dbfile >$hmmsearch_dir/$hmmout` or die "Problem running hmmsearch\n";
  }
  else {
    print "an hmmresult $hmmout already exists. Using this one!\n";
  }
  
  ## 4) process the hmm search result
  my $hitcount = 0;
  ## 4a) loop through the individual results
  ## now the modified version for hmmer3 comes
  my ($query_name, @results) = &parseHmmer3($hmmout, $hmmsearch_dir);
  if (! @results) {
    print "no hit found for $query_name\n";
    next;
  }
  ### reduce the hit list to only the best entry if $bhh flag is set
#  if ($bhh) {
#	my $keep = shift @results;
#	@results = qw();
#	$results[0] = $keep;
#}
  chomp $query_name;
  print "Results for HMM $query_name (Hits: ";	#mp inserted "HMM"
	print join(", ", @results), ")\n" if $debug;	#mp extended report
  my ($check, $refspec_final) = &determineRefspecFinal($query_name, @refspec);
  if ($check == 0) {
    die "error in retrieving refspec data\n";
  }
	HMMER_RESULT:			#mp added label for loop
  for (my $k = 0; $k < @results; $k++) {	#mp for every hmmsearch result, do the following (perhaps this is where the addtl. cdna seqs come from)
    my $hitname = $results[$k];
    print "processing hit: $hitname\n";	#mp added report
    my $keep = 0;
    my $hitseq = '';
    $refseq = '';

    ## 4b) test for the reciprocity criterion fulfilled
    ($keep, $hitseq)  = &check4reciprocity($query_name, $hitname, $refspec_final, @refspec);
    if ($keep == 1) {
			## blast search with the hmm hit identifies the core-ortholog sequence of the reference species
			## check for the taxon from the taxon_file.txt. IN FACT, THIS ROUTINE IS OUTDATED. I ALWAYS USE THE
			## GLOBAL TAXON NAME FOR THE HAMSTERED SEQUENCES.
			my $taxon = '';
      if ($taxon_check){
				if ($taxon_check == 1) {
					$taxon = &getTaxon($hitname);
				}
				elsif ($taxon_check == 2) {
					$taxon = $taxon_global;
				}
      }
      ## put the info about the hits into an object for later post-processing
      ### HERE COMES THE NEW STUFF THAT DEALS WITH THE DIFFERENT POSSIBILITIES: STRICT, RELAXED OR WHATEVER...
			#mp the fileobj is defined here
      $fileobj = &determineReferences ($fileobj, $taxon, $refspec_final, $hitname, $hitseq, $hitcount, $query_name);
      $hitcount++;
			print "Reciprocity fulfilled. hitcount: $hitcount\n" if $debug;	#mp 
		}
    else {
      print "Reciprocity not fulfilled!\n";
    }
  }	#mp end HMMER_RESULT

  ## 5) do the rest only if at least one hit was obtained
	#mp i.e. if the reciprocity criterion is fulfilled
  if (defined $fileobj) {

    ## 5a) if the hits are derived from ESTs, get the best ORF
    if ($estflag) {
		#mp the fileobj is modified here, exonerate/genewise is run
      $fileobj =  &predictORF($query_name);	#mp added arg: refseq id that belongs to this queryseq
    }

    ## 5b) if the user has chosen to postprocess the results
    if ($rep) {
      &processHits($fileobj);
    }

    ## 6) prepare the output
    my @taxa = keys(%$fileobj);
		TAXON:	#mp added label for loop
    for (my $i = 0; $i< @taxa; $i++) {
      if ($rep) {
				push @newseqs, ">$query_name|$fileobj->{$taxa[$i]}->{refspec_final}|$taxa[$i]|$fileobj->{$taxa[$i]}->{refid}";
				push @newseqs, $fileobj->{$taxa[$i]}->{refprot};
				if ($estflag) {
					push @newcds, ">$query_name|$fileobj->{$taxa[$i]}->{refspec_final}|$taxa[$i]|$fileobj->{$taxa[$i]}->{refid}";
					push @newcds, $fileobj->{$taxa[$i]}->{refcds};
					#mp added cdna output
					push @newcdna, ">$query_name|$fileobj->{$taxa[$i]}->{refspec_final}|$taxa[$i]|$fileobj->{$taxa[$i]}->{refid}";
					push @newcdna, $fileobj->{$taxa[$i]}->{refcdna};
					#mp end added cdna output
				}
      }
      else {
				my $idobj = $fileobj->{$taxa[$i]}->{ids};	#mp ids are the seq IDs in the EST file
				my $protobj = $fileobj->{$taxa[$i]}->{prot};
				my $cdsobj  = $fileobj->{$taxa[$i]}->{cds};
				#--------------------------------------------------
				# my $cdnaobj = $fileobj->{$taxa[$i]}->{cdna};
				#-------------------------------------------------- 
				my $refspecobj = $fileobj->{$taxa[$i]}->{refspec};
				for (my $j = 0; $j < @$idobj; $j++) {
					push @newseqs, ">$query_name|$refspecobj->[$j]|$taxa[$i]|$idobj->[$j]";
					push @newseqs, $protobj->[$j];
					if ($estflag) {
						push @newcds, ">$query_name|$taxa[$i]|$idobj->[$j]|$refspecobj->[$j]";
						push @newcds, $cdsobj->[$j];	#mp don't worry everything is fine
						#mp start compiling cdna data... I hope it's structured like everything else
						#--------------------------------------------------
						# push @cdna, ">$query_name|$taxa[$i]|$idobj->[$j]|$refspecobj->[$j]_cdna";	#mp compile cdna for output
						# push @cdna, $cdnaobj->[$j];	#mp compile cdna for output
						#-------------------------------------------------- 
					}
				}
      }
      my $refs = $co_seqs->{$query_name};
      for (keys %$refs) {	#mp apparently only print cds if -est option was selected
				my $line = ">$query_name|$_|" . $refs->{$_}->{seqid} . "\n" . $refs->{$_}->{seq};
				push @seqs, $line;
      }
      chomp @seqs;
      print "\n";
      @seqs = (@seqs, @newseqs);	#mp y u no use push() like a sane person

			#mp output
      open (OUT, ">$fa_dir_neu/$query_name.fa");
      print OUT join "\n", @seqs;
      print OUT "\n";
      close OUT;

      if ($estflag) {	#mp apparently only print cds if -est option was selected

				#mp more output
				open (OUT, ">$fa_dir_neu/$query_name.cds.fa");
				print OUT join "\n", @newcds;
				close OUT;

				#mp output cdna to cdna file
				open (OUT, ">$cdna_dir/$query_name.cdna.fa") or die "Fatal: Could not open $cdna_dir/$query_name\_cdna.fa: $!\n";
				print OUT join "\n", @newcdna;
				close OUT;
				#mp end cdna output
      }
      for (my $i = 0; $i < @newseqs; $i+= 2) {
				my $line = $newseqs[$i] . "|" . $newseqs[$i+1];
				$line =~ s/>//;
				push @seqs2store, $line;
				if ($estflag) {
					my $cdsline = $newcds[$i] . "|" . $newcds[$i+1];
					$cdsline =~ s/>//;
					push @cds2store, $cdsline;
				}
      }
    }	#mp end TAXON
  }
	print "\n";
}	#mp end of loop through the hmms

if (@seqs2store > 0) {
#mp added system error messages to die messages
  if ($append) {
    open (OUT, ">>$seqs2store_file") or die "failed to open output file $seqs2store_file: $!\n";
  }
  else {
    open (OUT, ">$seqs2store_file") or die "failed to open output file $seqs2store_file: $!\n";
  }
  print OUT join "\n", @seqs2store;
  print OUT "\n";
  close OUT;
  if ($estflag) {
    if ($append) {
      open (OUT, ">>$cds2store_file") or die "failed to open output file $cds2store_file: $!\n";
    }
    else {
    open (OUT, ">$cds2store_file") or die "failed to open output file $cds2store_file: $!\n";
    }
    print OUT join "\n", @cds2store;
    print OUT "\n";
    close OUT;
  }
}
else {
  print "no hits found\n";
}
print "Done! $skipcount couples skipped during exonerate post-processing.\n";
exit;

##################### start subs ###############

####### checkInput performs a number of checks whether sufficient information
### and all data are available to run HaMStR
sub checkInput {#{{{
	print join(" ", (caller(0))[0..3]), "\n" if $debug;
  my @log;
  my $check = 1;
  $dbfile_short = $dbfile;
  $dbfile_short =~ s/\..*//;
  ## 0) Check for presence of the file with the sequences that should be hamstered
  if (-e "$dbfile") {
		#the file exists
		print "Using EST file $dbfile\n";
		push @log, "INFILE PROCESSING\n";
		print "removing newlines from the infile $dbfile such that a sequence forms a consecutive string\n";
		`../bin/nentferner.pl -in=$dbfile -out=$dbfile.mod`;
		if (-e "$dbfile.mod") {
			print "Newlines from the infile have been removedpl succeeded.\n";
			push @log, "\tNewlines from the infile have been removed\n";
			$dbfile = $dbfile . '.mod';
			
			if (defined $longhead) {
				`sed -i -e "s/[[:space:]]\\+/$idsep/g" -e 's/\\(>.\\{20\\}\\).*/\\1/' $dbfile`;

				push @log, "\tOption -longhead was chosen. Replaced whitespaces in the sequence identifier with '$idsep'\n";
			}
		}
		else {
			push @log, "Problems running the script nentferner.pl ... who cares\n";	#mp who cares...
			$check = 1;	#mp made the check always successful, fuck nentferner
		}
  }
  else {
	#the file does not exist:
	push @log, "The specified infile $dbfile does not exist. PLEASE PROVIDE A VALID INFILE!\n";
	$check = 0;
	print join(" ", (caller(0))[0..3]), ", leaving\n" if $debug;
	return ($check, @log);
  }
  ## 1) check for filetype
  print "Checking for filetype:\t";
  if (!defined $estflag and !defined $proteinflag) {
    push @log, "please determine the sequence type. Choose between the options -EST or -protein";
    print "failed\n";
    $check = 0;
  }
  else {
		if ($estflag) {
			$estfile = $dbfile;
			$dbfile = "$dbfile.tc";
	 
			push @log, "HaMStR will run on the ESTs in $estfile";
			push @log, "Translating ESTs";
			if (!(-e "$dbfile")) {
				print "translating $estfile, this may take a while\n";
				`../bin/translate.pl -in=$estfile -out=$dbfile`;
				open (LOG, "hamstrsearch.log");
				my @info = <LOG>;
				@log = (@log, @info);
				close LOG;
			}
			else {
				push @log, "Translated file already exists, using this one\n";
			}
			if (! -e "$dbfile") {
				push @log, "The translation of $estfile failed. Check the script translate.pl";
				print "failed\n";
				$check = 0;
			}
		}
		else {
				## file type is protein
				print "succeeded\n";
		}
  }
  ## 2) Check for presence of blastall
  print "Checking for the blast program\t";
  if (`which $blast_prog` =~ / no /) {
    push @log, "could not execute $blast_prog. Please check if this program is installed and executable";
    print "failed\n";
    $check = 0;
  }
  else {
    push @log, "check for $blast_prog succeeded";
    print "succeeded\n";
  }
  ## 3) Check for presence of hmmsearch
  print "Checking for hmmsearch\t";
  if (! `$hmmsearchprog -h`) {
    push @log, "could not execute $hmmsearchprog. Please check if this program is installed and executable";
    print "failed\n";
    $check = 0;
  }
  else {
      push @log, "check for $hmmsearchprog succeeded\n";
      print "succeeded\n";
  }

  ## 4) Check for presence of the directory structure
  print "checking for presence of the hmm files ";
  if (!(defined $hmmset)) {
    $hmmpath = '.';
    $hmmset = 'manual';
		print "$hmmpath/$hmmset";
  }
  else {
    $hmmpath = "$hmmpath/$hmmset";
    $fafile = "$hmmpath/$hmmset" . '.fa';
		print "$hmmpath/$hmmset";
  }
  $hmm_dir = "$hmmpath/$hmm_dir";
  $hmmsearch_dir = 'hmm_search_' . $dbfile_short . '_' . $hmmset;
  ## 5) check for the presence of the hmm-files and the fasta-file
  if (!(-e "$hmm_dir")) {
    push @log, "Could not find $hmm_dir";
    print "\tfailed\n";
    $check = 0;
  }
  else {
    if (defined $hmm) {
	@hmms = split ',', $hmm;
	chomp @hmms;
	### check for the presence of all hmms
	for (my $k = 0; $k < @hmms; $k++) {
      		if (! -e "$hmm_dir/$hmms[$k]") {
			push @log, "$hmms[$k] has been defined but could not be found in $hmm_dir/$hmms[$k]";
			$check = 0;
			last;
      		}
      		else {
			push @log, "$hmms[$k] has been found\n";
     	 	}	
    	}
	}
	
    else {
      push @log, "running HaMStR with all hmms in $hmm_dir";
      @hmms = `ls $hmm_dir`;
    }
    chomp @hmms;
    print "\tsucceeded\n";
  }
  
  ## 6) Test for presence of the fasta file containing the sequences of the core-ortholog cluster
  print "checking for presence of the core-ortholog file $fafile".":\t";
  if (defined $fafile) {
    if (! -e "$fafile") {
      push @log, "Could not find the file $fafile";
      print "failed\n";
      $check = 0;
    }
    else {
      push @log, "check for $fafile succeeded";
      print "succeeded\n";
    }
  }
  else {
    push @log, "Please provide path and name of fasta file containing the core-ortholog sequences";
    $check = 0;
    print "failed\n";
  }
  ## 7) Checks for the taxon_file
  print "testing whether the taxon has been determined:\t";  
  if (!(defined $taxon_file) or (!(-e "$taxon_file"))) {
    if (defined $taxon_global) {
      push @log, "using default taxon $taxon_global for all sequences";
      print "succeeded\n";
      $taxon_check = 2;
    }
    else {
      push @log, "No taxon_file found. Please provide a global taxon name using the option -taxon";
      print "failed\n";
      $check = 0;
    }
  }
  else {
    push @log, "using the file $taxon_file as taxon_file";
    print "succeeded\n";
    $taxon_check = 1;
  }
  ## 8) Check for reference taxon
  print "Checking for reference species and blast-dbs\t";
  if (!(defined $refspec_string) and (! defined $strict and ! defined $relaxed)) {
      push @log, "Please provide a reference species for the reblast!";
      print "failed\n";
      $check = 0;
  }
  elsif (defined $strict or defined $relaxed) {
      if (! defined $refspec_string) {
	  ## The user has not provided a string of reference taxa. Chose all from the fasta file containing
	  ## the core orthologs.
	  @refspec = `grep '>'  $fafile |cut -d '|' -f 2 |sort |uniq`;
	  chomp @refspec;
	  $refspec_string = join ',', @refspec;
      }
      else {
	  @refspec = split /,/, $refspec_string;
      }
      if ($strict) {push @log, "Strict flag has been set. Reference species for the reblast: $refspec_string";}
      else {push @log, "Relaxed flag has been set. Reference species for the reblast: $refspec_string";}
      if (@refspec == 0) {
	  print "failed\n";
	  $check = 0;
      }
      else {
	  print "succeeded\n";
      }
  }
  else {
    push @log, "Reference species for the re-blast: $refspec_string";
    @refspec = split /,/, $refspec_string;
    $refspec_name = $refspec[0];
    print "succeeded\n";
  }
  ## 9) Check for presence of the required blast dbs
  print "checking for blast-dbs:\t";
  for (my $i = 0; $i < @refspec; $i++) {
    my $blastpathtmp = "$blastpath/$refspec[$i]/$refspec[$i]" . "_prot";
    if (! (-e "$blastpathtmp.pin")) {
      push @log, "please edit the blastpath. Could not find $blastpathtmp";
      print "$blastpathtmp failed\n";
      $check = 0;
    }
    else {
      push @log, "check for $blastpathtmp succeeded";
      print "succeeded\n";
    }
  }

  ## 10) Set the file where the matched seqs are found
  my $strictstring = '';
  if (defined $strict) {
        $strictstring = '.strict';
}
  $seqs2store_file = 'hamstrsearch_' . $dbfile_short . '_' . $hmmset . $strictstring . '.out';
  $cds2store_file = 'hamstrsearch_' . $dbfile_short . '_' . $hmmset . '_cds' . $strictstring . '.out';

  ## 11) check for filter setting for BLAST
  print "checking for low complexity filter setting:\t";
	$filter =~ tr/ft/FT/;
	if ($filter ne 'T' and $filter ne 'F') {
		push @log, "Filter is set to $filter. Please set the low complexity filter either to F or T.";
		print "low complexity filter check failed\n";
		$check = 0;
   	}
	else {
	push @log, "check for low complexity filter setting succeeded. Chosen value is $filter";
	print "succeeded\n";
}

  ## 12) apply the evalue-cut-off to the hmmsearch program
  $hmmsearchprog = $hmmsearchprog . " -E $eval";
  push @log, "hmmsearch: $hmmsearchprog";

  ## 13) setting up the directories where the output files will be put into.
  $fa_dir_neu = 'fa_dir_' . $dbfile_short . '_' . $hmmset . '_' . $refspec[0];
	$cdna_dir = 'cdna_dir_' . $dbfile_short . '_' . $hmmset . '_' . $refspec[0];	#mp add cdna output dir
  if ($strict) {
      $fa_dir_neu .= '_strict';
  }
  if ($relaxed) {
      $fa_dir_neu .= '_relaxed';
  }
  if ($check == 1) {
    if (!(-e "$hmmsearch_dir")) {
      `mkdir -p $hmmsearch_dir`;	#mp added -p flag to mkdir
    }
    if (!(-e "$fa_dir_neu")) {
      `mkdir -p $fa_dir_neu`;	#mp added -p flag to mkdir
    }
		#mp added creation of cdna output dir
		if (!(-e "$cdna_dir")) {
			`mkdir -p $cdna_dir`;	#mp with -p flag
		}
		#mp end add creation of cdna output dir
    if (!(-e "$tmpdir")) {
      `mkdir -p $tmpdir`;	#mp added -p flag to mkdir
    }
    if (!(-e "hamstrsearch_$dbfile_short")) {
			`mkdir -p hamstrsearch_$dbfile_short`; #mp added mkdir -p for this dir
		}
  }
  ## 14) determin whether or not the -representative flag has been set
  if (defined $rep) {
	push @log, "HaMStR will run with the -representative option";
	}
  else {
	push @log, "HaMStR was called without the -representative option. More than one ortholog may be identified per core-ortholog group!";
	} 
	if (defined $use_exonerate) {
		push @log, "Using exonerate instead of genewise";
		print "Using exonerate\n";
	}
	print "Finished checking input.\n" if $debug;
	print join(" ", (caller(0))[0..3]), ", leaving\n" if $debug;
  return ($check, @log);
}#}}}

#################
## check4reciprocity is the second major part of the program. It checks
## whether the protein sequence that has been identified by the hmmsearch
## identifies in turn the protein from the reference taxon that was used to
## build the hmm.
sub check4reciprocity {#{{{
	print join(" ", (caller(0))[0..3]), "\n" if $debug;
  my ($query_name, $hitname, $refspec_final, @refspec) = @_;
  my $searchdb;
  my $strict_suc = -1; # keeps track of success for all taxa
  my $relaxed_suc = 0; # keeps track of success for at least one taxon
  ## get the sequence that was identified as hit in the pHMM search from the db_file
  my $hitseq = `grep -m 1 -A 1 ">$hitname\$" $dbfile | tail -n 1`;
  if (!defined $hitseq) {
    print "could not retrieve a sequence for $hitname. Skipping...\n";
		print join(" ", (caller(0))[0..3]), ", leaving\n" if $debug;
    return(0, '', '', '');
  }
  
  ## continue with the blast
  chomp $hitseq;
  ## now run the blast
  open (OUT, ">$tmpdir/$$.fa") or die "could not open out for writing\n";
  print OUT ">$hitname\n$hitseq";
  close OUT;
  
  ## now comes the new part that does one to many blast searches. We need to iterate through all
  ## entries in the file $refspec_final and perform the Blast against each reftaxon. Note, unless
  ## $strict or $relaxed flags are set, there will be only a single reftaxon. If $relaxed is chosen
  ## then we can stop the blast searches as soon as the reciprocity is fulfilled.
  for (my $k = 0; $k < @$refspec_final; $k++) {
    my $orthocount = $refspec_final->[$k]->{orthocount};
    ## 1) Perform the blast search with the k-th reftaxon
    print "Reftaxon: $refspec_final->[$k]->{refspec}\n";
    if ($blast_prog =~ /blastp/) {
      !`$blast_prog -db $refspec_final->[$k]->{searchdb} -num_descriptions 10 -num_alignments 10 -query $tmpdir/$$.fa  -out $tmpdir/$$.blast` or die "Problem running blastp\n";
    ### postprocess the outfile
	`sed -i -e 's/Length=\\([0-9]*\\)/  (\\1 letters)/' -e 's/^\\(>*\\)lcl|/\\1/' $tmpdir/$$.blast`;	
    }
    else {
      !`blastall -p blastp -d $refspec_final->[$k]->{searchdb} -F $filter -i $tmpdir/$$.fa -o $tmpdir/$$.blast` or die "Problem running blast\n";
    }
    ## 2) now parse the best blast hit

    my @hits = &getBestBlasthit("$tmpdir/$$.blast");
    if (@hits > 0) {
      my $idsref = $refspec_final->[$k]->{refid};
      my @original_ids = @$idsref;
      print "core_orthologs: ", join "\t", @original_ids , "\n";
      ## now loop through the best hits with the same evalue and check whether
      ## among these I find the same seq as in $original
      my $i = 0;
      my $suc = 0; # keeps track of success for a single taxon
      while ($suc == 0 and $i <@hits) {
				print "blast-hit: $hits[$i]";
				## now loop through all the refspec-sequences in the hmm file
				my $j = 0;
				while ($suc == 0 and $j < @original_ids) {
					if ($original_ids[$j] eq $hits[$i]) {
						print "\thitting\n";
						$refspec_final->[$k]->{hit} = $j;
						$suc = 1;
						$relaxed_suc = 1;
					}
					else {
						print "\nnot hitting $original_ids[$j]\n";
						$j ++;
					}
					if ($suc == 1) {
						$relaxed_suc = 1;
						if ($strict_suc == -1) {
							$strict_suc = 1;
						}
					}
				}
			$i++;
      }
      if ($suc == 0) {
				$strict_suc = 0; # none of the blast hits matched against the reftaxon seq
      }
    }
    else {
      print "no hit obtained\n";
      $strict_suc = 0;
    }	
    ## when the user has chosen the strict flag, there is no reason to continue when $suc
    ## has remained 0 (reciprocity criterion not fulfilled). Thus, return to main.
    if ($strict and $strict_suc == 0) {
			print join(" ", (caller(0))[0..3]), ", leaving\n" if $debug;
      return (0, $hitseq);
    }
  }

  if ($relaxed_suc == 1) {
		print join(" ", (caller(0))[0..3]), ", leaving\n" if $debug;
    return (1, $hitseq);
  }
  else {
		print join(" ", (caller(0))[0..3]), ", leaving\n" if $debug;
    return (0, $hitseq);
  }
}#}}}
#############

sub getBestBlasthit {#{{{
	print join(" ", (caller(0))[0..3]), "\n" if $debug;
    my @hits;
    my ($file) = @_;
    my $searchio = Bio::SearchIO->new(-file        => $file,
				      -format      => 'blast',
				      -report_type => 'blastp') or die "parse failed";
    while(my $result = $searchio->next_result){
	my $count = 0;
	my $sig;
	my $sig_old;
	while( my $hit = $result->next_hit){
	    ## now I enter all top hits having the same evalue into the result
	    $sig = $hit->score;
	    if (!defined $sig_old) {
		$sig_old = $sig;
	    }
	    if ($sig == $sig_old) {
		push @hits, $hit->accession;
	    }
	    else {
		last;
	    }
	}
    }
		print join(" ", (caller(0))[0..3]), ", leaving\n" if $debug;
    return(@hits);
}#}}}
##################
sub getTaxon {#{{{
	print join(" ", (caller(0))[0..3]), "\n" if $debug;
    my ($hitname) = @_;
#    my $q = "select name from taxon t, est_project e, est_info i, annotation_neu a where a.id = $hitname and a.contig_id = i.contig_id and i.project_id = e.project_id and e.taxon_id = t.taxon_id";
    if ($hitname =~ /\D/) {
	$hitname =~ s/_.*//;
    }
    my $taxon = `grep -m 1 "^$hitname," $taxon_file | sed -e 's/^.*,//'`;
    chomp $taxon;
    $taxon =~ s/^[0-9]+,//;
    $taxon =~ s/\s*$//;
    $taxon =~ s/\s/_/g;
    if ($taxon) {
	print join(" ", (caller(0))[0..3]), ", leaving\n" if $debug;
	return ($taxon);
    }
    else {
	print join(" ", (caller(0))[0..3]), ", leaving\n" if $debug;
	return();
    }
}#}}}
###############
sub determineReferences {#{{{
	print join(" ", (caller(0))[0..3]), "\n" if $debug;
	my ($fileobj, $taxon, $refspec_final, $hitname, $hitseq, $hitcount, $query_name) = @_;
	my $refseq = '';
	my $refspec;
	## now we have to distinguish between three cases:
	## 1) hamstr is running in normal mode and one refspec has been determined. In this case, $refspec_final
	## contains data only from a single species.
	## 2) hamstr is running in normal mode and alternative refspecs have been determined by the user.
	## $refspec_final may contain results from more than one species, but we need to consider only the first
	## entry.
	## 3) hamstr is running in the strict mode. In this case $refspec_final contains data from several taxa and we need
	## to select the taxon and sequence that is most similar to the hamstered sequence.
	## 4) hamstr is running in the relaxed mode. In this case $refspec_final may contain data from several taxa and
	## we need to select the taxon and the sequence that is most similar to the hamstered sequence.
	if (defined $strict or defined $relaxed) {
	## more than one refspec. Now find the one that fits best
		my $max_score = 0;
		for (my $i = 0; $i < @$refspec_final; $i++) {
				## first, check whether the reciprocity criterion has been fulfilled
			if (defined $refspec_final->[$i]->{hit}) {
				my $rcn = $refspec_final->[$i]->{hit};
				my $refseq_cand = $refspec_final->[$i]->{sequence}->[$rcn];
				my $refspec_cand_id = $refspec_final->[$i]->{refid}->[$rcn];
				my $refspec_cand = $refspec_final->[$i]->{refspec};
				my $score = &getAlignmentScore($refseq_cand, $hitseq);
				if ($score > $max_score) {
					$refspec = $refspec_cand;
					$refseq = $refseq_cand;
					$max_score = $score;
				}
			}
		}
	}
	else { ## no choice, just one refspec
		my $rcn = $refspec_final->[0]->{hit};
		$refseq = $refspec_final->[0]->{sequence}->[$rcn];
		$refspec = $refspec_final->[0]->{refspec};
	}
	$fileobj->{$taxon}->{prot}->[$hitcount] = $hitseq;
	$fileobj->{$taxon}->{ids}->[$hitcount] = $hitname;
	$fileobj->{$taxon}->{refseq}->[$hitcount]= $refseq;
	$fileobj->{$taxon}->{refspec}->[$hitcount] = $refspec;
	#--------------------------------------------------
	# #mp added WARNING condition in case some bug prevents refspec from being set (most likely leads to genewise/exonerate crashing later on)
	#-------------------------------------------------- 
	if (!$refspec) {
		my $errorfile = $dbfile.'_errors.txt';
		warn "refspec for query ID $hitname is empty!\nThis is a bug and will probably lead to genewise/exonerate crashing later. Please report it to the developers.\n";
		open(my $errfh, ">>$errorfile") || die "Could not open $errorfile: $!\n";
		print $errfh scalar(localtime), " Missing refspec for $hitname from hmmsearching with $query_name\n";
		close $errfh || die "Could not close $errorfile: $!\n";
	}
	else {
		print "refspec for query ID \"$hitname\" is \"$refspec\"\n";
	}
	#--------------------------------------------------
	# #mp end add WARNING condition
	#-------------------------------------------------- 
	print join(" ", (caller(0))[0..3]), ", leaving\n" if $debug;
	return($fileobj);
}#}}}
###############
sub processHits {#{{{
	print join(" ", (caller(0))[0..3]), "\n" if $debug;
  my ($fileobj) = @_; 
  ## 1) align all hit sequences for a taxon against the reference species
  my @taxa = keys(%$fileobj);
#mp let's just test what happens if we only align ONE hit against the refspec
  for (my $i = 0; $i < @taxa; $i++) {
		if ($fileobj->{$taxa[$i]}->{prot}) {
			print "$taxa[$i] has prot defined, doing orfRanking...\n" if $debug;
			&orfRanking($taxa[$i]);
			last;
		}
  }
}  #}}}
  

################
sub predictORF {#{{{
	print join(" ", (caller(0))[0..3]), "\n" if $debug;
	my $query_name = shift;	#mp HMM name
  my $fileobj_new;
  my @taxa = keys(%$fileobj);	#mp the taxa that you are hamstring -> normally only one

	TAXON:	#mp labelled the loop
  for (my $i = 0; $i < @taxa; $i++) {
    my $protobj 		= $fileobj->{$taxa[$i]}->{prot};		#mp indented
    my $idobj 			= $fileobj->{$taxa[$i]}->{ids};			#mp indented
    my $refseqobj 	= $fileobj->{$taxa[$i]}->{refseq};	#mp indented
    my $refspecobj 	= $fileobj->{$taxa[$i]}->{refspec};	#mp indented
    my @ids = @$idobj;

		ID:	#mp labelled the loop. These are the individual sequences from the EST file that produced a hmmsearch hit. For each, do:
    for (my $j = 0; $j < @ids; $j++) {
			my $refseq = $refseqobj->[$j];
			my $refspec = $refspecobj->[$j];
			## determine the reading frame and remove the info from the ID
			my ($rf) = $ids[$j] =~ /.*_RF([0-9]+)/;
			print "rf is $rf\n";
			$ids[$j] =~ s/_RF.*//;
			#mp removed "\\b" from regex and added info line
			#mp this regex would almost never work because headers are not trimmed anymore
			print "running: \"grep -A 1 \">$ids[$j]\" $estfile |tail -n 1\"\n" if $debug;
			my $est = `grep -A 1 ">$ids[$j]" $estfile |tail -n 1`;	#mp get the sequence from the original EST file
			if (! $est) {
				die "error in retrieval of est sequence for $ids[$j] in subroutine processHits\n";
			}
			print "found $ids[$j] in $estfile\n" if $debug;
			## the EST is translated in rev complement
			if ($rf > 3) {
				$est = revComp($est);
			}
			#--------------------------------------------------
			# printf "running $wiseprog with:\n>\"%s|%s\"\n\"%s\"\nand:\n>\"%s|%s\"\n\"%s\"\n", $taxa[$i], $ids[$j], $est, $refseq_id, $refspec, $refseq;	#mp added info line
			#-------------------------------------------------- 
			# TODO why are $refspec and $refseq sometimes left empty?
			#mp run either exonerate or genewise, depending on invocation
			my $gw = ($use_exonerate) ? exonerate->new($est, $refseq, "$tmpdir") : run_genewise_hamstr->new($est, $refseq, "$tmpdir");

			#mp save genewise/exonerate output to file
			my $gwoutfile = $gw->{protname} . "-" . $gw->{dnaname} . '_' . $refspec . '_' . $query_name . "_" . "_$ids[$j]-" . '.' .$wiseprog . "out";
			my $gwoutput = $gw->{gw};
			open(my $gwresultfh, ">$tmpdir/$gwoutfile") or die "Could not open for writing: $!\n";
			print $gwresultfh join("\n", @$gwoutput);
			close $gwresultfh or die "Could not close file $gwoutfile: $!\n";
			print "Wrote exonerate output to $gwoutfile\n" if $debug;
			#mp end save genewise/exonerate output 

			#mp Skip a couple if exonerate doesn't return anything 
			#mp most likely if it has been handed an empty prot sequence (for whatever reason)
			if ($gw->{gw_count} == 0) {
				++$skipcount;
				print "$ids[$j] and $query_name|$refspec returned an empty exonerate result, skipping this couple ($skipcount skipped so far).\n";
				next TAXON;
			}
			#mp end skip couple
			my $translation = $gw->translation;
			my $cds = $gw->codons;
			$translation =~ s/[-!]//g;	#mp deletes gaps and stop codons
			$fileobj_new->{$taxa[$i]}->{ids}->[$j] = $ids[$j];
			$fileobj_new->{$taxa[$i]}->{prot}->[$j] = $translation;
			$fileobj_new->{$taxa[$i]}->{cds}->[$j] = $cds;
			$fileobj_new->{$taxa[$i]}->{cdna}->[$j] = $gw->extract_cdna;
			$fileobj_new->{$taxa[$i]}->{refseq}->[$j] = $refseq;
			$fileobj_new->{$taxa[$i]}->{refspec}->[$j] = $refspec;
		}
  }
	print join(" ", (caller(0))[0..3]), ", leaving\n" if $debug;	#mp added debug info
  return($fileobj_new);
}#}}}
############################
sub orfRanking {#{{{
	print join(" ", (caller(0))[0..3]), "\n" if $debug;
  my ($spec) = @_;
  my $result;
  my $refprot;
  my $refcds;
	my $refcdna;	#mp added $refcdna
  my @toalign;
  my $protobj = $fileobj->{$spec}->{prot};
  my $idobj = $fileobj->{$spec}->{ids};

  my $refcluster; ## variables to take the cluster and its id for later analysis
  my $refid;
  if (scalar @$protobj == 1) {	#mp added 'scalar' for clarity
    ## nothing to chose from
    $refprot = $protobj->[0];
    $refcds = $fileobj->{$spec}->{cds}->[0];
		$refcdna = $fileobj->{$spec}->{cdna}->[0];	#mp added refcdna
    my $length = length($refprot);
    $refid = $idobj->[0] . "-" . $length;
  }
  else {
    ## more than one cluster
		## note, I set the refseq fix to the first entry. This is to avoid that in this routine 
		## sequences from different taxa are used.  
		push @toalign, ">$fileobj->{$spec}->{refspec}->[0]";
		push @toalign, $fileobj->{$spec}->{refseq}->[0];
		## now walk through all the contigs
		for (my $i = 0; $i < @$protobj; $i++) {
			my @testseq = (">$idobj->[$i]", $protobj->[$i]);
			@testseq = (@testseq, @toalign);	#mp y u no use push() for readability

			open (OUT, ">$tmpdir/$pid.ref.fa") or die "could not open file for writing refseqs\n";
			print OUT join "\n", @testseq;
			close OUT;

			## run clustalw	#mp
			print "running '$alignmentprog $tmpdir/$pid.ref.fa -output=fasta -outfile=$tmpdir/$pid.ref.aln 2>&1 >$tmpdir/$pid.ref.log'\n" if $debug;	#mp added debug msg
			!(`$alignmentprog $tmpdir/$pid.ref.fa -output=fasta -outfile=$tmpdir/$pid.ref.aln 2>&1 >$tmpdir/$pid.ref.log`) or die "error running $alignmentprog\: $!\n";	#mp added error message
			## get the alignment score
			$result->[$i]->{score} =  `grep "Alignment Score" $tmpdir/$pid.ref.log |sed -e 's/[^0-9]//g'`;	#mp don't use external grep/sed :P
			if (!$result->[$i]->{score}) {
	      die "error in determining alignment score\n";
			}
			chomp $result->[$i]->{score};

			## get the aligned sequence
			open (ALN, "$tmpdir/$pid.ref.aln") 
				or die "failed to open alignment file\n";
			my @aln = <ALN>;
			close ALN;

			## remove the terminal gaps
			my $aseq = extractSeq($idobj->[$i], @aln);
			$aseq =~ s/-*$//;
			$result->[$i]->{aend} = length $aseq;
			my ($head) = $aseq =~ /^(-*).*/;
			($result->[$i]->{astart}) = length($head)+1;
		}
		### the results for all seqs has been gathered, now order them
		$result = &sortRef($result);

		print 'result: ', Dumper($result) if $debug;	#mp added debug output
		($refprot, $refcds, $refcdna, $refid) = &determineRef($result,$spec);	#mp added $refcdna -- THIS IS where the seq fragments got concatenated
  }
	print 'before: ', Dumper($fileobj) if $debug;	#mp added debug output
  $fileobj->{$spec}->{refprot} = $refprot;
  $fileobj->{$spec}->{refcds}  = $refcds;
  $fileobj->{$spec}->{refcdna}  = $refcdna;	#mp added refcdna
  $fileobj->{$spec}->{refid}   = $refid;
  $fileobj->{$spec}->{refspec_final} = $fileobj->{$spec}->{refspec}->[0];
	print 'after: ', Dumper($fileobj) if $debug;	#mp added debug output
	print join(" ", (caller(0))[0..3]), ", leaving\n" if $debug;
  return();
}#}}}
###########################
sub sortRef {#{{{
	print join(" ", (caller(0))[0..3]), "\n" if $debug;
	my $result = shift;
	my @sort;
	for (my $i = 0; $i < @$result; $i++) {
		push @sort, "$i,$result->[$i]->{astart},$result->[$i]->{aend},$result->[$i]->{score}";
	}
	open (OUT, ">$tmpdir/$pid.sort") or die "failed to write for sorting\n";
	print OUT join "\n", @sort;
	close OUT;
	`sort -n -t ',' -k 2 $tmpdir/$pid.sort >$tmpdir/$pid.sort.out`;	#mp how the f*ck can you use the external sort instead of the builtin sort function?!
	@sort = `less $tmpdir/$pid.sort`;	#mp WTF happened to reading a file with open()?
	chomp @sort;
	$result = undef;
	for (my $i = 0; $i < @sort; $i++) {
		($result->[$i]->{id}, $result->[$i]->{start}, $result->[$i]->{end}, $result->[$i]->{score}) = split ',', $sort[$i];
	}
	print join(" ", (caller(0))[0..3]), ", leaving\n" if $debug;
	return($result);
}#}}}
########################
sub determineRef {#{{{
	print join(" ", (caller(0))[0..3]), "\n" if $debug;
  my ($result, $spec) = @_;
  my $lastend = 0;
  my $lastscore = 0;
  my $final;
  my $count = 0;
  my $id = '';
  for (my $i = 0; $i < @$result; $i++) {
    if ($result->[$i]->{start} < $lastend or $lastend == 0) {
      if ($result->[$i]->{score} > $lastscore) {
				$lastend = $result->[$i]->{end};
				$lastscore = $result->[$i]->{score};
				$id = $result->[$i]->{id};
      }
    }
    elsif ($result->[$i]->{start} > $lastend) {
      ## a new part of the alignment is covered. Fix the results obtained so far
			#mp this seems to be the part where it is decided whether to concatenate or not -- INDEED :P
      $final->[$count]->{id} = $id;
      $lastend = $result->[$i]->{end};
      $id = $result->[$i]->{id};
      $count++;
    }
  }
  $final->[$count]->{id} = $id;
  ## now concatenate the results
  my $refprot = '';
  my $refid = '';
  my $refcds = '';
	my $refcdna = '';	#mp added refcdna, maybe this helps
  for (my $i = 0; $i < @$final; $i++) {
    my $seq = $fileobj->{$spec}->{prot}->[$final->[$i]->{id}];
    my $cdsseq = $fileobj->{$spec}->{cds}->[$final->[$i]->{id}];
		my $cdnaseq = $fileobj->{$spec}->{cdna}->[$final->[$i]->{id}];	#mp added cdnaseq
    my $length = length($seq);
    $refid .= "$fileobj->{$spec}->{ids}->[$final->[$i]->{id}]-$length" . "PP";
    $refprot .= $seq;
    if ($estflag) {
      $refcds .= $cdsseq;
			$refcdna .= $cdnaseq;	#mp added refcdna
    }
  }
	#mp remove trailing 'PP'
  $refid =~ s/PP$//;
	print join(" ", (caller(0))[0..3]), ", leaving\n" if $debug;
  return($refprot, $refcds, $refcdna, $refid);	#mp added returning $refcdna
}#}}}
#############################
sub extractSeq {#{{{
	print join(" ", (caller(0))[0..3]), "\n" if $debug;
  my ($id, @aln) = @_;
  my $seq = '';
  my $start = 0;
  for (my $i = 0; $i < @aln; $i++) {
    if ($aln[$i] =~ $id) {
      $start = 1;
    }
    elsif ($aln[$i] =~ />/ and $start == 1) {
      last;
    }
    elsif ($start == 1) {
      $seq .= $aln[$i];
    }
  }
  $seq =~ s/\s//g;
	print join(" ", (caller(0))[0..3]), ", leaving\n" if $debug;
  return ($seq);
}#}}}
##############################
sub revComp {#{{{
	print join(" ", (caller(0))[0..3]), "\n" if $debug;
	my ($seq) = @_;
	$seq =~ tr/AGCTYRKMWSagct/TCGARYMKWSTCGA/;
	$seq = reverse($seq);
	print join(" ", (caller(0))[0..3]), ", leaving\n" if $debug;
	return($seq);
}#}}}
##############################
sub parseHmmer3 {#{{{
	print join(" ", (caller(0))[0..3]), "\n" if $debug;
  my ($file, $path) = @_;
  if (!defined $path) {
    $path = '.';
  }
  open (IN, "$path/$file") or die "failed to open $file\n";
  my @data = <IN>;
  close IN;
  ### extract the hits
  my @hit;
  my $start = 0;
  my $stop = 0;
  my $i = 0;
  for (my $i = 0; $i < @data; $i++) {
    if (!($data[$i] =~ /\S/)) {
      next;
    }
    else {
      if ($data[$i] =~ /Scores for complete sequence/) {
				$start = 1;
				$i += 4;
      }
      elsif (($data[$i] =~ /inclusion threshold/) or ($data[$i] =~ /Domain/i)) {
				last;
      }
      if ($start == 1 and $stop == 0) {
				$data[$i] =~ s/^\s+//;
				my @list = split /\s+/, $data[$i];
				push @hit, $list[8];
				if (@hit == 1 and $bhh) {
					last;
				}
			}
		}
  }
  ### get the query_id
  my ($query) = grep /^Query:/, @data;
  $query =~ s/^Query:\s+//;	#mp over-use of regexes
  $query =~ s/\s.*//;
  if (defined $hit[0]) {
    chomp @hit;
		print join(" ", (caller(0))[0..3]), ", leaving\n" if $debug;
    return ($query, @hit);
  }
  else {
		print join(" ", (caller(0))[0..3]), ", leaving\n" if $debug;
    return ($query);
  }
}#}}}
#####################
#mp this sub reads the entire core-orthologs file into memory!
sub parseSeqfile {#{{{
	print join(" ", (caller(0))[0..3]), "\n" if $debug;
  my $seqref;
  my $id;
  my $spec;
  my $seqid;
  my $seq;
  my $file = shift;
  open (IN, "$file") or die "failed to open $file\n";
  my @seqs = <IN>;
  close IN;
  chomp @seqs;
  for (my $i = 0; $i < @seqs; $i++) {
    if ($seqs[$i] =~ />/) {
			$seqs[$i] =~ s/>//;
      if (defined $id and defined $seq) {
				$seqref->{$id}->{$spec}->{seqid} = $seqid;
				$seqref->{$id}->{$spec}->{seq} = $seq;
				$seq = undef;
      }
      ($id, $spec, $seqid) = split (/\|/, $seqs[$i]);
    }
    else {
      $seq .= $seqs[$i];
    }
  }
  if (defined  $id and defined $seq) {
		$seqref->{$id}->{$spec}->{seqid} = $seqid;
		$seqref->{$id}->{$spec}->{seq} = $seq;
		$seq = undef;
	}
	print join(" ", (caller(0))[0..3]), ", leaving\n" if $debug;
  return ($seqref);
}#}}}
##################
sub getAlignmentScore{ #{{{
	print join(" ", (caller(0))[0..3]), "\n" if $debug;
    my ($refseq_cand, $hitseq) = @_;
    my @testseq = ('>hitseq', $hitseq, '>refseq', $refseq_cand);
    open (OUT, ">$tmpdir/$pid.ref.fa") or die "could not open file for writing refseqs\n";
    print OUT join "\n", @testseq;
    close OUT;
    ## run clustalw
		print "running '$alignmentprog $tmpdir/$pid.ref.fa -output=fasta -outfile=$tmpdir/$pid.ref.aln 2>&1 >$tmpdir/$pid.ref.log'\n" if $debug;	#mp added debug msg
    !(`$alignmentprog $tmpdir/$pid.ref.fa -output=fasta -outfile=$tmpdir/$pid.ref.aln 2>&1 >$tmpdir/$pid.ref.log`) or die "error running clustalw\n";
    ## get the alignment score
    my $score =  `grep "Alignment Score" $tmpdir/$pid.ref.log |sed -e 's/[^0-9]//g'`;
    if (!$score) {
	die "error in determining alignment score! Problem with ClustalW\n";
    }
    chomp $score;
		print join(" ", (caller(0))[0..3]), ", leaving\n" if $debug;
    return ($score);
}#}}}
######################
sub determineRefspecFinal {#{{{
	print join(" ", (caller(0))[0..3]), "\n" if $debug;
  my ($query_name, @refspec) = @_;
  my $refspec_final;
  ## now get the id and the sequence used for building the hmm. Note, the latter will be
  ## needed at a later step to determine the best hit
  my @original;
  my $ac = 0;
  for (my $i = 0; $i < @refspec; $i++) {
		#mp $fafile = core-ortholog seq file
		#mp ($dbfile = EST file)
		#mp $query_name is the HMM ID, e.g., 411851
		#mp remove the core-ortholog taxon name and everything before it; remains: a number at the end of the header
    @original = `grep -A 1 "^>$query_name|$refspec[$i]" $fafile |sed -e "s/.*$refspec[$i]\|//"`;
    chomp @original;
    
    if (@original > 0) {
      $refspec_final->[$ac]->{refspec} = $refspec[$i];
      $refspec_final->[$ac]->{searchdb} = "$blastpath/$refspec[$i]/$refspec[$i]" . "_prot";
      ## now allow for more than one sequence per core-ortholog cluster and species
      $refspec_final->[$ac]->{orthocount} = 0;

      for (my $j = 0; $j < @original; $j+= 2) {
				$refspec_final->[$ac]->{refid}->[$refspec_final->[$ac]->{orthocount}] = $original[$j];
				$refspec_final->[$ac]->{sequence}->[$refspec_final->[$ac]->{orthocount}] = $original[$j+1];
				$refspec_final->[$ac]->{orthocount} += 1;
      }
      $ac++;
      @original = qw();
      if (!defined $strict and !defined $relaxed) {
				## one reftaxon is enough
				last;
      }
    }
    else {
      print "original sequence not be found with grepping for ^>$query_name|$refspec[$i]. Proceeding with next refspec\n";
    }
  }
  if (! defined $refspec_final->[0]->{refid}) {
    print "original sequence not found\n";
		print join(" ", (caller(0))[0..3]), ", leaving\n" if $debug;
    return (0, $refspec_final);
  } 
  ## now print some wordy information...
  if (!defined $strict and !defined $relaxed) {
    print "REFSPEC is $refspec_final->[0]->{refspec}\n";
  }
	print join(" ", (caller(0))[0..3]), ", leaving\n" if $debug;
  return(1, $refspec_final);
}#}}}
######################
#mp sub: save_cdna
#mp saves cdna in cdna output file 
#mp no longer called
sub save_cdna {#{{{
	print join(" ", (caller(0))[0..3]), "\n" if $debug;
	my $self = shift;
	my $hmm_id = shift;
	my $refspec_id = shift;
	my $est_id = shift;
	my $cdnafile = $cdna_dir . '/' . $hmm_id . '.cdna.fa';
	my @cdna_tmp;
	my $cdna_header;
	#mp search for the cdna part of the exonerate output, save in @cdna_tmp
	for (my $i = 0; $i < $self->{gw_count}; $i++) {
		if ($self->{gw}->[$i] =~ />.*_cdna$/) {	#mp corrected regex
			print '$self->{gw}->[$i] is: ' , $self->{gw}->[$i], "\n" if $debug;
			while ($self->{gw}->[$i] !~ m'//') {
				push @cdna_tmp, $self->{gw}->[$i];
				$i++;
			}
			#last; # end the for loop since nothing left to be done - wrong!
			# the whole thing is problematic, though, since it concatenates all the cdna seqs into the same list.
			# this leads to sometimes more than one cdna seq appearing in the cdna output file.
		}
	}
	
	#mp remove all fasta headers - this is equivalent to concatenating the two cdna fragments
	#mp let's just hope this doesn't introduce frameshift errors...
	
	my @cdna = grep( !/^>/, @cdna_tmp );

	#mp cdna fasta header
	$cdna_header = ">$hmm_id\|$refspec_id\|$taxon_global\|$est_id\_cdna";
	#mp open file or exit sub unsuccessfully
	open(my $CDNA_OUTFH, '>', $cdnafile) 
		or return 0;
	print $CDNA_OUTFH join("\n", ($cdna_header, @cdna));
	#mp close file or exit sub unsuccessfully
	close $CDNA_OUTFH
		or return 0;
	print "Saved cdna to $cdnafile\n" if $debug;
	print join(" ", (caller(0))[0..3]), " leaving\n" if $debug;
	return 1;
}#}}}
#mp end sub save_cdna

#mp Documentation moved to the bottom for improved readability

# DATE: Wed Dec 19 10:41:09 CEST 2007#{{{
# Date last modified: 
	##23. 07. 2010: found a bug in the extraction of the
	## hmm hit sequence from the sequence_file. A end-of-line char was missing.
	
	##09.08.2010: added the option to choose the new blastp program from ncbi. Just comment
	##out line 45 in the script and uncomment line 46. Note, in order to make this work I have
	##to slightly modify the blast output since otherwise it will not be parsed by the Bioperl
	##Blast parser. Currently this is a pretty dirty sed hack. It will also take care of removin
	##the string lcl| that is added in some instances to the hit id in the blast output.
	
	## I added the option that one can now provide a comma-separated string of phmm names as an
	## argument for the option -hmm 

	## 08.03.2011: 
	## 1) BUG-FIX: Hamstr will now remove automatically newlines from the input sequence file
	## 2) BUG-FIX: The sequence header remains now the same whether or not the flag -representative
	## has been chosen.
	
	## 10.04.2011
	## 1) added some information to the log file.

	## 20.05.2011
	## 1) BUG-FIX: The grep for the EST sequence in the sub-routine predictORF received also a hit when
	## the search pattern was only a substring of the EST sequence identifier. In some cases the wrong EST
	## was then used to predict the ORF. This has been fixed. 

	## 30.05.2011
	## 1) Extension: a command line option -longhead has been added. The user can now specify that the
	## full sequence id including whitespaces will considered throughout the hamstr search. Note, the
	## whitespaces will be replaced by the string specified in the variabel $idsep.
	## 2) Modification from the bug fix from 20.05.2011. In the grep for the original EST it is no longer
	## necessary that the search string and the EST id are identical over their entire length. Instead the
	## search string may be a prefix of the EST id ending with a whitespace.
	
	## 27.06.2011
	## 1) Extension: I added the option to run a true reciprocal best hit search. Only the best hit from the
	## hmmer search is used to check for reciprocity.#}}}

#mp converted helpmessage to POD

=pod#{{{
=head1 USAGE: hamstrsearch_local.pl -sequence_file=<> -hmmset=<> -taxon=<>  -refspec=<> [-est|-protein] [-hmm=<>] [-representative] [-h]


=head1 OPTIONS:

=head2 -sequence_file: 

path and name of the file containing the sequences hmmer is run against. Per default, this file should be in the data directory.
				
=head2 -est:

set this flag if you are searching in ESTs
				
=head2 -protein: 

set this flag if you are searching in protein sequences
				
=head2 -hmmset: 

specifies the name of the core-ortholog set.  The program will look for the files in the default directory 'core-orthologs' unless you specify a different one.
				
=head2 -taxon:

You need to specify a default taxon name from which your ESTs or protein sequences are derived.
				
=head2 -refspec: 

sets the reference species. Note, it has to be a species that contributed sequences to the hmms you are using. NO DEFAULT IS SET! For a list of possible reference taxa you can have a look at the speclist.txt file in the default core-ortholog sets that come with this distribution. Please use the 5 letter abreviations. If you choose to use core-orthologs were not every taxon is represented in all core-orthologs, you can provide a comma-separated list with the preferred refspec first. The lower-ranking reference species will only be used if a certain gene is not present in the preferred refspecies due to alternative paths in the transitive closure to define the core-orthologs.
				
        CURRENTLY NO CHECK IS IMPLEMENTED!
        NOTE: A BLAST-DB FOR THE REFERENCE SPECIES IS REQUIRED!
				
=head2 -eval_limit=<>

This options allows to set the e-value cut-off for the HMM search.  DEFAULT: 1
				
=head2 -hmm: 

option to provide only a single hmm to be used for the search.  Note, this file has to end with .hmm 
				
=head2 -strict

set this flag if the reciprocity criterion is only fulfilled when the re-blast against all primer taxa was successfull
				
=head2 -relaxed

set this flag if the reciprocity criterion is fulfilled when the re-blast against any of the primer taxa was successfull

=head2 -rbh

set this flag if you want to use a reciprocal best hit criterion. Only the highest scoring hit from the hmmer search will be used for re-blast.

=head2 -append

set this flag if the output should be appended to the files *.out and *_cds.out

=head2 -filter=<T|F>

set this flag to F if the re-blast should be performed without low-complexity filtering. Default is T.

=head2 -longhead

set this flag in the case your sequence identifier contain whitespaces and you whish to keep the entire sequence identifier throughout your analysis. HaMStR will then replace the whitespaces with a '__'. If this flag is not set, HaMStR will truncate the sequence Identifier at the first whitespace, however if and only if the sequence identifier then remain unique.


=head2 ### 

The following options should only be used when you chose to alter the default structure of the hamstrsearch_local directories. Currently, this has not been extensively tested.


=head2 -fasta_file: 

path and name of the file containing the core-ortholog sequences you don't have to use this option when you 

=head2 -hmmpath: 

sets the path to the hmm_dir. By default this is set to the current directory.

=head2 -blastpath: 

sets the path where the blast-dbs are located. Per default this is ../blast_dir. Note, the program expects the structure blastpath/refspec/refspec_prot.  To overrule this structure you will have to edit the script.

=head1 LICENSE

This program is freely distributed under a GPL. See -version for more info.
Copyright (c) GPL limited: portions of the code are from separate copyrights
=cut#}}}
