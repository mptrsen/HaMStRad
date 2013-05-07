#!/usr//bin/env perl
# PROGRAMNAME: hamstrad.pl	#mp

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
# use lib '/dmp-work/dmp/tools/hamstrad/bin';	#mp wtf?
use Getopt::Long;
use Bio::SearchIO;
use Bio::Search::Hit::BlastHit;
use run_genewise_hamstr;
use exonerate;	#mp
use Data::Dumper;	#mp
use File::Spec;	#mp
use File::Basename;	#mp
use Seqload::Fasta; #mp

# PROGRAM DESCRIPTION: See bottom of this file.
######################## start main #############################
my $version = '1.07';	#mp yay
### EDIT THE FOLLOWING LINES TO CUSTOMIZE YOUR SCRIPT
## note: all the variables can also be set via the command line
my $hmmsearchprog = 'hmmsearch'; #program for the hmm search
my $genewise = 'genewise'; #mp
my $exonerate = 'exonerate'; #mp
my $wiseprog = $genewise;	#mp change with -use_exonerate
my $alignmentprog = 'clustalw2';	#mp may be 'clustalw' or 'clustalw2'

########## EDIT THE FOLLOWING TWO LINES TO CHOOSE YOUR BLAST PROGRAM ##########
#my $blast_prog = 'blastall';
my $blast_prog = 'blastp';	#mp may be changed with -blast_prog=<blastprog>
###############################################################################

my $hmmpath = '../core_orthologs'; #path where the hmms are located
my $blastpath = '../blast_dir'; #path to the blast-dbs
my $tmpdir = 'tmp';
my $wiseconfigdir = '../wisecfg';
my $eval = 1; # eval cutoff for the hmm search
my $idsep = '__'; #character used to replace whitespaces in the sequence header with (flag -longhead)

my $hmm_dir = 'hmm_dir';
##############################
my $pid = $$;
my $help;
my $seq2store_file='';
my $cds2store_file='';
my $cds_dir;	#mp
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
my $output_dir;	#mp
my $aa_dir = '';
my $nt_dir = '';	#mp
my $cdnaobj;	#mp
my $gwrefprot;
my $seqtype;
my $align;
my $rep;
my $estflag;
my $proteinflag;
my $refseq;
my $strict;
my $strictstring = '';
my $relaxed;
my $relaxed_string;
my $refspec_final = '';
my $concat;
my $seqs2store_file;
my $append;
my $longhead;
my $check = 1;
my @log = qw();
my $log_dir;	#mp
my $logfile;	#mp
my $filter = 'T';
my $bhh;
my $exhaustive;	#mp
my $debug;	#mp
my $use_exonerate;	#mp
my $ncpu = 1; #mp number of CPU cores that hmmsearch can use
my $exonerate_dir;
my $genewise_dir;
my $skipcount = 0;	#mp
my $couplecount = 0;	#mp
my $showversion = 0;	#mp

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
						"d" => \$debug,	#mp
						"use_exonerate" => \$use_exonerate,	#mp
						"ncpu=i" => \$ncpu, #mp number of CPUs that can be used by hmmsearch
						"blast_prog=s" => \$blast_prog,	#mp
						"clustal_prog=s" => \$alignmentprog,	#mp
						"version|v" => \$showversion,	#mp
);	#mp


if ($help) {
  print $helpmessage;
  exit;
}
#####determine the hostname#######
#mp who the fuck cares for the hostname?
push @log, "$0 version $version\n";
my $hostname = `hostname`;
chomp $hostname;
print "$0 version $version\n";
if ($showversion) { exit }	#mp print only the version if requested
print "hostname is $hostname\n";
push @log, "HOSTNAME\n\t$hostname\n";
#################################
$wiseprog = $use_exonerate ? $exonerate : $genewise; #mp

## 1) check if all information is available to run HaMStR
($check, @log) = &checkInput();
if ($check == 0) {
  print join "\n", @log, "\n";	#mp moved this to the top
  print "$helpmessage";
  print "#######################################\n";
  print "There was an error running $0 $version\n";
  exit;
}
else {
  open (OUT, '>', File::Spec->catfile($log_dir, 'hamstrsearch.log')) or die "could not open logfile: $!\n";	#mp changed logfile
  print OUT join "\n", @log;
  close OUT;
}
### read in of the core-ortholog sequences
#mp literally, read in the entire core-orthologs fasta file.
my $co_seqs = parseSeqfile("$fafile");

#mp read in the original, untranslated sequences
#mp this will hold the entire sequence file in memory for easy retrieval of individual sequences
my $untranslated_sequence_of = slurp_fasta($estfile);


## 2) loop through the hmms
## process each hmm file separately
for (my $i = 0; $i < @hmms; $i++) {
  $fileobj = undef;
  my @seqs = qw();
  my @newseqs = qw();## var to contain the sequences to be added to the orthologous cluster
  my @newcds = qw();
	my @newcdna = qw() if $use_exonerate;	#mp added newcdna array: stores cdna for output
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
    print "no hit found for HMM $query_name\n";	#mp added 'HMM'
    next;
  }
  ### reduce the hit list to only the best entry if $bhh flag is set
#  if ($bhh) {
#	my $keep = shift @results;
#	@results = qw();
#	$results[0] = $keep;
#}
  chomp $query_name;
  print "Results for HMM $query_name\n";	#mp inserted "HMM"
	print "Hits: ", join(", ", @results), "\n" if $debug;	#mp extended report
  my ($check, $refspec_final) = &determineRefspecFinal($query_name, @refspec);
  if ($check == 0) {
    die "error in retrieving refspec data\n";
  }
	HMMER_RESULT:			#mp added label for loop
  for (my $k = 0; $k < @results; $k++) {	#mp for every hmmsearch result, do the following 
    my $hitname = $results[$k];
    print "Processing HMMsearch hit: $hitname\n";	#mp added report
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
					if ($use_exonerate) {
						push @newcdna, ">$query_name|$fileobj->{$taxa[$i]}->{refspec_final}|$taxa[$i]|$fileobj->{$taxa[$i]}->{refid}";
						push @newcdna, $fileobj->{$taxa[$i]}->{refcdna};
					}
					#mp end added cdna output
				}
      }
      else {
				my $idobj = $fileobj->{$taxa[$i]}->{ids};	#mp ids are the seq IDs in the EST file
				my $protobj = $fileobj->{$taxa[$i]}->{prot};
				my $cdsobj  = $fileobj->{$taxa[$i]}->{cds};
				my $refspecobj = $fileobj->{$taxa[$i]}->{refspec};
				for (my $j = 0; $j < @$idobj; $j++) {
					push @newseqs, ">$query_name|$refspecobj->[$j]|$taxa[$i]|$idobj->[$j]";
					push @newseqs, $protobj->[$j];
					if ($estflag) {
						push @newcds, ">$query_name|$taxa[$i]|$idobj->[$j]|$refspecobj->[$j]";
						push @newcds, $cdsobj->[$j];	#mp don't worry everything is fine
					}
				}
      }
      my $refs = $co_seqs->{$query_name};	#mp core ortholog seqs for this query (HMM)
      for (keys %$refs) {	#mp apparently only print cds if -est option was selected
				my $line = ">$query_name|$_|" . $refs->{$_}->{seqid} . "\n" . $refs->{$_}->{seq};
				push @seqs, $line;
      }
      chomp @seqs;
      print "\n";
      @seqs = (@seqs, @newseqs);	#mp y u no use push() like a sane person

			#mp output
      open (OUT, '>', File::Spec->catfile($aa_dir, "$query_name.aa.fa"));	#mp
      print OUT join "\n", @seqs;
      print OUT "\n";
      close OUT;

      if ($estflag) {	#mp apparently only print cds if -est option was selected

				#mp more output
				open (OUT, '>', File::Spec->catfile($cds_dir, $query_name . '.cds.fa')) or die "Fatal: Could not open " . File::Spec->catfile($cds_dir, "$query_name.cds.fa") . ": $!\n";
				print OUT join("\n", @newcds) . "\n";	#mp added newline
				close OUT;

				#mp output cdna to cdna file
				if ($use_exonerate) {
					open (OUT, '>', File::Spec->catfile($nt_dir, $query_name . '.nt.fa')) or die "Fatal: Could not open $nt_dir/$query_name.nt.fa: $!\n";
					print OUT join("\n", @newcdna) . "\n";	#mp added newline
					close OUT;
				}
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
    open (OUT, '>>', $seqs2store_file) or die "failed to open output file $seqs2store_file: $!\n";
  }
  else {
    open (OUT, '>', $seqs2store_file) or die "failed to open output file $seqs2store_file: $!\n";
  }
  print OUT join "\n", @seqs2store;
  print OUT "\n";
  close OUT;
  if ($estflag) {
    if ($append) {
      open (OUT, '>>', $cds2store_file) or die "failed to open output file $cds2store_file: $!\n";
    }
    else {
    open (OUT, '>', $cds2store_file) or die "failed to open output file $cds2store_file: $!\n";
    }
    print OUT join("\n", @cds2store) . "\n";
    print OUT "\n";
    close OUT;
  }
}
else {
  print "no hits found\n";
}
print "$0 $dbfile\: Done!\n";		#mp
print "$skipcount of $couplecount couples skipped during exonerate post-processing.\n"; #mp
exit;

##################### start subs ###############

####### checkInput performs a number of checks whether sufficient information
### and all data are available to run HaMStR
sub checkInput {#{{{
	print join(" ", (caller(0))[0..3]), "\n" if $debug;	#mp
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
			print "Newlines from the infile have been removed.\n";	#mp corrected typo
			push @log, "\tNewlines from the infile have been removed\n";
			$dbfile = $dbfile . '.mod';
			
			if (defined $longhead) {
				`sed -i -e "s/[[:space:]]\\+/$idsep/g" -e 's/\\(>.\\{20\\}\\).*/\\1/' $dbfile`;

				push @log, "\tOption -longhead was chosen. Replaced whitespaces in the sequence identifier with '$idsep' and truncated to 20 characters\n";	#mp added truncation info
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
	print join(" ", (caller(0))[0..3]), ", leaving\n" if $debug;	#mp
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
				#mp I dunno what this is good for, the @info data is never used again
				#--------------------------------------------------
				# open (LOG, '<', File::Spec->catfile($log_dir, 'hamstrsearch.log'));
				# my @info = <LOG>;
				# @log = (@log, @info);
				# close LOG;
				#-------------------------------------------------- 
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
  ## 2) Check for presence of blast program	#mp
  print "Checking for the blast program\t";
  if (system('which', $blast_prog)) {
    push @log, "could not find $blast_prog in PATH. Please check if this program is installed and executable";	#mp
    die "failed: could not find $blast_prog in PATH. Use a different BLAST program with the -blast_prog option.\n";	#mp
    $check = 0;
  }
  else {
    push @log, "check for $blast_prog succeeded";
    print "succeeded\n";
  }
	##mp 2a) Check for presence of alignment program too :P 
  print "Checking for the clustalw program\t";
  if (system('which', $alignmentprog)) {
    push @log, "could not find $alignmentprog in PATH. Please check if this program is installed and executable";	#mp
    die "failed: could not find $alignmentprog in PATH. Choose a different clustalw program with the -clustal_prog option.\n";	#mp
    $check = 0;
  }
  else {
    push @log, "check for $alignmentprog succeeded";
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
      push @log, "check for $hmmsearchprog succeeded";
      print "succeeded\n";
  }
	##mp 3a) How about we check for genewise/exonerate as well? Did you forget that?
	print "Checking for $wiseprog\t";
	if (! `$wiseprog -h`) {
		push @log, "could not execute $wiseprog. Please check if this program is installed and executable";
		print "failed\n";
		$check = 0;
	}
	else {
		push @log, "check for $wiseprog succeeded";
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
  if (defined $strict) {
        $strictstring = '.strict';
}

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
  $hmmsearchprog = $hmmsearchprog . " -E $eval --cpu $ncpu";
  push @log, "hmmsearch: $hmmsearchprog";

  ## 13) setting up the directories where the output files will be put into.
	#mp changed output dir structure
	$strictstring = $strict ? 'strict_' : '';
	$relaxed_string = defined($relaxed) ? 'relaxed_' : '';

  $output_dir = File::Spec->catdir($dbfile_short . '_' . $hmmset . '_' . $strictstring . $relaxed_string . join('_', @refspec));	#mp 
	push @log, "Using output directory $output_dir";	#mp
	print "Using output directory $output_dir\n";	#mp 

	$aa_dir = File::Spec->catdir($output_dir, 'aa');	#mp File::Spec
	$nt_dir = File::Spec->catdir($output_dir, 'nt');	#mp File::Spec
	$cds_dir = File::Spec->catdir($nt_dir, 'cds');	#mp File::Spec
	$log_dir = File::Spec->catdir($output_dir, 'log');	#mp File::Spec
	#mp overwrite the tmpdir variable with a per-run tmpdir to avoid confusion
	$tmpdir = File::Spec->catdir($output_dir, 'tmp'); #mp added per-run tempdir
  $hmmsearch_dir = File::Spec->catdir($log_dir, 'hmmsearch');	#mp File::Spec
	$exonerate_dir = File::Spec->catdir($log_dir, 'exonerate') if $use_exonerate;	#mp File::Spec
	$genewise_dir = File::Spec->catdir($log_dir, 'genewise') unless $use_exonerate;	#mp added genewise output dir
  $seqs2store_file = File::Spec->catfile($log_dir, 'hamstrsearch_' . basename($dbfile_short) . '_' . $hmmset . $strictstring . '.out');	#mp File::Spec
  $cds2store_file = File::Spec->catfile($log_dir, 'hamstrsearch_' . basename($dbfile_short) . '_' . $hmmset . '_cds' . $strictstring . '.out');	#mp File::Spec
  if ($check == 1) {
    if (!(-e "$hmmsearch_dir")) {
      `mkdir -p $hmmsearch_dir`;	#mp added -p flag to mkdir
    }
    if (!(-e "$aa_dir")) {
      `mkdir -p $aa_dir`;	#mp added -p flag to mkdir
    }
		if (!(-e "$nt_dir")) {
			`mkdir -p $nt_dir`;	#mp with -p flag
		}
		if (!(-e "$cds_dir")) {
			`mkdir -p $cds_dir`;	#mp with -p flag
		}
		#mp added creation of additional output dirs
		if ($use_exonerate) {
			if (!(-e $exonerate_dir)) {
				`mkdir -p $exonerate_dir`;	#mp with -p flag
			}
		}
		#mp added genewise output dir
		if (!$use_exonerate) {
			if (!(-e "$genewise_dir")) {
				`mkdir -p $genewise_dir`;	#mp with -p flag
			}
		}
		#mp end added genewise output dir
		if (!(-e "$log_dir")) {
			`mkdir -p $log_dir`;	#mp with -p flag
		}
		if (!(-e "$tmpdir")) {
			`mkdir -p $tmpdir`;	#mp with -p flag
		}
		#mp end add creation of additional output dirs
		#mp end changed output dir structure
    if (!(-e "$tmpdir")) {
      `mkdir -p $tmpdir`;	#mp added -p flag to mkdir
    }
		#--------------------------------------------------
		# #mp not used, candidate for removal
    # if (!(-e "hamstrsearch_$dbfile_short")) {
		# 	`mkdir -p hamstrsearch_$dbfile_short`; #mp added mkdir -p for this dir
		# }
		# #mp end candidate for removal
		#-------------------------------------------------- 
  }
  ## 14) determin whether or not the -representative flag has been set
  if (defined $rep) {
	push @log, "HaMStR will run with the -representative option";
	}
  else {
	push @log, "HaMStR was called without the -representative option. More than one ortholog may be identified per core-ortholog group!";
	} 
	if (defined $use_exonerate) {
		push @log, "using exonerate instead of genewise.";
		print "Using exonerate instead of genewise.\n";
	}
	print "Finished checking input.\n" if $debug;	#mp
	print join(" ", (caller(0))[0..3]), ", leaving\n" if $debug;	#mp
  return ($check, @log);
}#}}}

#################
## check4reciprocity is the second major part of the program. It checks
## whether the protein sequence that has been identified by the hmmsearch
## identifies in turn the protein from the reference taxon that was used to
## build the hmm.
sub check4reciprocity {#{{{
	print join(" ", (caller(0))[0..3]), "\n" if $debug;	#mp
  my ($query_name, $hitname, $refspec_final, @refspec) = @_;
  my $searchdb;
  my $strict_suc = -1; # keeps track of success for all taxa
  my $relaxed_suc = 0; # keeps track of success for at least one taxon
  ## get the sequence that was identified as hit in the pHMM search from the db_file
	print "Finding '$hitname' in $dbfile...\n" if $debug; #mp
  my $hitseq = `grep -F -m 1 -A 1 ">$hitname" $dbfile | tail -n 1`;	#mp added -F
  if (!defined $hitseq) {
    print "could not retrieve a sequence for $hitname. Skipping...\n";
		print join(" ", (caller(0))[0..3]), ", leaving\n" if $debug;	#mp
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
      !`$blast_prog -p blastp -d $refspec_final->[$k]->{searchdb} -F $filter -i $tmpdir/$$.fa -o $tmpdir/$$.blast` or die "Problem running blastall\n";	#mp
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
						print "\thitting...\n";
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
			print join(" ", (caller(0))[0..3]), ", leaving\n" if $debug;	#mp
      return (0, $hitseq);
    }
  }

  if ($relaxed_suc == 1) {
		print join(" ", (caller(0))[0..3]), ", leaving\n" if $debug;	#mp
    return (1, $hitseq);
  }
  else {
		print join(" ", (caller(0))[0..3]), ", leaving\n" if $debug;	#mp
    return (0, $hitseq);
  }
}#}}}
#############

sub getBestBlasthit {#{{{
	print join(" ", (caller(0))[0..3]), "\n" if $debug;	#mp
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
		print join(" ", (caller(0))[0..3]), ", leaving\n" if $debug;	#mp
    return(@hits);
}#}}}
##################
sub getTaxon {#{{{
	print join(" ", (caller(0))[0..3]), "\n" if $debug;	#mp
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
			print join(" ", (caller(0))[0..3]), ", leaving\n" if $debug;	#mp
			return ($taxon);
    }
    else {
			print join(" ", (caller(0))[0..3]), ", leaving\n" if $debug;	#mp
			return();
    }
}#}}}
###############
sub determineReferences {#{{{
	print join(" ", (caller(0))[0..3]), "\n" if $debug;	#mp
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
					print "Alignment score: $score: setting refspec to $refspec_cand\n";
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
		print 'fileobj: ' . Dumper($fileobj) if $debug;	#mp
		print "refspec for query ID $hitname is empty!\nThis means that this orthology assignment has an alignment score of zero or less (bad alignment). We are not sure whether this is a bug, but it will probably lead to genewise/exonerate crashing later.\n";
		print "\n";
		open(my $errfh, ">>$errorfile") || die "Could not open $errorfile: $!\n";
		print $errfh scalar(localtime), " Bad ClustalW alignment: Missing refspec for $hitname from hmmsearching with $query_name\n";
		close $errfh || die "Could not close $errorfile: $!\n";
	}
	else {
		print "refspec for query ID \"$hitname\" is \"$refspec\"\n";
	}
	#--------------------------------------------------
	# #mp end add WARNING condition
	#-------------------------------------------------- 
	print join(" ", (caller(0))[0..3]), ", leaving\n" if $debug;	#mp
	return($fileobj);
}#}}}
###############
sub processHits {#{{{
	print join(" ", (caller(0))[0..3]), "\n" if $debug;	#mp
  my ($fileobj) = @_; 
  ## 1) align all hit sequences for a taxon against the reference species
  my @taxa = keys(%$fileobj);
#mp let's just test what happens if we only align ONE hit against the refspec
  for (my $i = 0; $i < @taxa; $i++) {
		if ($fileobj->{$taxa[$i]}->{prot}) {
			print "$taxa[$i] has prot defined, doing orfRanking...\n" if $debug; #mp
			&orfRanking($taxa[$i]);
			#last;	#mp TODO is this right ?
		}
  }
}  #}}}
  

################
sub predictORF {#{{{
	print join(" ", (caller(0))[0..3]), "\n" if $debug;	#mp
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

		#print Dumper($fileobj->{$taxa[$i]}); exit;	#mp debugging

		ID:	#mp labelled the loop. These are the individual sequences from the EST file that produced a hmmsearch hit. For each, do:
    for (my $j = 0; $j < @ids; $j++) {
			my $refseq = $refseqobj->[$j];
			my $refspec = $refspecobj->[$j];
			## determine the reading frame and remove the info from the ID
			my ($rf) = $ids[$j] =~ /.*_RF([0-9]+)/;
			print "Predicting ORF for hit $ids[$j] (RF $rf)\n";
			$ids[$j] =~ s/_RF\d\.\d$//;	#mp modified regex to avoid problems with headers ending on '_RF'
			#mp removed "\\b" from regex and added info line
			#mp this regex would almost never work because headers are not trimmed anymore
			#mp
			#mp TODO Problem: Sonderzeichen vs. \b
			#mp added -F flag for fixed-string search. potentially dangerous :/
			#mp also this needs -m 1 because for multiple matches the last one would be used, 
			#mp which may not always be the correct one
			#my $est = `grep -F -m 1 -A 1 ">$ids[$j]" $estfile |tail -n 1`;	#mp get the sequence from the original EST file
			print "Fetching original sequence for $ids[$j] from EST file...\n" if $debug; #mp
			my $est = $untranslated_sequence_of->{$ids[$j]};
			if (! $est) {
				die "error in retrieval of est sequence for $ids[$j] in subroutine processHits\n";
			}
			print "found $ids[$j] in $estfile\n" if $debug;	#mp
			## the EST is translated in rev complement
			if ($rf > 3) {
				$est = revComp($est);
			}
			# TODO why are $refspec and $refseq sometimes left empty?
			#mp run either exonerate or genewise, depending on invocation
			++$couplecount;	#mp
			my $gw;
			if ($use_exonerate) {
				$gw = exonerate->new($est, $refseq, "$tmpdir");
				#mp save exonerate output to file
				my $gwoutfile = File::Spec->catfile($exonerate_dir, $query_name . "_" . $refspec . '_' . "$ids[$j]" . '.' . basename($wiseprog) . "out");	#mp File::Spec
				my $gwoutput = $gw->{gw};
				open(my $gwresultfh, '>', $gwoutfile) or die "Could not open $gwoutfile for writing: $!\n";
				print $gwresultfh join("\n", @$gwoutput) . "\n";
				close $gwresultfh or die "Could not close file $gwoutfile: $!\n";
				print "Wrote exonerate output to $gwoutfile\n" if $debug;	#mp
				#mp end save exonerate output 
				#mp Skip a couple if exonerate doesn't return anything 
				#mp most likely if it has been handed an empty prot sequence (for whatever reason)
				if ($gw->{gw_count} == 0) {
					++$skipcount;
					print "$ids[$j] and $query_name|$refspec returned an empty exonerate result, skipping this couple ($skipcount skipped so far).\n";
					next TAXON;
				}
				#mp end skip couple
			}
			else {
				$gw = run_genewise_hamstr->new($est, $refseq, "$tmpdir");
				#mp save genewise output to file
				my $gwoutfile = File::Spec->catfile($tmpdir, $query_name . "_" . $refspec . '_' . "$ids[$j]" . '.' . basename($wiseprog) . "out");	#mp File::Spec
				my $gwoutput = $gw->{gw};
				open(my $gwresultfh, '>', $gwoutfile) or die "Could not open $gwoutfile for writing: $!\n";
				print $gwresultfh join("\n", @$gwoutput) . "\n";
				close $gwresultfh or die "Could not close file $gwoutfile: $!\n";
				print "Wrote genewise output to $gwoutfile\n" if $debug;	#mp
				#mp end save genewise output 
			}

			my $translation = $gw->translation;
			my $cds = $gw->codons;
			$translation =~ s/[-!]//g;	#mp deletes gaps and stop codons (but not those coded with '*')
			$fileobj_new->{$taxa[$i]}->{ids}->[$j] = $ids[$j];
			$fileobj_new->{$taxa[$i]}->{prot}->[$j] = $translation;
			$fileobj_new->{$taxa[$i]}->{cds}->[$j] = $cds;
			$fileobj_new->{$taxa[$i]}->{cdna}->[$j] = $gw->extract_cdna if $use_exonerate;	#mp 
			$fileobj_new->{$taxa[$i]}->{refseq}->[$j] = $refseq;
			$fileobj_new->{$taxa[$i]}->{refspec}->[$j] = $refspec;
		}
  }
	print join(" ", (caller(0))[0..3]), ", leaving\n" if $debug;	#mp added debug info
  return($fileobj_new);
}#}}}
############################
sub orfRanking {#{{{
	print join(" ", (caller(0))[0..3]), "\n" if $debug;	#mp
  my ($spec) = @_;
  my $result;
  my $refprot;
  my $refcds;
	my $refcdna if $use_exonerate;	#mp added $refcdna
  my @toalign;
  my $protobj = $fileobj->{$spec}->{prot};
  my $idobj = $fileobj->{$spec}->{ids};

  my $refcluster; ## variables to take the cluster and its id for later analysis
  my $refid;
  if (scalar @$protobj == 1) {	#mp added 'scalar' for clarity
    ## nothing to chose from
    $refprot = $protobj->[0];
    $refcds = $fileobj->{$spec}->{cds}->[0];
		$refcdna = $fileobj->{$spec}->{cdna}->[0] if $use_exonerate;	#mp added refcdna
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
		for (my $i = 0; $i < scalar @$protobj; $i++) {	#mp added 'scalar' for clarity
			my @testseq = (">$idobj->[$i]", $protobj->[$i]);
			@testseq = (@testseq, @toalign);	#mp y u no use push() for readability

			open (OUT, ">$tmpdir/$pid.ref.fa") or die "could not open file for writing refseqs\n";
			print OUT join "\n", @testseq;
			close OUT;

			## run clustalw	#mp
			print "running '$alignmentprog -infile=$tmpdir/$pid.ref.fa -output=fasta -outfile=$tmpdir/$pid.ref.aln 2>&1 >$tmpdir/$pid.ref.log'\n" if $debug;	#mp added debug msg
			!(`$alignmentprog -infile=$tmpdir/$pid.ref.fa -output=fasta -outfile=$tmpdir/$pid.ref.aln 2>&1 >$tmpdir/$pid.ref.log`) or die "error running $alignmentprog\: $!\n";	#mp added error message
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

		#mp added $refcdna -- THIS IS where the seq fragments got concatenated
		if ($use_exonerate) {
			($refprot, $refcds, $refcdna, $refid) = &determineRef($result,$spec);	
		}
		else {
			($refprot, $refcds, $refid) = &determineRef($result,$spec);	#mp 
		}
		#mp end added $refcdna
  }
  $fileobj->{$spec}->{refprot} = $refprot;
  $fileobj->{$spec}->{refcds}  = $refcds;
  $fileobj->{$spec}->{refcdna}  = $refcdna if $use_exonerate;	#mp added refcdna
  $fileobj->{$spec}->{refid}   = $refid;
  $fileobj->{$spec}->{refspec_final} = $fileobj->{$spec}->{refspec}->[0];
	print join(" ", (caller(0))[0..3]), ", leaving\n" if $debug;	#mp
  return();
}#}}}
###########################
sub sortRef {#{{{
	print join(" ", (caller(0))[0..3]), "\n" if $debug;	#mp
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
	print join(" ", (caller(0))[0..3]), ", leaving\n" if $debug;	#mp
	return($result);
}#}}}
########################
sub determineRef {#{{{
	print join(" ", (caller(0))[0..3]), "\n" if $debug;	#mp
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
	my $refcdna = '' if $use_exonerate;	#mp added refcdna, maybe this helps
	my @headers;	#mp
  for (my $i = 0; $i < @$final; $i++) {

		#mp add check for duplicate transcripts
		my $header = $fileobj->{$spec}->{ids}->[$final->[$i]->{id}];
		if (grep(/$header/, @headers)) {
			printf("Skipping %s: already present in this sequence\n", $header);
			next;
		}
		push(@headers, $header);
		#mp end check for duplicate transcripts

    my $seq = $fileobj->{$spec}->{prot}->[$final->[$i]->{id}];
    my $cdsseq = $fileobj->{$spec}->{cds}->[$final->[$i]->{id}];
		my $cdnaseq = $fileobj->{$spec}->{cdna}->[$final->[$i]->{id}] if $use_exonerate;	#mp added cdnaseq
    my $length = length($seq);
    $refid .= "$fileobj->{$spec}->{ids}->[$final->[$i]->{id}]-$length" . "PP";
    $refprot .= $seq;
    if ($estflag) {
      $refcds .= $cdsseq;
			$refcdna .= $cdnaseq if $use_exonerate;	#mp added refcdna
    }
  }
	#mp remove trailing 'PP'
  $refid =~ s/PP$//;
	print join(" ", (caller(0))[0..3]), ", leaving\n" if $debug;	#mp
	#mp added returning $refcdna if using exonerate
	if ($use_exonerate) {
		return($refprot, $refcds, $refcdna, $refid);	
	}
	else {
		return($refprot, $refcds, $refid);	
	}
	#mp end added returning $refcdna
}#}}}
#############################
sub extractSeq {#{{{
	print join(" ", (caller(0))[0..3]), "\n" if $debug;	#mp
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
	print join(" ", (caller(0))[0..3]), ", leaving\n" if $debug;	#mp
  return ($seq);
}#}}}
##############################
sub revComp {#{{{
	print join(" ", (caller(0))[0..3]), "\n" if $debug;	#mp
	my ($seq) = @_;
	$seq =~ tr/AGCTYRKMWSagct/TCGARYMKWSTCGA/;
	$seq = reverse($seq);
	print join(" ", (caller(0))[0..3]), ", leaving\n" if $debug;	#mp
	return($seq);
}#}}}
##############################
sub parseHmmer3 {#{{{
	print join(" ", (caller(0))[0..3]), "\n" if $debug;	#mp
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
		print join(" ", (caller(0))[0..3]), ", leaving\n" if $debug;	#mp
    return ($query, @hit);
  }
  else {
		print join(" ", (caller(0))[0..3]), ", leaving\n" if $debug;	#mp
    return ($query);
  }
}#}}}
#####################
#mp this sub reads the entire core-orthologs file into memory!
sub parseSeqfile {#{{{
	print join(" ", (caller(0))[0..3]), "\n" if $debug;	#mp
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
	print join(" ", (caller(0))[0..3]), ", leaving\n" if $debug;	#mp
  return ($seqref);
}#}}}
##################
sub getAlignmentScore{ #{{{
	print join(" ", (caller(0))[0..3]), "\n" if $debug;	#mp
    my ($refseq_cand, $hitseq) = @_;
    my @testseq = ('>hitseq', $hitseq, '>refseq', $refseq_cand);
    open (OUT, ">$tmpdir/$pid.ref.fa") or die "could not open file for writing refseqs\n";
    print OUT join "\n", @testseq;
    close OUT;
    ## run clustalw
		print "running '$alignmentprog -infile=$tmpdir/$pid.ref.fa -output=fasta -outfile=$tmpdir/$pid.ref.aln 2>&1 >$tmpdir/$pid.ref.log'\n" if $debug;	#mp added debug msg 
    !(`$alignmentprog -infile=$tmpdir/$pid.ref.fa -output=fasta -outfile=$tmpdir/$pid.ref.aln 2>&1 >$tmpdir/$pid.ref.log`) or die "error running clustalw\n";	#mp added -infile=
    ## get the alignment score
    my $score =  `grep "Alignment Score" $tmpdir/$pid.ref.log |sed -e 's/[^0-9]//g'`;
    if (!$score) {
	die "error in determining alignment score! Problem with ClustalW\n";
    }
    chomp $score;
		print join(" ", (caller(0))[0..3]), ", leaving\n" if $debug;	#mp
    return ($score);
}#}}}
######################
sub determineRefspecFinal {#{{{
	print join(" ", (caller(0))[0..3]), "\n" if $debug;	#mp
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
		print join(" ", (caller(0))[0..3]), ", leaving\n" if $debug;	#mp
    return (0, $refspec_final);
  } 
  ## now print some wordy information...
  if (!defined $strict and !defined $relaxed) {
    print "REFSPEC is $refspec_final->[0]->{refspec}\n";
  }
	print join(" ", (caller(0))[0..3]), ", leaving\n" if $debug;	#mp
  return(1, $refspec_final);
}#}}}
######################
#mp sub: save_cdna
#mp saves cdna in cdna output file 
#mp no longer called, candidate for removal
sub save_cdna {#{{{
	print join(" ", (caller(0))[0..3]), "\n" if $debug;	#mp
	my $self = shift;
	my $hmm_id = shift;
	my $refspec_id = shift;
	my $est_id = shift;
	my $cdnafile = $nt_dir . '/' . $hmm_id . '.cdna.fa';
	my @cdna_tmp;
	my $cdna_header;
	#mp search for the cdna part of the exonerate output, save in @cdna_tmp
	for (my $i = 0; $i < $self->{gw_count}; $i++) {
		if ($self->{gw}->[$i] =~ />.*_cdna$/) {	#mp corrected regex
			print '$self->{gw}->[$i] is: ' , $self->{gw}->[$i], "\n" if $debug;	#mp
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
	print "Saved cdna to $cdnafile\n" if $debug;	#mp
	print join(" ", (caller(0))[0..3]), " leaving\n" if $debug;	#mp
	return 1;
}#}}}
#mp end sub save_cdna

#mp sub: slurp_fasta
#mp reads the content of a Fasta file into a hashref
sub slurp_fasta {
	my $fastafile = shift @_;
	my $data = { };
	my $fastafh = Seqload::Fasta->open($fastafile);
	while (my ($h, $s) = $fastafh->next_seq()) {
		$data->{$h} = $s;
	}
	return $data;
}

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

#{{{
=head1 NAME: 

HaMStRad - HaMStR advanced, a pipeline for orthology assessment based on EST data


=head1 SYNOPSIS: 

hamstrad.pl -sequence_file=<> -hmmset=<> -taxon=<>  -refspec=<> [-est|-protein] [-hmm=<>] [-representative] [-h]	


=head1 OPTIONS:

=head2 -sequence_file=NAME 

path and name of the file containing the sequences hmmer is run against. Per default, this file should be in the data directory.
				
=head2 -est

set this flag if you are searching in ESTs
				
=head2 -protein 

set this flag if you are searching in protein sequences
				
=head2 -hmmset=NAME 

specifies the name of the core-ortholog set.  The program will look for the files in the default directory 'core-orthologs' unless you specify a different one.
				
=head2 -taxon=NAME

You need to specify a default taxon name from which your ESTs or protein sequences are derived.
				
=head2 -refspec=NAME 

sets the reference species. Note, it has to be a species that contributed sequences to the hmms you are using. NO DEFAULT IS SET! For a list of possible reference taxa you can have a look at the speclist.txt file in the default core-ortholog sets that come with this distribution. Please use the 5 letter abreviations. If you choose to use core-orthologs were not every taxon is represented in all core-orthologs, you can provide a comma-separated list with the preferred refspec first. The lower-ranking reference species will only be used if a certain gene is not present in the preferred refspecies due to alternative paths in the transitive closure to define the core-orthologs.
				
        CURRENTLY NO CHECK IS IMPLEMENTED!
        NOTE: A BLAST-DB FOR THE REFERENCE SPECIES IS REQUIRED!
				
=head2 -eval_limit=<>

This options allows to set the e-value cut-off for the HMM search.  DEFAULT: 1
				
=head2 -hmm=HMMFILE 

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

=head2 -use_exonerate

use exonerate instead of genewise. This also enables corresponding nucleotide output.

=head2 -ncpu <N>

set number of CPU cores that hmmsearch can use. Defaults to 1 if left unspecified.

=head2 -blast_prog=NAME

sets the name of the BLAST program. May be 'blastp' or 'blastall'. Default: blastp

=head2 -clustal_prog=NAME

sets the name of the clustalw program. May be 'clustalw' or 'clustalw2'. Default: clustalw2


=head1 The following options should only be used when you chose to alter the default structure of the hamstrad directories. Currently, this has not been extensively tested.

=head2 -fasta_file=PATH

path and name of the file containing the core-ortholog sequences you don't have to use this option when you 

=head2 -hmmpath=PATH

sets the path to the hmm_dir. By default this is set to the current directory.

=head2 -blastpath=PATH

sets the path where the blast-dbs are located. Per default this is ../blast_dir. Note, the program expects the structure blastpath/refspec/refspec_prot.  To overrule this structure you will have to edit the script.

=head1 LICENSE

This program is freely distributed under a GPL. See -h for more info.
Copyright (c) GPL limited: portions of the code are from separate copyrights
=cut#}}}
