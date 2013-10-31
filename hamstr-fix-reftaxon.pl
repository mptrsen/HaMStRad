#!/usr/bin/perl
=head1 NAME

hamstr-fix-reftaxon.pl

=head1 SYNOPSIS

perl hamstr-fix-reftaxon.pl --inputdir INPUTDIR --outputdir OUTPUTDIR

=head1 DESCRIPTION

hamstr-fix-reftaxon attempts to fix the reference taxon assignment that HaMStR messed up. It does so by generating a pairwise clustalw alignment and selecting the reference taxon sequence with the highest score. It then rewrites the output files using the correct header.

=head1 OPTIONS

=head2 --inputdir DIR

Input directory. Is expected to contain HaMStR output files (*.aa.fa). If you do not specify an input directory, filenames are read from the argument list (command line).

=head2 --outputdir DIR

Output directory. Where the rewritten files should be placed. It is recommended not to use the same dir that you read from, but instead overwrite the affected files with their corrected copies afterwards. If you do not specify an output directory, the input files will be corrected in-place, i.e., safely overwritten.

=head2 --clustalw /path/to/clustalw

Path to clustalw executable. On some systems, this is called clustalw2, or is not in $PATH, or something else.

=head1 AUTHOR

Malte Petersen <mptrsen@uni-bonn.de>

=cut

use strict;
use warnings;
use autodie;

use File::Basename;
use File::Copy;
use File::Find;
use File::Path qw( make_path remove_tree ); # this also uses File::Spec
use File::Temp;
use Getopt::Long;
use IO::Dir;
use IO::File;
use List::Util qw(first);
use Data::Dumper;

my $clustalw = 'clustalw';
my $outdir = undef;
my $aaoutdir = undef;
my $ntoutdir = undef;
my $logoutdir = undef;
my $aaindir = undef;
my $ntindir = undef;
my $logindir = undef;
my $indir = undef;
my @infiles = ();

GetOptions(
	'clustalw=s' => \$clustalw,
	'outputdir=s' => \$outdir,
	'inputdir=s' => \$indir,
);

if ($indir) {
	# get list of input files
	$aaindir = File::Spec->catdir($indir, 'aa');
	$ntindir = File::Spec->catdir($indir, 'nt');
	$logindir = File::Spec->catdir($indir, 'log');
	unless (-d $aaindir) { die "Fatal: not a directory: $aaindir\n" }
	unless (-d $ntindir) { die "Fatal: not a directory: $ntindir\n" }
	@infiles = get_files($aaindir);
	if (scalar @infiles == 0) { die "Fatal: empty directory $aaindir\n" }
}
else {
	print "Usage: $0 --inputdir INPUTDIR --outputdir OUTPUTDIR\n" and exit;
}

if ($outdir) {
	$aaoutdir = File::Spec->catdir($outdir, 'aa');
	unless (-d $aaoutdir) { make_path($aaoutdir) or die "Fatal: could not create $aaoutdir: $!\n" }
	$ntoutdir = File::Spec->catdir($outdir, 'nt');
	unless (-d $ntoutdir) { make_path($ntoutdir) or die "Fatal: could not create $ntoutdir: $!\n" }
	$logoutdir = File::Spec->catdir($outdir, 'log'); 
	unless (-d $logoutdir) { make_path($logoutdir) or die "Fatal: could not create $logoutdir: $!\n" }
}
else {
	print "Usage: $0 --inputdir INPUTDIR --outputdir OUTPUTDIR\n" and exit;
}

my ($logfile, $cdslogfile) = get_logfiles($logindir);
my $logfh = IO::File->new(File::Spec->catfile($logoutdir, basename($logfile)), 'w');
my $cdslogfh = IO::File->new(File::Spec->catfile($logoutdir, basename($cdslogfile)), 'w');

foreach my $inf (@infiles) {
	# read sequences into memory
	my $seq_of = slurp_fasta($inf);

	# the taxon header
	my $header = first { /.*\|.*\|.*\|.*/ } keys %$seq_of;
	my ($geneid, $coretaxon, $taxon, $id) = split /\|/, $header;

	# take it out of the sequence pool
	my $sequence = delete $seq_of->{$header};

	# determine the correct reftaxon and take that out of the pool as well
	my ($hiscoretaxon, $hiscoreheader) = get_real_coretaxon($seq_of, $header, $sequence);
	my $hiscoresequence = delete $seq_of->{$hiscoreheader};
	# the hash now contains only all sequences without the refspec sequence and the taxon sequence
	# this is so we can write those back in order

	# format the header
	my $new_header = sprintf "%s|%s|%s|%s", $geneid, $hiscoretaxon, $taxon, $id;

	# output
	my $aaoutf = undef;
	my $ntoutf = undef;
	my $outfh = undef;

	$aaoutf = File::Spec->catfile($aaoutdir, basename($inf));

	$outfh = IO::File->new($aaoutf, 'w') or die "Fatal: could not open $aaoutf for writing: $!\n";
	printf $outfh ">%s\n%s\n", $_, $seq_of->{$_} foreach keys %$seq_of;
	printf $outfh ">%s\n%s\n", $hiscoreheader, $hiscoresequence;
	printf $outfh ">%s\n%s\n", $new_header, $sequence;
	undef $outfh;

	# rewrite the nucleotide file correspondingly
	rewrite_nucfile($inf, $ntoutdir, $header, $hiscoretaxon);

	printf $logfh "%s|%s\n", $new_header, $sequence;
	# report
	printf "%s: reftaxon for %s is %s (was: %s)\n", basename($aaoutf), $taxon, $hiscoretaxon, $coretaxon;
}


sub rewrite_nucfile {
	my $f = basename( shift @_ );
	my $outdir = shift @_;
	my ($header, $hiscoretaxon) = @_;
	$f =~ s/\.aa\./.nt./;
	my $outf = File::Spec->catfile($outdir, $f);

	# it's only a single sequence
	my $ntseqs = slurp_fasta( File::Spec->catfile($ntindir, $f) );
	my $sequence = delete $ntseqs->{$header};

	# generate the real header
	my ($geneid, $coretaxon, $taxon, $id) = split /\|/, $header;
	my $real_header = sprintf "%s|%s|%s|%s", $geneid, $hiscoretaxon, $taxon, $id;

	# write the sequence to file
	my $fh = IO::File->new(File::Spec->catfile($outf), 'w') or die "Fatal: could not open $outf for writing: $!\n";
	printf $fh ">%s\n%s\n", $real_header, $sequence;
	undef $fh;

	printf $cdslogfh "%s|%s\n", $real_header, $sequence;
}

sub get_logfiles {
	my $indir = shift @_;
	my @logfiles = ();
	find( sub { push @logfiles, $File::Find::name if /^hamstr.*\.out$/ }, $indir);
	if (scalar @logfiles != 2) {
		die "Fatal: missing log files in $indir\n";
	}
	my $cdslogfile = first { /(cds\.out|cdsstrict_\.out)$/ } @logfiles;
	my $logfile = first { ! /(cds\.out|cdsstrict_\.out)$/ } @logfiles;
	return ($logfile, $cdslogfile);
}

# get a list of (*.fa) files in the the dir
# call: get_files($dirname)
# returns: list of scalar string filenames
sub get_files {
	my $dirn = shift @_;
	my $dirh = IO::Dir->new($dirn);
	die "Fatal: could not open dir $dirn\: $!\n" unless defined $dirh;
	my @files = ();
	while (my $f = $dirh->read) {
		# skip stuff starting with a dot
		next if $f =~ /^\./;
		if (-f File::Spec->catfile($dirn, $f)) {
			push @files, File::Spec->catfile($dirn, $f);
		}
	}
	undef $dirh;
	@files = grep { /\.fa$/ } @files;
	return @files;
}

sub get_real_coretaxon {
	my $sequences = shift @_;
	my $header = shift @_;
	my $sequence = shift @_;
	my $hiscore = 0;
	my $hiscoretaxon = '';
	my $hiscoreheader = '';

	foreach (keys %$sequences) {
		my ($geneid, $taxon, $id) = split /\|/;
		my $fh = File::Temp->new();
		printf $fh ">%s\n%s\n", $_, $sequences->{$_};
		printf $fh ">%s\n%s\n", $header, $sequence;
		close $fh;
		my $outfh = File::Temp->new();
		my $result = [ `$clustalw -infile=$fh -outfile=$outfh` ];
		chomp @$result;
		my $score = first { $_ =~ /Alignment Score/ } @$result;
		$score =~ /Score\s*(-?\d+)/ and $score = $1;
		if ($score > $hiscore) {
			$hiscore = $score;
			$hiscoretaxon = $taxon;
			$hiscoreheader = $_;
		}
	}

	return ($hiscoretaxon, $hiscoreheader);
}

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


####################
# 
# Seqload::Fasta package
# for simple and error-free loading of fasta sequence data
# 
####################

package Seqload::Fasta;
use strict;
use warnings;
use Carp;
require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw( fasta2csv check_if_fasta );

# Constructor. Returns a sequence database object.
sub open {
  my ($class, $filename) = @_;
  open (my $fh, '<', $filename)
    or confess "Fatal: Could not open $filename\: $!\n";
  my $self = {
    'filename' => $filename,
    'fh'       => $fh
  };
  bless($self, $class);
  return $self;
}

# Returns the next sequence as an array (hdr, seq). 
# Useful for looping through a seq database.
sub next_seq {
  my $self = shift;
  my $fh = $self->{'fh'};
	# this is the trick that makes this work
  local $/ = "\n>"; # change the line separator
  return unless defined(my $item = readline($fh));  # read the line(s)
  chomp $item;
  
  if ($. == 1 and $item !~ /^>/) {  # first line is not a header
    croak "Fatal: " . $self->{'filename'} . "is not a FASTA file: Missing descriptor line\n";
  }

	# remove the '>'
  $item =~ s/^>//;

	# split to a maximum of two items (header, sequence)
  my ($hdr, $seq) = split(/\n/, $item, 2);
	$hdr =~ s/\s+$//;	# remove all trailing whitespace
  $seq =~ s/>//g if defined $seq;
  $seq =~ s/\s+//g if defined $seq; # remove all whitespace, including newlines

  return($hdr, $seq);
}

# Closes the file and undefs the database object.
sub close {
  my $self = shift;
  my $fh = $self->{'fh'};
  my $filename = $self->{'filename'};
  close($fh) or carp("Warning: Could not close $filename\: $!\n");
  undef($self);
}

# Destructor. This is called when you undef() an object
sub DESTROY {
  my $self = shift;
  $self->close;
}
