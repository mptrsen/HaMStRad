#!/usr/bin/env perl
use lib '../bin/';
use Filehandler;
use strict;
use Getopt::Long;

# PROGRAM NAME: NENTFERNER.PL

# AUTHOR: INGO EBERSBERGER, ebersber@eva.mpg.de
# DATE: 
# DESCRIPTION: THE PROGRAM ALL REMOVES NEWLINES FROM THE INPUT TEXT

# DATE LAST MODIFIED: 22/01/2001; 12.02.2004

#################### START MAIN PROGRAM #####################
my $p;
## parse the command line
my $infile;
my $space = '';
my $sep = '';
my $path;
my $ls;
my $rfq;
my $outname;
my $help;

my $message = "
options:
-in=<path/infile>\n
-out=<outfile>\n
-space=<yes||no>\tdefault:no\n
-sep=<tab||newline||no>\tdefault:newline. Set value to 'no' if you just want to join the lines\n
-leading_space: set this flag to remove one blank from the beginning of each line in the file\n
-reformat_qual: set this flag to reformat phred qual files such that every qual value occupies two bytes\n
";
if (@ARGV==0) {
	$help =1;
}

if (-e "nentferner.out") {
    `rm -f nentferner.out`;
}
GetOptions ("h" => \$help,
	    "in=s" => \$infile, 
	    "out=s"=> \$outname,
	    "space=s" => \$space,
	    "sep=s" => \$sep,
	    "leading_space" => \$ls,
	    "reformat_qual" => \$rfq);
 
###############
## help
if (defined $help) {
    die "$message\n";
}

my $seq;
my $crunch = 0;
my @head;
my $intcount;
my $join;
## interpretation of command line values
if ($infile =~ /\//) {
    ($path, $infile) = $infile =~ /(.*)\/(.*)?$/;
}
if (!(defined $path)) {
    $path = '.';
}
#
if ($space eq 'yes') {
	$space = ' ';
}
else {
	$space = '';
	print "\nnewlines will be removed without replacement\n";
}
#
if ($sep eq 'tab') {
    $join = "\t";
    $sep = "\t";
}
elsif ($sep eq 'no') {
    $join = '';
    $sep = '';
}
else {
    $join = ' ';
    $sep = "\n";
}
#i####### if reformat_qual is set, automatically set $sep to ' ' and $ls to 1
if ($rfq) {
	$space = ' ';
	$ls = 1;
}
if (!(defined $outname)) {
    $outname = 'nentferner.out';
}

die "The file $path/$infile does not exist!\n" unless -e "$path/$infile";
tie (*IN, Filehandler::, "$infile", "$path", "\n");
tie (*OUT, Filehandler::, "$outname", ">$path", "\n");
while (readline (IN)) {
    if ($_ =~ />/) {
	## the fasta header
	if ($crunch == 1) {
	    if (defined $rfq) {
            $seq =~ s/(\s{1}\d{1}\s{1})/ $1/g;
            $seq =~ s/(\s{1}\d{1})(\s{1}\d{1}\s{1})/$1 $2/g;
	    $seq =~ s/^(\d{1}\s{1})/ $1/;
	    $seq =~ s/(\s{1}\d{1})$/ $1/;
	    }
	    my $outline = (join "$join", @head) . "$sep" . $seq . "\n";
	    print OUT $outline;
	}
	chomp $_;
	my $head = $_;
#	$head =~ s/>.*\|//;
	@head = split /\s{1,}/, $head;
	$seq = '';
	$crunch = 1;
	$intcount ++;
	if ($intcount%1000 == 0) {
	    print "$intcount lines processed\n";
	}
    }
    else {
	if (defined $ls) {
		$_ =~ s/^\s{1}//;
	}
	$_ =~ s/\s{1,}$/$space/g;
	$seq .= $_;
    }
}
if (defined $rfq) {
            $seq =~ s/(\s{1}\d{1}\s{1})/ $1/g;
            $seq =~ s/(\s{1}\d{1})(\s{1}\d{1}\s{1})/$1 $2/g;
            $seq =~ s/^(\d{1}\s{1})/ $1/;
            $seq =~ s/(\s{1}\d{1})$/ $1/;
            }

my $outline = (join "$join", @head) . "$sep" . $seq . "\n";
print OUT $outline;
close (OUT);
close (IN);
