#!/usr/bin/perl

# check_for_duplicates.pl
# Checks for and removes multiple assignment of sequences in HaMStR-AD output
# VERSION 2.0

use strict;
use warnings;
use File::Spec;
use Carp;
use File::Copy qw(copy);


## Read command line
&usage if @ARGV % 2;
my %options = @ARGV;
&check_options(%options);

# Get names of files
my @HAMSTR_RESULT_FOLDERS = read_dir( $options{'-d'}, {'hide' => 1} );

foreach my $HAMSTR_RESULT_FOLDER ( @HAMSTR_RESULT_FOLDERS ) {

    print "Processing: ", $HAMSTR_RESULT_FOLDER , "\n";

    # Get names of all files in aa directory
    my $aa_path = File::Spec->catdir( $options{'-d'}, $HAMSTR_RESULT_FOLDER, 'aa' );

    my @filenames = read_dir( $aa_path, {'hide' => 1} );

    my %presence_of;
    foreach my $filename ( @filenames ) {
        my $file = File::Spec->catfile( $aa_path, $filename );

        open( my $input_fh, '<', $file) or die "Cannot open file \"$file\": $!\n";

        while (my $line = <$input_fh>) {
            next if $line !~ m/^>/;
            chomp $line;
            my @header_parts = split('\|', $line);
            next if @header_parts < 4;
            my @hits = split ('PP', $header_parts[-1]);
            foreach my $hit (@hits) {
                $hit =~ m/(.*)-\d+$/;
                my $sequence_id = $1;
                push ( @{$presence_of{$1}}, $filename );
            }
        }
        close $input_fh;
    }

    # Print results into log file
    my $logfile = $options{'-l'}.$HAMSTR_RESULT_FOLDER.'_log.txt';
    open( my $output_fh, '>', $logfile) or die "Cannot open file \"$logfile\": $!\n";

    my %is_file_affected_from_redundancy;

    foreach my $sequence_id ( keys %presence_of ) {
#	print $sequence_id, "\n";
        next if @{$presence_of{$sequence_id}} < 2;
        print {$output_fh} $sequence_id;
        foreach my $file ( @{$presence_of{$sequence_id}} ) {
            print {$output_fh} "\t", $file;
            $file =~ m/^(.*)\.aa.fa$/;
            $is_file_affected_from_redundancy{$1} = 1; # record name of problematic ortholog group
        }
        print {$output_fh} "\n";
    }
    close $output_fh;

    # Copy file with redundant and non-redundant contigs to separate folders

    # Create new folders
    my $aa_nr_path = File::Spec->catdir( $options{'-d'}, $HAMSTR_RESULT_FOLDER, 'aa_with_no_redundancy' );
    my $aa_ro_path = File::Spec->catdir( $options{'-d'}, $HAMSTR_RESULT_FOLDER, 'aa_with_redundancy_only' );
    my $nt_nr_path = File::Spec->catdir( $options{'-d'}, $HAMSTR_RESULT_FOLDER, 'nt_with_no_redundancy' );
    my $nt_ro_path = File::Spec->catdir( $options{'-d'}, $HAMSTR_RESULT_FOLDER, 'nt_with_redundancy_only' );

    mkdir $aa_nr_path;
    mkdir $aa_ro_path;
    mkdir $nt_nr_path;
    mkdir $nt_ro_path;

    # Copy aa files

    @filenames = read_dir( $aa_path, {'hide' => 1} );

    foreach my $filename ( @filenames ) {

        my $source_file = File::Spec->catfile( $aa_path, $filename );
        my $target_file;
        $filename =~ m/^(.*)\.aa.fa$/ ;
        if ( exists $is_file_affected_from_redundancy{$1} ) {
            $target_file = File::Spec->catfile( $aa_ro_path, $filename );
        }
        else {
            $target_file = File::Spec->catfile( $aa_nr_path, $filename );
        }
        copy $source_file, $target_file;
    }

    # Copy nt files
    my $nt_path = File::Spec->catdir( $options{'-d'}, $HAMSTR_RESULT_FOLDER, 'nt' );
    @filenames = read_dir( $nt_path, {'hide' => 1} );

    foreach my $filename ( @filenames ) {

        my $source_file = File::Spec->catfile( $nt_path, $filename );
        my $target_file;
        $filename =~ m/^(.*)\.nt.fa$/ ;
        if ( exists $is_file_affected_from_redundancy{$1} ) {
            $target_file = File::Spec->catfile( $nt_ro_path, $filename );
        }
        else {
            $target_file = File::Spec->catfile( $nt_nr_path, $filename );
        }
        copy $source_file, $target_file;
    }


}


## Subroutines

# Help for command line input
sub usage {
  print "Unknown option: @_\n" if ( @_ );
  print "usage: perl check_for_and_remove_duplicates.pl [-d PATH_TO_DIRECTORY_WITH_(ONLY!)_HAMSTR_OUTPUT_FOLDERS]\n";
  print "                                    [-l PATH_WHERE_TO_STORE_THE_LOG_FILES]\n";

  print "\nProgram reports sequence ids that HAMSTR assigned to different orthologs.\n";
  print "\nIt also moves all hit files that are not affected by redundency to the folders\n";
  print "\'aa_without_redundany\' and \'nt_without_redundany\', respectively, and it moves\n";
  print "all hits that are affected by redundany to the folders \'aa_with_redundany_only\' \n";
  print "and \'nt_with_redundany_only\', respectively\n";

  exit;
}

sub check_options {
    my (%options) = @_;

    &usage() if keys %options < 2 or keys %options > 2;

    my @required = ('-d', '-l' );
    foreach my $option ( @required ) {
        &usage() if !exists  $options{$option};
        &usage() if !defined $options{$option};
        delete $options{$option}
    }
    &usage() if keys %options
}

sub read_dir {

	# Unpack @_
	my ( $path, $arg_ref ) = @_;

	# Defaults
	my %DEFAULT_OF = ( 'hide' => 1,  # 1: hide system files; 0: don't
	);

	# Check provided arguments
	croak 'Missing or superfluous aruments'         if @_ < 1  || @_ > 2;
	croak 'Option(s) not passed via anonymous hash' if @_ == 2 && ref $arg_ref ne 'HASH';

	foreach my $provided_options ( keys %{ $arg_ref } ) {
		croak 'Unknown option(s)' if !exists $DEFAULT_OF{ $provided_options };
	}

	# Set defaults
	#          If option given...            Use option             Else default
	my $hide = exists $arg_ref->{'hide'}  ?  $arg_ref->{'hide'}  :  $DEFAULT_OF{'hide'};

	# Open directory handle
	opendir ( my $dir, $path ) or
		croak "Couldn't find path \"$path\": $!";

	# Read file names
	my @files = readdir( $dir ) or
		croak "Couldn't read directory \"$path\": $!";

	# Close directory handle
	closedir ( $dir ) or
		croak "Couldn't close directory \"$path\": $!";

	# Filter hidden system files out
	if ( $hide ) {
		@files = grep {! /^\./ } @files
	}

	# Return file names
	return @files;

}
