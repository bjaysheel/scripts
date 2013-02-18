#!/usr/bin/perl

=head1 NAME
   synthetic_fastq_update.pl

=head1 SYNOPSIS

    USAGE: synthetic_fastq_update.pl --input synthetic fastq file --output output dir --qual qual file | --real fastq

=head1 OPTIONS

B<--input, -i>
   Input file synthetic fastq file

B<--output, -o>
	Output directory

B<--qual, -q>
   quality file

B<--real, -r>
	real fastq file

B<--help,-h>
   This help message

=head1  DESCRIPTION
	Update quality value in synthetic fastq with data in qual file or from
	some other fastq file

	if both real and qual files are provided only qual file will be used.

=head1  INPUT
	synthetic fastq file is required,
	quality file or real fastq files is required

=head1  OUTPUT
	synthetic fastq file with updated quality values

	if real fastq value is provided a qual file is created that can be used
	again later.

=head1  CONTACT
  Jaysheel D. Bhavsar @ bjaysheel[at]gmail[dot]com

==head1 EXAMPLE
   synthentic_fastq_update.pl --input synthetic_file --qual qual.file

=cut

use strict;
use warnings;
use Data::Dumper;
use File::Basename;
use Pod::Usage;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);

my %options = ();
my $results = GetOptions (\%options,
                          'input|i=s',
						  'output|o=s',
						  'qual|q=s',
						  'real|r=s',
						  'help|h') || pod2usage();

## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}
#############################################################################
## make sure everything passed was peachy
&check_parameters(\%options);

open (FHD, "<", $options{input}) or die "Could not open file $!\n";

#get all qual values in a an array
open (QUL, "<", $options{qual}) or die "Could not open file $!\n";
my @qual = <QUL>;
close (QUL);

my $update = $options{output}."/".basename($options{input})."_updated.fastq";
open (OUT, ">", $update) or die "Coulld not open file $!\n";

my $lCounter = 0;
my $qCounter = 0;

while (<FHD>) {
	if ($lCounter != 3){
		#output first 3 lines of input file as is
		print OUT $_;
		$lCounter++;
	} else {
		#rest qual counter if qual files is smaller than input file.
		if ($qCounter >= scalar(@qual)){
			$qCounter = 0;
		}

		#print qual data.
		print OUT $qual[$qCounter];

		$qCounter++;
		$lCounter = 0;
	}
}

close(OUT);
close(FHD);

#############################################################################
sub check_parameters {
    my $options = shift;

	my @required = ("input", "output");

	foreach my $key (@required) {
		unless ($options{$key}) {
			print STDERR "ERROR: ARG $key is required\n";
			pod2usage({-exitval => 2,  -message => "error message", -verbose => 1, -output => \*STDERR});
			exit(-1);
		}
	}

	#error if output dir does not exist
	unless (-d $options{output}){
		print STDERR "ERROR: not a directory $options{output}\n";
		pod2usage({-exitval => 2,  -message => "error message", -verbose => 1, -output => \*STDERR});
		exit(-1);
	}

	#check if qual file or real fastq file exist
	if (!defined $options{qual} && !defined $options{real}){
		print STDERR "ERROR: real fastq or qual file is required\n";
		pod2usage({-exitval => 2,  -message => "error message", -verbose => 1, -output => \*STDERR});
		exit(-1);
	}

	#warn if both file passed
	if (defined $options{qual} && defined $options{real}){
		print STDERR "WARNING: qual file is used\n";
	}

	#convert real into qual file
	if (defined $options{real}){
		createQualFile($options{real});
	}
}

sub createQualFile {
	my $file = shift;

	$options{qual} = $options{output}."/".basename($file).".qual";

	#output very 4th line of the real fastq file.
	system("sed -n '0~4p' $file > $options{qual}");
}
