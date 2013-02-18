#!/usr/local/biotools/perl/5.14.2/bin/perl

=head1 NAME
   sort_bam_header.pl

=head1 SYNOPSIS

    USAGE: sort_bam_header.pl --input input header txt

=head1 OPTIONS

B<--input,-i>
   Input file

B<--help,-h>
   This help message

=head1  DESCRIPTION


=head1  INPUT


=head1  OUTPUT


=head1  CONTACT
  Jaysheel D. Bhavsar @ bjaysheel[at]gmail[dot]com


==head1 EXAMPLE
   sort_bam_header.pl

=cut

use strict;
use warnings;
use Data::Dumper;
use Pod::Usage;
use Scalar::Util qw(looks_like_number);
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);

my %options = ();
my $results = GetOptions (\%options,
                          'input|i=s',
						  'help|h') || pod2usage();

## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}
#############################################################################
## make sure everything passed was peachy
&check_parameters(\%options);

open (FHD, $options{input}) or die "Could not open file $!\n";
my @array=();

while(<FHD>) {
	chomp $_;

	my @data = split(/\t/, $_);

	my $number = $data[1];
	$number =~ s/SN:chr//g;
	print $number."\n";

	push @array, {"SQ" => $data[0],
				  "chr" => $data[1],
				  "ln" => $data[2],
				  "number" => $number};
}

my @sorted = sort { lc($a->{number}) cmp lc($b->{number}); } @array;

#print ("SN:chrGL000191.1" cmp "SN:chrGL000226.1")."\n";

foreach my $hash (@sorted){
	print $hash->{SQ}."\t".$hash->{chr}."\t".$hash->{ln}."\n";
}

exit(0);

#############################################################################
sub check_parameters {
    my $options = shift;

	my @required = qw(input);

	foreach my $key (@required) {
		unless ($options{$key}) {
			print STDERR "ARG: $key is required\n";
			pod2usage({-exitval => 2,  -message => "error message", -verbose => 1, -output => \*STDERR});
			exit(-1);
		}
	}
}
