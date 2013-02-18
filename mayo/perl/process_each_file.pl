#!/usr/local/biotools/perl/5.14.2/bin/perl

=head1 NAME
   process_each_file.pl

=head1 SYNOPSIS

    USAGE: process_each_file.pl --input input_dir

=head1 OPTIONS

B<--input,-i>
   Input dir

B<--help,-h>
   This help message

=head1  DESCRIPTION


=head1  INPUT


=head1  OUTPUT


=head1  CONTACT
  Jaysheel D. Bhavsar @ bjaysheel[at]gmail[dot]com


==head1 EXAMPLE
   process_each_file.pl

=cut

use strict;
use warnings;
use Data::Dumper;
use DBI;
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

opendir (DIR, $options{input}) or die "Could not open dir $!\n";

while(my $file = readdir(DIR)) {
	next if ($file =~ /^\./);

	my $cmd = "qsub -V -q 7-days -l h_vmem=4G -l h_stack=10M -m a -M bhavsar.jaysheel\@mayo.edu -b y -N INDEXBAM.$file";
	$cmd .= " samtools index $options{input}/$file";

	print $cmd."\n";
	system($cmd);
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
