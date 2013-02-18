#!/usr/local/biotools/perl/5.14.2/bin/perl

=head1 NAME
   tophatFusionMerge.pl

=head1 SYNOPSIS
    USAGE: tophatFusionMerge.pl i=input sam file -o=output file name

=head1 OPTIONS


B<--input, -i>
	Input sam file

B<--output, -o>
	Output sam file name

B<--help,-h>


=head1  DESCRIPTION
	Get complete read from tophat fusion candidates

=head1  INPUT
    Tophat alignmed Sam file filted with only fusion candidates

=head1  OUTPUT
    A sam file with all filter candidates stiched togeter.

=head1  CONTACT



==head1 EXAMPLE
	./tophatFusionMerge.pl -i=input.bam -o=output.bam

=cut

use strict;
use warnings;
use Data::Dumper;
use Pod::Usage;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);

my %options = ();
my $results = GetOptions (\%options,
                          'input|i=s',
						  'output|o=s',
						  'log|l=s',
			              'debug=s',
						  'help|h') || pod2usage();

## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

#############################################################################
## make sure everything passed was peachy
&check_parameters(\%options);

open(FHD, "<", $options{input}) or die("Could not open $options{input}");
open(OUT, ">", $options{output}) or die("Could not write file $options{output}");

my $cline=0;
while (<FHD>){
    chomp $_;

    if ($cline==0) {
        my @fields = split(/\t/, $_);
        my $i=0;
        my @nfields;

        for (my $i=0; $i <10; $i++) {
            $nfields[$i] = $fields[$i];
        }

        my $XF_field;
        foreach my $field (@fields) {
            if ($field =~ /^XF/) {
                $XF_field = $field;
                last;
            }
        }

        my @xfields = split (/\s/, $XF_field);
        my ($fchr, $schr) = split (/-/,$xfields[1]);

        (my $cigar = $xfields[3]) =~ s/F/N/;

        $nfields[2]=$fchr;
        $nfields[3]=$xfields[2];
        $nfields[5]=uc $cigar;
        $nfields[9]=$xfields[4];
        $nfields[10]=$xfields[5];

        foreach my $field (@nfields) {
            print OUT $field."\t";
        }
        print OUT "\n";
    }

    $cline = ($cline==1)?0:1;
}
close(FHD);
close(OUT);

exit();

#############################################################################
sub check_parameters {
    my $options = shift;

	my @required = qw(input output);

	foreach my $key (@required) {
		unless ($options{$key}) {
			print STDERR "ARG: $key is required\n";
			pod2usage({-exitval => 2,  -message => "error message", -verbose => 1, -output => \*STDERR});
			exit(-1);
		}
	}
}
