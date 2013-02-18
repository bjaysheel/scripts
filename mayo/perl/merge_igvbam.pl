#!/usr/local/biotools/perl/5.14.2/bin/perl

=head1 NAME
   merge_igvbam.pl

=head1 SYNOPSIS

    USAGE: merge_igvbam.pl --sample_info input_sample_info [--dir1....  --dir2 ....]

=head1 OPTIONS

B<--input,-i>
   Input file

B<--output,-o>
   output prefix

B<--help,-h>
   This help message

=head1  DESCRIPTION


=head1  INPUT


=head1  OUTPUT


=head1  CONTACT
  Jaysheel D. Bhavsar @ bjaysheel[at]gmail[dot]com


==head1 EXAMPLE
   merge_igvbam.pl

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
                          'sample_info|s=s',
						  'output_dir|o=s',
						  'help|h') || pod2usage();

## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}
#############################################################################
## make sure everything passed was peachy
&check_parameters(\%options);

open (DAT, "<", $options{sample_info}) or die "Could not open file $!\n";

while(<DAT>) {
	chomp $_;

	$_ =~ s/^(BAM)//;

	my @data = split(/=|\t/, $_);

	$data[0] =~ s/://;

	my $cmd = "qsub -V -wd $options{output_dir}/logs -q 7-days -l h_vmem=6G -l h_stack=10M -m ae -M bhavsar.jaysheel\@mayo.edu -b y -N MERGEBAM.$data[0]";
	$cmd .= " /usr/java/latest/bin/java -Xmx4g -jar /projects/bsi/bictools/apps/alignment/picard/latest/MergeSamFiles.jar";
	$cmd .= " SORT_ORDER=coordinate OUTPUT=$options{output_dir}/IGV_BAM/$data[0].igv-sorted.bam";

	if ($data[1] =~ /^s_/) {
		$cmd .= " INPUT=". $options{dir1} ."/".$data[1];
		$cmd .= " INPUT=". $options{dir2} ."/".$data[2];
	} else {
		$cmd .= " INPUT=". $options{dir1} ."/".$data[2];
		$cmd .= " INPUT=". $options{dir2} ."/".$data[1];
	}

	print $cmd."\n";
	system("$cmd");
}

exit(0);

#############################################################################
sub check_parameters {
    my $options = shift;

	my @required = qw(sample_info output_dir);

	foreach my $key (@required) {
		unless ($options{$key}) {
			print STDERR "ARG: $key is required\n";
			pod2usage({-exitval => 2,  -message => "error message", -verbose => 1, -output => \*STDERR});
			exit(-1);
		}
	}

	unless($options{dir1}) { $options{dir1} = "/data2/bsi/secondary/Beutler_Andreas_m068039/exome/120322_SN7001166_0060_AC0LA9ACXX/IGV_BAM"; }
	unless($options{dir2}) { $options{dir2} = "/data2/bsi/secondary/Beutler_Andreas_m068039/exome/120726_SN616_0199_BC117YACXX.IGV_BAM/IGV_BAM"; }
}
