#!/usr/local/biotools/perl/5.14.2/bin/perl

=head1 NAME
   updateRG.pl

=head1 SYNOPSIS

    USAGE: updateRG.pl --input_dir

=head1 OPTIONS

B<--input_dir,-i>
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
   updateRG.pl

=cut

use strict;
use warnings;
use Data::Dumper;
use Pod::Usage;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);

my %options = ();
my $results = GetOptions (\%options,
                          'input_dir|i=s',
						  'help|h') || pod2usage();

## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}
#############################################################################
## make sure everything passed was peachy
&check_parameters(\%options);

opendir (DIR, $options{input_dir}) or die "Could not open file $!\n";

while(my $file = readdir(DIR)) {
	next if ($file =~ /^\./);
	next if ($file =~ /\.bai$/);

	my $id = $file;
	$id =~ s/\.igv-sorted\.bam$//;
	$id =~ s/^s_//;
	$id =~ s/_|-//g;

	my $cmd = "qsub -q 7-days -l h_vmem=2G -l h_stack=10M -m ae -M bhavsar.jaysheel\@mayo.edu -b y -N UPDATERG.$id";
	$cmd .= " /usr/java/latest/bin/java -Xmx1g -jar /projects/bsi/bictools/apps/alignment/picard/latest/AddOrReplaceReadGroups.jar INPUT=$options{input_dir}/$file";
	$cmd .= " OUTPUT=$options{input_dir}/$file.tmp SORT_ORDER=coordinate RGLB=$options{RGLB} RGPL=$options{RGPL} RGCN=$options{RGCN}";
	$cmd .= " RGID=$id RGSM=$id RGPU=$id";

	system("$cmd");
}

close(DIR);
exit(0);

#############################################################################
sub check_parameters {
    my $options = shift;

	my @required = ("input_dir");

	foreach my $key (@required) {
		unless ($options{$key}) {
			print STDERR "ARG: $key is required\n";
			pod2usage({-exitval => 2,  -message => "error message", -verbose => 1, -output => \*STDERR});
			exit(-1);
		}
	}

	unless($options{RGLB}) { $options{RGLB} = "hg19"; }
	unless($options{RGPL}) { $options{RGPL} = "ILLUMINA"; }
	unless($options{RGCN}) { $options{RGCN} = "MAYO"; }
}
