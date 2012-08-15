#!/usr/local/biotools/perl/5.14.2/bin/perl

=head1 NAME
   create_sample_info.pl

=head1 SYNOPSIS

    USAGE: create_sample_info.pl --input_dir=result.txt [--token=BAM]

=head1 OPTIONS

B<--input_dir, -i>
   input directory

B<--token, -t>
   init token for same file.

B<--help,-h>
   This help message

=head1  DESCRIPTION

=head1  INPUT


=head1  OUTPUT

=head1  CONTACT
  Jaysheel D. Bhavsar @ bjaysheel[at]gmail[dot]com


==head1 EXAMPLE


=cut

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);

my %options = ();
my $results = GetOptions (\%options,
                          'input_dir|i=s',
						  'output_file|o=s',
						  'prefix|p=s',
						  'regex|r=s',
						  'help|h') || pod2usage();

## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}
#############################################################################
&check_parameters(\%options);

my @dirs = ("/data2/delivery/Beutler_Andreas_m068039/120321_SN730_0147_BD0UH5ACXX/secondary/IGV_BAM"
		   ,"/data2/delivery/Beutler_Andreas_m068039/120322_SN7001166_0060_AC0LA9ACXX/secondary/IGV_BAM"
		   ,"/data2/delivery/Beutler_Andreas_m068039/120413_SN7001166_0061_AC0MEPACXX/secondary/IGV_BAM"
		   ,"/data2/delivery/Beutler_Andreas_m068039/120413_SN7001166_0062_BC0MELACXX/secondary/IGV_BAM");

open(OUT, ">", $options{output_file}) or die "Could not open file $options{output_file}\n";

my $samplenames = "";
my $laneindex = "";
my $lanenum = "";

foreach my $dir (@dirs){
	opendir DIR, $dir or die "Could not open dir $dir\n";
	my @files = readdir DIR;
	close DIR;

	my @data = split(/\//, $dir);
	my $flowcell = $data[-3];

	my $reports = "/data2/bsi/reports/$flowcell/fastqc";

	foreach my $file (@files){
		next if ($file =~ m/^\.|^\.\./);
		next if ($file !~ m/$options{regex}$/);

		## create symbolic link into input_dir;
		#ln -s $dir/$file $options{input_dir}/$file

		my $sample = $file;
		$sample =~ s/$options{regex}$//;
		my ($lane, $index) = getLaneInfo($sample, $reports);

		$sample =~ s/_|-//g;
		print OUT $options{prefix}.":".$sample."=".$file."\n";

		$samplenames .= $sample.":";
		$laneindex .= "$index:";
		$lanenum .= "$lane:";
	}
}

$samplenames =~ s/:$//;
$laneindex =~ s/:$//;
$lanenum =~ s/:$//;

print $samplenames."\n\n";
print $laneindex."\n\n";
print $lanenum."\n\n";

close OUT;
exit(0);

#############################################################################
sub check_parameters {
    my $options = shift;

	my @required = ("input_dir", "output_file");

	foreach my $key (@required) {
		unless ($options{$key}) {
			print STDERR "ARG: $key is required\n";
			pod2usage({-exitval => 2,  -message => "error message", -verbose => 1, -output => \*STDERR});
			exit(-1);
		}
	}

	$options{regex} = "fastq\.gz|\.bam" unless (defined $options{regex});
}

#############################################################################
sub getLaneInfo {
	my ($sample, $reports) = @_;
	$sample = substr($sample, 2);

	opendir RPT, $reports or die "Could not open dir $reports\n";
	while ((my $filename = readdir(RPT))) {
		next if ($filename !~ /$sample.*_L(\d)_R\d_I(\w+)_fastqc.zip/);

		return ($1, $2);
	}
}
