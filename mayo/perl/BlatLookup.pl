#!/usr/local/biotools/perl/5.14.2/bin/perl

=head1 NAME
   BlatLookup.pl

=head1 SYNOPSIS

    USAGE: BlatLookup.pl --input input_bed_file --output_dir output_dir --ref reference [--window 10]

=head1 OPTIONS

B<--input,-i>
   Input file

B<--output_dir,-o>
   output dir

B<--window, -w>
    Window size

B<--help,-h>
   This help message

=head1  DESCRIPTION


=head1  INPUT


=head1  OUTPUT


=head1  CONTACT
  Jaysheel D. Bhavsar @ bjaysheel[at]gmail[dot]com


==head1 EXAMPLE
   BlatLookup.pl --input=OnTarget.bam --window=100 --output=OUTFILE --ref reference_file

=cut

use strict;
use warnings;
use Data::Dumper;
use Bio::SeqIO;
use Pod::Usage;
use Scalar::Util qw(looks_like_number);
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);

my %options = ();
my $results = GetOptions (\%options,
                          'input|i=s',
						  'output_dir|o=s',
                          'window|w=s',
						  'ref|r=s',
						  'help|h') || pod2usage();

## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}
#############################################################################
## make sure everything passed was peachy
&check_parameters(\%options);

open (FHD, "<", $options{input}) or die "Cound not open file $options{input}\n";

## create index
#my $ref_idx = Bio::Index::Fasta->new( '-filename' => $options{output}.'.faa.idx',
#									  '-write_flag' => 1 );
#$ref_idx->make_index($options{ref});

my $ref = Bio::SeqIO->new(-format => 'fasta',
                          -file   => $options{ref});
my $ref_seq = Bio::Seq->new(-id => "NA", -seq=>"ACTG");

my $prev = "";
while (<FHD>) {
	chomp $_;
	my @data = split(/\t/, $_);

	if (($prev =~ /$data[0]/) && ($ref_seq->display_id =~ /NA/)) {
		next;
	}

	unless ($ref_seq->display_id =~ /$data[0]$/) {
		$ref_seq = getRefSeq($data[0]);
		print "opening new file for $data[0]\n";
		open (FA, ">", $options{output_dir}. "/lookup_" .$data[0]. ".fa") or die "Cound not open file $options{output_dir}/lookup_$data[0].fa\n";
	}

	if ($ref_seq =~ /NA/) {
		print "WARNING: Cannot find $data[0] in reference file skipping all instance of $data[0]\n";
		$prev = $data[0];
		next;
	}

	foreach my $loc ($data[1]..$data[2]) {
		print FA ">".$data[0]."_".$loc."\n";
		print FA $ref_seq->subseq($loc-($options{window}/2), $loc+($options{window}/2))."\n";
	}
}

close(FHD);
close(FA);

exit(0);

#############################################################################
sub check_parameters {
    my $options = shift;

	my @required = qw(input output_dir ref);

	foreach my $key (@required) {
		unless ($options{$key}) {
			print STDERR "ARG: $key is required\n";
			pod2usage({-exitval => 2,  -message => "error message", -verbose => 1, -output => \*STDERR});
			exit(-1);
		}
	}

	unless($options{window}) { $options{window} = 100; }

	## ontarget: /data4/bsi/refdata/exomeCapture/AgilentV4_ucsc_refflat_hg19_2011-01-24.10bp.merge.bed
	## ref: /data2/bsi/reference/sequence/human/ncbi/37.1/allchr.fa
}

#############################################################################
sub getRefSeq {
	my $seqId = shift;

	while (my $seq = $ref->next_seq) {
		if ($seq->display_id =~ /$seqId$/) {
			return $seq;
		}
	}

	return "NA";
}
