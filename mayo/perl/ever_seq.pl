#!/usr/bin/perl

=head1 NAME
   even_seq.pl

=head1 SYNOPSIS

    USAGE: ever_seq.pl --input input_txt file --output_dir

=head1 OPTIONS

B<--input,-i>
   A tab delimited input file with sameple name and base dir where sorted bam
   files for the same are located.

	e.g:
		6027_S  /data2/delivery/Poland_Gregory_gap01/120221_SN414_0174_BD0RLFACXX/secondary/IGV_BAM
		6027_U  /data2/delivery/Poland_Gregory_gap01/120221_SN414_0174_BD0RLFACXX/secondary/IGV_BAM
		6132_S  /data2/delivery/Poland_Gregory_gap01/120221_SN414_0174_BD0RLFACXX/secondary/IGV_BAM
		6132_U  /data2/delivery/Poland_Gregory_gap01/120221_SN414_0174_BD0RLFACXX/secondary/IGV_BAM
		6093_S  /data2/delivery/Poland_Gregory_gap01/120221_SN414_0174_BD0RLFACXX/secondary/IGV_BAM
		6093_U  /data2/delivery/Poland_Gregory_gap01/120221_SN414_0174_BD0RLFACXX/secondary/IGV_BAM
		6223_S  /data2/delivery/Poland_Gregory_gap01/120221_SN414_0174_BD0RLFACXX/secondary/IGV_BAM
		6223_U  /data2/delivery/Poland_Gregory_gap01/120221_SN414_0174_BD0RLFACXX/secondary/IGV_BAM

B<--output_dir, -d>
	Output directory location.

B<--help,-h>
   This help message

=head1  DESCRIPTION

=head1  INPUT


=head1  OUTPUT


=head1  CONTACT
  Jaysheel D. Bhavsar @ bjaysheel[at]gmail[dot]com


==head1 EXAMPLE
   ever_seq.pl --input filename.txt

=cut

use strict;
use warnings;
use Data::Dumper;
use Cwd;
use File::Basename;
use Pod::Usage;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);

my %options = ();
my $results = GetOptions (\%options,
                          'input|i=s',
						  'output_dir|d=s',
						  'email|e=s',
						  'queue|q=s',
						  'help|h') || pod2usage();

## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

#############################################################################
## make sure everything passed was peachy
&check_parameters(\%options);

open (FHD, "<", $options{input}) or die "Cound not open file $options{input}\n";
my $tmp_sh = $options{output_dir}."/ever_seq_cmd.sh";

while (<FHD>){
	chomp $_;

	my @data = split(/\t/,$_);

	if (defined $ENV{PAYTHONPATH}){
		$ENV{PYTHONPATH}="/data4/bsi/investigator_projects/couch/Progams/EVER-seq-1.0.7/usr/local/biotools/python/2.7/lib/python2.7/site-packages:".$ENV{PYTHONPATH};
	} else {
		$ENV{PYTHONPATH}="/data4/bsi/investigator_projects/couch/Progams/EVER-seq-1.0.7/usr/local/biotools/python/2.7/lib/python2.7/site-packages";
	}

	$ENV{PATH} = "/data4/bsi/investigator_projects/couch/Progams/EVER-seq-1.0.7/usr/local/biotools/python/2.7/bin:".$ENV{PATH};

	my $ref_Bed = "/data2/bsi/RandD/Synthetic_Genomes/Tool_Tester/RNASEQ/ref_exons.bed"; # Replace this with your ref genome in bed format

	my $filename = $data[1]."/s_".$data[0]."-sorted.bam";
	my $output_prefix = $options{output_dir}."/".$data[0]."/".$data[0];
	my $tmp_sam = $options{output_dir}."/".$data[0].".sam";

	#create output dir for each sample
	system("mkdir -p $options{output_dir}/$data[0]");

	#wirte new command for each input.
	open (OUT, ">", $tmp_sh) or die "Counld not write to temp file $tmp_sh\n$!\n";
	#print OUT "echo \"Creating sam file $tmp_sam\" \n";
	#print OUT "/projects/bsi/bictools/apps/alignment/samtools/samtools-0.1.18/samtools view $filename > $tmp_sam \n";
	#print OUT "\n";
	print OUT "echo \"Run RPKM saturation\" \n";
	print OUT "python /data4/bsi/investigator_projects/couch/Progams/EVER-seq-1.0.7/scripts/RPKM_saturation.py ";
	print OUT "-i $tmp_sam -o $output_prefix -r $ref_Bed -l 10 -u 100 -s 10 \n";
	#print OUT "echo \"Removing temporary SAM file $tmp_sam\" \n";
	#print OUT "rm $tmp_sam\n";
	close(OUT);

	system("qsub -V -cwd -m sae -M $options{email} -l h_vmem=150G -l h_stack=200M -q $options{queue} $tmp_sh");
}
#remove tmp file.
`rm $tmp_sh`;

close(FHD);

#############################################################################
sub check_parameters {
    my $options = shift;

	my @required = ("input", "output_dir");

	foreach my $key (@required) {
		unless ($options{$key}) {
			print STDERR "ARG: $key is required\n";
			pod2usage({-exitval => 2,  -message => "error message", -verbose => 1, -output => \*STDERR});
			exit(-1);
		}
	}

	if (-f $options{output_dir}){
		print STDERR "File passed where Directory expected\n";
		pod2usage({-exitval => 2,  -message => "error message", -verbose => 1, -output => \*STDERR});
		exit(-1);
	}

	unless (-d $options{output_dir}){
		mkdir $options{output_dir} or die "Could not create output directory $options{output_dir}\n$!\n";
		pod2usage({-exitval => 2,  -message => "error message", -verbose => 1, -output => \*STDERR});
		exit(-1);
	}

	unless (defined $options{email}) { $options{email} = "bhavsar.jaysheel\@mayo.edu"; }
	unless (defined $options{queue}) { $options{queue} = "lg-mem"; }
}
