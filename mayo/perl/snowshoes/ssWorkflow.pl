#!/usr/local/biotools/perl/5.14.2/bin/perl

=head1 NAME
   ssWorkflow.pl

=head1 SYNOPSIS

    USAGE:

=head1 OPTIONS

B<--input,-i>


B<--output, -o>


B<--reference,-r>


B<--window, -w>


B<--blat_path, -b>


B<--samtools_path, -sam>


B<--blat_ref, -br>


B<--minScore, -m>

B<--minIdentity, -t>


B<--help,-h>


=head1  DESCRIPTION


=head1  INPUT


=head1  OUTPUT


=head1  CONTACT
  bjaysheel@gmail.com


==head1 EXAMPLE


=cut

use strict;
use warnings;
use Data::Dumper;
use Cwd;
use Pod::Usage;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use POSIX;

my %options = ();
my $results = GetOptions (\%options,
                          'sample_info|s=s',
						  'config_info|c=s',
						  'output_dir|o=s',
						  'help|h') || pod2usage();

## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}
#############################################################################
## make sure everything passed was peachy
&check_parameters(\%options);


timer(); #call timer to see when process ended.

my $cmd = "";

create_dir_struct();

##link input fasta files to fastq_temp dir.

## Align the two FASTQ files to indexed BWA genome
$cmd = "$options{bwa_directory}/bwa aln -l 32 -t 4 $options{ref_data}/BWA_indexed_genome_build37/allchr_hg19.fa";
$cmd .= " $options{output_dir}/fastq_temp/s_4_1_sequence > $options{output_dir}/BWA_aln_genome/s_4_aln_1.sai";
#system("$cmd");

$cmd = "$options{bwa_directory}/bwa aln -l 32 -t 4 $options{ref_data}/BWA_indexed_genome_build37/allchr_hg19.fa";
$cmd .= " $options{output_dir}/fastq_temp/s_4_2_sequence > $options{output_dir}/BWA_aln_genome/s_4_aln_2.sai";
#system("$cmd");

## Align the two fastq files to indexed BWA junction
$cmd = "$options{bwa_directory}/bwa aln -l 32 -t 4";
$cmd .= " $options{ref_data}/BWA_indexed_Junction_build37/hg19_refFlat_filtered.50.junction";
$cmd .= " $options{output_dir}/fastq_temp/s_4_1_sequence > $options{output_dir}/BWA_aln_Junction/s_4_aln_1.sai";
#system("$cmd");

$cmd = "$options{bwa_directory}/bwa aln -l 32 -t 4";
$cmd .= " $options{ref_data}/BWA_indexed_Junction_build37/hg19_refFlat_filtered.50.junction";
$cmd .= " $options{output_dir}/fastq_temp/s_4_2_sequence > $options{output_dir}/BWA_aln_Junction/s_4_aln_2.sai";
#system("$cmd");

## Create SAM files from the alignment .sai files for the genome alignment
$cmd = "$options{bwa_directory}/bwa samse $options{ref_data}/BWA_indexed_genome_build37/allchr_hg19.fa";
$cmd .= " $options{output_dir}/BWA_aln_genome/s_4_aln_1.sai $options{output_dir}/fastq_temp/s_4_1_sequence >";
$cmd .= " $options{output_dir}/BWA_aln_genome/s_4_1.sam";
#system("$cmd");

$cmd = "$options{bwa_directory}/bwa samse $options{ref_data}/BWA_indexed_genome_build37/allchr_hg19.fa";
$cmd .= " $options{output_dir}/BWA_aln_genome/s_4_aln_2.sai $options{output_dir}/fastq_temp/s_4_2_sequence >";
$cmd .= " $options{output_dir}/BWA_aln_genome/s_4_2.sam";
#system("$cmd");

## Create SAM files from the alignment .sai files for the junction alignment
$cmd = "$options{bwa_directory}/bwa samse $options{ref_data}/BWA_indexed_Junction_build37/";
$cmd .= "hg19_refFlat_filtered.50.junction $options{output_dir}/BWA_aln_Junction/s_4_aln_1.sai";
$cmd .= " $options{output_dir}/fastq_temp/s_4_1_sequence > $options{output_dir}/BWA_aln_Junction/s_4_1.sam";
#system("$cmd");

$cmd = "$options{bwa_directory}/bwa samse $options{ref_data}/BWA_indexed_Junction_build37/";
$cmd .= "hg19_refFlat_filtered.50.junction $options{output_dir}/BWA_aln_Junction/s_4_aln_2.sai";
$cmd .= " $options{output_dir}/fastq_temp/s_4_2_sequence > $options{output_dir}/BWA_aln_Junction/s_4_2.sam";
#system("$cmd");

## Convert SAM files into ID sorted SAM files using SAMTools for genome alignment
$cmd = "$options{samtool_dir}/samtools view -bS $options{output_dir}/BWA_aln_genome/s_4_1.sam -o";
$cmd .= " $options{output_dir}/BWA_aln_genome/s_4_1.bam";
#system("$cmd");

$cmd = "$options{samtool_dir}/samtools view -bS $options{output_dir}/BWA_aln_genome/s_4_2.sam -o";
$cmd .= " $options{output_dir}/BWA_aln_genome/s_4_2.bam";
#system("$cmd");

$cmd = "$options{samtool_dir}/samtools sort -n -m 4000000000 $options{output_dir}/BWA_aln_genome/s_4_1.bam";
$cmd .= " $options{output_dir}/BWA_aln_genome/s_4_1_ID_sorted";
#system("$cmd");

$cmd = "$options{samtool_dir}/samtools sort -n -m 4000000000 $options{output_dir}/BWA_aln_genome/s_4_2.bam";
$cmd .= " $options{output_dir}/BWA_aln_genome/s_4_2_ID_sorted";
#system("$cmd");

$cmd = "$options{samtool_dir}/samtools view -X $options{output_dir}/BWA_aln_genome/s_4_1_ID_sorted.bam -o";
$cmd .= " $options{output_dir}/BWA_aln_genome/sorted_sam/s_4_1_ID_sorted.sam";
#system("$cmd");

$cmd = "$options{samtool_dir}/samtools view -X $options{output_dir}/BWA_aln_genome/s_4_2_ID_sorted.bam -o";
$cmd .= " $options{output_dir}/BWA_aln_genome/sorted_sam/s_4_2_ID_sorted.sam";
#system("$cmd");

## remove unwanted files here.
#system ("rm $options{outpout_dir}/BWA_aln_genome/*.sam*");
#system ("rm $options{outpout_dir}/BWA_aln_genome/*.bam*");


## Convert SAM files into ID sorted SAM files using SAMTools for Junction alignment
$cmd = "$options{samtool_dir}/samtools view -bS $options{output_dir}/BWA_aln_Junction/s_4_1.sam -o";
$cmd .= " $options{output_dir}/BWA_aln_Junction/s_4_1.bam";
#system("$cmd");

$cmd = "$options{samtool_dir}/samtools view -bS $options{output_dir}/BWA_aln_Junction/s_4_2.sam -o";
$cmd .= " $options{output_dir}/BWA_aln_Junction/s_4_2.bam";
#system("$cmd");

$cmd = "$options{samtool_dir}/samtools sort -n -m 4000000000 $options{output_dir}/BWA_aln_Junction/s_4_1.bam";
$cmd .= " $options{output_dir}/BWA_aln_Junction/s_4_1_ID_sorted";
#system("$cmd");

$cmd = "$options{samtool_dir}/samtools sort -n -m 4000000000 $options{output_dir}/BWA_aln_Junction/s_4_2.bam";
$cmd .= " $options{output_dir}/BWA_aln_Junction/s_4_2_ID_sorted";
#system("$cmd");

$cmd = "$options{samtool_dir}/samtools view -X $options{output_dir}/BWA_aln_Junction/s_4_1_ID_sorted.bam -o";
$cmd .= " $options{output_dir}/BWA_aln_Junction/sorted_sam/s_4_1_ID_sorted.sam";
#system("$cmd");

$cmd = "$options{samtool_dir}/samtools view -X $options{output_dir}/BWA_aln_Junction/s_4_2_ID_sorted.bam -o";
$cmd .= " $options{output_dir}/BWA_aln_Junction/sorted_sam/s_4_2_ID_sorted.sam";
#system("$cmd");

## remove unwanted files here.
#system ("rm $options{outpout_dir}/BWA_aln_Junction/*.sam*");
#system ("rm $options{outpout_dir}/BWA_aln_Junction/*.bam*");

## call script 1
$cmd = "perl 1_process_sorted_SAM_to_mapped_reads $options{output_dir}";
#system("$cmd");

## call script 2
$cmd = "perl 2_mapped_Reads_processed $options{output_dir}";
#system("$cmd");

## call script 3
$cmd = "perl 3_get_fusion_transcript_and_protein_results_NEW $options{output_dir}";
#system("$cmd");

timer(); #call timer to see when process ended.
exit();

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
}

#############################################################################
sub create_dir_struct {
	my $options = shift;

	system ("mkdir $options{output_dir}/BWA_aln_Junction");
	system ("mkdir $options{output_dir}/BWA_aln_Junction/sorted_sam");
	system ("mkdir $options{output_dir}/BWA_aln_genome");
	system ("mkdir $options{output_dir}/BWA_aln_genome/sorted_sam");
	system ("mkdir $options{output_dir}/map");
	system ("mkdir $options{output_dir}/fastq_temp");
	system ("mkdir $options{output_dir}/fasta_fusion_one_end_mapped");
	system ("mkdir $options{output_dir}/fusion_candidates");
	system ("mkdir $options{output_dir}/results");
	system ("mkdir $options{output_dir}/small_megablast_files");

}

#############################################################################
sub timer {
    my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
    my @weekDays = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
    my ($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
    my $year = 1900 + $yearOffset;
    my $theTime = "$hour:$minute:$second, $weekDays[$dayOfWeek] $months[$month] $dayOfMonth, $year";
    print "Time now: " . $theTime."\n";
}
