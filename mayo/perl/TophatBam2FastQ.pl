#!/usr/local/biotools/perl/5.14.2/bin/perl

=head1 NAME
   TophatBam2FastQ.pl

=head1 SYNOPSIS
    USAGE: TophatBam2FastQ.pl -u=unmapped.bam -a=accepted_hits.bam -o=output_dir -s=sample_name

=head1 OPTIONS

B<--unmapped, -u>
	unmapped bam file

B<--accepted, -a>
	accpeted hits bam file

B<--output_dir, -o>
	output dir

B<--help,-h>


=head1  DESCRIPTION
	Convert Tophat 2 aligment output to fastq files

=head1  INPUT
	Tophat accepted_hits.bam and unmapped_hits.bam along with sample_name to be used
	as prefix for R1 and R2 fastq files.

=head1  OUTPUT
	SAMPLENAME.R1.fastq and SAMPLENAME.R2.fastq

=head1  CONTACT
  bjaysheel@gmail.com


==head1 EXAMPLE
	./TophatBam2FastQ.pl -u=unmapped.bam -a=accepted_hits.bam -o=output -s=ABCD

=cut

use strict;
use warnings;
use Data::Dumper;
use Pod::Usage;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);

my %options = ();
my $results = GetOptions (\%options,
                          'unmapped|u=s',
						  'accepted|a=s',
						  'output_dir|o=s',
						  'sample_name|s=s',
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

my $this = {};
my $cmd = "";

## Unmapped sequence ids have /1 and /2 suffix which need to be removed
## from read id as bam2fastx will add that suffix based on sam bit flag.

## get only R1 reads from unmapped.bam file
$cmd = "$options{samtools}/samtools view -f 0x40 $options{unmapped}";
$cmd .= " | sed -e 's/\\/[0-9]//' > $options{output_dir}/unmapped.edited.R1.sam";
execute_cmd($cmd);

## get only R2 reads from unmapped.bam file
$cmd = "$options{samtools}/samtools view -f 0x80 $options{unmapped}";
$cmd .= " | sed -e 's/\\/[0-9]//' > $options{output_dir}/unmapped.edited.R2.sam";
execute_cmd($cmd);

## get accepted hits header to be used later while creating bam file
$cmd = "$options{samtools}/samtools view -H $options{accepted} > $options{output_dir}/accepted.header.sam";
execute_cmd($cmd);

## get all primary R1 non fusion reads.
$cmd = "$options{samtools}/samtools view -f 0x40 -F 0x100 $options{accepted}";
$cmd .= " | awk '{ if (\$0 !~ \"XF:Z\") {print} }'";
$cmd .= " > $options{output_dir}/accepted_hits.nofusion.R1.sam";
execute_cmd($cmd);

## get all primary R2 non fusion reads.
$cmd = "$options{samtools}/samtools view -f 0x80 -F 0x100 $options{accepted}";
$cmd .= " | awk '{ if (\$0 !~ \"XF:Z\") {print} }'";
$cmd .= " > $options{output_dir}/accepted_hits.nofusion.R2.sam";
execute_cmd($cmd);

## get all R1 fusion reads
$cmd = "$options{samtools}/samtools view -f 0x40 $options{accepted}";
$cmd .= " | awk '{ if (\$0 ~ \"XF:Z\") {print} }'";
$cmd .= " > $options{output_dir}/accepted_hits.fusion.R1.sam";
execute_cmd($cmd);

## get all R2 fusion reads
$cmd = "$options{samtools}/samtools view -f 0x80 $options{accepted}";
$cmd .= " | awk '{ if (\$0 ~ \"XF:Z\") {print} }'";
$cmd .= " > $options{output_dir}/accepted_hits.fusion.R2.sam";
execute_cmd($cmd);

## sort all fusion reads by read id so that durion fusionmerge all duplicate entries are
## ignored
$cmd = "sort -k1 $options{output_dir}/accepted_hits.fusion.R1.sam > $options{output_dir}/accepted_hits.fusion.R1.sam.tmp";
$cmd .= "; mv $options{output_dir}/accepted_hits.fusion.R1.sam.tmp $options{output_dir}/accepted_hits.fusion.R1.sam";
execute_cmd($cmd);

## sort all fusion reads by read id so that durion fusionmerge all duplicate entries are
## ignored
$cmd = "sort -k1 $options{output_dir}/accepted_hits.fusion.R2.sam > $options{output_dir}/accepted_hits.fusion.R2.sam.tmp";
$cmd .= "; mv $options{output_dir}/accepted_hits.fusion.R2.sam.tmp $options{output_dir}/accepted_hits.fusion.R2.sam";
execute_cmd($cmd);

## merge fusion reads R1
tophatFusionMerge("$options{output_dir}/accepted_hits.fusion.R1.sam", "$options{output_dir}/accepted_hits.fusion.R1.edited.sam");

## merge fusion reads R2
tophatFusionMerge("$options{output_dir}/accepted_hits.fusion.R2.sam", "$options{output_dir}/accepted_hits.fusion.R2.edited.sam");

## stitch together all R1 sam files.
$cmd = "cat $options{output_dir}/accepted.header.sam $options{output_dir}/accepted_hits.fusion.R1.edited.sam";
$cmd .= " $options{output_dir}/accepted_hits.nofusion.R1.sam";
$cmd .= " $options{output_dir}/unmapped.edited.R1.sam > $options{output_dir}/stiched.R1.sam";
execute_cmd($cmd);

## stich together all R2 sam files.
$cmd = "cat $options{output_dir}/accepted.header.sam $options{output_dir}/accepted_hits.fusion.R2.edited.sam";
$cmd .= " $options{output_dir}/accepted_hits.nofusion.R2.sam";
$cmd .= " $options{output_dir}/unmapped.edited.R2.sam > $options{output_dir}/stiched.R2.sam";
execute_cmd($cmd);

## convert R1 sam to bam
$cmd = "$options{samtools}/samtools view -bS $options{output_dir}/stiched.R1.sam > $options{output_dir}/stiched.R1.bam";
execute_cmd($cmd);

## convert R2 sam to bam
$cmd = "$options{samtools}/samtools view -bS $options{output_dir}/stiched.R2.sam > $options{output_dir}/stiched.R2.bam";
execute_cmd($cmd);

## sort R1 bam by id.
$cmd = "/usr/java/latest/bin/java -Xmx6g -Xms512m";
$cmd .= " -jar $options{picard}/SortSam.jar";
$cmd .= " INPUT=$options{output_dir}/stiched.R1.bam";
$cmd .= " OUTPUT=$options{output_dir}/stiched.R1.sorted.id.bam";
$cmd .= " SO=queryname MAX_RECORDS_IN_RAM=1000000";
$cmd .= " TMP_DIR=$options{output_dir}/tmp/ VALIDATION_STRINGENCY=SILENT";
execute_cmd($cmd);

## sort R2 bam by id.
$cmd = "/usr/java/latest/bin/java -Xmx6g -Xms512m";
$cmd .= " -jar $options{picard}/SortSam.jar";
$cmd .= " INPUT=$options{output_dir}/stiched.R2.bam";
$cmd .= " OUTPUT=$options{output_dir}/stiched.R2.sorted.id.bam";
$cmd .= " SO=queryname MAX_RECORDS_IN_RAM=1000000";
$cmd .= " TMP_DIR=$options{output_dir}/tmp/ VALIDATION_STRINGENCY=SILENT";
execute_cmd($cmd);

## convert bam to fastq
$cmd = "$options{tophat}/bam2fastx -q -A -Q -o $options{output_dir}/$options{sample_name}.R1.fastq $options{output_dir}/stiched.R1.sorted.id.bam";
execute_cmd($cmd);

## convert bam to fastq
$cmd = "$options{tophat}/bam2fastx -q -A -Q -o $options{output_dir}/$options{sample_name}.R2.fastq $options{output_dir}/stiched.R2.sorted.id.bam";
execute_cmd($cmd);

$cmd = "rm $options{output_dir}/*.sam $options{output_dir}/*.bam";
execute_cmd($cmd);

my $num_r1 = `wc -l $options{output_dir}/$options{sample_name}.R1.fastq`;
my $num_r2 = `wc -l $options{output_dir}/$options{sample_name}.R2.fastq`;

print "Number of reads in $options{sample_name}.R1.fastq: $num_r1\n";
print "Number of reads in $options{sample_name}.R2.fastq: $num_r2\n";

exit();

#############################################################################
sub check_parameters {
    my $options = shift;

	my @required = qw(unmapped accepted output_dir sample_name);

	foreach my $key (@required) {
		unless ($options{$key}) {
			print STDERR "ARG: $key is required\n";
			pod2usage({-exitval => 2,  -message => "error message", -verbose => 1, -output => \*STDERR});
			exit(-1);
		}
	}

	$options{'samtools'} = "/projects/bsi/bictools/apps/alignment/samtools/latest" unless ($options{'samtools'});
	$options{'picard'} = "/projects/bsi/bictools/apps/alignment/picard/latest" unless ($options{'picard'});
	$options{'tophat'} = "/projects/bsi/bictools/apps/alignment/tophat/latest" unless ($options{'tophat'});
}

#############################################################################
sub execute_cmd {
	my $cmd = shift;

	print STDOUT $cmd."\n";
	system($cmd);
}

#############################################################################
sub create_dir_struct {
	my $options = shift;

	my $dir = "$options{output_dir}/tmp";
	if ( -d $dir) {
		print STDERR "Directory $dir exist";
	} else {
		execute_cmd("mkdir -p $dir");
	}
}

#############################################################################
sub tophatFusionMerge {
	my $input = shift;
	my $output = shift;

	open(FHD, "<", $input) or die("Could not open $input");
	open(OUT, ">", $output) or die("Could not write file $output");

	my $prev="nnnn";  ## cannot be empty string or else regex will be invalid.
	while (<FHD>){
		chomp $_;
		my @fields = split(/\t/, $_);

		unless($fields[0] =~ /$prev/i) {
			my $i=0;
			my @nfields;

			## temp store of first 10 fileds
			for (my $i=0; $i <10; $i++) {
				$nfields[$i] = $fields[$i];
			}

			## find the field that contatins fusion info and original read
			## i.e field starting with XF
			my $XF_field;
			foreach my $f (@fields) {
				if ($f =~ /^XF/) {
					$XF_field = $f;
					last;
				}
			}

			my @xfields = split (/\s/, $XF_field);
			my ($fchr, $schr) = split (/-/,$xfields[1]);

			## update CIGAR string by removing F flag
			(my $cigar = $xfields[3]) =~ s/F/N/;

			## update fileds with appropriate chr and CIGAR string based
			## on read value from fusion field
			$nfields[2]=$fchr;
			$nfields[3]=$xfields[2];
			$nfields[5]=uc $cigar;
			$nfields[9]=$xfields[4];
			$nfields[10]=$xfields[5];

			## print new read to sam file.
			print OUT join("\t", @nfields)."\n";
		}
		$prev = $fields[0];
	}
	close(FHD);
	close(OUT);
}
