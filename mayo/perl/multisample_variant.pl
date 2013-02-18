#!/usr/local/biotools/perl/5.14.2/bin/perl

=head1 NAME
   multisample_variant.pl

=head1 SYNOPSIS

    USAGE: multisample_variant.pl --input input_dir --output output_dir

=head1 OPTIONS

B<--input,-i>
   Input dir of all bam files

B<--output_dir,-o>
   output dir

B<--help,-h>
   This help message

=head1  DESCRIPTION


=head1  INPUT


=head1  OUTPUT


=head1  CONTACT
  Jaysheel D. Bhavsar @ bjaysheel[at]gmail[dot]com


==head1 EXAMPLE
   multisample_variant.pl --input=input_dir --output=output_dir

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
						  'ref|r=s',
                          'help|h') || pod2usage();

## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}
#############################################################################
## make sure everything passed was peachy
&check_parameters(\%options);

#my @chrs = qw(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM);
#my @chrs = qw(chr1 chr2 chr3 chr8 chr9 chr10 chr12 chr15 chr20 chr21 chrX);
my @chrs = qw(chrX);
my $bam = "";
my $hold = "";
my $cmd = "";

my $qsub = "qsub -V -wd $options{output}/logs -l h_vmem=100G -l h_stack=10M -b y -q ngs-rand -m ae -M bhavsar.jaysheel\@mayo.edu";

opendir (DIR, $options{input}) or die "Could not open dir $!\n";
while (my $file = readdir(DIR)){
	next if ($file =~ /^\./);
	next if ($file =~ /\.bai$/);

	## sample exception
	next if ($file =~ /(ABP9In25|ABP9In27|ABP9In09|ABP9In03|ABP8In22|ABP8In23).*/);
	next if ($file =~ /(ABP8In11|ABP8In20|ABP8In21|ABP8In10).*/);

	$bam .= "-I $options{input}/$file ";
}

foreach my $chr (@chrs) {
	my $c = 0;

	opendir(DIR, $options{target_bed}) or die "Could not open dir $options{target_dir}\n";

	while (my $file = readdir(DIR)) {

		#### skip . and .. file;
		next if ($file =~ /^\.$|^\.\.$/);

		#### only process bed files from unix split command (a special case);
		next if ($file !~ /^$chr[\w\w\w\w]/);

		#### rename split bed file as unix split command does not add extension
		if ($file !~ /\.bed$/) {
			print "Adding .bed ext to file: $file\n";

			system("mv $options{target_bed}/$file $options{target_bed}/$file.bed");
			$file .= ".bed";
		}

		#### skip files that have already been ran.
		next if ((-e "$options{output}/$chr/MultiSample_$file.vcf.idx") && (-s "$options{output}/$chr/MultiSample_$file.vcf.idx"));

		#### get job name
		my $job_name = "MULTISAMPLE.VARIANT.$file";

		$hold .= "$job_name,";


		##### setup unified genotyper command.
		$cmd = "$options{java}/java -Xmx98g -Xms512m -Djava.io.tmpdir=$options{output}/temp/$chr";
		$cmd .= " -jar $options{gatk}/GenomeAnalysisTK.jar";
		$cmd .= " -K $options{gatk}/Hossain.Asif_mayo.edu.key";
		$cmd .= " -R $options{ref} -et NO_ET";
		$cmd .= " -T UnifiedGenotyper";
		$cmd .= " --max_alternate_alleles 5";
		$cmd .= " --min_base_quality_score 20 -glm BOTH";
		#$cmd .= " -L $options{target_bed}/${chr}.prob.target.bed";
		$cmd .= " -L $options{target_bed}/$file";
		$cmd .= " $bam";
		$cmd .= " --out $options{output}/$chr/MultiSample_$file.vcf";

		$cmd = $qsub ." -N $job_name ". $cmd;

		sleep 5;
		execute_cmd($cmd);
	}

	close(DIR);
}

$cmd = "perl /projects/bsi/bictools/scripts/dev/jbhavsar/scripts/mayo/perl/multisample_vqsr.pl -i=$options{output} -o=$options{output}";
$cmd = $qsub . " -hold_jid $hold -N MULTISAMPLE_VQSR " . $cmd;
#execute_cmd($cmd);

exit(0);

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

	system("mkdir -p $options{output}/temp/");
	system("mkdir -p $options{output}/logs/");

	unless($options{java}) { $options{java} = "/usr/java/latest/bin"; }
	unless($options{gatk}) { $options{gatk} = "/projects/bsi/bictools/apps/alignment/GenomeAnalysisTK/latest"; }
	unless($options{ref}) { $options{ref} = "/data2/bsi/reference/sequence/human/ncbi/37.1/allchr.fa"; }
	unless($options{target_bed}) { $options{target_bed} = "/data2/bsi/secondary/Beutler_Andreas_m068039/exome/multisample/TruSeq/target_bed"; }
}

#############################################################################
sub execute_cmd {
	my $cmd = shift;

	print $cmd."\n";
	system($cmd);
}
