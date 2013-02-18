#!/usr/local/biotools/perl/5.14.2/bin/perl

=head1 NAME
   multisample_vqsr.pl

=head1 SYNOPSIS

    USAGE: multisample_vqsr.pl --input input_dir --output output_dir

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


=head1 EXAMPLE
   multisample_vqsr.pl --input=input_dir --output=output_dir

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

$ENV{'PERL5LIB'}="/projects/bsi/bictools/apps/variant_detection/vcftools/latest/perl";

my $cmd = "";
my $concat = "";

opendir (DIR, $options{input}) or die "Could not open dir $!\n";

my $count=0;
my @arr;

while (my $file = readdir(DIR)){
	next if ($file !~ /\.vcf$|\.vcf\.gz$/);

	## zip input vcf files
	if ($file !~ /\.gz$/) {
		$cmd = "$options{tabix}/bgzip $options{input}/$file";
		execute_cmd($cmd);
	}

	## create index of each vcf file.
	if ($file !~ /\.gz$/) {
		#### add .gz to file name.
		$cmd = "$options{tabix}/tabix -p vcf $options{input}/$file.gz";
		execute_cmd($cmd);

		$concat .= " $options{input}/$file.gz";
	} else {
		$cmd = "$options{tabix}/tabix -p vcf $options{input}/$file";
		execute_cmd($cmd);

		$concat .= " $options{input}/$file";
	}

	$count++;

	if($count > 500) {
		push(@arr, $concat);
		$count = 0;
		$concat = "";
	}
}
close(DIR);
my $i=1;

foreach my $c (@arr){
	## concatinate all indexed vcf files.
	$cmd = "$options{vcftools}/bin/vcf-concat $c > $options{output}/Merged.variant.raw.$i.vcf";
	execute_cmd($cmd);
	$i++;
}

exit();

## sort vcf file.
$cmd = "perl $options{gps_workflow}/vcfsort.pl $options{ref}.fai $options{output}/Merged.variant.raw.vcf";
$cmd .= " > $options{output}/Merged.variant.raw.vcf.temp";
execute_cmd($cmd);

$cmd = "mv $options{output}/Merged.variant.raw.vcf.temp $options{output}/Merged.variant.raw.vcf";
execute_cmd($cmd);

## remove multi allele
$cmd = "cat Merged.variant.raw.vcf | awk '\$0 ~ /^#/ || \$5 ~ /,/' > Merged.variant.raw.muilt.vcf";
execute_cmd($cmd);

$cmd = "cat Merged.variant.raw.vcf | awk '\$0 ~ /^#/ || \$5 !~ /,/' > Merged.variant.raw.vcf.tmp";
execute_cmd($cmd);

$cmd = "mv Merged.variant.raw.vcf.tmp Merged.variant.raw.vcf";
execute_cmd($cmd);

## split vcf file into indels and svn
$cmd = "perl $options{gps_workflow}/vcf_to_variant_vcf.pl -i=$options{output}/Merged.variant.raw.vcf";
$cmd .= " -v=$options{output}/Merged.variant.SNV.vcf -l=$options{output}/Merged.variant.INDEL.vcf -t=both";
execute_cmd($cmd);

## run vqsr training on SNV
my $check = 0;
$count = 0;
while ($check == 0 && $count < 10){
	$cmd = "$options{java}/java -Xmx3g -Xms512m -jar $options{gatk}/GenomeAnalysisTK.jar";
	$cmd .= " -R $options{ref} -et NO_ET";
	$cmd .= " -K $options{gatk}/Hossain.Asif_mayo.edu.key";
	$cmd .= " -T VariantRecalibrator";
	$cmd .= " -mode SNP";
	$cmd .= " -input $options{output}/Merged.variant.SNV.vcf";
	$cmd .= " -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $options{hapmap_vcf}";
	$cmd .= " -resource:omni,known=false,training=true,truth=false,prior=12.0 $options{omni_vcf}";
	$cmd .= " -resource:dbsnp,known=true,training=false,truth=false,prior=8.0 $options{dbSNP_vcf}";
	$cmd .= " -an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an DP";
	$cmd .= " -recalFile $options{output}/temp/Merged.variant.SNV.recal";
	$cmd .= " -tranchesFile $options{output}/temp/Merged.variant.SNV.tranches";
	$cmd .= " --maxGaussians 4";
	$cmd .= " --percentBadVariants 0.05";
	$cmd .= " -rscriptFile $options{output}/plot/Merged.variant.SNV.plots.R";
	execute_cmd($cmd);

	sleep 60;
	if (! -s "$options{output}/temp/Merged.variant.SNV.tranches.pdf" ){
		$cmd = "rm `grep -l $options{output}/plot/Merged.variant.SNV.plots.R *.log`";
		execute_cmd($cmd);

		$cmd = "rm core.*";
		execute_cmd($cmd);
	} else {
		$check = 1;
	}
	$count++;
}

## Apply Recalibrator
my $recal_count = `cat $options{output}/temp/Merged.variant.SNV.recal | wc -l`;
if (($recal_count > 2) && (-s "$options{output}/temp/Merged.variant.SNV.tranches")){
	$cmd = "$options{java}/java -Xmx6g -XX:-UseGCOverheadLimit -Xms512m";
	$cmd .= " -jar $options{gatk}/GenomeAnalysisTK.jar";
	$cmd .= " -R $options{ref} -et NO_ET";
	$cmd .= " -K $options{gatk}/Hossain.Asif_mayo.edu.key";
	$cmd .= " -mode SNP -T ApplyRecalibration";
	$cmd .= " -input $options{output}/Merged.variant.SNV.vcf";
	$cmd .= " -recalFile $options{output}/temp/Merged.variant.SNV.recal";
	$cmd .= " -tranchesFile $options{output}/temp/Merged.variant.SNV.tranches";
	$cmd .= " -o $options{output}/Merged.variant.filter.SNV.vcf";
	execute_cmd($cmd);

	sleep 60;
} else {
	print "ERROR: Failed to create SNV recal or tranches\n";
}

## run vqsr training on INDEL
$check = 0;
$count = 0;
while ($check == 0 && $count < 10){
	$cmd = "$options{java}/java -Xmx3g -Xms512m -jar $options{gatk}/GenomeAnalysisTK.jar";
	$cmd .= " -R $options{ref} -et NO_ET";
	$cmd .= " -K $options{gatk}/Hossain.Asif_mayo.edu.key";
	$cmd .= " -T VariantRecalibrator";
	$cmd .= " -mode INDEL";
	$cmd .= " -input $options{output}/Merged.variant.INDEL.vcf";
	$cmd .= " -resource:mills,VCF,known=true,training=true,truth=true,prior=12.0 $options{mills}";
	$cmd .= " -an QD -an FS -an HaplotypeScore -an ReadPosRankSum";
	$cmd .= " -recalFile $options{output}/temp/Merged.variant.INDEL.recal";
	$cmd .= " -tranchesFile $options{output}/temp/Merged.variant.INDEL.tranches";
	$cmd .= " --maxGaussians 4";
	$cmd .= " --percentBadVariants 0.05";
	$cmd .= " -rscriptFile $options{output}/plot/Merged.variant.INDEL.plots.R";
	execute_cmd($cmd);

	sleep 60;
	if (! -s "$options{output}/temp/Merged.variant.INDEL.tranches.pdf" ){
		$cmd = "rm `grep -l $options{output}/plot/Merged.variant.INDEL.plots.R *.log`";
		execute_cmd($cmd);

		$cmd = "rm core.*";
		execute_cmd($cmd);
	} else {
		$check = 1;
	}
	$count++;
}

## Apply Recalibrator
$recal_count = `cat $options{output}/temp/Merged.variant.INDEL.recal | wc -l`;
if (($recal_count > 2) && (-s "$options{output}/temp/Merged.variant.INDEL.tranches")){
	$cmd = "$options{java}/java -Xmx6g -XX:-UseGCOverheadLimit -Xms512m";
	$cmd .= " -jar $options{gatk}/GenomeAnalysisTK.jar";
	$cmd .= " -R $options{ref} -et NO_ET";
	$cmd .= " -K $options{gatk}/Hossain.Asif_mayo.edu.key";
	$cmd .= " -mode INDEL -T ApplyRecalibration";
	$cmd .= " -input $options{output}/Merged.variant.INDEL.vcf";
	$cmd .= " -recalFile $options{output}/temp/Merged.variant.INDEL.recal";
	$cmd .= " -tranchesFile $options{output}/temp/Merged.variant.INDEL.tranches";
	$cmd .= " -o $options{output}/Merged.variant.filter.INDEL.vcf";
	execute_cmd($cmd);

	sleep 60;
} else {
	print "ERROR: Failed to create INDEL recal or tranches\n";
}


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

	system("mkdir -p $options{output}/temp");
	system("mkdir -p $options{output}/plot");

	unless($options{java}) { $options{java} = "/usr/java/latest/bin"; }
	unless($options{gatk}) { $options{gatk} = "/projects/bsi/bictools/apps/alignment/GenomeAnalysisTK/latest"; }
	unless($options{tabix}) { $options{tabix} = "/projects/bsi/bictools/apps/alignment/tabix/latest"; }
	unless($options{vcftools}) { $options{vcftools} = "/projects/bsi/bictools/apps/variant_detection/vcftools/latest"; }
	unless($options{ref}) { $options{ref} = "/data2/bsi/reference/sequence/human/ncbi/37.1/allchr.fa"; }
	unless($options{gps_workflow}) { $options{gps_workflow} = "/projects/bsi/bictools/scripts/dnaseq/GENOME_GPS/trunk"; }
	unless($options{hapmap_vcf}) { $options{hapmap_vcf} = "/data2/bsi/reference/genetics/hapmap/hg19/hapmap_3.3.hg19.sites.vcf"; }
	unless($options{omni_vcf}) { $options{omni_vcf} = "/data2/bsi/reference/genetics/omni/hg19/1000G_omni2.5.hg19.sites.vcf"; }
	unless($options{dbSNP_vcf}) { $options{dbSNP_vcf} = "/data2/bsi/reference/annotation/dbSNP/hg19/dbsnp_135.hg19.vcf.gz"; }
	unless($options{mills}) { $options{mills} = "/data2/bsi/reference/misc/ALL.wgs.indels_mills_devine_hg19_leftAligned_collapsed.sites.vcf"; }
}


#############################################################################
sub execute_cmd {
	my $cmd = shift;

	print $cmd."\n";
	system($cmd);
}
