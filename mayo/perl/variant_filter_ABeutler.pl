#!/usr/local/biotools/perl/5.14.2/bin/perl

use strict;
use warnings;

my $ref = "/data2/bsi/reference/sequence/human/ncbi/37.1/allchr.fa";
my $gatk = "/projects/bsi/bictools/apps/alignment/GenomeAnalysisTK/1.6-7-g2be5704/";
my $out = "/data2/bsi/secondary/Beutler_Andreas_m068039/exome/120627_SN616_0194_AC0Y5DACXX_combine_run/Reports_per_Sample";

opendir (DIR, "$out") or die "Could not open dir\n";

my $c=0;
while (my $file = readdir(DIR)) {
	#exit if ($c > 1);
	next if ($file !~ /variants.filter.vcf$/);

	print "Processing $out/$file\n";

	sleep 5;
	my $cmd = "qsub -V -wd /data2/bsi/secondary/Beutler_Andreas_m068039/exome/TRUSEQ_MULTICELL/custom_filters/logs ";
	$cmd .= "-q 1-day -l h_vmem=3G -l h_stack=10M -m ae -M bhavsar.jaysheel\@mayo.edu /home/m101236/jbhavsar/scripts/amit_var_filter.sh $gatk $ref $out $file";

	print $cmd."\n";
	system("$cmd");

	$c++;
}

exit(0);
