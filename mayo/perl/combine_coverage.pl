#!/usr/local/biotools/perl/5.14.2/bin/perl

use strict;
use warnings;

my $out = "/data2/bsi/secondary/Beutler_Andreas_m068039/exome/120726_SN616_0199_BC117YACXX.combine-run/numbers";

opendir (DIR, "$out") or die "Could not open dir\n";
open (OUT, ">", "$out/all_coverage_data.tsv") or die "Could not create output file\n";
my $c=0;
while (my $file = readdir(DIR)) {
	#exit if ($c > 1);
	next if ($file !~ /coverage.out$/);

	print "Processing $out/$file\n";

	#sleep 5;
	open (FHD, "<", "$out/$file") or die "Could not open file $out/$file\n";
	my @data = <FHD>;
	close FHD;

	$file =~ s/.coverage.out$//;
	my $str = "$file\t";
	foreach (@data){
		chomp $_;
		$str .= $_."\t";
	}
	$str =~ s/\t$/\n/;

	print OUT $str;
	$c++;
}

exit(0);
