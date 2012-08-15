#!/usr/local/biotools/perl/5.14.2/bin/perl

use strict;
use warnings;
use Data::Dumper;
use File::Basename;

my $dir = $ARGV[0];
my $ext = $ARGV[1];

opendir (DIR, $dir) or die "Could not open file\n";

my $sample = "";
my $lane = "";
my $index = "";
my %sample_hash = ();

while (my $file = readdir(DIR)) {
	next if ($file !~ /$ext$/);

	# expecting file in following format
	# some_unique_id.FLOWCELLID_LANE_INDEX.EXT
	my @bits = split(/\./, $file);
	my @flow = split(/_/, $bits[1]);

	# in case a sample sample was run multiple times
	# use id and lane number as sample name
	$bits[0] =~ s/^s_//;
	$bits[0] =~ s/L\d$//;
	my $key = $bits[0];#. "" .$flow[1];
	$key =~ s/_|-//g;
	$sample .= $key. ":";
	$lane .= substr($flow[1], 1). ":";
	$index .= substr($flow[2], 1). ":";

	$sample_hash{$key} .= $file ."\t";
}

# remove all -/_ from sample
$sample =~ s/:$//; # remove last :
$lane =~ s/:$//;
$index =~ s/:$//;

#print sample info for run_info file
print "SAMPLENAMES=" .$sample. "\n";
print "LANEINDEX=" .$lane. "\n";
print "LABINDEXES=" .$index. "\n";

# print sample info for sample_info file
foreach my $key (keys %sample_hash){
	$sample_hash{$key} =~ s/\t$//;
	print uc($ext). ":" .$key. "=" . $sample_hash{$key}."\n";
}
