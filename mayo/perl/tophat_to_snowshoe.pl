#!/usr/local/biotools/perl/5.14.2/bin/perl

=head1 NAME
   tophat_to_snowshoe.pl

=head1 SYNOPSIS

    USAGE: tophat_to_snowshow.pl --result_txt=result.txt --snowshoe_xls=snowshoe.xls

=head1 OPTIONS

B<--result_txt,-rt>
   TopHat result.txt file

B<--outdir,-od>
   output directory


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
use MyBlat;
use Data::Dumper;
use Pod::Usage;
use Spreadsheet::ParseExcel::SaveParser;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);

my %options = ();
my $results = GetOptions (\%options,
                          'result_txt|rt=s',
						  'snowshoe_xls|sx=s',
						  'potential_fusion|pf=s',
						  'output_dir|o=s',
						  'threshold|t=s',
						  'help|h') || pod2usage();

## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

my $THRESHOLD = 10;
#############################################################################
## make sure everything passed was peachy
&check_parameters(\%options);

my $myblat = new MyBlat;
$myblat->init();

my $parser = Spreadsheet::ParseExcel::SaveParser->new();
my $workbook_orig = $parser->Parse($options{snowshoe_xls});
my $worksheet = $workbook_orig->worksheet(0);

my ($row_min, $row_max) = $worksheet->row_range();
$row_min +=1; #skip header.

my @ss_array_of_hash = ();

#create hash of xls file content.
for my $r ($row_min .. $row_max){
	push (@ss_array_of_hash, {SAMPLE=>$worksheet->get_cell($r,0)->value,
						   FUSION_PAIR_ALPHA=>(defined $worksheet->get_cell($r,1)) ? $worksheet->get_cell($r,1)->value : "",
						   FUSION_GENE_DIR=>(defined $worksheet->get_cell($r,2)) ? $worksheet->get_cell($r,2)->value : "",
						   TYPE=>(defined $worksheet->get_cell($r,3)) ? $worksheet->get_cell($r,3)->value : "",
						   POTENTIAL_FUSION_MECH=>(defined $worksheet->get_cell($r,4)) ? $worksheet->get_cell($r,4)->value : "",
						   FUSION_STRAND=>(defined $worksheet->get_cell($r,5)) ? $worksheet->get_cell($r,5)->value : "",
						   TOTAL_ENCOMPASSING=>(defined $worksheet->get_cell($r,6)) ? $worksheet->get_cell($r,6)->value : "",
						   TOTAL_SPLIT_READS=>(defined $worksheet->get_cell($r,7)) ? $worksheet->get_cell($r,7)->value : "",
						   EXON_BOUNDARY_FUSION=>(defined $worksheet->get_cell($r,9)) ? $worksheet->get_cell($r,9)->value : "",
						   EXON1=>(defined $worksheet->get_cell($r,10)) ? $worksheet->get_cell($r,9)->value : "",
						   EXON2=>(defined $worksheet->get_cell($r,11)) ? $worksheet->get_cell($r,11)->value : "",
						   PRIMER=>(defined $worksheet->get_cell($r,12)) ? $worksheet->get_cell($r,12)->value : "",
						   GENE_ORDER_ACCORDING_TO_PRIMER=>(defined $worksheet->get_cell($r,13)) ? $worksheet->get_cell($r,13)->value : "",
						   JUNCTION_MATCH=>"",
						   FOUND_IN_TOPHAT=>"",
						   FOUND_IN_POTENTIAL_TOPHAT=>"",
						   BLAT_HITS=>0});
}

open(RSTXT, "<", $options{result_txt}) or die "Could not open file to read $!\n";
open(TH_OUT, ">", $options{output_dir}."/tophat.txt") or die "Could not open file to write $!\n";
open(SS_OUT, ">", $options{output_dir}."/snowshoe.txt") or die "Could not open file to write $!\n";

my $header = join("\t","#Sample_name", "Fusion_pair_alphabetical", "Fusion_gene_direction",
					   "Type", "Potential_fusion_mech", "Fusion_strand", "Total_encompassing_reads",
					   "Total_split_reads", "Exon_boundary_fusion", "Exon1", "Exon2", "Primer");

print TH_OUT $header ."\tFound_in_SnowShoes\tFOUND_IN_POTENTIAL_SnowShoes\tJunction_match\tBLAT_HITS\n";

while(<RSTXT>) {
	chomp $_;

	my @data = split(/\t/, $_);
	my ($direction, $primer) = potentialFusion($data[0], $data[1], $data[4], $data[3], $data[6]);

	print TH_OUT $data[0]."\t"; #sample

	#fusion pair alphabetical
	if ($data[1] lt $data[4]){
		print TH_OUT $data[1]."_".$data[4]."\t";
	} else {
		print TH_OUT $data[4]."_".$data[1]."\t";
	}

	print TH_OUT $data[1]."->".$data[4]."\t"; #fusion gene directional

	#type
	if ($data[2] =~ $data[5]){
		print TH_OUT "intra\t";

		#potential fusion mechanism
		print TH_OUT "T";

		if ($direction =~ /fr|rf/i){
			print TH_OUT " and I"
		} elsif (($data[3] < $data[6]) && ($direction =~ /ff/i)){
			print TH_OUT " and D";
		} elsif (($data[3] > $data[6]) && ($direction =~ /rr/i)){
			print TH_OUT " and D";
		}

		print TH_OUT "\t";
	} else {
		print TH_OUT "inter\t";

		#potential fusion mechanism
		print TH_OUT "T\t";
	}

	#strand
	if (substr($direction,0,1) =~ /f/i){
		print TH_OUT "+\t";
	} else { print TH_OUT "-\t"; }

	my $split_reads = $data[7];
	my $encompassing = $data[8]+$data[9];

	print TH_OUT ($split_reads + $encompassing) ."\t"; #total
	print TH_OUT $encompassing ."\t"; #total_encompassing
	print TH_OUT $split_reads ."\t"; #split_reads

	print TH_OUT "NA\t"; #exon boundary fusion

	my ($inSS, $junc, $exon1, $exon2) = checkSnowShoe($data[0], $data[1], $data[4], $data[3], $data[6]);
	print TH_OUT "$exon1\t"; #exon 1
	print TH_OUT "$exon2\t"; #exon 2
	print TH_OUT $primer."\t"; # primer

	#seen in Snowshoe
	print TH_OUT "$inSS\t";
	print TH_OUT "\t"; #BLANK column
	print TH_OUT "$junc\t"; #junction match

	$myblat->execute($primer, $options{output_dir});
	print TH_OUT $myblat->numOfHits(90, $options{output_dir})."\n"; #blat verificaiton.
}

close(RSTXT);

#look for snowshoe results in potential fusion file.
for (my $i=0; $i<scalar(@ss_array_of_hash); $i++){
	#look for snowshoe results that have not been see in results.txt file.

	if ((!length($ss_array_of_hash[$i]->{FOUND_IN_TOPHAT})) && (length($ss_array_of_hash[$i]->{FUSION_GENE_DIR}))){
		my ($exon1,$exon2) = split(/->/, $ss_array_of_hash[$i]->{FUSION_GENE_DIR});

		#get all potential fusion entires for given exons.
		my $potential_fusion = `grep -A 1 -B 4 "^$exon1.*$exon2" $options{potential_fusion}`;
		chomp $potential_fusion;

		if (length($potential_fusion)){
			#get sample name
			my ($sample,$rest) = split(/\s/, $potential_fusion, 2);

			#split potential fusion based on sample name if there are
			#multiple entries for the same two exons
			my @fusions = split(/$sample/, $potential_fusion);

			#loop through all fusion to check coordinates within threshold
			#first array index "0" is always empty.
			for (my $j=1; $j<scalar(@fusions); $j++) {
				my ($line, $rest) = split (/\n/, $fusions[$j],2);

				$line =~ s/^\s+//;
				$line =~ s/\s+$//;
				my ($chroms, $coord1, $coord2, $direction) = split(/\s+/, $line, 5);

				my ($in_potential, $junction_match) = checkFusionCoordinates($ss_array_of_hash[$i]->{EXON1}, $ss_array_of_hash[$i]->{EXON2}, $exon1, $coord1, $coord2);

				$ss_array_of_hash[$i]->{FOUND_IN_POTENTIAL_TOPHAT} = ($in_potential) ? "YES" : "NO";
				$ss_array_of_hash[$i]->{JUNCTION_MATHC} = ($junction_match) ? "YES" : "NO";
				$ss_array_of_hash[$i]->{FOUND_IN_TOPHAT} = 0;

				$myblat->execute($ss_array_of_hash[$i]->{PRIMER}, $options{output_dir});
				$ss_array_of_hash[$i]->{BLAT_HITS} = $myblat->numOfHits(90, $options{output_dir});
			}
		}
	}
}

for (my $i=0; $i<scalar(@ss_array_of_hash); $i++){
	print SS_OUT $ss_array_of_hash[$i]->{SAMPLE}."\t";
	print SS_OUT $ss_array_of_hash[$i]->{FUSION_PAIR_ALPHA}."\t";
	print SS_OUT $ss_array_of_hash[$i]->{FUSION_GENE_DIR}."\t";
	print SS_OUT $ss_array_of_hash[$i]->{TYPE}."\t";
	print SS_OUT $ss_array_of_hash[$i]->{POTENTIAL_FUSION_MECH}."\t";
	print SS_OUT $ss_array_of_hash[$i]->{FUSION_STRAND}."\t";
	print SS_OUT $ss_array_of_hash[$i]->{TOTAL_ENCOMPASSING}."\t";
	print SS_OUT $ss_array_of_hash[$i]->{TOTAL_SPLIT_READS}."\t";
	print SS_OUT $ss_array_of_hash[$i]->{EXON_BOUNDARY_FUSION}."\t";
	print SS_OUT $ss_array_of_hash[$i]->{EXON1}."\t";
	print SS_OUT $ss_array_of_hash[$i]->{EXON2}."\t";
	print SS_OUT $ss_array_of_hash[$i]->{PRIMER}."\t";
	print SS_OUT $ss_array_of_hash[$i]->{FOUND_IN_TOPHAT}."\t";
	print SS_OUT $ss_array_of_hash[$i]->{FOUND_IN_POTENTIAL_TOPHAT}."\t";
	print SS_OUT $ss_array_of_hash[$i]->{JUNCTION_MATCH}."\t";
	print SS_OUT $ss_array_of_hash[$i]->{BLAT_HITS}."\n";
}

close(SS_OUT);
close(TH_OUT);
#############################################################################
sub check_parameters {
    my $options = shift;

	my @required = ("result_txt", "snowshoe_xls", "potential_fusion", "output_dir");

	foreach my $key (@required) {
		unless ($options{$key}) {
			print STDERR "ARG: $key is required\n";
			pod2usage({-exitval => 2,  -message => "error message", -verbose => 1, -output => \*STDERR});
			exit(-1);
		}
	}

	if (defined $options{threshold}) {
		$THRESHOLD = $options{threshold};
	}
}

#############################################################################
# get primer and direction from potential_fusion.txt
# given exons and coordinates from fusion results.txt file
sub potentialFusion{
	my ($sample, $exon1, $exon2, $coord1, $coord2) = @_;
	my $direction = "";
	my $primer = "";

	#get primer from potential fusion file.
	my $fusion_detail = `grep -A 1 -B 4 "^$exon1.*$exon2" $options{potential_fusion}`;

	if ($fusion_detail =~ /.*($coord1)\s($coord2)\s(\w\w).*\n(\w+\s\w+)\n(\w+\s\w+)\n.*/){
		$direction = $3;

		#primer is 1 and 4 of the sequence set in potentail fusion
		# ----------1----------- ------------2-----------
		# ----------3----------- ------------4-----------

		my @p = split(/\s/,$4);
		$primer = $p[0];

		@p = split(/\s/,$5);
		$primer .= $p[1];
	}

	return ($direction, $primer);
}

#############################################################################
# given tophat exons and coordinates
# check in snowshoe hash for same fusion, if found check if its with in defined
# threshold.
sub checkSnowShoe{
	my ($sample, $exon1, $exon2, $coord1, $coord2) = @_;

	my $f_gene_dir = $exon1."->".$exon2;

	for (my $i=0; $i<scalar(@ss_array_of_hash); $i++){

		# checking correct sample name. results.txt and snowshoes fusion.xls both
		# can contain fusion report for multiple samples.
		if (($ss_array_of_hash[$i]->{SAMPLE} =~ /$sample/i) && (length($ss_array_of_hash[$i]->{FUSION_GENE_DIR}))) {
			$myblat->execute($ss_array_of_hash[$i]->{PRIMER},$options{output_dir});
			$ss_array_of_hash[$i]->{BLAT_HITS} = $myblat->numOfHits(90, $options{output_dir});

			if ($ss_array_of_hash[$i]->{FUSION_GENE_DIR} =~ /$f_gene_dir/i){
				my ($in_tophat, $junction_match) = checkFusionCoordinates($ss_array_of_hash[$i]->{EXON1}, $ss_array_of_hash[$i]->{EXON2}, $exon1, $coord1, $coord2);

				$ss_array_of_hash[$i]->{FOUND_IN_TOPHAT} = ($in_tophat) ? "YES" : "NO";
				$ss_array_of_hash[$i]->{FOUND_IN_POTENTIAL_TOPHAT} = ($in_tophat) ? "YES" : "NO";
				$ss_array_of_hash[$i]->{JUNCTION_MATCH} = ($junction_match) ? "YES" : "NO";

				my @e1 = split(/:/, $ss_array_of_hash[$i]->{EXON1});
				my @e2 = split(/:/, $ss_array_of_hash[$i]->{EXON2});

				if ($exon1 =~ /$e2[2]/){
					@e1 = split(/:/, $ss_array_of_hash[$i]->{EXON2});
					@e2 = split(/:/, $ss_array_of_hash[$i]->{EXON1});
				}

				if ($in_tophat && $junction_match) {
					return ("YES", "YES", $ss_array_of_hash[$i]->{EXON1}, $ss_array_of_hash[$i]->{EXON2});
				} elsif ($in_tophat && ! $junction_match) {
					return ("YES", "NO", $e1[1].":".$e1[2].".:".$e1[-1], $e2[1].":".$e2[2].".:".$e2[-1]);
				} else {
					return ("NO", "NO", $e1[1].":".$e1[2].".:".$e1[-1], $e2[1].":".$e2[2].".:".$e2[-1]);
				}
			}
		}
	}

	return ("","","",""); # no match in snowshoe, novel to tophat.
}

sub checkFusionCoordinates{
	my ($exon1, $exon2, $th_exon1, $coord1, $coord2) = @_;

	my @e1 = split(/:/, $exon1);
	my @e2 = split(/:/, $exon2);

	if ($th_exon1 =~ /$e2[2]/){
		@e1 = split(/:/, $exon2);
		@e2 = split(/:/, $exon1);
	}

	#init to any number greater than 10;
	my $diff1 = 100;
	my $diff2 = 100;

	if (($e1[-1] =~ /\+/) && ($e2[-1] =~ /\+/)){
		# when --------> and -------->
		$diff1 = $e1[-2] - ($coord1+1);
		$diff2 = $e2[-3] - $coord2;
	} elsif (($e1[-1] =~ /\+/) && ($e2[-1] =~ /-/)){
		# when --------> and <--------
		$diff1 = $e1[-2] - ($coord1+1);
		$diff2 = $e2[-2] - ($coord2+1);
	} elsif (($e1[-1] =~ /-/) && ($e2[-1] =~ /-/)){
		# when <-------- and <--------
		$diff1 = $e1[-3] - $coord1;
		$diff2 = $e2[-2] - ($coord2+1);
	} elsif (($e1[-1] =~ /-/) && ($e2[-1] =~ /\+/)){
		# when <-------- and -------->
		$diff1 = $e1[-3] - $coord1;
		$diff2 = $e2[-3] - $coord2;
	}

	if (($diff1 == 0) && ($diff2 == 0)) {
		return(1, 1); # exact match so jucntion match is true
	} elsif (($diff1 < $THRESHOLD) && ($diff2 < $THRESHOLD)) {
		return(1, 0); # not an exact match
	} else {
		return(0, 0); # no match within threshold
	}
}
