#!/usr/bin/perl

=head1 NAME
   bam_breaker.pl

=head1 SYNOPSIS

    USAGE: bam_breaker.pl --input input_bam_file [--fragment 50]

=head1 OPTIONS

B<--input,-i>
   Input file

B<--fragment,-f>
   Fragment size

B<--help,-h>
   This help message

=head1  DESCRIPTION

=head1  INPUT


=head1  OUTPUT


=head1  CONTACT
  Jaysheel D. Bhavsar @ bjaysheel[at]gmail[dot]com


==head1 EXAMPLE
   bam_breaker.pl --input filename.bam

=cut

use strict;
use warnings;
use Data::Dumper;
use Cwd;
use Pod::Usage;
use Scalar::Util qw(looks_like_number);
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);

my %options = ();
my $results = GetOptions (\%options,
                          'input|i=s',
						  'output|o=s',
						  'fragment|f=s',
						  'help|h') || pod2usage();

## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}
#############################################################################
## make sure everything passed was peachy
&check_parameters(\%options);

my $tmp_sam = getcwd()."/tmp_sam_pos.sam";

#execute samtools to convert bam file into sam.
my $cmd = "samtools view -h $options{input} > $tmp_sam";
system($cmd);

open (FHD, "<", $tmp_sam) or die "Cound not open file $tmp_sam\n";
open (OUT, ">", $options{output}) or die "Cound not open file $options{output}\n";

my $line = 0;

my $rev_max = 0;
my $for_max = 0;
my @max_array = ();

while (<FHD>){
	$line++;

	#skip processing lines starting with @ just print to output file.
	if ($_ =~ /^@/){
		print OUT $_;
		#print $_;
		next;
	}

	chomp $_;

	my @data = split(/\t/, $_);

#	if (($rev_max < length($data[9])) && ($data[1] == 16)){
#		$max_array[0] = $_;
#		$rev_max = length($data[9]);
#	}
#
#	if (($for_max < length($data[9])) && ($data[1] != 16)){
#		$max_array[1] = $_;
#		$for_max = length($data[9]);
#	}
#}

#foreach my $max(@max_array){

#	my @data = split(/\t/,$max);

#	print $data[0]."\t".length($data[9])."\n";

	#counter and local vars
	my $start = 0;
	my $frag_count = 0;
	my $op_count = 0;
	my $prev_delta = 0;
	my $curr_delta = 0;
	my $operations = "";
	my @arr = ();

	#re-write operations with a delimiter and remove last comma
	$data[5] =~ s/([A-Za-z])/$1,/g;
	$data[5] =~ s/,$//;

	#create array of hashes of operations
	map { push @arr, {type=>$2, value=>$1} if (/(\d+)([A-Za-z])/g) } split(/,/,$data[5]);

	#loop through each operation
	foreach my $op (@arr){

		if ($op->{type} =~ /D/i) {
			$curr_delta += $op->{value};
		} elsif ($op->{type} =~ /I/i){
			$curr_delta -= $op->{value};
		}

		#skip any delete operations
		if ($op->{type} =~ /D/i){
			$operations .= $op->{value}."".$op->{type};
		} else {
			$op_count += $op->{value};

			#if total number of operations equal desired fragment size,
			#print and reset.
			if ($op_count == $options{fragment}){
				$operations .= $op->{value}."".$op->{type};
				print_fragment($operations, \@data, $frag_count, $start, $prev_delta);

				#adjust delta if last operation is insert or delete
				#since entire opration value was added to current delta
				#remove count that make it grater than 50.
				#delta cannot be greater than 50.
				if ($op->{type} =~ /D/i){
					$curr_delta -= ($op_count-$options{fragment});
				} elsif ($op->{type} =~ /I/i){
					$curr_delta += ($op_count-$options{fragment});
				}

				#rest
				$operations = "";
				$op_count=0;
				$frag_count+=1;
				$start += $options{fragment} + $prev_delta;
				$prev_delta = $curr_delta;
				$curr_delta = 0;

			#if total number of operations is greater than desired fragment size
			#update last operations, print and reset
			} elsif ($op_count > $options{fragment}){
				$operations .= $op->{value} - ($op_count - $options{fragment}) . "$op->{type}";
				print_fragment($operations, \@data, $frag_count, $start, $prev_delta);

				#adjust delta if last operation is insert or delete
				#since entire opration value was added to current delta
				#remove count that make it grater than 50.
				#delta cannot be greater than 50.
				if ($op->{type} =~ /D/i){
					$curr_delta -= ($op_count-$options{fragment});
				} elsif ($op->{type} =~ /I/i){
					$curr_delta += ($op_count-$options{fragment});
				}

				#reset
				$op->{value} = $op_count - $options{fragment};
				$operations = "";
				$op_count = 0;
				$frag_count+=1;
				$start += $options{fragment} + $prev_delta;
				$prev_delta = $curr_delta;
				$curr_delta = 0;
				redo;

			} else {
				$operations .= $op->{value}."".$op->{type};
			}
		}
	}
}

close(FHD);
close(OUT);

#remove temp sam file.
`rm $tmp_sam`;

#############################################################################
sub check_parameters {
    my $options = shift;

	my @required = ("input", "output");

	foreach my $key (@required) {
		unless ($options{$key}) {
			print STDERR "ARG: $key is required\n";
			pod2usage({-exitval => 2,  -message => "error message", -verbose => 1, -output => \*STDERR});
			exit(-1);
		}
	}

	unless($options{fragment}) { $options{fragment} = 50; }
}

#############################################################################
sub print_fragment{
	my ($op, $arr, $num, $start, $delta, $strand) = @_;

	my @tmp_arr = @{$arr};

	$tmp_arr[0] .= "_".($num+1);
	$tmp_arr[5] = $op;
	$tmp_arr[3] += ($start + $delta);
	$tmp_arr[9] = substr($tmp_arr[9],($num*$options{fragment}),$options{fragment});
	$tmp_arr[10] = substr($tmp_arr[10],($num*$options{fragment}),$options{fragment});

	print OUT join("\t",@tmp_arr)."\n";
}
