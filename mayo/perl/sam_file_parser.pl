#!/usr/bin/perl

=head1 NAME
   sam_file_parser.pl

=head1 SYNOPSIS

    USAGE: sam_file_parser.pl --input input_sam_file --output output_sam_file [--threshold 10000 --neighbour 10 --minMatch 4]

=head1 OPTIONS

B<--input,-i>
   Input file

B<--output,-o>
   output file

B<--threshold,-t>
    Map threshold to match with history

B<--neighbour,-n>
    Neighbours to search within

B<--minMatch, -m>
	Min number of matches within neighbour

B<--help,-h>
   This help message

=head1  DESCRIPTION
    Parse sam file and only keep sequencing reads that maps within THRESHOLD of at
	least on sequence withing the HISTORY

=head1  INPUT


=head1  OUTPUT
	A SAM file where each sequence meet the THRESHOLD specified within given
	NEIGHBOURHOOD

=head1  CONTACT
  Jaysheel D. Bhavsar @ bjaysheel[at]gmail[dot]com


==head1 EXAMPLE
   sam_file_parser.pl --input filename.sam --threshold 10000 --neighbour 9

=cut

use strict;
use warnings;
use Data::Dumper;
use DBI;
use Pod::Usage;
use Scalar::Util qw(looks_like_number);
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);

my %options = ();
my $results = GetOptions (\%options,
                          'input|i=s',
						  'output|o=s',
                          'threshold|t=s',
						  'neighbour|n=s',
						  'count|c=s',
						  'help|h') || pod2usage();

## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}
#############################################################################
## make sure everything passed was peachy
&check_parameters(\%options);

my @history=();
my @seen=();
my $x = 0;

open (FHD, "<", $options{input}) or die "Cound not open file $options{input}\n";
open (OUT, ">", $options{output}) or die "Cound not open file $options{output}\n";

while (<FHD>){
	chomp $_;

	#skip processing lines starting with @ just print to output file.
	if ($_ =~ /^@/){
		print OUT $_."\n";
		next;
	}

	check_sequence($_);

	$x++;

	#exit if ($x > 81);
}

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

	unless($options{threshold}) { $options{threshold} = 10000; }
	unless($options{neighbour}) { $options{neighbour} = 10; }
	unless($options{minMatch}) { $options{minMatch} = 4; }
}

sub check_sequence {
	my $line = shift;

	my @data = split(/\t/,$line);
	my $count = 0;

	#look at above neighbours to see if threshold is meet
	#if so print to output and mark it to avoid duplicate entry
	#in output file.
	foreach my $h (@history){
		my @hline = split(/\t/,$h);

		#threshold check
		if ((abs($hline[7]-$data[7]) < $options{threshold}) && ($hline[6] eq $data[6])){

			$count++;

			#check at neighbours below if threhold meets and haven't been
			#printed ot output file print it and mark it
			if (seen($hline[0], $hline[1],0)){
				print OUT $h."\n";
			}
		}
	}


	if ($count >= $options{minMatch}){
		#print current line if threshold meet and mark it.
		print OUT $line ."\n";
	}
	seen($data[0],$data[1],$count);

	# keep history of lines upto threshold
	# this is a first in first out queue.
	if (scalar (@history) >= $options{neighbour}){
		pop @history;
	}
	unshift @history, $line;
}

#check if this sequence has already passed and has been printed.
#if seen and printed return true, else add to seen as it will be
#printed and return false.
sub seen {
	my ($name, $pos, $c) = @_;

	#check if we have already see this sequence.
	foreach my $s (@seen){
		if (($s->{name} =~ /$name/) && ($s->{pos} =~ /$pos/) && ($s->{count} >= $options{minMatch})){
			return 1;
		} elsif (($s->{name} eq /$name/) && ($s->{pos} =~ /$pos/) && ($s->{count} < $options{minMatch})){
			$s->{count}++;
			return 0;
		}
	}

	#seen array only needs to be as large as history
	#because we are not going to be comparing elments outside
	#the history window.
	if (scalar (@seen) >= $options{neighbour}){
		pop @seen;
	}
	unshift @seen, {name=>$name, pos=>$pos, count=>$c};

	return 0;
}
