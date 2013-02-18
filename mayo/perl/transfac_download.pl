#!/usr/bin/perl

=head1 NAME
   transfac_download.pl

=head1 SYNOPSIS

    USAGE: transfac_download.pl --output output dir [--url base url --range 1..1700]

=head1 OPTIONS

B<--url,-u>
   Base url of transfac pages to download

B<--range,-r>
	transfac matrix identifier number

B<--output, -o>
	Output dir to store html pages.

B<--help,-h>
   This help message

=head1  DESCRIPTION
    Download matrix html pages from transfac for a given range.

=head1  INPUT


=head1  OUTPUT
	tranfac html pages

=head1  CONTACT
  Jaysheel D. Bhavsar @ bjaysheel[at]gmail[dot]com


==head1 EXAMPLE
   transfac_download.pl

=cut

use strict;
use warnings;
use Data::Dumper;
use Pod::Usage;
use LWP;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);

my %options = ();
my $results = GetOptions (\%options,
                          'url|u=s',
						  'range|r=s',
						  'output|o=s',
                          'help|h') || pod2usage();

## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}
#############################################################################
## make sure everything passed was peachy
&check_parameters(\%options);

my $browser = LWP::UserAgent->new();

my $request = HTTP::Request->new( GET => 'https://portal.biobase-international.com/cgi-bin/build_t/idb/1.0/pageview.cgi?view=MatrixReport&matrix_acc=M00795' );
my $response = $browser->request( $request );
print $response->dump()."\n";
exit();

foreach my $val (@{$options{range_arr}}){
	my $response = $browser->get($options{url}."/".$val);

	open (OUT, ">", $options{output}."/".$val.".html") or die "Couldn not open file $!\n";
	binmode(OUT);
	print OUT $response->content;
	close OUT;
}

#############################################################################
sub check_parameters {
    my $options = shift;

	unless ($options{output}) {
		pod2usage({-exitval => 2, -message => "error message", -verbose => 1, -output => \*STDERR});
		exit(1);
	}

	if (defined $options{range}){
		parse_range($options{range});
	} else {
		@{$options{range_arr}} = (1..1700);
	}

	#transfac url.
	#unless($options{url}) { $options{url} = "http://en.wikipedia.org/wiki"; };

	unless($options{url}) { $options{url} = "http://en.wikipedia.org/wiki"; };
	@{$options{range_arr}} = ("Theseus", "Perseus", "Orpheus", "Odysseus", "Meleage");
}

sub parse_range{
	my $range = shift;

	my @top = split(/,/,$range);

	foreach my $i (@top){
		if ($i =~ /\.\./){
			my($start,$stop) = ($i =~ /^(\d+)\.\.(\d+)$/);
			push (@{$options{range_arr}}, ($start..$stop));
		} else { push (@{$options{range_arr}}, $i); }
	}
}
