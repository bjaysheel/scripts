#!/usr/local/biotools/perl/5.14.2/bin/perl

=head1 NAME
   process_qjob.pl

=head1 SYNOPSIS

    USAGE: process_qjob.pl -i=input -f=filter -a=action -e=1

=head1 OPTIONS


B<--input, -i>
	Input qstat xml output file

B<--filter, -f>
	Filter xml file based on value(s) in comma separeated list.
	Currently filter only supported on job name.

B<--action, -a>
	Action to perform on each job. Currently only support delete action (qdel: default)

B<--except, -e>
	Perform ACTION on all jobs EXCEPT those that pass filter or
	Perform ACTION on all jobs that pass filter

B<--help,-h>


=head1 DESCRIPTION
	Given qstat xml file, filter each job based on FILTER values, and perform operation ACION

=head1 INPUT
	qstat -xml output file

=head1 OUTPUT
	Perform action on each job that meets FILTER criteria

=head1 VERSION
	0.0.1

=head1  CONTACT
  bjaysheel@gmail.com


==head1 EXAMPLE
	./process_qjob.pl -i=input.xml -f=SOME_VALUE -a=qdel -e=1

=cut

use strict;
use warnings;
use Data::Dumper;
use Pod::Usage;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use XML::Twig;

my %options = ();
my $results = GetOptions (\%options,
                          'input|i=s',
						  'filter|f=s',
						  'action|a=s',
						  'except|e=s',
						  'prt|p=s',
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

my %hash=();
my $ids = "";

my $twig = XML::Twig->new(twig_handlers =>
							{$options{tag} => \&process_job},
						);

$twig->parsefile($options{input});

$ids =~ s/,,/,/;
$ids =~ s/,$//;

my $cmd = "$options{action} $ids";

print "\n".$cmd."\n";

if (! $options{prt}) {
	system($cmd);
}


exit();

#############################################################################
sub check_parameters {
    my $options = shift;

	my @required = qw(input filter);

	foreach my $key (@required) {
		unless (defined $options{$key}) {
			print STDERR "ARG: $key is required\n";
			pod2usage({-exitval => 2,  -message => "error message", -verbose => 1, -output => \*STDERR});
			exit(-1);
		}
	}

	$options{action} = "qdel";
	$options{except} = 1 unless(defined $options{except});
	$options{tag} = "JB_name" unless (defined $options{tag});
	$options{prt} = 0 unless (defined $options{prt});
}

#############################################################################
sub process_job {
	my ($t, $elt) = @_;

	$options{filter} =~ s/,/|/;

	if ($options{except}) {
		if ($elt->text !~ /$options{filter}/) {
			print "Processing ". $elt->text . "\t";
			print $elt->parent->first_child->text."\t";
			print $elt->parent->{'att'}->{'state'}."\n";

			#### make id list unique
			unless (exists $hash{$elt->parent->first_child->text}) {
				$ids .= $elt->parent->first_child->text.",";
				$hash{$elt->parent->first_child->text} = 1;
			}

		}
	} else {
		if ($elt->text =~ /$options{filter}/) {
			print "Processing ". $elt->text . "\t";
			print $elt->parent->first_child->text."\t";
			print $elt->parent->{'att'}->{'state'}."\n";

			#### make id list unique
			unless (exists $hash{$elt->parent->first_child->text}) {
				$ids .= $elt->parent->first_child->text.",";
				$hash{$elt->parent->first_child->text} = 1;
			}
		}
	}

	return $ids;
}
