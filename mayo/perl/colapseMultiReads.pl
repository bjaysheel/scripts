#!/usr/local/biotools/perl/5.14.2/bin/perl

=head1 NAME
   colapseMultiReads.pl

=head1 SYNOPSIS
    USAGE: colapseMultiReads.pl -i=input file -o=output dir

=head1 OPTIONS

B<--input, -i>
	Sam/Bam file

B<--output_dir, -o>
	output dir

B<--help,-h>


=head1  DESCRIPTION
	Colapse multi reads into one

=head1  INPUT

=head1  OUTPUT


=head1  CONTACT
  bjaysheel@gmail.com


==head1 EXAMPLE
	./colapseMultiReads.pl -i=input_file.bam -o=output_dir

=cut

use strict;
use warnings;
use Data::Dumper;
use Pod::Usage;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);

my %options = ();
my $results = GetOptions (\%options,
                          'input|i=s',
						  'output_dir|o=s',
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

my $cmd="";
my $name="";

if ($options{input} =~ /\.bam/){
	$name = $options{input};
	$name =~ s/bam/sam/;
	$name = $options{output_dir}."/".$name;

	## get sam file
	$cmd = "$options{samtools}/samtools view $options{input} > $name";
	execute_cmd($cmd);
} else {

}

my $name = $options{in}



## update unmapped sequence id and remove ZT:X:X tag
$cmd = "$options{samtools}/samtools view $options{unmapped}";
$cmd .= " | sed -e 's/\\/[0-9]//p' | sed -e 's/ZT:[a-zA-Z]:[a-zA-Z]//p'";
$cmd .= " > unmapped.edited.sam";
execute_cmd($cmd);

## get accepted hits header
$cmd = "$options{samtools}/samtools view -H $options{accepted} > $options{output_dir}/accepted.header.sam";
execute_cmd($cmd);

## sort accepted hits by id.
$cmd = "/usr/java/latest/bin/java -Xmx6g -Xms512m";
$cmd .= " -jar $options{picard}/SortSam.jar";
$cmd .= " INPUT=$options{accepted}";
$cmd .= " OUTPUT=$options{output_dir}/accepted_hits-sorted.id.bam";
$cmd .= " SO=queryname MAX_RECORDS_IN_RAM=1000000";
$cmd .= " TMP_DIR=$options{output_dir}/tmp/ VALIDATION_STRINGENCY=SILENT";
execute_cmd($cmd);

## get only unique non fusion reads
$cmd = "$options{samtools}/samtools view $options{output_dir}/accepted_hits-sorted.id.bam";
$cmd .= " | awk -F '\\t' '{ for(i=12; i<=NF; i++){ if (\$i ~ \"NH:i:1\$\"){print}} }'";
$cmd .= " | awk '{ if (\$0 !~ \"XF:Z\") {print} }'";
$cmd .= " > $options{output_dir}/accepted_hits.unique.nofusion.sam";
execute_cmd($cmd);

## get fusion reads
$cmd = "$options{samtools}/samtools view $options{output_dir}/accepted_hits-sorted.id.bam";
$cmd .= " | awk '{ if (\$0 ~ \"XF:Z\") {print} }'";
$cmd .= " > $options{output_dir}/accepted_hits.fusion.sam";
execute_cmd($cmd);

## get multi reads
$cmd = "$options{samtools}/samtools view $options{output_dir}/accepted_hits-sorted.id.bam";
$cmd .= " | awk '{ if (\$0 !~ \"NH:i:1\"){print} }'";
$cmd .= " > $options{output_dir}/accepted_hits.multiread.sam";
execute_cmd($cmd);

exit();

#############################################################################
sub check_parameters {
    my $options = shift;

	my @required = qw(input output_dir);

	foreach my $key (@required) {
		unless ($options{$key}) {
			print STDERR "ARG: $key is required\n";
			pod2usage({-exitval => 2,  -message => "error message", -verbose => 1, -output => \*STDERR});
			exit(-1);
		}
	}

	$options{'samtools'} = "/projects/bsi/bictools/apps/alignment/samtools/latest" unless ($options{'samtools'});
	$options{'picard'} = "/projects/bsi/bictools/apps/alignment/picard/latest" unless ($options{'picard'});
}

#############################################################################
sub execute_cmd {
	my $cmd = shift;

	print STDOUT $cmd."\n";
	system($cmd);
}

#############################################################################
sub create_dir_struct {
	my $options = shift;

	my $dir = "$options{output_dir}/tmp";
	if ( -d $dir) {
		print STDERR "Directory $dir exist";
	} else {
		execute_cmd("mkdir -p $dir");
	}
}
