package MyBlat;

use strict;
use warnings;
use Data::Dumper;

sub new {
    my ($class) = @_;
    my $self = {};

    $self->{blat_ref} = "";
	$self->{blat_path} = "";

 	bless($self,$class);
	$self->init();

    return $self;
}

sub init {
	my $self = shift;

	$self->{blat_ref} = "/data2/bsi/reference/db/hg19.2bit";
	$self->{blat_lookup} = "/data2/bsi/staff_analysis/b.lookup/primer/lookup.tbx";
	$self->{blat_path} = "/projects/bsi/bictools/apps/alignment/blat/34_x64";
}

sub execute {
	my $self = shift;
	my ($sequence, $output_dir) = @_;

	my $tmp_fsa = "tmp_seq.fsa";
	my $tmp_psl = "tmp_seq.psl";

	#write FASTA file to BLAT
	open (TMP, ">", "$output_dir/$tmp_fsa") or die "Could not write temp sequence file\n$!\n";
	print TMP ">Blat_seq\n";
	print TMP "$sequence\n";
	close(TMP);

	#execute blat search.
	`$self->{blat_path}/blat $self->{blat_ref} $output_dir/$tmp_fsa -noHead $output_dir/$tmp_psl`;

	if (($? >> 8) < 0){
		print STDERR "ERROR: Problem executing blat\n$!\n";
		exit(-1);
	}
}

sub numOfHits{
	my $self = shift;
	my ($identity, $output_dir) = @_;

	my $tmp_psl = "tmp_seq.psl";
	my $hit_count = 0;

	#read BLAT results.
	open (PSL, "<", "$output_dir/$tmp_psl") or die "Could not open psl file\n$!\n";

	while (<PSL>){
		chomp $_;
		my @data = split(/\t/, $_);

		if ((($data[0]-$data[7])/$data[10])*100 >= $identity){
			$hit_count++;
		}
	}
	close(PSL);

	return $hit_count;
}

sub clean {
	my $self = shift;
	my ($output_dir) = @_;

	`rm -rf $output_dir/tmp_seq.fsa`;
	`rm -rf $output_dir/tmp_seq.psl`;
}

1;
