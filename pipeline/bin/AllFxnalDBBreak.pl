#!/usr/bin/perl -w

=head1 NAME
   libraryHistorgram.pl

=head1 SYNOPSIS

    USAGE: libraryHistogram.pl --server server-name --env dbi [--library libraryId]

=head1 OPTIONS

B<--server,-s>
   Server name from where MGOL blastp records are updated

B<--library,-l>
    Specific libraryId whoes MGOL hit_names to updates

B<--env,-e>
    Specific environment where this script is executed.  Based on these values
    db connection and file locations are set.  Possible values are
    igs, dbi, ageek or test

B<--help,-h>
   This help message

=head1  DESCRIPTION
    Create XML document that contaions information to draw all Fxnal
    count chart. ie: number of hits in keeg, seed, uniref100p, cog, alcame

=head1  INPUT
    The input is defined with --server,  --library.

=head1  OUTPUT
   Updated blastp tables for all/specifed library.

=head1  CONTACT
  Jaysheel D. Bhavsar @ bjaysheel[at]gmail[dot]com


==head1 EXAMPLE
   libraryHistogram.pl --server calliope --env dbi --library 31

=cut

use IO::File;
use strict;
use DBI;
use LIBInfo;
use XML::Writer;
use Pod::Usage;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
BEGIN {
  use Ergatis::Logger;
}

my %options = ();
my $results = GetOptions (\%options,
                          'server|s=s',
                          'library|b=s',
						  'env|e=s',
                          'input|i=s',
						  'outdir|o=s',
                          'log|l=s',
                          'debug|d=s',
                          'help|h') || pod2usage();

## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
                                  'LOG_LEVEL'=>$options{'debug'});
$logger = $logger->get_logger();
#############################################################################
#### DEFINE GLOBAL VAIRABLES.
##############################################################################
my $db_user;
my $db_pass;
my $dbname;
my $db_host;
my $host;

my $dbh0;
my $dbh1;
my $dbh;

my $libinfo = LIBInfo->new();
my $libObject;

my $file_loc = $options{outdir};

## make sure everything passed was peachy
&check_parameters(\%options);
##############################################################################
timer(); #call timer to see when process started.

my $lib_sel = $dbh0->prepare(q{SELECT id FROM library WHERE deleted=0 and server=?});
my $blst_sel = $dbh->prepare(q{SELECT DISTINCT b.sequenceId FROM blastp b RIGHT JOIN sequence s on b.sequenceId=s.id
			       WHERE b.e_value < 0.001 and b.fxn_topHit=1 and b.database_name = ? and s.libraryId=?
					and b.deleted=0 and s.deleted=0
			       ORDER BY sequenceId});

my $rslt = '';
my @libArray;

if ($options{library} <= 0){
    $lib_sel->execute($options{server});
    $rslt = $lib_sel->fetchall_arrayref({});

    foreach my $lib (@$rslt){
	push @libArray, $lib->{'id'};
    }
} else {
    push @libArray, $options{library};
}

foreach my $lib (@libArray){
    print "Processing library $lib\n";

    my $xml_file = "FXNAL_OVERVIEW_XMLDOC_".$lib.".xml";
    my $id_file = "FXNAL_OVERVIEW_IDDOC_".$lib.".xml";
    my $tag_num = 1;

    my $xml_out = new IO::File(">".$file_loc."/xDocs/".$xml_file)
		    or die "Could not open file ".$file_loc."/xDocs/".$xml_file." to write\n";
    my $id_out = new IO::File(">".$file_loc."/xDocs/".$id_file)
		    or die "Could not open file ".$file_loc."/xDocs/".$id_file." to write\n";

    my $xml_writer = new XML::Writer(OUTPUT=>$xml_out);
    my $id_writer = new XML::Writer(OUTPUT=>$id_out);

    $xml_writer->xmlDecl("UTF-8");
    $id_writer->xmlDecl("UTF-8");

    $xml_writer->startTag("root");
    $id_writer->startTag("root");

    ($xml_writer, $id_writer, $tag_num) = getNodeFor('ACLAME',$lib,$xml_writer,$id_writer,$id_file,$tag_num);
    ($xml_writer, $id_writer, $tag_num) = getNodeFor('COG',$lib,$xml_writer,$id_writer,$id_file,$tag_num);
    ($xml_writer, $id_writer, $tag_num) = getNodeFor('UNIREF100P',$lib,$xml_writer,$id_writer,$id_file,$tag_num);
    ($xml_writer, $id_writer, $tag_num) = getNodeFor('KEGG',$lib,$xml_writer,$id_writer,$id_file,$tag_num);
    ($xml_writer, $id_writer, $tag_num) = getNodeFor('SEED',$lib,$xml_writer,$id_writer,$id_file,$tag_num);

    $xml_writer->endTag("root");
    $id_writer->endTag("root");

    $xml_writer->end();
    $id_writer->end();

    $xml_out->close();
    $id_out->close();
}

$dbh->disconnect;
$dbh0->disconnect;
$dbh1->disconnect;

timer(); #call timer to see when process ended.
exit(0);

###############################################################################
####  SUBS
###############################################################################
sub check_parameters {
    my $options = shift;

    my $flag = 0;

    # if library list file or library file has been specified
    # get library info. server, id and library name.
    if ((defined $options{input}) && (length($options{input}))){
      $libObject = $libinfo->getLibFileInfo($options{input});
      $flag = 1;
    }

    # if server is not specifed and library file is not specifed show error
    if (!$options{server} && !$flag){
      pod2usage({-exitval => 2,  -message => "error message", -verbose => 1, -output => \*STDERR});
      exit(-1);
    }

    # if exec env is not specified show error
    unless ($options{env}) {
      pod2usage({-exitval => 2,  -message => "error message", -verbose => 1, -output => \*STDERR});
      exit(-1);
    }

    # if no library info set library to -1;
    unless ($options{library}){
        $options{library} = -1;
    }

    # if getting info from library file set server and library info.
    if ($flag){
        $options{library} = $libObject->{id};
        $options{server} = $libObject->{server};
    }

    if ($options{env} eq 'dbi'){
	  $db_user = q|bhavsar|;
	  $db_pass = q|P3^seus|;
	  $dbname = q|VIROME|;
	  $db_host = $options{server}.q|.dbi.udel.edu|;
	  $host = q|virome.dbi.udel.edu|;
    }elsif ($options{env} eq 'camera') {
      $db_user = q|virome_app|;
          $db_pass = q|camera123|;
          $dbname = q|virome_stage|;
          $db_host = q|dory.internal.crbs.ucsd.edu|;
          $host = q|dory.internal.crbs.ucsd.edu|;
   }elsif ($options{env} eq 'igs'){
	  $db_user = q|dnasko|;
	  $db_pass = q|dnas_76|;
	  $dbname = q|virome_processing|;
	  $db_host = q|dnode001.igs.umaryland.edu|;
	  $host = q|dnode001.igs.umaryland.edu|;
    }elsif ($options{env} eq 'ageek') {
	  $db_user = q|bhavsar|;
	  $db_pass = q|Application99|;
	  $dbname = $options{server};
	  $db_host = q|10.254.0.1|;
	  $host = q|10.254.0.1|;
    }else {
	  $db_user = q|kingquattro|;
	  $db_pass = q|Un!c0rn|;
	  $dbname = q|VIROME|;
	  $db_host = q|localhost|;
	  $host = q|localhost|;
    }

    $dbh0 = DBI->connect("DBI:mysql:database=virome_stage;host=$host",
	"$db_user", "$db_pass",{PrintError=>1, RaiseError =>1, AutoCommit =>1});

    $dbh1 = DBI->connect("DBI:mysql:database=uniref_lookup2;host=$host",
	"$db_user", "$db_pass",{PrintError=>1, RaiseError =>1, AutoCommit =>1});

    $dbh = DBI->connect("DBI:mysql:database=$dbname;host=$db_host",
	"$db_user", "$db_pass",{PrintError=>1, RaiseError =>1, AutoCommit =>1});
}

###############################################################################
sub timer {
    my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
    my @weekDays = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
    my ($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
    my $year = 1900 + $yearOffset;
    my $theTime = "$hour:$minute:$second, $weekDays[$dayOfWeek] $months[$month] $dayOfMonth, $year";
    print "Time now: " . $theTime."\n";
}

###############################################################################
sub getNodeFor {
    my ($db, $lib, $xw, $iw, $fname, $tag) = @_;

    $blst_sel->execute($db,$lib);
    my $rslt = $blst_sel->fetchall_arrayref({});

    if (@{$rslt} > 0){
	my $idList = '';
	my $count = 0+@{$rslt};

	foreach my $row (@$rslt) {
	    $idList .= $row->{sequenceId} . ",";
	}

	$idList =~ s/,$//;

	if ($db eq 'UNIREF100P'){
	    $db = 'GO';
	}

	$xw->emptyTag("CATEGORY", 'LABEL'=>$db, 'VALUE'=>$count, 'TAG'=>'TAG_'.$tag, 'IDFNAME'=>$fname);
	$iw->emptyTag("TAG_".$tag, 'IDLIST'=>$idList);
    }

    return $xw, $iw, $tag+1;
}
