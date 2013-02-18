#!/usr/local/biotools/perl/5.10.0/bin/perl

##
# Usage: SoftSearch.pl [-lqrmsd] -b <BAM> -f <Genome.fa> -sam <samtools path> -bed <bedtools path>
# Created 1-30-2012 by Steven Hart, PhD
# hart.steven@mayo.edu
# Required bedtools & samtools to be in path

use Getopt::Long;
use warnings;
use strict;

my ($INPUT_BAM,$INPUT_FASTA,$OUTPUT_FILE,$minSoft,$minSoftReads,$dist_To_Soft,$bedtools,$samtools);

my ($minRP, $temp_output, $num_sd, $MapQ, $chrom, $unmated_pairs, $minBQ, $pair_only, $soft_only, $command);

#Declare variables
GetOptions(
	'b=s' => \$INPUT_BAM,
	'f=s' => \$INPUT_FASTA,
	'o:s' => \$OUTPUT_FILE,
	'm:i' => \$minRP,
	'l:i' => \$minSoft,
	'r:i' => \$minSoftReads,
	't:i' => \$temp_output,
	's:s' => \$num_sd,
	'd:i' => \$dist_To_Soft,
	'q:i' => \$MapQ,
	'c:s' => \$chrom,
	'u:s' => \$unmated_pairs,
	'x:s' => \$minBQ,
	'p:s' => \$pair_only,
	'g:s' => \$soft_only,
	"help|h|?|"	=> \&usage);


#if(defined($INPUT_BAM)){$INPUT_BAM=$INPUT_BAM} else {print usage(); die "Where is the BAM file?\n\n" ;}
#if(defined($INPUT_FASTA)){$INPUT_FASTA=$INPUT_FASTA} else {print usage(); die "Where is the fasta file?\n\n"; }

unless (defined($INPUT_BAM)) {
	usage();
	die "Where is the BAM file?\n\n";
}

unless (defined($INPUT_FASTA)){
	usage();
	die "Where is the fasta file?\n\n";
}

my $index = `ls $INPUT_BAM.bai`;
if($index eq "") {
	die "\n\nERROR: you need index your BAM file\n\n"
}

#### get current time
print "Start Time : " . spGetCurDateTime() . "\n";
my $now = time;

#if(defined($OUTPUT_FILE)) {$OUTPUT_FILE=$OUTPUT_FILE} else {$OUTPUT_FILE="output.vcf"; print "\nNo outfile specified.  Using output.vcf as default\n\n"; }
#if(defined($minSoft)){$minSoft=$minSoft} else {$minSoft=5; }
#if(defined($minRP)){$minRP=$minRP} else {$minRP=5; }
#if(defined($minSoftReads)){$minSoftReads=$minSoftReads} else {$minSoftReads=5; }
#if(defined($dist_To_Soft)){$dist_To_Soft=$dist_To_Soft} else {$dist_To_Soft=300; }
#if(defined($num_sd)){$num_sd=$num_sd} else {$num_sd=4; }
#if(defined($MapQ)){$MapQ=$MapQ; } else {$MapQ=35; }
#if(defined($chrom)){$chrom=join("","chr",$chrom); } else {$chrom=""; }
#if(!defined($unmated_pairs)){$unmated_pairs=0; }

unless (defined($OUTPUT_FILE)){
	$OUTPUT_FILE = "output.vcf";
	print "\nNo outfile specified.  Using output.vcf as default\n\n";
}

unless (defined($minSoft)) { $minSoft = 5; }
unless (defined($minRP)) { $minRP = 5; }
unless (defined($minSoftReads)) { $minSoftReads = 5; }
unless (defined($dist_To_Soft)) { $dist_To_Soft = 300; }
unless (defined($num_sd)) { $num_sd = 4; }
unless (defined($MapQ)) { $MapQ = 35; }
if (defined($chrom)) {
	$chrom = join("","chr",$chrom);
} else { $chrom=""; }
unless (defined($unmated_pairs)) { $unmated_pairs=0; }

my $badQualValue = chr($MapQ);
if(defined($minBQ)) { $badQualValue = chr($minBQ); }
if($badQualValue eq "#") { $badQualValue = "\#"; }

# adding and cheking for samtools and bedtools in the PATh
## check for bedtools and samtools in the path
$bedtools=`which intersectBed` ;
unless (defined($bedtools) && length($bedtools)) {
	die "\nError:\n\tno bedtools. Please install a bedtools and add to the path\n";
}

#$samtools=`samtools 2>&1`;
$samtools = `which samtools`;
print "$samtools"."\n";
#unless ($samtools !~ /(samtools)/i) {
unless (defined($samtools) && length($samtools)) {
	die "\nError:\n\tno samtools. Please install a samtools and add to the path\n";
}

print "Usage = SoftSearch.pl -l $minSoft -q $MapQ -r $minSoftReads -d $dist_To_Soft -m $minRP -s $num_sd -c $chrom -b $INPUT_BAM -f $INPUT_FASTA -o $OUTPUT_FILE \n\n";

sub usage {
	print "\nusage: SoftSearch.pl [-cqlrmsd] -b <BAM> -f <Genome.fa> \n";
	print "-q\tMinimum mapping quality [35]\n";
	print "-l\tMinimum length of soft-clipped segment [5]\n";
	print "-r\tMinimum depth of soft-clipped reads at position [5]\n";
	print "-m\tMinimum number of discordant read pairs [5]\n";
	print "-s\tNumber of sd away from mean to be considered discordant [4]\n";
	print "-u\tNumber of unmated pairs [0]\n";
	print "-d\tMinimum distance between soft-clipped segments and \n\tdiscordant read pairs [300]\n";
	print "-o\tOutput file name [output.vcf]\n";
	print "-t\tPrint temp files for debugging [no|yes]\n";
    print "-c\tuse only this chrom\n";
	print "-p\tuse paired-end mode only [no|yes]\n";
	print "-g\tdisable paired-only seach [no|yes]\n\n";
	exit(1);
}


#############################################################
# create temporary variable name
#############################################################
srand (time ^ $$ ^ unpack "%L*", `ps axww | gzip -f`);
our $random_name=join "", map { ("a".."z")[rand 26] } 1..8;

#############################################################
# Calcualte insert size distribution of properly mated reads
#############################################################
my $index_input_bam= $INPUT_BAM  . ".bai";
if (! -e $index_input_bam)	{
	my $command=join("","samtools index ${INPUT_BAM}");
	print "Indexing the input BAM file\n";
	print "$command\n\n";
	system("$command");
}
my $metrics=`samtools view -q $MapQ -f2 $INPUT_BAM $chrom|cut -f9|head -10000|awk '{if (\$1<0){\$1=-\$1}else{\$1=\$1} sum+=\$1; sumsq+=\$1*\$1} END {print sum/NR, sqrt(sumsq/NR - (sum/NR)**2)}'`;
my ($mean,$stdev)=split(/ /,$metrics);
$stdev=~s/\n//;
my $upper_limit=int($mean+($num_sd*$stdev));
my $lower_limit=int($mean-($num_sd*$stdev));
print "The mean insert size is $mean +/- $stdev (sd)\n
The upper limit = $upper_limit\n
The lower limit = $lower_limit\n";


my $tmp_name = join ("",$random_name,".tmp.bam");

#############################################################
# Make sam file that has soft clipped reads
#############################################################
#give file a name
my $random_file_sc = "";
if(!defined($pair_only)){
	$random_file_sc = join ("",$random_name,".sc.sam");
	$command = join ("","samtools view -q $MapQ -f 1 -F 1036  $INPUT_BAM $chrom |awk '{OFS=\"\\t\"}(\$6~/S/)'>",$random_file_sc);

	print "Making SAM file of soft-clipped reads\n\n";
	print "\n\n$command\n\n";
	system("$command");

	#############################################################
	# Find areas that have deep enough soft-clip coverage
	print "Identifying soft-clipped regions that are at least $minSoft bp long and Base Qual is not equal to $badQualValue\n\n";
	open (FILE,"$random_file_sc")||die "Can't open soft-clipped sam file\n";
	my $tmpfile = join("",$random_file_sc,".sc.passfilter");
	open (OUT,">$tmpfile")||die "Can't write files here!\n";
	while(<FILE>){
		@_=split(/\s/,$_);
		my @CIGAR = split(/(S|M|I|D|N|H|X|P)/, $_[5]);
		## find the positions of S (soft clipping field) field
		my @softclip_pos;
		for (my $i=0;$i<=$#CIGAR;$i++)	{
			if ($CIGAR[$i] eq 'S')	{
				push (@softclip_pos,$i-1);
			}
		}

		for (my $i=0;$i<=$#softclip_pos;$i++) {
			if	($CIGAR[$softclip_pos[$i]] > $minSoft){
				###############Make sure base qualities don't have more than 2 bad marks
				my $qual = $_[10];
				my $TR = ($qual =~ tr/$badQualValue//);
				if($badQualValue eq "#") { $TR = ($qual=~tr/\#//); }
				#Skip the soft clip if there is more than 2 bad qual values
				next if($TR > 2);
				print OUT;
				last;
			}
		}
	}

	close FILE;
	close OUT;
	$command=join(" ","mv",$tmpfile,$random_file_sc);
	print "$command\n\n";
	system("$command");
}

#############################################################
#  Make discordant read BAM files
#############################################################
#make bam file that has discordant read paird
print "Making BAM file that has discordant read pairs\n\n";
my $random_file_disc=join("",$random_name,".disc.bam");
$command=join ("","samtools view -q $MapQ -h -f 1 -F 1024 $INPUT_BAM $chrom | awk '{OFS=\"\\t\"}{if (\$1~/^@/){print \$0}else{if(\$7 ~/=/){\$7=\$3};if(\$9<0){\$9=-\$9}if(((\$9<$lower_limit)||(\$9>$upper_limit))&&((\$7~/chr/)&&(\$3~/chr/))){print \$0}}}'| samtools view -hSb - ",">",$random_file_disc);
print "$command\n";
system("$command");
unlink("$tmp_name","$index");
print "Indexing BAM file\n\n";
$command=join (" ","samtools index ",$random_file_disc);
print "$command\n\n";

system("$command");
my $index2=join("",$random_file_disc,".bai");


#########################################################
#Stack up Sof$random_file"tClips
#########################################################
my $random_file=join("",$random_name,".sc.direction");
if(!defined($pair_only)){
	open (FILE,"$random_file_sc")|| die "Can't open sam file\n";
#	$random_file=join("",$random_name,".sc.direction");
	print "Calling sides of soft-clips\n\nTMPOUT=$random_file\tINPUT=$random_file_sc\n\n";
	open (TMPOUT,">$random_file")|| die "Can't create tmp file\n";
	while (<FILE>){
		@_=split(/\s/,$_);
		my @CIGAR=split(/(S|M|I|D|N|H|X|P)/, $_[5]);
		my @softclip_pos;

		for (my $i=0; $i<=$#CIGAR; $i++) {
			if ($CIGAR[$i] eq "S") {
				push (@softclip_pos, $i-1)
			}
		}

		#If there are two parts that are soft clipped, set the lower one to 0
		if($#softclip_pos > 1) {
			if ($CIGAR[$softclip_pos[0]]<$CIGAR[$softclip_pos[1]]) {
				$CIGAR[$softclip_pos[0]]=0
			} else {
				$CIGAR[$softclip_pos[1]]=0
			}
		}

		#If the soft clip occurs at beginning of read and its on the plus strand, then it's a left clip
        if(($CIGAR[$softclip_pos[0]]>$minSoft)&&($_[5]!~/^.*S$/)){
            my $softBases=substr($_[9],0,$CIGAR[$softclip_pos[0]]);
            my $end1=$_[3]-1;
			my $end2=$end1+1;

			if(($_[1]==97)||($_[1]==99)||($_[1]==65)||($_[1]==163)||($_[1]==161)||($_[1]==129)){
                print  TMPOUT "$_[2]\t$end1\t$end2\t$softBases|left\n"; #chr breakpoint num_of_bases_clipped\
			} else{
                $end1=$_[3]-1; $end2=$end1+1;
                print  TMPOUT "$_[2]\t$end1\t$end2\t$softBases|left\n"; #chr breakpoint num_of_bases_clipped\
            }
		}
		#If the soft clip occurs at end of read and its on the minus strand, then it's a left clip
        else {

			unless (defined $softclip_pos[1]) {
				print "File being read:: $random_file\n";
				print "Line being processed:: $_\n\n";
				print Dumper(\@CIGAR);
				print "\n\n";
				print Dumper(\@softclip_pos);

				exit();
			}


            my $softBases=substr($_[9],0,$CIGAR[$softclip_pos[1]]);
            if(($_[1]==83)||($_[1]==81)||($_[1]==177)||($_[1]==147)||($_[1]==145)||($_[1]==113)){
                my $end1=$_[3]+length($softBases)-1;
				my $end2=$end1+1;
                print TMPOUT "$_[2]\t$end1\t$end2\t$softBases|right\n"; #chr breakpoint num_of_bases_clipped
            } else {
                my $end1=$_[3]+length($softBases)-1;
				my $end2=$end1+1;
                print TMPOUT "$_[2]\t$end1\t$end2\t$softBases|right\n"; #chr breakpoint num_of_bases_clipped
			}
		}
	}

close FILE;
close TMPOUT;
}
#
#
#
#
######################################################
# Transform Read pair groups into softclip equivalents
######################################################
#
#
#
if(defined($soft_only)){
	use FindBin qw($Bin);
	my $path=$Bin;
	print"\n\nPATH=$path\n\n";
	my $tmp_out=join("",$random_name,"RP.out");
	$command=join("","perl ",$path,"/Bam2pair.pl -b ",$random_file_disc," -o ",$tmp_out," -isize 100 -winsize ",$upper_limit," -min ",$minRP);
	print "$command\n";
	system("$command");
	$command=join("","perl -ane '\$end1=\@F[1];\$end2=\@F[3];print join(\"\\t\",\@F[0..1],\$end1,\"unknown|left\");print \"\\n\";print join(\"\\t\",\@F[2..3],\$end2,\"unknown|left\");print \"\\n\"' ", $tmp_out," >> ",$random_file);
	print "$command\n";
	system($command);
	unlink($tmp_out);
}
#
#
#
######################################################
# EXPERIMENTAL
######################################################
#
#
#

#die "look at $random_file for left and right clips\n";
#######################################################
unlink("$random_file","$tmp_name","$random_file","$random_file_sc","$random_file_disc","$index","$index2","$random_name") if (-z $random_file || ! -e $random_file) ;
die "Softclipped file is empty($random_file).\nNo soft clipping found using desired paramters\n\n" if (-z $random_file || ! -e $random_file) ;

#############################################################
#  Make sure there are enough soft-clippped supporting reads
#############################################################
my $outfile=join("",$random_file,".sc.mergeBed");
#sortbed -i .sc.direction | mergeBed -nms -d 25 -i stdin > .sc.mergeBed
$command=join(" ","sortBed -i",$random_file," | mergeBed  -nms -i stdin","|grep \";\"","|awk '{OFS=\"\t\"}(NF==4)'",">",$outfile);
print "$command\n";
system("$command");
unlink("$tmp_name","$random_file","$outfile","$random_file_sc","$random_file_disc","$index","$index2","$random_name") if (-z $outfile || !-e $outfile) ;
die "mergeBed file is empty.\nNo strucutral variants found\n\n" if (-z $outfile || ! -e $outfile) ;
print "completed mergeBed\n\n";

###############################################################
# If left and right are on the same line, make into 2 lines
###############################################################
open (INFILE,$outfile)||die "couldn't open temp file : $. \n\n";
my $tmp2=join("",$random_name,".sc.fixed.mergeBed");
print "INFILE=$outfile\tOUTFILE=$tmp2\n\n";
#INPUT FORMAT=chr9\t131467\t131473\tATGCTTATTAAAA|left;TTATTAAAAGCATA|left
open (OUTFILE,">$tmp2")||die "couldn't create temp file : $. \n\n";
while(<INFILE>){
	my $l=$_;chomp $l;
	my @a=split(/\t/,$l);
	my $inf = $a[3];
	my @info = split(/\;/, $inf);
	my @left=();
	my @right=();
	@left=grep(/left/, @info);
	@right=grep(/right/, @info);
	#New
	my $lft = join(";",@left);
	my $rht = join(";",@right);
	$inf = join(";",@info);
	if((@left)&&(@right)){
		print OUTFILE "$a[0]\t$a[1]\t$a[2]\t$lft\n$a[0]\t$a[1]\t$a[2]\t$rht\n";
	} else {
		my $all=join("\t",@a[0..2],$inf);
		print OUTFILE "$all\n";
	}
}
# make sure output file name is $outfile
$command=join(" ","sed -e '/ /s//\t/g'", $tmp2,"|awk 'BEGIN{OFS=\"\\t\"}(NF==4)'", "|perl -pne 's/ /\t/g'>",$outfile);
system("$command");
print "$command\n";
unlink("$tmp_name","$random_file","$tmp2","$outfile","$random_file_sc","$random_file_disc","$index","random_name") if (-z $outfile || ! -e $outfile) ;
die "Fixed mergeBed file is empty($outfile).\nNo strucutral variants found\n\n" if (-z $outfile || ! -e $outfile) ;
print "completed fixing mergeBed\n\n";

###############################################################
# Seperate directions of soft clips
###############################################################
my $left_sc=join("","left",$tmp2);
my $right_sc=join("","right",$tmp2);
$command=join(" ","grep left ", $tmp2, " |sed -e '/left /s//left\;/g' | sed -e '/ /s//\t/g' >",$left_sc);
system("$command");
print "$command\n";
$command=join(" ","grep right ", $tmp2, " |sed -e '/right /s//right\;/g' | sed -e '/ /s//\t/g' >",$right_sc);
system("$command");
print "$command\n";

###############################################################
# Count the number and identify directions of soft clips
###############################################################
print "\nCount the number and identify directions of soft clips\nINFILE=$outfile\n\n";
print "looking in $outfile\n";
$outfile=join("",$random_name,".sc.fixed.mergeBed");

open (INFILE,$outfile)||die "couldn't open temp file\n\n";
my $tmp3=join("",$random_file,"predSV");
open (OUTFILE,">$tmp3")||die "couldn't create temp file\n\n";
while(<INFILE>){
	@_=split(/\t/,$_);
	my $count=tr/\;//;
	$count=$count+1;
	my $left=0;
	my $right=0;
	while ($_ =~ /left/g) { $left++ } # count number of right clips
	while ($_ =~ /right/g) { $right++ } # count number of left clips
	###############################################################
	if ($count >= $minSoftReads){
		my @clips=split(/\;|\|/,$_[3]); #get longets soft-clipped read
		my $max=0;
		my $temp=0;
		my $temp2=0;
		my $temp3=0;
		my $dir=0;
		my $maxSclip=0;
		for (my $i=0;$i<$count;$i++) {
			my $plus1=$i+1;
			$temp=length($clips[$i]);
			$temp2=$clips[$plus1];
			$temp3=$clips[$i];
			if ($temp > $max){$maxSclip=$temp3;$max =$temp;$dir=$temp2}
			else {$max=$max;$dir=$dir;$maxSclip=$maxSclip}
			$i++;
		}
		my $order2=join("|",$left,$right);
		print OUTFILE join ("\t",@_[0..2],$maxSclip,$max,$dir,$count,$order2) . "\n";
	} elsif($_=~/unknown/) {
		print OUTFILE join ("\t",@_[0..2],"NA","NA","left","NA","NA|NA") . "\n";
        print OUTFILE join ("\t",@_[0..2],"NA","NA","right","NA","NA|NA") . "\n";
	}
#Format is Chrom,start, end,longest Soft-clip,length of longest Soft-clip, direction of longest soft-clip,#supporting softclips,#right Sclips|#left Sclips
}
close INFILE;
close OUTFILE;

unlink("$tmp_name","$random_file","$tmp3","$outfile","$random_file_sc","$random_file_disc","$index","$random_name","$right_sc","$left_sc") if (-z $tmp3 || !-e $tmp3) ;
die "No structural variants found while Counting the number and identify directions of soft clips.\n" if (-z $tmp3 || !-e $tmp3) ;

###############################################################
# Print header information
###############################################################
open (OUT,">$OUTPUT_FILE")||die "Can't write files here!\n";
&print_header();
close OUT;

print "\n random_file=$random_file\n \n tmp3=$tmp3\n  outfile=$outfile\n  random_file_sc=$random_file_sc\n  $random_file_disc\n  $index\n  random_name=$random_name\n ";
###############################################################
# DO the bulk of the work
###############################################################
open (FILE,"$tmp3")|| die "Can't open file\n";
open (OUT,">>$OUTPUT_FILE")|| die "Can't open file\n";

my $tmp_hits=join("",$random_file,"hits");

print "\nusing $tmp3 and writing to $OUTPUT_FILE \n";
while (<FILE>){
	#$_=~s/\n//;
	@_=split(/\t/,$_);
	#If left clip {+- or -- or -+ }{+- are uninformative b/c they go upstream}
	#If right clip {++ or -- or +-}
	if($_[5] eq "left"){
		#print "LEFT###########################################\n$_\n";
		#my $plus_Reads=$_;
		my @plus_Reads=split(/\s/,$_);
		$plus_Reads[7]=~s/\n//g;
		#Get all types of compatible reads
		#Get improperly paired reads (@ max distance)
		my $start=$_[1];
		my $end=$_[2]+$dist_To_Soft;
		my $target=join("",$_[0],":",$start,"-",$end);


		#print "LEFT TARGET=$target\t";
		#get reads on minus strand
		#my $hits=`samtools view -q $MapQ -f 16 -c $random_file_disc $target`;
		#print "NUM HITS=$hits\n";

		my $hits=`samtools view -q $MapQ -f 16 $random_file_disc $target | awk -F'\t' '{print \$3"\t"\$7"\t"\$8}' > $tmp_hits`;
		$hits=`awk -F '\t' '{if(\$2~/=/){\$2=\$1}else{\$2=\$2};end2=\$3+1;print \$2"\t"\$3"\t"end2;}' $tmp_hits | awk 'NF==3' | sortBed -i stdin | mergeBed -n -d $dist_To_Soft -i stdin | sort -k4,4nr | head -1`;

		#Count number of unmated reads on minus strand
		my $num_unmapped_pairs=`samtools view -q $MapQ -f24 -c $INPUT_BAM $target`;
		$num_unmapped_pairs=~s/\n//;


		my $m_info=$hits;
		next if (! defined ($m_info));
		#$mate_info=~s/\n//;
		my @mate_info=split(/\t/,$m_info);
		my $nRP=$mate_info[3];
		next if((!$nRP)||($minRP > $nRP));
		next if (! defined $nRP);
		next if (($nRP<$minRP)&&($unmated_pairs > ($num_unmapped_pairs+$nRP))); 	#end if there is no mate info or nRP+uRP<minRP

		my @mate_soft;
		my $mate_chr = "";
		my $mate_start;

		##################################################################################
		# Find out if mates have nearby soft-clips (to refine the breakpoints)
		##################################################################################
		#Look for evidence of soft-clipping near mate
		if($m_info){
			@mate_info=split(/\t/,$m_info);
			$start=int($mate_info[2]+$dist_To_Soft);
			#print "LOOKING FOR $mate_info[0]\t$mate_info[1]\t$start in $left_sc\n";
			my $m_soft=`echo -e \"$mate_info[0]\t$mate_info[1]\t$start\"| awk -F'\t' 'NR==3' | intersectBed -a stdin -b $left_sc|head -1`;
			$m_soft=~s/\n//g;
			@mate_soft=split(/\s/, $m_soft);
		}
		if(@mate_soft && scalar(@mate_soft)){
			$mate_chr=$mate_soft[0];
			$mate_start=$mate_soft[1];
		}
		else{
			@mate_info=split(/\s/, $m_info);
			$mate_chr=$mate_info[0];
			$mate_start=$mate_info[1];
		}
		next if ($mate_chr eq ""); 	#end if there is no mate info
		#########################################
		# Get orientation info
		#########################################
		#if sam diff chr{CTX}, elsif INS, or DEL
		#LEFT CLIP;
		my $CTX=`samtools view -q $MapQ -F1024 -f 16 $random_file_disc $target|awk '(\$7 !~/=/)'|wc -l|awk '{print "CTX,"\$0}'`;$CTX=~s/\n//;
		my $DEL=`samtools view -q $MapQ -f 16 -F 1056 $random_file_disc $target|awk '{if((\$7 ~ /=/)&&(\$9>$upper_limit)){print }}'|wc -l|awk '{print "DEL,"\$1}'`;$DEL=~s/\n//;
		my $INS=`samtools view -q $MapQ -f 16 -F 1056 $random_file_disc $target|awk '{if((\$7 ~ /=/)&&(\$9<$lower_limit && \$9 > 0 )){print "INS"}}'|wc -l|awk '{print "INS,"\$1}'`;$INS=~s/\n//;
		my $INV=`samtools view -q $MapQ -f 48 -F 1024 $random_file_disc $target|awk '(\$8 !~/0/)'|wc -l|awk '{print "INV,"\$1}'`;$INV=~s/\n//;
		my $TDUP=`samtools view -q $MapQ  -F 1056 -f 16 $random_file_disc $target |awk '(\$8>\$4)'|wc -l|awk '{print "TDUP,"\$1}'`;$TDUP=~s/\n//;
        my $NOV_INS=`samtools view  -q $MapQ  -F 1024 -f 24  $random_file_disc $target|wc -l|awk '{print "NOVEL_INS,"\$1}'`;$NOV_INS=~s/\n//;
		#print "$CTX,$DEL,$INS,$INV,$TDUP,$NOV_INS\t";
		my $TYPE="";
		#If values are not defined, set them to "0"
		if(!$CTX){$CTX="CTX,0"}
        if(!$DEL){$DEL="DEL,0"}
        if(!$INS){$INS="INS,0"}
        if(!$INV){$INV="INV,0"}
		if(!$TDUP){$TDUP="TDUP,0"}
        if(!$NOV_INS){$NOV_INS="NOV_INS,0"}
		my @data=("$CTX","$DEL","$INS","$INV","$TDUP");#print "@data\n";
		my $max=0;
		my $type="UNKNOWN";
		for (my $i=0;$i<@data;$i++){
			my @tmp=split(",",$data[$i]);
			if($tmp[1]>=$max){$max=$tmp[1];$type=$tmp[0]}
		}
		$TYPE=$type;
		#print "TYPE=$TYPE\n";

		#########################################
		# Get DNA info
		#########################################
		#print "PLUS_READS=$plus_Reads[0],$plus_Reads[1]\nMATE=$mate_chr,$mate_start,$INPUT_FASTA\n";
		my $REF1_base=getSeq($plus_Reads[0],$plus_Reads[1],$INPUT_FASTA);
		my $REF2_base=getSeq($mate_chr,$mate_start,$INPUT_FASTA);

		#########################################
		# print in VCF format
		#########################################
		my $BND1_name=join "", map { ("a".."z")[rand 26] } 1..3;$BND1_name=join "_","BND",$BND1_name;
		my $BND2_name=join "", map { ("a".."z")[rand 26] } 1..3;$BND2_name=join "_","BND",$BND2_name;
		my $lSC=$plus_Reads[4];
		my $nSC=$plus_Reads[6];
		$nRP=~s/\n//;
		my $isize=$plus_Reads[1]-$mate_start;
		my $INFO_1=join(";","EV=LEFTCLIP,SVTYPE=BND","EVENT=$TYPE","ISIZE=$isize","lSC=$lSC","nSC=$nSC","nRP=$nRP","uRP=$num_unmapped_pairs","MATE_ID=$BND2_name");
		my $INFO_2=join(";","EV=LEFTCLIP,SVTYPE=BND","EVENT=$TYPE","ISIZE=$isize","lSC=$lSC","nSC=$nSC","nRP=$nRP","uRP=$num_unmapped_pairs","MATE_ID=$BND1_name");

		my $ALT_1 = "";
		my $ALT_2 = "";

		if($TYPE !~ /cCTX|INV/){
			$ALT_1=join("","]",$mate_chr,":",$mate_start,"]",$REF1_base);
			$ALT_2=join("",$REF2_base,"[",$plus_Reads[0],":",$plus_Reads[1],"[");
			#		2      321682 bnd_V  T   ]13:123456]T  6    PASS SVTYPE=BND
			#		13     123456 bnd_U  C   C[2:321682[   6    PASS SVTYPE=BND
		} else {
			$ALT_1=join("","]",$plus_Reads[0],":",$plus_Reads[1],"]",$REF2_base);
			$ALT_2=join("",$REF1_base,"[",$mate_chr,":",$mate_start,"[");
		}
		if(($mate_chr=~/chr/)&&($plus_Reads[0]=~/chr/)){
			print OUT join ("\t",$plus_Reads[0],$plus_Reads[1],".",$BND1_name,$ALT_1,".","PASS",$INFO_1,"\n");
			print OUT join ("\t",$mate_chr,$mate_start,".",$BND2_name,$ALT_2,".","PASS",$INFO_2,"\n");
		}
	}

	if($_[5] eq "right"){
		#print "RIGHT############################################\n$_\n";
		#my $p_Reads=$_;
		my @plus_Reads = split(/\s/,$_);
		$plus_Reads[7]=~s/\n//g;

		#############
		#Get all types of compatible reads
		#Get improperly paired reads (@ max distance)
		my $start=$_[1]-$dist_To_Soft;
		my $end=$_[2];
		my $target=join("",$_[0],":",$start,"-",$end);


		#Count number of unmated reads on plus strand
		my $num_unmapped_pairs=`samtools view -q $MapQ -f8 -c $INPUT_BAM $target`;
		next if($unmated_pairs > $num_unmapped_pairs);

		#print "RIGHT TARGET=$target\t";
		#get reads on plus strand
		#my $hits=`samtools view -q $MapQ -F 16 -c $random_file_disc $target`;
		#print "NUM HITS=$hits\n";
		my $hits=`samtools view -q $MapQ -F 16 $random_file_disc $target | awk -F'\t' '{print \$3"\t"\$7"\t"\$8}' > $tmp_hits`;
		$hits=`awk -F '\t' '{if(\$2~/=/){\$2=\$1}else{\$2=\$2};end2=\$3+1;print \$2"\t"\$3"\t"end2;}' $tmp_hits | awk 'NF==3' | sortBed -i stdin | mergeBed -n -d $dist_To_Soft -i stdin | sort -k 4,4nr | head -1`;

		#Count number of unmated reads on plus strand
		$num_unmapped_pairs=`samtools view -q $MapQ -f8 -c $INPUT_BAM $target`;
		$num_unmapped_pairs=~s/\n//;

		#if (defined($hits))	{
		#	print "HITS=\n$hits\n";
		#}else	{ #print "HITS=\n no hit\n";
		#}

		my $m_info=$hits;
		next if (! defined ($m_info));
		#$mate_info=~s/\n//;
		my @mate_info=split(/\t/,$m_info);
		my $nRP=$mate_info[3];
		next if (! defined $nRP);next if($minRP > $nRP);
	    next if (($nRP<$minRP)&&($unmated_pairs > ($num_unmapped_pairs+$nRP)));  #end if there is no mate info or nRP+uRP<minRP

		my @mate_soft;
		my $mate_chr = "";
		my $mate_start;

		##################################################################################
		# Find out if mates have nearby soft-clips (to refine the breakpoints)
		##################################################################################
		#Look for evidence of soft-clipping near mate
		if($m_info){
			@mate_info=split(/\t/, $m_info);
			$start=int($mate_info[2]+$dist_To_Soft);
#			print "LOOKING FOR $mate_info[0]\t$mate_info[1]\t$start in $left_sc\n";
			my $m_soft=`echo -e \"$mate_info[0]\t$mate_info[1]\t$start\"| awk -F'\t' 'NR==3' | intersectBed -a stdin -b $left_sc|head -1`;
			$m_soft=~s/\n//g;
			@mate_soft = split(/\s/, $m_soft);
		}
		if(@mate_soft){
			$mate_chr=$mate_soft[0];
			$mate_start=$mate_soft[1];
		}
		else{
			@mate_info = split(/\s/, $m_info);
			$mate_chr=$mate_info[0];
			$mate_start=$mate_info[1];
		}
		next if ($mate_chr eq ""); 	#end if there is no mate info
		#########################################
		# Get orientation info
		#########################################
		#if sam diff chr{CTX}, elsif INS, or DEL
		###RIGHT
        my $CTX=`samtools view -q $MapQ -F1040 -f 32 $random_file_disc $target|awk '(\$7 !~/=/)'|wc -l|awk '{print "CTX,"\$0}'`;$CTX=~s/\n//;
        my $DEL=`samtools view -q $MapQ -f 32 -F 1040 $random_file_disc $target|awk '{if((\$7 ~ /=/)&&(\$9>$upper_limit)){print "DEL"}}'|wc -l|awk '{print "DEL,"\$1}'`;$DEL=~s/\n//;
        my $INS=`samtools view -q $MapQ -f 32 -F 1040 $random_file_disc $target|awk '{if((\$7 ~ /=/)&&(\$9<$lower_limit && \$9 > 0)){print "INS"}}'|wc -l|awk '{print "INS,"\$1}'`;$INS=~s/\n//;
		my $INV=`samtools view -q $MapQ -F 1071 $random_file_disc $target|awk '(\$8 !~/0/)'|wc -l|awk '{print "INV,"\$1}'`;$INV=~s/\n//;
        my $TDUP=`samtools view -q $MapQ -f 32 -F1024  $random_file_disc $target |awk '(\$8<\$4)'|wc -l|awk '{print "TDUP,"\$1}'`;$TDUP=~s/\n//;
		my $NOV_INS=`samtools view  -q $MapQ  -F 1040 -f 8  $random_file_disc $target|wc -l|awk '{print "NOVEL_INS,"\$1}'`;$NOV_INS=~s/\n//;
		my $TYPE="";

		#print "$CTX,$DEL,$INS,$INV,$TDUP,$NOV_INS\t";

		#If values are not defined, set them to "0"
		if(!$CTX){$CTX="CTX,0"}
		if(!$DEL){$DEL="DEL,0"}
		if(!$INS){$INS="INS,0"}
		if(!$INV){$INV="INV,0"}
		if(!$TDUP){$TDUP="TDUP,0"}
		if(!$NOV_INS){$NOV_INS="NOVEL_INS,0"}
        my @data=("$CTX","$DEL","$INS","$INV","$TDUP","$NOV_INS");#print "@data\n";
        my $max=0;
		my $type="UNKNOWN";
        for (my $i=0;$i<@data;$i++){
            my @tmp=split(",",$data[$i]);
			if($tmp[1]>=$max){$max=$tmp[1];$type=$tmp[0]}
		}

		$TYPE=$type;
		#print "TYPE=$TYPE\n";

		#########################################
		# Get DNA info
		#########################################

		my $REF1_base=getSeq($plus_Reads[0],$plus_Reads[1],$INPUT_FASTA);
		my $REF2_base=getSeq($mate_chr,$mate_start,$INPUT_FASTA);

		#########################################
		# print in VCF format
		#########################################
		my $BND1_name=join "", map { ("a".."z")[rand 26] } 1..3;$BND1_name=join "_","BND",$BND1_name;
		my $BND2_name=join "", map { ("a".."z")[rand 26] } 1..3;$BND2_name=join "_","BND",$BND2_name;
		my $lSC=$plus_Reads[4];
		my $nSC=$plus_Reads[6];
		$nRP=~s/\n//;
		my $isize=$mate_start-$plus_Reads[1];
		my $INFO_1=join(";","EV=RIGHTCLIP,SVTYPE=BND","EVENT=$TYPE","ISIZE=$isize","lSC=$lSC","nSC=$nSC","nRP=$nRP","uRP=$num_unmapped_pairs","MATE_ID=$BND2_name");
		my $INFO_2=join(";","EV=RIGHT_CLIP,SVTYPE=BND","EVENT=$TYPE","ISIZE=$isize","lSC=$lSC","nSC=$nSC","nRP=$nRP","uRP=$num_unmapped_pairs","MATE_ID=$BND1_name");

		my $ALT_1;
		my $ALT_2;

		if($TYPE !~ /cCTX|INV/){
			$ALT_1=join("","]",$mate_chr,":",$mate_start,"]",$REF1_base);
			$ALT_2=join("",$REF2_base,"[",$plus_Reads[0],":",$plus_Reads[1],"[");
			#		2      321682 bnd_V  T   ]13:123456]T  6    PASS SVTYPE=BND
			#		13     123456 bnd_U  C   C[2:321682[   6    PASS SVTYPE=BND
		} else {
			$ALT_1=join("","]",$plus_Reads[0],":",$plus_Reads[1],"]",$REF2_base);
			$ALT_2=join("",$REF1_base,"[",$mate_chr,":",$mate_start,"[");
		}

		if(($mate_chr=~/chr/)&&($plus_Reads[0]=~/chr/)){
			print OUT join ("\t",$plus_Reads[0],$plus_Reads[1],".",$BND1_name,$ALT_1,".","PASS",$INFO_1,"\n");
			print OUT join ("\t",$mate_chr,$mate_start,".",$BND2_name,$ALT_2,".","PASS",$INFO_2,"\n");
		}
	}
		################################
		#
}
close FILE;
close OUT;

#############################################################
#############################################################
# Delete temp files
my $meregedBed=join("",$random_name,".sc.direction.sc.mergeBed");
if(defined($temp_output)){$temp_output=$temp_output} else {$temp_output="no"}
if ($temp_output eq "no"){
	unlink("$tmp_name","$random_file","$tmp2",,"$tmp3","$outfile","$random_file_sc","$random_file_disc","$index","$random_name","$index2","$right_sc","$left_sc","$meregedBed");
}
print "Analysis Completed\n\nYou did it!!!\n";
print "Finish Time : " . &spGetCurDateTime() . "\n";
$now = time - $now;
printf("\n\nTotal running time: %02d:%02d:%02d\n\n", int($now / 3600), int(($now % 3600) / 60),
int($now % 60));
exit;

sub rev_comp {
  my $dna = @_;
  my $revcomp = reverse($dna);
  $revcomp =~ tr/ACGTacgt/TGCAtgca/;
  return $revcomp;
}


## to get reference base
sub getSeq{
	my ($chr,$pos,$fasta)=@_;
	if($chr !~ /^chr/){die "$chr is not correct\n";}
	if ($pos <0){die "$pos is not a number\n"}
	my @result = `samtools faidx $fasta $chr:$pos-$pos`;
	chomp($result[1]);
	return uc($result[1]);
}

## to get time
sub spGetCurDateTime {
	my ($sec, $min, $hour, $mday, $mon, $year) = localtime();
	my $curDateTime = sprintf "%4d-%02d-%02d %02d:%02d:%02d",
	$year+1900, $mon+1, $mday, $hour, $min, $sec;
	return ($curDateTime);
}

### print header
sub print_header {
	my $date=&spGetCurDateTime();
	my $header = qq{##fileformat=VCFv4.1
##fileDate=$date
##source=SoftSearch.pl
##reference=$INPUT_FASTA
##Usage= SoftSearch.pl -l $minSoft -q 20 -r $minSoftReads -d $dist_To_Soft -m $minRP -u $unmated_pairs -s $num_sd -b $INPUT_BAM -f $INPUT_FASTA -o $OUTPUT_FILE
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=EVENT,Number=1,Type=String,Description="ID of event associated to breakend">
##INFO=<ID=lSC,Number=1,Type=Integer,Description="Length of the longest soft clips supporting the BND">
##INFO=<ID=nSC,Number=1,Type=Integer,Description="Number of supporting soft-clips\">
##INFO=<ID=nRP,Number=1,Type=Integer,Description="Number of read pairs overlapping Soft-Clips">
##INFO=<ID=uRP,Number=1,Type=Integer,Description="Number of unmated read pairs nearby Soft-Clips">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n};
	print OUT $header;
}
