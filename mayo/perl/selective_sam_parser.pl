#!/usr/bin/perl

print STDERR "here\n";

open(MYINPUTFILE,"<", $ARGV[0]); # open stdin for input

print STDERR "file open, reading file...\n";
my(@lines) = <MYINPUTFILE>; # read file into list


$bt=10;
$at=9;

print STDERR "File read\n";
$linetot=scalar(@lines);

for ($i=0; $i <= $linetot; $i++) {

    ($jk1,$jk2,$jk3,$jk4,$jk5,$jk6,$col7,$compvalue,$jkrest)=split(/\t/,$lines[$i]);
    #print $compvalue,"\n";
   $printme="false";
    for ( $bv=1; $bv <= $bt; $bv++) {
       if ( (($i+1-$bv) >= 0)  )   {
            ($jk1,$jk2,$jk3,$jk4,$jk5,$jk6,$prevcol7,$prevcompvalue,$jkrest)=split(/\t/,$lines[$i-$bv]);
            if ( (abs($prevcompvalue-$compvalue) < 10000) && ( $prevcol7 == $col7) ) {
		if ($lines[$i] =~ /R0211812_0073:4:1:10020:173193#0/){
			print "LOOK UP\n";
			print $lines[$i]."\n".$lines[$i-$bv]."\n";
			print $prevcompvalue."\t".$compvalue."\t".$prevcol7."\t".$col7."\t".$i-$bv."\n";
		}
               $printme="true";
            }
       }
    }
    for ( $av=1; $av <= $at; $av++) {
       if ( (($i+$av) <= $linetot)  )   {
            ($jk1,$jk2,$jk3,$jk4,$jk5,$jk6,$avcol7,$avcompvalue,$jkrest)=split(/\t/,$lines[$i+$av]);
            if ( (abs($avcompvalue-$compvalue) < 10000) && ( $avcol7 == $col7) ) {
		if ($lines[$i] =~ /R0211812_0073:4:1:10020:173193#0/){
			print "LOOK DOWN\n";
                        print $lines[$i]."\n".$lines[$i-$av]."\n";
			print "a: ".$avcompvalue."\t";
			print "b: ".$compvalue."\t";
			print "c: ".$avcol7."\t";
			print "d: ".$col7."\t";
			print "e: ".($i+$av)."\n";
                }
               $printme="true";
            }
       }

    }

    if ( $printme eq "true" ) {
        #print "$lines[$i]"; # print in order
   }
}

close(MYINPUTFILE);
