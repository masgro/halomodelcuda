#!/usr/bin/perl
use strict;
use warnings;

my $argc = @ARGV;
if(($argc < 2) || ($ARGV[0] > $ARGV[1])){
	print "Usage: ./run_chi.pl Mass_min Mass_max\n";
	exit();
}

my $eje = 'hm.cono';
my $m1f = $ARGV[0];
my $m1 = sprintf "%4.2f", $m1f;
my $m2f = $ARGV[1];
my $m2 = sprintf "%4.2f", $m2f;

my $nchi  =   41;
my $bcmin =   0.20;
my $bcmax =   1.00;
my $abmin =   0.20;
my $abmax =   1.00;
my $dbc = ($bcmax - $bcmin)/($nchi - 1);
my $dab = ($abmax - $abmin)/($nchi - 1);

my $angulof = 45;
my $angulo  = sprintf "%02d", $angulof;

my $term = "1h";
my $nrun = 1;
my $kind = "-DCG";
my $mag  = 19; 

######### BEGIN 1H TERM ######
 for(my $i = 0; $i < $nchi; $i++) {
 	for(my $j = 0; $j < $nchi; $j++) {
 
 		my $filename = 'stop';
 		if (-e $filename) {
       system("rm -rf stop");
 			die("Stopping...");
 		}
 
 		#SI QUEREMOS HACER CHI2 EN T1h
 		my $BCMEDIO = $bcmin + $dbc * $i; 
 		my $ABMEDIO = $abmin + $dab * $j; 
 		my $bcmedio = sprintf "%.2f", $BCMEDIO;
 		my $abmedio = sprintf "%.2f", $ABMEDIO;
 
 		##SI NO QUEREMOS ALINEAMIENTO EN T2h
 		my $align_b = "1.90";
 		my $align_c = "1.90";
 
 		#SI QUEREMOS HACER CHI2 EN T2h
 		#my $ALIGN_B = $bcmin + $dbc * $i; 
 		#my $ALIGN_C = $abmin + $dab * $j; 
 
     #SI QUEREMOS QUE VARIE EXPONENCIALMENTE
     #$ALIGN_B = exp($ALIGN_B);
     #$ALIGN_C = exp($ALIGN_C);
 	
 		#my $align_b = sprintf "%.2f", $ALIGN_B;
 		#my $align_c = sprintf "%.2f", $ALIGN_C;
     
     #SI QUEREMOS QUE LOS HALOS SEAN REDONDOS
     #my $bcmedio = 0.99;
     #my $abmedio = 0.99;
 
 		system("make clean");
   	my $cmd = "make kind=$kind term=$term angulo=$angulo bcmedio=$bcmedio abmedio=$abmedio align_b=$align_b align_c=$align_c m1=$m1 m2=$m2 mag=$mag $eje";
 		system($cmd);
 
 		for(my $k = 0; $k < $nrun; $k++){
   		my $cmd = "$eje $k";
   		system($cmd);
 		}
 
 		system("make clean");
   	$cmd = "make kind=$kind term=$term angulo=$angulo bcmedio=$bcmedio abmedio=$abmedio align_b=$align_b align_c=$align_c	m1=$m1 m2=$m2 mag=$mag reescribe.x";
 		system($cmd);
 
   	$cmd = "reescribe.x $nrun";
 		system($cmd);
 
     $cmd = "rm -rf funcorr_?_?h.??";
     system($cmd);
 	}
 	my $num = sprintf "%02d", $i;
 	my $cmd="nvidia-smi -f log.$num";
 	system($cmd);
 }
## ####### END 1H TERM ######

#die("Stopping...");

$nchi = 41;
$bcmin = 0.4;
$bcmax = 1.2;
$abmin = 0.4;
$abmax = 1.2;
$dbc = ($bcmax - $bcmin)/($nchi - 1);
$dab = ($abmax - $abmin)/($nchi - 1);

$term = "2h";

for(my $i = 0; $i < $nchi; $i++) {
	for(my $j = 0; $j < $nchi; $j++) {

		my $filename = 'stop';
		if (-e $filename) {
			die("Stopping...");
		}

		#SI QUEREMOS HACER CHI2 EN T2h
		my $ALIGN_B = $bcmin + $dbc * $i; 
		my $ALIGN_C = $abmin + $dab * $j; 

    #SI QUEREMOS QUE VARIE EXPONENCIALMENTE
    #$ALIGN_B = exp($ALIGN_B);
    #$ALIGN_C = exp($ALIGN_C);
	
		my $align_b = sprintf "%.2f", $ALIGN_B;
		my $align_c = sprintf "%.2f", $ALIGN_C;
    
    #SI QUEREMOS QUE LOS HALOS SEAN REDONDOS
    my $bcmedio = "1.00";
    my $abmedio = "1.00";

		system("make clean");
  	my $cmd = "make kind=$kind term=$term angulo=$angulo bcmedio=$bcmedio abmedio=$abmedio align_b=$align_b align_c=$align_c m1=$m1 m2=$m2 mag=$mag $eje";
    print $cmd;
		system($cmd);

		for(my $k = 0; $k < $nrun; $k++){
  		my $cmd = "$eje $k";
  		system($cmd);
		}

		system("make clean");
  	$cmd = "make kind=$kind term=$term angulo=$angulo bcmedio=$bcmedio abmedio=$abmedio align_b=$align_b align_c=$align_c	m1=$m1 m2=$m2 mag=$mag reescribe.x";
		system($cmd);

  	$cmd = "reescribe.x $nrun";
		system($cmd);

    $cmd = "rm -rf funcorr_?_?h.??";
    system($cmd);
	}
	my $num = sprintf "%02d", $i;
	my $cmd="nvidia-smi -f log.$num";
	system($cmd);
}



