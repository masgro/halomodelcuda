#!/usr/bin/perl
use strict;
use warnings;
use diagnostics;

my $nchi  =    9;
my $bcmin = 0.20;
my $bcmax = 1.00;
my $abmin = 0.20;
my $abmax = 1.00;
my $dbc = ($bcmax - $bcmin)/($nchi - 1);
my $dab = ($abmax - $abmin)/($nchi - 1);

my $nrun = 1;

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
		my $align_b = "1.000";
		my $align_c = "1.000";

  	my $cmd = "hm.cono.x $abmedio $bcmedio $align_b $align_c";
  	system($cmd);

  	$cmd = "reescribe.x $nrun $abmedio $bcmedio $align_b $align_c";
		system($cmd);
	}
	my $num = sprintf "%02d", $i;
	my $cmd="nvidia-smi -f log.$num";
	system($cmd);
}

#die("Stopping...");

$nchi  =    9;
$bcmin = 0.20;
$bcmax = 1.00;
$abmin = 0.20;
$abmax = 1.00;
$dbc = ($bcmax - $bcmin)/($nchi - 1);
$dab = ($abmax - $abmin)/($nchi - 1);

for(my $i = 0; $i < $nchi; $i++) {
	for(my $j = 0; $j < $nchi; $j++) {

		my $filename = 'stop';
		if (-e $filename) {
			die("Stopping...");
		}

		#SI QUEREMOS HACER CHI2 EN T2h
		my $ALIGN_B = $bcmin + $dbc * $i; 
		my $ALIGN_C = $abmin + $dab * $j; 

		my $align_b = sprintf "%.2f", $ALIGN_B;
		my $align_c = sprintf "%.2f", $ALIGN_C;
    
    #SI QUEREMOS QUE LOS HALOS SEAN REDONDOS
    my $bcmedio = "1.000";
    my $abmedio = "1.000";

  	my $cmd = "hm.cono.2h.x $abmedio $bcmedio $align_b $align_c";
  	system($cmd);

  	$cmd = "reescribe.2h.x $nrun $abmedio $bcmedio $align_b $align_c";
		system($cmd);
	}
	my $num = sprintf "%02d", $i;
	my $cmd = "nvidia-smi -f log.$num";
	system($cmd);
}



