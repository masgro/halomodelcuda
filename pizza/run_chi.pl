#!/usr/bin/perl
use strict;
use warnings;
use diagnostics;

my $argc = @ARGV;
if(($argc < 2) || ($ARGV[0] > $ARGV[1])){
	print "$ARGV[0] $ARGV[1]\n";
  print "Usage: ./run_chi.pl Mass_min Mass_max\n";
  exit();
}

my $folder = "/home/msgro/halomodelcuda/pizza";

my $mmin  = $ARGV[0];
my $mmax  = $ARGV[1];

my $nchi  =  40+1;
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

  	my $cmd = "$folder/hm.cono.x $mmin $mmax $abmedio $bcmedio $align_b $align_c";
  	system($cmd);

  	$cmd = "$folder/reescribe.x $nrun $mmin $mmax $abmedio $bcmedio $align_b $align_c";
		system($cmd);
	}
	#my $num = sprintf "%02d", $i;
	#my $cmd="nvidia-smi -f log.$num";
	#system($cmd);
}

#die("Stopping...");

$nchi =   24+1;
$bcmin = 0.10;
$bcmax = 1.30;
$abmin = 0.10;
$abmax = 1.30;
$dbc = ($bcmax - $bcmin)/($nchi - 1);
$dab = ($abmax - $abmin)/($nchi - 1);
$nrun = 10;

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
    my $bcmedio = "0.500";
    my $abmedio = "0.500";

    for(my $k = 0; $k < $nrun; $k++) {
  	  my $cmd = "$folder/hm.cono.2h.x $mmin $mmax $abmedio $bcmedio $align_b $align_c $k";
  	  system($cmd);
	  }

  	my $cmd = "$folder/reescribe.2h.x $nrun $mmin $mmax $abmedio $bcmedio $align_b $align_c";
		system($cmd);

    $cmd = "rm -rf funcorr_?_?h.??";
    system($cmd);
	}
	#my $num = sprintf "%02d", $i;
	#my $cmd = "nvidia-smi -f log.$num";
	#system($cmd);
}
