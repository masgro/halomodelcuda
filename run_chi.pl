#!/usr/bin/perl
use strict;
use warnings;

my $eje = 'hm.cono';
my $nchi = 31;
my $m1f = $ARGV[0];
my $m1 = sprintf "%4.2f", $m1f;
my $m2f = $ARGV[1];
my $m2 = sprintf "%4.2f", $m2f;
my $bcmin =   0.2;
my $bcmax =   1.4;
my $abmin =   0.2;
my $abmax =   1.4;
my $dbc = ($bcmax - $bcmin)/($nchi - 1);
my $dab = ($abmax - $abmin)/($nchi - 1);

my $angulof = 45.0;
my $angulo  = sprintf "%02d", $angulof;

my $term = "2h";
my $nrun = 1;

for(my $i = 0; $i < $nchi; $i++) {
	for(my $j = 0; $j < $nchi; $j++) {

		#SI QUEREMOS HACER CHI2 EN T1h
		#my $BCMEDIO = $bcmin + $dbc * $i; 
		#my $ABMEDIO = $abmin + $dab * $j; 
		#my $bcmedio = sprintf "%.2f", $BCMEDIO;
		#my $abmedio = sprintf "%.2f", $ABMEDIO;

		##SI NO QUEREMOS ALINEAMIENTO EN T2h
		#my $align_b = "0.00";
		#my $align_c = "0.00";

		#SI QUEREMOS HACER CHI2 EN T2h
		my $ALIGN_B = $bcmin + $dbc * $i; 
		my $ALIGN_C = $abmin + $dab * $j; 

    #SI QUEREMOS QUE VARIE EXPONENCIALMENTE
    #$ALIGN_B = exp($ALIGN_B);
    #$ALIGN_C = exp($ALIGN_C);
	
		my $align_b = sprintf "%.2f", $ALIGN_B;
		my $align_c = sprintf "%.2f", $ALIGN_C;
    
    #SI QUEREMOS QUE LOS HALOS SEAN REDONDOS
    my $bcmedio = 0.99;
    my $abmedio = 0.99;

		system("make clean");
  	my $cmd = "make term=$term angulo=$angulo bcmedio=$bcmedio abmedio=$abmedio align_b=$align_b align_c=$align_c m1=$m1 m2=$m2 $eje";
		system($cmd);

		for(my $k = 0; $k < $nrun; $k++)
		{
  		my $cmd = "$eje $k";
  		system($cmd);
		}

		system("make clean");
  	$cmd = "make term=$term angulo=$angulo bcmedio=$bcmedio abmedio=$abmedio align_b=$align_b align_c=$align_c	m1=$m1 m2=$m2 reescribe.x";
		system($cmd);

  	$cmd = "reescribe.x $nrun";
		system($cmd);
	}
	my $num = sprintf "%02d", $i;
	my $cmd="nvidia-smi -f log.$num";
	system($cmd);
}
