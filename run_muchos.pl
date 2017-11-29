#!/usr/bin/perl
use strict;
use warnings;

my $exec = 'hm';

my $nrun = 5;

my $term = '2h';

my $angulof = 1.0;
my $angulo  = sprintf "%02d", $angulof;

my @m1v = (12.00,12.00,12.50,13.00,13.50,14.00,14.50);
my @m2v = (16.00,12.50,13.00,13.50,14.00,14.50,15.00);

my @bcmediov = (0.78,0.79,0.68,0.56,0.58,0.56);
my @abmediov = (0.89,0.84,0.89,0.89,0.74,0.78);

my @align_bv = (0.91, 0.00, 0.00, 0.62, 0.00, 0.57);
my @align_cv = (0.35, 0.00, 0.00, 0.28, 0.00, 0.32);

my $bin_mass_max = 1;
my $tipo_max = 3;

my $i = 0;
do {

	my $m1 = sprintf "%.2f", $m1v[$i];
	my $m2 = sprintf "%.2f", $m2v[$i];
	
	#my $bcmedio = sprintf "%.2f", $bcmediov[$i];
	#my $abmedio = sprintf "%.2f", $abmediov[$i];

	my $bcmedio = "0.80";
	my $abmedio = "0.80";

	my $tipo = 0;
	do {
    #my $align_b = sprintf "%.2f", $align_bv[$i];
    #my $align_c = sprintf "%.2f", $align_cv[$i];

    my $align_b = "0.50";
    my $align_c = "0.35";
  
    system("make clean");
    my $cmd = "make tipo=$tipo term=$term angulo=$angulo bcmedio=$bcmedio abmedio=$abmedio align_b=$align_b align_c=$align_c m1=$m1 m2=$m2 $exec";
	  print $cmd;
    system($cmd);
    
    for(my $k = 0; $k < $nrun; $k++)
    {
    	my $cmd = "$exec $k";
    	system($cmd);
    }
    
    system("make clean");
    $cmd = "make tipo=$tipo term=$term angulo=$angulo bcmedio=$bcmedio abmedio=$abmedio align_b=$align_b align_c=$align_c m1=$m1 m2=$m2 reescribe.x";
    system($cmd);
    
    $cmd = "reescribe.x $nrun";
    system($cmd);

		$tipo++;
	}until($tipo == $tipo_max);
  
  $i++;
} until ($i == $bin_mass_max);

