#!/usr/bin/perl
use strict;
use warnings;

my $exec     = 'hm.cono';
my $nrun     =  1;
my $term     = "1h";
my $kind     = "-DCG";
my $mag      = 21;

my $m1f = 12.00;
my $m2f = 15.00;
my $bcmediof = 1.0;
my $abmediof = 1.0;
my $align_bf = 1.4;
my $align_cf = 1.4;
my $angulof  = 45.0;

my $m1 = sprintf "%4.2f", $m1f;
my $m2 = sprintf "%4.2f", $m2f;
my $bcmedio = sprintf "%.2f", $bcmediof;
my $abmedio = sprintf "%.2f", $abmediof;
my $align_b = sprintf "%.2f", $align_bf;
my $align_c = sprintf "%.2f", $align_cf;
my $angulo  = sprintf "%02d", $angulof;

system("make clean");
my $cmd = "make kind=$kind term=$term angulo=$angulo bcmedio=$bcmedio abmedio=$abmedio align_b=$align_b align_c=$align_c m1=$m1 m2=$m2 mag=$mag $exec";
system($cmd);

for(my $k = 0; $k < $nrun; $k++){
	my $cmd = "$exec $k";
	system($cmd);
}

system("make clean");
$cmd = "make kind=$kind term=$term angulo=$angulo bcmedio=$bcmedio abmedio=$abmedio align_b=$align_b align_c=$align_c	m1=$m1 m2=$m2 mag=$mag reescribe.x";
system($cmd);

$cmd = "reescribe.x $nrun";
system($cmd);

$cmd = "rm -rf funcorr_?_?h.??";
system($cmd);

###################
#TERMINO DE 2-HALOS
###################
$term     = "2h";
$bcmediof = 1.00;
$abmediof = 1.00;

$m1 = sprintf "%4.2f", $m1f;
$m2 = sprintf "%4.2f", $m2f;
$bcmedio = sprintf "%.2f", $bcmediof;
$abmedio = sprintf "%.2f", $abmediof;
$align_b = sprintf "%.2f", $align_bf;
$align_c = sprintf "%.2f", $align_cf;
$angulo  = sprintf "%02d", $angulof;

system("make clean");
$cmd = "make kind=$kind term=$term angulo=$angulo bcmedio=$bcmedio abmedio=$abmedio align_b=$align_b align_c=$align_c m1=$m1 m2=$m2 mag=$mag $exec";
system($cmd);

for(my $k = 0; $k < $nrun; $k++){
	my $cmd = "$exec $k";
	system($cmd);
}

system("make clean");
$cmd = "make kind=$kind term=$term angulo=$angulo bcmedio=$bcmedio abmedio=$abmedio align_b=$align_b align_c=$align_c	m1=$m1 m2=$m2 mag=$mag reescribe.x";
system($cmd);

$cmd = "reescribe.x $nrun";
system($cmd);

system("make clean");
$cmd = "rm -rf funcorr_?_?h.??";
system($cmd);
