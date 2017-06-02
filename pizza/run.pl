#!/usr/bin/perl
use strict;
use warnings;

my $argc = @ARGV;
if(($argc < 3) || ($ARGV[1] > $ARGV[2])){
	print "Usage: ./run_chi.pl Mag Mass_min Mass_max\n";
	exit();
}

my $exec = 'hm.cono';
my $nrun =         1;
my $term =      "1h";
my $kind =    "-DCG";
#my $mag =      "17";
#my $m1f = 13.74;
#my $m1  = sprintf "%4.2f", $m1f;
#my $m2f = 15.00;
#my $m2  = sprintf "%4.2f", $m2f;

my $mag = $ARGV[0];
my $m1  = $ARGV[1];
my $m2  = $ARGV[2];

my $bcmediof = 1.00;
my $bcmedio = sprintf "%.2f", $bcmediof;
my $abmediof = 1.00;
my $abmedio = sprintf "%.2f", $abmediof;

my $align_bf = 1.50;
my $align_b = sprintf "%.2f", $align_bf;
my $align_cf = 1.50;
my $align_c = sprintf "%.2f", $align_cf;

my $angulof = 45;
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

system("make clean");
$cmd = "rm -rf funcorr_?_1h.??";
system($cmd);

#die("Stopping...");

$term =  "2h";
$nrun =    1;

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
$cmd = "rm -rf funcorr_?_2h.??";
system($cmd);
