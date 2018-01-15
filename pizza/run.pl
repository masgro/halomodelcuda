#!/usr/bin/perl
use strict;
use warnings;

#my $argc = @ARGV;
#if(($argc < 3) || ($ARGV[1] > $ARGV[2])){
#	print "Usage: ./run.pl Mag Mass_min Mass_max\n";
#	exit();
#}

my $exec = 'hm.cono';
my $nrun =         1;
my $term =      "1h";
my $kind =    "-DCG";
#my $mag =      "19";
#my $m1f = 13.37;
#my $m1  = sprintf "%4.2f", $m1f;
#my $m2f = 13.47;
#my $m2  = sprintf "%4.2f", $m2f;

my $mag = $ARGV[0];
my $m1  = $ARGV[1];
my $m2  = $ARGV[2];

my $bcmediof = 0.75;
my $bcmedio = sprintf "%.2f", $bcmediof;
my $abmediof = 0.75;
my $abmedio = sprintf "%.2f", $abmediof;

my $align_bf = 1.90;
my $align_b = sprintf "%.2f", $align_bf;
my $align_cf = 1.90;
my $align_c = sprintf "%.2f", $align_cf;

#my $bcmedio = $ARGV[3];
#my $abmedio = $ARGV[4];
#my $align_b = $ARGV[5];
#my $align_c = $ARGV[6];

my $angulof = 45;
my $angulo  = sprintf "%02d", $angulof;

system("make --quiet clean");
my $cmd = "make --quiet kind=$kind term=$term angulo=$angulo bcmedio=$bcmedio abmedio=$abmedio align_b=$align_b align_c=$align_c m1=$m1 m2=$m2 mag=$mag $exec";
system($cmd);

for(my $k = 0; $k < $nrun; $k++){
	my $cmd = "$exec $k";
	system($cmd);
}

system("make --quiet clean");
$cmd = "make --quiet kind=$kind term=$term angulo=$angulo bcmedio=$bcmedio abmedio=$abmedio align_b=$align_b align_c=$align_c	m1=$m1 m2=$m2 mag=$mag reescribe.x";
system($cmd);

$cmd = "reescribe.x $nrun";
system($cmd);

system("make --quiet clean");
$cmd = "rm -rf funcorr_?_1h.??";
system($cmd);

#die("Stopping...");

$term =  "2h";
$nrun =    1;

$abmediof = 1.00;
$abmedio = sprintf "%.2f", $abmediof;
$bcmediof = 1.00;
$bcmedio = sprintf "%.2f", $bcmediof;

$align_bf = 0.40;
$align_b = sprintf "%.2f", $align_bf;
$align_cf = 0.45;
$align_c = sprintf "%.2f", $align_cf;


#my $num = $ARGV[7];
my $num = 0;

system("make --quiet clean");
$cmd = "make --quiet num=$num kind=$kind term=$term angulo=$angulo bcmedio=$bcmedio abmedio=$abmedio align_b=$align_b align_c=$align_c m1=$m1 m2=$m2 mag=$mag $exec";
system($cmd);

for(my $k = 0; $k < $nrun; $k++){
	my $cmd = "$exec $k";
	system($cmd);
}

system("make --quiet clean");
$cmd = "make --quiet num=$num kind=$kind term=$term angulo=$angulo bcmedio=$bcmedio abmedio=$abmedio align_b=$align_b align_c=$align_c	m1=$m1 m2=$m2 mag=$mag reescribe.x";
system($cmd);

$cmd = "reescribe.x $nrun";
system($cmd);

system("make --quiet clean");
$cmd = "rm -rf funcorr_?_2h.??";
system($cmd);
