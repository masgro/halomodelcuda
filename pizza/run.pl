#!/usr/bin/perl
use strict;
use warnings;

my $argc = @ARGV;
if(($argc < 3) || ($ARGV[1] > $ARGV[2])){
  print "Usage: ./run.pl Mag Mass_min Mass_max\n";
  exit();
}

my $mag = $ARGV[0];
my $m1  = $ARGV[1];
my $m2  = $ARGV[2];

my $exec = 'hm.cono.x';
my $nrun =      1;
my $term =        "1h";
my $kind =      "-DCG";

my $bcmediof = 0.75;
my $abmediof = 0.75;
my $align_bf = 1.90;
my $align_cf = 1.90;
my $angulof  = 90.0;

my $bcmedio = sprintf "%.2f", $bcmediof;
my $abmedio = sprintf "%.2f", $abmediof;
my $align_b = sprintf "%.2f", $align_bf;
my $align_c = sprintf "%.2f", $align_cf;
my $angulo  = sprintf "%02d", $angulof;

my $cmd="hm.cono.x $m1 $m2 $abmedio $bcmedio $align_b $align_c";
system($cmd);

$cmd = "reescribe.x $nrun $m1 $m2 $abmedio $bcmedio $align_b $align_c";
system($cmd);

#die("Stopping...");

###################
#TERMINO DE 2-HALOS
###################
$term     = "2h";
#$bcmediof = 0.80;
#$abmediof = 0.80;
$nrun     = 1;

$bcmedio = sprintf "%.2f", $bcmediof;
$abmedio = sprintf "%.2f", $abmediof;
$align_b = sprintf "%.2f", $align_bf;
$align_c = sprintf "%.2f", $align_cf;
$angulo  = sprintf "%02d", $angulof;

for(my $k = 0; $k < $nrun; $k++){
  $cmd="hm.cono.2h.x $m1 $m2 $abmedio $bcmedio $align_b $align_c $k";
  system($cmd);
}

$cmd = "reescribe.2h.x $nrun $m1 $m2 $abmedio $bcmedio $align_b $align_c";
system($cmd);
