#!/usr/bin/perl
use strict;
use warnings;

my $eje = 'hm.cono';
#my $m1f = 14.0;
#my $m1 = sprintf "%4.2f", $m1f;
#my $m2f = 14.5;
#my $m2 = sprintf "%4.2f", $m2f;

my $nrun = 20;

##my $bcmediof = 0.7;
##my $bcmedio = sprintf "%3.1f", $bcmediof;
##my $abmediof = 0.8;
##my $abmedio = sprintf "%3.1f", $abmediof;

my $angulof = 45.0;
my $angulo  = sprintf "%02d", $angulof;

my @m1v = (12.00,12.50,13.00,13.50,14.00,14.50);
my @m2v = (12.50,13.00,13.50,14.00,14.50,15.00);

my @bcmediov = (0.79,0.79,0.68,0.58,0.58,0.53);
my @abmediov = (0.89,0.84,0.89,0.84,0.74,0.68);

my $i = 0;
do {

	my $m1 = sprintf "%4.2f", $m1v[$i];
	my $m2 = sprintf "%4.2f", $m2v[$i];
	
	my $bcmedio = sprintf "%.2f", $bcmediov[$i];
	my $abmedio = sprintf "%.2f", $abmediov[$i];

  ####NOALIGN
  my $align_bf = 0.0;
  my $align_b = sprintf "%5.3f", $align_bf;
  my $align_cf = 0.0;
  my $align_c = sprintf "%5.3f", $align_cf;
  
  system("make clean");
  my $cmd = "make tipo=0 term=2h angulo=$angulo bcmedio=$bcmedio abmedio=$abmedio align_b=$align_b align_c=$align_c m1=$m1 m2=$m2 $eje";
	print $cmd;
  system($cmd);
  
  for(my $k = 0; $k < $nrun; $k++)
  {
  	my $cmd = "$eje $k";
  	system($cmd);
  }
  
  system("make clean");
  $cmd = "make term=2h angulo=$angulo bcmedio=$bcmedio abmedio=$abmedio align_b=$align_b align_c=$align_c	m1=$m1 m2=$m2 reescribe.x";
  system($cmd);
  
  $cmd = "reescribe.x $nrun";
  system($cmd);
  
  
  
  ####ALIGNB
  my $align_bf = 1.0;
  my $align_b = sprintf "%5.3f", $align_bf;
  my $align_cf = 0.0;
  my $align_c = sprintf "%5.3f", $align_cf;
  
  system("make clean");
  my $cmd = "make tipo=1 term=2h angulo=$angulo bcmedio=$bcmedio abmedio=$abmedio align_b=$align_b align_c=$align_c m1=$m1 m2=$m2 $eje";
  system($cmd);
  
  for(my $k = 0; $k < $nrun; $k++)
  {
  	my $cmd = "$eje $k";
  	system($cmd);
  }
  
  system("make clean");
  $cmd = "make term=2h angulo=$angulo bcmedio=$bcmedio abmedio=$abmedio align_b=$align_b align_c=$align_c	m1=$m1 m2=$m2 reescribe.x";
  system($cmd);
  
  $cmd = "reescribe.x $nrun";
  system($cmd);
  
  
  ####ALIGNC
  my $align_bf = 0.0;
  my $align_b = sprintf "%5.3f", $align_bf;
  my $align_cf = 1.0;
  my $align_c = sprintf "%5.3f", $align_cf;
  
  system("make clean");
  my $cmd = "make tipo=2 term=2h angulo=$angulo bcmedio=$bcmedio abmedio=$abmedio align_b=$align_b align_c=$align_c m1=$m1 m2=$m2 $eje";
  system($cmd);
  
  for(my $k = 0; $k < $nrun; $k++)
  {
  	my $cmd = "$eje $k";
  	system($cmd);
  }
  
  system("make clean");
  $cmd = "make term=2h angulo=$angulo bcmedio=$bcmedio abmedio=$abmedio align_b=$align_b align_c=$align_c	m1=$m1 m2=$m2 reescribe.x";
  system($cmd);
  
  $cmd = "reescribe.x $nrun";
  system($cmd);

  $i++;
} until ($i == 6);

