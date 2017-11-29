#!/bin/sh

eje='hm.cono'
nf=5

make clean
make

i=0

while [ $i -lt $nf ]
do
  $eje $i GA
  i=`expr $i + 1`
done

make clean
make term=2h

i=0
while [ $i -lt $nf ]
do
  $eje $i GA
  i=`expr $i + 1`
done
