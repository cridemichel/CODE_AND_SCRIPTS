#!/bin/sh
rm  eg_vs_d2eg.dat
for f in Tgamma*
do
stddevtail+mean.sh $f/'rcmz.dat_R*' 100 >> eg_vs_d2eg.dat
done
