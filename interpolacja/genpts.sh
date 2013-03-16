#!/bin/sh
gcc genpts.c -o genpts -lgslcblas -lgsl -lm -O2
./genpts $1
gnuplot genpts.p
rm out.tmp
rm outPts.tmp
