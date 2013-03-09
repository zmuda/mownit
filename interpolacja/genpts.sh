#!/bin/sh
gcc genpts.c -o genpts -lgslcblas -lgsl -lm
./genpts $1
# >> out.tmp
#gnuplot genpts.p
#rm out.tmp
