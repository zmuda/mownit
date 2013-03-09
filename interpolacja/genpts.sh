#!/bin/sh
gcc genpts.c -o genpts
./genpts $1
#gnuplot genpts.p
#rm out.tmp
