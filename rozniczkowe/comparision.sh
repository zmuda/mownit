#!/bin/sh
gcc comparision.c -o comparision -lm -lgsl -lgslcblas
./comparision $1 $2 $3
gnuplot comparision.p
rm outEuler.tmp
rm outGSL.tmp
