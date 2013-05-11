#!/bin/sh
gcc animation.c -o animation -lm -lgsl -lgslcblas
./animation $1 $2 $3
gnuplot animation.gpi
