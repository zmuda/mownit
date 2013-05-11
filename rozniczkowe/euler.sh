#!/bin/sh
gcc euler.c -o euler -lm
./euler
gnuplot euler.p
rm outEuler.tmp
