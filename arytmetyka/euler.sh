#!/bin/sh
gcc euler.c -o euler
./euler $1
gnuplot euler.p
rm outEuler.tmp
