#!/bin/sh
gcc ciag.c -o ciag
./ciag $1
gnuplot ciag.p
rm out.tmp
