#!/bin/sh
gcc mrs.c -o mrs -lgslcblas -lgsl -lm
# cu' - u" = x^2 - rownanie dyfuzji-konwekcji
# xE(0,1); u(0)=u(1)=0;
# $1 = liczba podzialow; $2 = stala konwekcji c
./mrs $1 $2
echo TODOgnuplot
