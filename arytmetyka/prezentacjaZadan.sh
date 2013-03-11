#!/bin/sh
echo ''
echo '----------------------------'
echo ''

echo 'Porownanie kalkulacji wyrazow ciagu x{n+1}= x{n} + 3.0 * x{n} (1 - x{n})'
echo 'Podaj ile pierwszych wyrazow zaprezentowac'
read wyrazow 

gcc ciag.c -o ciag
./ciag $wyrazow
echo '(enter przerywa)'
gnuplot ciag.p
rm out.tmp

echo ''
echo '----------------------------'
echo ''

echo 'Prezentacja reprezentacji floatow' 
gcc reprezentacja.c -o rep -lgslcblas -lgsl -lm
./rep

echo ''
echo '----------------------------'
echo ''

echo 'Prezentacja rzutowania'
gcc jednatrzecia.c -o jt -lgslcblas -lgsl -lm
./jt

echo ''
echo '----------------------------'
echo ''


gcc mrs.c -o mrs -lgslcblas -lgsl -lm
echo 'Zaimlementowano rozwiazanie rownania rozniczkowego za pomoca mrs'
echo "cu'"' - u" = x^2 - rownanie dyfuzji-konwekcji'
echo 'x e(0,1); u(0)=u(1)=0;'
echo 'Podaj liczbe podzialow'
read num
echo 'Podaj stala konwekcji c'
read stala
./mrs $num $stala
gnuplot mrs.p
rm out_0.tmp
