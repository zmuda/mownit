CC=gcc
CFLAGS=-lm -lgsl -lgslcblas
E=echo

all: testing benching

test: test.c
	$(CC) $@.c -o $@.out $(CFLAGS) 
bench: bench.c
	$(CC) $@.c -o $@.out $(CFLAGS) 

testing: test
	$E "cylindryczny pocisk obracajaca sie w jednorodnym osrodku - symulacja - dyskretyzacja dziedziny"
	$E "spelnia wymagania rozkladu Cholesky'ego"
	./test.out s3rmt3m3.mtx
	$E "wyboczanie(sciskanie->zginanie) wlotu slnika Boeing'a"
	$E "spelnia wymagania rozkladu Cholesky'ego"
	./test.out bcsstk27.mtx
	$E "uogolnione wartosci wlasne przeplywu plynow"
	$E "spelnia wymagania rozkladu Cholesky'ego"
	./test.out bcsstk13.mtx

benching: bench
	./bench.out
	gnuplot plot.p

clean:
	rm -rf *.out
	rm -rf *.o
	rm -rf *.tmp
	rm -rf *~
