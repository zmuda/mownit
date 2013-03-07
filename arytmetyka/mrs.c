#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_vector.h>

/*
# mrs jest niestabilna numerycznie
# rozwiązanie przykładowego problemu
# cu' - u" = x^2 - rownanie dyfuzji-konwekcji
# xe(0,1); u(0)=u(1)=0;
# arg1 = liczba podzialow; arg2 = stala konwekcji c
*/
int main(int argc, char *argv[]){
    unsigned n;
    double c;
    if(argc<3){
        printf("\nWymagane dwa argumenty: liczba podzialow; stala konwekcji c\n");
        return 1;
    } else {
        int i=atoi(argv[1]);
        c=atof(argv[2]);
        if( i<2 ){
            printf("\nNie mozna dzielic na mniej niz dwa elementy\n");
            return 2;
        }
        n=i;
    }

    /*
    macierz wspolczynników w mrs jest trójdiagonalna
    przy obranych warunkach poczatkowych jest ponadto
     symetryczna wzgledem obu diagonal
    */
    const unsigned N =n+1;
    gsl_vector* diag = gsl_vector_alloc( N-2 );
    //analitycznie uzyskana wartosc
    gsl_vector_set_all (diag, 2*n*n);
    gsl_vector* offdiag = gsl_vector_alloc( N-3 );
    //analitycznie uzyskana wartosc
    gsl_vector_set_all (offdiag, (c/2)*n-n*n);
    gsl_vector* f = gsl_vector_alloc( N-2 );
    //wypelnienie probkowanymi wartosciami funkcji x^2
    int i=N-2;
    double h = 1.0/n;
    while(i--){
        gsl_vector_set (f,i, (i*h+h) * (i*h+h) );
    }
    gsl_vector* x = gsl_vector_alloc( N-2 );

    //solver dedykowany do trojdiagonalej symetrycznej
    //macierzy wspolczynnikow
    // - mniejsza zlozonosc O(n^2)
    gsl_linalg_solve_symm_tridiag (diag,offdiag,f,x);

    //plik tymczasowy dla gnuplot
    FILE* out = fopen("out_0.tmp","w");
    i=N-2;
    //warunek brzegowy #1
    fprintf(out,"%f\t%f\n",1.0,0.0);
    while(i--){
        fprintf(out,"%f\t%f\n",i*h+h,gsl_vector_get (x,i));
    }
    //warunek brzegowy #2
    fprintf(out,"%f\t%f",0.0,0.0);

    gsl_vector_free(diag);
    gsl_vector_free(offdiag);
    return 0;
}
