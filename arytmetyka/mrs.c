#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_vector.h>

/*
nie ma takiej sily, ktora pozwoli domyslic sie, ze kiedys to byla 1/3
(przed przyporzadkowaniem do float'a)
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

    const unsigned N =n+1;
    gsl_vector* diag = gsl_vector_alloc( sizeof(double)*(N-2) );
    gsl_vector_set_all (diag, 2*n*n);
    gsl_vector* offdiag = gsl_vector_alloc( sizeof(double)*(N-3) );
    gsl_vector_set_all (offdiag, (c/2)*n-n*n);
    gsl_vector* f = gsl_vector_alloc( sizeof(double)*(N-2) );
    int i=N-2;
    double h = 1.0/n;
    while(i--){
        gsl_vector_set (f, sizeof(double)*i, i*h * i*h);
    }
    gsl_vector* x = gsl_vector_alloc( sizeof(double)*(N-2) );
    //symetrycznie
    //gsl_vector* upper = gsl_vector_alloc( sizeof(double)*(N-3) );

    gsl_linalg_solve_symm_tridiag (diag,offdiag,f,x);

    gsl_vector_free(diag);
    gsl_vector_free(offdiag);
    //symetrycznie
    //gsl_vector_free (gsl_vector * upper);
    return 0;
}
