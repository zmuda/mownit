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
    gsl_vector* diag = gsl_vector_alloc( N-2 );
    gsl_vector_set_all (diag, 2*n*n);
    gsl_vector* offdiag = gsl_vector_alloc( N-3 );
    gsl_vector_set_all (offdiag, (c/2)*n-n*n);
    gsl_vector* f = gsl_vector_alloc( N-2 );
    int i=N-2;
    double h = 1.0/n;
    while(i--){
 //     printf("%i",i);
        gsl_vector_set (f,i, (i*h+h) * (i*h+h) );
    }
    gsl_vector* x = gsl_vector_alloc( N-2 );
    //symetrycznie
    //gsl_vector* upper = gsl_vector_alloc( sizeof(double)*(N-3) );
/*
    int j=N-1;
    while(i--){
        printf();
    }
*/

    gsl_linalg_solve_symm_tridiag (diag,offdiag,f,x);

    FILE* out = fopen("out_0.tmp","w");
    i=N-2;
    fprintf(out,"%f\t%f\n",1.0,0.0);
    while(i--){
        fprintf(out,"%f\t%f\n",i*h+h,gsl_vector_get (x,i));
    }
    fprintf(out,"%f\t%f",0.0,0.0);

    gsl_vector_free(diag);
    gsl_vector_free(offdiag);
    //symetrycznie
    //gsl_vector_free (gsl_vector * upper);
    return 0;
}
