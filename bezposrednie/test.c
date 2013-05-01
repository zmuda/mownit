//TEN KOD NIE JEST DO CZYTANIA
#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_linalg.h>
#include "mmio.h"
#include "mmio.c"

void displayVector(gsl_vector *x){
    int i=0;
    while(i<x->size){
        printf("%f\n",x->data[i]);
        i++;
    }
}
void LU(gsl_matrix_view m, gsl_vector_view b,int height){
    gsl_permutation * p = gsl_permutation_alloc (height);
    gsl_vector *x = gsl_vector_alloc (height);
    int s;
    gsl_linalg_LU_decomp (&m.matrix, p, &s);
    gsl_linalg_LU_solve (&m.matrix, p, &b.vector, x);
    gsl_permutation_free (p);
#ifdef VERBOSE
    displayVector(x);
#endif
    gsl_vector_free (x);
}
void Cholesky(gsl_matrix_view m, gsl_vector_view b,int height){
    gsl_vector *x = gsl_vector_alloc (height);
    int s;
    gsl_linalg_cholesky_decomp (&m.matrix);
    gsl_linalg_cholesky_solve (&m.matrix, &b.vector, x);
#ifdef VERBOSE
    displayVector(x);
#endif
    gsl_vector_free (x);
}
void gen(double* m, int width, int height) {
    srand(time(NULL));
    int i;
    for(i=0; i<width*height; ++i) {
        m[i]=((double)rand())/(rand()+1);
    }
}
void display(double* m,int width, int height) {
    int i;
    for(i=0; i<width*height; i++) {
        if(!(i%width))printf("\n");
        printf("%.2f\t",m[i]);

    }
    printf("\n");
}
double* readMM(char* name, int* width, int* height){

    FILE *f;
    MM_typecode matcode;
    if ((f = fopen(name, "r")) == NULL)
            exit(1);
    if (mm_read_banner(f, &matcode) != 0){
        printf("Could not process Matrix Market banner.\n");
        exit(2);
    }
    int nz;
    if ((mm_read_mtx_crd_size(f, height, width, &nz)) !=0)
        exit(3);
    int width2=*width;
    int height2=*height;
    double* ret= malloc(sizeof(double)*width2*height2);
    int x,y;
    double val;
    int i;
    for(i=0;i<width2*height2;i++)ret[i]=0;
    while(nz--){
        fscanf(f, "%d %d %lg\n", &y, &x, &val);
        ret[(y-1)*width2+(x-1)]=val;
    }
    return ret;
}
int main(int argc, char *argv[]) {
    int width,height;
    if(argc<2)return 99;
    double* example=readMM(argv[1],&width,&height);
    double* f=malloc(sizeof(double)*height);
    gen(f,height,1);
#ifdef VERBOSE
    display(f,1,height);
    printf("\n\nas a \"left side\" gives:\n\n");
#endif
    gsl_matrix_view AA = gsl_matrix_view_array(example, height, width);
    gsl_vector_view BB = gsl_vector_view_array(f, height);
    Cholesky(AA,BB,height);
#ifndef VERBOSE
    printf("DONE\n");
#endif
    return 0;
}
