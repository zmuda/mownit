//TEN KOD NIE JEST DO CZYTANIA
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_linalg.h>
#include "parametry.h"
void LU(gsl_matrix_view m, gsl_vector_view b,int height){
    gsl_permutation * p = gsl_permutation_alloc (height);
    gsl_vector *x = gsl_vector_alloc (height);
    int s;
    gsl_linalg_LU_decomp (&m.matrix, p, &s);
    gsl_linalg_LU_solve (&m.matrix, p, &b.vector, x);
    gsl_permutation_free (p);
    gsl_vector_free (x);
}
void Cholesky(gsl_matrix_view m, gsl_vector_view b,int height){
    gsl_vector *x = gsl_vector_alloc (height);
    int s;
    gsl_linalg_cholesky_decomp (&m.matrix);
    gsl_linalg_cholesky_solve (&m.matrix, &b.vector, x);
    gsl_vector_free (x);
}
void gensym(double* m, int dim) {
    srand(time(NULL));
    int i,j;
    for(i=0;i<dim;i++){
        for(j=i;j<dim;j++){
            m[i*dim+j]=m[j*dim+i]=((double)rand())/(rand()+1);
        }
    }
}
void genPositiveDefinite(double* m, int dim) {
    srand(time(NULL));
    int i,j;
    for(i=0;i<dim;i++){
        for(j=i;j<dim;j++){
            m[i*dim+j]=m[j*dim+i]=1.0/(i+j+1);
            if(i==j)m[i*dim+j]+=1;
        }
    }
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
void print_stats(gsl_matrix *cov,double chisq) {
#define COV(i,j) (gsl_matrix_get(cov,(i),(j)))
    printf ("# covariance matrix:\n");
    printf ("[ %+.5e, %+.5e, %+.5e, %+.5e  \n",
            COV(0,0), COV(0,1), COV(0,2), COV(0,3));
    printf ("  %+.5e, %+.5e, %+.5e, %+.5e\n",
            COV(1,0), COV(1,1), COV(1,2), COV(1,3));
    printf ("  %+.5e, %+.5e, %+.5e, %+.5e\n",
            COV(2,0), COV(2,1), COV(2,2), COV(2,3));
    printf ("  %+.5e, %+.5e, %+.5e, %+.5e ]",
            COV(3,0), COV(3,1), COV(3,2), COV(3,3));
    printf ("# chisq = %g\n", chisq);
}
//////////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[]) {
    double times[TIMES];
    struct timeval stop, start;
    int l;

    int width = INIT;
    int height = INIT;

    double *adata, *bdata;
    adata=malloc(sizeof(double)*width*height);
    bdata=malloc(sizeof(double)*1*height);
    genPositiveDefinite(adata,width);
    gen(bdata,height,1);

    int n =(1+(INIT-MIN)/STEP);
    gsl_matrix *X,*XB, *cov;
    gsl_vector *yN,*yU, *wN,*wU, *cU, *cN;
    X = gsl_matrix_alloc (n, 4);
    XB = gsl_matrix_alloc (n, 4);
    yN = gsl_vector_alloc (n);
    yU = gsl_vector_alloc (n);
    wN = gsl_vector_alloc (n);
    wU = gsl_vector_alloc (n);
#define CN(i) (gsl_vector_get(cN,(i)))
#define CU(i) (gsl_vector_get(cU,(i)))
    cN = gsl_vector_alloc (4);
    cU = gsl_vector_alloc (4);
    cov = gsl_matrix_alloc (4, 4);
    FILE* fd = fopen("tmp.out","w");
    FILE* fd2 = fopen("plot.p","w");
    int i=n;
    double xi, yi, ei, chisq;
    while(width>MIN) {
        i--;
#ifdef VERBOSE
        printf("\t<%dx%d>\n",width,height);
#endif
        gsl_matrix_set (X, i, 0, 1.0);
        gsl_matrix_set (X, i, 1, width);
        gsl_matrix_set (X, i, 2, width*width);
        gsl_matrix_set (X, i, 3, width*width*width);
        gsl_matrix_set (XB, i, 0, 1.0);
        gsl_matrix_set (XB, i, 1, pow(width,COMPL/3));
        gsl_matrix_set (XB, i, 2, pow(width,COMPL*2/3));
        gsl_matrix_set (XB, i, 3, pow(width,COMPL));
        genPositiveDefinite(adata,width);
        gsl_matrix_view AA = gsl_matrix_view_array(adata, height, width);
        gsl_vector_view BB = gsl_vector_view_array(bdata, height);


        l=TIMES;
        while(l--) {
            gettimeofday(&start, NULL);
            LU(AA,BB,height);
            gettimeofday(&stop, NULL);
            times[l]=(double)stop.tv_usec+(double)stop.tv_sec*1000000 - (double)start.tv_usec-(double)start.tv_sec*1000000;
        }
        yi=gsl_stats_mean(times,1,TIMES);
        ei=gsl_stats_sd(times,1,TIMES);
        gsl_vector_set (yN, i, yi);
        gsl_vector_set (wN, i, 1.0/(ei*ei));
#ifdef VERBOSE
        printf("LU: srednio zajelo %f usec\n",yi);
        printf("\tOdchylenie standardowe %f usec\n",ei);
#endif
        fprintf(fd,"%f\t%f\t%f\t",(double)height,yi,ei);
        l=TIMES;
        while(l--) {
            gettimeofday(&start, NULL);
            Cholesky(AA,BB,height);
            gettimeofday(&stop, NULL);
            times[l]=(double)stop.tv_usec+(double)stop.tv_sec*1000000 - (double)start.tv_usec-(double)start.tv_sec*1000000;
        }
        yi=gsl_stats_mean(times,1,TIMES);
        ei=gsl_stats_sd(times,1,TIMES);
        gsl_vector_set (yU, i, yi);
        gsl_vector_set (wU, i, 1.0/(ei*ei));
#ifdef VERBOSE
        printf("Cholesky: srednio zajelo %f usec\n",yi);
        printf("\tOdchylenie standardowe %f usec\n",ei);
#endif
        fprintf(fd,"%f\t%f\n",yi,ei);
        width-=STEP;
        height-=STEP;
    }
    gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (n, 4);
    gsl_multifit_wlinear (X, wN, yN, cN, cov, &chisq, work);
    fprintf (fd2,"LU(x) = %g + %g*x + %g*x**2 + %g*x**3\n",
             CN(0), CN(1), CN(2), CN(3));
#ifdef VERBOSE
    printf ("# best fit: Y = %g + %g X + %g X^2 + %g X^3\n",
            CN(0), CN(1), CN(2), CN(3));
    print_stats(cov,chisq);
#endif
    gsl_multifit_wlinear (X, wU, yU, cU, cov, &chisq, work);
    fprintf (fd2,"Cholesky(x)= %g + %g*x + %g*x**2 + %g*x**3\n",
             CU(0), CU(1), CU(2), CU(3));
#ifdef VERBOSE
    printf ("# best fit: Y = %g + %g X + %g X^2 + %g X^3\n",
            CU(0), CU(1), CU(2), CU(3));
    print_stats(cov,chisq);
#endif
    fprintf(fd2,"plot 'tmp.out' using 1:2:3 title 'LU' linecolor rgb 'blue' with yerrorbars");
    fprintf(fd2,",'tmp.out' using 1:4:5 title 'Cholesky' linecolor rgb 'green' with yerrorbars");
    fprintf(fd2,", Cholesky(x) with lines linecolor rgb 'green'");
    fprintf(fd2,",LU(x) with lines linecolor rgb 'blue'");
    fprintf(fd2,"\n\npause -1\nset terminal svg\nset output\"output.svg\"\nreplot");

    fclose(fd);
    fclose(fd2);
    gsl_matrix_free (X);
    gsl_matrix_free (XB);
    gsl_vector_free (yN);
    gsl_vector_free (yU);
    gsl_vector_free (wN);
    gsl_vector_free (wU);
    gsl_vector_free (cN);
    gsl_vector_free (cU);
    gsl_matrix_free (cov);
    return 0;
}
