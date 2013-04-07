#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit.h>
#include "parametry.h"
void mul(double a[], double b[], double res[], int width,int height) {
    int i,j,k;
    for(i=0; i<height; i++) {
        for(j=0; j<width; j++) {
            for(k=0; k<height; k++) {
                res[k+i*height]+=a[i*width+j]*b[k+j*height];
                //printf("(%d %d %d)\n",k+i*height,i*width+j,k+j*height);
            }
        }
    }
}
void mulNaive(double a[], double b[], double res[], int width,int height) {
    int i,j,k;
    for(k=0; k<height; k++) {
        for(j=0; j<height; j++) {
            for(i=0; i<width; i++) {
                res[j+k*height]+=a[i+k*width]*b[j+i*height];
                //printf("%d %d %d\n",j+k*height,i+k*width,j+i*height);
            }
        }
    }
}
void gen(double* m, int width, int height) {
    //srand(time(NULL));
    int i;
    for(i=0; i<width*height; ++i) {
        m[i]=((double)rand())/(rand()+1);
    }

}
void zeros(double* m, int width, int height) {
    int i;
    for(i=0; i<width*height; ++i) {
        m[i]=0;
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
    double *c,*b,*a;
    double times[TIMES];
    struct timeval stop, start;
    int l;

    int width = INIT;
    int height = INIT;
    c=malloc(sizeof(double)*height*height);
    b=malloc(sizeof(double)*width*height);
    a=malloc(sizeof(double)*width*height);
    gen(a,width,height);
    gen(b,height,width);
    zeros(c,height,height);

    int n =(1+(INIT-MIN)/STEP);
    gsl_matrix *X,*XB, *cov;
    gsl_vector *yN,*yB,*yU, *wN,*wU,*wB, *cU, *cB, *cN;
    X = gsl_matrix_alloc (n, 4);
    XB = gsl_matrix_alloc (n, 4);
    yN = gsl_vector_alloc (n);
    yU = gsl_vector_alloc (n);
    yB = gsl_vector_alloc (n);
    wN = gsl_vector_alloc (n);
    wU = gsl_vector_alloc (n);
    wB = gsl_vector_alloc (n);
#define CN(i) (gsl_vector_get(cN,(i)))
#define CU(i) (gsl_vector_get(cU,(i)))
#define CB(i) (gsl_vector_get(cB,(i)))
    cN = gsl_vector_alloc (4);
    cU = gsl_vector_alloc (4);
    cB = gsl_vector_alloc (4);
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
        gsl_matrix_view A = gsl_matrix_view_array(a, height, width);
        gsl_matrix_view B = gsl_matrix_view_array(b, width, height);
        gsl_matrix_view C = gsl_matrix_view_array(c, height, height);
        l=TIMES;
        while(l--) {
            gettimeofday(&start, NULL);
            mulNaive(a,b,c,width,height);
            gettimeofday(&stop, NULL);
            times[l]=(double)stop.tv_usec+(double)stop.tv_sec*1000000 - (double)start.tv_usec-(double)start.tv_sec*1000000;
        }
        yi=gsl_stats_mean(times,1,TIMES);
        ei=gsl_stats_sd(times,1,TIMES);
        gsl_vector_set (yN, i, yi);
        gsl_vector_set (wN, i, 1.0/(ei*ei));
#ifdef VERBOSE
        printf("Proste mnozenie srednio zajelo %f usec\n",yi);
        printf("\tOdchylenie standardowe %f usec\n",ei);
#endif
        fprintf(fd,"%f\t%f\t%f\t",(double)height,yi,ei);
        l=TIMES;
        while(l--) {
            gettimeofday(&start, NULL);
            mul(a,b,c,width,height);
            gettimeofday(&stop, NULL);
            times[l]=(double)stop.tv_usec+(double)stop.tv_sec*1000000 - (double)start.tv_usec-(double)start.tv_sec*1000000;
        }
        yi=gsl_stats_mean(times,1,TIMES);
        ei=gsl_stats_sd(times,1,TIMES);
        gsl_vector_set (yU, i, yi);
        gsl_vector_set (wU, i, 1.0/(ei*ei));
#ifdef VERBOSE
        printf("Ulepszone mnozenie srednio zajelo %f usec\n",yi);
        printf("\tOdchylenie standardowe %f usec\n",ei);
#endif
        fprintf(fd,"%f\t%f\t",yi,ei);
        l=TIMES;
        while(l--) {
            gettimeofday(&start, NULL);
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, &A.matrix, &B.matrix,0.0, &C.matrix);
            gettimeofday(&stop, NULL);
            times[l]=(double)stop.tv_usec+(double)stop.tv_sec*1000000 - (double)start.tv_usec-(double)start.tv_sec*1000000;
        }
        yi=gsl_stats_mean(times,1,TIMES);
        ei=gsl_stats_sd(times,1,TIMES);
        gsl_vector_set (yB, i, yi);
        gsl_vector_set (wB, i, 1.0/(ei*ei));
#ifdef VERBOSE
        printf("BLAS mnozenie srednio zajelo %f usec\n",yi);
        printf("\tOdchylenie standardowe %f usec\n",ei);
#endif
        fprintf(fd,"%f\t%f\n",yi,ei);
        width-=STEP;
        height-=STEP;
    }
    gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (n, 4);
    gsl_multifit_wlinear (X, wN, yN, cN, cov, &chisq, work);
    fprintf (fd2,"naive(x) = %g + %g*x + %g*x**2 + %g*x**3\n",
             CN(0), CN(1), CN(2), CN(3));
#ifdef VERBOSE
    printf ("# best fit: Y = %g + %g X + %g X^2 + %g X^3\n",
            CN(0), CN(1), CN(2), CN(3));
    print_stats(cov,chisq);
#endif
    gsl_multifit_wlinear (X, wU, yU, cU, cov, &chisq, work);
    fprintf (fd2,"improved(x)= %g + %g*x + %g*x**2 + %g*x**3\n",
             CU(0), CU(1), CU(2), CU(3));
#ifdef VERBOSE
    printf ("# best fit: Y = %g + %g X + %g X^2 + %g X^3\n",
            CU(0), CU(1), CU(2), CU(3));
    print_stats(cov,chisq);
#endif
    gsl_multifit_wlinear (XB, wB, yB, cB, cov, &chisq, work);
    gsl_multifit_linear_free (work);
    fprintf (fd2,"blas(x) = %g + %g*x**%g + %g*x**%g + %g*x**%g\n",
             CB(0), CB(1),COMPL/3, CB(2),COMPL*2/3, CB(3),COMPL);
#ifdef VERBOSE
    printf ("# best fit: Y = %g + %g X^%g + %g X^%g + %g X^%g\n",
            CB(0), CB(1),COMPL/3, CB(2),COMPL*2/3, CB(3),COMPL);
    print_stats(cov,chisq);
#endif
    fprintf(fd2,"plot 'tmp.out' using 1:2 title 'naive' linecolor rgb 'blue'");
    fprintf(fd2,",'tmp.out' using 1:4 title 'improved' linecolor rgb 'green'");
    fprintf(fd2,",'tmp.out' using 1:6 title 'blas' linecolor rgb 'red'");
    fprintf(fd2,",blas(x) with lines linecolor rgb 'red'");
    fprintf(fd2,", improved(x) with lines linecolor rgb 'green'");
    fprintf(fd2,",naive(x) with lines linecolor rgb 'blue'");
    fprintf(fd2,"\n\npause -1\nset terminal svg\nset output\"output.svg\"\nreplot");

    fclose(fd);
    fclose(fd2);
    free(a);
    free(b);
    free(c);
    gsl_matrix_free (X);
    gsl_matrix_free (XB);
    gsl_vector_free (yN);
    gsl_vector_free (yU);
    gsl_vector_free (yB);
    gsl_vector_free (wN);
    gsl_vector_free (wB);
    gsl_vector_free (wU);
    gsl_vector_free (cB);
    gsl_vector_free (cN);
    gsl_vector_free (cU);
    gsl_matrix_free (cov);
    return 0;
}
