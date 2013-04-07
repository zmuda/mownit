#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit.h>

//#define WIDTH 300
//#define HEIGHT 300
#define TIMES 10


void mul(double a[], double b[], double res[], int width,int height){
    int i,j,k;
    for(i=0;i<height;i++){
        for(j=0;j<width;j++){
            for(k=0;k<height;k++){
                res[k+i*height]+=a[i*width+j]*b[k+j*height];
                //printf("(%d %d %d)\n",k+i*height,i*width+j,k+j*height);
            }
        }
    }
}
void mulNaive(double a[], double b[], double res[], int width,int height){
    int i,j,k;
    for(k=0;k<height;k++){
        for(j=0;j<height;j++){
            for(i=0;i<width;i++){
                res[j+k*height]+=a[i+k*width]*b[j+i*height];
                //printf("%d %d %d\n",j+k*height,i+k*width,j+i*height);
            }
        }
    }
}

void gen(double* m, int width, int height){
    srand(time(NULL));
    int i;
    for(i=0;i<width*height;++i){
        m[i]=((double)rand())/(rand()+1);
    }

}
void zeros(double* m, int width, int height){
    int i;
    for(i=0;i<width*height;++i){
        m[i]=0;
    }

}
void display(double* m,int width, int height){
    int i;
    for(i=0;i<width*height;i++){
        if(!(i%width))printf("\n");
        printf("%.2f\t",m[i]);

    }
    printf("\n");
}

//////////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[]){
    double *c,*b,*a;
    double times[TIMES];
    struct timeval stop, start;
    int l;

#define INIT 400
#define STEP 12
#define MIN 12
    int WIDTH = INIT;
    int HEIGHT = INIT;
    c=malloc(sizeof(double)*HEIGHT*HEIGHT);
    b=malloc(sizeof(double)*WIDTH*HEIGHT);
    a=malloc(sizeof(double)*WIDTH*HEIGHT);
    gen(a,WIDTH,HEIGHT);
    gen(b,HEIGHT,WIDTH);
    zeros(c,HEIGHT,HEIGHT);


/*double *yN,*yU,*yB,*x,*eN,*eU,*eB;
yN=malloc(sizeof(double)*(1+(INIT-MIN)/35));
yU=malloc(sizeof(double)*(1+(INIT-MIN)/35));
yB=malloc(sizeof(double)*(1+(INIT-MIN)/35));
x=malloc(sizeof(double)*(1+(INIT-MIN)/35));
eN=malloc(sizeof(double)*(1+(INIT-MIN)/35));
eU=malloc(sizeof(double)*(1+(INIT-MIN)/35));
eB=malloc(sizeof(double)*(1+(INIT-MIN)/35));
for(l=0;l<(1+(INIT-MIN)/35);l++){
    x[l]=yN[l]=yU[l]=yB[l]=eN[l]=eU[l]=eB[l]=0;
}*/

    int n =(1+(INIT-MIN)/STEP);
    gsl_matrix *X, *cov;
    gsl_vector *yN,*yB,*yU, *wN,*wU,*wB, *cU, *cB, *cN;
    X = gsl_matrix_alloc (n, 4);
    yN = gsl_vector_alloc (n);
    yU = gsl_vector_alloc (n);
    yB = gsl_vector_alloc (n);
    wN = gsl_vector_alloc (n);
    wU = gsl_vector_alloc (n);
    wB = gsl_vector_alloc (n);
    #define CN(i) (gsl_vector_get(cN,(i)))
    #define CU(i) (gsl_vector_get(cU,(i)))
    #define CB(i) (gsl_vector_get(cB,(i)))
    #define COV(i,j) (gsl_matrix_get(cov,(i),(j)))
    cN = gsl_vector_alloc (4);
    cU = gsl_vector_alloc (4);
    cB = gsl_vector_alloc (4);
    cov = gsl_matrix_alloc (4, 4);

    FILE* fd = fopen("tmp.out","w");
    int i=n;
    double xi, yi, ei, chisq;
    while(WIDTH>MIN){
        i--;
        gsl_matrix_set (X, i, 0, 1.0);
        gsl_matrix_set (X, i, 1, WIDTH);
        gsl_matrix_set (X, i, 2, WIDTH*WIDTH);
        gsl_matrix_set (X, i, 3, WIDTH*WIDTH*WIDTH);

        gsl_matrix_view A = gsl_matrix_view_array(a, HEIGHT, WIDTH);
        gsl_matrix_view B = gsl_matrix_view_array(b, WIDTH, HEIGHT);
        gsl_matrix_view C = gsl_matrix_view_array(c, HEIGHT, HEIGHT);

        l=TIMES;
        while(l--){
            gettimeofday(&start, NULL);
            mulNaive(a,b,c,WIDTH,HEIGHT);
            gettimeofday(&stop, NULL);
            times[l]=(double)stop.tv_usec+(double)stop.tv_sec*1000000 - (double)start.tv_usec-(double)start.tv_sec*1000000;
            //printf("zajela %.1f usec\n",times[l]);
        }
        yi=gsl_stats_mean(times,1,TIMES);
        ei=gsl_stats_sd(times,1,TIMES);
        gsl_vector_set (yN, i, yi);
        gsl_vector_set (wN, i, 1.0/(ei*ei));
        printf("Proste mnozenie srednio zajelo %f usec\n",yi);
        printf("\tOdchylenie standardowe %f usec\n",ei);
        fprintf(fd,"%f\t%f\t%f\t",(double)HEIGHT,yi,ei);


        l=TIMES;
        while(l--){
            gettimeofday(&start, NULL);
            mul(a,b,c,WIDTH,HEIGHT);
            gettimeofday(&stop, NULL);
            times[l]=(double)stop.tv_usec+(double)stop.tv_sec*1000000 - (double)start.tv_usec-(double)start.tv_sec*1000000;
            //printf("zajela %.1f usec\n",times[l]);
        }
        yi=gsl_stats_mean(times,1,TIMES);
        ei=gsl_stats_sd(times,1,TIMES);
        gsl_vector_set (yU, i, yi);
        gsl_vector_set (wU, i, 1.0/(ei*ei));
        printf("Ulepszone mnozenie srednio zajelo %f usec\n",yi);
        printf("\tOdchylenie standardowe %f usec\n",ei);
        fprintf(fd,"%f\t%f\t",yi,ei);

        l=TIMES;
        while(l--){
            gettimeofday(&start, NULL);
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, &A.matrix, &B.matrix,0.0, &C.matrix);
            gettimeofday(&stop, NULL);
            times[l]=(double)stop.tv_usec+(double)stop.tv_sec*1000000 - (double)start.tv_usec-(double)start.tv_sec*1000000;
            //printf("zajela %.1f usec\n",times[l]);
        }
        yi=gsl_stats_mean(times,1,TIMES);
        ei=gsl_stats_sd(times,1,TIMES);
        gsl_vector_set (yB, i, yi);
        gsl_vector_set (wB, i, 1.0/(ei*ei));
        printf("BLAS mnozenie srednio zajelo %f usec\n",yi);
        printf("\tOdchylenie standardowe %f usec\n",ei);
        fprintf(fd,"%f\t%f\n",yi,ei);

        WIDTH-=STEP;
        HEIGHT-=STEP;
    }
    {
         gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (n, 4);
         gsl_multifit_wlinear (X, wN, yN, cN, cov, &chisq, work);
         gsl_multifit_linear_free (work);
    }
    {
         printf ("# best fit: Y = %g + %g X + %g X^2 + %g X^3\n",
                 CN(0), CN(1), CN(2), CN(3));

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
    {
         gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (n, 4);
         gsl_multifit_wlinear (X, wU, yU, cU, cov, &chisq, work);
         gsl_multifit_linear_free (work);
    }
    {
         printf ("# best fit: Y = %g + %g X + %g X^2 + %g X^3\n",
                 CU(0), CU(1), CU(2), CU(3));

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
        {
         gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (n, 4);
         gsl_multifit_wlinear (X, wB, yB, cB, cov, &chisq, work);
         gsl_multifit_linear_free (work);
    }
    {
         printf ("# best fit: Y = %g + %g X + %g X^2 + %g X^3\n",
                 CB(0), CB(1), CB(2), CB(3));

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
    return 0;
}
