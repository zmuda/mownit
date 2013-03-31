#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_blas.h>


void mul(double a[], double b[], double res[], int width,int height){
    int i,j,k;
    for(i=0;i<height;i++){
        for(j=0;j<width;j++){
            for(k=0;k<height;k++){
                res[k+i*height]+=a[i*width+j]*b[k+j*height];
            }
        }
    }
    /**/
}
void kick(double* m, int width,int height){
    double* tmp = malloc(sizeof(double)*width*height);
    int i,j;
    for(i=0;i<height;i++){
        for(j=0;j<width;j++){
            tmp[i*width+j]=m[j*height+i];
        }
    }
    memcpy(m,tmp,sizeof(double)*width*height);
}

void mul_upgraded2(double a[], double b[], double res[], int width,int height){
    //kick(b,width,height);
    int i,j,k;
    for(i=0;i<height;i++){
        for(j=0;j<height;j++){
            for(k=0;k<width;k++){
                //res[i][j]+=a[i][k]*b[k][j];//ord
                //res[i][j]+=a[i][k]*b[j][k];//kicked
                res[i*height+j]+=a[i*width+k]*b[j*width+k];
            }
        }
    }
}

void mul_upgraded(double a[], double b[], double res[], int width,int height){
    //kick(b,width,height);
    int i,j,k;
    for(i=0;i<width;i++){
        for(j=0;j<width;j++){
            for(k=0;k<height;k++){
                //res[i][j]+=a[i][k]*b[k][j];//ord
                //res[i][j]+=a[i][k]*b[j][k];//kicked
                res[i*width+j]+=a[i*height+k]*b[j*height+k];
            }
        }
    }
}

#define WIDTH 5
#define HEIGHT 4

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

}
//////////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[]){

/*
    double a[] = { 0.11, 0.12, 0.13,
                      0.21, 0.22, 0.23 };

    double b[] = { 1011, 1012,
                      1021, 1022,
                      1031, 1032 };

    double c[] = { 0.00, 0.00,
                      0.00, 0.00 };
/**/
    double *c,*b,*a;

    c=malloc(sizeof(double)*HEIGHT*HEIGHT);
    b=malloc(sizeof(double)*WIDTH*HEIGHT);
    a=malloc(sizeof(double)*WIDTH*HEIGHT);
    gen(a,WIDTH,HEIGHT);
    gen(b,HEIGHT,WIDTH);
    zeros(c,HEIGHT,HEIGHT);

/*
    mul(a,b,c,WIDTH,HEIGHT);
    display(c,HEIGHT,HEIGHT);
    zeros(c,HEIGHT,HEIGHT);
    mul_upgraded(a,b,c,WIDTH,HEIGHT);
    display(c,HEIGHT,HEIGHT);
*/


    gsl_matrix_view A = gsl_matrix_view_array(a, HEIGHT, WIDTH);
    gsl_matrix_view B = gsl_matrix_view_array(b, WIDTH, HEIGHT);
    gsl_matrix_view C = gsl_matrix_view_array(c, HEIGHT, HEIGHT);

    double times[10];
    struct timeval stop, start;
    int l;


    l=10;
    while(l--){
        gettimeofday(&start, NULL);
        mul_upgraded(a,b,c,WIDTH,HEIGHT);
        gettimeofday(&stop, NULL);
        times[l]=(double)stop.tv_usec+(double)stop.tv_sec*1000000 - (double)start.tv_usec-(double)start.tv_sec*1000000;
        printf("zajela %.1f usec\n",times[l]);
    }
    printf("Ulepszone mnozenie srednio zajelo %f usec\n",gsl_stats_mean(times,1,10));
    printf("\tOdchylenie standardowe %f usec\n",gsl_stats_sd(times,1,10));

    display(c,HEIGHT,HEIGHT);
    zeros(c,HEIGHT,HEIGHT);

    kick(b,WIDTH,HEIGHT);
    //sleep(1);
    l=10;
    while(l--){
        gettimeofday(&start, NULL);
        mul(a,b,c,WIDTH,HEIGHT);
        gettimeofday(&stop, NULL);
        times[l]=(double)stop.tv_usec+(double)stop.tv_sec*1000000 - (double)start.tv_usec-(double)start.tv_sec*1000000;
        printf("zajela %.1f usec\n",times[l]);
    }
    printf("Proste mnozenie srednio zajelo %f usec\n",gsl_stats_mean(times,1,10));
    printf("\tOdchylenie standardowe %f usec\n",gsl_stats_sd(times,1,10));

    //sleep(1);
    display(c,HEIGHT,HEIGHT);



    //sleep(1);
    l=10;
    while(l--){
        gettimeofday(&start, NULL);
        gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, &A.matrix, &B.matrix,0.0, &C.matrix);
        gettimeofday(&stop, NULL);
        times[l]=(double)stop.tv_usec+(double)stop.tv_sec*1000000 - (double)start.tv_usec-(double)start.tv_sec*1000000;
        printf("zajela %.1f usec\n",times[l]);
    }
    printf("BLAS mnozenie srednio zajelo %f usec\n",gsl_stats_mean(times,1,10));
    printf("\tOdchylenie standardowe %f usec\n",gsl_stats_sd(times,1,10));


/*
 Function: int gsl_blas_dgemm (CBLAS_TRANSPOSE_t TransA,
                               CBLAS_TRANSPOSE_t TransB,
                               double alpha, const gsl_matrix * A,
                               const gsl_matrix * B, double beta,
                               gsl_matrix * C)
 */
    //gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, &A.matrix, &B.matrix,0.0, &C.matrix);
    //mul(a,b,c,3,2);

    //printf ("[ %g, %g\n", c[0], c[1]);
    //printf ("  %g, %g ]\n", c[2], c[3]);

/*
    l=10;
    while(l--){
        len =lenBak;

        gettimeofday(&start, NULL);
        while(len--){
            nwtn[len]=polynomial_eval(newtonInterp,x[len]);
        }
        gettimeofday(&stop, NULL);
        //printf("Mojej met. Newton'a ewaluacja zajęła %lu usec\n", stop.tv_usec - start.tv_usec);
        times[l]=stop.tv_usec - start.tv_usec;
    }
    printf("Mojej met. Newton'a ewaluacja zajęła średnio %f usec\n",gsl_stats_mean(times,1,10));
    printf("\tOdchylenie standardowe %f usec\n",gsl_stats_sd(times,1,10));
*/

    /**/

    return 0;
}
