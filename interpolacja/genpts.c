#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>


const int WIDTH = 800;
const int HEIGHT = 600;
double *ya,*xa;
int stepsBak;

void gen(int steps){
    srand(time(NULL));
    int x,y;
    double spare = WIDTH;
    x=0;
    xa = (double*) malloc(steps*sizeof(double));
    ya = (double*) malloc(steps*sizeof(double));
    stepsBak=steps;
    while(steps--){
        spare = WIDTH-x;
        y= rand()%(HEIGHT);
        int tmp = (spare/steps)+1;
        if(tmp<2)tmp=spare;
        x+= rand()%tmp+1;
        xa[stepsBak-steps-1]=x;
        ya[stepsBak-steps-1]=y;
        //printf("%f\t%f\n",xa[stepsBak-steps-1],ya[stepsBak-steps-1]);
    }
}
void gsl_polynomial(){
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_interp* inter = gsl_interp_alloc(gsl_interp_polynomial, stepsBak);
    gsl_interp_init (inter,xa,ya, stepsBak);

    double xi, yi;
    for (xi = xa[0]; xi < xa[stepsBak-1]; xi += 1.0 ) {
        yi = gsl_interp_eval(inter, xa, ya, xi, acc);
        printf("%lf %lf\n", xi, yi);
    }
    gsl_interp_free(inter);
    gsl_interp_accel_free(acc);
    free(xa);
    free(ya);
}
//////////////////////////////////////////////////////////////////////////////////
/*
stuct polynomial {
    int n;
    double *a;
};
*/
//////////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[]){
    int steps = -1;
    if(argc<2){
        printf("\nNie podano liczby punktow\n");
        return 1;
    } else {
        steps = atoi(argv[1]);
        if(steps<1){
            printf("\nLiczba punktow musi byc niemniejza niz 1\n");
            return 2;
        }
    }
    gen(steps);
    gsl_polynomial();
    return 0;
}
