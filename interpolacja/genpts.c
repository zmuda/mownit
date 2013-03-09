#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>


void work(int steps){
    const int WIDTH = 800;
    const int HEIGHT = 600;

    FILE* out = fopen("out.tmp","w");
    srand(time(NULL));
    int x,y;
    double spare = WIDTH;
    x=0;
    while(steps--){
        spare = WIDTH-x;
        y= rand()%(HEIGHT);
        int tmp = (spare/steps);
        if(tmp<2)tmp=2;
        x+= rand()%tmp;
        fprintf(out,"%i\t%i\n",x,y);
    }
    //fprintf(out,"%f\t%i\t%lf\n",xf,i,xd);
    fclose(out);
}


void gen(int steps){
    const int WIDTH = 800;
    const int HEIGHT = 600;

    FILE* out = fopen("out.tmp","w");
    srand(time(NULL));
    int x,y;
    double spare = WIDTH;
    x=0;
    double* xa = (double*) malloc(steps*sizeof(double));
    double* ya = (double*) malloc(steps*sizeof(double));
    int stepsBak=steps;
    while(steps--){
        spare = WIDTH-x;
        y= rand()%(HEIGHT);
        int tmp = (spare/steps)+1;
        if(tmp<2)tmp=2;
        x+= rand()%tmp;
        xa[stepsBak-steps-1]=x;
        ya[stepsBak-steps-1]=y;
        printf("%f\t%f\n",xa[stepsBak-steps-1],ya[stepsBak-steps-1]);
    }
}



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
    return 0;
}
