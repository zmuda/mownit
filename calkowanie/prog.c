#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_statistics.h>
#include "funkcja1.h"
#define f(x) f(x,NULL)
typedef double funkcyjny( double x, void* p);

double trapestry(funkcyjny f,double begin, double end, unsigned probes) {
    if(probes<2)exit(2);
    double step = (end-begin)/(probes-1);
    double ret=0;
    int i;
    double x=begin;
    for(i=0; i<probes-1; i++) {
        ret+=(f(x)+f(x+step))/2*step;
        x+=step;
    }
    return ret;
}
double trapestryNonfixed(funkcyjny f,double begin, double end, double estimatedError, int* stepsTaken) {
    double last = trapestry(f,begin,end,2);
    double res = trapestry(f,begin,end,3);
    int i=3;
    while(fabs(last-res)>estimatedError) {
        last=res;
        res=trapestry(f,begin,end,++i);
    }
    *stepsTaken=i;
    return res;
}

struct timeval stop, start;
int main(int argc, char *argv[]) {
    double ERROR=1e-12;
    int steps;
    double res=trapestryNonfixed(f,0,1,ERROR,&steps);
    printf("Dla tolerancji:\t %.36lf\n",ERROR);
    printf("Trapezy:\t%.36lf\tkroków: %d\n",res,steps);

    double times[10];
    int l = 10;
    while(l--) {
        gettimeofday(&start, NULL);
        trapestry(f,0,1,steps);
        gettimeofday(&stop, NULL);
        times[l]=(double)(stop.tv_usec-start.tv_usec);
    }
    printf("średnio %f usec\n",gsl_stats_mean(times,1,10));

    gsl_function F;
    F.function = &f;
    size_t neval;
    double result, error;
    gsl_integration_qng(&F, 0, 1, 0, ERROR, &result, &error, &neval);
    printf("gsl(QNG):\t%.36lf\tkroków: %d\n",result,neval);

    l=10;
    while(l--) {
        gettimeofday(&start, NULL);
        gsl_integration_qng(&F, 0, 1, ERROR, 0, &result, &error, &neval);
        gettimeofday(&stop, NULL);
        times[l]=(double)(stop.tv_usec-start.tv_usec);
    }
    l=10;
    printf("średnio %f usec\n",gsl_stats_mean(times,1,10));
    printf("wolfram\t\t%f",-2.0);

    ERROR=1e-4;
    printf("\n\t--\t--\t--\t--\n\n");
    //QAG
    F.function = &k;
    printf("ADAPTACYJNIE\n");
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (POINTS);
    gsl_integration_qag(&F, 0.01, 100, ERROR, 0, POINTS, 1, w,&result, &error);
    printf("Dla ilosci krokow:\t%d\n",POINTS);
    printf("gsl(QAG):\t%.36lf\tkroków: %d\n",result,POINTS);
    printf("trapestry:\t%.36lf\tkroków: %d\n",trapestry(k,0.01,100,POINTS),POINTS);
    printf("wolfram\t\t%f",1.54838);

    printf("\n\t--\t--\t--\t--\n\n");
    //QAGS
    F.function = &h;
    printf("OSOBLIWOŚCI\n");
    gsl_integration_qags(&F, 0, 1, ERROR, 0, POINTS, w,&result, &error);
    printf("Dla ilosci krokow:\t%d\n",POINTS);
    printf("gsl(QAGS):\t%.36lf\tkroków: %d\n",result,POINTS);
    printf("trapestry:\t%.36lf\tkroków: %d\n",trapestry(h,0,1,POINTS),POINTS);
    printf("wolfram\t\t%f",-4.0);

    printf("\n\t--\t--\t--\t--\n\n");
    //QAWO
    F.function = &g;
    printf("OSCYLACYJNE\n");
    gsl_integration_qawo_table * wf = gsl_integration_qawo_table_alloc(1,3.14,GSL_INTEG_COSINE,POINTS);
    gsl_integration_qawo(&F,0.0,ERROR,ERROR,POINTS,w,wf,&result,&error);
    printf("Dla ilosci krokow:\t%d\n",POINTS);
    printf("gsl(QAGS):\t%.36lf\tkroków: %d\n",result,POINTS);
    printf("trapestry:\t%.36lf\tkroków: %d\n",trapestry(g,0,3.14,POINTS),POINTS);
    printf("wolfram\t\t%f",-.104717);

    gsl_integration_qawo_table_free(wf);
    gsl_integration_workspace_free(w);
    return 0;
}
