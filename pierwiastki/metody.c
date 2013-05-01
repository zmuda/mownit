#include<stdio.h>
#include<sys/times.h>
#include<gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include<unistd.h>

#define ELO -1
#define EPS 1e-15
#define TIMES 1000000


void pr_time(clock_t real, struct tms *start, struct tms *end)
{
    static long clktck = 0;
    if(clktck == 0)
        clktck = sysconf(_SC_CLK_TCK);
    fprintf(stderr,"\treal: %7.2f\n",real/(double)clktck);
    fprintf(stderr,"\tuser: %7.2f\n",(end->tms_utime-start->tms_utime)/(double)clktck);
    fprintf(stderr,"\tsys: %7.2f\n",(end->tms_stime-start->tms_stime)/(double)clktck);
}
int do_bucket(gsl_root_fsolver* s, int max_it, double* value){
    int it = 0, status;
    double x_lo, x_hi;
    double r;
    do{
        it++;
        gsl_root_fsolver_iterate (s);
        r = gsl_root_fsolver_root (s);
        x_lo = gsl_root_fsolver_x_lower (s);
        x_hi = gsl_root_fsolver_x_upper (s);
        status = gsl_root_test_interval (x_lo,x_hi,EPS,0);
    }while (status == GSL_CONTINUE && it < max_it);
    *value = r;
    return it;
}
int do_differential(gsl_root_fdfsolver* s, int max_it, double x_st, double* value){
    int it = 0, status;
    double x0, x = x_st;
    do{
        it++;
        gsl_root_fdfsolver_iterate (s);
        x0 = x;
        x = gsl_root_fdfsolver_root (s);
        status = gsl_root_test_delta (x,x0,EPS,0);
    }while (status == GSL_CONTINUE && it < max_it);
    *value=x;
    return it;
}
#define b1 12.0
#define c1 3.0
#define b2 4.0
#define c2 4.0
double f1(double x, void* param){
    return (x + b1)*x + c1;
}
double fnew(double x, void* param){
    if(x>2)return -1.0;
    return (x -4)*x + 3;
}
double dfnew(double x, void* param){
    if(x>2)return .0;
    return 2*x-4;
}

void fdfnew(double x, void* param, double* f, double* df)
{
    *f =  fnew(x,NULL);
    *df = dfnew(x,NULL);
}

double f2(double x, void* param){
    return (x + b2)*x + c2;
}

double df1(double x, void* param){
    return 2*x+b1;
}

double df2(double x, void* param){
    return 2*x+b2;
}

void fdf1(double x, void* param, double* f, double* df)
{
    *f =  (x + b1)*x + c1;
    *df = 2*x + b1;
}
void fdf2(double x, void* param, double* f, double* df)
{
    *f =  (x + b2)*x + c2;
    *df = 2*x + b2;
}


void bucket(int max_it, const gsl_root_fsolver_type * B){
    gsl_function F1,F2;
    F1.function=&f1;
    F2.function=&f2;

    printf("f.kwadratowa o jednokrotnych pierwiastkach\n");
    struct tms str, end;
    int i;
    int iters;
    double result;
    clock_t strt = times(&str);
    gsl_root_fsolver* slvr = gsl_root_fsolver_alloc(B);

    for(i=0;i<TIMES;i++){
        gsl_root_fsolver_set (slvr,&F1,-6,50);
        iters = do_bucket(slvr,max_it,&result);
    }
    clock_t endt = times(&end);
    pr_time(endt-strt,&str,&end);
    printf("\titer:\t%d\n",iters);

#ifndef BENCH
    printf("tutaj zawiedzie (podwojny pierwiastek)\n");
    strt = times(&str);
    for(i=0;i<TIMES;i++){
        gsl_root_fsolver_set (slvr,&F2,0,4);
        iters = do_bucket(slvr,max_it,&result);
    }
    endt = times(&end);
    pr_time(endt-strt,&str,&end);
    printf("\titer:\t%d\n",iters);
#endif
    gsl_root_fsolver_free (slvr);

}
void differential(int max_it,const gsl_root_fdfsolver_type * D)
{
    gsl_function_fdf F1,F2;
    F1.f = &f1;
    F1.df = &df1;
    F1.fdf = &fdf1;
    F2.f = &fnew;
    F2.df = &dfnew;
    F2.fdf = &fdfnew;

    gsl_root_fdfsolver* slvr = gsl_root_fdfsolver_alloc(D);
    struct tms str, end;
    int iters,i;
    double result;
    printf("f.kwadratowa o jednokrotnych pierwiastkach\n");
    clock_t strt = times(&str);
    for(i=0;i<TIMES;i++){
        gsl_root_fdfsolver_set (slvr,&F1,20);
        iters = do_differential(slvr,max_it,20,&result);
    }
    clock_t endt = times(&end);
    pr_time(endt-strt,&str,&end);
    printf("\titer:\t%d\n",iters);
#ifndef BENCH
    printf("tutaj zawiedzie (pochodna sie zeruje na przedziale)\n");
    strt = times(&str);
    for(i=0;i<TIMES;i++){
        gsl_root_fdfsolver_set (slvr,&F2,20);
        iters = do_differential(slvr,max_it,20,&result);
    }
    endt = times(&end);
    pr_time(endt-strt,&str,&end);
    printf("\titer:\t%d\n",iters);
#endif
    gsl_root_fdfsolver_free (slvr);
}

int main(int argc, char* argv[]){
int it =200;
#ifdef BENCH
        printf("<Regulafalsi>\n");
        bucket(it,gsl_root_fsolver_falsepos);
        printf("<Met. Brent'a>\n");
        bucket(it,gsl_root_fsolver_brent);
        printf("<Bisekcja>\n");
        bucket(it,gsl_root_fsolver_bisection);
        printf("<Met. Newtona>\n");
        differential(it,gsl_root_fdfsolver_newton);
        printf("<Met. Steffensona>\n");
        differential(it,gsl_root_fdfsolver_steffenson);
        printf("<Met. secant>\n");
        differential(it,gsl_root_fdfsolver_secant);
    return 0;
#endif
    if(argc!=2){
        printf("\nprovide opt-code\n");
        return 1;
    }
    int num = atoi(argv[1]);
    if(num>0 && num<7){
        printf("\n\tKazde podejscie ma limit %d prob\n\ti powtarzane jest %d razy (pomiar czsu)\n\n",it,TIMES);
        if(num==1){
        printf("<Regulafalsi>\n");
        bucket(it,gsl_root_fsolver_falsepos);
        }
        if(num==2){
        printf("<Met. Brent'a>\n");
        bucket(it,gsl_root_fsolver_brent);
        }
        if(num==3){
        printf("<Bisekcja>\n");
        bucket(it,gsl_root_fsolver_bisection);
        }
        if(num==4){
        printf("<Met. Newtona>\n");
        differential(it,gsl_root_fdfsolver_newton);
        }
        if(num==5){
        printf("<Met. Steffensona>\n");
        differential(it,gsl_root_fdfsolver_steffenson);
        }
        if(num==6){
        printf("<Met. secant>\n");
        differential(it,gsl_root_fdfsolver_secant);
        }
    }else{
        printf("\nwrong opt-code\n");
        return 2;
    }
    return 0;
}
