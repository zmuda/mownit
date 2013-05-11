#include <stdlib.h>
#include <stdio.h>
#include <math.h>
/*
y' +7y +1 = 0
y(0)=1
na przedziale (0,1)
*/
double fun(double x){
    return sin(x);
}
#define G 10.0
#define L 1.0
double alpha;
double v;
double t;
double dt;
double euler(double* tt, double stepLen, double v0){
    if(stepLen>0){
        alpha=0;
        v=v0;
        t==0;
        dt=stepLen;
    }
    {
		t+=dt;
		*tt=t;
		v+=-dt*G/L*sin(alpha);
        alpha+=dt*v;
        //fprintf(out,"%lf\t%lf\n",t,alpha);
	}
	return alpha;
}


int main(int argc, char *argv[]){
	FILE* out = fopen("outEuler.tmp","w");
	//euler(out,1, 1.0);
	double t;
	double alpha = euler(&t, .1, 0.1);
    do{
        //sleep(1);
        fprintf(out,"%lf\t%lf\n",t,alpha);
        //printf("%lf\t%lf\n",t,alpha);
        alpha = euler(&t,-1,-1);
    }while(t<10);
    return 0;
}
