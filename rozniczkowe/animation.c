#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

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
    }
	return alpha;
}
int func(double t, const double y[], double f[], void *params) {
	f[0] = y[1];
	f[1] = -G/L*sin(y[0]);
	return GSL_SUCCESS;
}
int jac(double t, const double y[], double *dfdy, double dfdt[], void *params) {
	gsl_matrix_view dfdy_mat  = gsl_matrix_view_array (dfdy, 2, 2);
	gsl_matrix * m = &dfdy_mat.matrix;
	gsl_matrix_set(m, 0, 0, 0.0);
	gsl_matrix_set(m, 0, 1, 1.0);
	gsl_matrix_set(m, 1, 0, -2.0*G/L*y[0]*y[1] - 1.0);
	gsl_matrix_set(m, 1, 1, -G/L*(y[0]*y[0] - 1.0));
	dfdt[0] = 0.0;
	dfdt[1] = 0.0;
	return GSL_SUCCESS;
}
double step=0.03;
double velo=3;
double duration=10.0;
int main(int argc, char *argv[]){
    if(argc<3){
        printf("usage: [step] [initial velocity] [duration]\n");
        return 1;
    }
    step=atof(argv[1]);
    velo=atof(argv[2]);
    duration=atof(argv[3]);
	FILE* out = fopen("animation.gpi","w");
    gsl_odeiv2_system sys = {func, jac, 2, NULL};
    const gsl_odeiv2_step_type * T = gsl_odeiv2_step_rk8pd;
    gsl_odeiv2_step * s = gsl_odeiv2_step_alloc (T, 2);
    gsl_odeiv2_control * c = gsl_odeiv2_control_y_new (1, 0.0);
    gsl_odeiv2_evolve * e = gsl_odeiv2_evolve_alloc (2);
    double y[2] = { 0.0, velo };
	double t;
	double alpha = euler(&t, step, velo);
	//fprintf(out,"clear\nreset\nset terminal gif animate delay 10\nset output \"animate.gif\"\nset isosample 40\n");
    fprintf(out,"clear\nreset\nset terminal gif animate delay 10\nset output \"animate.gif\"\n");
    fprintf(out,"set xrange [-%lf:%lf]\nset yrange [-%lf:%lf]\n",1.1*L,1.1*L,1.1*L,1.1*L);

    do{
        //fprintf(out,"plot \"<echo '%lf %lf'\" with points\n",t,alpha);
        fprintf(out,"plot \"<echo '0 0 %lf %lf'\" with vectors\n",-L*sin(alpha),-L*cos(alpha));
        //fprintf(out,"set arrow from 0, 0 to %lf, %lf as 1\nshow arrow 1\nplot\n",-L*sin(alpha),-L*cos(alpha));
        //printf("%lf\t%lf\n",t,alpha);
        alpha = euler(&t,-1,-1);
    }while(t<duration-step);
    fprintf(out,"plot \"<echo '%lf %lf'\" with points",t,alpha);
    fclose(out);

    return 0;
}
