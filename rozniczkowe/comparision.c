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
double step=0.1;
double velo=0.1;
double duration=100.0;
int main(int argc, char *argv[]){
    if(argc<3){
        printf("usage: [step] [initial velocity] [duration]\n");
        return 1;
    }
    step=atof(argv[1]);
    velo=atof(argv[2]);
    duration=atof(argv[3]);
	FILE* out = fopen("outEuler.tmp","w");
    gsl_odeiv2_system sys = {func, jac, 2, NULL};
    const gsl_odeiv2_step_type * T = gsl_odeiv2_step_rk8pd;
    gsl_odeiv2_step * s = gsl_odeiv2_step_alloc (T, 2);
    gsl_odeiv2_control * c = gsl_odeiv2_control_y_new (1, 0.0);
    gsl_odeiv2_evolve * e = gsl_odeiv2_evolve_alloc (2);
    double y[2] = { 0.0, velo };
	double t;
	double alpha = euler(&t, step, velo);
    do{
        fprintf(out,"%lf\t%lf\n",t,alpha);
        //printf("%lf\t%lf\n",t,alpha);
        alpha = euler(&t,-1,-1);
    }while(t<duration);

    fclose(out);
    out = fopen("outGSL.tmp","w");
    t=0;
    while(t<duration-step){
        //gsl_odeiv2_evolve_apply (e, c, s,&sys,&t, 10,&step, y);
        gsl_odeiv2_evolve_apply_fixed_step (e, c, s,&sys,&t,step, y);
        fprintf(out,"%lf\t%lf\n",t,y[0]);
        //printf("%lf\t%lf\n",t,y[0]);
    }
    fclose(out);

    return 0;
}
