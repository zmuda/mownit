#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_linalg.h>

#define DIMM 23

#define H 1.0
#define DT 1.0
#define VAL 10.0
#define ERR .0001
#define STEP 30
#define TIMES 10

double T[DIMM][DIMM];
double Tmatrix[DIMM*DIMM*DIMM*DIMM];
double Tlast[DIMM*DIMM];
double Tfound[DIMM*DIMM];


gsl_matrix_view m;
gsl_vector_view b;
gsl_vector *x;
gsl_permutation * p;

int updateVal_implicit(int i, int j, int dimm){
    if(i==0 || j==0 || i==dimm-1 || j==dimm-1){
        T[i][j]=VAL;
        return 1;
    }
    double tmp=T[i][j];
    T[i][j]=( H*H*tmp + DT*(T[i][j+1]+T[i][j-1]+T[i+1][j]+T[i-1][j]) ) / (H*H+4*DT);
    if(tmp-T[i][j]<0){
        return (tmp-T[i][j]>-ERR);
    } else {
        return (tmp-T[i][j]<ERR);
    }
}
void print_implicit(int dimm){
    int i,j;
    for(i=0;i<dimm;i++){
        for(j=0;j<dimm;j++){
            if(i==0 || j==0 || i==dimm-1 || j==dimm-1) printf("  %.0f\t",T[i][j]);
            else printf("%.3f\t",T[i][j]);
        }
        printf("\n");
    }
}
void print_explicit(int dimm){
    int i,j;
    for(i=0;i<dimm*dimm;i++){
        for(j=0;j<dimm*dimm;j++){
            printf("%.1f ",Tmatrix[i*dimm*dimm+j]);
        }
        printf("\t%.3f\n",Tfound[i*dimm+j]);
    }printf("\n");
}
void print_explicit_solution(int dimm){
    int i,j;
    for(i=0;i<dimm;i++){
        for(j=0;j<dimm;j++){
            printf("%.2f\t",Tfound[i*dimm+j]);
        }
        printf("\n");
    }printf("\n");
}

void init_implicit(int dimm){
    int i,j;
    for(i=0;i<dimm;i++){
        for(j=0;j<dimm;j++){
            T[i][j]=0;
        }
    }
}
int k(int x, int y, int dimm){
    return x*dimm+y;
}
void init_explicit(int dimm){
    int i,j;
    for(i=0;i<dimm;i++){
        for(j=0;j<dimm;j++){
            if(i==0 || j==0 || i==dimm-1 || j==dimm-1){
                Tmatrix[(i*dimm+j)*(dimm*dimm)+(i*dimm+j)]=1;
                Tfound[i*dimm+j]=VAL;
            } else {
                Tfound[i*dimm+j]=0;
                Tmatrix[(i*dimm+j)*(dimm*dimm)+(i*dimm+j)]=(H*H+4*DT)/H/H;
                Tmatrix[(i*dimm+j)*(dimm*dimm)+(i*dimm+j-dimm)]=-DT/H/H;
                Tmatrix[(i*dimm+j)*(dimm*dimm)+(i*dimm+j-1)]=-DT/H/H;
                Tmatrix[(i*dimm+j)*(dimm*dimm)+(i*dimm+j+1)]=-DT/H/H;
                Tmatrix[(i*dimm+j)*(dimm*dimm)+(i*dimm+j+dimm)]=-DT/H/H;
            }
        }
    }
    m = gsl_matrix_view_array (Tmatrix, dimm*dimm, dimm*dimm);
    x = gsl_vector_alloc (dimm*dimm);
    p = gsl_permutation_alloc (dimm*dimm);
    int s;
    gsl_linalg_LU_decomp (&m.matrix, p, &s);
    //print_explicit(dimm);
    //print_explicit_solution(dimm);
}
int iteration_implicit(int dimm){
    int i,j;
    int res=1;
    for(i=0;i<dimm;i++){
        for(j=0;j<dimm;j++){
            res*=updateVal_implicit(i,j,dimm);
        }
    }
    return res;
}
int iteration_explicit(int dimm){
    int i,j;
    i=dimm*dimm;
    while(i--){
        Tlast[i]=Tfound[i];
    }
    b = gsl_vector_view_array (Tlast, dimm*dimm);
    gsl_linalg_LU_solve (&m.matrix, p, &b.vector, x);
    i=dimm*dimm;
    int ret =1;
    while(i--){
        Tfound[i]=x->data[i];
        if(Tfound[i]-Tlast[i]>ERR)ret=0;
    }
    return ret;
}
int main3(int argc, char** argv){
    init_implicit(DIMM);
    while(!iteration_implicit(DIMM)){
        print_implicit(DIMM);
    }
    return 0;
}
int main4(int argc, char** argv){
    init_explicit(DIMM);
    while(!iteration_explicit(DIMM)){
        print_explicit_solution(DIMM);
        sleep(1);
    }
    gsl_permutation_free (p);
    gsl_vector_free (x);
    return 0;
}

int main(int argc, char** argv){
    double times[TIMES];
    struct timeval stop, start;
    int l;

    int dimm=DIMM;
    //while(dimm>0){
        printf("\n\tIMPLICIT\trozmiar:\t%d\n",dimm);
        l=TIMES;
        while(l--) {
            gettimeofday(&start, NULL);
            init_implicit(dimm);
            while(!iteration_implicit(DIMM));
            gettimeofday(&stop, NULL);
            times[l]=(double)stop.tv_usec+(double)stop.tv_sec*1000000 - (double)start.tv_usec-(double)start.tv_sec*1000000;
            //printf("Zajelo %f\n",times[l]);
        }
        double yi,ei;
        yi=gsl_stats_mean(times,1,TIMES);
        ei=gsl_stats_sd(times,1,TIMES);
        printf("Proste mnozenie srednio zajelo %f usec\n",yi);
        printf("\tOdchylenie standardowe %f usec\n",ei);
     //   dimm-=STEP;
    //}
    dimm=DIMM;
    //while(dimm>0){
        printf("\n\tEXPLICIT\trozmiar:\t%d\n",dimm);
        l=TIMES;
        while(l--) {
            gettimeofday(&start, NULL);
            init_explicit(dimm);
            while(!iteration_explicit(DIMM));
            gettimeofday(&stop, NULL);
            times[l]=(double)stop.tv_usec+(double)stop.tv_sec*1000000 - (double)start.tv_usec-(double)start.tv_sec*1000000;
            //printf("Zajelo %f\n",times[l]);
        }
        yi=gsl_stats_mean(times,1,TIMES);
        ei=gsl_stats_sd(times,1,TIMES);
        printf("Proste mnozenie srednio zajelo %f usec\n",yi);
        printf("\tOdchylenie standardowe %f usec\n",ei);
     //   dimm-=STEP;
    //}

    return 0;
}
