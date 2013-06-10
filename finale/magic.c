#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_linalg.h>

#define DIMM 20
//#define IMPLICIT

#define H 1.0
#define DT 1.0
#define VAL 10.0
#define ERR .000001
#define STEP 30
#define TIMES 10

#define ZMIN -1.
#define ZMAX 1.

double T[DIMM][DIMM];
double Tmatrix[DIMM*DIMM*DIMM*DIMM];
double Tlast[DIMM*DIMM];
double Tfound[DIMM*DIMM];

double f(double x, double y){
    x-=.5;
    y-=.5;
    return y*y*8*y;
}

gsl_matrix_view m;
gsl_vector_view b;
gsl_vector *x;
gsl_permutation * p;

int updateVal_implicit(int i, int j, int dimm){
    if(i==0 || j==0 || i==dimm-1 || j==dimm-1){
        T[i][j]=f(((double)i)/(dimm-1),((double)j)/(dimm-1));
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
            //if(i==0 || j==0 || i==dimm-1 || j==dimm-1) printf("  %.0f\t",T[i][j]);
            //else
            printf("%.3f\t",T[i][j]);
        }
        printf("\n");
    }
}
void outprint(char* name, int dimm){
    int i,j;
    FILE* out = fopen(name,"w");
    for(i=0;i<DIMM;i++){
        for(j=0;j<DIMM;j++){
            #ifdef IMPLICIT
            fprintf(out,"%d\t%d\t%f\n",i,j,T[i][j]);
            #else
            fprintf(out,"%d\t%d\t%f\n",i,j,Tfound[i*dimm+j]);
            #endif
        }
        fprintf(out,"\n");
    }
    fclose(out);
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
                Tfound[i*dimm+j]=f(((double)i)/(dimm-1),((double)j)/(dimm-1));
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
    int i=0;
    char buff[8];
    FILE* f = fopen("plot.gpi","w");
    fprintf(f,"set terminal gif animate delay 10\nset output \"animate.gif\"\n");
    fprintf(f,"set xrange[0:%d]\n",DIMM-1);
    fprintf(f,"set yrange[0:%d]\n",DIMM-1);
    fprintf(f,"set zrange[%f:%f]\n",ZMIN,ZMAX);
    //fprintf(f,"set pm3d\n",VAL);
    //fprintf(f,"set palette model XYZ rgbformulae 7,5,15\n",VAL);

    #ifdef IMPLICIT
    init_implicit(DIMM);
    while(!iteration_implicit(DIMM)){
    #else
    init_explicit(DIMM);
    while(!iteration_explicit(DIMM)){
    #endif
        //printf(".");
        sprintf(buff,"tmp/%d",i++);
        outprint(buff,DIMM);
        fprintf(f,"splot \"tmp/%d\" with lines\n",i);
    }
    fclose(f);
    system("gnuplot plot.gpi");
    #ifdef IMPLICIT
    print_implicit(DIMM);
    #else
    print_explicit_solution(DIMM);
    #endif
    return 0;
}
