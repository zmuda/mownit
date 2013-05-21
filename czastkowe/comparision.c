#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_statistics.h>

#define DIMM 200
#define VAL 10
#define ERR .0001
#define STEP 30
#define TIMES 10

double V[DIMM][DIMM];


int updateVal(int i, int j, int dimm){
    if(i==0 || j==0 || i==dimm-1 || j==dimm-1){
        V[i][j]=0;
        return 1;
    }
    if(i==dimm/2 || i==(int)((dimm-.5)/2)){
        if(j==dimm/2 || j==(int)((dimm-.5)/2)){
            V[i][j]=VAL;
            return 1;
        }
    }
    double tmp=V[i][j];
    V[i][j]=(V[i][j+1]+V[i][j-1]+V[i+1][j]+V[i-1][j])/4;
    if(tmp-V[i][j]<0){
        return (tmp-V[i][j]>-ERR);
    } else {
        return (tmp-V[i][j]<ERR);
    }
}
void print(int dimm){
    int i,j;
    for(i=0;i<dimm;i++){
        for(j=0;j<dimm;j++){
            if(i==0 || j==0 || i==dimm-1 || j==dimm-1) printf("  %.0f\t",V[i][j]);
            else printf("%.3f\t",V[i][j]);
        }
        printf("\n");
    }
}
void outprint(char* name){
    int i,j;
    FILE* out = fopen(name,"w");
    for(i=0;i<DIMM;i++){
        for(j=0;j<DIMM;j++){
            fprintf(out,"%d\t%d\t%f\n",i,j,V[i][j]);
        }
        fprintf(out,"\n");
    }
    fclose(out);
}
void fix(int dimm){
    V[dimm/2][dimm/2]=VAL;
    V[dimm/2][(int)((dimm-.5)/2)]=VAL;
    V[(int)((dimm-.5)/2)][(int)((dimm-.5)/2)]=VAL;
    V[(int)((dimm-.5)/2)][dimm/2]=VAL;
}
void init(int dimm){
    int i,j;
    for(i=0;i<dimm;i++){
        for(j=0;j<dimm;j++){
            V[i][j]=0;
        }
    }
    fix(dimm);
}
int iteration(int dimm){
    int i,j;
    int res=1;
    for(i=0;i<dimm;i++){
        for(j=0;j<dimm;j++){
            res*=updateVal(i,j,dimm);
        }
    }
    return res;
}

int main(int argc, char** argv){
    double times[TIMES];
    struct timeval stop, start;
    int l;

    int dimm=DIMM;
    while(dimm>0){
        printf("\t\trozmiar:\t%d\n",dimm);
        l=TIMES;
        while(l--) {
            gettimeofday(&start, NULL);
            init(dimm);
            while(!iteration(DIMM));
            gettimeofday(&stop, NULL);
            times[l]=(double)stop.tv_usec+(double)stop.tv_sec*1000000 - (double)start.tv_usec-(double)start.tv_sec*1000000;
            printf("Zajelo %f\n",times[l]);
        }
        double yi,ei;
        yi=gsl_stats_mean(times,1,TIMES);
        ei=gsl_stats_sd(times,1,TIMES);
        printf("Proste mnozenie srednio zajelo %f usec\n",yi);
        printf("\tOdchylenie standardowe %f usec\n",ei);
        dimm-=STEP;
    }
    return 0;
}
