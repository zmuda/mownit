#include <stdio.h>
#include <stdlib.h>

#define DIMM 70
#define VAL 10
#define ERR .0001

double V[DIMM][DIMM];


int updateVal(int i, int j){
    if(i==0 || j==0 || i==DIMM-1 || j==DIMM-1){
        V[i][j]=0;
        return 1;
    }
    if(i==DIMM/2 || i==(int)((DIMM-.5)/2)){
        if(j==DIMM/2 || j==(int)((DIMM-.5)/2)){
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
void print(){
    int i,j;
    for(i=0;i<DIMM;i++){
        for(j=0;j<DIMM;j++){
            if(i==0 || j==0 || i==DIMM-1 || j==DIMM-1) printf("  %.0f\t",V[i][j]);
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
void fix(){
    V[DIMM/2][DIMM/2]=VAL;
    V[DIMM/2][(int)((DIMM-.5)/2)]=VAL;
    V[(int)((DIMM-.5)/2)][(int)((DIMM-.5)/2)]=VAL;
    V[(int)((DIMM-.5)/2)][DIMM/2]=VAL;
}
void init(){
    int i,j;
    for(i=0;i<DIMM;i++){
        for(j=0;j<DIMM;j++){
            V[i][j]=0;
        }
    }
    fix();
}
int iteration(){
    int i,j;
    int res=1;
    for(i=0;i<DIMM;i++){
        for(j=0;j<DIMM;j++){
            res*=updateVal(i,j);
        }
    }
    return res;
}

int main(int argc, char** argv){
    init();
    int i=0;
    char buff[8];
    FILE* f = fopen("plot.gpi","w");
    fprintf(f,"set terminal gif animate delay 10\nset output \"animate.gif\"\n");
    fprintf(f,"set xrange[0:%d]\n",DIMM-1);
    fprintf(f,"set yrange[0:%d]\n",DIMM-1);
    fprintf(f,"set zrange[0:%d]\n",VAL);
    fprintf(f,"set pm3d\n",VAL);
    fprintf(f,"set palette model XYZ rgbformulae 7,5,15\n",VAL);

    while(!iteration()){
        //printf(".");
        sprintf(buff,"tmp/%d",i++);
        outprint(buff);
        fprintf(f,"splot \"tmp/%d\" with lines\n",i);
    }
    fclose(f);
    system("gnuplot plot.gpi");

    print();
    return 0;
}
