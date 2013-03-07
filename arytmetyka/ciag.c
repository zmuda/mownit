#include <stdlib.h>
#include <stdio.h>

/*
 x{n+1}= x{n} + 3.0 * x{n} (1 - x{n})
 w tym ciągu błąd reprezentacji jest wzmacniany z każdą iteracją (rekurencja)
 dodatkowo jest źle uwarunkowany dla x->0 albo x->4/3
 no a dla x->1 mamy blad anulowania
*/

void work(int steps){
    float xf = 0.01f;
    float tmpf;
    double xd = 0.01;
    double tmpd;
    //printf("\n\t:double:\t:float:\n\n");
    FILE* out = fopen("out.tmp","w");
    fprintf(out,"%f\t%i\t%lf\n",xf,0,xd);
    int i=0;
    while(++i<steps){
        tmpf=1-xf;
        tmpf*=xf;
        tmpf*=3.0f;
        xf+=tmpf;

        tmpd=1-xd;
        tmpd*=xd;
        tmpd*=3.0;
        xd+=tmpd;

        fprintf(out,"%f\t%i\t%lf\n",xf,i,xd);
        //printf("\t%f\t%lf\n",xf,xd);
    }
    fclose(out);
}


int main(int argc, char *argv[]){
    int steps = -1;
    if(argc<2){
        printf("\nNie podano liczby poczatkowych wyrazow\n");
        return 1;
    } else {
        steps = atoi(argv[1]);
        if(steps<1){
            printf("\nLiczba wyrazow musi byc niemniejza niz 1\n");
            return 2;
        }
    }
    work(steps);
    return 0;
}

