#include <stdlib.h>
#include <stdio.h>

void work(int steps){
    float xf = 0.01f;
    float tmpf;
    double xd = 0.01;
    double tmpd;
    while(steps--){

    }
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

    return 0;
}
