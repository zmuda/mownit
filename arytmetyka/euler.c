#include <stdlib.h>
#include <stdio.h>

/*
y' +7y +1 = 0
y(0)=1
na przedziale (0,1)
*/

void euler(FILE* out, unsigned steps, double y0){
	double stepLen = 1.0/steps;
	double x = 0;
	double y = y0;
	while(x<1){
		fprintf(out,"%lf\t%lf\n",x,y);
		y+= stepLen*(-7y-1);		
		x+=stepLen;
	}
}


int main(int argc, char *argv[]){
	double x0=1.0;
	unsigned steps=-1;
	if(argc<2){
        	printf("\nPodaj ilosc krokow\n");
       		return 1;
   	}
	steps = atoi(argv[1]);
	FILE* out = fopen("outEuler.tmp","w");
	euler(out,steps,1.0);
}
