#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_ieee_utils.h>

/*
nie ma takiej sily, ktora pozwoli domyslic sie, ze kiedys to byla 1/3
(przed przyporzadkowaniem do float'a)
*/
int main(){
    printf("Jedna trzecia float\n");
    float a=1.0/3;
    gsl_ieee_printf_float(&a);
    printf("\nJedna trzecia zrzutowana na double\n");
    double b=a;
    gsl_ieee_printf_double  (&b);
    printf("\n");

    return 0;
}
