#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_ieee_utils.h>

/*
mantysa nie jest znormalizowana, kiedy nie nalezy do przedzialu [1,2)
*/
int main(){
    printf("Kilka coraz mniejszych float'ów\n");
    float a=13;
    //while(1+a>1){
    while(a>0){
        a/=2;
        gsl_ieee_printf_float(&a);
        printf("\n");
    }
    return 0;
}
