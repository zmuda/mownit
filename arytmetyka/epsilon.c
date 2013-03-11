#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_ieee_utils.h>

int main(){
    float a=1;
    float tmp;
    float lasta;
    tmp=1+a;
    while(tmp>1){
        lasta=a;
        a/=2;
        tmp=1+a;
    }
    printf("\nMaszynewe epsilon: %.40f\n",lasta);
    printf("(zakładając, że pomieściło się na 40 mscach dziesiętnych dla float)\n",lasta);
    printf("\nMaszynewe epsilon float w reprezentacji cecha-mantysa\n");
    gsl_ieee_printf_float(&lasta);
    printf("\n");
    return 0;
}

