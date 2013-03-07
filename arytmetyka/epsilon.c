#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_ieee_utils.h>

int main(){
    float a=1;
    while(1+a>1){
        a/=2;
    }
    printf("\nMaszynewe epsilon: %.40f\n",a);
    printf("(zakładając, że pomieściło się na 40 mscach dziesiętnych dla float)\n",a);
    printf("\nMaszynewe epsilon w reprezentacji cecha-mantysa\n");
    gsl_ieee_printf_float(&a);
    printf("\n");
    return 0;
}

