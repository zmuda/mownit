#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>


const int WIDTH = 800;
const int HEIGHT = 600;
double *ya,*xa;
int stepsBak;

void gen(int steps){
    srand(time(NULL));
    int x,y;
    double spare = WIDTH;
    x=0;
    xa = (double*) malloc(steps*sizeof(double));
    ya = (double*) malloc(steps*sizeof(double));
    stepsBak=steps;
    while(steps--){
        spare = WIDTH-x;
        y= rand()%(HEIGHT);
        int tmp = (spare/steps)+1;
        if(tmp<2)tmp=spare;
        x+= rand()%tmp+1;
        xa[stepsBak-steps-1]=x;
        ya[stepsBak-steps-1]=y;
        xa[stepsBak-steps-1]=x;
        ya[stepsBak-steps-1]=x*x*x+x*x+x+1;

        //printf("%f\t%f\n",xa[stepsBak-steps-1],ya[stepsBak-steps-1]);
    }
}
void gsl_polynomial(){
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_interp* inter = gsl_interp_alloc(gsl_interp_polynomial, stepsBak);
    gsl_interp_init (inter,xa,ya, stepsBak);

    double xi, yi;
    for (xi = xa[0]; xi < xa[stepsBak-1]; xi += 1.0 ) {
        yi = gsl_interp_eval(inter, xa, ya, xi, acc);
        printf("%lf %lf\n", xi, yi);
    }
    gsl_interp_free(inter);
    gsl_interp_accel_free(acc);
    free(xa);
    free(ya);
}
//////////////////////////////////////////////////////////////////////////////////
struct polynomialStruct {
    int n;
    int type;
    double *a;
};
typedef struct polynomialStruct polynomial;

polynomial* polynomial_alloc(int mode, int rootsNum){
    polynomial* ret = (polynomial*)malloc( sizeof(polynomial) );
    ret->a=malloc(rootsNum*sizeof(double));
    ret->n=rootsNum;
    ret->type=mode;
    int i=rootsNum;
    while(i--){
        (ret->a)[i]=0.0;
    }
    return ret;
    return 0;
}
void polynomial_free(polynomial* pol){
    free(pol->a);
    free(pol);
}
int aisaplusb(polynomial* a,polynomial* b){
    if(a->n==b->n){
        int i=a->n;
        while(i--){
            a->a[i]+=b->a[i];
        }
        return 0;
    } else return 1;
}
polynomial* polycopy(polynomial* pol){
    polynomial* ret = polynomial_alloc(pol->type,pol->n);
    int i=pol->n;
    while(i--){
        ret->a[i]=pol->a[i];
    }
    return ret;
}
int aisadotmonomial(polynomial* a, double x){
    int i= a->n;
    if(a->a[i-1]!=0.0){
        return 1;
    }
    polynomial* copy = polycopy(a);
    while(i-->1){
        a->a[i] = x* copy->a[i];
        a->a[i] += copy->a[i-1];
    }
    a->a[0]=x*copy->a[0];
}
int aisadotscalar(polynomial* a, double x){
    int i= a->n;
    while(i--){
        a->a[i] *=x;
    }
    return 0;
}

double polynomial_eval(polynomial* pol, double x){
    double tmp=0    ;
    int i=pol->n;
    while(i--){
        double ii=i;
        tmp+=(pol->a[i])*pow(x,ii);
    }
    return tmp;
}

void Lagrange(polynomial* pol, double* xa, double* ya, int num){
    int roots = pol->n;
    int i = roots;
    int j;
    double under;
    polynomial* tmp;
    while(i--){
        tmp=polynomial_alloc(-1,pol->n);
        tmp->a[0]=1;
        under=ya[i];
        j=roots;
        while(j--){
            if(j!=i){
                under/=( xa[i]-xa[j] );
                aisadotmonomial(tmp,xa[j]);
            }
        }
        aisadotscalar(tmp,under);
        aisaplusb(pol,tmp);
    }
    i= roots;
    while(i--){
        if( i%2 ){
            pol->a[i]=-pol->a[i];
        }
        if( !(roots%2) ){
            pol->a[i]=-pol->a[i];
        }
    }
}
void Newton(polynomial* pol, double* xa, double* ya, int num){
    int roots = pol->n;
    int i = roots;
    int j = roots;
    double** tab = malloc(roots * sizeof(double*) );
    while(j--){
        tab[j]=malloc( (j+1)* sizeof(double) );
    }
    j=roots;
    while(j--){
        tab[j][0]=ya[j];
    }
    //tab[1][1]=(tab[0][1]+tab[0][0])/(ya[1]-ya[0]);

    j=1;
    while(j<roots){
        int k=1;
        while(k<=j){
            tab[j][k]=(tab[j][k-1]-tab[j-1][k-1])/(ya[j]-ya[j-k]);
            //printf("(%i,%i)=((%i,%i)-(%i,%i))/(%i-%i)\n",j,k,j,k-1,j-1,k-1,j,j-k);
            printf("(%i,%i)=%f\t",j,k,tab[j][k]);
            k++;
        }
        j++;
    }
    /*
    j=roots;
    while(j--){
        pol->a[j]=tab[j][j];
        printf("<%f>\n",tab[j][j]);
    }
    */
}

void polynomial_init(polynomial* pol, double* xa, double* ya, int num){
    //if(pol->type==0){
    //Lagrange(pol,xa,ya,num);
    //} else {
        Newton(pol,xa,ya,num);
    //}
}



//////////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[]){
    int steps = -1;
    if(argc<2){
        printf("\nNie podano liczby punktow\n");
        return 1;
    } else {
        steps = atoi(argv[1]);
        if(steps<1){
            printf("\nLiczba punktow musi byc niemniejza niz 1\n");
            return 2;
        }
    }
    //printf("FUUUUUUUUUUUUUUUUU");
    gen(steps);
    //gsl_polynomial();
    polynomial* tmp=polynomial_alloc(-1,stepsBak);
    polynomial_init(tmp,xa,ya,stepsBak);

    //tmp->a[0]=1;
    //aisadotscalar(tmp,5);

    int i= stepsBak;
    while(i--){
        printf("%f\n",tmp->a[i]);
    }

    i= stepsBak;
    while(i--){
        printf("%f\t%f\t%f\n",xa[i],ya[i],polynomial_eval(tmp,xa[i]));
    }

    polynomial_free(tmp);

    return 0;
}
