#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <sys/time.h>
#include <gsl/gsl_statistics.h>

#define LAGRANGE 0
#define NEWTON 1

/**
* w takim prostokącie "zamknięte" będą wygenerowane punkty
*/
const int WIDTH = 800;
const int HEIGHT = 600;
/**
* wygenerowane punkty
*/
double *ya,*xa;
/**
* ilosc punktów
*/
int pointsNum;

/**
* alokuje tablice ze współrzędnymi punktami
* i generuje losowe (dość - zamknięte w prostokącie
* i dość równomiernie rozmieszczone na osi OX) punkty
*/
void gen(int num){
    srand(time(NULL));
    int x,y;
    //double spare = WIDTH;
    x=0;
    xa = (double*) malloc(num*sizeof(double));
    ya = (double*) malloc(num*sizeof(double));
    pointsNum=num;
    while(num--){
        y= rand()%(HEIGHT);
        x+= rand()%17+1;
        xa[pointsNum-num-1]=x;
        ya[pointsNum-num-1]=y;
    }
}

//////////////////////////////////////////////////////////////////////////////////
struct polynomialStruct {
    int n;
    int type;
    double *a;
};
typedef struct polynomialStruct polynomial;
//////////////////////////////////////////////////////////////////////////////////

/**
* alokuje wielomian, który będzie interpolował
* mode      - wybierz LAGRANGE albo NEWTON co oznacza metodę
* rootsNum  - stopień wielomianu (tj. ilość punktów)
*/
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
/**
* zwalnia zasoby wielomianu
*/
void polynomial_free(polynomial* pol){
    free(pol->a);
    free(pol);
}
/**
* do wielomianu a dodaje wielomian b
*/
int addPolynomial(polynomial* a,polynomial* b){
    if(a->n==b->n){
        int i=a->n;
        while(i--){
            a->a[i]+=b->a[i];
        }
        return 0;
    } else return 1;
}
/**
* towrzy kopię wielomianu pol
*/
polynomial* polycopy(polynomial* pol){
    polynomial* ret = polynomial_alloc(pol->type,pol->n);
    int i=pol->n;
    while(i--){
        ret->a[i]=pol->a[i];
    }
    return ret;
}
/**
* mnoży wielomian a razy jednomian (x + 'x') gdzie 'x' to argument
* zwraca 0 w razie sukcesu, wpp kod błędu
*/
int multiplyByMonomial(polynomial* a, double x){
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
/**
* mnoży wielomian razy skalar 'x'; zwraca 0 albo kod błędu
*/
int multiplyByScalar(polynomial* a, double x){
    int i= a->n;
    while(i--){
        a->a[i] *=x;
    }
    return 0;
}
/**
* oblicza wartoś wielomianu interpretacyjnego pol w punkcie x
* !wielomian musi być zainicjalizowany polynomial_init
*/
double polynomial_eval(polynomial* pol, double x){
    double tmp=0    ;
    int i=pol->n;
    while(i--){
        double ii=i;
        tmp+=(pol->a[i])*pow(x,ii);
    }
    return tmp;
}
/**
* pomocnicza funkcja ustalająca wielomian interpolacyjny
*/
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
                multiplyByMonomial(tmp,xa[j]);
            }
        }
        multiplyByScalar(tmp,under);
        addPolynomial(pol,tmp);
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
/**
* pomocnicza funkcja ustalająca wielomian interpolacyjny
*/
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

    j=1;
    while(j<roots){
        int k=1;
        while(k<=j){
            tab[j][k]=(tab[j][k-1]-tab[j-1][k-1])/(xa[j]-xa[j-k]);
            k++;
        }
        j++;
    }

    polynomial* tmp = polynomial_alloc(-1,pol->n);
    polynomial* n = polynomial_alloc(-1,pol->n);
    n->a[0]=1;
    i=0;
    while(i<roots){
        polynomial* m = polycopy(n);
        multiplyByScalar(n,tab[i][i]);
        addPolynomial(tmp,n);
        j=roots;
        n=m;
        multiplyByMonomial(n,-xa[i]);
        i++;
    }
    i=roots;
    while(i--){
        pol->a[i]=tmp->a[i];
    }
    polynomial_free(tmp);
    polynomial_free(n);
}
/**
* inicjalizuje wielomian; wybiera odpowiednią metodę
* !wilomian musi być zaalokowany polynomial_aloc
*/
void polynomial_init(polynomial* pol, double* xa, double* ya, int num){
    if(pol->type==LAGRANGE){
        Lagrange(pol,xa,ya,num);
    } else {
        Newton(pol,xa,ya,num);
    }
}

//////////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[]){
    //"parsowanie" argumentu
    int num = -1;
    if(argc<2){
        printf("\nNie podano liczby punktow\n");
        return 1;
    } else {
        num = atoi(argv[1]);
        if(num<1){
            printf("\nLiczba punktow musi byc niemniejza niz 1\n");
            return 2;
        }
        if(num>70){
            printf("\nLiczba punktow musi byc niewiększa niz 70\n");
            return 3;
        }
    }

    gen(num);
    struct timeval stop, start;

    printf("\n////////////////////////////////////////////////////////");
    //      TESTING

    int l = 10;
    double times[10];
    while(l--){
        gettimeofday(&start, NULL);
        gsl_interp_accel *aTmp = gsl_interp_accel_alloc();
        gsl_interp* iTmp = gsl_interp_alloc(gsl_interp_polynomial, pointsNum);
        gsl_interp_init (iTmp,xa,ya, pointsNum);
        gsl_interp_free(iTmp);
        gsl_interp_accel_free(aTmp);

        gettimeofday(&stop, NULL);
        times[l]=stop.tv_usec - start.tv_usec;
    }
    printf("\nGSL polynomial inicjalizacja zajęła średnio %lu usec\n",gsl_stats_mean(times,1,10));
    printf("\tOdchylenie standardowe %lu usec\n",gsl_stats_sd(times,1,10));

    l=10;
    while(l--){
        gettimeofday(&start, NULL);
        gsl_interp_accel *aTmp = gsl_interp_accel_alloc();
        gsl_interp* iTmp = gsl_interp_alloc(gsl_interp_cspline, pointsNum);
        gsl_interp_init (iTmp,xa,ya, pointsNum);
        gsl_interp_free(iTmp);
        gsl_interp_accel_free(aTmp);

        gettimeofday(&stop, NULL);
        times[l]=stop.tv_usec - start.tv_usec;
    }
    printf("\nGSL cspline inicjalizacja zajęła średnio %lu usec\n",gsl_stats_mean(times,1,10));
    printf("\tOdchylenie standardowe %lu usec\n",gsl_stats_sd(times,1,10));

    l=10;
    while(l--){
        gettimeofday(&start, NULL);
        gsl_interp_accel *aTmp = gsl_interp_accel_alloc();
        gsl_interp* iTmp = gsl_interp_alloc(gsl_interp_akima, pointsNum);
        gsl_interp_init (iTmp,xa,ya, pointsNum);
        gsl_interp_free(iTmp);
        gsl_interp_accel_free(aTmp);

        gettimeofday(&stop, NULL);
        times[l]=stop.tv_usec - start.tv_usec;
    }
    printf("\nGSL akima inicjalizacja zajęła średnio %lu usec\n",gsl_stats_mean(times,1,10));
    printf("\tOdchylenie standardowe %lu usec\n",gsl_stats_sd(times,1,10));

    l=10;
    while(l--){
        gettimeofday(&start, NULL);
        polynomial* lagrangeInterpT =polynomial_alloc(LAGRANGE,pointsNum);
        polynomial_init(lagrangeInterpT,xa,ya, pointsNum);
        polynomial_free(lagrangeInterpT);

        gettimeofday(&stop, NULL);
        times[l]=stop.tv_usec - start.tv_usec;
    }
    printf("\nMojej met. Lagrange'a inicjalizacja zajęła średnio %lu usec\n",gsl_stats_mean(times,1,10));
    printf("\tOdchylenie standardowe %lu usec\n",gsl_stats_sd(times,1,10));

    l=10;
    while(l--){
        gettimeofday(&start, NULL);
        polynomial* newtonInterpT =polynomial_alloc(NEWTON,pointsNum);
        polynomial_init(newtonInterpT,xa,ya, pointsNum);
        polynomial_free(newtonInterpT);

        gettimeofday(&stop, NULL);
        times[l]=stop.tv_usec - start.tv_usec;
    }
    printf("\nMojej met. Newton'a inicjalizacja zajęła średnio %lu usec\n",gsl_stats_mean(times,1,10));
    printf("\tOdchylenie standardowe %lu usec\n",gsl_stats_sd(times,1,10));

    printf("////////////////////////////////////////////////////////\n\n");

gettimeofday(&start, NULL);
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_interp* interp = gsl_interp_alloc(gsl_interp_polynomial, pointsNum);
    gsl_interp_init (interp,xa,ya, pointsNum);
    gsl_interp* cspline = gsl_interp_alloc(gsl_interp_cspline , pointsNum);
    gsl_interp_init (cspline,xa,ya, pointsNum);
    gsl_interp* akima = gsl_interp_alloc(gsl_interp_akima , pointsNum);
    gsl_interp_init (akima,xa,ya, pointsNum);
    polynomial* lagrangeInterp =polynomial_alloc(LAGRANGE,pointsNum);
    polynomial* newtonInterp =polynomial_alloc(NEWTON,pointsNum);
    polynomial_init(lagrangeInterp,xa,ya, pointsNum);
    polynomial_init(newtonInterp,xa,ya, pointsNum);
gettimeofday(&stop, NULL);
printf("\n<< %lu >>\n",stop.tv_usec - start.tv_usec);

    int res = pointsNum/3;
    int lenBak= (xa[pointsNum-1] - xa[0]+1)*res;
    int len =lenBak;
    //printf("\n\n\n<<%i>>\n\n\n",len);
    double* x = malloc(len*sizeof(double));
    double* gsl = malloc(len*sizeof(double));
    double* lgrng = malloc(len*sizeof(double));
    double* nwtn = malloc(len*sizeof(double));
    double* akma = malloc(len*sizeof(double));
    double* cspn = malloc(len*sizeof(double));

    double step = 1.0/res;

    // oblicznie położeń punktów siatki
    double xi= xa[pointsNum-1];
    while(len--){
        if(xi<xa[0])xi=xa[0];// gsl się buntuje
        x[len]=xi;
        xi-=step;
    }



    // polynomial interp     - mierzymy czas
    len =lenBak;
    gettimeofday(&start, NULL);
    while(len--){
        gsl[len]=gsl_interp_eval(interp, xa, ya, x[len], acc);
    }
    gettimeofday(&stop, NULL);
    printf("GSL polynomial ewaluacja zajęła %lu usec\n", stop.tv_usec - start.tv_usec);

    // akima - mierzymy czas
    len =lenBak;
    gettimeofday(&start, NULL);
    while(len--){
        akma[len]=gsl_interp_eval(akima, xa, ya, x[len], acc);
    }
    gettimeofday(&stop, NULL);
    printf("GSL akima ewaluacja zajęła %lu usec\n", stop.tv_usec - start.tv_usec);

    // cspline - mierzymy czas
    len =lenBak;
    gettimeofday(&start, NULL);
    while(len--){
        cspn[len]=gsl_interp_eval(cspline, xa, ya, x[len], acc);
    }
    gettimeofday(&stop, NULL);
    printf("GSL cspline ewaluacja zajęła %lu usec\n", stop.tv_usec - start.tv_usec);


    // lagrange - mierzymy czas
    len =lenBak;
    gettimeofday(&start, NULL);
    while(len--){
        lgrng[len]=polynomial_eval(lagrangeInterp,x[len]);
    }
    gettimeofday(&stop, NULL);
    printf("Mojej met. Lagrange'a ewaluacja zajęła %lu usec\n", stop.tv_usec - start.tv_usec);

    // newton -mierzymy czas
    len =lenBak;
    gettimeofday(&start, NULL);
    while(len--){
        nwtn[len]=polynomial_eval(newtonInterp,x[len]);
    }
    gettimeofday(&stop, NULL);
    printf("Mojej met. Newton'a ewaluacja zajęła %lu usec\n", stop.tv_usec - start.tv_usec);

    FILE* out = fopen("out.tmp","w");
    len = (xa[pointsNum-1] - xa[0]+1)*res;
    while(len--){
        fprintf(out,"%f\t%f\t%f\t%f\t%f\t%f\n",x[len],gsl[len],lgrng[len],nwtn[len],cspn[len],akma[len]);
    }
    polynomial_free(lagrangeInterp);
    polynomial_free(newtonInterp);
    gsl_interp_free(interp);
    gsl_interp_free(cspline);
    gsl_interp_free(akima);
    gsl_interp_accel_free(acc);
    fclose(out);

    FILE* out2 = fopen("outPts.tmp","w");
    int i =pointsNum;
    while(i--){
        fprintf(out2,"%f\t%f\n",xa[i],ya[i]);
    }
    fclose(out2);
    free(xa);free(ya);free(x);free(lgrng);free(nwtn);free(gsl);
    return 0;
}
