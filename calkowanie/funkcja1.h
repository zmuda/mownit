//zwykle
double f(double x,void* null){
    return 3*(x*x-2*x);
}
//oscylacja
double g(double x,void* null){
    return  x * sin(30 * x) * cos(x);
}
//osobliwosc
double h(double x,void* null){
    //return  1/(x-1);
    return  log(x) / sqrt(x);
}
//adaptacyjne
double k(double x,void* null){
    return  sin(x*2)/x;
}

//const double ERROR = 1e-5;
const int POINTS =  30;





