#include<algorithm>
#include<gsl/gsl_math.h>
#include<armadillo>
#include"../include/utils.h"
#include<gsl/gsl_sf.h>
using namespace std;

void utils::w_from_logW(const double *src, double* out, unsigned int n){
    double maxLogWeight = *max_element(src,src+n);
    for(unsigned int i=0;i<n;i++) out[i] = exp(src[i]-maxLogWeight);
}

double utils::Z_from_logW(const double *src, unsigned int n){
    double maxLogWeight = *max_element(src,src+n);
    double Z=0;
    for(unsigned int i=0; i<n; i++){
        Z+=exp(src[i]-maxLogWeight);
    }
    return(Z*exp(maxLogWeight));
}

double utils:: Z_factor(double x,double xprop,double sigma,double cutoff){
    return(log(gsl_sf_erf_Q((cutoff-x)/sigma)/gsl_sf_erf_Q((cutoff-xprop)/sigma)));
}


