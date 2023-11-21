#include "../include/reparam.h"
#include <cmath>
#include <iostream>

using namespace std;

reparam::reparam(int N)
{
    x = new double[N];
    y = new double[N];

    x[0] = 0; x[N - 1] = 1;
    y[0] = 0; y[N - 1] = 1;
    for (int i = 1; i < N - 1; i++)
    {
        x[i] = 1.0 * i / (N - 1);
        y[i] = taur_over_taud(x[i]);
    }

    acc    = gsl_interp_accel_alloc();
    spline = gsl_spline_alloc(gsl_interp_cspline, N);
    poly = gsl_interp_alloc(gsl_interp_polynomial, N);

    gsl_spline_init(spline, y, x, N);
    gsl_interp_init(poly, y, x, N);
}

reparam::~reparam()
{
    gsl_spline_free(spline);
    gsl_interp_free(poly);
    gsl_interp_accel_free(acc);
    delete[] x;
    delete[] y;
}

double reparam::taur_over_taud(double x)
{
    return(x * log(1 / x) / (1 - x));
}

void reparam::map(double Amax, double decay_time, double raising_time,
        double& A, double& gamma, double& omega )
{
    double lambdap = exp(-1 / decay_time);
    double val;
    //poly->type->eval(poly->state,y,x,poly->size, (raising_time / decay_time), acc, &val);
//    double lambdam = exp(-1 / val / decay_time);
    double lambdam;

    lambdam = exp(-1 / gsl_spline_eval(spline, (raising_time / decay_time), acc) / decay_time);

    gamma = lambdap + lambdam;
    omega = -lambdap * lambdam;
    A     = Amax * (lambdap - lambdam) / (pow(lambdap, raising_time) - pow(lambdam, raising_time));
}


