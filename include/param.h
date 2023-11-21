#ifndef PARAM_H
#define PARAM_H

#include<fstream>
#include "constants.h"
#include "../include/reparam.h"

class param {
public:
    double A;
    double Amax;
    double c0;
    double r0, r1;
    double* wbb;
    double sigma2;
    double gamma;
    double omega;
    double rise_time;
    double decay_time;
    void print();
    void write(ofstream&,double);
    param& operator=(const param&);
    string filename;
    param();
    param(string);

};

#endif
