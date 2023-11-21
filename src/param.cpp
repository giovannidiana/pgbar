#include "../include/param.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <iostream>

using namespace std;

param::param(){
    wbb = new double[2];
}

param::param(string fname){
    wbb = new double[2];
    filename=fname;
}

void param::print(){
    cout<<"Amax   = "<<Amax<<endl
        <<"c0     = "<<c0<<endl
        <<"taud   = "<<decay_time<<endl
        <<"taur   = "<<rise_time<<endl
        <<"sigma2 = "<<sigma2<<endl
        <<"r0     = "<<r0<<endl
        <<"r1     = "<<r1<<endl;
}

void param::write(ofstream &out,double sf){
    if(!out.is_open()){
        out.open(filename);
        out<<"Amax,c0,decay_time,rise_time,sigma2,r0,r1,w01,w10"<<endl;
    }

    out<<Amax<<","
        <<c0<<","
        <<decay_time/sf*1000<<","
        <<rise_time/sf*1000<<","
        <<sigma2<<","
        <<r0<<","
        <<r1<<","
        <<wbb[0]<<","
        <<wbb[1]<<endl;
}

param& param::operator=(const param& p){
    Amax = p.Amax;
    A    = p.A;
    c0   = p.c0;
    r0   = p.r0;
    r1   = p.r1;
    wbb[0] = p.wbb[0];
    wbb[1] = p.wbb[1];
    sigma2 = p.sigma2;
    gamma = p.gamma;
    omega = p.omega;
    decay_time = p.decay_time;
    rise_time = p.rise_time;
    return *this;
}
