#include <iostream>
#include "../include/particle.h"
#include "../include/constants.h"
#include <fstream>
#include <sys/stat.h>
#include "../include/utils.h"

using namespace std;

int main(int argc, char *argv[])
{
    if(argc<2){
        cout<<"$1: constants json file"<<endl;
        cout<<"$2: gt values (X only. Y is saved in trace.dat)"<<endl;
        cout<<"$3: PGAS samples"<<endl;
        cout<<"$4: param samples"<<endl;
        exit(0);
    }

    ofstream gtfile;
    ofstream parfile;
    ofstream trajsamples;
    ofstream correlation_file("correlation.dat");

    gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);
    constpar constants(argv[1]);

    param testpar;
    testpar.Amax=1;
    testpar.c0=0;
    testpar.r0=0.1; testpar.r1=4;
    testpar.wbb[0] = .1;
    testpar.wbb[1] = .2;
    testpar.sigma2=pow(0.2,2);
    testpar.decay_time = 1000/1000.0 * constants.sampling_frequency;
    testpar.rise_time = 80/1000.0 * constants.sampling_frequency;

    testpar.print();

    reparam tau2gamma(1000);
    tau2gamma.map(testpar.Amax, testpar.decay_time,testpar.rise_time,testpar.A,testpar.gamma,testpar.omega);

    int length=4000;
    Trajectory T(length, argv[2]);
    Trajectory traj_sam1(length,""), traj_sam2(length,argv[3]);

    for(unsigned int t=0;t<T.size;++t){
        traj_sam1.B(t)=0;
        traj_sam1.burst(t)=0;
        traj_sam1.C(t)=0;
        traj_sam1.S(t)=0;
        traj_sam1.Y(t)=0;
    }

    T.simulate(rng, testpar, &constants);
    T.write(gtfile,0);
    T.Y.save("trace.dat",arma::raw_ascii);

    SMC sampler(T.Y,constants,true);
    param testpar2(argv[4]);
    
    for(unsigned int i=0;i<constants.niter;i++){
        double sscor = utils::subsampled_correlation(traj_sam1.S,T.S,30,20);
        correlation_file<<sscor<<endl;
        sampler.PGAS(testpar,traj_sam1,traj_sam2);
        sampler.sampleParameters(testpar,testpar2,traj_sam2, &tau2gamma);
        traj_sam1=traj_sam2;
        testpar=testpar2;
        testpar2.write(parfile,constants.sampling_frequency);
        traj_sam2.write(trajsamples,i);
    }

    parfile.close();
    trajsamples.close();
    gtfile.close();


    return 0;
}
