#include <iostream>
#include "../include/particle.h"
#include "../include/constants.h"
#include <fstream>
#include <sys/stat.h>
#include "../include/utils.h"
#include<sstream>

using namespace std;

int main(int argc, char *argv[])
{
    if(argc<2){
        cout<<"$1: folder"<<endl;
        cout<<"$2: ISI"<<endl;
        cout<<"$3: SNR"<<endl;
        cout<<"$4: increment"<<endl;
        cout<<"$5: constants"<<endl;
        exit(0);
    }

    stringstream ss;
    string output_folder=argv[1];
    double ISI = atoi(argv[2]);
    double sigma = 1./atof(argv[3]);
    int increment = atoi(argv[4]);

    ss<< argv[2]<<'_'<<argv[3]<<'_'<<argv[4];
    string label = ss.str();

    gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng,(int)(100*atof(argv[3]))+increment);

    constpar constants(argv[5]);
    constants.sampling_frequency = 600;

    param testpar(output_folder+"/"+label+"-params.csv");
    testpar.Amax=constants.Amax_prior_mean;
    testpar.c0=0;
    testpar.r0=0.1; testpar.r1=1;
    testpar.wbb[0] = .1;
    testpar.wbb[1] = .2;
    testpar.sigma2=pow(sigma,2);
    testpar.decay_time = 40/1000.0 * constants.sampling_frequency;
    testpar.rise_time = 3.7/1000.0 * constants.sampling_frequency;

    reparam tau2gamma(1000);
    tau2gamma.map(testpar.Amax, testpar.decay_time,testpar.rise_time,testpar.A,testpar.gamma,testpar.omega);

    ofstream paramfile;
    testpar.write(paramfile,constants.sampling_frequency);

    int length=floor(constants.sampling_frequency*0.2);
    cout<<output_folder+"/"+label+"-gt.csv"<<endl;
    Trajectory T(length, output_folder+"/"+label+"-gt.csv");

    T.simulate_doubles2(rng, testpar, &constants, ISI);
    
    ofstream gtfile;
    T.write(gtfile,0);

    arma::mat datamat(length,2);
    datamat.col(0) = arma::linspace(0,length-1,length)/constants.sampling_frequency;
    datamat.col(1) = T.Y;
    
    datamat.save(output_folder+"/"+label+"-data.csv",arma::raw_ascii);
    T.S.save(output_folder+"/"+label+"-spikes.csv",arma::raw_ascii);

    gtfile.close();


    return 0;
}
