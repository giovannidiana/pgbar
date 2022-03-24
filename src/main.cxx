#include <iostream>
#include<iomanip>
#include "../include/particle.h"
#include "../include/constants.h"
#include <fstream>
#include <sys/stat.h>
#include "../include/utils.h"

using namespace std;

int main(int argc, char *argv[])
{

    if (argc < 2)
    {
        cout << "usage:" << endl;
        cout << "./bin/main <data> <constants> <output_folder> <column> <gtSpikes> <tag> <initParamFile>" << endl;
        return(0);
    }

    string data_file = argv[1];
    string constants_file = argv[2];
    string output_folder = argv[3];
    unsigned int column=atoi(argv[4]);
    string gtSpike_file = argv[5];
    string tag = argv[6];

    string initParamFile;
    bool has_init_param=false;
    if(argc>7){
        has_init_param=true;
        initParamFile = argv[7];
    }

    constpar constants(constants_file);

    constants.output_folder = output_folder;

    cout << "constants: " << constants_file << endl;
    cout << "output_folder: " << output_folder << endl;
    cout << "using column " << column << endl;

    // Create folder if missing
    //
    struct stat sb;

    if (stat(argv[3], &sb) != 0)
    {
        const int dir_err = mkdir(argv[3], S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        if (-1 == dir_err)
        {
            printf("Error creating directory!");
            exit(1);
        }
    }

    ofstream parsamples;
    ofstream trajsamples;
    ofstream correlation_file(output_folder+"/correlations_"+tag+".dat");
    arma::vec gtSpikes; gtSpikes.load(gtSpike_file,arma::raw_ascii);

    SMC sampler(data_file,column,constants,false);

    reparam tau2gamma(1000);

    int length=sampler.data_y.n_elem;
    Trajectory traj_sam1(length,""), traj_sam2(length,output_folder+"/traj_samples_"+tag+".dat");

    for(unsigned int t=0;t<length;++t){
        traj_sam1.B(t)=0;
        traj_sam1.burst(t)=0;
        traj_sam1.C(t)=0;
        traj_sam1.S(t)=0;
        traj_sam1.Y(t)=0;
    }
    
    param testpar(output_folder+"/param_samples_"+tag+".dat");
    param testpar2(output_folder+"/param_samples_"+tag+".dat");

    if(has_init_param){
        cout<<"Using parameters:"<<endl;
        ifstream initParam(initParamFile);
        initParam>>testpar.c0>>testpar.r0>>testpar.r1>>testpar.wbb[0]>>testpar.wbb[1]>>testpar.sigma2>>testpar.decay_time>>testpar.rise_time;
        testpar.print();
        initParam.close();
    } else {
        cout<<"Draw new parameters"<<endl;
        testpar.c0=0;
        testpar.r0=0.1;
        testpar.r1=1;
        testpar.wbb[0]=0.1;
        testpar.wbb[1]=0.1;
        testpar.sigma2=pow(0.1,2);
        testpar.decay_time = constants.decay_time_mean;
        testpar.rise_time = constants.rise_time_mean;
        testpar.Amax=constants.Amax_prior_mean;
    }

    tau2gamma.map(testpar.Amax, testpar.decay_time,testpar.rise_time,testpar.A,testpar.gamma,testpar.omega);

    for(unsigned int i=0;i<constants.niter ;i++){
        double sscor = utils::subsampled_correlation(traj_sam1.S,gtSpikes,100,25);
        correlation_file<<sscor<<' '<<arma::accu(traj_sam1.S)<<endl;
        cout<<"iteration:"<<setw(5)<<i<<" correlation with GT: "<<setw(10)<<sscor<<", spikes: "<<setw(5)<<arma::accu(traj_sam1.S)<<setw(5)<<arma::accu(gtSpikes)<<"\r"<<flush;
        sampler.PGAS(testpar,traj_sam1,traj_sam2);
        traj_sam1=traj_sam2;
        sampler.sampleParameters(testpar,testpar2,traj_sam2, &tau2gamma);
        testpar=testpar2;
        testpar.write(parsamples,constants.sampling_frequency);
        traj_sam2.write(trajsamples,i);
    }

    parsamples.close();
    trajsamples.close();

    // Save last param set
    //
    ofstream lastParams(output_folder+"/last_params_"+tag+".dat");
    lastParams<<testpar.c0<<" "<<testpar.r0<<" "<<testpar.r1<<" "<<testpar.wbb[0]<<" "<<testpar.wbb[1]<<" "<<testpar.sigma2<<" "<<testpar.decay_time<<" "<<testpar.rise_time<<endl;


    return 0;
}
