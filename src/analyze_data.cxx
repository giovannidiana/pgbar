#include <iostream>
#include<iomanip>
#include "../include/particle.h"
#include "../include/constants.h"
#include <fstream>
#include <sys/stat.h>
#include "../include/utils.h"
#include <getopt.h>

using namespace std;

int main(int argc, char *argv[])
{

    string data_file;
    string constants_file;
    string output_folder;
    unsigned int column;
    string gtSpike_file;
    string tag;

    string trainedPriorFile;
    bool has_trained_priors=false;
    bool has_gtspikes=false;
    bool append=false;
    int existing_samples=1;
    unsigned int niter=0;
    unsigned int trim=1;
    unsigned int verbose=1;
    int seed=0;
    int nparticles=-1;

    int options_required[] = {1,1,1,1,1,0,0,0,0,0,0,0};
    int options_provided[] = {0,0,0,0,0,0,0,0,0,0,0,0};

    int c;

    static struct option long_options[] = {
        {"data_file",required_argument,NULL,'d'},
        {"constants_file",required_argument,NULL,'c'},
        {"output_folder",required_argument,NULL,'o'},
        {"column",required_argument,NULL,'C'},
        {"tag",required_argument,NULL,'t'},
        {"niter",required_argument,NULL,'n'},
        {"prior",required_argument,NULL,'p'},
        {"continue",no_argument,NULL,'a'},
        {"trim",required_argument,NULL,'T'},
        {"quiet",no_argument,NULL,'q'},
        {"ground_truth",required_argument,NULL,'g'},
        {"rng_seed",required_argument,NULL,'s'},
        {"particles",required_argument,NULL,'j'},
        {NULL,0,NULL,0}
    };

    while (true){
        int option_index=0;
        c=getopt_long(argc,argv,"",long_options,&option_index);
        options_provided[option_index]=1;

        if(c==-1) break;

        switch(c){
            case 'd':
                data_file.assign(optarg);
                 break;
            case 'c':
                 constants_file.assign(optarg);
                 break;
            case 'o':
                 output_folder.assign(optarg);
                 // Create folder if missing
                 //
                 struct stat sb;

                 if (stat(output_folder.c_str(), &sb) != 0)
                 {
                    const int dir_err = mkdir(output_folder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
                    if (-1 == dir_err)
                    {
                        printf("Error creating directory!");
                        exit(1);
                    }
                 }

                 break;
            case 'C':
                 column=stoi(optarg);
                 break;
            case 't':
                 tag.assign(optarg);
                 break;
            case 'n':
                 niter=stoi(optarg);
                 break;
            case 'p':
                 has_trained_priors=true;
                 trainedPriorFile.assign(optarg);
                 break;
            case 'a':
                 append=true;
                 break;

            case 'q':
                 verbose=0;
                 break;

            case 'T':
                 trim=stoi(optarg);
                 break;
            
            case 'g':
                 gtSpike_file.assign(optarg);
                 has_gtspikes=true;
                 break;
            
            case 's':
                 seed=stoi(optarg);
                 break;

            case 'j':
                 nparticles = stoi(optarg);
                 break;
                 
            default:
                 cerr << "usage:" << endl;
                 cerr << "./bin/main <data> <constants> <output_folder> <column> <gtSpikes> <tag> <trainedPriorFile>" << endl;
                 return 1;
        }
    }

	// check all required options are there
	for(unsigned int i=0;i<8;i++){
		if(options_required[i]==1 && options_provided[i]==0){
			cerr<<"! missing "<<long_options[i].name<<endl;
			return 1;
		}
	}
    

    constpar constants(constants_file);

    constants.output_folder = output_folder;

    cout << "constants: " << constants_file << endl;
    cout << "output_folder: " << output_folder << endl;
    cout << "using column " << column << endl;

    ofstream parsamples;
    ofstream noiseout;
    ofstream trajsamples;
    istringstream last_params;

    if(append){
        ifstream ps(output_folder+"/param_samples_"+tag+".dat");
        string line;        
        while(ps>>std::ws && std::getline(ps,line)) existing_samples++;

        if(existing_samples==0){
            cerr << "empty existing file!!"<<endl;
            return 1;
        }
        
        last_params.str(line); 
        ps.close();

        parsamples.open(output_folder+"/param_samples_"+tag+".dat",std::ios_base::app);
        trajsamples.open(output_folder+"/traj_samples_"+tag+".dat",std::ios_base::app);
    } 

    reparam tau2gamma(5000);
    
    param testpar(output_folder+"/param_samples_"+tag+".dat");
    param testpar2(output_folder+"/param_samples_"+tag+".dat");
    noiseout.open(output_folder+"/noise_"+tag+".dat");

    arma::vec gtSpikes;
    if(has_gtspikes){
        constants.KNOWN_SPIKES=true;
        gtSpikes.load(gtSpike_file,arma::raw_ascii);
    }

    if(has_trained_priors){
        string dum;
        cout<<"Using trained priors:"<<endl;
        // update the constants
        ifstream trainedPrior(trainedPriorFile);
        trainedPrior>>dum>>constants.Amax_prior_mean>>constants.Amax_prior_sd;
        trainedPrior>>dum>>dum>>constants.c0_prior_sd;
        trainedPrior>>dum>>constants.decay_time_mean_ms>>constants.decay_time_sd;
        trainedPrior>>dum>>constants.rise_time_mean_ms>>constants.rise_time_sd;
        trainedPrior>>dum>>constants.alpha_sigma2>>constants.beta_sigma2;
        trainedPrior>>dum>>constants.alpha_rate_b0>>constants.beta_rate_b0;
        trainedPrior>>dum>>constants.alpha_rate_b1>>constants.beta_rate_b1;
        trainedPrior>>dum>>constants.alpha_w01>>constants.beta_w01;
        trainedPrior>>dum>>constants.alpha_w10>>constants.beta_w10;
        trainedPrior.close();
    } 

    if(niter>0){
        constants.niter=niter;
    }

    if(nparticles>0){
        constants.nparticles=nparticles;
    }

    // build scheduler for bm_sigma
    arma::vec log10_bm_sigma_schedule1 = arma::linspace(-3,log10(constants.bm_sigma),constants.niter/2);
    arma::vec log10_bm_sigma_schedule2(constants.niter/2); 
    log10_bm_sigma_schedule2.fill( log10(constants.bm_sigma) );
    arma::vec log10_bm_sigma_schedule = arma::join_cols(log10_bm_sigma_schedule1, log10_bm_sigma_schedule2);

    // Initialize the sampler ( this will also reset the scales, that's why we need to initialize after we update the constants )
    SMC sampler(data_file,column,constants,false,seed);

    // Initialize the trajectory
    Trajectory traj_sam1(sampler.TIME,""), traj_sam2(sampler.TIME,output_folder+"/traj_samples_"+tag+".dat");

    for(unsigned int t=0;t<sampler.TIME;++t){
        traj_sam1.B(t)=0;
        traj_sam1.burst(t)=0;
        traj_sam1.C(t)=0;
        traj_sam1.S(t)=0;
        if(has_gtspikes) traj_sam1.S(t)=gtSpikes(t);
        traj_sam1.Y(t)=0;
    }

    // set initial parameters 

    if(append){
        cout<<"Use parameters from previous analysis"<<endl;
        vector<string> parse_params;
        while(last_params.good()){
            string substr;
            getline(last_params,substr,',');
            parse_params.push_back(substr);
        }

        testpar.Amax       = stod(parse_params[0]);
        testpar.c0         = stod(parse_params[1]);
        testpar.decay_time = stod(parse_params[2]);
        testpar.rise_time  = stod(parse_params[3]);
        testpar.sigma2     = stod(parse_params[4]);
        testpar.r0         = stod(parse_params[5]);
        testpar.r1         = stod(parse_params[6]);
        testpar.wbb[0]     = stod(parse_params[7]);
        testpar.wbb[1]     = stod(parse_params[8]);

        testpar.decay_time*=constants.sampling_frequency/1000.0;
        testpar.rise_time*=constants.sampling_frequency/1000.0;
    
    } else {
        cout<<"Draw new parameters"<<endl;
        testpar.c0=0;
        testpar.r0=constants.alpha_rate_b0/constants.beta_rate_b0;
        testpar.r1=constants.alpha_rate_b1/constants.beta_rate_b1;
        testpar.wbb[0]=constants.alpha_w01/constants.beta_w01;
        testpar.wbb[1]=constants.alpha_w10/constants.beta_w10;
        testpar.sigma2=constants.beta_sigma2/(constants.alpha_sigma2+1); // the mode of the inverse gamma
        testpar.decay_time = constants.decay_time_mean;
        testpar.rise_time = constants.rise_time_mean;
        testpar.Amax=constants.Amax_prior_mean;
    }

    tau2gamma.map(testpar.Amax, testpar.decay_time,testpar.rise_time,testpar.A,testpar.gamma,testpar.omega);

    for(unsigned int i=existing_samples-1;i<constants.niter ;i++){

        // set bm_sigma from the scheduler
        //
        //constants.bm_sigma = pow(10.0,log10_bm_sigma_schedule[i]);

        sampler.PGAS(testpar,traj_sam1,traj_sam2);
        traj_sam1=traj_sam2;
        for(unsigned int k=0;k<10;k++){
            sampler.sampleParameters(testpar,testpar2,traj_sam2, &tau2gamma);
            testpar=testpar2;
        }
        if(i % trim == 0){
            testpar.write(parsamples,constants.sampling_frequency);
            traj_sam2.write(trajsamples,i/trim);

            // write noise hyperparameters
            double res=pow(arma::norm(sampler.data_y - traj_sam2.B-traj_sam2.C),2);
            noiseout << constants.alpha_sigma2 + traj_sam2.size / 2.0 << ' ' << constants.beta_sigma2 + 0.5 * res << endl;
        }
        if(verbose>0) cout<<"iteration:"<<setw(5)<<i<<", spikes: "<<setw(5)<<arma::accu(traj_sam1.S)<<'\r'<<flush;
    }

    parsamples.close();
    trajsamples.close();

    // Save last param set
    //
    ofstream lastParams(output_folder+"/last_params_"+tag+".dat");
    lastParams<<testpar.c0<<" "<<testpar.r0<<" "<<testpar.r1<<" "<<testpar.wbb[0]<<" "<<testpar.wbb[1]<<" "<<testpar.sigma2<<" "<<testpar.decay_time<<" "<<testpar.rise_time<<endl;


    return 0;
}
