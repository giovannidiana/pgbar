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
    bool append=false;
    int existing_samples=1;
    unsigned int niter=0;

	int options_required[]={1,1,1,1,1,1,0,0,0,0};
	int options_provided[]={0,0,0,0,0,0,0,0,0,0};
    
    int c;
    static struct option long_options[] = {
        {"data_file",required_argument,NULL,'d'},
        {"constants_file",required_argument,NULL,'c'},
        {"output_folder",required_argument,NULL,'o'},
        {"column",required_argument,NULL,'C'},
        {"ground_truth",required_argument,NULL,'g'},
        {"tag",required_argument,NULL,'t'},
        {"niter",required_argument,NULL,'n'},
        {"prior",required_argument,NULL,'p'},
        {"continue",no_argument,NULL,'a'},
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
            case 'g':
                 gtSpike_file.assign(optarg);
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
    ofstream trajsamples;
    ofstream correlation_file;
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
        correlation_file.open(output_folder+"/correlations_"+tag+".dat",std::ios_base::app);
    } else {
        correlation_file.open(output_folder+"/correlations_"+tag+".dat");
    }


    arma::vec gtSpikes; gtSpikes.load(gtSpike_file,arma::raw_ascii);

    reparam tau2gamma(1000);

    int length=gtSpikes.n_elem;
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

    // Initialize the sampler ( this will also reset the scales, that's why we need to initialize after we update the constants )
    SMC sampler(data_file,column,constants,false);

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

    for(unsigned int i=existing_samples-1;i<constants.niter;i++){
        //double sscor = utils::subsampled_correlation(traj_sam1.S,gtSpikes,constants.sampling_frequency,7.5);
        sampler.PGAS(testpar,traj_sam1,traj_sam2);
        traj_sam1=traj_sam2;
        sampler.sampleParameters(testpar,testpar2,traj_sam2, &tau2gamma);
        testpar=testpar2;
        testpar.write(parsamples,constants.sampling_frequency);
        traj_sam2.write(trajsamples,i);
        double sscor = utils::subsampled_and_filtered_correlation(traj_sam2.S,gtSpikes,constants.sampling_frequency,7.5,200);
        correlation_file<<sscor<<' '<<arma::accu(traj_sam1.S)<<endl;
        cout<<"iteration:"<<setw(5)<<i<<" correlation with GT: "<<setw(10)<<sscor<<", spikes: "<<setw(5)<<arma::accu(traj_sam1.S)<<setw(10)<<arma::accu(gtSpikes)<<'\r'<<flush;
    }

    parsamples.close();
    trajsamples.close();

    // Save last param set
    //
    ofstream lastParams(output_folder+"/last_params_"+tag+".dat");
    lastParams<<testpar.c0<<" "<<testpar.r0<<" "<<testpar.r1<<" "<<testpar.wbb[0]<<" "<<testpar.wbb[1]<<" "<<testpar.sigma2<<" "<<testpar.decay_time<<" "<<testpar.rise_time<<endl;


    return 0;
}
