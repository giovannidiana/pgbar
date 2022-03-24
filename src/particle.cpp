#include"../include/particle.h"
#include"../include/utils.h"
#include<gsl/gsl_sf_gamma.h>
#include<gsl/gsl_math.h>
#include<fstream>
#include<string>
#include<ctime>

using namespace std;

// PARTICLE CLASS
// -----------------------------------------------------------
//
Particle::Particle(){
    C.set_size(2);
    B=0;
    burst=0;
    S=0;
    ancestor=-1;
    logWeight=0;
}

Particle& Particle::operator=(const Particle& p){
    B=p.B;
    burst=p.burst;
    C=p.C;
    S=p.S;
    logWeight=p.logWeight;
    return *this;
}

void Particle::print(){
    cout<<"    C    : "<<C(0)<<endl
        <<"    B    : "<<B<<endl
        <<"    burst: "<<burst<<endl
        <<"    S    : "<<S<<endl;
}

// Trajectory class
// ------------------------------------------------------------
//
Trajectory::Trajectory(unsigned int s, string fname){

    filename=fname;
    size=s;
    B.resize(s);
    burst.resize(s);
    C.resize(s);
    S.resize(s);
    Y.resize(s);

}

void Trajectory::simulate(gsl_rng* rng, const param& par, const constpar *constants){

    double dt=1.0/constants->sampling_frequency;
    double rate[2] = {par.r0*dt,par.r1*dt};

    B(0) = gsl_ran_gaussian(rng,constants->bm_sigma);
    burst(0) = 0;
    S(0) = gsl_ran_poisson(rng,par.r0*dt);
    C(0) = par.c0 + par.A*S(0);
    Y(0) = C(0)+B(0)+gsl_ran_gaussian(rng,sqrt(par.sigma2));
    
    double W[2][2] = {{1-par.wbb[0]*dt, par.wbb[0]*dt},
        {par.wbb[1]*dt, 1-par.wbb[1]*dt}};

    for(unsigned int t=1;t<size;t++){
        burst(t) = (gsl_rng_uniform(rng) < W[burst(t-1)][0]) ? 0 : 1;

        B(t)     = B(t-1)+gsl_ran_gaussian(rng,constants->bm_sigma*sqrt(dt)); 
        S(t)     = gsl_ran_poisson(rng,rate[burst(t)]);
        C(t)     = (t>1) ? C(t-2)*par.omega + C(t-1)*par.gamma + par.A*S(t) : C(t-1)*par.gamma + par.A*S(t);
        Y(t)     = C(t)+B(t)+gsl_ran_gaussian(rng,sqrt(par.sigma2)); 
    }
}

void Trajectory::simulate_doubles(gsl_rng* rng, const param& par, const constpar *constants){

    double dt=1.0/constants->sampling_frequency;
    double rate[2] = {par.r0*dt,par.r1*dt};

    B(0) = gsl_ran_gaussian(rng,constants->bm_sigma);
    burst(0) = 0;
    S(0) = 0;
    C(0) = par.c0 + par.A*S(0);
    Y(0) = C(0)+B(0)+gsl_ran_gaussian(rng,sqrt(par.sigma2));
    
    double W[2][2] = {{1-par.wbb[0]*dt, par.wbb[0]*dt},
        {par.wbb[1]*dt, 1-par.wbb[1]*dt}};

    bool stim1=false;
    bool stim2=false;

    for(unsigned int t=1;t<size;t++){
        burst(t) = (t*dt>=0.09 && t*dt<0.111) ? 1 : 0;

        B(t)     = B(t-1)+gsl_ran_gaussian(rng,constants->bm_sigma*sqrt(dt)); 
        S(t)     = 0;
        
        if(!stim1 && t*dt>=0.1) {
            S(t) = 1;
            stim1=true;
        } 

        if(!stim2 && t*dt>=0.11) {
            S(t) = 1;
            stim2=true;
        }

        C(t)     = (t>1) ? C(t-2)*par.omega + C(t-1)*par.gamma + par.A*S(t) : C(t-1)*par.gamma + par.A*S(t);
        Y(t)     = C(t)+B(t)+gsl_ran_gaussian(rng,sqrt(par.sigma2)); 
    }
}

// simulate double stim with variable inter spike interval
void Trajectory::simulate_doubles2(gsl_rng* rng, const param& par, const constpar *constants, int isi){

    double dt=1.0/constants->sampling_frequency;
    double rate[2] = {par.r0*dt,par.r1*dt};

    B(0) = gsl_ran_gaussian(rng,constants->bm_sigma);
    burst(0) = 0;
    S(0) = 0;
    C(0) = par.c0 + par.A*S(0);
    Y(0) = C(0)+B(0)+gsl_ran_gaussian(rng,sqrt(par.sigma2));
    
    double W[2][2] = {{1-par.wbb[0]*dt, par.wbb[0]*dt},
        {par.wbb[1]*dt, 1-par.wbb[1]*dt}};

    bool stim1=false;
    bool stim2=false;

    for(unsigned int t=1;t<size;t++){
        burst(t) = (t*dt>=0.09 && t*dt<0.111) ? 1 : 0;

        B(t)     = B(t-1)+gsl_ran_gaussian(rng,constants->bm_sigma*sqrt(dt)); 
        S(t)     = 0;
        
        if(!stim1 & t>=20) {
            S(t) = 1;
            stim1=true;
        } 

        if(!stim2 && t>=20+isi) {
            S(t) = 1;
            stim2=true;
        }

        C(t)     = (t>1) ? C(t-2)*par.omega + C(t-1)*par.gamma + par.A*S(t) : C(t-1)*par.gamma + par.A*S(t);
        Y(t)     = C(t)+B(t)+gsl_ran_gaussian(rng,sqrt(par.sigma2)); 
    }
}

void Trajectory::write(ofstream &outfile, unsigned int index){
    
    if(!outfile.is_open()){
        outfile.open(filename);
        outfile<<"index,burst,B,S,C"<<endl;
    }

    for(unsigned int i=0;i<size;++i){
        outfile<<index<<','
            <<burst(i)<<','
            <<B(i)<<','
            <<S(i)<<','
            <<C(i)<<endl;
    }

}

Trajectory& Trajectory::operator=(const Trajectory& traj){
        this->burst = traj.burst;
        this->B     = traj.B;
        this->S     = traj.S;
        this->C     = traj.C;
        this->Y     = traj.Y;

        return *this ;
}

// SMC class
// -------------------------------------------------------------------
//
SMC::SMC(string filename, int index, constpar& cst, bool has_header, int seed){

    rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng, (seed==0) ? cst.seed : seed);

    arma::field <string> header;

    if(has_header){
        tracemat.load(arma::csv_name(filename, header));
    } else {
        tracemat.load(filename);
    }

    data_time  = tracemat.col(0);
    data_y    = tracemat.col(index);
    constants = &cst;

    // The time units here are assumed to be seconds
    constants->sampling_frequency = 1.0 / (tracemat(1, 0)-tracemat(0,0));
    cout << "setting sampling frequency to: "<<constants->sampling_frequency << endl;
    constants->set_time_scales();

    // set number of particles
    nparticles=constants->nparticles;
    TIME = tracemat.n_rows;
    cout<<"nparticles: "<<nparticles<<endl;
    cout<<"TIME      : "<<TIME<<endl;

    // set particle system
    particleSystem.resize(TIME);
    for(unsigned int t=0;t<TIME;++t){
        particleSystem[t] = new Particle[nparticles];
        for(unsigned int j=0;j<nparticles;j++) particleSystem[t][j].index=j;
    }
}

SMC::SMC(arma::vec &Y, constpar& cst, bool v){

    verbose=v;

    constants = &cst;

    rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng,constants->seed);

    data_time  = arma::regspace(0,Y.n_elem,1.0/constants->sampling_frequency);
    data_y    = Y;

    // The time units here are assumed to be seconds
    if(verbose) cout << "setting sampling frequency to: "<<constants->sampling_frequency << endl;
    constants->set_time_scales();

    // set number of particles
    nparticles=constants->nparticles;
    TIME = Y.n_rows;

    // set particle system
    particleSystem.resize(TIME);
    for(unsigned int t=0;t<TIME;++t){
        particleSystem[t] = new Particle[nparticles];
        for(unsigned int j=0;j<nparticles;j++) particleSystem[t][j].index=j;
    }
}


void SMC::rmu(Particle &p, double y, const param &par, bool set=false){

    if(!set){
        // draw state if not set

        p.burst = (gsl_rng_uniform(rng) > 0.5 ) ? 1 : 0;

        if(!constants->KNOWN_SPIKES){
            double rate=(p.burst==0) ? par.r0 : par.r1;
            p.S     = gsl_ran_poisson(rng, rate/constants->sampling_frequency);
        }

        p.C     = {par.c0 + par.A*p.S, 0};
        p.B     = gsl_ran_gaussian(rng,1);
    }

    p.logWeight=-0.5*log(2*M_PI*par.sigma2)-0.5/par.sigma2*pow(y-p.C(0)-p.B,2); // from g(y|x) only
                
}
double SMC::logf(const Particle &pin, const Particle &pout, const param &par){
    double lf = 0;
    double dt = 1.0/constants->sampling_frequency;
    double W[2][2] = {{1-par.wbb[0]*dt, par.wbb[0]*dt},
        {par.wbb[1]*dt, 1-par.wbb[1]*dt}};
    double rate[2] = {par.r0*dt,par.r1*dt};


    lf += log(W[pin.burst][pout.burst]) +
        pout.S*log(rate[pout.burst])-rate[pout.burst] - log(gsl_sf_gamma(pout.S+1)) +
        -0.5*pow((pout.B-pin.B)/(constants->bm_sigma*sqrt(dt)),2);

    return(lf);
}

// Here part.S is supposed to be given.
//
void SMC::move_and_weight_GTS(Particle &part, const Particle& parent, double y, const param &par, bool set=false){

    double probs[2];
    double log_probs[2];
    double Z;
    double ct;
    unsigned int burst;
    unsigned int ns=part.S;

    double dt = 1.0/constants->sampling_frequency;

    double rate[2] = {par.r0*dt,par.r1*dt};
    double W[2][2] = {{1-par.wbb[0]*dt, par.wbb[0]*dt},
        {par.wbb[1]*dt, 1-par.wbb[1]*dt}};

    double z[2] = {par.sigma2/(par.sigma2+dt*pow(constants->bm_sigma,2)),
                   dt*pow(constants->bm_sigma,2)/(par.sigma2+dt*pow(constants->bm_sigma,2))};
    double sigma_B_posterior = sqrt(dt*pow(constants->bm_sigma,2)*par.sigma2/(par.sigma2+dt*pow(constants->bm_sigma,2)));
    
    arma::mat M = { {par.gamma, par.omega},
                    {1,0}};

    arma::mat MC=M*parent.C;
    ct = MC(0) + par.A*ns;

    log_probs[0] = log(W[parent.burst][0]);
    log_probs[0] += ns*log(rate[0]) - log(gsl_sf_gamma(ns+1)) -rate[0];
    log_probs[0] += -0.5/(par.sigma2+pow(constants->bm_sigma,2))*pow(y-ct-parent.B,2);
    log_probs[1] = log(W[parent.burst][1]);
    log_probs[1] += ns*log(rate[1]) - log(gsl_sf_gamma(ns+1)) -rate[1];
    log_probs[1] += -0.5/(par.sigma2+pow(constants->bm_sigma,2))*pow(y-ct-parent.B,2);

    utils::w_from_logW(log_probs,probs,2);
    Z=utils::Z_from_logW(log_probs,2);

    part.logWeight = log(Z);

    if(!set){
        // move particle if not already set
        gsl_ran_discrete_t *rdisc = gsl_ran_discrete_preproc(2, probs);
        int idx = gsl_ran_discrete(rng,rdisc);

        part.burst = idx;
        part.B     = z[0]*parent.B+z[1]*(y-parent.C(0)) + gsl_ran_gaussian(rng,sigma_B_posterior);

        gsl_ran_discrete_free(rdisc);
    }
    
    arma::vec svec = {ns,0.0};
    part.C     = MC+par.A*svec;

}

void SMC::move_and_weight(Particle &part, const Particle& parent, double y, const param &par, bool set=false){

    const int maxspikes = 10;  // The number of spikes goes from 0 to maxspikes-1
    double probs[2*maxspikes];
    double log_probs[2*maxspikes];
    double Z;
    double ct;
    unsigned int burst, ns;
    double dt = 1.0/constants->sampling_frequency;

    double rate[2] = {par.r0*dt,par.r1*dt};
    double W[2][2] = {{1-par.wbb[0]*dt, par.wbb[0]*dt},
        {par.wbb[1]*dt, 1-par.wbb[1]*dt}};
    
    double z[2] = {par.sigma2/(par.sigma2+dt*pow(constants->bm_sigma,2)),
                   dt*pow(constants->bm_sigma,2)/(par.sigma2+dt*pow(constants->bm_sigma,2))};
    double sigma_B_posterior = sqrt(dt*pow(constants->bm_sigma,2)*par.sigma2/(par.sigma2+dt*pow(constants->bm_sigma,2)));
    
    arma::mat M = { {par.gamma, par.omega},
                    {1,0}};

    arma::mat MC=M*parent.C;

    for(unsigned int i=0;i<2*maxspikes;i++){

        burst = floor(i/maxspikes);
        ns    = i%maxspikes;
        ct    = MC(0)+par.A*ns;

        log_probs[i] = log(W[parent.burst][burst]);
        log_probs[i] += ns*log(rate[burst]) - log(gsl_sf_gamma(ns+1)) -rate[burst];
        log_probs[i] += -0.5/(par.sigma2+pow(constants->bm_sigma,2))*pow(y-ct-parent.B,2);

    } 

    utils::w_from_logW(log_probs,probs,2*maxspikes);
    Z=utils::Z_from_logW(log_probs,2*maxspikes);

    part.logWeight = log(Z);

    if(!set){
        // move particle if not already set
        gsl_ran_discrete_t *rdisc = gsl_ran_discrete_preproc(2*maxspikes, probs);
        int idx = gsl_ran_discrete(rng,rdisc);

        part.burst = floor(idx/maxspikes);
        part.S     = idx%maxspikes;
        part.B     = z[0]*parent.B+z[1]*(y-parent.C(0)) + gsl_ran_gaussian(rng,sigma_B_posterior);

        gsl_ran_discrete_free(rdisc);
    }
    
    arma::vec svec = {part.S,0.0};
    part.C     = MC+par.A*svec;

}

void SMC::sampleParameters(const param &pin, param &pout, Trajectory &traj, reparam *tau2gamma){

    // Sampling the posterior burst transition rates
    //
    int counts[2][2] = {{0,0},{0,0}};

    for(unsigned int t=1;t<traj.size;t++){
        counts[traj.burst(t-1)][traj.burst(t)]++;
    }
    double dt = 1.0/constants->sampling_frequency;

    double alpha_w01_post = constants->alpha_w01 + counts[0][1];
    double beta_w01_post  = constants->beta_w01  + counts[0][0]*dt;
    double alpha_w10_post = constants->alpha_w10 + counts[1][0];
    double beta_w10_post  = constants->beta_w10  + counts[1][1]*dt;

    pout.wbb[0] = gsl_ran_gamma(rng,alpha_w01_post,1.0/beta_w01_post);
    pout.wbb[1] = gsl_ran_gamma(rng,alpha_w10_post,1.0/beta_w10_post);

    // Sampling the posterior firing rates
    //
    arma::uvec which_b0 = arma::find(traj.burst==0);
    double alpha_rate_b0_post = constants->alpha_rate_b0 + arma::accu(traj.S(which_b0));
    double beta_rate_b0_post  = constants->beta_rate_b0  + which_b0.n_elem*dt;

    arma::uvec which_b1 = arma::find(traj.burst==1);
    double alpha_rate_b1_post = constants->alpha_rate_b1 + arma::accu(traj.S(which_b1));
    double beta_rate_b1_post  = constants->beta_rate_b1  + which_b1.n_elem*dt;

    pout.r0 = gsl_ran_gamma(rng,alpha_rate_b0_post,1.0/beta_rate_b0_post);
    pout.r1 = gsl_ran_gamma(rng,alpha_rate_b1_post,1.0/beta_rate_b1_post);

    // Sampling standard deviation
    //
    double res = pow(arma::norm(data_y - traj.B - traj.C),2);
    pout.sigma2 = 1.0 / gsl_ran_gamma(rng, constants->alpha_sigma2 + traj.size / 2.0, 1.0 / (constants->beta_sigma2 + 0.5 * res));

    // Sampling Amax, rise and decay times

    unsigned int sampled_variable = gsl_rng_uniform_int(rng,4);
    if(!constants->SAMPLE_KINETICS) sampled_variable=sampled_variable % 2;
    unsigned int var_selec[4] = {0,0,0,0}; 
    var_selec[sampled_variable] = 1;

    if(constants->KNOWN_SPIKES){
        var_selec[0]=1;
        var_selec[1]=1;
        if(constants->SAMPLE_KINETICS){
            var_selec[2]=1;
            var_selec[3]=1;
        }
    }

    // Sampling c0 with positivity constraint
    double c0_test   = pin.c0         + var_selec[0]*gsl_ran_gaussian(rng,constants->c0_proposal_sd);
    while(var_selec[0]==1 && c0_test<0){
        c0_test   = pin.c0         + var_selec[0]*gsl_ran_gaussian(rng,constants->c0_proposal_sd);
    }
    
    // Sampling Amax
    double Amax_test = pin.Amax + var_selec[1]*gsl_ran_gaussian(rng,constants->Amax_proposal_sd);
    while(Amax_test<=0)
        Amax_test = pin.Amax + var_selec[1]*gsl_ran_gaussian(rng,constants->Amax_proposal_sd);

    double taud_test = pin.decay_time + var_selec[2]*gsl_ran_gaussian(rng,constants->decay_time_proposal_sd);

    double taur_test = pin.rise_time  + var_selec[3]*gsl_ran_gaussian(rng,constants->rise_time_proposal_sd);

    while(taud_test<=taur_test || taur_test<=0){
        taud_test = pin.decay_time + var_selec[2]*gsl_ran_gaussian(rng,constants->decay_time_proposal_sd);
        taur_test = pin.rise_time  + var_selec[3]*gsl_ran_gaussian(rng,constants->rise_time_proposal_sd);
    }


    double res_test;

    //cout<<"proposed values:"<<endl
    //   <<"    "<<Amax_test<<' '<<taud_test/constants->sampling_frequency*1000.<<' '<<taur_test/constants->sampling_frequency*1000<<' '<<c0_test<<endl;
    double A_test,gamma_test,omega_test;
    arma::vec C_test(traj.size);
    
    tau2gamma->map(Amax_test, taud_test, taur_test, A_test, gamma_test, omega_test);

    C_test(0) = c0_test + A_test*traj.S(0);
    C_test(1) = gamma_test*C_test(0) + A_test*traj.S(1);
    for(unsigned int t=2;t<traj.size;t++) C_test(t) = C_test(t-2)*omega_test + C_test(t-1)*gamma_test + A_test*traj.S(t);

    res_test = pow(arma::norm(data_y - traj.B - C_test),2);
    //cout<<"    res test: "<<res_test<<endl;

    double log_alpha_MH = -0.5*pow((Amax_test-constants->Amax_prior_mean)/constants->Amax_prior_sd,2) + 0.5*pow((pin.Amax-constants->Amax_prior_mean)/constants->Amax_prior_sd,2) 
                          -0.5*pow((taud_test-constants->decay_time_mean)/constants->decay_time_sd,2) + 0.5*pow((pin.decay_time-constants->decay_time_mean)/constants->decay_time_sd,2)
                          -0.5*pow((taur_test-constants->rise_time_mean)/constants->rise_time_sd,2)   + 0.5*pow((pin.rise_time-constants->rise_time_mean)/constants->rise_time_sd,2)
                          -0.5*pow((c0_test)/constants->c0_prior_sd,2)+ 0.5*pow((pin.c0)/constants->c0_prior_sd,2) 
                          -0.5*(res_test)/pout.sigma2 + 0.5*(res/pout.sigma2);
                          //+ utils::Z_factor(pin.Amax,Amax_test,constants->Amax_prior_sd,sqrt(pout.sigma2));  // this term would be needed for joint proposals in A and sigma2

    if(gsl_rng_uniform(rng) < exp(log_alpha_MH)) {
        pout.Amax  = Amax_test;
        pout.A     = A_test;
        pout.c0    = c0_test;
        pout.rise_time  = taur_test;
        pout.decay_time  = taud_test;
        pout.gamma = gamma_test;
        pout.omega = omega_test;
        traj.C     = C_test;
    } else {
        pout.Amax  = pin.Amax;
        pout.A     = pin.A;
        pout.c0    = pin.c0;
        pout.rise_time = pin.rise_time;
        pout.decay_time = pin.decay_time;
        pout.gamma = pin.gamma;
        pout.omega = pin.omega;
    }


}

void SMC::PF(const param &par){

    // define particle system
    cout<<"TIME = "<<TIME<<endl;
    cout<<"nparticles = "<<nparticles<<endl;

    // define weights
    double w[nparticles],logW[nparticles];

    // initialize all particles at time 0
    for(unsigned int i=0;i<nparticles;++i){
        rmu(particleSystem[0][i],data_y(0),par);
    }

    for(unsigned int t=1;t<TIME;++t){
        cout<<"time = "<<t<<"\r";
        // Multinomial resampling
        // // // collect weights
        for(unsigned int i=0;i<nparticles;i++){
            logW[i] = particleSystem[t-1][i].logWeight;
        }

        utils::w_from_logW(logW,w,nparticles);

        gsl_ran_discrete_t *rdisc = gsl_ran_discrete_preproc(nparticles, w);
        
        for(unsigned int i=0;i<nparticles;i++){
            
            //set ancestor of particle i
            unsigned int a = gsl_ran_discrete(rng,rdisc);
            particleSystem[t][i].ancestor = a; 
            particleSystem[t][i].logWeight = -log(nparticles);

            //move particle
            move_and_weight(particleSystem[t][i], particleSystem[t-1][a], data_y(t), par);
        }

        gsl_ran_discrete_free(rdisc);

    }

    ofstream outfile_C("PF_output_C.dat");
    ofstream outfile_B("PF_output_B.dat");
    ofstream outfile_A("PF_output_ancestors.dat");
    ofstream outfile_test("PF_all_part.dat");
    ofstream outfile_S("PF_output_S.dat");


    for(unsigned int i=0;i<nparticles;i++){
        unsigned int a=i;
        unsigned int t=TIME;

        while(t>0){
            t--;
            outfile_A<<a<<' ';
            outfile_C<<particleSystem[t][a].C(0)<<' ';
            outfile_S<<particleSystem[t][a].S<<' ';
            outfile_B<<particleSystem[t][a].B<<' ';
            a=particleSystem[t][a].ancestor;
        }
        outfile_C<<endl;
        outfile_B<<endl;
        outfile_A<<endl;
        outfile_S<<endl;
    }

  
}
 
void SMC::PGAS(const param &par, const Trajectory &traj_in, Trajectory &traj_out){

    unsigned int t,a;
    int i;

    // set particle 0 from input trajectory
    for(t=0;t<TIME;t++){
        particleSystem[t][0].B=traj_in.B(t);
        particleSystem[t][0].burst=traj_in.burst(t);
        if(t==0) {
            particleSystem[t][0].C={traj_in.C(t),0};
        } else {
            particleSystem[t][0].C={traj_in.C(t),traj_in.C(t-1)};
        }   
        particleSystem[t][0].S=traj_in.S(t);
        particleSystem[t][0].ancestor=0;
    }

    // set spike count to all particles if training
    if(constants->KNOWN_SPIKES){
        for(t=0;t<TIME;t++){
            for(i=0;i<nparticles;i++){
                particleSystem[t][i].S=traj_in.S(t);
            }
        }
    }

    // define weights
    double w[nparticles],logW[nparticles];
    double ar_w[nparticles], ar_logW[nparticles];

    // initialize all particles at time 0

    for(i=0;i<nparticles;++i){
        rmu(particleSystem[0][i],data_y(0),par,i==0);
    }

    for(t=1;t<TIME;++t){
        // Retrieve weights and calculate ancestor resampling weights for particle 0
        for(i=0;i<nparticles;i++){
            logW[i] = particleSystem[t-1][i].logWeight;
            ar_logW[i] = logW[i]+logf(particleSystem[t-1][i],particleSystem[t][0],par);
        }

        utils::w_from_logW(logW,w,nparticles);
        utils::w_from_logW(ar_logW,ar_w,nparticles);

        gsl_ran_discrete_t *rdisc = gsl_ran_discrete_preproc(nparticles, w);
        gsl_ran_discrete_t *ar_rdisc = gsl_ran_discrete_preproc(nparticles, ar_w);

        // Ancestor resampling of particle 0
        a = gsl_ran_discrete(rng,ar_rdisc);
        particleSystem[t][0].ancestor = a;
        
        // Resampling particles 1:nparticles
        for(i=1;i<nparticles;i++){
            a = gsl_ran_discrete(rng,rdisc);
            particleSystem[t][i].ancestor = a;
        }

        // Move and weight particles
        for(i=0;i<nparticles;i++){
            a = particleSystem[t][i].ancestor;
            if(constants->KNOWN_SPIKES){
                move_and_weight_GTS(particleSystem[t][i], particleSystem[t-1][a], data_y(t), par, i==0);
            } else {
                move_and_weight(particleSystem[t][i], particleSystem[t-1][a], data_y(t), par, i==0);
            }
        }

        gsl_ran_discrete_free(rdisc);
        gsl_ran_discrete_free(ar_rdisc);
        
    }

    // Now use the last particle set to resample the new trajectory
    for(i=0;i<nparticles;i++){
        logW[i] = particleSystem[TIME-1][i].logWeight;
    }

    utils::w_from_logW(logW,w,nparticles);

    gsl_ran_discrete_t *rdisc = gsl_ran_discrete_preproc(nparticles, w);
    i = gsl_ran_discrete(rng,rdisc);

    t=TIME;

    while(t>0){
        t--;
        traj_out.B(t)=particleSystem[t][i].B;
        traj_out.burst(t)=particleSystem[t][i].burst;
        traj_out.C(t)=particleSystem[t][i].C(0);
        traj_out.S(t)=particleSystem[t][i].S;
        traj_out.Y(t)=data_y(t);
        i=particleSystem[t][i].ancestor;
    }
 
}
 
