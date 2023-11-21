#include "../include/constants.h"
#include <cmath>
#include <jsoncpp/json/json.h>
#include <fstream>
#include <iostream>

using namespace std;

constpar::constpar(string filename)
{
    Json::Reader reader;
    Json::Value  cfg;
    ifstream     paramfile(filename);

    paramfile >> cfg;

    // Integrated Brownian motion 
    bm_sigma  = cfg["BM"]["bm_sigma"].asDouble();   

    // data
    sampling_frequency = cfg["data"]["sampling_frequency"].asDouble();

    // AR constants
    AR = cfg["AR"]["p"].asInt();

    // Preprocessing
    baseline_frames     = cfg["preproc"]["baseline_frames"].asInt();
    nospike_before      = cfg["preproc"]["nospike_before"].asInt();
    normalization_index = cfg["preproc"]["normalization_index"].asInt();

    // Priors
    rise_time_mean_ms = cfg["priors"]["rise time"].asDouble();
    decay_time_mean_ms   = cfg["priors"]["decay time"].asDouble();
    alpha_rate_b0        = cfg["priors"]["alpha rate 0"].asDouble();
    beta_rate_b0         = cfg["priors"]["beta rate 0"].asDouble();
    alpha_rate_b1        = cfg["priors"]["alpha rate 1"].asDouble();
    beta_rate_b1         = cfg["priors"]["beta rate 1"].asDouble();
    Amax_prior_mean = cfg["priors"]["Amax_prior_mean"].asDouble();
    Amax_prior_sd   = cfg["priors"]["Amax_prior_sd"].asDouble();
    c0_prior_sd     = cfg["priors"]["c0_prior_sd"].asDouble();
    alpha_w01       = cfg["priors"]["alpha w01"].asDouble();
    beta_w01        = cfg["priors"]["beta w01"].asDouble();
    alpha_w10       = cfg["priors"]["alpha w10"].asDouble();
    beta_w10        = cfg["priors"]["beta w10"].asDouble();
    alpha_sigma2    = cfg["priors"]["alpha sigma2"].asDouble();
    beta_sigma2     = cfg["priors"]["beta sigma2"].asDouble();
    decay_time_sd   = cfg["priors"]["decay_time_sd"].asDouble();
    rise_time_sd    = cfg["priors"]["rise_time_sd"].asDouble();

    seed                = cfg["MCMC"]["seed"].asInt();
    niter               = cfg["MCMC"]["niter"].asInt();
    nparticles          = cfg["MCMC"]["nparticles"].asInt();
    SAMPLE_KINETICS     = cfg["MCMC"]["SAMPLE_KINETICS"].asBool();
    SAMPLE_SPIKES       = cfg["MCMC"]["SAMPLE_SPIKES"].asBool();
    MOVE_SPIKES         = cfg["MCMC"]["MOVE_SPIKES"].asBool();
    SAMPLE_PARAMETERS   = cfg["MCMC"]["SAMPLE_PARAMETERS"].asBool();
    KNOWN_SPIKES        = cfg["MCMC"]["KNOWN_SPIKES"].asBool();

    // proposals
    decay_time_proposal_sd   = cfg["proposals"]["decay_time_proposal_sd"].asDouble();
    rise_time_proposal_sd    = cfg["proposals"]["rise_time_proposal_sd"].asDouble();
    Amax_proposal_sd         = cfg["proposals"]["Amax_proposal_sd"].asDouble();
    c0_proposal_sd           = cfg["proposals"]["c0_proposal_sd"].asDouble();

}

void constpar::set_time_scales()
{
    decay_time_mean         = sampling_frequency * decay_time_mean_ms / 1000.0;
    rise_time_mean          = sampling_frequency * rise_time_mean_ms / 1000.0;
    decay_time_sd           = sampling_frequency * decay_time_sd / 1000;
    rise_time_sd            = sampling_frequency * rise_time_sd / 1000;
    decay_time_proposal_sd *= sampling_frequency/1000;
    rise_time_proposal_sd  *= sampling_frequency/1000;
}

void constpar::print()
{
    cout << "________SETTINGS_____________" << endl;
    cout << "_____________________________" << endl;
    cout << "Brownian motion" << endl;
    cout << "    sigma = " << bm_sigma << endl;
    cout << "AR" << endl;
    cout << "    p = " << AR << endl;
    cout << "    sampling frequency = " << sampling_frequency << endl;
    cout << "priors" << endl;
    cout << "    decay time = " << decay_time_mean << " (frames)" << endl;
    cout << "    rise time = " << rise_time_mean << " (frames)" << endl;
    cout << "    alpha_sigma2 = " << alpha_sigma2 << endl;
    cout << "    beta_sigma2 = " << beta_sigma2 << endl;
    cout << "    prior_mean_Amax = " << Amax_prior_mean << endl;
    cout << "    prior_sd_Amax = " << Amax_prior_sd << endl;
    cout << "MCMC" << endl;
    cout << "    niter = " << niter << endl;
    cout << "    SAMPLE_KINETICS   = " << SAMPLE_KINETICS << endl;
    cout << "    SAMPLE_PARAMETERS = " << SAMPLE_PARAMETERS << endl;
    cout << "    CROP_TRACE          = " << CROP_TRACE << endl;
    cout << "_____________________________" << endl;
}
