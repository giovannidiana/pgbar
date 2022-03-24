#ifndef CONSTANTS_H
#define CONSTANTS_H
#include <string>
using namespace std;

class constpar {
public:
    constpar();
    constpar(string);

    void print();

    double bm_sigma = 1e-4;

    double alpha_sigma2;
    double beta_sigma2;      // Note that beta in GSL is actually the SCALE, not the rate.
    double alpha_rate_b0, beta_rate_b0;
    double alpha_rate_b1, beta_rate_b1;
    double alpha_w01, beta_w01;
    double alpha_w10, beta_w10;
    double Amax_prior_mean;
    double Amax_prior_sd;
    double c0_prior_sd;

    double sampling_frequency;
    double decay_time_mean_ms, rise_time_mean_ms;
    double decay_time_mean;
    double decay_time_sd;
    double rise_time_mean;
    double rise_time_sd;
    double decay_time_proposal_sd;
    double rise_time_proposal_sd;
    double Amax_proposal_sd;
    double c0_proposal_sd;

    int AR;
    bool MOVE_SPIKES         = true;
    bool SAMPLE_KINETICS     = true;
    bool SAMPLE_SPIKES       = true;
    bool SAMPLE_PARAMETERS   = true;
    bool CROP_TRACE          = false;
    bool KNOWN_SPIKES        = false;

    int seed;
    int niter;
    int nparticles;

    // Preprocessing
    int t_min, t_max;
    int baseline_frames;
    int nospike_before;
    int normalization_index;

    string output_folder;
    string tag;

    void set_time_scales();
};
#endif
