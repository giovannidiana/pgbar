# PGBAR: Accurate spike inference using sequential Monte Carlo on bursting generative model

This repository contains:

- The PGBAR software to infer spike times from fluorescence time series data. The method is described in our [preprint](https://doi.org/10.1101/2022.04.05.487201). We use a second order autoregressive model to describe calcium-dependent fluorescence with the addition of a time-dependent baseline to correct for low frequency modulation. Bayesian inference of spike times and autoregressive model parameters is done using the particle Gibbs method, which employs sequential Monte Carlo to estimate time-dependent variables of the model. PGBAR generates samples from the full posterior distribution of spike times and model parameters. 
- Linescan data from cerebellar granule cells (soma and bouton recordings)
- Visualization tools running on the browser

## Dependencies
* [GNU Scientific Libraries](https://www.gnu.org/software/gsl/)
* [armadillo](http://arma.sourceforge.net/) (version >11)
* [jsoncpp](https://github.com/open-source-parsers/jsoncpp)

## Install
The main software for data analysis is written in `C++`. To compile it you need a C++ compiler in your system (e.g. g++ or clang).
Once the necessary dependencies are present, the software can then be installed by simply run `make` from the command line.


## Getting started

### Settings 
In order to run the inference, the user must provide a JSON file specifying the settings required by the algorithm such as the prior distributions of model parameters such as peak response, rise and decay time constants. The file [constants/constants_template.json](constants/constants_template.json) provides a template:

```
{
    "data":{
        "sampling_frequency":20
    },

    "AR":{
        "p":2
    },

    "BM":{
        "bm_sigma":1e-3
    },

    "priors":{
        "rise time":50,
        "rise_time_sd":5,
        "decay time":1500.0,
        "decay_time_sd":2,
        "Amax_prior_mean": 1.0,
        "Amax_prior_sd": 0.2,
        "alpha rate 0": 1,
        "beta rate 0": 20,
        "alpha rate 1":10,
        "beta rate 1": 10,
        "c0_prior_sd": 1,
        "alpha w01": 10,
        "beta w01": 100,
        "alpha w10": 10,
        "beta w10": 100,
        "alpha sigma2": 2,
        "beta sigma2": 2
    },

    "proposals":{
        "c0_proposal_sd":0.05,
        "Amax_proposal_sd":0.05,
        "decay_time_proposal_sd":10,
        "rise_time_proposal_sd":1
    },

    "MCMC":{
        "seed":2,
        "niter":200,
        "nparticles":500,
        "SAMPLE_KINETICS":1,
        "SAMPLE_SPIKES":1,
        "MOVE_SPIKES":1,
        "SAMPLE_PARAMETERS":1
    }
}
```

### Analysis of $\Delta F/F$

The command to run the main analysis software is
```

```

