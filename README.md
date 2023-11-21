# BARSMC: Accurate spike inference using sequential Monte Carlo on bursting generative model

This repository contains:

- The PGBAR software to infer spike times from fluorescence time series data. The method is described in our [preprint](https://doi.org/10.1101/2022.04.05.487201). We use a second order autoregressive model to describe calcium-dependent fluorescence with the addition of a time-dependent baseline to correct for low frequency modulation. Bayesian inference of spike times and autoregressive model parameters is done using the particle Gibbs method, which employs sequential Monte Carlo to estimate time-dependent variables of the model. PGBAR generates samples from the full posterior distribution of spike times and model parameters. 
- Linescan data from cerebellar granule cells (soma and bouton recordings)
- Visualization tools running on the browser

## Dependencies
* [GNU Scientific Libraries](https://www.gnu.org/software/gsl/)
* [armadillo](http://arma.sourceforge.net/) (version >11)
* [jsoncpp](https://github.com/open-source-parsers/jsoncpp)



