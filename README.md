# Companion Code for "The Empirical Implications for the Interest-Rate Lower Bound"
### by Chris Gust,  Ed Herbst, David Lopez-Salido, and Matt Smith

The paper can be found here: https://www.aeaweb.org/articles?id=10.1257/aer.20121437&&from=f

contact: Ed Herbst [ed.herbst@gmail.com]

up to date version of this code: http://github.com/eph/empirical-implications-of-zlb

## Software Requirements
This code was written for Linux.  It may be possible to run it on OSX or Windows.

To run this code, you need:

1. Intel Fortran (with MKL libraries)
2. MPICH2 (i.e., a message passing interface with support for `ifort`.)  Infiniband is *highly* recommended.
3. Python (+ associated packages, for figures, tables, and estimation of the linear model.  In particular, the DSGE package is required.  Install this using anaconda: `conda install dsge -c eherbst`.)
4. Matlab (for plotting + analying impulse responses and simulated moments)

## Additional Fortran Libraries
1. SparseAMA (https://github.com/es335mathwiz/sparseAMA)
   Sparse matrix version of the Anderson Moore Algorithm.

2. json-fortran (https://github.com/jacobwilliams/json-fortran)
   A Fortran 2008 JSON API.

3. FLAP (https://github.com/szaghi/FLAP)
   Fortran command line arguments parser.

4. fortress (https://github.com/eph/fortress)
   Fortran utilities for Bayesian estimation of time series models.  Only needed to estimate the linearized model.

## Installation + Running
1. Install the `Additional Fortran Libraries`
2. Edit the makefile to reflect these libraries locations (and possibly ```LD_LIBRARY_PATH```)
3. To install individual programs, at the prompt (from this directory) type:
   ```make PROGRAM```
   where the possible programs are described below:

   | Program Name           | Description                                                                           |
   | ------------           | -----------                                                                           |
   | `driver_irfs`          | Computes the impulse reponse functions in Figures 2 and 3.                            |
   | `driver_simdata`       | Simulates data from the model.                                                        |
   | `driver_selectmoments` | Simulates data from the model, and computes ZLB statistics for histogram in Figure 6. | 
   | `driver_prwmh`         | Runs the particle-filter-based MCMC algorithm.                                        |
   | `driver_smoother`      | Runs the particle filter and smoother for a given set of parameter values.            |
   | `driver_equitypremium` | Given a set of smoothed estimates, computes the equity premium.                       |
   | `driver_altsim`        | Given a set of smoothed estimates, run a counterfactual.                              |
4. Each program takes command line arguments.  To see descriptions of the arguments type.
   ```
   ./PROGRAM --help
   ```

    For example:
	```
	eherbst@thnkpd:$ ./driver_prwmh --help
usage: driver_prwmh  [--p0 value] [--covariance value] [--prior-file value] [--scaling value] [--nsim value] [--npart value] [--zlb value] [--output-file value] [--seed value] [--help] [--version]

Program for estimating model via particle MCMC
 

Optional switches:
   --p0 value, -p0 value
    default value input/mean.txt
    starting value
   --covariance value, -H value
    default value input/cholM.txt
    Cholesky of covariance for proposal innvoations
   --scaling value, -c value
    default value 0.15
    scaling of innovation
   --nsim value, -n value
    default value 50000
    Length of MCMC chain
   --npart value
    default value 1500000
    Number of Particles for PF
   --zlb value
    default value .true.
    Impose ZLB
   --output-file value
    default value output.json
    Output file
   --seed value
    default value 1848
    Seed for RNG
   --help, -h
    Print this help message
   --version, -v
    Print version
	```

    Finally, remember when actually running the program to invoke it using MPI (except `driver_smoother`)
    ```
    mpirun -n NPROCS ./PROGRAM --some arguments 
    ```
    
Figures + Table Generation
--------------------------
## Figures
1. Comparison of the posterior
    * Estimate the linear model: `driver_linear`
    * Estimate the nonlinear model: `driver_prwmh`
    * Run `python/fig_linear_vs_nonlinear.py`

2. Response to an Exogenous Increase in the Risk Premium
    * Driver to run: `driver_irfs`
    * Figure script: `matlab -r "run('matlab/plot_nonlinear_irf.m')"`

3. Response to a Fall in Investment Efficiency
    * Driver to run: `driver_irfs -s 2`
    * Figure script: `matlab -r "run('matlab/plot_nonlinear_irf.m')"`
4. Smoothed Estimates of Model Objects
    * Driver to run: `driver_prwmh`
    * Get 1000 draws from the posterior: `python python/sample_posterior.py`
    * For each of 1000 draws run: `driver_smoother -p0 results/thinned_posterior/NNNNpara.txt --output-file outputNNNN.json`
    * Make plot `python python/fig_observable_fit.py`
    
5. Model objects for different values of Measurement Error

    * Estimate the model with a different value for the measurement errors (i.e,
   change the last values of the parameter vector). 
   *  Then follow the steps in [4]. (You may want to change the directories to not overwrite earlier
   estimations.)

6. Distribution of the Probability and Duration of Being at the ZLB

   * After estimating the model and thinning the posterior, run `./driver_zlbstats`
   * Matlab script: `matlab -r "run('matlab/generate_zlb_freq_duration.m')"`

7. The Path of the Estimated Shocks and Equitiy Premium During the Great Recession

    * Get the smoothed objects (See [4])
    * Estimate the equity premium for each set of smoothed estimates using `driver_equitypremium`
    * Run the python script `python/fig_shocks_combo.py`

8. The Path of the Estimated Technology Shocks
    * Get the smoothed objects (See [4] and [5])
    * Run the python script `python/fig_output_inf_tech.py`

9. The Contribution of the Estimated Shocks to the Great Recession

    * Get the smoothed objects (see [4]).
    * For each of set of smoothed estimates, run the 3 alt sims using `driver_altsim`
    * Run the python script `python/fig_great_recession.py`

10. The Contribution of the Zero Lower Bound to the Great Recession
    * Get the smoothed objects (see [4])
    * For each of set of smoothed estimates, run the nozlb_2008 altsim using `driver_altsim`
    * Run the python script `python/fig_effect_of_zlb_levels.py`

11. The Effect of the Lower Bound on Aggregate Spending
    * Get the smoothed objects (see [4])
    * For each of set of smoothed estimates, run the nozlb_2008 altsim using `driver_altsim`
    * Run the python script `python/fig_effect_of_zlb.py`

## Tables
1. Posterior Distribution of the Parameters
   For `N=0,1,2,3`

   ```
   mpirun -n NPROC ./driver_prwmh --p0 STARTING_VALUE.txt --seed ASEED --output-file results/mcmc/outputN.json
   ```

   Then

   ```
   python python/tab_posterior.py
   ```

2. Standard Deviations of Aggregate Variables at Business Cycle Frequencies
   Run the posterior estimation, then thin the posterior with:

   ```
   python python/thin_posterior.py
   ```
   Then run the matlab script:

   ```
   matlab -r "run('matlab/generate_momentstable.m')"
   ```
