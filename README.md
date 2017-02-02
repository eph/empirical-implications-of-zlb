Companion Code for "The Empirical Implications for the Interest-Rate Lower Bound"
---------------------------------------------------------------------------------
by Christopher Gust,  Matthew Smith, and Ed Herbst

contact: Ed Herbst [edward.p.herbst@frb.gov]

up to date version of this code: http://github.com/eph/empirical-implications-of-zlb

Software Requirements
---------------------
This code was written for linux.  It may be possible to run it on OSX or Windows, I have no idea...

1. Intel Fortran (with MKL libraries)
2. MPICH2 (i.e., a message passing interface with support for `ifort`.)  Infinibad is *highly* recommended.
3. Python (+ associated packages, for figures, tables, and estimation of the linear model.)
4. Matlab (for plotting + analying impulse responses and simulated moments)

Additional Fortran Libraries
----------------------------
1. SparseAMA (https://github.com/es335mathwiz/sparseAMA)
   Sparse matrix version of the Anderson Moore Algorithm.

2. json-fortran (https://github.com/jacobwilliams/json-fortran)
   A Fortran 2008 JSON API.

3. FLAP (https://github.com/szaghi/FLAP)
   Fortran command line arguments parser.

Installation + Running
----------------------
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


    Finally, remember when actually running the program to invoke it using MPI.
    ```
    mpirun -n NPROCS ./PROGRAM --some arguments 
    ```
