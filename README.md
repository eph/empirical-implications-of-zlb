Companion Code for "The Empirical Implications for the Interest-Rate Lower Bound"
---------------------------------------------------------------------------------
by Christopher Gust,  Matthew Smith, and Ed Herbst

contact: Ed Herbst [edward.p.herbst@frb.gov]

Software Requirements
---------------------
This code was written for linux.  It may be possible to run it on OSX or Windows, I have no idea...

1. Intel Fortran (with MKL libraries)
2. MPICH2 (i.e., a message passing interface with support for `ifort`.)
3. Python (+ associated packages, for figures, tables, and estimation of the linear model.)

Additional Fortran Libraries
----------------------------
1. SparseAMA (https://github.com/es335mathwiz/sparseAMA)
   Sparse matrix version of the Anderson Moore Algorithm.

2. json-fortran (https://github.com/jacobwilliams/json-fortran)
   A Fortran 2008 JSON API

Installation
------------
1. Install the `Additional Fortran Libraries`
2. Edit the makefile to reflect these libraries locations (and possibly ```LD_LIBRARY_PATH```)
3. To install individual programs, type:
   * `make rundriver_irfs`
     Installs the fortran program for simulating impulse responses.  The program is executed with
     ```make
     mpirun -n NPROC ./rundrivers_irfs --nsim NSIM --shockindex SHOCK --dset DSET > irfs.out
     ```
     where ```DSET=1``` for posterior mean and 2 for the posterior median.
     ```SHOCK=1``` for the risk premium and ```SHOCK=2``` for the MEI shock.
     More details in the code.

   * `make ru



