# ['gfortran','ifort']
COMPILER := gfortran

#------------------------------------------------------------
# external fortran libraries
#------------------------------------------------------------ 
CONDA=$(HOME)/anaconda3
JSON=-I$(CONDA)/include/json-fortran -L$(CONDA)/lib/json-fortran -ljsonfortran
FLAP=-I$(CONDA)/include/flap -L$(CONDA)/lib -lflap
FORTRESS=-I$(CONDA)/include/fortress -L$(CONDA)/lib/ -lfortress


#------------------------------------------------------------
# compilers, linkers, and flags
#------------------------------------------------------------
ifeq ($(COMPILER),gfortran)
FC = mpif90 
LINK = mpif90 
FC2 = -O3 -ffree-line-length-1000
LAPACK = -lopenblas
else ifeq ($(COMPILER),ifort)
FC = mpifort
LINK = mpifort
# for debug FC2 = -g -check all -warn all
# use -fp-model precise for value-safe optimization of fp calculations (SLOW)
FC2 = -O3 -nocheck -inline-level=2 -shared-intel -mcmodel=medium -xSSE4.2 -ipo
LAPACK = -mkl
endif 

VPATH=src/model:src/temp:src/fortress:\
      src/drivers:src/linear_model:\
      generated_linear_model/base:generated_linear_model/base/amiller


LOBJS = rng_serial.o utils.o linear_solution.o polydef.o \
	model_details.o get_decisionrule.o get_decisionrule_parallel.o \
	simulate_model.o class_model.o  \
	pdf_fcns.o pdfs.o priorfcn.o inbounds.o




%.o: %.f90
	$(FC) $(FORTRESS) $(FC2) -c $< #-mkl -c  $<

driver_irfs: $(LOBJS) driver_irfs.f90 
	$(FC) $(FC2) $^ -o driver_irfs $(LAPACK) $(FLAP) $(FORTRESS) $(JSON)

driver_simdata: $(LOBJS) driver_simdata.f90
	$(FC) $(FC2) $^ -o driver_simdata  $(LAPACK) $(FLAP) $(FORTRESS) $(JSON)

driver_selectmoments: $(LOBJS) driver_selectmoments.f90
	$(FC) $(FC2)  $^ -o driver_selectmoments $(LAPACK) $(FLAP) $(FORTRESS) $(JSON)

driver_prwmh: $(LOBJS) RandomNumber.o particles.o ParallelParticleFilter.o driver_prwmh.f90
	$(FC) $(FC2) $^ -o driver_prwmh $(LAPACK) $(FLAP) $(FORTRESS) $(JSON)

driver_smoother: $(LOBJS) RandomNumber.o particles.o class_ParticleSmoother.o driver_smoother.f90
	$(FC) $(FC2) $^ -o driver_smoother $(LAPACK) $(FLAP) $(FORTRESS) $(JSON)

driver_equity_premium : $(LOBJS) driver_equity_premium.f90
	$(FC) $(FC2) $^ -o driver_equity_premium $(LAPACK) $(FLAP) $(FORTRESS) $(JSON)

driver_altsim: $(LOBJS) driver_altsim.f90
	$(FC) $(FC2) $^ -o driver_altsim $(LAPACK) $(FLAP) $(FORTRESS) $(JSON)

driver_nomr: driver_nomr.f90 rng_serial.o linear_solution.o polydef.o model_details_laggedactual.o  utils.o get_decisionrule.o get_decisionrule_parallel.o simulate_model.o class_model.o  $(MODNAME)_AMA_matrices.o RandomNumber.o pdf_fcns.o pdfs.o priorfcn.o inbounds.o 
	$(FC) $(FC2)  $^ -o driver_nomr $(LAPACK) $(FLAP) $(FORTRESS) $(JSON)

clean:	
	rm -f $(LOBJS) 
	rm -f *.o
	rm -f *.mod
	rm -f *~

