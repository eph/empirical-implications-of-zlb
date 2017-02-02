### makefile --- 
## 
## Description: 
## 
## Author: Ed Herbst [edward.p.herbst@frb.gov]
## Last-Updated: 02/02/17
## 

#------------------------------------------------------------
# compilers and linkers
#------------------------------------------------------------
CC = icc
FC = mpif90 -f90=ifort 
LINK = mpif90 -f90=ifort


#------------------------------------------------------------
# external fortran libraries
#------------------------------------------------------------ 
CONDA=$(HOME)/miniconda2/envs/ifort
# sparseAMA
SPAMADIR=/msu/home/m1cjg01/research/aer_revision_code/Fortran/sparseAMAViaModelEZ/sparseAMA
FOBJS = $(patsubst %.f,%.o,$(wildcard $(SPAMADIR)/src/main/fortran/*.f))
COBJS = $(patsubst %.c,%.o,$(wildcard $(SPAMADIR)/src/main/c/*.c))
# AIM FILE
MODNAME=modelv5_2



# json-fortran
JSON=-I$(CONDA)/include/json-fortran -L$(CONDA)/lib/json-fortran -ljsonfortran

# FLAP

FLAP=-I$(CONDA)/include/flap -L$(CONDA)/lib -lflap

# SLICOT
SLICOT = -L/mq/home/m1eph00/lib -lslicot_sequential


#------------------------------------------------------------
# compiler flags
#------------------------------------------------------------
# for debug FC2 = -g -check all -warn all
# use -fp-model precise for value-safe optimization of fp calculations (SLOW)
FC2 = -O3 -nocheck -inline-level=2 -shared-intel -mcmodel=medium -xSSE4.2 -ipo 
FFLAGS = -I$(SPAMADIR)/src/main/include -I/opt/intel/mkl/include
CFLAGS = -c -I$(SPAMADIR)/src/main/include -I/opt/intel/mkl/include
VPATH=src/model:src/temp:src/fortress:src/drivers:src/linear_model:generated_linear_model/base:generated_linear_model/base/amiller


LOBJS = rng_serial.o utils.o linear_solution.o polydef.o \
	model_details.o get_decisionrule.o get_decisionrule_parallel.o \
	simulate_model.o class_model.o $(MODNAME)_AMA_matrices.o \
	pdf_fcns.o pdfs.o priorfcn.o inbounds.o


%.o: %.f90
	$(FC) $(FC2) -mkl -c  $<

$(MODNAME)_AMA_matrices.o : $(MODNAME)_AMA_matrices.c
	$(CC) -shared -fPIC -ipo  $(CFLAGS) -c $< -I$(SPAMADIR)/src/main/include 

driver_irfs: $(LOBJS) driver_irfs.f90 
	$(FC) $(FC2) $(FOBJS) $(COBJS) $(FLAP) -mkl $^ -o driver_irfs

driver_simdata: $(LOBJS) driver_simdata.f90
	$(FC) $(FC2) $(FOBJS) $(COBJS) $(FLAP) -mkl $^ -o driver_simdata

driver_selectmoments: $(LOBJS) driver_selectmoments.f90
	$(FC) $(FC2) $(FOBJS) $(COBJS) $(FLAP) -mkl $^ -o driver_selectmoments

driver_prwmh: $(LOBJS) RandomNumber.o particles.o ParallelParticleFilter.o driver_prwmh.f90
	$(FC) $(FC2) $(FOBJS) $(COBJS) $(FLAP)  -mkl $^ -o driver_prwmh $(JSON)

driver_smoother: $(LOBJS) RandomNumber.o particles.o class_ParticleSmoother.o driver_smoother.f90
	$(FC) $(FC2) $(FOBJS) $(COBJS) $(FLAP)  -mkl $^ -o driver_smoother $(JSON)

driver_equity_premium : $(LOBJS) driver_equity_premium.f90
	$(FC) $(FC2) $(FOBJS) $(COBJS) $(FLAP) -mkl $^ -o driver_equity_premium $(JSON)

driver_altsim: $(LOBJS) driver_altsim.f90
	$(FC) $(FC2) $(FOBJS) $(COBJS) $(FLAP) -mkl $^ -o driver_altsim $(JSON)

driver_nomr: driver_nomr.f90 rng_serial.o linear_solution.o polydef.o model_details_laggedactual.o  utils.o get_decisionrule.o get_decisionrule_parallel.o simulate_model.o class_model.o  $(MODNAME)_AMA_matrices.o RandomNumber.o pdf_fcns.o pdfs.o priorfcn.o inbounds.o 
	$(FC) $(FC2) $(FOBJS) $(COBJS) -mkl $^ -o driver_nomr $(JSON)


.PHONY : rwmh_driver 


driver_rhoeta: driver_rhoeta.f90 rng_serial.o linear_solution.o polydef.o model_details.o  utils.o get_decisionrule.o get_decisionrule_parallel.o simulate_model.o class_model.o  $(MODNAME)_AMA_matrices.o RandomNumber.o particles.o ParallelParticleFilter.o pdf_fcns.o pdfs.o priorfcn.o inbounds.o TemperedParticleFilter.o
	$(FC) $(FC2) $(FOBJS) $(COBJS) -mkl $^ -o driver_rhoeta $(JSON)

rwmh_driver_linear_model: rwmh_driver_linear_model.f90 filter.o gensys.o as63.o prior.o model_linear.o RandomNumber.o particles.o ParallelParticleFilter.o pdf_fcns.o 
	$(FC) $(FC2) -mkl $^ -o rwmh_driver_linear_model $(SLICOT) $(JSON)

test_linear_model: test_linear.f90 gensys.o filter.o model_linear.o RandomNumber.o particles.o TemperedParticleFilter.o ParallelParticleFilter.o
	$(FC) $(FC2) -openmp -mkl  $^ -o test_linear_model $(SLICOT) $(JSON)

smoother_driver_linear: smoother_driver.f90 gensys.o filter.o as63.o prior.o model_linear.o RandomNumber.o particles.o class_ParticleSmoother.o 
	$(FC) $(FC2) -openmp -mkl  $^ -o driver_smoother_linear $(SLICOT) $(JSON)

clean:	
	rm -f $(LOBJS) 
	rm -f *.o
	rm -f *.mod
	rm -f *~

