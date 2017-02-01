### makefile --- 
## 
## Description: 
## 
## Author: Ed Herbst [edward.p.herbst@frb.gov]
## Last-Updated: 02/01/17
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
# sparseAMA
SPAMADIR=/msu/home/m1cjg01/research/aer_revision_code/Fortran/sparseAMAViaModelEZ/sparseAMA
FOBJS = $(patsubst %.f,%.o,$(wildcard $(SPAMADIR)/src/main/fortran/*.f))
COBJS = $(patsubst %.c,%.o,$(wildcard $(SPAMADIR)/src/main/c/*.c))
# AIM FILE
MODNAME=modelv5_2


# json-fortran
JSON =  -I/mq/home/m1eph00/tmp/json-fortran/lib -L/mq/home/m1eph00/lib/ -ljsonfortran

# SLICOT
SLICOT = -L/mq/home/m1eph00/lib -lslicot_sequential




LOBJS = rng_serial.o utils.o linear_solution.o polydef.o model_details.o get_decisionrule.o get_decisionrule_parallel.o  class_model.o driver.o $(MODNAME)_AMA_matrices.o
LOBJS2 = rng_serial.o utils.o linear_solution.o polydef.o model_details.o get_decisionrule.o get_decisionrule_parallel.o simulate_model.o class_model.o driver2.o $(MODNAME)_AMA_matrices.o


#------------------------------------------------------------
# compiler flags
#------------------------------------------------------------
# for debug FC2 = -g -check all -warn all
# use -fp-model precise for value-safe optimization of fp calculations (SLOW)
FC2 = -openmp -O3 -nocheck -inline-level=2 -shared-intel -mcmodel=medium -xSSE4.2 -ipo 
FFLAGS = -I$(SPAMADIR)/src/main/include -I/opt/intel/mkl/include
CFLAGS = -c -I$(SPAMADIR)/src/main/include -I/opt/intel/mkl/include
VPATH=src/model:src/temp:src/fortress:src/drivers:src/linear_model:generated_linear_model/base:generated_linear_model/base/amiller



%.o: %.f90
	$(FC) $(FC2) -mkl -c  $<

$(MODNAME)_AMA_matrices.o : $(MODNAME)_AMA_matrices.c
	$(CC)  $(CFLAGS) -c $< -I$(SPAMADIR)/src/main/include -shared -fPIC -ipo 

# DRIVERS
rundriver: $(LOBJS) $(FOBJS) $(COBJS)
	$(LINK) $(FC2) -mkl $(LOBJS) $(FOBJS) $(COBJS)  -lm $(LAPACKLIBS) -o rundriver

rundriver_irfs: driver_irfs.f90 rng_serial.o linear_solution.o polydef.o model_details.o utils.o get_decisionrule.o get_decisionrule_parallel.o simulate_model.o pdf_fcns.o pdfs.o priorfcn.o inbounds.o class_model.o $(MODNAME)_AMA_matrices.o 
	$(FC) $(FC2) $(FOBJS) $(COBJS) -mkl $^ -o rundriver_irfs

rundriver_zlbstats: driver_zlbstats.f90 rng_serial.o linear_solution.o polydef.o model_details.o utils.o get_decisionrule.o get_decisionrule_parallel.o simulate_model.o pdf_fcns.o pdfs.o priorfcn.o inbounds.o class_model.o $(MODNAME)_AMA_matrices.o 
	$(FC) $(FC2) $(FOBJS) $(COBJS) -mkl $^ -o rundriver_zlbstats

rundriver_simdata: driver_simdata.f90 rng_serial.o linear_solution.o polydef.o model_details.o utils.o get_decisionrule.o get_decisionrule_parallel.o simulate_model.o pdf_fcns.o pdfs.o priorfcn.o inbounds.o class_model.o $(MODNAME)_AMA_matrices.o 
	$(FC) $(FC2) $(FOBJS) $(COBJS) -mkl $^ -o rundriver_simdata

driver_altsim: driver_altsim.f90 rng_serial.o linear_solution.o polydef.o model_details.o  utils.o get_decisionrule.o get_decisionrule_parallel.o class_model.o  $(MODNAME)_AMA_matrices.o RandomNumber.o pdf_fcns.o pdfs.o priorfcn.o inbounds.o simulate_model.o
	$(FC) $(FC2) $(FOBJS) $(COBJS) -mkl $^ -o driver_altsim $(JSON)

driver_nomr: driver_nomr.f90 rng_serial.o linear_solution.o polydef.o model_details_laggedactual.o  utils.o get_decisionrule.o get_decisionrule_parallel.o simulate_model.o class_model.o  $(MODNAME)_AMA_matrices.o RandomNumber.o pdf_fcns.o pdfs.o priorfcn.o inbounds.o 
	$(FC) $(FC2) $(FOBJS) $(COBJS) -mkl $^ -o driver_nomr $(JSON)


.PHONY : rwmh_driver 

rwmh_driver: prwmh_driver.f90 rng_serial.o linear_solution.o polydef.o model_details.o  utils.o get_decisionrule.o get_decisionrule_parallel.o class_model.o  simulate_model.o $(MODNAME)_AMA_matrices.o RandomNumber.o particles.o ParallelParticleFilter.o pdf_fcns.o pdfs.o priorfcn.o inbounds.o TemperedParticleFilter.o
	$(FC) $(FC2) $(FOBJS) $(COBJS) -mkl $^ -o rwmh_driver $(JSON)

driver_rhoeta: driver_rhoeta.f90 rng_serial.o linear_solution.o polydef.o model_details.o  utils.o get_decisionrule.o get_decisionrule_parallel.o simulate_model.o class_model.o  $(MODNAME)_AMA_matrices.o RandomNumber.o particles.o ParallelParticleFilter.o pdf_fcns.o pdfs.o priorfcn.o inbounds.o TemperedParticleFilter.o
	$(FC) $(FC2) $(FOBJS) $(COBJS) -mkl $^ -o driver_rhoeta $(JSON)

rwmh_driver_linear_model: rwmh_driver_linear_model.f90 filter.o gensys.o as63.o prior.o model_linear.o RandomNumber.o particles.o ParallelParticleFilter.o pdf_fcns.o 
	$(FC) $(FC2) -mkl $^ -o rwmh_driver_linear_model $(SLICOT) $(JSON)

test_linear_model: test_linear.f90 gensys.o filter.o model_linear.o RandomNumber.o particles.o TemperedParticleFilter.o ParallelParticleFilter.o
	$(FC) $(FC2) -openmp -mkl  $^ -o test_linear_model $(SLICOT) $(JSON)

 smoother_driver: smoother_driver.f90 rng_serial.o linear_solution.o polydef.o model_details.o  utils.o get_decisionrule.o get_decisionrule_parallel.o class_model.o  $(MODNAME)_AMA_matrices.o RandomNumber.o particles.o class_ParticleSmoother.o pdf_fcns.o pdfs.o priorfcn.o inbounds.o 
	$(FC) $(FC2) $(FOBJS) $(COBJS) -openmp -mkl  $^ -o driver_smoother $(SLICOT) $(JSON)


smoother_driver_linear: smoother_driver.f90 gensys.o filter.o as63.o prior.o model_linear.o RandomNumber.o particles.o class_ParticleSmoother.o 
	$(FC) $(FC2) -openmp -mkl  $^ -o driver_smoother_linear $(SLICOT) $(JSON)


.PHONY : linear_model 

linear_model:
	rm -rf generated_linear_model
	python python/create_linear_model.py
	cd generated_linear_model && make smc_driver_mpi

.PHONY: estimate_linear_model
estimate_linear_model: linear_model
	echo "cd generated_linear_model && tmpiexec -n 12 ./smc_model_v5_2_linear" > ~/tmp/run.txt
	qsub -l procs=12 ~/tmp/run.txt


clean:	
	rm -f $(LOBJS) 
	rm -f *.o
	rm -f *.mod
	rm -f *~

