files=cart_mpi_3d.f90 initialization_3d.f90 lbm_3d.f90 main.f90
FC_FLAGS=-g -O3
OFFLOAD_FLAGS=#-mp=gpu -gpu=cc80
#FC=/soft/thetagpu/hpc-sdk/Linux_x86_64/21.3/compilers/bin/nvfortran
MPIFC=mpif90
run:$(files)	
	$(MPIFC) $(files) -o run ${FC_FLAGS} ${OFFLOAD_FLAGS}
	rm -f *.mod
	rm -f *.o
	rm -f *~ 
