files=cart_mpi.f90 initialization.f90 lbm.f90 main.f90


run:$(files)
	mpiifort -qopenmp $(files) -o run -O3
	rm -f *.mod
	rm -f *~
