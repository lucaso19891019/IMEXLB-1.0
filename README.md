# IMEXLB-1
This code package simulates the fuild with lattice Boltzmann method. It can solve 2D and 3D problem. The code is accelerated by MPI and OpenMP GPU offloading directives.
The example problem is cylinder flow (2D) and sphere flow (3D).
Use nvhpc-mpi compiler.
For 3D applications, add "_3d" to the file names in Makefile.
Run the code with command: mpirun -np (# of processors) run
