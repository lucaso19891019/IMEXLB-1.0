# IMEXLB-1

## Purpose
The IMEXLB project aims to develop a Lattice Boltzmann Method (LBM) proxy application code-suite for heterogeneous platforms (such as ThetaGPU). A ProxyApp, by definition, is a proxy for a full-fledged application code that simulates a wider array of problems. The IMEXLB ProxyApp, therefore, is a small self-contained code unit, with minimal dependencies, that will demonstrate both the science as well as the route to achieving parallel performance on a heterogeneous computing platform for one or more exemplar problems.

For this project, in particular, the IMEXLB ProxyApp suite will consist of Fortran and C++ codes that will simulate two-phase canonical flows using the MPI+X parallelization paradigm, where X denotes the SIMD parallelization paradigms: OpenMP-4.5+/Kokkos/DPC++/RAJA. This ProxyApp will serve as a reference for other LB code developers who may wish to either use IMEXLB as a building block, or follow the SIMD parallelization strategy to modify their own code, to do their simulations on platforms like Polaris, Aurora, & Frontier.

## Characteristics
This code package simulates the fuild with lattice Boltzmann method. It can solve 2D and 3D problem. The code is accelerated by MPI and OpenMP GPU offloading directives.
The example problem is cylinder flow (2D) and sphere flow (3D).
Use nvhpc-mpi compiler.
For 3D applications, add "_3d" to the file names in Makefile.
Run the code with command: mpirun -np (# of processors) run
