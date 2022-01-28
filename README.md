# IMEXLB-1.0

## Purpose
The IMEXLB project aims to develop a Lattice Boltzmann Method (LBM) proxy application code-suite for heterogeneous platforms (such as ThetaGPU). A ProxyApp, by definition, is a proxy for a full-fledged application code that simulates a wider array of problems. The IMEXLB ProxyApp, therefore, is a small self-contained code unit, with minimal dependencies, that will demonstrate both the science as well as the route to achieving parallel performance on a heterogeneous computing platform for one or more exemplar problems.

For this project, in particular, the IMEXLB ProxyApp suite will consist of Fortran and C++ codes that will simulate two-phase canonical flows using the MPI+X parallelization paradigm, where X denotes the SIMD parallelization paradigms: OpenMP-4.5+/Kokkos/DPC++/RAJA. This ProxyApp will serve as a reference for other LB code developers who may wish to either use IMEXLB as a building block, or follow the SIMD parallelization strategy to modify their own code, to do their simulations on platforms like Polaris, Aurora, & Frontier.

## Background
LBM is a relatively novel approach to solve the Navier-Stokes equations (NSE) in the low-Mach number regime. The governing equation can be derived from the Boltzmann equation after discretizing the phase space with constant microscopic lattice velocities. One major drive behind the use of LBM in the CFD community is the ease of parallelization, but the increasing popularity of LBM can also be attributed to its unique feature: LBM solves a simplified Boltzmann equation that is essentially a set of 1st-order hyperbolic PDEs with constant microscopic lattice velocities, for which a plethora of simple yet excellent discretization schemes are available. Furthermore, all the complex non-linear effects are modeled locally at the grid points.

## Code Characteristics 
* Written in Fortran 90
* MPI+OpenMP hybrid parallelism with GPU offloading directives 
* 2D (D2Q9) and 3D (D3Q27) problems 
* Example problem is flow past a circle (2D) and sphere(3D) 

## Building
The code has currently been tested on ALCF's [ThetaGPU](https://www.alcf.anl.gov/support-center/theta/theta-thetagpu-overview), as such, the following instructions assume you are building on ThetaGPU. In the future, we will provide more general instructions for non-ThetaGPU architectures. 

From the ThetaGPU login, please request a node. For example:

```
qsub -I --attrs pubnet=true -A <YOUR_PROJECT_NAME> -n 1 -q single-gpu -t 30
```

Please clone this code and use the available nvhpc-mpi compiler by doing:

```
module load nvhpc-mpi
make
```

Note: For 3D applications, add "_3d" to the file names in Makefile except for "main.f90". Run the code with command: 

```
mpirun -np (# of processors) run
```

## License
IMEXLBM is released under the BSD 3-clause license (see LICENSE file). All new contributions must be made under the BSD 3-clause license.

## Acknowledgments
This research was supported by the Exascale Computing Project (17-SC-20-SC), a joint project of the U.S. Department of Energy's Office of Science and National Nuclear Security Administration, responsible for delivering a capable exascale ecosystem, including software, applications, and hardware technology, to support the nation's exascale computing imperative.

If you want to cite IMEXLBM please use: 
>David F. Richards, et al., “FY21 Proxy App Suite Release: Report for ECP Proxy App Project Milestone ADCD504-12” LLNL-TR-827482, Sept. 2021. https://proxyapps.exascaleproject.org/wp-content/uploads/2021/10/proxyReport21.pdf
>Geng Liu, Saumil Patel, Ramesh Balakrishnan and Taehun Lee, “IMEXLBM 1.0: A Proxy Application based on the Lattice Boltzmann Method for solving Computational Fluid Dynamic problems on GPUs”, Jan. 2022. https://arXiv:2201.11330

## Development Team
Points of contact for further assistance - Geng Liu (gliu@anl.gov), Taehun Lee (thlee@ccny.cuny.edu), Saumil Patel (spatel@anl.gov), Ramesh Balakrishnan (bramesh@anl.gov).

