# IMEXLB-1

## Purpose
The IMEXLB project aims to develop a Lattice Boltzmann Method (LBM) proxy application code-suite for heterogeneous platforms (such as ThetaGPU). A ProxyApp, by definition, is a proxy for a full-fledged application code that simulates a wider array of problems. The IMEXLB ProxyApp, therefore, is a small self-contained code unit, with minimal dependencies, that will demonstrate both the science as well as the route to achieving parallel performance on a heterogeneous computing platform for one or more exemplar problems.

For this project, in particular, the IMEXLB ProxyApp suite will consist of Fortran and C++ codes that will simulate two-phase canonical flows using the MPI+X parallelization paradigm, where X denotes the SIMD parallelization paradigms: OpenMP-4.5+/Kokkos/DPC++/RAJA. This ProxyApp will serve as a reference for other LB code developers who may wish to either use IMEXLB as a building block, or follow the SIMD parallelization strategy to modify their own code, to do their simulations on platforms like Polaris, Aurora, & Frontier.

## Background
LBM is a relatively novel approach to solve the Navier-Stokes equations (NSE) in the low-Mach number regime. The governing equation can be derived from the Boltzmann equation after discretizing the phase space with constant microscopic lattice velocities. One major drive behind the use of LBM in the CFD community is the ease of parallelization, but the increasing popularity of LBM can also be attributed to its unique feature: LBM solves a simplified Boltzmann equation that is essentially a set of 1st-order hyperbolic PDEs with constant microscopic lattice velocities, for which a plethora of simple yet excellent discretization schemes are available. Furthermore, all the complex non-linear effects are modeled locally at the grid points. Thus, the overall numerical scheme becomes extremely simple and efficient. LBM can deliver better solution quality than a comparably discretized NSE solver at substantially lower computational cost. Furthermore, the linear advection part of LBM can be solved exactly due to unity CFL property. The absence of dispersive and dissipative errors in LBM makes it an ideal choice for the simulation of turbulent flows, interfacial flow, and acoustics, in which conservation of kinetic energy and isotropy play a crucial role. The highly scalable IMEXLB ProxyApp suite is a first step towards the development of other full-fledged LB codes that will complement NSE codes that are being developed at Argonne (e.g. Nek5000, nekRS, PHASTA), and is a much needed building block for simulating multiphase flows on exascale computing platforms at Argonne and the other DOE laboratories.

## Code Characteristics
This code package simulates the fuild with lattice Boltzmann method. It can solve 2D and 3D problem. The code is accelerated by MPI and OpenMP GPU offloading directives. The example problem is cylinder flow (2D) and sphere flow (3D). Use nvhpc-mpi compiler.
For 3D applications, add "_3d" to the file names in Makefile. Run the code with command: mpirun -np (# of processors) run.


## Citing IMEXLB-1
>David F. Richards, et al., “,FY21 Proxy App Suite Release: Report for ECP Proxy App Project Milestone ADCD504-12” LLNL-TR-827482, Sept. 2021.https://proxyapps.exascaleproject.org/wp-content/uploads/2021/10/proxyReport21.pdf

## Development Team
Points of contact for further assistance - Geng Liu (gliu@anl.gov), Taehun Lee (thlee@ccny.cuny.edu), Saumil Patel (spatel@anl.gov), Ramesh Balakrishnan (bramesh@anl.gov).

