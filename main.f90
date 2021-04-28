!Geng Liu, April 2nd, 2021
program main
  !The main program is where LBM is implemented.
  use lbm

  implicit none
  
  !Initialize MPI and get MPI rank and # of processors.
  call MPI_INIT (ierr)
  call MPI_COMM_RANK (MPI_COMM_WORLD, rank, ierr)
  call MPI_COMM_SIZE (MPI_COMM_WORLD, num_procs, ierr)

  !Read data for initialization
  call DataInput
  !Set up cartesian MPI
  call MPISetup
  
  !Define lattice model
  call InitLattice
  !Allocate arrays
  call AllocateArrays
  !Define geometry
  call SetGeometry
  
  !OpenMP GPU parallelization: global offloading
  !$OMP TARGET DATA map(to:t,e,f,fb,p,u,geo,bl,br,bu,bd,b_user,fluid_id) device(0)
  !Device number to be changed for cross node implementations.
  
  !Initialize u,p
  call InitUP

  !Initialization of PDFs
  call InitPDF

  !LBM loop
  do iter=0,max_step

     !print out maximum velocity magnitude every "interv" steps
     if (mod(iter,interv).eq.0)then
        if(rank.eq.0)then         
          write(*,*)"T=",iter/t_intv
        endif
        call Monitor
     endif

     call Collision    
     call BoundaryCondition
     call Propagation    
     call PostProcessing     

  enddo
 
  !Write results to file
  call Write

  !$OMP END TARGET DATA

  !Deallocate arrays
  call DeAllocateArrays

  !MPI end here
  call MPI_FINALIZE (ierr)

endprogram main
