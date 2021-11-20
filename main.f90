!Geng Liu, April 2nd, 2021
program main 
  !The main program is where LBM is implemented.
  use lbm

  implicit none
  integer gpuid
  
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
  gpuid=mod(rank,8)+1!# of GPU per node: 8
  !$OMP TARGET ENTER DATA map(to:t,e,local_length,local_start,f,fb,p,u,geo,bl,br,bu,bd,b_user,fluid_id) device(gpuid)
  !Device number to be changed for cross node implementations.
  
  !Initialize u,p
  call InitUP

  !Initialization of PDFs
  call InitPDF

  count=0!Output count
  
  !LBM loop 
  do iter=0,max_step

     !print out maximum velocity magnitude every "interv" steps
     if (mod(iter,interv).eq.0)then        
        call Monitor
!        call WriteBinary
        count=count+1
     endif

     call Collision
     
     call BoundaryCondition
     call Propagation    
     call PostProcessing     

  enddo
 
  !Write results to file
  call WriteBinary

  !$OMP TARGET EXIT DATA map(delete:t,e,local_length,local_start,f,fb,p,u,geo,bl,br,bu,bd,b_user,fluid_id)

  !Deallocate arrays
  call DeAllocateArrays

  !MPI end here
  call MPI_FINALIZE (ierr)

endprogram main
