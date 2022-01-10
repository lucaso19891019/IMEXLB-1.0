!Geng Liu, April 2nd, 2021
program main 
  !The main program is where LBM is implemented.
  use lbm

  implicit none
  !gpuid: ID of GPU device
  integer gpuid
  !start: Starting time mark.
  !Finish: Ending time mark.
  double precision start,finish
  
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
  
  !OpenMP GPU parallelization: global offloading (Remove when GPU offloading is turned off.)
  gpuid=mod(rank,omp_get_num_devices())!Obtain the GPU ID corresponding to current rank.
  call omp_set_default_device(gpuid)!Set Default GPU for current rank.
  !$OMP TARGET ENTER DATA map(to:t,e,local_length,local_start,f,fb,rho,u,ome,lap,dcrho,dccp,dmrho,dmcp,fluid_id) device(gpuid)
  !Device number to be changed for cross node implementations.
  
  !Initialize u,p
  call InitUP

  !Initialization of PDFs
  call InitPDF
  
  count=0!Output count
  
  !LBM loop 
  call cpu_time(start)
  do iter=0,max_step

     !print out maximum velocity magnitude every "interv" steps
     if (mod(iter,interv).eq.0)then        
        call Monitor
        call WriteBinary
        count=count+1
     endif

     call Collision
     
     call BoundaryCondition
     call Propagation    
     call PostProcessing     

  enddo
  call cpu_time(finish)
 
  !Print elapsed time
  if(rank.eq.0)then
     write(*,*)'Elapsed time: ',finish-start
  endif
  
  !Write results to file
  !call WriteBinary

  !$OMP TARGET EXIT DATA map(delete:t,e,local_length,local_start,f,fb,rho,u,ome,lap,dcrho,dccp,dmrho,dmcp,fluid_id)

  !Deallocate arrays
  call DeAllocateArrays

  !MPI end here
  call MPI_FINALIZE (ierr)

endprogram main
