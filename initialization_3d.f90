!Geng Liu, April 2nd, 2021
module initialization
  !This module includes the variables and subroutines for fluid initialization.
  use cart_mpi
  use omp_lib
  implicit none

  !dy: The index increase in 1D arrays reshaped from 3D arrays per 1 y-coordinate increment.
  !dz: The index increase in 1D arrays reshaped from 3D arrays per 1 z-coordinate increment.
  integer dy,dz

  !charlength: characteristic length
  !pos1: center position of droplet 1
  !pos2: center position of droplet 2
  !beta: a constant that is related to the compressibility of bulk phases
  !radius: cylinder flow geometry related: cylinder radius, half of charlength in this case
  double precision charlength,pos1,pos2,beta,radius
  
  !rhol: liquid density
  !rhov: vapor density
  !dt: time increment
  !taul: liquid relaxation time
  !tauv: vapor relaxation time
  !ep: interface thickness
  !t_intv: dimensionless time increment
  !inter: time interval size(dimensionless) for monitor/output
  !max_t: maximum time(dimensionless)
  double precision rhol,rhov,dt,taul,tauv,ep,t_intv,inter,max_t
  !max_step: maximum number of time steps
  !iter: iteration index
  !interv: time interval size(in lattice unit) for monitor/output
  integer max_step,iter,interv
  !count: output count
  integer count

  !t: lattice weights
  !e: lattice velocities
  double precision,dimension(0:nq-1)::t
  integer,dimension(0:nq*dim-1)::e

  !f,fb: PDF and backup PDF arrays
  double precision,dimension(:),allocatable::f,fb

  !rho: density array
  !ome: relaxation frequency array
  !cp: chemical potential array
  !lap: laplacian array
  double precision,dimension(:),allocatable::rho,ome,cp,lap
  
  !Vector that store fluid bulk points' indices:
  integer,dimension(:),allocatable::fluid_id
  !Size of index vectors
  integer size_fluid

  !u: velocity array
  !dcrho: central gradient array of density
  !dmrho: mixed gradient array of density
  !dccp: central gradient array of chemical potential
  !dmcp: mixed gradient array of chemiical potential
  double precision,dimension(:),allocatable::u,dcrho,dmrho,dccp,dmcp

  
       
contains
  !DataInput: This subroutine reads data from "input_3d.in", and generates derived initial data.
  subroutine DataInput
    
    open(100,file='input_3d.in',status='old')
    read(100,*)nx,ny,nz     
    read(100,*)charlength,pos1,pos2,beta
    !charlength & pos are given in the unit of nx, need to be converted to lattice unit.
    charlength=charlength*nx
    pos1=pos1*nx
    pos2=pos2*nx
    !Read relaxation time, density, time increment and interface thickness
    read(100,*)taul,tauv,dt,rhol,rhov,ep
 
    !Read max time and output interval
    read(100,*)max_t,inter
    
    close(100)

    !Derived phyical property
 
    
    radius=charlength*0.5d0

    !Convert time to number of steps:
    t_intv=sqrt(rhol*radius**3*12.d0/((rhol-rhov)**4*beta*ep))
    interv=inter*t_intv
    max_step=max_t*t_intv

    !Display
    if(rank.eq.0)then

       write(*,'(A24,F6.2)')   "Characteristic Length :",charlength
       write(*,'(A24,F6.2)')   "Interface Thickness   :",ep
       write(*,'(A24,F6.2)')   "Surface Tension Sigma :",(rhol-rhov)**4*ep*beta/12.d0
  
       write(*,'(A24,F6.2)')   "Relaxation Time Liquid:",taul
       write(*,'(A24,F6.2)')   "Relaxation Time Vapor :",tauv

       write(*,'(A24,F4.1)')   "Density Liquid        :",rhol
       write(*,'(A24,F4.1)')   "Density Vapor         :",rhov

       write(*,'(A24,I8)')     "Maximum # of Steps    :",max_step
       write(*,'(A24,I8)')     "Interval Size         :",interv

       write(*,*)"============================="

    endif

  endsubroutine DataInput

  !---------------------------------------------
  !InitLattice: This subroutine defines the lattice velocites and weights
  !Numbering of D3Q27 lattice directions: (see order of directions in paper.)
  subroutine InitLattice
    integer i,j,k,iq   
    !Lattice velocities and weights(3D)
    do j=-1,1
       do i=-1,1
          do k=-1,1
             iq=(k+1)*9+(i+1)*3+(j+1)
             e(iq*dim)=i
             e(iq*dim+1)=j
             e(iq*dim+2)=k
             t(iq)=8.d0/4.d0**(abs(i)+abs(j)+abs(k))/27.d0
          enddo
       enddo
    enddo
  endsubroutine InitLattice

  !---------------------------------------------
  !AllocateArrays: This subroutine allocates local arrays.
  subroutine AllocateArrays
    
    allocate(f(0:array_size*nq-1))
    allocate(fb(0:array_size*nq-1))


    allocate(lap(0:array_size-1))
    allocate(cp(0:array_size-1))
    allocate(rho(0:array_size-1))
    allocate(ome(0:array_size-1))

    allocate(fluid_id(0:array_size-1))
 

    allocate(u(0:array_size*dim-1))
    allocate(dcrho(0:array_size*dim-1))
    allocate(dmrho(0:array_size*dim-1))
    allocate(dccp(0:array_size*dim-1))
    allocate(dmcp(0:array_size*dim-1))

  endsubroutine AllocateArrays
  
  !-----------------------------------------------------
  !DeAllocateArrays: This subroutine deallocates local arrays.
  subroutine DeAllocateArrays
    deallocate(f)
    deallocate(fb)

    deallocate(lap)
    
    deallocate(cp)
    deallocate(rho)
    deallocate(ome)

    deallocate(fluid_id)
    
    deallocate(u)
    deallocate(dcrho)
    deallocate(dmrho)
    deallocate(dccp)
    deallocate(dmcp)
    
    
  endsubroutine DeAllocateArrays  

  !---------------------------------------------
  !SetGeometry: This subroutine initializes the flag array according to the given geometry
  subroutine SetGeometry 
    integer i,j,k,x,y,z,id,idn
    ! dy and dz are global variables that are generated here, and should not be modified in other subroutines or functiions, dx is 1 and thus not needed.
    dy=(local_length(1)+2*ghost)
    dz=(local_length(1)+2*ghost)*(local_length(2)+2*ghost)
    idn=0
    
    do k=0,local_length(3)-1
       do j=0,local_length(2)-1
          do i=0,local_length(1)-1
             x=i+local_start(1)
             y=j+local_start(2)
             z=k+local_start(3)
             !Coordinates obtained.

             !Convert 3D coordinates to 1D index
             id=(i+ghost)+dy*(j+ghost)+dz*(k+ghost)

             !Store fluid bulk IDs in vector fluid_id 
             fluid_id(idn)=id
             idn=idn+1                                
          enddo
       enddo
    enddo

    size_fluid=idn

    
    
  endsubroutine SetGeometry

  !--------------------------------------
  subroutine InitUP
    !This subroutine initializes velocity and pressure fields. This suboutine can be offloaded to devices.
    integer i,j,k,x,y,z,id
    double precision dis1, dis2, phi, tau
    
    
    
    !Initilize density, relaxation frequency and chemical potential
    !$OMP TARGET TEAMS DISTRIBUTE collapse(3) private(x,y,z,id,dis1,dis2,phi)
    do k=0,local_length(3)-1
       do j=0,local_length(2)-1
          do i=0,local_length(1)-1
             
             id=(i+ghost)+dy*(j+ghost)+dz*(k+ghost)
             
             x=i+local_start(1)
             y=j+local_start(2)
             z=k+local_start(3)
             !Coordinates obtained.

             !Sphere
             dis1=sqrt(dble((x-nx/2.d0)**2+(y-pos1+ep/2)**2+(z-nz/2.d0)**2))-radius
             dis2=sqrt(dble((x-nx/2.d0)**2+(y-pos2-ep/2)**2+(z-nz/2.d0)**2))-radius
             !phi: phase field
             phi=(tanh(2.d0*min(dis1,dis2)/ep)+1.d0)*0.5d0

             
             rho(id) = phi*rhov+(1.d0-phi)*rhol

             tau=phi*tauv+(1.d0-phi)*taul
             ome(id) = 1.d0/(tau+0.5d0)
             
             cp(id)=0.d0

          enddo
       enddo
    enddo
    !$OMP END TARGET TEAMS DISTRIBUTE
  endsubroutine InitUP

  

endmodule initialization
