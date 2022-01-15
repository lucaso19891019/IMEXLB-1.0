!Geng Liu, April 2nd, 2021
module initialization
  !This module includes the variables and subroutines for fluid initialization.
  use cart_mpi
  use omp_lib
  implicit none

  !charlength: characteristic length
  !pos: cylinder flow geometry related: position of cylinder in x-direction
  !re: Reynolds number
  !uu: characteristic velocity magnitude
  !ma: Mach number
  !radius: cylinder flow geometry related: cylinder radius, half of charlength in this case
  double precision charlength,pos,re,uu,ma,radius
  !nu: kinematic viscosity
  !rho0: fluid density
  !dt: time increment
  !ome: relaxation parameter 1/(tau+0.5)
  !tau: relaxation time
  !t_intv: dimensionless time increment
  !inter: time interval size(dimensionless) for monitor/output
  !max_t: maximum time(dimensionless)
  double precision nu,rho0,dt,ome,tau,t_intv,inter,max_t
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

  !p: pressure array
  double precision,dimension(:),allocatable::p
  !geo: geometry flag array
  integer,dimension(:),allocatable::geo
  !Vectors that store boundary points' indices: l-left, r-right, u-top, d-bottom, _user-user defined
  integer,dimension(:),allocatable::bl,br,bu,bd,b_user
  !Vector that store fluid bulk points' indices:
  integer,dimension(:),allocatable::fluid_id
  !dy: The index increase in 1D arrays reshaped from 3D arrays per 1 y-coordinate increment.
  integer dy
  !Size of index vectors
  integer size_fluid,l_size,r_size,u_size,d_size,user_size

  !u: velocity array
  double precision,dimension(:),allocatable::u

  
       
contains
  !DataInput: This subroutine reads data from "input.in", and generates derived initial data.
  subroutine DataInput
    
    open(100,file='input.in',status='old')
    read(100,*)nx,ny     
    read(100,*)charlength,pos,re
    !charlength & pos are given in the unit of nx, need to be converted to lattice unit.
    charlength=charlength*nx
    pos=pos*nx
    !Read relaxation parameter
    read(100,*)ome,dt,rho0
    nu=(2.d0-ome)*dt/(6.d0*ome)    
    !Read max time
    read(100,*)max_t,inter
    t_intv=charlength**2/nu
    interv=inter*t_intv
    max_step=max_t*t_intv
    close(100)

    !Derived phyical property
    uu=re*nu/charlength
    tau=3*nu
    ma=uu*sqrt(3.d0)
    radius=charlength*0.5d0  

    !Display
    if(rank.eq.0)then

       write(*,'(A24,F6.2)')   "Reynolds Number      :",re
       write(*,'(A24,F6.2)')   "Characteristic Length:",charlength
       write(*,'(A24,F6.2)')   "Ma                   :",ma
       write(*,'(A24,F6.2)')   "Kinematic Viscosity  :",nu
       write(*,'(A24,F6.2)')   "Relaxation Parameter :",ome
       write(*,'(A24,F6.2)')   "Relaxation Time      :",tau

       write(*,'(A24,F4.1)')   "Density              :",rho0

       write(*,'(A24,I8)')     "Maximum # of Steps   :",max_step
       write(*,'(A24,I8)')     "Interval Size        :",interv

       write(*,*)"============================="

    endif

  endsubroutine DataInput

  !---------------------------------------------
  !InitLattice: This subroutine defines the lattice velocites and weights
  !Numbering of D2Q9 lattice directions:
  !   2  5  8
  !    \ | /
  !   1--4--7
  !    / | \
  !   0  3  6
  subroutine InitLattice
    integer i,j,iq   
    !Lattice velocities and weights(2D)
    do j=-1,1
       do i=-1,1
          iq=(i+1)*3+(j+1)
          e(iq*dim)=i
          e(iq*dim+1)=j
          t(iq)=4*(1-0.75*i**2)*(1-0.75*j**2)/9.d0
       enddo
    enddo
  endsubroutine InitLattice

  !---------------------------------------------
  !AllocateArrays: This subroutine allocates local arrays.
  subroutine AllocateArrays
    
    allocate(f(0:array_size*nq-1))
    allocate(fb(0:array_size*nq-1))


    allocate(p(0:array_size-1))
    allocate(geo(0:array_size-1))
    allocate(fluid_id(0:array_size-1))
    allocate(b_user(0:array_size-1))

    allocate(u(0:array_size*dim-1))

    allocate(br(0:local_length(2)+2))
    allocate(bl(0:local_length(2)+2))
    allocate(bu(0:local_length(1)+2))
    allocate(bd(0:local_length(1)+2))
  endsubroutine AllocateArrays
  
  !-----------------------------------------------------
  !DeAllocateArrays: This subroutine deallocates local arrays.
  subroutine DeAllocateArrays
    deallocate(f)
    deallocate(fb)

    deallocate(p)
    deallocate(geo)
    deallocate(fluid_id)
    deallocate(b_user)
    
    deallocate(u)
    
    deallocate(br)
    deallocate(bl)
    deallocate(bu)
    deallocate(bd)
  endsubroutine DeAllocateArrays  

  !---------------------------------------------
  !SetGeometry: This subroutine initializes the flag array according to the given geometry
  subroutine SetGeometry 
    integer i,j,id,idn,x,y,idl,idr,idu,idd,id_user,iq
    logical flag
    ! dy is a global variable that is generated here, and should not be modified in other subroutines or functiions, dx is 1 and thus not needed.
    dy=(local_length(1)+2*ghost)
    
    idn=0
    idl=0
    idr=0
    idu=0
    idd=0

    do j=0,local_length(2)+2
       do i=0,local_length(1)+2
          x=i-1
          y=j-1
          id=(i-1+ghost)+dy*(j-1+ghost)
          
          if(x>=0.and.x<local_length(1).and.y>=0.and.y<local_length(2))then
             x=x+local_start(1)
             y=y+local_start(2)
             !Coordinates obtained.
             
             geo(id)=0
             !Cylinder
             if(dble((x-pos)**2+(y-ny/2.d0)**2).gt.radius**2)then
                geo(id)=1
                fluid_id(idn)=id
                idn=idn+1

             endif
          else
             x=x+local_start(1)
             y=y+local_start(2)
             !Walls
             if(x.eq.-1)then
                bl(idl)=id
                idl=idl+1               
             endif
             if(x.eq.nx+1)then
                br(idr)=id
                idr=idr+1
             endif
             if(y.eq.-1)then
                bd(idd)=id
                idd=idd+1
             endif
             if(y.eq.ny+1)then
                bu(idu)=id
                idu=idu+1
             endif
          endif
       enddo
    enddo

    !Flag communication
    call PassInt(geo)
    
    size_fluid=idn
    l_size=idl
    r_size=idr
    u_size=idu
    d_size=idd
    
    !User defined boundary (cylinder)
    flag=.false.
    id_user=0
    do j=0,local_length(2)+2
       do i=0,local_length(1)+2
          x=i-1
          y=j-1
          id=(i-1+ghost)+dy*(j-1+ghost)
          do iq=0,nq-1
             if(geo(id)==0.and.geo(id+e(iq*dim)+dy*(iq*dim+1))==1)then
                b_user(id_user)=id
                id_user=id_user+1
                flag=.true.
                exit
             endif
          enddo
          if(flag)then
             flag=.false.
             exit
          endif
       enddo
    enddo
    user_size=id_user
    
  endsubroutine SetGeometry

  !--------------------------------------
  !This subroutine initializes velocity and pressure fields. This suboutine can be offloaded to devices.
  subroutine InitUP
    integer i,j,id
    double precision r,y
    double precision r_vec(0:99)
    !Generate random number vector on CPU
    call random_number(r_vec) 
    
    !Initilize momentum, pressure
    !$OMP TARGET TEAMS DISTRIBUTE collapse(2) map(to:r_vec) private(r,y,id)
    do j=0,local_length(2)-1
       do i=0,local_length(1)-1

          id=i+ghost+dy*(j+ghost)
          !Perturbation
          r=1.d0+(r_vec(mod(id,100))-0.5d0)*2.d0*0.2d0
          y=j+local_start(2)
          
          p(id)=0.d0
          u(id*dim)=uu*4.d0*y*(ny-y)/ny**2*geo(id)*r 
          u(id*dim+1)=0.d0
       enddo
    enddo
    !$OMP END TARGET TEAMS DISTRIBUTE
  endsubroutine InitUP

endmodule initialization
