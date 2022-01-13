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
  !flag: geometry flag array
  integer,dimension(:),allocatable::geo
  !Vectors that store boundary points' indices: l-left, r-right, u-top, d-bottom, _user-user defined
  integer,dimension(:),allocatable::bl,br,bu,bd,bf,bb,b_user
  !Vector that store fluid bulk points' indices:
  integer,dimension(:),allocatable::fluid_id
  !dy: The index increase in 1D arrays reshaped from 3D arrays per 1 y-coordinate increment.
  !dy: The index increase in 1D arrays reshaped from 3D arrays per 1 z-coordinate increment.
  integer dy,dz
  !Size of index vectors
  integer size_fluid,l_size,r_size,u_size,d_size,f_size,b_size,user_size

  !u: velocity array
  double precision,dimension(:),allocatable::u

  
       
contains
  !DataInput: This subroutine reads data from "input.in", and generates derived initial data.
  subroutine DataInput
    
    open(100,file='input_3d.in',status='old')
    read(100,*)nx,ny,nz     
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


    allocate(p(0:array_size-1))
    allocate(geo(0:array_size-1))
    allocate(fluid_id(0:array_size-1))
    allocate(b_user(0:array_size-1))

    allocate(u(0:array_size*dim-1))

    allocate(br(0:(local_length(2)+2*ghost)*(local_length(3)+2*ghost)))
    allocate(bl(0:(local_length(2)+2*ghost)*(local_length(3)+2*ghost)))
    allocate(bu(0:(local_length(2)+2*ghost)*(local_length(1)+2*ghost)))
    allocate(bd(0:(local_length(2)+2*ghost)*(local_length(1)+2*ghost)))
    allocate(bf(0:(local_length(1)+2*ghost)*(local_length(3)+2*ghost)))
    allocate(bb(0:(local_length(1)+2*ghost)*(local_length(3)+2*ghost)))
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
    deallocate(bf)
    deallocate(bb)
  endsubroutine DeAllocateArrays  

  !---------------------------------------------
  !SetGeometry: This subroutine initializes the flag array according to the given geometry
  subroutine SetGeometry 
    integer i,j,k,id,idn,x,y,z,idl,idr,idu,idf,idb,idd,id_user,iq
    logical flag
    idn=0
    idl=0
    idr=0
    idu=0
    idd=0
    idf=0
    idb=0

    do k=0,local_length(3)+2
    do j=0,local_length(2)+2
       do i=0,local_length(1)+2
          x=i-1
          y=j-1
          z=k-1
          id=(i-1+ghost)+(local_length(1)+2*ghost)*(j-1+ghost)+(local_length(1)+2*ghost)*(local_length(2)+2*ghost)*(k-1+ghost)
          
          if(x>=0.and.x<local_length(1).and.y>=0.and.y<local_length(2).and.z>=0.and.z<local_length(3))then
             x=x+local_start(1)
             y=y+local_start(2)
             z=z+local_start(3)
             !Coordinates obtained.

             geo(id)=0
             !Sphere
             if(dble((x-pos)**2+(y-ny/2.d0)**2+(z-nz/2.d0)**2).gt.radius**2)then
                geo(id)=1
                fluid_id(idn)=id
                idn=idn+1
             endif
          else
             x=x+local_start(1)
             y=y+local_start(2)
             z=z+local_start(3)
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
                bf(idf)=id
                idf=idf+1
             endif
             if(y.eq.ny+1)then
                bb(idb)=id
                idb=idb+1
             endif
             if(z.eq.-1)then
                bd(idd)=id
                idd=idd+1
             endif
             if(z.eq.nz+1)then
                bu(idu)=id
                idu=idu+1
             endif
          endif
       enddo
       enddo
    enddo

    !Flag communication
    call PassInt(geo)
    
    size_fluid=idn
    l_size=idl
    r_size=idr
    u_size=idu
    d_size=idd
    f_size=idf
    b_size=idb
    
    !User defined boundary (sphere)
    flag=.false.
    id_user=0
    do idn=0,size_fluid-1
       id=fluid_id(idn)
       do iq=0,nq-1
          if(geo(id)==0.and.geo(id+e(iq*dim)+(local_length(1)+2*ghost)*(iq*dim+1)+(local_length(1)+2*ghost)*(local_length(2)+2*ghost)*(iq*dim+2))==1)then
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
    user_size=id_user
    
  endsubroutine SetGeometry

  !--------------------------------------
  subroutine InitUP
    !This subroutine initializes velocity and pressure fields. This suboutine can be offloaded to devices.
    integer i,j,k,id
    double precision r,z
    double precision r_vec(0:99)
    !Generate random number vector on CPU
    call random_number(r_vec)
    
    !Initilize momentum, pressure
    !$OMP TARGET TEAMS DISTRIBUTE collapse(3) map(to:r_vec) private(r,z,id)
    do k=0,local_length(3)-1
       do j=0,local_length(2)-1
          do i=0,local_length(1)-1
             id=i+ghost+(local_length(1)+2*ghost)*(j+ghost)+(local_length(1)+2*ghost)*(local_length(2)+2*ghost)*(k+ghost)
             !Perturbation
             r=1.d0+(r_vec(mod(id,100))-0.5d0)*2.d0*0.2d0
             z=k+local_start(3)
             
             p(id)=0.d0
             u(id*dim)=uu*4.d0*z*(nz-z)/nz**2*r*geo(id) 
             u(id*dim+1)=0.d0
             u(id*dim+2)=0.d0
          enddo
       enddo
    enddo
    !$OMP END TARGET TEAMS DISTRIBUTE
  endsubroutine InitUP 

endmodule initialization
