!Geng Liu, April 2nd, 2021
module lbm
  !This module includes the subroutines for LBM algorithms.
  use initialization
  use omp_lib

  !MPI local maximum velocity
  double precision um

  !Writing offset
  integer(kind=MPI_OFFSET_KIND)offset
contains
  !Definition to some temporary variables in subroutines:
  !idn: index for id vectors
  !id: index for arrays
  !iq: index for lattice directions
  !io: lattice direction opposite to iq
  !ind: index for dimensions
  !udu: magnitude of velocity squared
  !edu: inner product of physical velocity and lattice velocity
  !feq: equalibrium distributions

  !---------------------------------------------
  !InitPDF: This subroutine initializes the momentum, pressure and PDFs.
  subroutine InitPDF
    integer idn,id,iq,ind
    double precision udu,edu

    !Initialize PDF
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO private(id,ind,iq,udu,edu) schedule(static,1)
    do idn=0,size_fluid-1
       id=fluid_id(idn)
       udu=0.d0
       do ind=0,dim-1
          udu=udu+u(id*dim+ind)**2
       enddo

       do iq=0,nq-1
          edu=0.d0
          do ind=0,dim-1
             edu=edu+e(iq*dim+ind)*u(id*dim+ind)
          enddo

          !Initialize the PDFs with equilibrium values (incompressible model)
          f(id*nq+iq)=t(iq)*(p(id)*3.d0+rho0*(3*edu+4.5*edu**2-1.5*udu))
       enddo
    enddo
    !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO         
  endsubroutine InitPDF

  !------------------------------------------------------------------------
  !Collision: This subroutine is the collision process of LBM. The collision result of array f is written in array fb.
  subroutine Collision
    integer idn,id,iq,ind
    double precision udu,edu,feq

    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO private(id,ind,iq,udu,edu,feq) schedule(static,1)
    do idn=0,size_fluid-1
       id=fluid_id(idn)
       udu=0.d0
       do ind=0,dim-1
          udu=udu+u(id*dim+ind)**2
       enddo

       do iq=0,nq-1
          edu=0.d0
          do ind=0,dim-1
             edu=edu+e(iq*dim+ind)*u(id*dim+ind)
          enddo
          
          feq=t(iq)*(p(id)*3.d0+rho0*(3*edu+4.5*edu**2-1.5*udu))
          fb(id*nq+iq)=(1-ome)*f(id*nq+iq)+ome*feq
       enddo
    enddo
    !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

  endsubroutine Collision

  !---------------------------------------------
  !Propagation: This subroutine includes the propagation process of LBM, as well as the bounce back process. The propagation result of array fb is written in array f.
  subroutine Propagation
    integer idn,id,iq

    !Perfect Shift
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO private(id,iq) schedule(static,1)
    do idn=0,size_fluid-1
       id=fluid_id(idn)      
       do iq=0,nq-1
          f(id*nq+iq)=fb((id-e(iq*dim)-(local_length(1)+2*ghost)*e(iq*dim+1)-(local_length(1)+2*ghost)*(local_length(2)+2*ghost)*e(iq*dim+2))*nq+iq)
       enddo
    enddo
    !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
  endsubroutine Propagation

  !---------------------------------------------
  !BoundaryCondition: This subroutine applies the boundary conditions other than no-slip conditions
  subroutine BoundaryCondition
    integer idn,id,iq

    !OpenMP and MPI synchronization before propagation:
    !$OMP TARGET UPDATE from(fb)
    call PassF(fb)
    !$OMP TARGET UPDATE to(fb)
    

    !Up Wall
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO private(id,iq) nowait schedule(static,1)
    do idn=0,u_size-1
       id=bu(idn)
       do iq=0,nq-1
          fb(id*nq+iq) = fb((id+2*e(iq*dim)+2*(local_length(1)+2*ghost)*e(iq*dim+1)+2*(local_length(1)+2*ghost)*(local_length(2)+2*ghost)*e(iq*dim+2))*nq+nq-1-iq)        
       enddo
    enddo
    !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

    !Down Wall
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO private(id,iq) nowait schedule(static,1)
    do idn=0,d_size-1
       id=bd(idn)
       do iq=0,nq-1
          fb(id*nq+iq) = fb((id+2*e(iq*dim)+2*(local_length(1)+2*ghost)*e(iq*dim+1)+2*(local_length(1)+2*ghost)*(local_length(2)+2*ghost)*e(iq*dim+2))*nq+nq-1-iq)        
       enddo
    enddo
    !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

    !User Boundary
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO private(id,iq) schedule(static,1)
    do idn=0,user_size-1
       id=b_user(idn)
       do iq=0,nq-1
          fb(id*nq+iq) = fb((id+2*e(iq*dim)+2*(local_length(1)+2*ghost)*e(iq*dim+1)+2*(local_length(1)+2*ghost)*(local_length(2)+2*ghost)*e(iq*dim+2))*nq+nq-1-iq)        
       enddo
    enddo
    !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
!--------------------------------------------------------------------------
    !Inlet condition: gradient free condition in boundary normal direction
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO private(id,iq) nowait schedule(static,1)
    do idn=0,l_size-1
       id=bl(idn)
       do iq=0,nq-1
          fb(id*nq+iq) = fb((id+1)*nq+iq)        
       enddo
    enddo
    !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

    !Outlet condition: gradient free condition along characteristics
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO private(id,iq) nowait schedule(static,1)
    do idn=0,r_size-1
       id=br(idn)
       do iq=0,nq-1
          fb(id*nq+iq) = fb((id+e(iq*dim)+(local_length(1)+2*ghost)*e(iq*dim+1)+(local_length(1)+2*ghost)*(local_length(2)+2*ghost)*e(iq*dim+2))*nq+iq)        
       enddo
    enddo
    !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

  endsubroutine BoundaryCondition

  !---------------------------------------------
  !PostProcessing: This subroutine evalutes the physical properties, including pressure and momentum.
  subroutine PostProcessing

    integer idn,id,iq,z

    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO private(id,ind,iq) schedule(static,1)
    do idn=0,size_fluid-1
       id=fluid_id(idn)
       p(id)=0
       do ind=0,dim-1
          u(id*dim+ind)=0
       enddo
       do iq=0,nq-1
          p(id)=p(id)+f(id*nq+iq)
          do ind=0,dim-1
             u(id*dim+ind)=u(id*dim+ind)+f(id*nq+iq)*e(iq*dim+ind)
          enddo
       enddo
       do ind=0,dim-1
          u(id*dim+ind)=u(id*dim+ind)/rho0*geo(id)
       enddo
       p(id)=p(id)/3.d0*geo(id)
    enddo
    !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO 
!----------------------------------------------------------------------------------
    !Boundary Condition
    !Dirichlet boundary condition for inlet
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO private(id,y) nowait schedule(static,1)
    do idn=0,l_size-1
       id=bl(idn)+1
       z=id/((local_length(1)+2*ghost)*(local_length(2)+2*ghost))-ghost+local_start(3)
       u(id*dim)=uu*rho0*4.d0*z*(nz-z)/nz**2
       u(id*dim+1)=0.d0
    enddo
    !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
    !0 pressure condition for outlet
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO  private(id) schedule(static,1)
    do idn=0,r_size-1
       id=br(idn)-1
       p(id)=0
    enddo
    !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

  endsubroutine PostProcessing

  !---------------------------------------------
  !Output: This subourtine writes the result to output files.
  subroutine Write
    character filename*20,format*100
    integer i,j,k,id
    double precision x,y,z
    integer file,request,contig_type,contig_type_2d,write_2d_type,write_3d_type
    character buffer*150
    integer print_size,num_var,num_digit
    integer(kind=MPI_OFFSET_KIND)offset
    integer(kind=MPI_ADDRESS_KIND)extent

    !$OMP TARGET UPDATE from(u,p)
    
    write(filename,'(A)')"output.dat"
    call MPI_FILE_OPEN(CART_COMM,filename,MPI_MODE_CREATE+MPI_MODE_EXCL+MPI_MODE_WRONLY,mpi_INFO_NULL,file,ierr)
    if(ierr.ne.MPI_SUCCESS)then
       if(rank==0)then
          call MPI_FILE_DELETE(filename,MPI_INFO_MULL,ierr)
       endif
       call MPI_FILE_OPEN(CART_COMM,filename,MPI_MODE_CREATE+MPI_MODE_WRONLY,mpi_INFO_NULL,file,ierr)
    endif

    write(buffer,'(A,I3,A,I3,A,I3,A)')"TITLE = ""Cylinder"""//NEW_LINE('A')//"VARIABLES = ""X"", ""Y"", ""Z"", ""p"", ""ux"", ""uy"", ""uz"""//NEW_LINE('A')//"ZONE I = ",nx+1,", J = ",ny+1,", K = ",nz+1,NEW_LINE('A')
    print_size=LEN_TRIM(buffer)
    offset=0
    if(rank==0)then
       call MPI_FILE_SEEK(file,offset,MPI_SEEK_SET,ierr)
       call MPI_FILE_IWRITE(file,buffer,print_size,MPI_CHARACTER,request,ierr)
    endif
    call MPI_BARRIER(CART_COMM,ierr)
    offset=offset+print_size

    num_var=4
    num_digit=14
    call MPI_TYPE_CONTIGUOUS(local_length(1)*(num_var+dim)*num_digit,MPI_CHARACTER,contig_type,ierr)
    extent=(nx+1)*(num_var+dim)*num_digit
    call MPI_TYPE_CREATE_RESIZED(contig_type,0,extent,write_2d_type,ierr)
    call MPI_TYPE_COMMIT(write_2d_type,ierr)
    call MPI_TYPE_CONTIGUOUS(local_length(2),write_2d_type,contig_type_2d,ierr)
    extent=(nx+1)*(ny+1)*(num_var+dim)*num_digit
    call MPI_TYPE_CREATE_RESIZED(contig_type_2d,0,extent,write_3d_type,ierr)
    call MPI_TYPE_COMMIT(write_3d_type,ierr)


    
    offset=offset+(local_start(1)+local_start(2)*(nx+1)+local_start(3)*(ny+1)*(nx+1))*(num_digit*(num_var+dim))
    call MPI_FILE_SET_VIEW(file,offset,MPI_CHARACTER,write_3d_type,"native",MPI_INFO_NULL,ierr)
    write(format,'(A)')"(sp,es13.6e2,X,es13.6e2,X,es13.6e2,4(X,es13.6e2),A)"
    do k=0,local_length(3)-1
     do j=0,local_length(2)-1
       do i=0,local_length(1)-1
          
          x=local_start(1)+i
          y=local_start(2)+j
          z=local_start(3)+k
          id=(i+ghost)+(local_length(1)+2*ghost)*(j+ghost)+(local_length(1)+2*ghost)*(local_length(2)+2*ghost)*(k+ghost)
          
          write(buffer,format)x,y,z,p(id),u(id*dim),u(id*dim+1),u(id*dim+2),NEW_LINE('A')
          print_size=LEN_TRIM(buffer)
          call MPI_FILE_IWRITE(file,TRIM(buffer),print_size,MPI_CHARACTER,request,ierr)

       enddo
    enddo
    enddo
    
    call MPI_BARRIER(CART_COMM,ierr)
    call MPI_FILE_CLOSE(file,ierr)

  endsubroutine Write

  !------------------------------------------
  !Monitor: This subroutine prints the maximum magnitude of velocity
  subroutine Monitor
    integer idn,id
    double precision um_global

    !Find the MPI local maximum magnitude of velocity
    um=0
    !OpenMP reduction
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO map(tofrom:um) private(id) reduction(max:um) schedule(static,1)
    do idn=1,size_fluid-1
       id=fluid_id(idn)
       um=max(um,sqrt(u(id*dim)**2+u(id*dim+1)**2))
    enddo
    !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
    call MPI_BARRIER(CART_COMM,ierr)
    !MPI reduction
    call MPI_REDUCE(um,um_global,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,CART_COMM,ierr)
    !Print
    if(rank.eq.0)then
       write(*,*)"T=",iter/t_intv
       write(*,*)"Umax=",um_global/(charlength*t_intv)
       write(*,*)"-------------------"              
    endif
    call MPI_BARRIER(CART_COMM,ierr)
  end subroutine Monitor

  !---------------------------------------------
  !Output: This subourtine writes the result to binary output files.
  subroutine WriteBinary
    character filename*20,num*3
    integer i,j,k,ind
    integer numvars !Number of Variables, including the coordinates
    double precision uxmax,       uymax,       uzmax,       pmax
    double precision uxmin,       uymin,       uzmin,       pmin
    double precision uxmax_global,uymax_global,uzmax_global,pmax_global
    double precision uxmin_global,uymin_global,uzmin_global,pmin_global
    integer file,request,write_3d_type,write_type
    integer print_size
    
    double precision, dimension(:),allocatable:: buffer

    !$OMP TARGET UPDATE from(u,p)

    write(num,'(I3.3)')count
    write(filename,'(A)')"output_"//num//".plt"
    call MPI_FILE_OPEN(CART_COMM,filename,MPI_MODE_CREATE+MPI_MODE_EXCL+MPI_MODE_WRONLY,mpi_INFO_NULL,file,ierr)
    if(ierr.ne.MPI_SUCCESS)then
       if(rank==0)then
          call MPI_FILE_DELETE(filename,MPI_INFO_MULL,ierr)
       endif
       call MPI_FILE_OPEN(CART_COMM,filename,MPI_MODE_CREATE+MPI_MODE_WRONLY,mpi_INFO_NULL,file,ierr)
    endif

    numvars=3+4 !3 coordinates (x,y,z), 4 fields(p,ux,uy,uz)
    
    if(rank==0)then
       offset=0
       call MPI_FILE_SEEK(file,offset,MPI_SEEK_SET,ierr)
       
       print_size=8       
       call MPI_FILE_IWRITE(file,'#!TDV112',print_size,MPI_CHARACTER,request,ierr)
       offset=offset+print_size
       
       print_size=4+numvars  +6           +1         +1         +1 +1 +2 +2 +2
            !                 Len(Title)|  len(var1)| len(var2)|...        
       call MPI_FILE_IWRITE(file,(/1,0,        str2ascii('Sphere',6),0,  numvars,&
            !                        FullData| Title
            ichar('x'),0,ichar('y'),0,ichar('z'),0,ichar('p'),0,str2ascii('ux',2),0,&
            str2ascii('uy',2),0,str2ascii('uz',2),0/)&
            !Variable Names            
            ,print_size,MPI_INTEGER,request,ierr)
       !Note that integer 0 comes after any string
       offset=offset+print_size*4
       
       call MPI_FILE_IWRITE(file,299.0,1,MPI_REAL,request,ierr)
       offset=offset+4
       
       print_size=8+3
       call MPI_FILE_IWRITE(file,(/str2ascii('ZONE '//num,8),0,-1,0/)&
            !                      ZoneName|                      StrandID
            ,print_size,MPI_INTEGER,request,ierr)
       offset=offset+print_size*4
       
       call MPI_FILE_IWRITE(file,iter/t_intv,1,MPI_DOUBLE_PRECISION,request,ierr)
            !                    Solution Time
       offset=offset+8

       print_size=9
       call MPI_FILE_IWRITE(file,(/-1,     0,0,0,0,nx+1,ny+1,nz+1,0/)&
            !                      Ordered|        IMax|JMax|KMax
            ,print_size,MPI_INTEGER,request,ierr)
       offset=offset+print_size*4

       call MPI_FILE_IWRITE(file,(/357.0,299.0/),2,MPI_REAL,request,ierr)
       offset=offset+8

       print_size=numvars+3
       call MPI_FILE_IWRITE(file,(/2, 2, 2, 2, 2, 2, 2,       0,0,-1/)&
            !            VariableFormat:1.Float,2.Double
            ,print_size,MPI_INTEGER,request,ierr)
       offset=offset+print_size*4      
    endif
    call MPI_BARRIER(CART_COMM,ierr)

    !Find the MPI local max/min of output variables
    uxmax=0.d0
    uxmin=0.d0    
    uymax=0.d0
    uymin=0.d0
    uzmax=0.d0
    uzmin=0.d0
    pmax=0.d0        
    pmin=0.d0
    !OpenMP reduction
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO map(tofrom:uxmax,uxmin,uymax,uymin,uzmax,uzmin,pmax,pmin) private(id) reduction(max:uxmax,uymax,uzmax,pmax) reduction(min:uxmmin,uymin,uzmin,pmin) schedule(static,1)
    do idn=1,size_fluid-1
       id=fluid_id(idn)
       uxmax=max(uxmax,u(id*dim))
       uxmin=min(uxmin,u(id*dim))
       uymax=max(uymax,u(id*dim+1))
       uymin=min(uymin,u(id*dim+1))
       uzmax=max(uzmax,u(id*dim+2))
       uzmin=min(uzmin,u(id*dim+2))
       pmax=max(pmax,p(id))      
       pmin=min(pmin,p(id))
    enddo
    !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
    call MPI_BARRIER(CART_COMM,ierr)
    !MPI reduction
    call MPI_REDUCE(uxmax,uxmax_global,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,CART_COMM,ierr)
    call MPI_REDUCE(uxmin,uxmin_global,1,MPI_DOUBLE_PRECISION,MPI_MIN,0,CART_COMM,ierr)
    call MPI_REDUCE(uymax,uymax_global,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,CART_COMM,ierr)
    call MPI_REDUCE(uymin,uymin_global,1,MPI_DOUBLE_PRECISION,MPI_MIN,0,CART_COMM,ierr)
    call MPI_REDUCE(uzmax,uzmax_global,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,CART_COMM,ierr)
    call MPI_REDUCE(uzmin,uzmin_global,1,MPI_DOUBLE_PRECISION,MPI_MIN,0,CART_COMM,ierr)
    call MPI_REDUCE(pmax,pmax_global,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,CART_COMM,ierr)
    call MPI_REDUCE(pmin,pmin_global,1,MPI_DOUBLE_PRECISION,MPI_MIN,0,CART_COMM,ierr)
    
    !Write variable limits
    if(rank==0)then
       call MPI_FILE_IWRITE(file,(/0.d0,(nx/charlength),0.d0,(ny/charlength),0.d0,(nz/charlength),&
       !                           xmin,xmax            ymin,ymax            zmin,zmax
            pmin_global*charlength/t_intv**2,pmax_global*charlength/t_intv**2,&
            uxmin_global/(charlength*t_intv),uxmax_global/(charlength*t_intv),&
            uymin_global/(charlength*t_intv),uymax_global/(charlength*t_intv),&
            uzmin_global/(charlength*t_intv),uzmax_global/(charlength*t_intv)/),&   
            numvars*2,       MPI_DOUBLE_PRECISION,request,ierr)
       offset=offset+numvars*2*8
    endif

    !Pass offset
    call MPI_Bcast(offset, 1, MPI_INTEGER, 0, CART_COMM,ierr)
    
    !SetWritingView    
    call MPI_TYPE_CREATE_SUBARRAY(dim,global_length,local_length,local_start,MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,write_3d_type,ierr)
    call MPI_TYPE_CONTIGUOUS(numvars*1,write_3d_type,write_type,ierr)
    !numvars*1 works with compiling option -fc=ifx, numvars does not work
    !Both numvars*1 and numvars does not work with -fc=ifort
    call MPI_TYPE_COMMIT(write_type,ierr)
    call MPI_FILE_SET_VIEW(file,offset,MPI_DOUBLE_PRECISION,write_type,"native",MPI_INFO_NULL,ierr)
    
    !Allocate writing buffer
    print_size=1
    do i=1,dim
       print_size=print_size*local_length(i)
    enddo     
    ALLOCATE(buffer(print_size))
    
    !Writing DATA
    do k=0,local_length(3)-1
       do j=0,local_length(2)-1
          do i=0,local_length(1)-1
             buffer(i+j*local_length(1)+k*local_length(1)*local_length(2)+1)=&
                  (local_start(1)+i)/charlength           
          enddo
       enddo
    enddo
    call MPI_FILE_WRITE_ALL(file,buffer,print_size,MPI_DOUBLE_PRECISION,request,ierr)

    do k=0,local_length(3)-1
       do j=0,local_length(2)-1
          do i=0,local_length(1)-1      
             buffer(i+j*local_length(1)+k*local_length(1)*local_length(2)+1)=&
                  (local_start(2)+j)/charlength
          enddo
       enddo
    enddo
    call MPI_FILE_WRITE_ALL(file,buffer,print_size,MPI_DOUBLE_PRECISION,request,ierr)

    do k=0,local_length(3)-1
       do j=0,local_length(2)-1
          do i=0,local_length(1)-1      
             buffer(i+j*local_length(1)+k*local_length(1)*local_length(2)+1)=&
                  (local_start(3)+k)/charlength
          enddo
       enddo
    enddo
    call MPI_FILE_WRITE_ALL(file,buffer,print_size,MPI_DOUBLE_PRECISION,request,ierr)

    do k=0,local_length(3)-1
       do j=0,local_length(2)-1
          do i=0,local_length(1)-1      
             buffer(i+j*local_length(1)+k*local_length(1)*local_length(2)+1)=&
                  p((i+ghost)+(local_length(1)+2*ghost)*(j+ghost)&
                  +(local_length(1)+2*ghost)*(local_length(2)+2*ghost)*(k+ghost))&
                  *charlength/t_intv**2!Unit conversion
          enddo
       enddo
    enddo
    call MPI_FILE_WRITE_ALL(file,buffer,print_size,MPI_DOUBLE_PRECISION,request,ierr)

    do ind=0,dim-1
       do k=0,local_length(3)-1
          do j=0,local_length(2)-1
             do i=0,local_length(1)-1      
                buffer(i+j*local_length(1)+k*local_length(1)*local_length(2)+1)=&
                     u(((i+ghost)+(local_length(1)+2*ghost)*(j+ghost)&
                     +(local_length(1)+2*ghost)*(local_length(2)+2*ghost)*(k+ghost))*dim+ind)&
                     /(charlength*t_intv)!Unit conversion
             enddo
          enddo
       enddo
       call MPI_FILE_WRITE_ALL(file,buffer,print_size,MPI_DOUBLE_PRECISION,request,ierr)
    enddo

    !Deallocate writing buffer
    DEALLOCATE(buffer)    
     
    call MPI_BARRIER(CART_COMM,ierr)
    call MPI_FILE_CLOSE(file,ierr)

  contains
    function str2ascii(s,n)result(o)
      character,dimension(30)::s
      integer n
      integer,dimension(n)::o
      integer i
      do i=1,n
         o(i)=ichar(s(i))
      enddo
    endfunction str2ascii
    
  endsubroutine WriteBinary

endmodule lbm
