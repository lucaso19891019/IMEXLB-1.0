!Geng Liu, April 2nd, 2021
module lbm
  !This module includes the subroutines for LBM algorithms.
  use initialization
  use omp_lib

  !MPI local maximum velocity
  double precision um
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
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD private(id,ind,iq,udu,edu) schedule(static,1)
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
          f(id*nq+iq)=t(iq)*(p(id)+rho0*(3*edu+4.5*edu**2-1.5*udu))
       enddo
    enddo
    !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD          
  endsubroutine InitPDF

  !------------------------------------------------------------------------
  !Collision: This subroutine is the collision process of LBM. The collision result of array f is written in array fb.
  subroutine Collision
    integer idn,id,iq,ind
    double precision udu,edu,feq

    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD private(id,ind,iq,udu,edu,feq) schedule(static,1)
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
          
          feq=t(iq)*(p(id)+rho0*(3*edu+4.5*edu**2-1.5*udu))
          fb(id*nq+iq)=(1-ome)*f(id*nq+iq)+ome*feq
       enddo
    enddo
    !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD

  endsubroutine Collision

  !---------------------------------------------
  !Propagation: This subroutine includes the propagation process of LBM, as well as the bounce back process. The propagation result of array fb is written in array f.
  subroutine Propagation
    integer idn,id,iq

    !Perfect Shift
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD private(id,iq) schedule(static,1)
    do idn=0,size_fluid-1
       id=fluid_id(idn)      
       do iq=0,nq-1
          f(id*nq+iq)=fb((id-e(iq*dim)-(local_length(1)+2*ghost)*e(iq*dim+1))*nq+iq)
       enddo
    enddo
    !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
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
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD private(id,iq) nowait schedule(static,1)
    do idn=0,u_size-1
       id=bu(idn)
       do iq=0,nq-1
          fb(id*nq+iq) = fb((id+2*e(iq*dim)+2*(local_length(1)+2*ghost)*e(iq*dim+1))*nq+nq-1-iq)        
       enddo
    enddo
    !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD

    !Down Wall
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD private(id,iq) nowait schedule(static,1)
    do idn=0,d_size-1
       id=bd(idn)
       do iq=0,nq-1
          fb(id*nq+iq) = fb((id+2*e(iq*dim)+2*(local_length(1)+2*ghost)*e(iq*dim+1))*nq+nq-1-iq)        
       enddo
    enddo
    !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD

    !User Boundary
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD private(id,iq) schedule(static,1)
    do idn=0,user_size-1
       id=b_user(idn)
       do iq=0,nq-1
          fb(id*nq+iq) = fb((id+2*e(iq*dim)+2*(local_length(1)+2*ghost)*e(iq*dim+1))*nq+nq-1-iq)        
       enddo
    enddo
    !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
!--------------------------------------------------------------------------
    !Inlet condition: gradient free condition in boundary normal direction
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD private(id,iq) nowait schedule(static,1)
    do idn=0,l_size-1
       id=bl(idn)
       do iq=0,nq-1
          fb(id*nq+iq) = fb((id+1)*nq+iq)        
       enddo
    enddo
    !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD

    !Outlet condition: gradient free condition along characteristics
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD private(id,iq) nowait schedule(static,1)
    do idn=0,r_size-1
       id=br(idn)
       do iq=0,nq-1
          fb(id*nq+iq) = fb((id+e(iq*dim)+(local_length(1)+2*ghost)*e(iq*dim+1))*nq+iq)        
       enddo
    enddo
    !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD

  endsubroutine BoundaryCondition

  !---------------------------------------------
  !PostProcessing: This subroutine evalutes the physical properties, including pressure and momentum.
  subroutine PostProcessing

    integer idn,id,iq,y

    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD private(id,ind,iq) schedule(static,1)
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
    !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD    
!----------------------------------------------------------------------------------
    !Boundary Condition
    !Dirichlet boundary condition for inlet
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD private(id,y) nowait schedule(static,1)
    do idn=0,l_size-1
       id=bl(idn)
       y=id/(local_length(1)+2*ghost)-ghost+local_start(2)
       u(id*dim)=uu*rho0*4.d0*y*(ny-y)/ny**2
       u(id*dim+1)=0.d0
    enddo
    !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
    !0 pressure condition for outlet
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD private(id) schedule(static,1)
    do idn=0,r_size-1
       id=br(idn)
       p(id)=0
    enddo
    !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD

  endsubroutine PostProcessing

  !---------------------------------------------
  !Output: This subourtine writes the result to output files.
  subroutine Write
    character filename*20,format*100
    integer i,j,id
    double precision x,y
    integer file,request,contig_type,write_2d_type
    character buffer*120
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

    write(buffer,'(A,I3,A,I3,A)')"TITLE = ""Cylinder"""//NEW_LINE('A')//"VARIABLES = ""X"", ""Y"", ""p"", ""ux"", ""uy"""//NEW_LINE('A')//"ZONE I = ",nx+1,", J = ",ny+1,NEW_LINE('A')
    print_size=LEN_TRIM(buffer)
    offset=0
    if(rank==0)then
       call MPI_FILE_SEEK(file,offset,MPI_SEEK_SET,ierr)
       call MPI_FILE_IWRITE(file,buffer,print_size,MPI_CHARACTER,request,ierr)
    endif
    call MPI_BARRIER(CART_COMM,ierr)
    offset=offset+print_size

    num_var=3
    num_digit=14
    call MPI_TYPE_CONTIGUOUS(local_length(1)*(num_var+dim)*num_digit,MPI_CHARACTER,contig_type,ierr)
    extent=(nx+1)*(num_var+dim)*num_digit
    call MPI_TYPE_CREATE_RESIZED(contig_type,0,extent,write_2d_type,ierr)
    call MPI_TYPE_COMMIT(write_2d_type,ierr)
    offset=offset+(local_start(1)+local_start(2)*(nx+1))*(num_digit*(num_var+dim))
    call MPI_FILE_SET_VIEW(file,offset,MPI_CHARACTER,write_2d_type,"native",MPI_INFO_NULL,ierr)
    
    write(format,'(A)')"(sp,es13.6e2,X,es13.6e2,3(X,es13.6e2),A)"
    
    do j=0,local_length(2)-1 
       do i=0,local_length(1)-1          
          x=local_start(1)+i
          y=local_start(2)+j
          id=(i+ghost)+(local_length(1)+2*ghost)*(j+ghost)
          
          write(buffer,format)x,y,p(id),u(id*dim),u(id*dim+1),NEW_LINE('A')
          print_size=LEN_TRIM(buffer)
          call MPI_FILE_IWRITE(file,TRIM(buffer),print_size,MPI_CHARACTER,request,ierr)

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
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD map(tofrom:um) private(id) reduction(max:um) schedule(static,1)
    do idn=1,size_fluid-1
       id=fluid_id(idn)
       um=max(um,sqrt(u(id*dim)**2+u(id*dim+1)**2))
    enddo
    !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
    
    !MPI reduction
    call MPI_REDUCE(um,um_global,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,CART_COMM,ierr)
    !Print
    if(rank.eq.0)then         
       write(*,*)"Umax=",um_global
       write(*,*)"-------------------"              
    endif
    call MPI_BARRIER(CART_COMM,ierr)
  end subroutine Monitor

endmodule lbm
