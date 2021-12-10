!Geng Liu, April 2nd, 2021
module lbm
  !This module includes the subroutines for LBM algorithms.
  use initialization

  !Maximum velocity on current MPI processor
  double precision um

  !Writing offset
  integer(kind=MPI_OFFSET_KIND)offset
contains
  !Definition to some temporary variables in subroutines:
  !idn: index for index vectors
  !id: index for array elements
  !iq: index for lattice directions
  !io: lattice direction opposite to iq
  !ind: index for dimensions
  !udu: magnitude of velocity squared
  !edu: inner product of physical velocity and lattice velocity
  !feq: equalibrium distributions

  !---------------------------------------------
  !InitPDF: This subroutine initializes the PDFs. This subroutine can be offloaded to devices.
  subroutine InitPDF
    integer idn,id,iq,ind
    double precision udu,edu,gamma,eddrhoc,eddcpc,uddrhoc,uddcpc

    call Laplacian

    do idn=0,size_fluid-1
       id=fluid_id(idn)
       cp(id)=2.d0*beta*((rho(id)-rhol)*(rho(id)-rhov)*(2.d0*rho(id)-rhol-rhov)&
            -(ep*(rhol-rhov)/4.d0)**2*lap(id))
    enddo
    call PassD(cp)
    
    call GradC(dcrho,rho)
    call GradC(dccp,cp)
    
    !Initialize PDF
    !$OMP TARGET TEAMS DISTRIBUTE private(idn,id,ind,iq,udu,edu) 
    do idn=0,size_fluid-1
       id=fluid_id(idn)
       udu=0.d0
       uddrhoc=0.d0
       uddcpc=0.d0
       do ind=0,dim-1
          udu=udu+u(id*dim+ind)**2
          uddrhoc=uddrhoc+u(id*dim+ind)*dcrho(id*dim+ind)
          uddcpc=uddcpc+u(id*dim+ind)*dccp(id*dim+ind)
       enddo

       
       do iq=0,nq-1
          edu=0.d0
          do ind=0,dim-1
             edu=edu+e(iq*dim+ind)*u(id*dim+ind)
          enddo

          !Initialize the PDFs with equilibrium values
          gamma=t(iq)*(1.d0+(3*edu+4.5*edu**2-1.5*udu))

          eddrhoc=0.5d0*(rho(id+e(iq*dim)+dy*e(iq*dim+1)+dz*e(iq*dim+2))&
               -rho(id-e(iq*dim)-dy*e(iq*dim+1)-dz*e(iq*dim+2)))
          eddcpc=0.5d0*(cp(id+e(iq*dim)+dy*e(iq*dim+1)+dz*e(iq*dim+2))&
               -cp(id-e(iq*dim)-dy*e(iq*dim+1)-dz*e(iq*dim+2)))

          f(id*nq+iq)=gamma*rho(id)-0.5*gamma*(eddrhoc-uddrhoc-3.d0*rho(id)*(eddcpc-uddcpc))
       enddo
    enddo
    !$OMP END TARGET TEAMS DISTRIBUTE    
  endsubroutine InitPDF

  !------------------------------------------------------------------------
  !Collision: This subroutine is the collision process of LBM. The collision result of array f is stored in array fb. This subroutine can be offloaded to devices.
  subroutine Collision
    integer idn,id,iq,ind
    double precision udu,edu,feq,gamma,eddrhoc,eddcpc,uddrhoc,uddcpc,eddrhom,eddcpm,uddrhom,uddcpm

    call GradM(dmrho,rho)
    call GradM(dmcp,cp)
    
    !$OMP TARGET TEAMS DISTRIBUTE private(idn,id,ind,iq,udu,edu,feq)
    do idn=0,size_fluid-1
       id=fluid_id(idn)
       udu=0.d0
       uddrhoc=0.d0
       uddcpc=0.d0
       uddrhom=0.d0
       uddcpm=0.d0
       do ind=0,dim-1
          udu=udu+u(id*dim+ind)**2
          uddrhoc=uddrhoc+u(id*dim+ind)*dcrho(id*dim+ind)
          uddrhom=uddrhom+u(id*dim+ind)*dmrho(id*dim+ind)
          uddcpc=uddcpc+u(id*dim+ind)*dccp(id*dim+ind)
          uddcpm=uddcpm+u(id*dim+ind)*dmcp(id*dim+ind)
       enddo

       do iq=0,nq-1
          edu=0.d0
          do ind=0,dim-1
             edu=edu+e(iq*dim+ind)*u(id*dim+ind)
          enddo

          gamma=t(iq)*(1.d0+(3*edu+4.5*edu**2-1.5*udu))
          
          eddrhoc=0.5d0*(rho(id+e(iq*dim)+dy*e(iq*dim+1)+dz*e(iq*dim+2))&
               -rho(id-e(iq*dim)-dy*e(iq*dim+1)-dz*e(iq*dim+2)))
          eddcpc=0.5d0*(cp(id+e(iq*dim)+dy*e(iq*dim+1)+dz*e(iq*dim+2))&
               -cp(id-e(iq*dim)-dy*e(iq*dim+1)-dz*e(iq*dim+2)))

          feq=gamma*rho(id)-0.5*gamma*(eddrhoc-uddrhoc-3.d0*rho(id)*(eddcpc-uddcpc))

          eddrhom=0.25d0*(5.d0*rho(id+e(iq*dim)+dy*e(iq*dim+1)+dz*e(iq*dim+2))&
               -3.d0*rho(id)&
               -rho(id-e(iq*dim)-dy*e(iq*dim+1)-dz*e(iq*dim+2))&
               -rho(id+2*e(iq*dim)+2*dy*e(iq*dim+1)+2*dz*e(iq*dim+2)))
          eddcpm=0.25d0*(5.d0*cp(id+e(iq*dim)+dy*e(iq*dim+1)+dz*e(iq*dim+2))&
               -3.d0*cp(id)&
               -cp(id-e(iq*dim)-dy*e(iq*dim+1)-dz*e(iq*dim+2))&
               -cp(id+2*e(iq*dim)+2*dy*e(iq*dim+1)+2*dz*e(iq*dim+2)))

          fb(id*nq+iq)=(1-ome(id))*f(id*nq+iq)+ome(id)*feq&
               +gamma*(eddrhom-uddrhom-3.d0*rho(id)*(eddcpm-uddcpm))
       enddo
    enddo
    !$OMP END TARGET TEAMS DISTRIBUTE

  endsubroutine Collision

  !---------------------------------------------
  !Propagation: This subroutine is the propagation process of LBM. The propagation result of array fb is stored in array f. This subroutine can be offloaded to devices.
  subroutine Propagation
    integer idn,id,iq

    !Perfect Shift
    !$OMP TARGET TEAMS DISTRIBUTE private(idn,id,iq)
    do idn=0,size_fluid-1
       id=fluid_id(idn)      
       do iq=0,nq-1
          f(id*nq+iq)=fb((id-e(iq*dim)-dy*e(iq*dim+1)-dz*e(iq*dim+2))&
               *nq+iq)
       enddo
    enddo
    !$OMP END TARGET TEAMS DISTRIBUTE
  endsubroutine Propagation

  !---------------------------------------------
  !BoundaryCondition: This subroutine applies the boundary conditions. The boundary conditions should be pre-propagation conditions, post-propagation conditions should either be transformed to pre-propagation conditions, or be applied in the PostProcessing subroutine. Periodic condition is the default condition and doesn't require specification here. This subroutine can be offloaded to devices.
  subroutine BoundaryCondition
    integer idn,id,iq

    !OpenMP and MPI synchronization before propagation:
    !$OMP TARGET UPDATE from(fb)
    call PassF(fb)
    !$OMP TARGET UPDATE to(fb)
    

  endsubroutine BoundaryCondition

  !---------------------------------------------
  !PostProcessing: This subroutine evalutes the physical properties, including pressure and velocity. It can also include some post-propagation boundary conditions. This subroutine can be offloaded to devices. 
  subroutine PostProcessing

    integer idn,id,iq,ind
    double precision phi,tau

    !$OMP TARGET TEAMS DISTRIBUTE private(idn,id,iq)
    do idn=0,size_fluid-1
       id=fluid_id(idn)
       rho(id)=0.d0
       do iq=0,nq-1
          rho(id)=rho(id)+f(id*nq+iq)
       enddo
       phi=(rhol-rho(id))/(rhol-rhov)
       tau=phi*tauv+(1.d0-phi)*taul
       ome(id) = 1.d0/(tau+0.5d0)
    enddo
    !$OMP END TARGET TEAMS DISTRIBUTE

    call Laplacian

    do idn=0,size_fluid-1
       id=fluid_id(idn)
       cp(id)=2.d0*beta*((rho(id)-rhol)*(rho(id)-rhov)*(2.d0*rho(id)-rhol-rhov)&
            -(ep*(rhol-rhov)/4.d0)**2*lap(id))
    enddo
    call PassD(cp)

    call GradC(dcrho,rho)
    call GradC(dccp,cp)
    
    
    !$OMP TARGET TEAMS DISTRIBUTE private(idn,id,ind,iq)
    do idn=0,size_fluid-1
       id=fluid_id(idn)      
       do ind=0,dim-1
          u(id*dim+ind)=0
       enddo
       do iq=0,nq-1
          do ind=0,dim-1
             u(id*dim+ind)=u(id*dim+ind)+f(id*nq+iq)*e(iq*dim+ind)
          enddo
       enddo
       do ind=0,dim-1
          u(id*dim+ind)=(u(id*dim+ind)+&
               0.5*(dcrho(id*dim+ind)/3.d0-rho(id)*dccp(id*dim+ind)))/rho(id)
       enddo
    enddo
    !$OMP END TARGET TEAMS DISTRIBUTE


  endsubroutine PostProcessing

  !------------------------------------------
  !Monitor: This subroutine prints the user defined runtime global values, in this case maximum magnitude of velocity. This subroutine can be partly offloaded to devices.
  subroutine Monitor
    integer idn,id
    double precision um_global

    !Find the MPI local maximum magnitude of velocity
    um=0
    !OpenMP reduction
    !$OMP TARGET TEAMS DISTRIBUTE map(tofrom:um) private(idn,id) reduction(max:um)
    do idn=1,size_fluid-1
       id=fluid_id(idn)
       um=max(um,sqrt(u(id*dim)**2+u(id*dim+1)**2))
    enddo
    !$OMP END TARGET TEAMS DISTRIBUTE
    call MPI_BARRIER(CART_COMM,ierr)
    !MPI reduction
    call MPI_REDUCE(um,um_global,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,CART_COMM,ierr)
    !Print
    if(rank.eq.0)then
       write(*,*)"T=",iter/t_intv
       write(*,*)"Umax=",um_global/(charlength*t_intv)
       write(*,*)"-------------------"              
    endif
  end subroutine Monitor

  !---------------------------------------------
  !WriteBinary: This subourtine writes the result to a single binary file in PLT format.
  subroutine WriteBinary
    character filename*20,num*3
    integer i,j,k,ind
    integer numvars !Number of Variables, including the coordinates
    double precision uxmax,       uymax,       uzmax,       cpmax,       rhomax
    double precision uxmin,       uymin,       uzmin,       cpmin,       rhomin
    double precision uxmax_global,uymax_global,uzmax_global,cpmax_global,rhomax_global
    double precision uxmin_global,uymin_global,uzmin_global,cpmin_global,rhomin_global
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

    numvars=3+5 !3 coordinates (x,y,z), 5 fields(rho,cp,ux,uy,uz)
    
    if(rank==0)then
       offset=0
       call MPI_FILE_SEEK(file,offset,MPI_SEEK_SET,ierr)
       
       print_size=8       
       call MPI_FILE_IWRITE(file,'#!TDV112',print_size,MPI_CHARACTER,request,ierr)
       offset=offset+print_size
       
       print_size=4+numvars  +6           +1         +1         +1 +3 +2 +2 +2 +2
            !                 Len(Title)|  len(var1)| len(var2)|...        
       call MPI_FILE_IWRITE(file,(/1,0,        str2ascii('Sphere',6),0,  numvars,&
            !                        FullData| Title
            ichar('x'),0,ichar('y'),0,ichar('z'),0,str2ascii('rho',3),0,str2ascii('cp',2),0,&
            str2ascii('ux',2),0,str2ascii('uy',2),0,str2ascii('uz',2),0/)&
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
       call MPI_FILE_IWRITE(file,(/2, 2, 2, 2, 2, 2, 2, 2,       0,0,-1/)&
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
    cpmax=0.d0        
    cpmin=0.d0
    rhomax=0.5d0
    rhomin=0.5d0
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
       cpmax=max(cpmax,cp(id))      
       cpmin=min(cpmin,cp(id))
       rhomax=max(rhomax,rho(id))      
       rhomin=min(rhomin,rho(id))
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
    call MPI_REDUCE(cpmax,cpmax_global,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,CART_COMM,ierr)
    call MPI_REDUCE(cpmin,cpmin_global,1,MPI_DOUBLE_PRECISION,MPI_MIN,0,CART_COMM,ierr)
    call MPI_REDUCE(rhomax,rhomax_global,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,CART_COMM,ierr)
    call MPI_REDUCE(rhomin,rhomin_global,1,MPI_DOUBLE_PRECISION,MPI_MIN,0,CART_COMM,ierr)
    
    !Write variable limits
    if(rank==0)then
       call MPI_FILE_IWRITE(file,(/0.d0,(nx/charlength),0.d0,(ny/charlength),0.d0,(nz/charlength),&
       !                           xmin,xmax            ymin,ymax            zmin,zmax
            rhomin_global,rhomax_global,cpmin_global,cpmax_global,&
            uxmin_global/(charlength*t_intv),uxmax_global/(charlength*t_intv),&
            uymin_global/(charlength*t_intv),uymax_global/(charlength*t_intv),&
            uzmin_global/(charlength*t_intv),uzmax_global/(charlength*t_intv)/),&   
            numvars*2,       MPI_DOUBLE_PRECISION,request,ierr)
       offset=offset+numvars*2*8
    endif

    !Pass offset
    call MPI_Bcast(offset, 1, MPI_INTEGER, 0, CART_COMM,ierr)
    
    !SetWritingView    
    call MPI_TYPE_CREATE_SUBARRAY(dim,global_length,local_length,local_start,&
         MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,write_3d_type,ierr)
    call MPI_TYPE_CONTIGUOUS(numvars*1,write_3d_type,write_type,ierr)
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
    call MPI_FILE_WRITE_ALL(file,buffer,print_size,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)

    do k=0,local_length(3)-1
       do j=0,local_length(2)-1
          do i=0,local_length(1)-1      
             buffer(i+j*local_length(1)+k*local_length(1)*local_length(2)+1)=&
                  (local_start(2)+j)/charlength
          enddo
       enddo
    enddo
    call MPI_FILE_WRITE_ALL(file,buffer,print_size,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)

    do k=0,local_length(3)-1
       do j=0,local_length(2)-1
          do i=0,local_length(1)-1      
             buffer(i+j*local_length(1)+k*local_length(1)*local_length(2)+1)=&
                  (local_start(3)+k)/charlength
          enddo
       enddo
    enddo
    call MPI_FILE_WRITE_ALL(file,buffer,print_size,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)

    do k=0,local_length(3)-1
       do j=0,local_length(2)-1
          do i=0,local_length(1)-1      
             buffer(i+j*local_length(1)+k*local_length(1)*local_length(2)+1)=&
                  rho((i+ghost)+dy*(j+ghost)+dz*(k+ghost))             
          enddo
       enddo
    enddo
    call MPI_FILE_WRITE_ALL(file,buffer,print_size,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
    
    do k=0,local_length(3)-1
       do j=0,local_length(2)-1
          do i=0,local_length(1)-1      
             buffer(i+j*local_length(1)+k*local_length(1)*local_length(2)+1)=&
                  cp((i+ghost)+dy*(j+ghost)+dz*(k+ghost))
          enddo
       enddo
    enddo
    call MPI_FILE_WRITE_ALL(file,buffer,print_size,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)

    do ind=0,dim-1
       do k=0,local_length(3)-1
          do j=0,local_length(2)-1
             do i=0,local_length(1)-1      
                buffer(i+j*local_length(1)+k*local_length(1)*local_length(2)+1)=&
                     u(((i+ghost)+dy*(j+ghost)+dz*(k+ghost))*dim+ind)&
                     /(charlength*t_intv)!Unit conversion
             enddo
          enddo
       enddo
       call MPI_FILE_WRITE_ALL(file,buffer,print_size,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
    enddo

    !Deallocate writing buffer
    DEALLOCATE(buffer)    
     
    call MPI_BARRIER(CART_COMM,ierr)
    call MPI_FILE_CLOSE(file,ierr)

  contains
    !A function that transforms a string to ASCII code
    !s: String
    !n: Size of string
    !o: Generated ASCII code
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

  !------------------------------------------------------
  !Laplacian
  subroutine Laplacian
    integer id,idn,iq
    double precision diff
    call PassD(rho)
    do idn=0,size_fluid-1
       id=fluid_id(idn)
       lap(id)=0.d0
       do iq=0,nq-1
          diff=rho(id+e(iq*dim)+dy*e(iq*dim+1)+dz*e(iq*dim+2))&
               -2.d0*rho(id)&
               +rho(id-e(iq*dim)-dy*e(iq*dim+1)-dz*e(iq*dim+2))
          lap(id)=lap(id)+t(iq)*diff
       enddo
       lap(id)=lap(id)*3.d0
    enddo
  endsubroutine Laplacian

  !------------------------------------------------------
  !Central Gradient
  subroutine GradC(dval,val)
    double precision,dimension(0:array_size-1),intent(in)::val
    double precision,dimension(0:array_size*dim-1),intent(out)::dval
    integer id,idn,iq,ind
    double precision diff
    do idn=0,size_fluid-1
       id=fluid_id(idn)
       do ind=0,dim-1
          dval(id*dim+ind)=0.d0
       enddo
       do iq=0,nq-1
          diff=0.5d0*(val(id+e(iq*dim)+dy*e(iq*dim+1)+dz*e(iq*dim+2))&
               -val(id-e(iq*dim)-dy*e(iq*dim+1)-dz*e(iq*dim+2)))
          do ind=0,dim-1
             dval(id*dim+ind)=dval(id*dim+ind)+t(iq)*diff*e(iq*dim+ind)*3.d0
          enddo
       enddo
       
    enddo
  endsubroutine GradC

  !------------------------------------------------------
  !Mixed Gradient
  subroutine GradM(dval,val)
    double precision,dimension(0:array_size-1),intent(in)::val
    double precision,dimension(0:array_size*dim-1),intent(out)::dval
    integer id,idn,iq,ind
    double precision diff
    do idn=0,size_fluid-1
       id=fluid_id(idn)
       do ind=0,dim-1
          dval(id*dim+ind)=0.d0
       enddo
       do iq=0,nq-1
          diff=0.25d0*(5.d0*val(id+e(iq*dim)+dy*e(iq*dim+1)+dz*e(iq*dim+2))&
               -3.d0*val(id)&
               -val(id-e(iq*dim)-dy*e(iq*dim+1)-dz*e(iq*dim+2))&
               -val(id+2*e(iq*dim)+2*dy*e(iq*dim+1)+2*dz*e(iq*dim+2)))
          do ind=0,dim-1
             dval(id*dim+ind)=dval(id*dim+ind)+t(iq)*diff*e(iq*dim+ind)*3.d0
          enddo
       enddo
       
    enddo

    
  endsubroutine GradM
  
endmodule lbm
