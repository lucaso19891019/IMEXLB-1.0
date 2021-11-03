!Geng Liu, April 2nd, 2021.
module cart_mpi
  !This module generates a 3D cartesian MPI communicator, and defines the arrays distributed to each MPI processor. The module also defines the data types for sending and receiving subarrays.
  use mpi
  implicit none
  !dim: dimension of the problem
  !nq: number of lattice velocities
  !ghost: thickness of ghost nodes
  !num_neighbots: number of neighbors for one processor
  integer,parameter:: dim=3,nq=27,ghost=3,num_neighbors=26
  
  !nx: number of points in x direction
  !ny: number of points in y direction
  !nz: number of points in z direction
  integer nx,ny,nz

  !ierr: MPI error signal
  !num_procs: number of MPI processors
  !rank: 1D ID of this MPI processor
  integer ierr, num_procs, rank

  !CART_COMM: newly generated MPI cartesian communicator
  integer CART_COMM

  !cart_coords: cartesian coordinates of current MPI processor
  !cart_num_procs: number of MPI processors in different cartesian directions
  integer,dimension(dim)::cart_coords,cart_num_procs
  !cart_periods: flags for periodic communication in each direction
  logical,dimension(dim)::cart_periods
  !global_length: directional length of whole domain
  !local_length: directional length of array distributed to current MPI processor
  !local_start: local starting point's position in the whole domain
  !local_end: local ending point's position in the whole domain
  integer,dimension(dim)::global_length,local_length,local_start,local_end
  !size of scalor field array in current processor (reshaped to 1D)
  integer array_size

  !neighbor_procs: ID of each neighbor MPI processor
  integer,dimension(0:num_neighbors)::neighbor_procs
  !Subarray datatypes for receiving data from each neighbor processor
  !Note: INT-integer, F-PDF (or variable of same dimension), U-Velocity (or variable of same dimension), DPR-double precision, LOG-logical
  integer,dimension(0:num_neighbors)::INT_RECV_TYPE,F_RECV_TYPE,U_RECV_TYPE,DPR_RECV_TYPE,LOG_RECV_TYPE
  !Subarray datatypes for sending data to each neighbor processor
  integer,dimension(0:num_neighbors)::INT_SEND_TYPE,F_SEND_TYPE,U_SEND_TYPE,DPR_SEND_TYPE,LOG_SEND_TYPE

contains
  !MPISetup: By calling other subroutines, this subroutine generates a cartesian MPI communicator, and defines the arrays distributed to each MPI processor as well as the data types for sending and receiving subarrays.
  subroutine MPISetup
     
    call SetCartesian
    call SetLocalLength
    call SetNeighbors

  endsubroutine MPISetup
  
  !------------------------------------------
  !SetCartesian: This subroutine generates a cartesian MPI communicator. It optimizes the number of processors in each direction.
  subroutine SetCartesian   
    double precision max_e, config_e
    integer px,py,pyz,pz
    logical reorder_flag
    integer i,j
    
    global_length(1)=nx+1
    global_length(2)=ny+1
    global_length(3)=nz+1
    
    !Optimization of the number of processors in each direction.
    max_e=0.d0
    do i=1,num_procs
       if(mod(num_procs,i).eq.0)then
          px=i
          pyz=num_procs/i
          do j=1,pyz
             if(mod(pyz,j).eq.0)then
                py=j
                pz=pyz/j
                if(.not.((global_length(3).eq.2).and.(pz>1)))then
                   config_e=1.d0/(2*(global_length(2)*(px-1)*global_length(3)/py/pz&
                        &+global_length(1)*(py-1)*global_length(3)/px/pz&
                        &+global_length(1)*global_length(2)*(pz-1)/px/py))
                   if(config_e.ge.max_e)then
                      max_e=config_e
                      cart_num_procs(1)=px
                      cart_num_procs(2)=py
                      cart_num_procs(3)=pz
                   endif
                endif
             endif
          enddo
       endif
    enddo

    do i=1,dim
       cart_periods(i)=.true.
    enddo
    !reorder_flag: whether to reorder the MPI processors or not after generating the cartesian communicator.
    reorder_flag=.true.
    !Generate cartesian communicator
    call MPI_CART_CREATE(MPI_COMM_WORLD,dim,cart_num_procs,cart_periods,reorder_flag,CART_COMM,ierr)
    !Get reordered ID of current MPI processor
    call MPI_COMM_RANK(CART_COMM,rank,ierr)
    !Get current processor's coordinates
    call MPI_CART_COORDS(CART_COMM,rank,dim,cart_coords,ierr)

    call MPI_BARRIER(CART_COMM,ierr)    
  endsubroutine SetCartesian
  
  !-------------------------------------------------------------
  !SetLocalLength: Based on the communicator generated in SetCartesian subroutine, this subroutine generates the local array information, including the local array length in each direction, and the starting and ending points' positions in the whole domain.
  subroutine SetLocalLength
    integer i,type1,ind
    integer,dimension(dim,num_procs)::length_public
    integer,dimension(dim)::this_cart_coords
    integer,dimension(:,:,:,:),allocatable::local_length_db
    !Get the directional length of local array
    do i=1,dim
       local_length(i)=global_length(i)/cart_num_procs(i)+(1-sign(1,(cart_coords(i)-mod(global_length(i),cart_num_procs(i)))))/2
       if(local_length(i)<ghost)then
          print*, 'Dimension = ',i,' in Cartesian proc id =',rank,' smaller than GhostCells'
          call MPI_ABORT(MPI_COMM_WORLD,0,ierr)
       endif
    enddo

    !Collect the local length information of each MPI processor in a temporary database.
    allocate(local_length_db(dim,0:cart_num_procs(1)-1,0:cart_num_procs(2)-1,0:cart_num_procs(3)-1))
    
    call MPI_TYPE_VECTOR(dim,1,1,MPI_INTEGER,type1,ierr)
    call MPI_TYPE_COMMIT(type1,ierr)
    call MPI_ALLGATHER(local_length(:),1,type1,length_public,1,type1,CART_COMM,ierr)
    call MPI_TYPE_FREE(type1,ierr)
    do i=1,num_procs
       call MPI_CART_COORDS(CART_COMM,i-1,dim,this_cart_coords,ierr)
       local_length_db(:,this_cart_coords(1),this_cart_coords(2),this_cart_coords(3))=length_public(:,i)
    enddo

    !Evaluate the starting and ending points' positions
    local_start(1)=sum(local_length_db(1,0:cart_coords(1)-1,cart_coords(2),cart_coords(3)))
    local_start(2)=sum(local_length_db(2,cart_coords(1),0:cart_coords(2)-1,cart_coords(3)))
    local_start(3)=sum(local_length_db(3,cart_coords(1),cart_coords(2),0:cart_coords(3)-1))

    local_end(1)=sum(local_length_db(1,0:cart_coords(1),cart_coords(2),cart_coords(3)))-1
    local_end(2)=sum(local_length_db(2,cart_coords(1),0:cart_coords(2),cart_coords(3)))-1
    local_end(3)=sum(local_length_db(3,cart_coords(1),cart_coords(2),0:cart_coords(3)))-1

    deallocate(local_length_db)

    array_size=1;
    do ind=1,dim
       array_size=array_size*(local_length(ind)+2*ghost)
    enddo
    
  endsubroutine SetLocalLength

  !--------------------------------------------
  !SetNeighbors: This subroutine defines the neighbors of this MPI processor, and generates sending and receiving subarray datatypes for communication with these neighbors.
  subroutine SetNeighbors
    !neighbor_local_rank: rank of neighbor reordered in the world of neighbors of current MPI processor (local relative rank)
    !neighbor_pos_p(scalar): local relative coordinate in a certain direction.
    !neighbor_pos(vector): local relative coordinate vector.
    integer i,j,k,ind, neighbor_local_rank,neighbor_pos_p
    integer,dimension(dim)::neighbor_pos

    !neighbor_cart_coords: global coordinates of neighbor
    integer,dimension(dim)::neighbor_cart_coords
    !Subarray information for communication: f-PDF, u-Velocity, no prefix-scalor 
    integer,dimension(dim+1)::fsizes,fsubsize,frecv_start,fsend_start,usizes,usubsize,urecv_start,usend_start
    integer,dimension(dim)::sizes,subsize,recv_start,send_start
    
    do i=-1,1
       do j=-1,1
          do k=-1,1
             neighbor_pos(1)=i
             neighbor_pos(2)=j
             neighbor_pos(3)=k
             neighbor_local_rank=(k+1)*9+(i+1)*3+(j+1)

          do ind=1,dim
             neighbor_cart_coords(ind)=cart_coords(ind)+neighbor_pos(ind)
          enddo

          call MPI_CART_RANK(CART_COMM,neighbor_cart_coords,neighbor_procs(neighbor_local_rank),ierr)

          fsizes(1)=nq
          fsubsize(1)=nq
          frecv_start(1)=0
          fsend_start(1)=0

          usizes(1)=dim
          usubsize(1)=dim
          urecv_start(1)=0
          usend_start(1)=0
          
          do ind=1,dim
             neighbor_pos_p=neighbor_pos(ind)
             
             sizes(ind)=local_length(ind)+2*ghost             
             subsize(ind)=local_length(ind)*(1-abs(neighbor_pos_p))+ghost*abs(neighbor_pos_p)
             recv_start(ind)=neighbor_pos_p*(1-neighbor_pos_p)/2*ghost+neighbor_pos_p*(1+neighbor_pos_p)/2*local_length(ind)+ghost
             send_start(ind)=neighbor_pos_p*(1+neighbor_pos_p)/2*(local_length(ind)-ghost)+ghost

             fsizes(ind+1)=sizes(ind)
             fsubsize(ind+1)=subsize(ind)
             frecv_start(ind+1)=recv_start(ind)
             fsend_start(ind+1)=send_start(ind)

             usizes(ind+1)=sizes(ind)
             usubsize(ind+1)=subsize(ind)
             urecv_start(ind+1)=recv_start(ind)
             usend_start(ind+1)=send_start(ind)
          enddo

          !Generating subarray datatypes
          call MPI_TYPE_CREATE_SUBARRAY(dim,sizes,subsize,recv_start,MPI_ORDER_FORTRAN,&
               &MPI_INTEGER,INT_RECV_TYPE(neighbor_local_rank),ierr)
          call MPI_TYPE_COMMIT(INT_RECV_TYPE(neighbor_local_rank),ierr)
          call MPI_TYPE_CREATE_SUBARRAY(dim,sizes,subsize,send_start,MPI_ORDER_FORTRAN,&
               &MPI_INTEGER,INT_SEND_TYPE(neighbor_local_rank),ierr)
          call MPI_TYPE_COMMIT(INT_SEND_TYPE(neighbor_local_rank),ierr)

          call MPI_TYPE_CREATE_SUBARRAY(dim+1,fsizes,fsubsize,frecv_start,MPI_ORDER_FORTRAN,&
               &MPI_DOUBLE_PRECISION,F_RECV_TYPE(neighbor_local_rank),ierr)
          call MPI_TYPE_COMMIT(F_RECV_TYPE(neighbor_local_rank),ierr)
          call MPI_TYPE_CREATE_SUBARRAY(dim+1,fsizes,fsubsize,fsend_start,MPI_ORDER_FORTRAN,&
               &MPI_DOUBLE_PRECISION,F_SEND_TYPE(neighbor_local_rank),ierr)
          call MPI_TYPE_COMMIT(F_SEND_TYPE(neighbor_local_rank),ierr)

          call MPI_TYPE_CREATE_SUBARRAY(dim+1,usizes,usubsize,urecv_start,MPI_ORDER_FORTRAN,&
               &MPI_DOUBLE_PRECISION,U_RECV_TYPE(neighbor_local_rank),ierr)
          call MPI_TYPE_COMMIT(U_RECV_TYPE(neighbor_local_rank),ierr)
          call MPI_TYPE_CREATE_SUBARRAY(dim+1,usizes,usubsize,usend_start,MPI_ORDER_FORTRAN,&
               &MPI_DOUBLE_PRECISION,U_SEND_TYPE(neighbor_local_rank),ierr)
          call MPI_TYPE_COMMIT(U_SEND_TYPE(neighbor_local_rank),ierr)

          call MPI_TYPE_CREATE_SUBARRAY(dim,sizes,subsize,recv_start,MPI_ORDER_FORTRAN,&
               &MPI_DOUBLE_PRECISION,DPR_RECV_TYPE(neighbor_local_rank),ierr)
          call MPI_TYPE_COMMIT(DPR_RECV_TYPE(neighbor_local_rank),ierr)
          call MPI_TYPE_CREATE_SUBARRAY(dim,sizes,subsize,send_start,MPI_ORDER_FORTRAN,&
               &MPI_DOUBLE_PRECISION,DPR_SEND_TYPE(neighbor_local_rank),ierr)
          call MPI_TYPE_COMMIT(DPR_SEND_TYPE(neighbor_local_rank),ierr)

          call MPI_TYPE_CREATE_SUBARRAY(dim,sizes,subsize,recv_start,MPI_ORDER_FORTRAN,&
               &MPI_LOGICAL,LOG_RECV_TYPE(neighbor_local_rank),ierr)
          call MPI_TYPE_COMMIT(LOG_RECV_TYPE(neighbor_local_rank),ierr)
          call MPI_TYPE_CREATE_SUBARRAY(dim,sizes,subsize,send_start,MPI_ORDER_FORTRAN,&
               &MPI_LOGICAL,LOG_SEND_TYPE(neighbor_local_rank),ierr)
          call MPI_TYPE_COMMIT(LOG_SEND_TYPE(neighbor_local_rank),ierr)          
       enddo
    enddo
    enddo
    
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    
  endsubroutine SetNeighbors
  
  !-------------------------------------------------
  !PassF: This subroutine is called when communication of the Array (PDF or variable of the same dimension) is needed. The communication is based on the sending and receiving subarray datatypes generated in SetNeighbors
subroutine PassF(Array)

  double precision,dimension(0:nq*array_size-1)::array
  !neighbor_local_rank:local(relative) rank of neighbor
  !neighbor_oppos_rank: local(relative) rank of neighbor opposite to neighbor_local_rank
  integer neighbor_local_rank,neighbor_oppos_rank
  
   !Requests and communication status.
   integer,dimension(0:num_neighbors,2):: req
   integer,dimension(MPI_STATUS_SIZE,0:num_neighbors,2) :: communication_status

   call MPI_BARRIER(CART_COMM,ierr)

   !Communication
   do neighbor_local_rank=0,num_neighbors/2-1

      neighbor_oppos_rank=num_neighbors-neighbor_local_rank

      call MPI_ISEND(Array,1,F_SEND_TYPE(neighbor_local_rank),neighbor_procs(neighbor_local_rank),&
           &neighbor_local_rank,CART_COMM,req(neighbor_local_rank,1),ierr)
      call MPI_IRECV(Array,1,F_RECV_TYPE(neighbor_oppos_rank),neighbor_procs(neighbor_oppos_rank),&
           &neighbor_local_rank,CART_COMM,req(neighbor_local_rank,2),ierr)
      call MPI_WAIT(req(neighbor_local_rank,2),communication_status(:,neighbor_local_rank,2),ierr)
      call MPI_WAIT(req(neighbor_local_rank,1),communication_status(:,neighbor_local_rank,1),ierr)

      call MPI_BARRIER(CART_COMM,ierr)

      call MPI_ISEND(Array,1,F_SEND_TYPE(neighbor_oppos_rank),neighbor_procs(neighbor_oppos_rank),&
           &neighbor_oppos_rank,CART_COMM,req(neighbor_oppos_rank,1),ierr) 
      call MPI_IRECV(Array,1,F_RECV_TYPE(neighbor_local_rank),neighbor_procs(neighbor_local_rank),&
           &neighbor_oppos_rank,CART_COMM,req(neighbor_oppos_rank,2),ierr)
      call MPI_WAIT(req(neighbor_oppos_rank,2),communication_status(:,neighbor_oppos_rank,2),ierr)
      call MPI_WAIT(req(neighbor_oppos_rank,1),communication_status(:,neighbor_oppos_rank,1),ierr)

      call MPI_BARRIER(CART_COMM,ierr)

   enddo
 endsubroutine PassF

 !-------------------------------------------
 !PassInt: This subroutine is communication for integer arrays, communication for arrays of other types can be similar to PassF and PassInt (array type and sending and receiving types should be replaced).
 subroutine PassInt(Array)
  integer,dimension(0:array_size-1)::array
  !neighbor_local_rank:local(relative) rank of neighbor
  !neighbor_oppos_rank: local(relative) rank of neighbor opposite to neighbor_local_rank
  integer neighbor_local_rank,neighbor_oppos_rank
  
   !Requests and communication status.
   integer,dimension(0:num_neighbors,2):: req
   integer,dimension(MPI_STATUS_SIZE,0:num_neighbors,2) :: communication_status

   call MPI_BARRIER(CART_COMM,ierr)

   !Communication
   do neighbor_local_rank=0,num_neighbors/2-1

      neighbor_oppos_rank=num_neighbors-neighbor_local_rank

      call MPI_ISEND(Array,1,INT_SEND_TYPE(neighbor_local_rank),neighbor_procs(neighbor_local_rank),&
           &neighbor_local_rank,CART_COMM,req(neighbor_local_rank,1),ierr)
      call MPI_IRECV(Array,1,INT_RECV_TYPE(neighbor_oppos_rank),neighbor_procs(neighbor_oppos_rank),&
           &neighbor_local_rank,CART_COMM,req(neighbor_local_rank,2),ierr)
      call MPI_WAIT(req(neighbor_local_rank,2),communication_status(:,neighbor_local_rank,2),ierr)
      call MPI_WAIT(req(neighbor_local_rank,1),communication_status(:,neighbor_local_rank,1),ierr)

      call MPI_BARRIER(CART_COMM,ierr)

      call MPI_ISEND(Array,1,INT_SEND_TYPE(neighbor_oppos_rank),neighbor_procs(neighbor_oppos_rank),&
           &neighbor_oppos_rank,CART_COMM,req(neighbor_oppos_rank,1),ierr) 
      call MPI_IRECV(Array,1,INT_RECV_TYPE(neighbor_local_rank),neighbor_procs(neighbor_local_rank),&
           &neighbor_oppos_rank,CART_COMM,req(neighbor_oppos_rank,2),ierr)
      call MPI_WAIT(req(neighbor_oppos_rank,2),communication_status(:,neighbor_oppos_rank,2),ierr)
      call MPI_WAIT(req(neighbor_oppos_rank,1),communication_status(:,neighbor_oppos_rank,1),ierr)

      call MPI_BARRIER(CART_COMM,ierr)

   enddo
endsubroutine PassInt
 
endmodule cart_mpi
