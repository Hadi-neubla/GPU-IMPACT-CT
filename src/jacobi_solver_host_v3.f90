!**************************************************************************************************                                           
!* by Hadi Zolfaghari, ARTORG Center, University of Bern, (zolfaghari.haadi@gmail.com)            *      
!* October 2015 - March 2020                                                                      *            
!**************************************************************************************************
subroutine device_jacobi_solver(grid_id,omega,SS1,SS2,SS3,NN1,NN2,NN3,rel,bb)

USE mod_dims
USE mod_vars
USE mod_vars_GPU
USE mod_exchange
USE mod_diff
USE mod_laplace
USE mod_inout
USE mod_cudafor



use iso_c_binding

! Module where we define CUDA functions to use them from Fortran

implicit none


   integer, target, intent(in)  ::  grid_id  ! the grid level in multigrid --> corresponding to 'g = 1,2,3 ..' in IMPACT
   
   real(8), target, intent(in)     ::  omega

   INTEGER, target, INTENT(in)    ::  SS1      !< lower index bound dimension 1
   INTEGER, target, INTENT(in)    ::  SS2      !< lower index bound dimension 2
   INTEGER, target, INTENT(in)    ::  SS3      !< lower index bound dimension 3
 
   INTEGER, target, INTENT(in)    ::  NN1      !< upper index bound dimension 1
   INTEGER, target, INTENT(in)    ::  NN2      !< upper index bound dimension 2
   INTEGER, target, INTENT(in)    ::  NN3      !< upper index bound dimension 3



   REAL(8)   , target, INTENT(inout) ::  rel (b1L:(NN(1,grid_id)+b1U),b2L:(NN(2,grid_id)+b2U),b3L:(NN(3,grid_id)+b3U))!< ? 
   REAL(8)  , target, INTENT(inout) ::  bb (b1L:(NN(1,grid_id)+b1U),b2L:(NN(2,grid_id)+b2U),b3L:(NN(3,grid_id)+b3U))!< ? 
!   REAL(8)   , target, INTENT(  out) ::  comp(b1L:(NN(1,grid_id)+b1U),b2L:(NN(2,grid_id)+b2U),b3L:(NN(3,grid_id)+b3U))!< gradient operator


   real(8), target    :: err_tot, err_tot_global,err_tot_old, rel_err_tot
   integer, target :: mode_j, iter_level

   REAL(8)   , target ::  rel_old(b1L:(NN(1,grid_id)+b1U),b2L:(NN(2,grid_id)+b2U),b3L:(NN(3,grid_id)+b3U))!< ? 
   REAL(8)   , target ::  GPUres(1:N1p,1:N2p,1:N3p)!< ? 

!   REAL(8)   , target ::  west_ghost(1:N2-1,1:N3-1)
!   REAL(8)   , target ::  east_ghost(1:N2-1,1:N3-1)
!   REAL(8)   , target ::  front_ghost(1:N1-1,1:N2-1)
!   REAL(8)   , target ::  rear_ghost(1:N1-1,1:N2-1)
!   REAL(8)   , target ::  lower_ghost(1:N1-1,1:N3-1)
!   REAL(8)   , target ::  upper_ghost(1:N1-1,1:N3-1)
!
!
!   REAL(8)   , target ::  west_chunk(1:N2-1,1:N3-1)
!   REAL(8)   , target ::  east_chunk(1:N2-1,1:N3-1)
!   REAL(8)   , target ::  front_chunk(1:N1-1,1:N2-1)
!   REAL(8)   , target ::  rear_chunk(1:N1-1,1:N2-1)
!   REAL(8)   , target ::  lower_chunk(1:N1-1,1:N3-1)
!   REAL(8)   , target ::  upper_chunk(1:N1-1,1:N3-1)





!   REAL(8)   , target ::  slide_south(b1L:(NN(1,grid_id)+b1U),b2L:(NN(2,grid_id)+b2U),b3L:(NN(3,grid_id)+b3U))!< ? 
!   REAL(8)   , target ::  slide_north(b1L:(NN(1,grid_id)+b1U),b2L:(NN(2,grid_id)+b2U),b3L:(NN(3,grid_id)+b3U))!< ? 



  ! Declare the size of rel and other arguments passing to jacobi relaxation -- ref : mod_laplace

   integer   ::   rel_size(3)

   integer   ::   N1R, N2R, N3R

  ! Declare sincos kernel arguments with "target" -- for c_loc
  integer, target     :: nx, ny, nz ! the rel, bb and comp size components



!   integer   ::   cdg1_size(3)   
!   integer   ::   cdg2_size(3)
!   integer   ::   cdg3_size(3)

!  integer, target     :: nx_cdg1, ny_cdg1, nz_cdg1 ! the grad and phi size components
!  integer, target     :: nx_cdg2, ny_cdg2, nz_cdg2 ! the grad and phi size components
!  integer, target     :: nx_cdg3, ny_cdg3, nz_cdg3 ! the grad and phi size components


  !type(c_ptr) :: rel_dev, bb_dev, comp_dev, cdg1_dev,  cdg2_dev, cdg3_dev
!  integer :: sum_host


!  integer(c_int) :: istat
!
!  integer, parameter :: szblock = 1
!
!  ! CUDA constants
!  integer(c_int), parameter :: cudaMemcpyHostToDevice = 1
!  integer(c_int), parameter :: cudaMemcpyDeviceToHost = 2
!
!  integer(c_size_t) :: szreal, szint, szptr, offset = 0

  type(dim3) :: blocks, threads


!  type(cudaStream_t) :: stream1, stream2, stream3

  integer :: i, j, k, step
  real(8), volatile :: start, finish

  integer  ::   iter

  real(8)     ::  norm_inf, norm_inf2, norm_inf_global, norm_inf_global2, norm_zero
  real(8)     :: norm_rhs, norm_rhs_global, norm_inf_rel
  real(8)     ::  t1, t2, t3
  include 'mpif.h'

!! Aplying the exchange !!TODO this part is CPU only and a subject for shared memory optimization

!   CALL exchange(dimension_id,0,phi) 


!  szreal = sizeof(rel(1, 1, 1))
!  szint = sizeof(nx)
!  szptr = sizeof(rel_dev)
 
 
!    b1L, b1U, NN(1,1) they are -3 3 65 for 64^3 problem
 ! We repeat GPU logic twice to show that the first kernel invocation includes
  ! the time of GPU device initialization, which is done one time per
  ! application run.

!    call cpu_time(start)



    ! Get the size of arrays passing to jacobi relaxation -- > ref : mod_laplace
!    rel_size = shape(rel)
!
!    nx = rel_size(1)  
!    ny = rel_size(2) 
!    nz = rel_size(3) 
!    


!    cdg1_size = shape(cdg1)
!    cdg2_size = shape(cdg2)
!    cdg3_size = shape(cdg3)
!
!
!    nx_cdg1 = cdg1_size(1)
!    ny_cdg1 = cdg1_size(2)
!    nz_cdg1 = cdg1_size(3)
!
!    nx_cdg2 = cdg2_size(1)
!    ny_cdg2 = cdg2_size(2)
!    nz_cdg2 = cdg2_size(3)
!
!    nx_cdg3 = cdg3_size(1)
!    ny_cdg3 = cdg3_size(2)
!    nz_cdg3 = cdg3_size(3)



    S1R  = SNB(1,1,grid_id)
    S2R  = SNB(1,2,grid_id)
    S3R  = SNB(1,3,grid_id)

    N1R  = SNB(2,1,grid_id)
    N2R  = SNB(2,2,grid_id)
    N3R  = SNB(2,3,grid_id)

!    print*, S1R, N1R, 'SNB'
!
!    S1R  = SNF(1,1,grid_id)
!    N1R  = SNF(2,1,grid_id)
!    
!    print*, S1R, N1R, 'SNF'

!    print*, nx, ny, nz, nx_cdg1, nx_cdg2, nx_cdg3, ny_cdg1, ny_cdg2, ny_cdg3, nz_cdg1, nz_cdg2, nz_cdg3, 'nx, ny, nz, nx_cdg1, nx_cdg2, nx_cdg3, ny_cdg1, ny_cdg2, ny_cdg3, nz_cdg1, nz_cdg2, nz_cdg3'
    ! Allocate memory on GPU.

!if (gpu_step == 1) then

!    istat = cudaMalloc(cdg1_dev, nx_cdg1 * ny_cdg1 * nz_cdg1 * szreal)
!    call cuda_err_check(istat)
!
!    istat = cudaMalloc(cdg2_dev, nx_cdg2 * ny_cdg2 * nz_cdg2 * szreal)
!    call cuda_err_check(istat)
!
!    istat = cudaMalloc(cdg3_dev, nx_cdg3 * ny_cdg3 * nz_cdg3 * szreal)
!    call cuda_err_check(istat)  
! 
!
!

!    istat = cudaMalloc(rel_dev, nx_dev * ny_dev * nz_dev * szreal)
!    call cuda_err_check(istat)
!    istat = cudaMalloc(bb_dev, nx_dev * ny_dev * nz_dev *  szreal)
!    call cuda_err_check(istat)



!!    istat = cudaMalloc(comp_dev, nx * ny * nz *  szreal)
!!    call cuda_err_check(istat)
!
!
!    ! defining ghost_cells to be copied to GPU
!    istat = cudaMalloc(west_ghost_dev, (N2-1) * (N3-1) * szreal)
!    call cuda_err_check(istat)
!    istat = cudaMalloc(east_ghost_dev, (N2-1) * (N3-1) * szreal)
!    call cuda_err_check(istat)
!    istat = cudaMalloc(front_ghost_dev, (N1-1) * (N2-1) * szreal)
!    call cuda_err_check(istat)
!    istat = cudaMalloc(rear_ghost_dev, (N1-1) * (N2-1) * szreal)
!    call cuda_err_check(istat)
!    istat = cudaMalloc(upper_ghost_dev, (N1-1) * (N3-1) * szreal)
!    call cuda_err_check(istat)
!    istat = cudaMalloc(lower_ghost_dev, (N1-1) * (N3-1) * szreal)
!    call cuda_err_check(istat)
!
!	! XXX Several blocks of code below are equivalent to sum_kernel<<<...>>> in C.
!
!    ! Create CUDA compute grid configuration.
!    
!     
!!     print*, 'grid size',  (nx - 8)/8, (ny - 8)/8, (nz - 8)/8
!!    blocks = dim3(1, 1, 1)
!!    threads = dim3(1, 1, 1)
!    ! Copy input data from CPU to GPU memory.
!    istat = cudaMemcpy(cdg1_dev, c_loc(cdg1), nx_cdg1 * ny_cdg1 * nz_cdg1 * szreal, cudaMemcpyHostToDevice)
!    call cuda_err_check(istat)
!    istat = cudaMemcpy(cdg2_dev, c_loc(cdg2), nx_cdg2 * ny_cdg2 * nz_cdg2 * szreal, cudaMemcpyHostToDevice)
!    call cuda_err_check(istat)
!    istat = cudaMemcpy(cdg3_dev, c_loc(cdg3), nx_cdg3 * ny_cdg3 * nz_cdg3 * szreal, cudaMemcpyHostToDevice)
!    call cuda_err_check(istat)
    
!    endif


    istat = cudaMemcpy(bb_dev, c_loc(bb), nx_dev * ny_dev * nz_dev * szreal, cudaMemcpyHostToDevice)
    call cuda_err_check(istat)

    istat = cudaMemcpy(rel_dev, c_loc(rel), nx_dev * ny_dev * nz_dev * szreal, cudaMemcpyHostToDevice)
    call cuda_err_check(istat)

    !----------- initialize the ghosts
    west_ghost  = rel(S1R-1,1:N2-1,1:N3-1)
    east_ghost  = rel(N1R+1,1:N2-1,1:N3-1)
    lower_ghost = rel(1:N1-1,S2R-1,1:N3-1)
    upper_ghost = rel(1:N1-1,N2R+1,1:N3-1) 
    rear_ghost  = rel(1:N1-1,1:N2-1,S3R-1)
    front_ghost = rel(1:N1-1,1:N2-1,N3R+1)

    west_chunk  = rel(S1R,1:N2-1,1:N3-1)
    east_chunk  = rel(N1R,1:N2-1,1:N3-1)
    lower_chunk = rel(1:N1-1,S2R,1:N3-1)
    upper_chunk = rel(1:N1-1,N2R,1:N3-1) 
    rear_chunk  = rel(1:N1-1,1:N2-1,S3R)
    front_chunk = rel(1:N1-1,1:N2-1,N3R)


#ifdef CUDAVERSION9
    blocks = dim3((nx_dev - 8)/8, (ny_dev - 8)/8, (nz_dev - 8)/8)
    if (mod(nx_Dev, szblock) .ne. 0) blocks%x = blocks%x + 1
    threads = dim3((nx_dev - 8)/blocks%x, (ny_dev - 8)/blocks%y, (nz_dev - 8)/blocks%z)
#endif


    iter=0
    mode_j=10
    iter_level=0
    norm_zero =0.
    iterate: do

    ! do the MPI exchange before copying to GPU

    iter  = iter+1


    ! print*, iter


!----------------------------exchange----------------------------
!----------------------------------------------------------------

!     CALL exchange_relax_GPU(1,0,0,0,0,.FALSE.,rel)

      
     CALL exchange_relax_GPU2

!----------------------------------------------------------------



    if ( (mod(iter,mode_j)==1))  then  !--------------------------------------------------------------------------------------------------------

!    iter_level=iter_level+mode_j/10
!    if (iter_level == 10) then
!        iter_level=0
!        mode_j=10*mode_j
!    endif 
    


    istat = cudaMemcpy(c_loc(rel), rel_dev, nx_dev * ny_dev * nz_dev * szreal, cudaMemcpyDeviceToHost)
    call cuda_err_check(istat)


     CALL exchange_relax_GPU(1,0,0,0,0,.FALSE.,rel)


      norm_inf = 0.
      norm_inf2= 0.
      !norm_rhs=0.
      !norm_rhs_global=0.
      err_tot = 0.
      err_tot_global = 0. 
       !print*, N1p
       DO k = SS3, NN3
          DO j = SS2, NN2
             DO i = SS1, NN1
!                GPUres(i,j,k) = bb(i,j,k)                                     &
!                            &      - cdg1(-1,i,1)*rel(i-1,j,k) - cdg1(1,i,1)*rel(i+1,j,k)     &
!                            &      - cdg2(-1,j,1)*rel(i,j-1,k) - cdg2(1,j,1)*rel(i,j+1,k)     &
!                            &      - cdg3(-1,k,1)*rel(i,j,k-1) - cdg3(1,k,1)*rel(i,j,k+1)    &
!                            &      - (cdg1(0,i,1) + cdg2(0,j,1) + cdg3(0,k,1))*rel(i,j,k)
!                norm_inf = MAX(ABS(GPUres(i,j,k)),norm_inf)
!                err_tot = err_tot + abs((rel(i,j,k)-rel_old(i,j,k)))/((NN3-SS3+1)*(NN2-SS2+1)*(NN1-SS1+1))
                err_tot = err_tot + abs((rel(i,j,k)-rel_old(i,j,k)))
!              ! norm_rhs = MAX(ABS(bb(i,j,k)),norm_rhs)
             END DO
          END DO
       END DO
       !err_tot = sqrt(sum((rel-rel_old)**2))
       !err_tot = sqrt(((rel(5,5,5)-rel_old(5,5,5))**2))

    !--------------------------------------------------------------------------------------------------------
    !--For infinity norm of the residual
!       CALL MPI_ALLREDUCE(norm_inf,norm_inf_global,1,MPI_REAL8,MPI_MAX,COMM_CART,merror) 
!       norm_inf = norm_inf_global
    !--For infinity norm of the 
       ! for automatic relative norm
       !CALL MPI_ALLREDUCE(norm_rhs,norm_rhs_global,1,MPI_REAL8,MPI_MAX,COMM_CART,merror) 
       !norm_rhs = norm_rhs_global
       !norm_inf_rel=norm_inf/norm_rhs
   ! --For total error
       CALL MPI_ALLREDUCE(err_tot,err_tot_global,1,MPI_REAL8,MPI_MAX,COMM_CART,merror) 
       err_tot = err_tot_global

       rel_err_tot = abs(log(err_tot_old/err_tot))


  !     CALL MPI_ALLREDUCE(norm_inf2,norm_inf_global2,1,MPI_REAL8,MPI_MAX,COMM_CART,merror)
  !     norm_inf2 = norm_inf_global2
     !-------------------------------------------------------------------------------------------------------
     !--- write out stats

   ! if ((substep ==1) .and. (rank==0) .and. (timestep .eq. 2))  OPEN(UNIT=139,FILE="rel_err_file.txt",FORM="FORMATTED",position="append",STATUS="OLD",ACTION="READWRITE")
   ! if ((substep ==1) .and. (rank==0) .and. (timestep .eq. 2)) WRITE(UNIT=139, FMT=*) (err_tot_old/err_tot)


!    if ((substep ==1) .and. (rank==0) .and. (timestep .eq. 2))  OPEN(UNIT=139,FILE="err_file.txt",FORM="FORMATTED",position="append",STATUS="OLD",ACTION="READWRITE")
!    if ((substep ==1) .and. (rank==0) .and. (timestep .eq. 2)) WRITE(UNIT=139, FMT=*) err_tot





!        if (rank ==0 ) print*, iter, norm_inf
        !if (norm_inf .le. 2.0e-5) EXIT ITERATE 
       if (rank ==0 ) print*, iter, mode_j*8.*0.8*3.1415*3.1415/(M1*M1), err_tot, rel_err_tot
       !if (rel_err_tot .le. (1*10.0e-5)) EXIT ITERATE
       if (rel_err_tot .le. (20.*mode_j*8.*1./3.0*0.8*3.1415*3.1415*(1./(M1*M1)+1./(M2*M2)+1./(M3*M3)))) EXIT ITERATE
      !  if (iter .eq. 2) EXIT ITERATE 
       rel_old = rel
       err_tot_old=err_tot

       endif
!    if (iter == 1) then
!       istat = cudaMemcpy(rel_dev, c_loc(rel), nx * ny * nz * szreal, cudaMemcpyHostToDevice)
!       call cuda_err_check(istat)
!    endif

    ! copying ghost cells


    istat = cudaMemcpy(west_ghost_dev, c_loc(west_ghost), (N2-1) * (N3-1) * szreal, cudaMemcpyHostToDevice)
    call cuda_err_check(istat)
    istat = cudaMemcpy(east_ghost_dev, c_loc(east_ghost), (N2-1) * (N3-1) * szreal, cudaMemcpyHostToDevice)
    call cuda_err_check(istat)
    istat = cudaMemcpy(lower_ghost_dev, c_loc(lower_ghost), (N1-1) * (N3-1) * szreal, cudaMemcpyHostToDevice)
    call cuda_err_check(istat)
    istat = cudaMemcpy(upper_ghost_dev, c_loc(upper_ghost), (N1-1) * (N3-1) * szreal, cudaMemcpyHostToDevice)
    call cuda_err_check(istat)
    istat = cudaMemcpy(rear_ghost_dev, c_loc(rear_ghost), (N1-1) * (N2-1) * szreal, cudaMemcpyHostToDevice)
    call cuda_err_check(istat)
    istat = cudaMemcpy(front_ghost_dev, c_loc(front_ghost), (N1-1) * (N2-1) * szreal, cudaMemcpyHostToDevice)
    call cuda_err_check(istat)
   ! print*, 'done'


!    print*, maxval(front_ghost)

!!    istat = cudaMemcpy(comp_dev, c_loc(comp), nx * ny * nz * szreal, cudaMemcpyHostToDevice)
!!    call cuda_err_check(istat)



! Creating CUDA Streams  !TEST -- wrapper is buggy -->  Maybe c_loc(stream1) can help
!    istat = cudaStreamCreate(stream1)
!    call cuda_err_check(istat)
!    istat = cudaStreamCreate(stream2)
!    call cuda_err_check(istat)
!    istat = cudaStreamCreate(stream3)
!    call cuda_err_check(istat)




   ! Setup CUDA compute grid configuration.

#ifdef CUDAVERSION9

        istat = cudaConfigureCall(blocks, threads, int8(0), c_null_ptr)
        call cuda_err_check(istat)
    
    !!! Stream 1 : the cuda API break down
    
        ! Setup CUDA kernel arguments (individually)
        offset = 0
        istat = cudaSetupScalarArgument(c_loc(grid_id), szint, offset)
        call cuda_err_check(istat)
        offset = offset + sizeof(grid_id)
        istat = cudaSetupScalarArgument(c_loc(nx_dev), szint, offset)
        call cuda_err_check(istat)
        offset = offset + sizeof(nx_dev)
        istat = cudaSetupScalarArgument(c_loc(ny_dev), szint, offset)
        call cuda_err_check(istat)
        offset = offset + sizeof(ny_dev)
        istat = cudaSetupScalarArgument(c_loc(nz_dev), szint, offset)
        call cuda_err_check(istat)
        offset = offset + sizeof(nz_dev)
    
        ! only the second index of the cdg matrices are not fixed and have to updated on the GPU
    
        istat = cudaSetupScalarArgument(c_loc(ny_cdg1), szint, offset)
        call cuda_err_check(istat)
        offset = offset + sizeof(ny_cdg1)
        istat = cudaSetupScalarArgument(c_loc(ny_cdg2), szint, offset)
        call cuda_err_check(istat)
        offset = offset + sizeof(ny_cdg2)
        istat = cudaSetupScalarArgument(c_loc(ny_cdg3), szint, offset)
        call cuda_err_check(istat)
        offset = offset + sizeof(ny_cdg3)
        istat = cudaSetupScalarArgument(c_loc(nx_cdg1), szint, offset)
        call cuda_err_check(istat)
        offset = offset + sizeof(nx_cdg1)
   
    
    
    
       ! Previous arguments were 4 bytes in size, total offset is 12 bytes.
        ! Since the next 3 arguments are 8 bytes each and have to be aligned
        ! by 8 bytes boundary, it's necessary to insert a 4-byte spacing here
        ! (from 12 to 16).
    !    offset = offset + 4
        istat = cudaSetupScalarArgument(c_loc(omega), szreal, offset)
        call cuda_err_check(istat)
        offset = offset + sizeof(omega)
        istat = cudaSetupArrayArgument(cdg1_dev, szptr, offset)
        call cuda_err_check(istat)
        offset = offset + sizeof(cdg1_dev)
        istat = cudaSetupArrayArgument(cdg2_dev, szptr, offset)
        call cuda_err_check(istat)
        offset = offset + sizeof(cdg2_dev)
        istat = cudaSetupArrayArgument(cdg3_dev, szptr, offset)
        call cuda_err_check(istat)
        offset = offset + sizeof(cdg3_dev)
        istat = cudaSetupArrayArgument(rel_dev, szptr, offset)
        call cuda_err_check(istat)
        offset = offset + sizeof(rel_dev)
        istat = cudaSetupArrayArgument(bb_dev, szptr, offset)
        call cuda_err_check(istat)
        offset = offset + sizeof(bb_dev)
  !      istat = cudaSetupArrayArgument(comp_dev, szptr, offset)
  !      call cuda_err_check(istat)
  !      offset = offset + sizeof(comp_dev)
        istat = cudaSetupArrayArgument(rear_ghost_dev, szptr, offset)
        call cuda_err_check(istat)
        offset = offset + sizeof(rear_ghost_dev)   
        istat = cudaSetupArrayArgument(front_ghost_dev, szptr, offset)
        call cuda_err_check(istat)
        offset = offset + sizeof(front_ghost_dev)
        istat = cudaSetupArrayArgument(upper_ghost_dev, szptr, offset)
        call cuda_err_check(istat)
        offset = offset + sizeof(upper_ghost_dev)
        istat = cudaSetupArrayArgument(lower_ghost_dev, szptr, offset)
        call cuda_err_check(istat)
        offset = offset + sizeof(lower_ghost_dev)
        istat = cudaSetupArrayArgument(west_ghost_dev, szptr, offset)
        call cuda_err_check(istat)
        offset = offset + sizeof(west_ghost_dev)
        istat = cudaSetupArrayArgument(east_ghost_dev, szptr, offset)
        call cuda_err_check(istat)
        offset = offset + sizeof(east_ghost_dev)
        ! Finally, launch CUDA kernel after we configured everything.
        ! Kernel is identified by C-string with its name (must be
        ! undecorated, i.e. with extern "C")


!        call MPI_Barrier(MPI_COMM_WORLD, merror)
!        t1 =   MPI_Wtime()

        istat = cudaLaunch('jor_kernel_v2' // c_null_char)

!        call MPI_Barrier(MPI_COMM_WORLD, merror)
 !       t2 =   MPI_Wtime()
!        print*, nx, ny, nz, 'compute time ', t2-t1

        call cuda_err_check(istat)




!--------------------------------- get  chunks------------------------------
!---------------------------------------------------------------------------
        istat = cudaConfigureCall(blocks, threads, int8(0), c_null_ptr)
        call cuda_err_check(istat)
    
    !!! Stream 1 : the cuda API break down
    
        ! Setup CUDA kernel arguments (individually)
        offset = 0
        istat = cudaSetupScalarArgument(c_loc(nx_dev), szint, offset)
        call cuda_err_check(istat)
        offset = offset + sizeof(nx_dev)
        istat = cudaSetupScalarArgument(c_loc(ny_dev), szint, offset)
        call cuda_err_check(istat)
        offset = offset + sizeof(ny_dev)
        istat = cudaSetupScalarArgument(c_loc(nz_dev), szint, offset)
        call cuda_err_check(istat)
        offset = offset + sizeof(nz_dev)
   
    
       ! Previous arguments were 4 bytes in size, total offset is 12 bytes.
        ! Since the next 3 arguments are 8 bytes each and have to be aligned
        ! by 8 bytes boundary, it's necessary to insert a 4-byte spacing here
        ! (from 12 to 16).
        offset = offset + 4
        istat = cudaSetupArrayArgument(rel_dev, szptr, offset)
        call cuda_err_check(istat)
        offset = offset + sizeof(rel_dev)
        istat = cudaSetupArrayArgument(rear_ghost_dev, szptr, offset)
        call cuda_err_check(istat)
        offset = offset + sizeof(rear_ghost_dev)   
        istat = cudaSetupArrayArgument(front_ghost_dev, szptr, offset)
        call cuda_err_check(istat)
        offset = offset + sizeof(front_ghost_dev)
        istat = cudaSetupArrayArgument(upper_ghost_dev, szptr, offset)
        call cuda_err_check(istat)
        offset = offset + sizeof(upper_ghost_dev)
        istat = cudaSetupArrayArgument(lower_ghost_dev, szptr, offset)
        call cuda_err_check(istat)
        offset = offset + sizeof(lower_ghost_dev)
        istat = cudaSetupArrayArgument(west_ghost_dev, szptr, offset)
        call cuda_err_check(istat)
        offset = offset + sizeof(west_ghost_dev)
        istat = cudaSetupArrayArgument(east_ghost_dev, szptr, offset)
        call cuda_err_check(istat)
        offset = offset + sizeof(east_ghost_dev)
        ! Finally, launch CUDA kernel after we configured everything.
        ! Kernel is identified by C-string with its name (must be
        ! undecorated, i.e. with extern "C")


!        call MPI_Barrier(MPI_COMM_WORLD, merror)
!        t1 =   MPI_Wtime()

        istat = cudaLaunch('get_ghost' // c_null_char)

!        call MPI_Barrier(MPI_COMM_WORLD, merror)
 !       t2 =   MPI_Wtime()
!        print*, nx, ny, nz, 'compute time ', t2-t1

        call cuda_err_check(istat)


#endif

        call jor_launcher((grid_id), (nx_dev), (ny_dev),(nz_dev), (ny_cdg1), (ny_cdg2), ny_cdg3, nx_cdg1,omega,cdg1_dev,cdg2_dev,cdg3_dev, (rel_dev), (bb_dev), (rear_ghost_dev), front_ghost_dev, upper_ghost_dev, lower_ghost_dev, west_ghost_dev, east_ghost_dev)  
        call get_ghost_launcher((nx_dev), (ny_dev),(nz_dev), (rel_dev), (rear_ghost_dev), front_ghost_dev, upper_ghost_dev, lower_ghost_dev, west_ghost_dev, east_ghost_dev)  


        istat = cudaMemcpy(c_loc(west_chunk), west_ghost_dev, (N2-1) * (N3-1) * szreal, cudaMemcpyDeviceToHost)
        call cuda_err_check(istat)
        istat = cudaMemcpy(c_loc(east_chunk), east_ghost_dev, (N2-1) * (N3-1) * szreal,  cudaMemcpyDeviceToHost)
        call cuda_err_check(istat)
        istat = cudaMemcpy(c_loc(lower_chunk), lower_ghost_dev, (N1-1) * (N3-1) * szreal,  cudaMemcpyDeviceToHost)
        call cuda_err_check(istat)
        istat = cudaMemcpy(c_loc(upper_chunk), upper_ghost_dev, (N1-1) * (N3-1) * szreal,  cudaMemcpyDeviceToHost)
        call cuda_err_check(istat)
        istat = cudaMemcpy(c_loc(rear_chunk), rear_ghost_dev, (N1-1) * (N2-1) * szreal,  cudaMemcpyDeviceToHost)
        call cuda_err_check(istat)
        istat = cudaMemcpy(c_loc(front_chunk), front_ghost_dev, (N1-1) * (N2-1) * szreal,  cudaMemcpyDeviceToHost)
        call cuda_err_check(istat)






!=============================================================================================================
!=============================================================================================================



    ! Copy back results from GPU to CPU memory.
!----for the jacobi version
!    istat = cudaMemcpy(c_loc(rel), comp_dev, nx * ny * nz * szreal, cudaMemcpyDeviceToHost)
!    call cuda_err_check(istat)
!----for the jor version


!    istat = cudaMemcpy(c_loc(rel), rel_dev, nx * ny * nz * szreal, cudaMemcpyDeviceToHost)
!    call cuda_err_check(istat)



    ! compute the global error
 !   if ((mod(iter, 2) == 0) .or. (mod(iter, 2) == 1)) then
  !     if  (mod(iterate, 2) == 0) then
!          istat = cudaMemcpy(c_loc(rel), comp_dev, nx * ny * nz * szreal, cudaMemcpyDeviceToHost)
!          call cuda_err_check(istat)   
 !      else  
 !         istat = cudaMemcpy(c_loc(rel_new), comp_dev, nx * ny * nz * szreal, cudaMemcpyDeviceToHost)
 !         call cuda_err_check(istat) 
 !      end if

!   endif



    end do iterate

    ! Free allocated GPU memory.
!if (((timestep ==n_timesteps) .and. (substep == 3))) then

!    istat = cudaFree(rel_dev)
!    call cuda_err_check(istat)
!    istat = cudaFree(bb_dev)
!    call cuda_err_check(istat)
!    istat = cudaFree(cdg1_dev)
!    call cuda_err_check(istat)
!    istat = cudaFree(cdg2_dev)
!    call cuda_err_check(istat)
!    istat = cudaFree(cdg3_dev)
!    call cuda_err_check(istat)
!
!    istat = cudaFree(rear_ghost_dev)
!    call cuda_err_check(istat)
!    istat = cudaFree(front_ghost_dev)
!    call cuda_err_check(istat)
!    istat = cudaFree(west_ghost_dev)
!    call cuda_err_check(istat)
!    istat = cudaFree(east_ghost_dev)
!    call cuda_err_check(istat)
!    istat = cudaFree(upper_ghost_dev)
!    call cuda_err_check(istat)
!    istat = cudaFree(lower_ghost_dev)
!    call cuda_err_check(istat)

!endif

! Applying boundary conditions

end subroutine device_jacobi_solver

