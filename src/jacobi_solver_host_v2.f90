!**************************************************************************************************                                           
!* by Hadi Zolfaghari, ARTORG Center, University of Bern, (zolfaghari.haadi@gmail.com)            *      
!* October 2015 - March 2020                                                                      *            
!**************************************************************************************************
subroutine device_jacobi_solver(grid_id,omega,SS1,SS2,SS3,NN1,NN2,NN3,rel,bb, comp)

USE mod_dims
USE mod_vars
USE mod_exchange
USE mod_diff
USE mod_laplace
USE mod_inout




use iso_c_binding

! Module where we define CUDA functions to use them from Fortran
use mod_cudafor

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
   REAL(8)   , target, INTENT(inout) ::  bb (b1L:(NN(1,grid_id)+b1U),b2L:(NN(2,grid_id)+b2U),b3L:(NN(3,grid_id)+b3U))!< ? 
   REAL(8)   , target, INTENT(  out) ::  comp(b1L:(NN(1,grid_id)+b1U),b2L:(NN(2,grid_id)+b2U),b3L:(NN(3,grid_id)+b3U))!< gradient operator


   REAL(8)   , target ::  rel_old(b1L:(NN(1,grid_id)+b1U),b2L:(NN(2,grid_id)+b2U),b3L:(NN(3,grid_id)+b3U))!< ? 
   REAL(8)   , target ::  GPUres(1:N1p,1:N2p,1:N3p)!< ? 

   REAL(8)   , target ::  west_ghost(1:N2-1,1:N3-1)
   REAL(8)   , target ::  east_ghost(1:N2-1,1:N3-1)
   REAL(8)   , target ::  front_ghost(1:N1-1,1:N2-1)
   REAL(8)   , target ::  rear_ghost(1:N1-1,1:N2-1)
   REAL(8)   , target ::  lower_ghost(1:N1-1,1:N3-1)
   REAL(8)   , target ::  upper_ghost(1:N1-1,1:N3-1)






!   REAL(8)   , target ::  slide_south(b1L:(NN(1,grid_id)+b1U),b2L:(NN(2,grid_id)+b2U),b3L:(NN(3,grid_id)+b3U))!< ? 
!   REAL(8)   , target ::  slide_north(b1L:(NN(1,grid_id)+b1U),b2L:(NN(2,grid_id)+b2U),b3L:(NN(3,grid_id)+b3U))!< ? 



  ! Declare the size of rel and other arguments passing to jacobi relaxation -- ref : mod_laplace

   integer   ::   rel_size(3)



   integer   ::   cdg1_size(3)   
   integer   ::   cdg2_size(3)
   integer   ::   cdg3_size(3)

   integer   ::   N1R, N2R, N3R

  ! Declare sincos kernel arguments with "target" -- for c_loc
  integer, target     :: nx, ny, nz ! the rel, bb and comp size components

  integer, target     :: nx_cdg1, ny_cdg1, nz_cdg1 ! the grad and phi size components
  integer, target     :: nx_cdg2, ny_cdg2, nz_cdg2 ! the grad and phi size components
  integer, target     :: nx_cdg3, ny_cdg3, nz_cdg3 ! the grad and phi size components

!  real(8), target, allocatable, dimension(:,:,:) :: h_x, h_y, h_xy1, h_xy2
!  real(8), target, dimension(b1L:(N1+b1U), b2L:(N2+b2U), b3L:(N3+b3U)) :: grad 
!  real(8), target, allocatable, dimension(:,:,:) :: phi_h, cGp1_h

!  type(c_ptr) :: d_x, d_y, d_xy

  !type(c_ptr) :: rel_dev, bb_dev, comp_dev, cdg1_dev,  cdg2_dev, cdg3_dev
  integer :: sum_host
  integer(c_int) :: istat

  integer, parameter :: szblock = 1

  ! CUDA constants
  integer(c_int), parameter :: cudaMemcpyHostToDevice = 1
  integer(c_int), parameter :: cudaMemcpyDeviceToHost = 2

  integer(c_size_t) :: szreal, szint, szptr, offset = 0

  type(dim3) :: blocks, threads


!  type(cudaStream_t) :: stream1, stream2, stream3

  integer :: i, j, k, step
  real(8), volatile :: start, finish

  integer  ::   iter

  real(8)     ::  norm_inf, norm_inf2, norm_inf_global, norm_inf_global2, norm_zero
  real(8)     ::  t1, t2, t3
  include 'mpif.h'

!! Aplying the exchange !!TODO this part is CPU only and a subject for shared memory optimization

!   CALL exchange(dimension_id,0,phi) 



  szreal = sizeof(rel(1, 1, 1))
  szint = sizeof(nx)
  szptr = sizeof(rel_dev)
 
 
!    b1L, b1U, NN(1,1) they are -3 3 65 for 64^3 problem
 ! We repeat GPU logic twice to show that the first kernel invocation includes
  ! the time of GPU device initialization, which is done one time per
  ! application run.

!    call cpu_time(start)



    ! Get the size of arrays passing to jacobi relaxation -- > ref : mod_laplace
    rel_size = shape(rel)

     


    cdg1_size = shape(cdg1)
    cdg2_size = shape(cdg2)
    cdg3_size = shape(cdg3)


    nx = rel_size(1)  
    ny = rel_size(2) 
    nz = rel_size(3) 
!    print*, nx, ny,nz

    nx_cdg1 = cdg1_size(1)
    ny_cdg1 = cdg1_size(2)
    nz_cdg1 = cdg1_size(3)

    nx_cdg2 = cdg2_size(1)
    ny_cdg2 = cdg2_size(2)
    nz_cdg2 = cdg2_size(3)

    nx_cdg3 = cdg3_size(1)
    ny_cdg3 = cdg3_size(2)
    nz_cdg3 = cdg3_size(3)



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

if (timestep == 0) then

    istat = cudaMalloc(cdg1_dev, nx_cdg1 * ny_cdg1 * nz_cdg1 * szreal)
    call cuda_err_check(istat)

    istat = cudaMalloc(cdg2_dev, nx_cdg2 * ny_cdg2 * nz_cdg2 * szreal)
    call cuda_err_check(istat)

    istat = cudaMalloc(cdg3_dev, nx_cdg3 * ny_cdg3 * nz_cdg3 * szreal)
    call cuda_err_check(istat)  
 


    istat = cudaMalloc(rel_dev, nx * ny * nz * szreal)
    call cuda_err_check(istat)
    ! rel, bb and comp come in same size -- ref to mod_laplace
    istat = cudaMalloc(bb_dev, nx * ny * nz *  szreal)
    call cuda_err_check(istat)
    istat = cudaMalloc(comp_dev, nx * ny * nz *  szreal)
    call cuda_err_check(istat)


    ! defining ghost_cells to be copied to GPU
    istat = cudaMalloc(west_ghost_dev, (N2-1) * (N3-1) * szreal)
    call cuda_err_check(istat)
    istat = cudaMalloc(east_ghost_dev, (N2-1) * (N3-1) * szreal)
    call cuda_err_check(istat)
    istat = cudaMalloc(front_ghost_dev, (N1-1) * (N2-1) * szreal)
    call cuda_err_check(istat)
    istat = cudaMalloc(rear_ghost_dev, (N1-1) * (N2-1) * szreal)
    call cuda_err_check(istat)
    istat = cudaMalloc(upper_ghost_dev, (N1-1) * (N3-1) * szreal)
    call cuda_err_check(istat)
    istat = cudaMalloc(lower_ghost_dev, (N1-1) * (N3-1) * szreal)
    call cuda_err_check(istat)

	! XXX Several blocks of code below are equivalent to sum_kernel<<<...>>> in C.

    ! Create CUDA compute grid configuration.
    
     
!     print*, 'grid size',  (nx - 8)/8, (ny - 8)/8, (nz - 8)/8
!    blocks = dim3(1, 1, 1)
!    threads = dim3(1, 1, 1)
    ! Copy input data from CPU to GPU memory.
    istat = cudaMemcpy(cdg1_dev, c_loc(cdg1), nx_cdg1 * ny_cdg1 * nz_cdg1 * szreal, cudaMemcpyHostToDevice)
    call cuda_err_check(istat)
    istat = cudaMemcpy(cdg2_dev, c_loc(cdg2), nx_cdg2 * ny_cdg2 * nz_cdg2 * szreal, cudaMemcpyHostToDevice)
    call cuda_err_check(istat)
    istat = cudaMemcpy(cdg3_dev, c_loc(cdg3), nx_cdg3 * ny_cdg3 * nz_cdg3 * szreal, cudaMemcpyHostToDevice)
    call cuda_err_check(istat)
    
    endif


    istat = cudaMemcpy(bb_dev, c_loc(bb), nx * ny * nz * szreal, cudaMemcpyHostToDevice)
    call cuda_err_check(istat)

    istat = cudaMemcpy(rel_dev, c_loc(rel), nx * ny * nz * szreal, cudaMemcpyHostToDevice)
    call cuda_err_check(istat)
 

    blocks = dim3((nx - 8)/8, (ny - 8)/8, (nz - 8)/8)
    if (mod(nx, szblock) .ne. 0) blocks%x = blocks%x + 1
    threads = dim3((nx - 8)/blocks%x, (ny - 8)/blocks%y, (nz - 8)/blocks%z)



    iter=0
    norm_zero =0.
    iterate: do

    ! do the MPI exchange before copying to GPU

    iter  = iter+1


!    CALL exchange(1,1,rel(:,:,:))
!    CALL exchange(2,2,rel(:,:,:))
!    CALL exchange(3,3,rel(:,:,:))
!     print*, iter

     CALL exchange_relax_GPU(1,0,0,0,0,.FALSE.,rel)
!     print*, N1, N2, N3

       norm_inf = 0.
       norm_inf2= 0.
       if ( (mod(iter,20)==1))  then  !--------------------------------------------------------------------------------------------------------
       !print*, N1p
       DO k = SS3, NN3
          DO j = SS2, NN2
             DO i = SS1, NN1
                GPUres(i,j,k) = bb(i,j,k)                                     &
                            &      - cdg1(-1,i,1)*rel(i-1,j,k) - cdg1(1,i,1)*rel(i+1,j,k)     &
                            &      - cdg2(-1,j,1)*rel(i,j-1,k) - cdg2(1,j,1)*rel(i,j+1,k)     &
                            &      - cdg3(-1,k,1)*rel(i,j,k-1) - cdg3(1,k,1)*rel(i,j,k+1)    &
                            &      - (cdg1(0,i,1) + cdg2(0,j,1) + cdg3(0,k,1))*rel(i,j,k)
                norm_inf = MAX(ABS(GPUres(i,j,k)),norm_inf)
              !  norm_inf2 = MAX(ABS(rel(i,j,k)-rel_old(i,j,k)),norm_inf2)
              !  rel_old(i,j,k) = rel(i,j,k)
             END DO
          END DO
       END DO

      ! print*, 'rank ', rank, ' has the maximum value of ', maxval(abs(GPUres)), ' located in ', maxloc(abs(GPUres))
 
!    end if

    !--------------------------------------------------------------------------------------------------------
       CALL MPI_ALLREDUCE(norm_inf,norm_inf_global,1,MPI_REAL8,MPI_MAX,COMM_CART,merror) ! MPI_REDUCE bringt nichts, weil exit_yes dann mit MPI_BCAST verteilt werden müsste ...
       norm_inf = norm_inf_global
    !   if (iter == 1) norm_zero = norm_inf 


  !     CALL MPI_ALLREDUCE(norm_inf2,norm_inf_global2,1,MPI_REAL8,MPI_MAX,COMM_CART,merror) ! MPI_REDUCE bringt nichts, weil exit_yes dann mit MPI_BCAST verteilt werden müsste ...
  !     norm_inf2 = norm_inf_global2
     !========================================================================================================
        if (norm_inf .le. 1.e-2) EXIT ITERATE 
      !  if (iter .eq. 2) EXIT ITERATE 
      !  if (rank ==0 ) print*, iter, norm_inf, norm_inf2
      ! if (rank ==0 ) print*, iter, norm_inf, norm_inf2
       endif


!    if (iter == 1) then
!       istat = cudaMemcpy(rel_dev, c_loc(rel), nx * ny * nz * szreal, cudaMemcpyHostToDevice)
!       call cuda_err_check(istat)
!    endif

    ! copying ghost cells

   ! print*, S1R, N1R

    west_ghost  = rel(S1R-1,1:N2-1,1:N3-1)
    east_ghost  = rel(N1R+1,1:N2-1,1:N3-1)
    lower_ghost = rel(1:N1-1,S2R-1,1:N3-1)
    upper_ghost = rel(1:N1-1,N2R+1,1:N3-1) 
    rear_ghost  = rel(1:N1-1,1:N2-1,S3R-1)
    front_ghost = rel(1:N1-1,1:N2-1,N3R+1)

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




!!    istat = cudaMemcpy(comp_dev, c_loc(comp), nx * ny * nz * szreal, cudaMemcpyHostToDevice)
!!    call cuda_err_check(istat)



! Creating CUDA Streams  !TEST -- wrapper is buggy -->  Maybe c_loc(stream1) can help
!    istat = cudaStreamCreate(stream1)
!    call cuda_err_check(istat)
!    istat = cudaStreamCreate(stream2)
!    call cuda_err_check(istat)
!    istat = cudaStreamCreate(stream3)
!    call cuda_err_check(istat)




!    print*, threads, blocks
!    print*, nx, ny, nz
   ! Setup CUDA compute grid configuration.

        istat = cudaConfigureCall(blocks, threads, int8(0), c_null_ptr)
        call cuda_err_check(istat)
    
    !!! Stream 1 : the cuda API break down
    
        ! Setup CUDA kernel arguments (individually)
        offset = 0
        istat = cudaSetupScalarArgument(c_loc(grid_id), szint, offset)
        call cuda_err_check(istat)
        offset = offset + sizeof(grid_id)
        istat = cudaSetupScalarArgument(c_loc(nx), szint, offset)
        call cuda_err_check(istat)
        offset = offset + sizeof(nx)
        istat = cudaSetupScalarArgument(c_loc(ny), szint, offset)
        call cuda_err_check(istat)
        offset = offset + sizeof(ny)
        istat = cudaSetupScalarArgument(c_loc(nz), szint, offset)
        call cuda_err_check(istat)
        offset = offset + sizeof(nz)
    
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
        istat = cudaSetupArrayArgument(comp_dev, szptr, offset)
        call cuda_err_check(istat)
        offset = offset + sizeof(comp_dev)
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


    ! Copy back results from GPU to CPU memory.
!----for the jacobi version
!    istat = cudaMemcpy(c_loc(rel), comp_dev, nx * ny * nz * szreal, cudaMemcpyDeviceToHost)
!    call cuda_err_check(istat)
!----for the jor version


    istat = cudaMemcpy(c_loc(rel), rel_dev, nx * ny * nz * szreal, cudaMemcpyDeviceToHost)
    call cuda_err_check(istat)



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
if (timestep ==n_timesteps) then

    istat = cudaFree(rel_dev)
    call cuda_err_check(istat)
    istat = cudaFree(bb_dev)
    call cuda_err_check(istat)
    istat = cudaFree(comp_dev)
    call cuda_err_check(istat)
    istat = cudaFree(cdg1_dev)
    call cuda_err_check(istat)
    istat = cudaFree(cdg2_dev)
    call cuda_err_check(istat)
    istat = cudaFree(cdg3_dev)
    call cuda_err_check(istat)

    istat = cudaFree(rear_ghost_dev)
    call cuda_err_check(istat)
    istat = cudaFree(front_ghost_dev)
    call cuda_err_check(istat)
    istat = cudaFree(west_ghost_dev)
    call cuda_err_check(istat)
    istat = cudaFree(east_ghost_dev)
    call cuda_err_check(istat)
    istat = cudaFree(upper_ghost_dev)
    call cuda_err_check(istat)
    istat = cudaFree(lower_ghost_dev)
    call cuda_err_check(istat)

endif

! Applying boundary conditions

end subroutine device_jacobi_solver

