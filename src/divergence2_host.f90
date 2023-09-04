!**************************************************************************************************                                    // !* by Hadi Zolfaghari, ARTORG Center, University of Bern, (zolfaghari.haadi@gmail.com)            *      
!* October 2015 - March 2020                                                                      *            
!* Modified by Hadi Zolfaghari, University of Cambridge (hz382@damtp.cam.ac.uk)                   *      
!* April 2020 -                                                                                   *            
!**************************************************************************************************
subroutine device_divergence2( phi, div)

USE mod_dims
USE mod_vars
USE mod_vars_GPU
USE mod_exchange
USE mod_diff
USE mod_laplace
USE mod_inout

use iso_c_binding

! Module where we define CUDA functions to use them from Fortran
use mod_cudafor

implicit none





   real(8),    target, intent(inout)       ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),3)         
   real(8),    target, intent(out)    ::  div(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))         

  ! Declare the size of phi and other arguments passing to stencil operator -- ref : mod_diff

 !  integer   ::   phi_size(4)

!   integer   ::   cS1_size(2)
!   integer   ::   cS2_size(2)
!   integer   ::   cS3_size(2)


!  integer, target     :: nx, ny, nz ! the output and phi size components
!  integer, target     :: nx_div_coef1, ny_div_coef1, nx_div_coef2, ny_div_coef2, nx_div_coef3, ny_div_coef3 !the cS1 size components

!  type(c_ptr) :: d_x, d_y, d_xy
!  type(c_ptr) :: div2_dev, phi_div2_dev, cDu1_dev, cDv2_dev, cDw3_dev
    
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


  real(8)            ::  t1, t2
  include 'mpif.h'

!! Aplying the exchange !!TODO check this for every stencil operation !! On the CPU

!   CALL exchange(dimension_id,0,phi) 

!  szreal = sizeof(nl(1, 1, 1,1))
!  szint = sizeof(nx)
!  szptr = sizeof(phi_div2_dev)
 
 
 ! We repeat GPU logic twice to show that the first kernel invocation includes
  ! the time of GPU device initialization, which is done one time per
  ! application run.

!    call cpu_time(start)

    ! Get the size of arrays passing to stencil subroutine
!    phi_size = shape(phi)
!    nx = phi_size(1)  
!    ny = phi_size(2) 
!    nz = phi_size(3) 



    if (GPU_verbose == 1) print*, 'On GPU divergence2 ...'
    

!    cS1_size = shape(cDu1)
!    nx_div_coef1 = cS1_size(1)
!    ny_div_coef1 = cS1_size(2)
!    cS2_size = shape(cDv2)
!    nx_div_coef2 = cS2_size(1)
!    ny_div_coef2 = cS2_size(2)
!    cS3_size = shape(cDw3)
!    nx_div_coef3 = cS3_size(1)
!    ny_div_coef3 = cS3_size(2)
!    istat = cudaMalloc(cDu1_dev, nx_div_coef1 * ny_div_coef1 * szreal)
!    call cuda_err_check(istat)
!    istat = cudaMemcpy(cDu1_dev, c_loc(cDu1), nx_div_coef1 * ny_div_coef1 * szreal, cudaMemcpyHostToDevice)
!    call cuda_err_check(istat)
!    istat = cudaMalloc(cDv2_dev, nx_div_coef2 * ny_div_coef2 * szreal)
!    call cuda_err_check(istat)
!    istat = cudaMemcpy(cDv2_dev, c_loc(cDv2), nx_div_coef2 * ny_div_coef2 * szreal, cudaMemcpyHostToDevice)
!    call cuda_err_check(istat)
!    istat = cudaMalloc(cDw3_dev, nx_div_coef3 * ny_div_coef3 * szreal)
!    call cuda_err_check(istat)
!    istat = cudaMemcpy(cDw3_dev, c_loc(cDw3), nx_div_coef3 * ny_div_coef3 * szreal, cudaMemcpyHostToDevice)
!    call cuda_err_check(istat)



!    print*, nx, ny, nz, nxs, nys, 'nx, ny, nz, nxs, nys', 'this is in output'
    ! Allocate memory on GPU.

   ! phi and output come in size size -- ref to mod_diff
!    istat = cudaMalloc(div2_dev, nx_dev * ny_dev * nz_dev *  szreal)
!    call cuda_err_check(istat)
!
!    istat = cudaMalloc(phi_div2_dev, nx_dev * ny_dev * nz_dev * 3 * szreal)
!    call cuda_err_check(istat)

! XXX Several blocks of code below are equivalent to sum_kernel<<<...>>> in C.


#ifdef CUDAVERSION9
    ! Create CUDA compute grid configuration.
    blocks = dim3((nx_dev-8)/8, (ny_dev-8)/8, (nz_dev-8)/8)
    if (mod(nx_dev, szblock) .ne. 0) blocks%x = blocks%x + 1
    threads = dim3((nx_dev-8)/blocks%x, (ny_dev-8)/blocks%y, (nz_dev-8)/blocks%z)
#endif 

!!    blocks = dim3(1, 1, 1)
!!    threads = dim3(1, 1, 1)
    ! Copy input data from CPU to GPU memory.
   

   !!  TODO velocity should be transferred component-wise
    istat = cudaMemcpy(phi_div2_dev, c_loc(phi), nx_dev * ny_dev * nz_dev * 3 * szreal, cudaMemcpyHostToDevice)
    call cuda_err_check(istat)
!    istat = cudaMemcpy(div2_dev, c_loc(div), nx * ny * nz * szreal, cudaMemcpyHostToDevice)
!    call cuda_err_check(istat)


! Creating CUDA Streams  !TEST -- wrapper is buggy -->  Maybe c_loc(stream1) can help
!    istat = cudaStreamCreate(stream1)
!    call cuda_err_check(istat)
!    istat = cudaStreamCreate(stream2)
!    call cuda_err_check(istat)
!    istat = cudaStreamCreate(stream3)
!    call cuda_err_check(istat)


!    istat = cudaMemcpy(output_dev, c_loc(output), nx * ny * nz * szreal, cudaMemcpyHostToDevice)
!    call cuda_err_check(istat)

#ifdef CUDAVERSION9
   ! Setup CUDA compute grid configuration.
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
    istat = cudaSetupScalarArgument(c_loc(nx_div_coef1), szint, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(nx_div_coef1)
     istat = cudaSetupScalarArgument(c_loc(nx_div_coef2), szint, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(nx_div_coef2)
    istat = cudaSetupScalarArgument(c_loc(nx_div_coef3), szint, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(nx_div_coef3)



   ! Previous arguments were 4 bytes in size, total offset is 12 bytes.
    ! Since the next 3 arguments are 8 bytes each and have to be aligned
    ! by 8 bytes boundary, it's necessary to insert a 4-byte spacing here
    ! (from 12 to 16).
!    offset = offset + 4

    istat = cudaSetupArrayArgument(phi_div2_dev, szptr, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(phi_div2_dev)

    istat = cudaSetupArrayArgument(div2_dev, szptr, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(div2_dev)

    istat = cudaSetupArrayArgument(cDu1_dev, szptr, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(cDu1_dev)

    istat = cudaSetupArrayArgument(cDv2_dev, szptr, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(cDv2_dev)

    istat = cudaSetupArrayArgument(cDw3_dev, szptr, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(cDw3_dev)

    ! Finally, launch CUDA kernel after we configured everything.
    ! Kernel is identified by C-string with its name (must be
    ! undecorated, i.e. with extern "C")
   !  call MPI_Barrier(MPI_COMM_WORLD, merror)
   !  t1 =   MPI_Wtime()


    istat = cudaLaunch('divergence2_kernel' // c_null_char)
!    call cuda_err_check(istat)

#endif

    
  call divergence2_launcher( (nx_dev), (ny_dev),(nz_dev), (nx_div_coef1), (nx_div_coef2), nx_div_coef3,  (phi_div2_dev), (div2_dev), cDu1_dev, cDv2_dev, cDw3_dev)

 !    call MPI_Barrier(MPI_COMM_WORLD, merror)
  !   t2 =   MPI_Wtime()


!     print*,  'compute time helmholtz' , t2-t1

    ! Copy back results from GPU to CPU memory.
    istat = cudaMemcpy(c_loc(div), div2_dev, nx_dev * ny_dev * nz_dev * szreal, cudaMemcpyDeviceToHost)
    call cuda_err_check(istat)

!    istat = cudaMemcpy(c_loc(phi), phi_div2_dev, nx * ny * nz * szreal, cudaMemcpyDeviceToHost)
!    call cuda_err_check(istat)



    ! Free allocated GPU memory.
!    istat = cudaFree(div2_dev)
!    call cuda_err_check(istat)
!    istat = cudaFree(phi_div2_dev)
!    call cuda_err_check(istat)
!    istat = cudaFree(cDu1_dev)
!    call cuda_err_check(istat)
!    istat = cudaFree(cDv2_dev)
!    call cuda_err_check(istat)
!    istat = cudaFree(cDw3_dev)
!    call cuda_err_check(istat)


! Applying boundary conditions


end subroutine device_divergence2

