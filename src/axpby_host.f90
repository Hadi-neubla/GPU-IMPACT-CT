!// !**************************************************************************************************                             
!// !* by Hadi Zolfaghari, ARTORG Center, University of Bern, (zolfaghari.haadi@gmail.com)            *      
!// !* October 2015 - March 2020                                                                      *            
!// !* Modified by Hadi Zolfaghari, University of Cambridge (hz382@damtp.cam.ac.uk)                   *      
!// !* April 2020 -                                                                                   *            
!// !**************************************************************************************************
subroutine device_axpby(a, b, x, y, z)

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


   REAL(8),    target, intent(in)    ::  a,b  !z = ax + by coefficients
   

   REAL(8)   , target, INTENT(in)      ::  x (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)) 
   REAL(8)   , target, INTENT(in)      ::  y (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
   REAL(8)   , target, INTENT(inout)   ::  z (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))

  ! Declare the size of x,y

!  integer   ::   phi_size(3)

!  integer, target     :: nx, ny, nz ! the grad and phi size components
!  type(c_ptr) :: d_x, d_y, d_z
    
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


!! Aplying the exchange !!TODO this part is CPU only and a subject for shared memory optimization

!  szreal = sizeof(x(1, 1, 1))
!  szint  = sizeof(nx)
!  szptr  = sizeof(d_x)
 
 
  ! the time of GPU device initialization, which is done one time per
  ! application run.

!    call cpu_time(start)



    ! Get the size of arrays passing to the kernel
!    phi_size = shape(x)

!    nx = phi_size(1)  
!    ny = phi_size(2) 
!    nz = phi_size(3) 
!    print*, nx, ny, nz, nxg, nyg, 'nx, ny, nz, nxg, nyg', 'this is in gradient'
    ! Allocate memory on GPU.

!    istat = cudaMalloc(d_x, nx_dev * ny_dev * nz_dev * szreal)
!    call cuda_err_check(istat)
!    istat = cudaMalloc(d_y, nx_dev * ny_dev * nz_dev * szreal)
!    call cuda_err_check(istat)
!    istat = cudaMalloc(d_z, nx_dev * ny_dev * nz_dev * szreal)
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

    istat = cudaMemcpy(d_x, c_loc(x), nx_dev * ny_dev * nz_dev * szreal, cudaMemcpyHostToDevice)
    call cuda_err_check(istat)
    istat = cudaMemcpy(d_y, c_loc(y), nx_dev * ny_dev * nz_dev * szreal, cudaMemcpyHostToDevice)
    call cuda_err_check(istat)
    istat = cudaMemcpy(d_z, c_loc(z), nx_dev * ny_dev * nz_dev * szreal, cudaMemcpyHostToDevice)
    call cuda_err_check(istat)



! Creating CUDA Streams  !TEST -- wrapper is buggy -->  Maybe c_loc(stream1) can help
!    istat = cudaStreamCreate(stream1)
!    call cuda_err_check(istat)
!    istat = cudaStreamCreate(stream2)
!    call cuda_err_check(istat)
!    istat = cudaStreamCreate(stream3)
!    call cuda_err_check(istat)




!    istat = cudaMemcpy(grad_dev, c_loc(grad), nx * ny * nz * szreal, cudaMemcpyHostToDevice)
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
   ! Previous arguments were 4 bytes in size, total offset is 12 bytes.
    ! Since the next 3 arguments are 8 bytes each and have to be aligned
    ! by 8 bytes boundary, it's necessary to insert a 4-byte spacing here
    ! (from 12 to 16).
    offset = offset + 4
    istat = cudaSetupScalarArgument(c_loc(a), szreal, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(a)
    istat = cudaSetupScalarArgument(c_loc(b), szreal, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(b)
    istat = cudaSetupArrayArgument(d_x, szptr, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(d_x)
    istat = cudaSetupArrayArgument(d_y, szptr, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(d_y)
    istat = cudaSetupArrayArgument(d_z, szptr, offset)
    call cuda_err_check(istat)

    ! Finally, launch CUDA kernel after we configured everything.
    ! Kernel is identified by C-string with its name (must be
    ! undecorated, i.e. with extern "C")
    istat = cudaLaunch('axpby_kernel' // c_null_char)
    call cuda_err_check(istat)

#endif

    call axpby_launcher((nx_dev), (ny_dev), (nz_dev), a, b, (d_x), (d_y), (d_z))   

    ! Copy back results from GPU to CPU memory.
    istat = cudaMemcpy(c_loc(z), d_z, nx_dev * ny_dev * nz_dev * szreal, cudaMemcpyDeviceToHost)
!    istat = cudaMemcpy(c_loc(cGp1), cGp1_dev, nxg* nyg* szreal ,cudaMemcpyDeviceToHost)
    call cuda_err_check(istat)



    ! Free allocated GPU memory.
!    istat = cudaFree(d_x)
!    call cuda_err_check(istat)
!    istat = cudaFree(d_y)
!    call cuda_err_check(istat)
!    istat = cudaFree(d_z)
!    call cuda_err_check(istat)
!
end subroutine device_axpby

