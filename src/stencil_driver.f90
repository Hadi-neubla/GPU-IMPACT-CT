!**************************************************************************************************                                           
!* by Hadi Zolfaghari, ARTORG Center, University of Bern, (zolfaghari.haadi@gmail.com)            *      
!* October 2015 - March 2020                                                                      *            
!**************************************************************************************************
subroutine cuda_err_check(istat)



use iso_c_binding

! Module where we define CUDA functions to use them from Fortran
use mod_cudafor

implicit none

  integer, intent(in) :: istat

  ! Variables for error handling
  type(c_ptr) :: msg_ptr
  integer, parameter :: msg_len = 1024
  character(len=msg_len), pointer :: msg

  if (istat .ne. 0) then
      ! Get void* pointer to the error string
      msg_ptr = cudaGetErrorString(istat)
      ! Cast void* error string pointer to string
      call c_f_pointer(msg_ptr, msg)
      write(*, *) 'CUDA Error: ', &
          ! Take entire message string until first occurence of '\0' symbol
          ! -- end-of-string marker in C
          msg(1:index(msg, c_null_char))
      return
  endif

end subroutine cuda_err_check



subroutine device_stencil(dimension_id, lower_ind, upper_ind, cS1, cS2, cS3, phi, output)

USE mod_dims
USE mod_vars
USE mod_exchange
USE mod_diff
USE mod_laplace
USE mod_helmholtz
USE mod_inout




use iso_c_binding

! Module where we define CUDA functions to use them from Fortran
use mod_cudafor

implicit none


   integer, target, intent(in)  ::  dimension_id  ! the dimension_id of output --> corresponding to 'm=1,2,3' in IMPACT
   integer, target, intent(in)  ::  lower_ind, upper_ind

   REAL(8)   , target, INTENT(inout) ::  phi (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))!< ? 
   REAL(8)   , target, INTENT(  out) ::  output(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))!< stencil operator output

   REAL(8)   , target, INTENT(in)    ::  cS1(lower_ind:upper_ind,0:N1) 
   REAL(8)   , target, INTENT(in)    ::  cS2(lower_ind:upper_ind,0:N2) 
   REAL(8)   , target, INTENT(in)    ::  cS3(lower_ind:upper_ind,0:N3) 


  ! Declare the size of phi and other arguments passing to stencil operator -- ref : mod_diff

   integer   ::   phi_size(3)
   integer   ::   cS1_size(2)
   integer   ::   cS2_size(2)
   integer   ::   cS3_size(2)



  integer, target     :: nx, ny, nz ! the output and phi size components
  integer, target     :: nxs, nys !the cS1 size components
!  real(8), target, dimension(b1L:(N1+b1U), b2L:(N2+b2U), b3L:(N3+b3U)) :: output 
!  real(8), target, allocatable, dimension(:,:,:) :: phi_h, cS1_h

  type(c_ptr) :: d_x, d_y, d_xy
  type(c_ptr) :: output_dev, phi_dev, cS1_dev, cS2_dev, cS3_dev
    
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


!! Aplying the exchange !!TODO check this for every stencil operation !! On the CPU

   CALL exchange(dimension_id,0,phi) 



  szreal = sizeof(output(1, 1, 1))
  szint = sizeof(nx)
  szptr = sizeof(d_x)
 
 
 ! We repeat GPU logic twice to show that the first kernel invocation includes
  ! the time of GPU device initialization, which is done one time per
  ! application run.

!    call cpu_time(start)



    ! Get the size of arrays passing to stencil subroutine
    phi_size = shape(phi)

     

    if (dimension_id .eq. 1)  cS1_size = shape(cS1)
    if (dimension_id .eq. 2)  cS2_size = shape(cS2)
    if (dimension_id .eq. 3)  cS3_size = shape(cS3)


    nx = phi_size(1)  
    ny = phi_size(2) 
    nz = phi_size(3) 
    if      (dimension_id .eq. 1) then
        nxs = cS1_size(1)
        nys = cS1_size(2)
    else if (dimension_id .eq. 2) then
        nxs = cS2_size(1)
        nys = cS2_size(2)
    else if (dimension_id .eq. 3) then
        nxs = cS3_size(1)
        nys = cS3_size(2)
    endif


!    print*, nx, ny, nz, nxs, nys, 'nx, ny, nz, nxs, nys', 'this is in output'
    ! Allocate memory on GPU.

    if      (dimension_id .eq. 1) then
        istat = cudaMalloc(cS1_dev, nxs * nys * szreal)
        call cuda_err_check(istat)
    else if (dimension_id .eq. 2) then
        istat = cudaMalloc(cS2_dev, nxs * nys * szreal)
        call cuda_err_check(istat)
    else if (dimension_id .eq. 3) then
        istat = cudaMalloc(cS3_dev, nxs * nys * szreal)
        call cuda_err_check(istat)  
    end if


    istat = cudaMalloc(phi_dev, nx * ny * nz * szreal)
    call cuda_err_check(istat)
    ! phi and output come in size size -- ref to mod_diff
    istat = cudaMalloc(output_dev, nx * ny * nz *  szreal)
    call cuda_err_check(istat)

	! XXX Several blocks of code below are equivalent to sum_kernel<<<...>>> in C.


!    print*, nx/8, ny/8, nz/8, 'block dimmension on gpu'

    ! Create CUDA compute grid configuration.
    blocks = dim3((nx-8)/8, (ny-8)/8, (nz-8)/8)
    if (mod(nx, szblock) .ne. 0) blocks%x = blocks%x + 1
    threads = dim3((nx-8)/blocks%x, (ny-8)/blocks%y, (nz-8)/blocks%z)

!!    blocks = dim3(1, 1, 1)
!!    threads = dim3(1, 1, 1)
    ! Copy input data from CPU to GPU memory.
    if      (dimension_id .eq. 1) then
        istat = cudaMemcpy(cS1_dev, c_loc(cS1), nxs * nys * szreal, cudaMemcpyHostToDevice)
        call cuda_err_check(istat)
    else if (dimension_id .eq. 2) then
        istat = cudaMemcpy(cS2_dev, c_loc(cS2), nxs * nys * szreal, cudaMemcpyHostToDevice)
        call cuda_err_check(istat)
    else if (dimension_id .eq. 3) then
        istat = cudaMemcpy(cS3_dev, c_loc(cS3), nxs * nys * szreal, cudaMemcpyHostToDevice)
        call cuda_err_check(istat)
    end if
    



    istat = cudaMemcpy(phi_dev, c_loc(phi), nx * ny * nz * szreal, cudaMemcpyHostToDevice)
    call cuda_err_check(istat)

! Creating CUDA Streams  !TEST -- wrapper is buggy -->  Maybe c_loc(stream1) can help
!    istat = cudaStreamCreate(stream1)
!    call cuda_err_check(istat)
!    istat = cudaStreamCreate(stream2)
!    call cuda_err_check(istat)
!    istat = cudaStreamCreate(stream3)
!    call cuda_err_check(istat)




!    istat = cudaMemcpy(output_dev, c_loc(output), nx * ny * nz * szreal, cudaMemcpyHostToDevice)
!    call cuda_err_check(istat)


!    print*, threads, blocks
!    print*, nx, ny, nz
   ! Setup CUDA compute grid configuration.
    istat = cudaConfigureCall(blocks, threads, int8(0), c_null_ptr)
    call cuda_err_check(istat)

!!! Stream 1 : the cuda API break down

    ! Setup CUDA kernel arguments (individually)
    offset = 0
    istat = cudaSetupScalarArgument(c_loc(dimension_id), szint, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(dimension_id)
    istat = cudaSetupScalarArgument(c_loc(lower_ind), szint, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(lower_ind)
    istat = cudaSetupScalarArgument(c_loc(upper_ind), szint, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(upper_ind)
    istat = cudaSetupScalarArgument(c_loc(nx), szint, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(nx)
    istat = cudaSetupScalarArgument(c_loc(ny), szint, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(ny)
    istat = cudaSetupScalarArgument(c_loc(nz), szint, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(nz)
    istat = cudaSetupScalarArgument(c_loc(nxs), szint, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(nxs)
    istat = cudaSetupScalarArgument(c_loc(nys), szint, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(nys)

   ! Previous arguments were 4 bytes in size, total offset is 12 bytes.
    ! Since the next 3 arguments are 8 bytes each and have to be aligned
    ! by 8 bytes boundary, it's necessary to insert a 4-byte spacing here
    ! (from 12 to 16).
   !  offset = offset + 4
    istat = cudaSetupArrayArgument(phi_dev, szptr, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(phi_dev)
    if      (dimension_id .eq. 1) then
        istat = cudaSetupArrayArgument(cS1_dev, szptr, offset)
        call cuda_err_check(istat)
        offset = offset + sizeof(cS1_dev)
    else if (dimension_id .eq. 2) then
        istat = cudaSetupArrayArgument(cS2_dev, szptr, offset)
        call cuda_err_check(istat)
        offset = offset + sizeof(cS2_dev)
    else if (dimension_id .eq. 3) then
        istat = cudaSetupArrayArgument(cS3_dev, szptr, offset)
        call cuda_err_check(istat)
        offset = offset + sizeof(cS3_dev)
    end if




    istat = cudaSetupArrayArgument(output_dev, szptr, offset)
    call cuda_err_check(istat)

    ! Finally, launch CUDA kernel after we configured everything.
    ! Kernel is identified by C-string with its name (must be
    ! undecorated, i.e. with extern "C")
    istat = cudaLaunch('stencil_kernel' // c_null_char)
    call cuda_err_check(istat)



    ! Copy back results from GPU to CPU memory.
    istat = cudaMemcpy(c_loc(output), output_dev, nx * ny * nz * szreal, cudaMemcpyDeviceToHost)
!    istat = cudaMemcpy(c_loc(cS1), cS1_dev, nxs* nys* szreal ,cudaMemcpyDeviceToHost)
    call cuda_err_check(istat)



    ! Free allocated GPU memory.
    istat = cudaFree(output_dev)
    call cuda_err_check(istat)
    istat = cudaFree(phi_dev)
    call cuda_err_check(istat)
    if      (dimension_id .eq. 1) then
        istat = cudaFree(cS1_dev)
        call cuda_err_check(istat)
    else if (dimension_id .eq. 2) then
        istat = cudaFree(cS2_dev)
        call cuda_err_check(istat)
    else if (dimension_id .eq. 3) then
        istat = cudaFree(cS3_dev)
        call cuda_err_check(istat)
    end if


! Applying boundary conditions
! This part is on the CPU and !!TODO: has to be shared-memory optimized


end subroutine device_stencil

