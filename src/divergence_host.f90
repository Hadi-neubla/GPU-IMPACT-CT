!**************************************************************************************************                                    // !* by Hadi Zolfaghari, ARTORG Center, University of Bern, (zolfaghari.haadi@gmail.com)            *      
!* October 2015 - March 2020                                                                      *            
!* Modified by Hadi Zolfaghari, University of Cambridge (hz382@damtp.cam.ac.uk)                   *      
!* April 2020 -                                                                                   *            
!**************************************************************************************************

subroutine device_divergence(dimension_id,phi,div)

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


   integer, target, intent(in)      ::  dimension_id  ! the dimension_id of gradient --> corresponding to 'm=1,2,3' in IMPACT

   REAL(8)   , target, INTENT(inout)   ::  phi (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))!< ?    
   REAL(8)   , target, INTENT( inout)  ::  div(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))!< divergence operator




! for coefficient copy test!
!   real(8), target                     ::  cDtest(d1L:d1U, 0:N1)

  ! Declare the size of phi and other arguments passing to gradient -- ref : mod_diff

   integer   ::   phi_size(3)
   integer, target     :: nx, ny, nz ! the grad and phi size components


!   integer   ::   cDu1_size(2)
!   integer   ::   cDv2_size(2)
!   integer   ::   cDw3_size(2)


  ! Declare divergence kernel arguments with "target" -- for c_loc
!  integer, target     :: nx_div_coef1, ny_div_coef1 !the cDu1 size components
!  integer, target     :: nx_div_coef2, ny_div_coef2 !the cDu1 size components
!  integer, target     :: nx_div_coef3, ny_div_coef3 !the cDu1 size components


!  type(c_ptr) :: div_dev, phi_div_dev, cDu1_dev, cDv2_dev, cDw3_dev, div_dev_copy ! div_dev_copy is to avoid overloading divergence
    
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

  integer :: i, j, k, step, nmatch
  real(8), volatile :: start, finish

  real(8)           :: t1, t2, t3

  include 'mpif.h'
!! Applying the exchange !!TODO this part is CPU only and a subject for shared memory optimization

!   CALL exchange(dimension_id,0,phi) 
  CALL exchange(dimension_id, dimension_id, phi)
 

!  szreal = sizeof(div(1, 1, 1))
!  szint = sizeof(nx)
!  szptr = sizeof(d_x)
 
 


!    call cpu_time(start)



    ! Get the size of arrays passing to gradient subroutine
!    phi_size = shape(phi)
!    nx = phi_size(1)  
!    ny = phi_size(2) 
!    nz = phi_size(3) 
 
     

!    if (dimension_id .eq. 1)  cDu1_size = shape(cDu1)
!    if (dimension_id .eq. 2)  cDv2_size = shape(cDv2)
!    if (dimension_id .eq. 3)  cDw3_size = shape(cDw3)
!
!
!   if      (dimension_id .eq. 1) then
!        nx_div_coef1 = cDu1_size(1)
!        ny_div_coef1 = cDu1_size(2)
!    else if (dimension_id .eq. 2) then
!        nx_div_coef2 = cDv2_size(1)
!        ny_div_coef2 = cDv2_size(2)
!
!    else if (dimension_id .eq. 3) then
!        nx_div_coef3 = cDw3_size(1)
!        ny_div_coef3 = cDw3_size(2)
!    endif


    ! Allocate memory on GPU.

!    if      (dimension_id .eq. 1) then
!        istat = cudaMalloc(cDu1_dev, nx_div_coef1 * ny_div_coef1 * szreal)
!        call cuda_err_check(istat)
!    else if (dimension_id .eq. 2) then
!        istat = cudaMalloc(cDv2_dev, nx_div_coef2 * ny_div_coef2 * szreal)
!        call cuda_err_check(istat)
!    else if (dimension_id .eq. 3) then
!        istat = cudaMalloc(cDw3_dev, nx_div_coef3 * ny_div_coef3 * szreal)
!        call cuda_err_check(istat)  
!    end if


!    istat = cudaMalloc(phi_div_dev, nx_dev * ny_dev * nz_dev * szreal)
!    call cuda_err_check(istat)
!    istat = cudaMalloc(div_dev, nx_dev * ny_dev * nz_dev *  szreal)
!    call cuda_err_check(istat)

!    istat = cudaMalloc(div_dev_copy, nx * ny * nz *  szreal)
!    call cuda_err_check(istat)



	! XXX Several blocks of code below are equivalent to sum_kernel<<<...>>> in C.


!    print*, nx/8, ny/8, nz/8, 'block dimmension on gpu'

    ! Create CUDA compute grid configuration.
    blocks = dim3((nx_dev-8)/8, (ny_dev-8)/8, (nz_dev-8)/8)
    if (mod(nx_dev, szblock) .ne. 0) blocks%x = blocks%x + 1
    threads = dim3((nx_dev-8)/blocks%x, (ny_dev-8)/blocks%y, (nz_dev-8)/blocks%z)

!    blocks = dim3(1, 1, 1)
!    threads = dim3(1, 1, 1)

!    istat = cudaHostRegister(c_loc(phi), nx * ny * nz *szreal, 1)


    ! Copy input data from CPU to GPU memory.
!    if      (dimension_id .eq. 1) then
!        istat = cudaMemcpy(cDu1_dev, c_loc(cDu1), nx_div_coef1 * ny_div_coef1 * szreal, cudaMemcpyHostToDevice)
!        call cuda_err_check(istat)
!    else if (dimension_id .eq. 2) then
!        istat = cudaMemcpy(cDv2_dev, c_loc(cDv2), nx_div_coef2 * ny_div_coef2 * szreal, cudaMemcpyHostToDevice)
!        call cuda_err_check(istat)
!    else if (dimension_id .eq. 3) then
!        istat = cudaMemcpy(cDw3_dev, c_loc(cDw3), nx_div_coef3 * ny_div_coef3 * szreal, cudaMemcpyHostToDevice)
!        call cuda_err_check(istat)
!    end if
    


!    call cpu_time(start) 
    istat = cudaMemcpy(phi_div_dev, c_loc(phi), nx_dev * ny_dev * nz_dev * szreal, cudaMemcpyHostToDevice)
    call cuda_err_check(istat)
 !   call cpu_time(finish)

    istat = cudaMemcpy(div_dev, c_loc(div), nx_dev * ny_dev * nz_dev * szreal, cudaMemcpyHostToDevice)
    call cuda_err_check(istat)




    !print*, 'the copy time without the pinned memory', finish- start

! Creating CUDA Streams  !TEST -- wrapper is buggy -->  Maybe c_loc(stream1) can help
!    istat = cudaStreamCreate(stream1)
!    call cuda_err_check(istat)
!    istat = cudaStreamCreate(stream2)
!    call cuda_err_check(istat)
!    istat = cudaStreamCreate(stream3)
!    call cuda_err_check(istat)




!    istat = cudaMemcpy(div_dev, c_loc(div), nx * ny * nz * szreal, cudaMemcpyHostToDevice)
!    call cuda_err_check(istat)

#ifdef CUDAVERSION9

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
    istat = cudaSetupScalarArgument(c_loc(nx_dev), szint, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(nx_dev)
    istat = cudaSetupScalarArgument(c_loc(ny_dev), szint, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(ny_dev)
    istat = cudaSetupScalarArgument(c_loc(nz_dev), szint, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(nz_dev)
    if      (dimension_id .eq. 1) then
    istat = cudaSetupScalarArgument(c_loc(nx_div_coef1), szint, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(nx_div_coef1)
    istat = cudaSetupScalarArgument(c_loc(ny_div_coef1), szint, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(ny_div_coef1)
    else if      (dimension_id .eq. 2) then
    istat = cudaSetupScalarArgument(c_loc(nx_div_coef2), szint, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(nx_div_coef2)
    istat = cudaSetupScalarArgument(c_loc(ny_div_coef2), szint, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(ny_div_coef2)
    else if      (dimension_id .eq. 3) then
    istat = cudaSetupScalarArgument(c_loc(nx_div_coef3), szint, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(nx_div_coef3)
    istat = cudaSetupScalarArgument(c_loc(ny_div_coef3), szint, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(ny_div_coef3)
    endif
   ! Previous arguments were 4 bytes in size, total offset is 12 bytes.
    ! Since the next 3 arguments are 8 bytes each and have to be aligned
    ! by 8 bytes boundary, it's necessary to insert a 4-byte spacing here
    ! (from 12 to 16).
!    offset = offset + 4
    istat = cudaSetupArrayArgument(phi_div_dev, szptr, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(phi_div_dev)

    istat = cudaSetupArrayArgument(div_dev, szptr, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(div_dev)


    if      (dimension_id .eq. 1) then
        istat = cudaSetupArrayArgument(cDu1_dev, szptr, offset)
        call cuda_err_check(istat)
        offset = offset + sizeof(cDu1_dev)
    else if (dimension_id .eq. 2) then
        istat = cudaSetupArrayArgument(cDv2_dev, szptr, offset)
        call cuda_err_check(istat)
        offset = offset + sizeof(cDv2_dev)
    else if (dimension_id .eq. 3) then
        istat = cudaSetupArrayArgument(cDw3_dev, szptr, offset)
        call cuda_err_check(istat)
        offset = offset + sizeof(cDw3_dev)
    end if




!    istat = cudaSetupArrayArgument(div_dev, szptr, offset)
!    call cuda_err_check(istat)

    ! Finally, launch CUDA kernel after we configured everything.
    ! Kernel is identified by C-string with its name (must be
    ! undecorated, i.e. with extern "C")


!     call MPI_Barrier(MPI_COMM_WORLD, merror)
!     t1 =   MPI_Wtime()

    istat = cudaLaunch('divergence_kernel' // c_null_char)
    call cuda_err_check(istat)
!     call MPI_Barrier(MPI_COMM_WORLD, merror)
!     t2 =   MPI_Wtime()
!     print*, 'compute time divergence', t2-t1


#endif 
    if (dimension_id == 1) then
       call divergence_launcher((dimension_id), (nx_dev), (ny_dev),(nz_dev), (nx_div_coef1), (ny_div_coef1), (phi_div_dev), (div_dev), cDu1_dev)
    elseif (dimension_id ==2) then
       call divergence_launcher((dimension_id), (nx_dev), (ny_dev),(nz_dev), (nx_div_coef2), (ny_div_coef2), (phi_div_dev), (div_dev), cDv2_dev)
    elseif (dimension_id ==3) then
       call divergence_launcher((dimension_id), (nx_dev), (ny_dev),(nz_dev), (nx_div_coef3), (ny_div_coef3), (phi_div_dev), (div_dev), (cDw3_dev))
    endif


    ! Copy back results from GPU to CPU memory.
    istat = cudaMemcpy(c_loc(div), div_dev, nx_dev * ny_dev * nz_dev * szreal, cudaMemcpyDeviceToHost)
!    istat = cudaMemcpy(c_loc(cDtest), cDu1_dev, nx_div_coef* ny_div_coef* szreal ,cudaMemcpyDeviceToHost)
    call cuda_err_check(istat)




! Testing the coefficent copy
!    nmatch = 0
!    do i =-3, 2
!       do j=0, N1
!         if (abs(cDu1(i,j)- cDtest(i,j)).le. 0.000000000001) nmatch = nmatch + 1
!       enddo
!    enddo
!    print*, nmatch, 'match number of copy test'
    ! Free allocated GPU memory.
!    istat = cudaFree(div_dev)
!    call cuda_err_check(istat)
!    istat = cudaFree(phi_div_dev)
!    call cuda_err_check(istat)

!    istat = cudaFree(div_dev_copy)
!    call cuda_err_check(istat)


!    if      (dimension_id .eq. 1) then
!        istat = cudaFree(cDu1_dev)
!        call cuda_err_check(istat)
!    else if (dimension_id .eq. 2) then
!        istat = cudaFree(cDv2_dev)
!        call cuda_err_check(istat)
!    else if (dimension_id .eq. 3) then
!        istat = cudaFree(cDw3_dev)
!        call cuda_err_check(istat)
!    end if



end subroutine device_divergence

