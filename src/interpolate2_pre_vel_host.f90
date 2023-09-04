!**************************************************************************************************                                    // !* by Hadi Zolfaghari, ARTORG Center, University of Bern, (zolfaghari.haadi@gmail.com)            *      
!* October 2015 - March 2020                                                                      *            
!* Modified by Hadi Zolfaghari, University of Cambridge (hz382@damtp.cam.ac.uk)                   *      
!* April 2020 -                                                                                   *            
!**************************************************************************************************

subroutine device_interpolate(dimension_id,phi_mat,inter_mat)

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


   integer, target, intent(in)  ::  dimension_id  ! the dimension_id of interpolate --> corresponding to 'm=1,2,3' in IMPACT
   

   REAL(8)   , target, INTENT(inout) ::  phi_mat(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))!< ? 
   REAL(8)   , target, INTENT(  out) ::  inter_mat(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))!< interpolate operator



  ! Declare the size of phi and other arguments passing to interpolate -- ref : mod_diff

   integer   ::   phi_size(3)
   integer, target     :: nx, ny, nz ! the inter and phi size components


!   integer   ::   cIpu_size(2)
!   integer   ::   cIpv_size(2)
!   integer   ::   cIpw_size(2)
!
!
!  integer, target     :: nxIpu, nyIpu !the cIpu size components

!  type(c_ptr) :: inter_pv_dev, phi_inter_pv_dev, cIpu_dev, cIpv_dev, cIpw_dev
    
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

!   CALL exchange(dimension_id,0,phi) 



!  szreal = sizeof(inter_mat(1, 1, 1))
!  szint = sizeof(nx)
!  szptr = sizeof(d_x)
 
 
 ! We repeat GPU logic twice to show that the first kernel invocation includes
  ! the time of GPU device initialization, which is done one time per
  ! application run.

!    call cpu_time(start)

    if (GPU_verbose == 1) print*, 'interpolation on GPU in direction',' ',':',' ', dimension_id

    ! Get the size of arrays passing to interpolate subroutine
!    phi_size = shape(phi_mat)
!
!    nx = phi_size(1)  
!    ny = phi_size(2) 
!    nz = phi_size(3) 
    

!    if (dimension_id .eq. 1)  cIpu_size = shape(cIpu)
!    if (dimension_id .eq. 2)  cIpv_size = shape(cIpv)
!    if (dimension_id .eq. 3)  cIpw_size = shape(cIpw)
!
!
!    if      (dimension_id .eq. 1) then
!        nxIpu1 = cIpu_size(1)
!        nyIpu1 = cIpu_size(2)
!    else if (dimension_id .eq. 2) then
!        nxIpu2 = cIpv_size(1)
!        nyIpu2 = cIpv_size(2)
!    else if (dimension_id .eq. 3) then
!        nxIpu3 = cIpw_size(1)
!        nyIpu3 = cIpw_size(2)
!    endif


!    print*, nx, ny, nz, nxIpu, nyIpu, 'nx, ny, nz, nxIpu, nyIpu', 'this is in interpolate'
    ! Allocate memory on GPU.

!    if      (dimension_id .eq. 1) then
!        istat = cudaMalloc(cIpu_dev, nxIpu1 * nyIpu1 * szreal)
!        call cuda_err_check(istat)
!    else if (dimension_id .eq. 2) then
!        istat = cudaMalloc(cIpv_dev, nxIpu2 * nyIpu2 * szreal)
!        call cuda_err_check(istat)
!    else if (dimension_id .eq. 3) then
!        istat = cudaMalloc(cIpw_dev, nxIpu3 * nyIpu3 * szreal)
!        call cuda_err_check(istat)  
!    end if


!    istat = cudaMalloc(phi_inter_pv_dev, nx_dev * ny_dev * nz_dev * szreal)
!    call cuda_err_check(istat)
!    istat = cudaMalloc(inter_pv_dev, nx_dev * ny_dev * nz_dev *  szreal)
!    call cuda_err_check(istat)

	! XXX Several blocks of code below are equivalent to sum_kernel<<<...>>> in C.


!    print*, nx/8, ny/8, nz/8, 'block dimmension on gpu'

    ! Create CUDA compute grid configuration.
    blocks = dim3((nx_dev-8)/8, (ny_dev-8)/8, (nz_dev-8)/8)
    if (mod(nx_dev, szblock) .ne. 0) blocks%x = blocks%x + 1
    threads = dim3((nx_dev-8)/blocks%x, (ny_dev-8)/blocks%y, (nz_dev-8)/blocks%z)

!!    blocks = dim3(1, 1, 1)
!!    threads = dim3(1, 1, 1)
    ! Copy input data from CPU to GPU memory.
!    if      (dimension_id .eq. 1) then
!        istat = cudaMemcpy(cIpu_dev, c_loc(cIpu), nxIpu1 * nyIpu1 * szreal, cudaMemcpyHostToDevice)
!        call cuda_err_check(istat)
!    else if (dimension_id .eq. 2) then
!        istat = cudaMemcpy(cIpv_dev, c_loc(cIpv), nxIpu2 * nyIpu2 * szreal, cudaMemcpyHostToDevice)
!        call cuda_err_check(istat)
!    else if (dimension_id .eq. 3) then
!        istat = cudaMemcpy(cIpw_dev, c_loc(cIpw), nxIpu3 * nyIpu3 * szreal, cudaMemcpyHostToDevice)
!        call cuda_err_check(istat)
!    end if
    



    istat = cudaMemcpy(phi_inter_pv_dev, c_loc(phi_mat), nx_dev * ny_dev * nz_dev * szreal, cudaMemcpyHostToDevice)
    call cuda_err_check(istat)

! Creating CUDA Streams  !TEST -- wrapper is buggy -->  Maybe c_loc(stream1) can help
!    istat = cudaStreamCreate(stream1)
!    call cuda_err_check(istat)
!    istat = cudaStreamCreate(stream2)
!    call cuda_err_check(istat)
!    istat = cudaStreamCreate(stream3)
!    call cuda_err_check(istat)




!    istat = cudaMemcpy(inter_pv_dev, c_loc(inter), nx * ny * nz * szreal, cudaMemcpyHostToDevice)
!    call cuda_err_check(istat)


#ifdef CUDAVERSION9

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
    istat = cudaSetupScalarArgument(c_loc(nxIpu1), szint, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(nxIpu1)
    istat = cudaSetupScalarArgument(c_loc(nyIpu1), szint, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(nyIpu1)
    else if      (dimension_id .eq. 2) then
    istat = cudaSetupScalarArgument(c_loc(nxIpu2), szint, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(nxIpu2)
    istat = cudaSetupScalarArgument(c_loc(nyIpu2), szint, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(nyIpu2)
    else if      (dimension_id .eq. 3) then
    istat = cudaSetupScalarArgument(c_loc(nxIpu3), szint, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(nxIpu3)
    istat = cudaSetupScalarArgument(c_loc(nyIpu3), szint, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(nyIpu3)
    end if


   ! Previous arguments were 4 bytes in size, total offset is 12 bytes.
    ! Since the next 3 arguments are 8 bytes each and have to be aligned
    ! by 8 bytes boundary, it's necessary to insert a 4-byte spacing here
    ! (from 12 to 16).
!    offset = offset + 4
    istat = cudaSetupArrayArgument(phi_inter_pv_dev, szptr, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(phi_inter_pv_dev)
    if      (dimension_id .eq. 1) then
        istat = cudaSetupArrayArgument(cIpu_dev, szptr, offset)
        call cuda_err_check(istat)
        offset = offset + sizeof(cIpu_dev)
    else if (dimension_id .eq. 2) then
        istat = cudaSetupArrayArgument(cIpv_dev, szptr, offset)
        call cuda_err_check(istat)
        offset = offset + sizeof(cIpv_dev)
    else if (dimension_id .eq. 3) then
        istat = cudaSetupArrayArgument(cIpw_dev, szptr, offset)
        call cuda_err_check(istat)
        offset = offset + sizeof(cIpw_dev)
    end if




    istat = cudaSetupArrayArgument(inter_pv_dev, szptr, offset)
    call cuda_err_check(istat)

    ! Finally, launch CUDA kernel after we configured everything.
    ! Kernel is identified by C-string with its name (must be
    ! undecorated, i.e. with extern "C")
    istat = cudaLaunch('interpolate_kernel' // c_null_char)
    call cuda_err_check(istat)

#endif

    if (dimension_id == 1) then
       call interpolate2_pre_vel_launcher((dimension_id), (nx_dev), (ny_dev),(nz_dev), nxIpu1, nyIpu1, phi_inter_pv_dev, cIpu_dev, inter_pv_dev)
    elseif (dimension_id ==2) then
       call interpolate2_pre_vel_launcher((dimension_id), (nx_dev), (ny_dev),(nz_dev), nxIpu2, nyIpu2, phi_inter_pv_dev, cIpv_dev, inter_pv_dev)
    elseif (dimension_id ==3) then
       call interpolate2_pre_vel_launcher((dimension_id), (nx_dev), (ny_dev),(nz_dev), nxIpu3, nyIpu3, phi_inter_pv_dev, cIpw_dev, inter_pv_dev)
    endif

    ! Copy back results from GPU to CPU memory.
    istat = cudaMemcpy(c_loc(inter_mat), inter_pv_dev, nx_dev * ny_dev * nz_dev * szreal, cudaMemcpyDeviceToHost)
!    istat = cudaMemcpy(c_loc(cIpu), cIpu_dev, nxIpu* nyIpu* szreal ,cudaMemcpyDeviceToHost)
    call cuda_err_check(istat)



    ! Free allocated GPU memory.
!    istat = cudaFree(inter_pv_dev)
!    call cuda_err_check(istat)
!    istat = cudaFree(phi_inter_pv_dev)
!    call cuda_err_check(istat)
!    if      (dimension_id .eq. 1) then
!        istat = cudaFree(cIpu_dev)
!        call cuda_err_check(istat)
!    else if (dimension_id .eq. 2) then
!        istat = cudaFree(cIpv_dev)
!        call cuda_err_check(istat)
!    else if (dimension_id .eq. 3) then
!        istat = cudaFree(cIpw_dev)
!        call cuda_err_check(istat)
!    end if


! Applying boundary conditions
! This part is on the CPU and !!TODO: has to be shared-memory optimized


end subroutine device_interpolate

