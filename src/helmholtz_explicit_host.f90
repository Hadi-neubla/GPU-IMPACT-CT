!**************************************************************************************************                                    // !* by Hadi Zolfaghari, ARTORG Center, University of Bern, (zolfaghari.haadi@gmail.com)            *      
!* October 2015 - March 2020                                                                      *            
!* Modified by Hadi Zolfaghari, University of Cambridge (hz382@damtp.cam.ac.uk)                   *      
!* April 2020 -                                                                                   *            
!**************************************************************************************************
subroutine device_helmholtz(component_id, vel_component, helm_component)

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


   integer, target, intent(in)    ::  component_id  ! row number in advection tensor u. del(u)



   real(8),    target, intent(in)       ::  vel_component(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))         
   real(8),    target, intent(inout)    ::  helm_component(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))         
!   real(8),    target ,intent(in)  ::       multLg 
  ! Declare the size of phi and other arguments passing to stencil operator -- ref : mod_diff

!   integer   ::   phi_size(3)
!   integer, target     :: nx, ny, nz ! the output and phi size components

!   integer   ::   cS1_size(2)
!   integer   ::   cS2_size(2)
!   integer   ::   cS3_size(2)
!   integer   ::   cS4_size(2)
!   integer   ::   cS5_size(2)
!   integer   ::   cS6_size(2)
!
!
!
!  integer, target     :: nxs1, nys1, nxs2, nys2, nxs3, nys3 !the cS1 size components
!  integer, target     :: nxs4, nys4, nxs5, nys5, nxs6, nys6 !the cS1 size components

!  type(c_ptr) :: helm_dev, vel_helm_dev, cS1_helm_dev, cS2_helm_dev, cS3_helm_dev
    
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

!  szreal = sizeof(nl(1, 1, 1, 1))
!  szint = sizeof(nx)
!  szptr = sizeof(helm_dev)
 
 
 ! We repeat GPU logic twice to show that the first kernel invocation includes
  ! the time of GPU device initialization, which is done one time per
  ! application run.

!    call cpu_time(start)

    ! Get the size of arrays passing to stencil subroutine
!    phi_size = shape(vel_component)
!    nx = phi_size(1)  
!    ny = phi_size(2) 
!    nz = phi_size(3) 



    if (GPU_verbose == 1) print*, 'On GPU Explicit Helmholtz Operation for ',' ',':',' ', component_id
    

!!    if (component_id .eq. 1) then
!          cS1_size = shape(cu11)
!          nxs1 = cS1_size(1)
!          nys1 = cS1_size(2)
!          cS2_size = shape(cp22)
!          nxs2 = cS2_size(1)
!          nys2 = cS2_size(2)
!          cS3_size = shape(cp33)
!          nxs3 = cS3_size(1)
!          nys3 = cS3_size(2)
!          istat = cudaMalloc(cS1_helm_dev, nxs1 * nys1 * szreal)
!          call cuda_err_check(istat)
!          istat = cudaMemcpy(cS1_helm_dev, c_loc(cu11), nxs1 * nys1 * szreal, cudaMemcpyHostToDevice)
!          call cuda_err_check(istat)
!          istat = cudaMalloc(cS2_helm_dev, nxs2 * nys2 * szreal)
!          call cuda_err_check(istat)
!          istat = cudaMemcpy(cS2_helm_dev, c_loc(cp22), nxs2 * nys2 * szreal, cudaMemcpyHostToDevice)
!          call cuda_err_check(istat)
!          istat = cudaMalloc(cS3_helm_dev, nxs3 * nys3 * szreal)
!          call cuda_err_check(istat)
!          istat = cudaMemcpy(cS3_helm_dev, c_loc(cp33), nxs3 * nys3 * szreal, cudaMemcpyHostToDevice)
!          call cuda_err_check(istat)
! ! endif 
!
! ! if (component_id .eq. 2) then
!!          cS1_size = shape(cp11)
!!          nxs1 = cS1_size(1)
!!          nys1 = cS1_size(2)
!          cS4_size = shape(cv22)
!          nxs4 = cS4_size(1)
!          nys4 = cS4_size(2)
!!          cS3_size = shape(cp33)
!!          nxs3 = cS3_size(1)
!!          nys3 = cS3_size(2)
!!          istat = cudaMalloc(cS1_helm_dev, nxs1 * nys1 * szreal)
!!          call cuda_err_check(istat)
!!          istat = cudaMemcpy(cS1_helm_dev, c_loc(cp11), nxs1 * nys1 * szreal, cudaMemcpyHostToDevice)
!!          call cuda_err_check(istat)
!          istat = cudaMalloc(cS4_helm_dev, nxs4 * nys4 * szreal)
!          call cuda_err_check(istat)
!          istat = cudaMemcpy(cS4_helm_dev, c_loc(cv22), nxs4 * nys4 * szreal, cudaMemcpyHostToDevice)
!          call cuda_err_check(istat)
! !         istat = cudaMalloc(cS3_helm_dev, nxs3 * nys3 * szreal)
! !         call cuda_err_check(istat)
! !         istat = cudaMemcpy(cS3_helm_dev, c_loc(cp33), nxs3 * nys3 * szreal, cudaMemcpyHostToDevice)
! !         call cuda_err_check(istat)
! !  endif
!
!  ! if (component_id .eq. 3) then
!
!          cS5_size = shape(cp11)
!          nxs5 = cS5_size(1)
!          nys5 = cS5_size(2)
! !         cS2_size = shape(cp22)
! !         nxs2 = cS2_size(1)
! !         nys2 = cS2_size(2)
!          cS6_size = shape(cw33)
!          nxs6 = cS6_size(1)
!          nys6 = cS6_size(2)
!          istat = cudaMalloc(cS5_helm_dev, nxs5 * nys5 * szreal)
!          call cuda_err_check(istat)
!          istat = cudaMemcpy(cS5_helm_dev, c_loc(cp11), nxs5 * nys5 * szreal, cudaMemcpyHostToDevice)
!          call cuda_err_check(istat)
!!          istat = cudaMalloc(cS2_helm_dev, nxs2 * nys2 * szreal)
!!          call cuda_err_check(istat)
!!          istat = cudaMemcpy(cS2_helm_dev, c_loc(cp22), nxs2 * nys2 * szreal, cudaMemcpyHostToDevice)
!!          call cuda_err_check(istat)
!          istat = cudaMalloc(cS6_helm_dev, nxs6 * nys6 * szreal)
!          call cuda_err_check(istat)
!          istat = cudaMemcpy(cS6_helm_dev, c_loc(cw33), nxs6 * nys6 * szreal, cudaMemcpyHostToDevice)
!          call cuda_err_check(istat) 

  ! endif



!    print*, nx, ny, nz, nxs, nys, 'nx, ny, nz, nxs, nys', 'this is in output'
    ! Allocate memory on GPU.

   ! phi and output come in size size -- ref to mod_diff
!    istat = cudaMalloc(helm_dev, nx_dev * ny_dev * nz_dev *  szreal)
!    call cuda_err_check(istat)
!
!    istat = cudaMalloc(vel_helm_dev, nx_dev * ny_dev * nz_dev *  szreal)
!    call cuda_err_check(istat)

! XXX Several blocks of code below are equivalent to sum_kernel<<<...>>> in C.


!    print*, nx/8, ny/8, nz/8, 'block dimmension on gpu'

    ! Create CUDA compute grid configuration.
#ifdef CUDAVERSION9
    blocks = dim3((nx_dev-8)/8, (ny_dev-8)/8, (nz_dev-8)/8)
    if (mod(nx_dev, szblock) .ne. 0) blocks%x = blocks%x + 1
    threads = dim3((nx_dev-8)/blocks%x, (ny_dev-8)/blocks%y, (nz_dev-8)/blocks%z)
#endif

!!    blocks = dim3(1, 1, 1)
!!    threads = dim3(1, 1, 1)
    ! Copy input data from CPU to GPU memory.
   

   !!  TODO velocity should be transferred component-wise

    if (component_id .eq. 1) then
    istat = cudaMemcpy(vel1_dev, c_loc(vel_component), nx_dev * ny_dev * nz_dev  * szreal, cudaMemcpyHostToDevice)
    call cuda_err_check(istat)
    else  if (component_id .eq. 2) then
    istat = cudaMemcpy(vel2_dev, c_loc(vel_component), nx_dev * ny_dev * nz_dev  * szreal, cudaMemcpyHostToDevice)
    call cuda_err_check(istat)
    else  if (component_id .eq. 3) then
    istat = cudaMemcpy(vel3_dev, c_loc(vel_component), nx_dev * ny_dev * nz_dev  * szreal, cudaMemcpyHostToDevice)
    call cuda_err_check(istat)
    end if
  
    istat = cudaMemcpy(helm_dev, c_loc(helm_component), nx_dev * ny_dev * nz_dev * szreal, cudaMemcpyHostToDevice)
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


#ifdef CUDAVERSION9

!    print*, threads, blocks
!    print*, nx, ny, nz
   ! Setup CUDA compute grid configuration.
    istat = cudaConfigureCall(blocks, threads, int8(0), c_null_ptr)
    call cuda_err_check(istat)

!!! Stream 1 : the cuda API break down

    ! Setup CUDA kernel arguments (individually)
    offset = 0
    istat = cudaSetupScalarArgument(c_loc(component_id), szint, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(component_id)
    istat = cudaSetupScalarArgument(c_loc(nx_dev), szint, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(nx_dev)
    istat = cudaSetupScalarArgument(c_loc(ny_dev), szint, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(ny_dev)
    istat = cudaSetupScalarArgument(c_loc(nz_dev), szint, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(nz_dev)

    if (component_id .eq. 1) then
    istat = cudaSetupScalarArgument(c_loc(nxs1), szint, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(nxs1)
    istat = cudaSetupScalarArgument(c_loc(nxs2), szint, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(nxs2)
    istat = cudaSetupScalarArgument(c_loc(nxs3), szint, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(nxs3)
    else if (component_id .eq. 2) then
    istat = cudaSetupScalarArgument(c_loc(nxs1), szint, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(nxs1)
    istat = cudaSetupScalarArgument(c_loc(nxs4), szint, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(nxs4)
    istat = cudaSetupScalarArgument(c_loc(nxs3), szint, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(nxs3)
    else if (component_id .eq. 3) then
    istat = cudaSetupScalarArgument(c_loc(nxs5), szint, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(nxs5)
    istat = cudaSetupScalarArgument(c_loc(nxs2), szint, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(nxs2)
    istat = cudaSetupScalarArgument(c_loc(nxs6), szint, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(nxs6)
    endif



   ! Previous arguments were 4 bytes in size, total offset is 12 bytes.
    ! Since the next 3 arguments are 8 bytes each and have to be aligned
    ! by 8 bytes boundary, it's necessary to insert a 4-byte spacing here
    ! (from 12 to 16).
    offset = offset + 4
    istat = cudaSetupScalarArgument(c_loc(multL), szreal, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(multL)
    if (component_id .eq. 1) then
    istat = cudaSetupArrayArgument(vel1_dev, szptr, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(vel1_dev)
    else if (component_id .eq. 2) then
    istat = cudaSetupArrayArgument(vel2_dev, szptr, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(vel2_dev)
    else if (component_id .eq. 3) then
    istat = cudaSetupArrayArgument(vel3_dev, szptr, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(vel3_dev)
    endif



    istat = cudaSetupArrayArgument(helm_dev, szptr, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(helm_dev)
    if (component_id .eq. 1) then
    istat = cudaSetupArrayArgument(cS1_helm_dev, szptr, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(cS1_helm_dev)

    istat = cudaSetupArrayArgument(cS2_helm_dev, szptr, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(cS2_helm_dev)

    istat = cudaSetupArrayArgument(cS3_helm_dev, szptr, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(cS3_helm_dev)
    else if (component_id .eq. 2) then
    istat = cudaSetupArrayArgument(cS5_helm_dev, szptr, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(cS5_helm_dev)

    istat = cudaSetupArrayArgument(cS4_helm_dev, szptr, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(cS4_helm_dev)

    istat = cudaSetupArrayArgument(cS3_helm_dev, szptr, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(cS3_helm_dev)
    else  if (component_id .eq. 3) then 
    istat = cudaSetupArrayArgument(cS5_helm_dev, szptr, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(cS5_helm_dev)

    istat = cudaSetupArrayArgument(cS2_helm_dev, szptr, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(cS2_helm_dev)

    istat = cudaSetupArrayArgument(cS6_helm_dev, szptr, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(cS6_helm_dev)

    end if

    ! Finally, launch CUDA kernel after we configured everything.
    ! Kernel is identified by C-string with its name (must be
    ! undecorated, i.e. with extern "C")
   !  call MPI_Barrier(MPI_COMM_WORLD, merror)
   !  t1 =   MPI_Wtime()


    istat = cudaLaunch('helmholtz_kernel' // c_null_char)
!    call cuda_err_check(istat)

 !    call MPI_Barrier(MPI_COMM_WORLD, merror)
  !   t2 =   MPI_Wtime()

#endif


    if (component_id == 2) then
       call helmholtz_explicit_launcher((component_id), (nx_dev), (ny_dev),(nz_dev), (nxs5), (nxs4), nxs3, multL, (vel2_dev), (helm_dev), (cS5_helm_dev), cS4_helm_dev, cS3_helm_dev)
     else if (component_id == 1) then
       call helmholtz_explicit_launcher((component_id), (nx_dev), (ny_dev),(nz_dev), (nxs1), (nxs2), nxs3, multL, (vel1_dev), (helm_dev), (cS1_helm_dev), cS2_helm_dev, cS3_helm_dev)
    else if (component_id == 3) then
       call helmholtz_explicit_launcher((component_id), (nx_dev), (ny_dev),(nz_dev), (nxs5), (nxs2), nxs6, multL, (vel3_dev), (helm_dev), (cS5_helm_dev), cS2_helm_dev, cS6_helm_dev)
    endif




!     print*,  'compute time helmholtz' , t2-t1

    ! Copy back results from GPU to CPU memory.
    istat = cudaMemcpy(c_loc(helm_component), helm_dev, nx_dev * ny_dev * nz_dev * szreal, cudaMemcpyDeviceToHost)

  !  istat = cudaMemcpy(c_loc(tempor), c_loc(multL), szreal, cudaMemcpyDeviceToHost)
    call cuda_err_check(istat)

  !  print*, tempor

    ! Free allocated GPU memory.
!    istat = cudaFree(helm_dev)
!    call cuda_err_check(istat)
!    istat = cudaFree(vel_helm_dev)
!    call cuda_err_check(istat)
!    istat = cudaFree(cS1_helm_dev)
!    call cuda_err_check(istat)
!    istat = cudaFree(cS2_helm_dev)
!    call cuda_err_check(istat)
!    istat = cudaFree(cS3_helm_dev)
!    call cuda_err_check(istat)
!    istat = cudaFree(cS4_helm_dev)
!    call cuda_err_check(istat)
!    istat = cudaFree(cS5_helm_dev)
!    call cuda_err_check(istat)
!    istat = cudaFree(cS6_helm_dev)
!    call cuda_err_check(istat)

! Applying boundary conditions


end subroutine device_helmholtz

