!**************************************************************************************************                                    // !* by Hadi Zolfaghari, ARTORG Center, University of Bern, (zolfaghari.haadi@gmail.com)            *      
!* October 2015 - March 2020                                                                      *            
!* Modified by Hadi Zolfaghari, University of Cambridge (hz382@damtp.cam.ac.uk)                   *      
!* April 2020 -                                                                                   *            
!**************************************************************************************************
subroutine device_Nonlinear_upwind(component_id, vel_component, nl_component)

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
   real(8),    target, intent(inout)    ::  nl_component(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))         

  ! Declare the size of phi and other arguments passing to stencil operator -- ref : mod_diff



   

  type(dim3) :: blocks, threads


!  type(cudaStream_t) :: stream1, stream2, stream3

  integer :: i, j, k, step
  real(8), volatile :: start, finish


!! Aplying the exchange !!TODO check this for every stencil operation !! On the CPU

!   CALL exchange(dimension_id,0,phi) 




 
 ! We repeat GPU logic twice to show that the first kernel invocation includes
  ! the time of GPU device initialization, which is done one time per
  ! application run.

!    call cpu_time(start)



    ! Get the size of arrays passing to stencil subroutine
!    phi_size = shape(vel_component)



    if (GPU_verbose == 1) print*, 'Advection on GPU for component',' ',':',' ', component_id
     


!    if (component_id .eq. 1) then
!          cS_size = shape(cu1)
!          nxs = cS_size(1)
!          nys = cS_size(2)
!          istat = cudamalloc(cs_dev, nxs * nys * szreal)
!          call cuda_err_check(istat)
!          istat = cudamalloc(cs2_dev, nxs * nys * szreal)
!          call cuda_err_check(istat)
!          istat = cudaMemcpy(cS_dev, c_loc(cNu1U), nxs * nys * szreal, cudaMemcpyHostToDevice)
!          call cuda_err_check(istat)
!          istat = cudaMemcpy(cS2_dev, c_loc(cNu1D), nxs * nys * szreal, cudaMemcpyHostToDevice)
!          call cuda_err_check(istat)
!   endif 
!
!    if (component_id .eq. 2) then
!          cS_size = shape(cp2)
!          nxs = cS_size(1)
!          nys = cS_size(2)
!          istat = cudaMalloc(cS_dev, nxs * nys * szreal)
!          call cuda_err_check(istat)
!          istat = cudaMalloc(cS2_dev, nxs * nys * szreal)
!          call cuda_err_check(istat)
!          istat = cudaMemcpy(cS_dev, c_loc(cNp2U), nxs * nys * szreal, cudaMemcpyHostToDevice)
!          call cuda_err_check(istat)
!          istat = cudaMemcpy(cS2_dev, c_loc(cNp2D), nxs * nys * szreal, cudaMemcpyHostToDevice)
!          call cuda_err_check(istat)
!    endif
!
!    if (component_id .eq. 3) then
!          cS_size = shape(cp3)
!          nxs = cS_size(1)
!          nys = cS_size(2)
!          istat = cudaMalloc(cS_dev, nxs * nys * szreal)
!          call cuda_err_check(istat)
!          istat = cudaMalloc(cS2_dev, nxs * nys * szreal)
!          call cuda_err_check(istat)
!          istat = cudaMemcpy(cS_dev, c_loc(cNp3U), nxs * nys * szreal, cudaMemcpyHostToDevice)
!          call cuda_err_check(istat) 
!          istat = cudaMemcpy(cS2_dev, c_loc(cNp3D), nxs * nys * szreal, cudaMemcpyHostToDevice)
!          call cuda_err_check(istat)
!    endif
!
!    if (component_id .eq. 4) then
!          cS_size = shape(cp1)
!          nxs = cS_size(1)
!          nys = cS_size(2) 
!          istat = cudaMalloc(cS_dev, nxs * nys * szreal)
!          call cuda_err_check(istat)
!          istat = cudaMalloc(cS2_dev, nxs * nys * szreal)
!          call cuda_err_check(istat)
!          istat = cudaMemcpy(cS_dev, c_loc(cNp1U), nxs * nys * szreal, cudaMemcpyHostToDevice)
!          call cuda_err_check(istat)
!          istat = cudaMemcpy(cS2_dev, c_loc(cNp1D), nxs * nys * szreal, cudaMemcpyHostToDevice)
!          call cuda_err_check(istat)
!    endif
!
!    if (component_id .eq. 5) then
!          cS_size = shape(cv2)
!          nxs = cS_size(1)
!          nys = cS_size(2) 
!          istat = cudaMalloc(cS_dev, nxs * nys * szreal)
!          call cuda_err_check(istat)
!          istat = cudaMalloc(cS2_dev, nxs * nys * szreal)
!          call cuda_err_check(istat)
!          istat = cudaMemcpy(cS_dev, c_loc(cNv2U), nxs * nys * szreal, cudaMemcpyHostToDevice)
!          call cuda_err_check(istat) 
!          istat = cudaMemcpy(cS2_dev, c_loc(cNv2D), nxs * nys * szreal, cudaMemcpyHostToDevice)
!          call cuda_err_check(istat)
!    endif
!
!    if (component_id .eq. 6) then  
!          cS_size = shape(cp3)
!          nxs = cS_size(1)
!          nys = cS_size(2)
!          istat = cudaMalloc(cS_dev, nxs * nys * szreal)
!          call cuda_err_check(istat)
!          istat = cudaMalloc(cS2_dev, nxs * nys * szreal)
!          call cuda_err_check(istat)
!          istat = cudaMemcpy(cS_dev, c_loc(cNp3U), nxs * nys * szreal, cudaMemcpyHostToDevice)
!          call cuda_err_check(istat)
!          istat = cudaMemcpy(cS2_dev, c_loc(cNp3D), nxs * nys * szreal, cudaMemcpyHostToDevice)
!          call cuda_err_check(istat)
!    endif
!
!    if (component_id .eq. 7) then
!          cS_size = shape(cp1)
!          nxs = cS_size(1)
!          nys = cS_size(2) 
!          istat = cudaMalloc(cS_dev, nxs * nys * szreal)
!          call cuda_err_check(istat)
!          istat = cudaMalloc(cS2_dev, nxs * nys * szreal)
!          call cuda_err_check(istat)
!          istat = cudaMemcpy(cS_dev, c_loc(cNp1U), nxs * nys * szreal, cudaMemcpyHostToDevice)
!          call cuda_err_check(istat)
!          istat = cudaMemcpy(cS2_dev, c_loc(cNp1D), nxs * nys * szreal, cudaMemcpyHostToDevice)
!          call cuda_err_check(istat)
!    endif
!
!    if (component_id .eq. 8) then
!          cS_size = shape(cp2)
!          nxs = cS_size(1)
!          nys = cS_size(2)
!          istat = cudaMalloc(cS_dev, nxs * nys * szreal)
!          call cuda_err_check(istat)
!          istat = cudaMalloc(cS2_dev, nxs * nys * szreal)
!          call cuda_err_check(istat)
!          istat = cudaMemcpy(cS_dev, c_loc(cNp2U), nxs * nys * szreal, cudaMemcpyHostToDevice)
!          call cuda_err_check(istat)
!          istat = cudaMemcpy(cS2_dev, c_loc(cNp2D), nxs * nys * szreal, cudaMemcpyHostToDevice)
!          call cuda_err_check(istat)
!    endif
!
!    if (component_id .eq. 9) then
!          cS_size = shape(cw3)
!          nxs = cS_size(1)
!          nys = cS_size(2) 
!          istat = cudaMalloc(cS_dev, nxs * nys * szreal)
!          call cuda_err_check(istat)
!          istat = cudaMalloc(cS2_dev, nxs * nys * szreal)
!          call cuda_err_check(istat)
!          istat = cudaMemcpy(cS_dev, c_loc(cNw3U), nxs * nys * szreal, cudaMemcpyHostToDevice)
!          call cuda_err_check(istat)
!          istat = cudaMemcpy(cS2_dev, c_loc(cNw3D), nxs * nys * szreal, cudaMemcpyHostToDevice)
!          call cuda_err_check(istat)
!    endif


!    nx = phi_size(1)  
!    ny = phi_size(2) 
!    nz = phi_size(3) 


    !print*, nx, ny, nz, nxs, nys, 'nx, ny, nz, nxs, nys', 'this is in output', rank
    ! Allocate memory on GPU.

!    istat = cudaMalloc(pp_dev, nx * ny * nz * szreal)
!    call cuda_err_check(istat)
!    ! phi and output come in size size -- ref to mod_diff
!    istat = cudaMalloc(nl_dev, nx * ny * nz *  szreal)
!    call cuda_err_check(istat)
!    istat = cudaMalloc(vel_dev, nx * ny * nz *  szreal)
!    call cuda_err_check(istat)

! XXX Several blocks of code below are equivalent to sum_kernel<<<...>>> in C.

#ifdef CUDAVERSION9

    ! Create CUDA compute grid configuration.
    blocks = dim3((nx_dev-8)/8, (ny_dev-8)/8, (nz_dev-8)/8)
    if (mod(nx_dev, szblock) .ne. 0) blocks%x = blocks%x + 1
    threads = dim3((nx_dev-8)/blocks%x, (ny_dev-8)/blocks%y, (nz_dev-8)/blocks%z)

#endif

    ! Copy input data from CPU to GPU memory.
   

    istat = cudaMemcpy(pp_dev, c_loc(pp), nx_dev * ny_dev * nz_dev * szreal, cudaMemcpyHostToDevice)
    call cuda_err_check(istat)
    !!  TODO velocity should be transferred component-wise
    istat = cudaMemcpy(vel_dev, c_loc(vel_component), nx_dev * ny_dev * nz_dev  * szreal, cudaMemcpyHostToDevice)
    call cuda_err_check(istat)
    istat = cudaMemcpy(nl_dev, c_loc(nl_component), nx_dev * ny_dev * nz_dev * szreal, cudaMemcpyHostToDevice)
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
    istat = cudaSetupScalarArgument(c_loc(nxn1), szint, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(nxn1)
    istat = cudaSetupScalarArgument(c_loc(nyn1), szint, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(nyn1)
    endif

    if (component_id .eq. 2) then
    istat = cudaSetupScalarArgument(c_loc(nxn2), szint, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(nxn2)
    istat = cudaSetupScalarArgument(c_loc(nyn2), szint, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(nyn2)
    endif

    if (component_id .eq. 3) then
    istat = cudaSetupScalarArgument(c_loc(nxn3), szint, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(nxn3)
    istat = cudaSetupScalarArgument(c_loc(nyn3), szint, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(nyn3)
    endif

    if (component_id .eq. 4) then
    istat = cudaSetupScalarArgument(c_loc(nxn4), szint, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(nxn4)
    istat = cudaSetupScalarArgument(c_loc(nyn4), szint, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(nyn4)
    endif

    if (component_id .eq. 5) then
    istat = cudaSetupScalarArgument(c_loc(nxn5), szint, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(nxn5)
    istat = cudaSetupScalarArgument(c_loc(nyn5), szint, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(nyn5)
    endif
    if (component_id .eq. 6) then
    istat = cudaSetupScalarArgument(c_loc(nxn3), szint, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(nxn3)
    istat = cudaSetupScalarArgument(c_loc(nyn3), szint, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(nyn3)
    endif

    if (component_id .eq. 7) then
    istat = cudaSetupScalarArgument(c_loc(nxn4), szint, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(nxn4)
    istat = cudaSetupScalarArgument(c_loc(nyn4), szint, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(nyn4)
    endif


    if (component_id .eq. 8) then
    istat = cudaSetupScalarArgument(c_loc(nxn2), szint, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(nxn2)
    istat = cudaSetupScalarArgument(c_loc(nyn2), szint, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(nyn2)
    endif

    if (component_id .eq. 9) then
    istat = cudaSetupScalarArgument(c_loc(nxn6), szint, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(nxn6)
    istat = cudaSetupScalarArgument(c_loc(nyn6), szint, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(nyn6)
    endif





   ! Previous arguments were 4 bytes in size, total offset is 12 bytes.
    ! Since the next 3 arguments are 8 bytes each and have to be aligned
    ! by 8 bytes boundary, it's necessary to insert a 4-byte spacing here
    ! (from 12 to 16).
   !  offset = offset + 4
    istat = cudaSetupArrayArgument(pp_dev, szptr, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(pp_dev)

    istat = cudaSetupArrayArgument(vel_dev, szptr, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(vel_dev)

    istat = cudaSetupArrayArgument(nl_dev, szptr, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(nl_dev)


    if (component_id .eq. 1) then
    istat = cudaSetupArrayArgument(cN1_dev, szptr, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(cN1_dev)

    istat = cudaSetupArrayArgument(cN21_dev, szptr, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(cN21_dev)
    endif

    if (component_id .eq. 2) then
    istat = cudaSetupArrayArgument(cN2_dev, szptr, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(cN2_dev)

    istat = cudaSetupArrayArgument(cN22_dev, szptr, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(cN22_dev)
    endif

    if (component_id .eq. 3) then
    istat = cudaSetupArrayArgument(cN3_dev, szptr, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(cN3_dev)

    istat = cudaSetupArrayArgument(cN23_dev, szptr, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(cN23_dev)
    endif

    if (component_id .eq. 4) then
    istat = cudaSetupArrayArgument(cN4_dev, szptr, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(cN4_dev)

    istat = cudaSetupArrayArgument(cN24_dev, szptr, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(cN24_dev)
    endif

    if (component_id .eq. 5) then
    istat = cudaSetupArrayArgument(cN5_dev, szptr, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(cN5_dev)

    istat = cudaSetupArrayArgument(cN25_dev, szptr, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(cN25_dev)
    endif

    if (component_id .eq. 6) then
    istat = cudaSetupArrayArgument(cN3_dev, szptr, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(cN3_dev)

    istat = cudaSetupArrayArgument(cN23_dev, szptr, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(cN23_dev)
    endif

    if (component_id .eq. 7) then
    istat = cudaSetupArrayArgument(cN4_dev, szptr, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(cN4_dev)

    istat = cudaSetupArrayArgument(cN24_dev, szptr, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(cN24_dev)
    endif

    if (component_id .eq. 8) then
    istat = cudaSetupArrayArgument(cN2_dev, szptr, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(cN2_dev)

    istat = cudaSetupArrayArgument(cN22_dev, szptr, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(cN22_dev)
    endif

    if (component_id .eq. 9) then
    istat = cudaSetupArrayArgument(cN6_dev, szptr, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(cN6_dev)

    istat = cudaSetupArrayArgument(cN26_dev, szptr, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(cN26_dev)
    endif


    ! Finally, launch CUDA kernel after we configured everything.
    ! Kernel is identified by C-string with its name (must be
    ! undecorated, i.e. with extern "C")
    istat = cudaLaunch('nonlinear_upwind_kernel' // c_null_char)
    call cuda_err_check(istat)


#endif


    if (component_id == 1) then
       call nonlinear_upwind_launcher((component_id), (nx_dev), (ny_dev),(nz_dev), (nxn1), (nyn1), (pp_dev), (vel_dev), (nl_dev), cN1_dev, cN21_dev)
    else if (component_id == 2) then
       call nonlinear_upwind_launcher((component_id), (nx_dev), (ny_dev),(nz_dev), (nxn2), (nyn2), (pp_dev), (vel_dev), (nl_dev), cN2_dev, cN22_dev)
    else if (component_id == 3) then
       call nonlinear_upwind_launcher((component_id), (nx_dev), (ny_dev),(nz_dev), (nxn3), (nyn3), (pp_dev), (vel_dev), (nl_dev), cN3_dev, cN23_dev)
    else if (component_id == 4) then
       call nonlinear_upwind_launcher((component_id), (nx_dev), (ny_dev),(nz_dev), (nxn4), (nyn4), (pp_dev), (vel_dev), (nl_dev), cN4_dev, cN24_dev)
    else if (component_id == 5) then
       call nonlinear_upwind_launcher((component_id), (nx_dev), (ny_dev),(nz_dev), (nxn5), (nyn5), (pp_dev), (vel_dev), (nl_dev), cN5_dev, cN25_dev)
    else if (component_id == 6) then
       call nonlinear_upwind_launcher((component_id), (nx_dev), (ny_dev),(nz_dev), (nxn3), (nyn3), (pp_dev), (vel_dev), (nl_dev), cN3_dev, cN23_dev)
    else if (component_id == 7) then
       call nonlinear_upwind_launcher((component_id), (nx_dev), (ny_dev),(nz_dev), (nxn4), (nyn4), (pp_dev), (vel_dev), (nl_dev), cN4_dev, cN24_dev)
    else if (component_id == 8) then
       call nonlinear_upwind_launcher((component_id), (nx_dev), (ny_dev),(nz_dev), (nxn2), (nyn2), (pp_dev), (vel_dev), (nl_dev), cN2_dev, cN22_dev)
    else if (component_id == 9) then
       call nonlinear_upwind_launcher((component_id), (nx_dev), (ny_dev),(nz_dev), (nxn6), (nyn6), (pp_dev), (vel_dev), (nl_dev), cN6_dev, cN26_dev)
    endif





    ! Copy back results from GPU to CPU memory.
    istat = cudaMemcpy(c_loc(nl_component), nl_dev, nx_dev * ny_dev * nz_dev * szreal, cudaMemcpyDeviceToHost)
!    istat = cudaMemcpy(c_loc(cS1), cS1_dev, nxs* nys* szreal ,cudaMemcpyDeviceToHost)
    call cuda_err_check(istat)



    ! Free allocated GPU memory.
!    istat = cudaFree(nl_dev)
!    call cuda_err_check(istat)
!    istat = cudaFree(vel_dev)
!    call cuda_err_check(istat)
!    istat = cudaFree(pp_dev)
!    call cuda_err_check(istat)
!    istat = cudaFree(cS_dev)
!    call cuda_err_check(istat)
!    istat = cudaFree(cS2_dev)
!    call cuda_err_check(istat)
!


! Applying boundary conditions


end subroutine device_Nonlinear_upwind

