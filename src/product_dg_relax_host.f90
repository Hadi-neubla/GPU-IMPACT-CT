!**************************************************************************************************                                           
!* by Hadi Zolfaghari, ARTORG Center, University of Bern, (zolfaghari.haadi@gmail.com)            *      
!* October 2015 - March 2020                                                                      *            
!**************************************************************************************************
subroutine device_product_dg_relax(grid_id,rel_g, comp_g)

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
   

   REAL(8)   , target, INTENT(inout) ::  rel_g (b1L:(NN(1,grid_id)+b1U),b2L:(NN(2,grid_id)+b2U),b3L:(NN(3,grid_id)+b3U))!< ? 
   REAL(8)   , target, INTENT(  out) ::  comp_g(b1L:(NN(1,grid_id)+b1U),b2L:(NN(2,grid_id)+b2U),b3L:(NN(3,grid_id)+b3U))!< gradient operator




  ! Declare the size of rel and other arguments passing to jacobi relaxation -- ref : mod_laplace

   integer   ::   rel_size_g(3)



   integer   ::   cdg1_size_g(3)   
   integer   ::   cdg2_size_g(3)
   integer   ::   cdg3_size_g(3)



  ! Declare sincos kernel arguments with "target" -- for c_loc
  integer, target     :: nx_dev_g, ny_dev_g, nz_dev_g ! the rel, bb and comp size components

  integer, target     :: nx_cdg1_g, ny_cdg1_g, nz_cdg1_g ! the grad and phi size components
  integer, target     :: nx_cdg2_g, ny_cdg2_g, nz_cdg2_g ! the grad and phi size components
  integer, target     :: nx_cdg3_g, ny_cdg3_g, nz_cdg3_g ! the grad and phi size components

!  real(8), target, allocatable, dimension(:,:,:) :: h_x, h_y, h_xy1, h_xy2
!  real(8), target, dimension(b1L:(N1+b1U), b2L:(N2+b2U), b3L:(N3+b3U)) :: grad 
!  real(8), target, allocatable, dimension(:,:,:) :: phi_h, cGp1_h

!  type(c_ptr) :: d_x, d_y, d_xy
  type(c_ptr) :: rel_dev2, comp_dev2, cdg1_dev2,  cdg2_dev2, cdg3_dev2
    
  integer :: sum_host
  integer(c_int) :: istat

  integer, parameter :: szblock = 1

  ! CUDA constants
  integer(c_int), parameter :: cudaMemcpyHostToDevice = 1
  integer(c_int), parameter :: cudaMemcpyDeviceToHost = 2

  integer(c_size_t) :: szreal_g, szint_g, szptr_g, offset = 0

  type(dim3) :: blocks, threads


!  type(cudaStream_t) :: stream1, stream2, stream3

  integer :: i, j, k, step
  real(8), volatile :: start, finish

  integer  ::  relax


!! Aplying the exchange !!TODO this part is CPU only and a subject for shared memory optimization

!   CALL exchange(dimension_id,0,phi) 



  szreal_g = sizeof(rel_g(1, 1, 1))
  szint_g = sizeof(nx_dev_g)
  szptr_g = sizeof(rel_dev2)
 
 

 ! We repeat GPU logic twice to show that the first kernel invocation includes
  ! the time of GPU device initialization, which is done one time per
  ! application run.

!    call cpu_time(start)



    ! Get the size of arrays passing to jacobi relaxation -- > ref : mod_laplace
    rel_size_g = shape(rel_g)

     


    cdg1_size_g = shape(cdg1)
    cdg2_size_g = shape(cdg2)
    cdg3_size_g = shape(cdg3)


    nx_dev_g = rel_size_g(1)  
    ny_dev_g = rel_size_g(2) 
    nz_dev_g = rel_size_g(3) 


    nx_cdg1_g = cdg1_size_g(1)
    ny_cdg1_g = cdg1_size_g(2)
    nz_cdg1_g = cdg1_size_g(3)

    nx_cdg2_g = cdg2_size_g(1)
    ny_cdg2_g = cdg2_size_g(2)
    nz_cdg2_g = cdg2_size_g(3)

    nx_cdg3_g = cdg3_size_g(1)
    ny_cdg3_g = cdg3_size_g(2)
    nz_cdg3_g = cdg3_size_g(3)


!    print*, nx, ny, nz, nx_cdg1, nx_cdg2, nx_cdg3, ny_cdg1, ny_cdg2, ny_cdg3, nz_cdg1, nz_cdg2, nz_cdg3, 'nx, ny, nz, nx_cdg1, nx_cdg2, nx_cdg3, ny_cdg1, ny_cdg2, ny_cdg3, nz_cdg1, nz_cdg2, nz_cdg3'
    ! Allocate memory on GPU.



    istat = cudaMalloc(cdg1_dev2, nx_cdg1_g * ny_cdg1_g * nz_cdg1_g * szreal_g)
    call cuda_err_check(istat)

    istat = cudaMalloc(cdg2_dev2, nx_cdg2_g * ny_cdg2_g * nz_cdg2_g * szreal_g)
    call cuda_err_check(istat)

    istat = cudaMalloc(cdg3_dev2, nx_cdg3_g * ny_cdg3_g * nz_cdg3_g * szreal_g)
    call cuda_err_check(istat)  
 


    istat = cudaMalloc(rel_dev2, nx_dev_g * ny_dev_g * nz_dev_g * szreal_g)
    call cuda_err_check(istat)
    ! rel, bb and comp come in same size -- ref to mod_laplace
   istat = cudaMalloc(comp_dev2, nx_dev_g * ny_dev_g * nz_dev_g *  szreal_g)
    call cuda_err_check(istat)

	! XXX Several blocks of code below are equivalent to sum_kernel<<<...>>> in C.

    ! Create CUDA compute grid configuration.
    
#ifdef CUDAVERSION9     
     blocks = dim3((nx_dev_g - 8)/8, (ny_dev_g - 8)/8, (nz_dev_g - 8)/8)
    if (mod(nx_dev_g, szblock) .ne. 0) blocks%x = blocks%x + 1
    threads = dim3((nx_dev_g - 8)/blocks%x, (ny_dev_g - 8)/blocks%y, (nz_dev_g - 8)/blocks%z)
#endif

!    blocks = dim3(1, 1, 1)
!    threads = dim3(1, 1, 1)
    ! Copy input data from CPU to GPU memory.
    istat = cudaMemcpy(cdg1_dev2, c_loc(cdg1), nx_cdg1_g * ny_cdg1_g * nz_cdg1_g * szreal_g, cudaMemcpyHostToDevice)
    call cuda_err_check(istat)
    istat = cudaMemcpy(cdg2_dev2, c_loc(cdg2), nx_cdg2_g * ny_cdg2_g * nz_cdg2_g * szreal_g, cudaMemcpyHostToDevice)
    call cuda_err_check(istat)
    istat = cudaMemcpy(cdg3_dev2, c_loc(cdg3), nx_cdg3_g * ny_cdg3_g * nz_cdg3_g * szreal_g, cudaMemcpyHostToDevice)
    call cuda_err_check(istat)
    

    istat = cudaMemcpy(rel_dev2, c_loc(rel_g), nx_dev_g * ny_dev_g * nz_dev_g * szreal_g, cudaMemcpyHostToDevice)
    call cuda_err_check(istat)

!!    istat = cudaMemcpy(comp_dev2, c_loc(comp), nx * ny * nz * szreal_g, cudaMemcpyHostToDevice)
!!    call cuda_err_check(istat)



! Creating CUDA Streams  !TEST -- wrapper is buggy -->  Maybe c_loc(stream1) can help
!    istat = cudaStreamCreate(stream1)
!    call cuda_err_check(istat)
!    istat = cudaStreamCreate(stream2)
!    call cuda_err_check(istat)
!    istat = cudaStreamCreate(stream3)
!    call cuda_err_check(istat)




!    istat = cudaMemcpy(grad_dev, c_loc(grad), nx * ny * nz * szreal_g, cudaMemcpyHostToDevice)
!    call cuda_err_check(istat)

#ifdef CUDAVERSION9
   ! Setup CUDA compute grid configuration.
    do relax =1, 1

        istat = cudaConfigureCall(blocks, threads, int8(0), c_null_ptr)
        call cuda_err_check(istat)
    
    !!! Stream 1 : the cuda API break down
    
        ! Setup CUDA kernel arguments (individually)
        offset = 0
        istat = cudaSetupScalarArgument(c_loc(grid_id), szint_g, offset)
        call cuda_err_check(istat)
        offset = offset + sizeof(grid_id)
        istat = cudaSetupScalarArgument(c_loc(nx_dev_g), szint_g, offset)
        call cuda_err_check(istat)
        offset = offset + sizeof(nx_dev_g)
        istat = cudaSetupScalarArgument(c_loc(ny_dev_g), szint_g, offset)
        call cuda_err_check(istat)
        offset = offset + sizeof(ny_dev_g)
        istat = cudaSetupScalarArgument(c_loc(nz_dev_g), szint_g, offset)
        call cuda_err_check(istat)
        offset = offset + sizeof(nz_dev_g)
    
        ! only the second index of the cdg matrices are not fixed and have to updated on the GPU
    
        istat = cudaSetupScalarArgument(c_loc(ny_cdg1_g), szint_g, offset)
        call cuda_err_check(istat)
        offset = offset + sizeof(ny_cdg1_g)
        istat = cudaSetupScalarArgument(c_loc(ny_cdg2_g), szint_g, offset)
        call cuda_err_check(istat)
        offset = offset + sizeof(ny_cdg2_g)
        istat = cudaSetupScalarArgument(c_loc(ny_cdg3_g), szint_g, offset)
        call cuda_err_check(istat)
        offset = offset + sizeof(ny_cdg3_g)
        istat = cudaSetupScalarArgument(c_loc(nx_cdg1_g), szint_g, offset)
        call cuda_err_check(istat)
        offset = offset + sizeof(nx_cdg1_g)
    
    
    
    
       ! Previous arguments were 4 bytes in size, total offset is 12 bytes.
        ! Since the next 3 arguments are 8 bytes each and have to be aligned
        ! by 8 bytes boundary, it's necessary to insert a 4-byte spacing here
        ! (from 12 to 16).
    !    offset = offset + 4
        istat = cudaSetupArrayArgument(cdg1_dev2, szptr_g, offset)
        call cuda_err_check(istat)
        offset = offset + sizeof(cdg1_dev2)
        istat = cudaSetupArrayArgument(cdg2_dev2, szptr_g, offset)
        call cuda_err_check(istat)
        offset = offset + sizeof(cdg2_dev2)
        istat = cudaSetupArrayArgument(cdg3_dev2, szptr_g, offset)
        call cuda_err_check(istat)
        offset = offset + sizeof(cdg3_dev2)
        istat = cudaSetupArrayArgument(rel_dev2, szptr_g, offset)
        call cuda_err_check(istat)
        offset = offset + sizeof(rel_dev2)
       istat = cudaSetupArrayArgument(comp_dev2, szptr_g, offset)
        call cuda_err_check(istat)
        offset = offset + sizeof(comp_dev2)
    
    
        ! Finally, launch CUDA kernel after we configured everything.
        ! Kernel is identified by C-string with its name (must be
        ! undecorated, i.e. with extern "C")
        istat = cudaLaunch('product_dg_relax_kernel' // c_null_char)
        call cuda_err_check(istat)

    enddo

#endif

  call product_dg_relax_launcher((grid_id), (nx_dev_g), (ny_dev_g),(nz_dev_g), (ny_cdg1_g), (ny_cdg2_g), ny_cdg3_g, nx_cdg1_g,cdg1_dev2,cdg2_dev2,cdg3_dev2, (rel_dev2), comp_dev2)


    ! Copy back results from GPU to CPU memory.
    istat = cudaMemcpy(c_loc(comp_g), comp_dev2, nx_dev_g * ny_dev_g * nz_dev_g * szreal_g, cudaMemcpyDeviceToHost)
    call cuda_err_check(istat)



    ! Free allocated GPU memory.
    istat = cudaFree(rel_dev2)
    call cuda_err_check(istat)
    istat = cudaFree(comp_dev2)
    call cuda_err_check(istat)
    istat = cudaFree(cdg1_dev2)
    call cuda_err_check(istat)
    istat = cudaFree(cdg2_dev2)
    call cuda_err_check(istat)
    istat = cudaFree(cdg3_dev2)
    call cuda_err_check(istat)


! Applying boundary conditions

end subroutine device_product_dg_relax

