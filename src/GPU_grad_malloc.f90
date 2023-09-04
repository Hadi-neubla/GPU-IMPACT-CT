!**************************************************************************************************                                           
!* by Hadi Zolfaghari, ARTORG Center, University of Bern, (zolfaghari.haadi@gmail.com)            *      
!* October 2015 - March 2020                                                                      *            
!**************************************************************************************************
!!!Subroutine for merged-up product div_grad
subroutine  grad_GPU_alloc(dimension_id)
USE mod_dims
USE mod_vars
USE mod_exchange
USE mod_diff
USE mod_laplace
USE mod_helmholtz
USE mod_inout
use GPU_vars    
use cudafor

integer, intent(in)    ::   dimension_id



    if (dimension_id .eq. 1)  cGp1_size = shape(cGp1)
    if (dimension_id .eq. 2)  cGp2_size = shape(cGp2)
    if (dimension_id .eq. 3)  cGp3_size = shape(cGp3)


    nx = phi_size(1)  
    ny = phi_size(2) 
    nz = phi_size(3) 
    if      (dimension_id .eq. 1) then
        nxg = cGp1_size(1)
        nyg = cGp1_size(2)
    else if (dimension_id .eq. 2) then
        nxg = cGp2_size(1)
        nyg = cGp2_size(2)
    else if (dimension_id .eq. 3) then
        nxg = cGp3_size(1)
        nyg = cGp3_size(2)
    endif


!    print*, nx, ny, nz, nxg, nyg, 'nx, ny, nz, nxg, nyg', 'this is in gradient'
    ! Allocate memory on GPU.

    if      (dimension_id .eq. 1) then
        istat = cudaMalloc(cGp1_dev, nxg * nyg * szreal)
        call cuda_err_check(istat)
    else if (dimension_id .eq. 2) then
        istat = cudaMalloc(cGp2_dev, nxg * nyg * szreal)
        call cuda_err_check(istat)
    else if (dimension_id .eq. 3) then
        istat = cudaMalloc(cGp3_dev, nxg * nyg * szreal)
        call cuda_err_check(istat)  
    end if


    istat = cudaMalloc(phi_dev, nx * ny * nz * szreal)
    call cuda_err_check(istat)
    ! phi and grad come in size size -- ref to mod_diff
    istat = cudaMalloc(grad_dev, nx * ny * nz *  szreal)
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
        istat = cudaMemcpy(cGp1_dev, c_loc(cGp1), nxg * nyg * szreal, cudaMemcpyHostToDevice)
        call cuda_err_check(istat)
    else if (dimension_id .eq. 2) then
        istat = cudaMemcpy(cGp2_dev, c_loc(cGp2), nxg * nyg * szreal, cudaMemcpyHostToDevice)
        call cuda_err_check(istat)
    else if (dimension_id .eq. 3) then
        istat = cudaMemcpy(cGp3_dev, c_loc(cGp3), nxg * nyg * szreal, cudaMemcpyHostToDevice)
        call cuda_err_check(istat)
    end if
    
    istat = cudaMemcpy(phi_dev, c_loc(phi), nx * ny * nz * szreal, cudaMemcpyHostToDevice)
    call cuda_err_check(istat)

end subroutine grad_GPU_alloc
