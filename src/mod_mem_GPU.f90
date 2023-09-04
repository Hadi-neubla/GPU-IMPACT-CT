!**************************************************************************************************                                           
!* by Hadi Zolfaghari, ARTORG Center, University of Bern, (zolfaghari.haadi@gmail.com)            *      
!* October 2015 - March 2020                                                                      *            
!**************************************************************************************************

MODULE mod_mem_GPU

  USE mod_dims
  USE mod_vars 
  USE mod_vars_GPU
  USE mod_cudafor 
  USE ISO_C_BINDING
  
  IMPLICIT NONE

  contains 
!> subroutine that sets up the variables for drivers  
subroutine GPU_setup


 !general

    szreal = sizeof(pre(1, 1, 1))
    szint = sizeof(N1)
    szptr = sizeof(grad_dev)

    phi_device_size = shape(pre)

    nx_dev = phi_device_size(1)
    ny_dev = phi_device_size(2)
    nz_dev = phi_device_size(3)



! Poisson solver
    cdg1_size = shape(cdg1)
    cdg2_size = shape(cdg2)
    cdg3_size = shape(cdg3)


    nx_cdg1 = cdg1_size(1)
    ny_cdg1 = cdg1_size(2)
    nz_cdg1 = cdg1_size(3)

    nx_cdg2 = cdg2_size(1)
    ny_cdg2 = cdg2_size(2)
    nz_cdg2 = cdg2_size(3)

    nx_cdg3 = cdg3_size(1)
    ny_cdg3 = cdg3_size(2)
    nz_cdg3 = cdg3_size(3)

 !gradient kernel
    cGp1_size = shape(cGp1)
    cGp2_size = shape(cGp2)
    cGp3_size = shape(cGp3)


    nx_grad_coef1 = cGp1_size(1)
    ny_grad_coef1 = cGp1_size(2)
    nx_grad_coef2 = cGp2_size(1)
    ny_grad_coef2 = cGp2_size(2)
    nx_grad_coef3 = cGp3_size(1)
    ny_grad_coef3 = cGp3_size(2)

! divergence kernel

    cDu1_size = shape(cDu1)
    cDv2_size = shape(cDv2)
    cDw3_size = shape(cDw3)


    nx_div_coef1 = cDu1_size(1)
    ny_div_coef1 = cDu1_size(2)
    nx_div_coef2 = cDv2_size(1)
    ny_div_coef2 = cDv2_size(2)
    nx_div_coef3 = cDw3_size(1)
    ny_div_coef3 = cDw3_size(2)

! interpolate pre_vel kernel
    cIpu_size = shape(cIpu)
    cIpv_size = shape(cIpv)
    cIpw_size = shape(cIpw)

    nxIpu1 = cIpu_size(1)
    nyIpu1 = cIpu_size(2)
    nxIpu2 = cIpv_size(1)
    nyIpu2 = cIpv_size(2)
    nxIpu3 = cIpw_size(1)
    nyIpu3 = cIpw_size(2)


! interpolate_vel_pre_kernel
    cIup_size = shape(cIup)
    cIvp_size = shape(cIvp)
    cIwp_size = shape(cIwp)



    nxIup1 = cIup_size(1)
    nyIup1 = cIup_size(2)
    nxIup2 = cIvp_size(1)
    nyIup2 = cIvp_size(2)
    nxIup3 = cIwp_size(1)
    nyIup3 = cIwp_size(2)

! helmholtz explicit kernel

    cS1_size = shape(cu11)
    nxs1 = cS1_size(1)
    nys1 = cS1_size(2)
    cS2_size = shape(cp22)
    nxs2 = cS2_size(1)
    nys2 = cS2_size(2)
    cS3_size = shape(cp33)
    nxs3 = cS3_size(1)
    nys3 = cS3_size(2)
    cS4_size = shape(cv22)
    nxs4 = cS4_size(1)
    nys4 = cS4_size(2)
    cS5_size = shape(cp11)
    nxs5 = cS5_size(1)
    nys5 = cS5_size(2)
    cS6_size = shape(cw33)
    nxs6 = cS6_size(1)
    nys6 = cS6_size(2)

! nonlinear_upwind kernel

    cN1_size = shape(cu1)
    nxn1 = cN1_size(1)
    nyn1 = cN1_size(2)
    cN2_size = shape(cp2)
    nxn2 = cN2_size(1)
    nyn2 = cN2_size(2)
    cN3_size = shape(cp3)
    nxn3 = cN3_size(1)
    nyn3 = cN3_size(2)
    cN4_size = shape(cp1)
    nxn4 = cN4_size(1)
    nyn4 = cN4_size(2)
    cN5_size = shape(cv2)
    nxn5 = cN5_size(1)
    nyn5 = cN5_size(2)
    cN6_size = shape(cw3)
    nxn6 = cN6_size(1)
    nyn6 = cN6_size(2)

end subroutine GPU_setup

!!> subroutine that mallocs all the necessary memory on the GPU
subroutine GPU_malloc

! Poisson solver

    istat = cudaMalloc(cdg1_dev, nx_cdg1 * ny_cdg1 * nz_cdg1 * szreal)
    call cuda_err_check(istat)

    istat = cudaMalloc(cdg2_dev, nx_cdg2 * ny_cdg2 * nz_cdg2 * szreal)
    call cuda_err_check(istat)

    istat = cudaMalloc(cdg3_dev, nx_cdg3 * ny_cdg3 * nz_cdg3 * szreal)
    call cuda_err_check(istat)  
 
    istat = cudaMalloc(rel_dev, nx_dev * ny_dev * nz_dev * szreal)
    call cuda_err_check(istat)
    istat = cudaMalloc(bb_dev, nx_dev * ny_dev * nz_dev *  szreal)
    call cuda_err_check(istat)


!    istat = cudaMalloc(rel_dev, nx * ny * nz * szreal)
!    call cuda_err_check(istat)
!    istat = cudaMalloc(bb_dev, nx * ny * nz *  szreal)
!    call cuda_err_check(istat)

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




! gradient kernel
    istat = cudaMalloc(cGp1_dev, nx_grad_coef1 * ny_grad_coef1 * szreal)
    call cuda_err_check(istat)
    istat = cudaMalloc(cGp2_dev, nx_grad_coef2 * ny_grad_coef2 * szreal)
    call cuda_err_check(istat)
    istat = cudaMalloc(cGp3_dev, nx_grad_coef3 * ny_grad_coef3 * szreal)
    call cuda_err_check(istat)

    istat = cudaMalloc(phi_grad_dev, nx_dev * ny_dev * nz_dev * szreal)
    call cuda_err_check(istat)
    istat = cudaMalloc(grad_dev, nx_dev * ny_dev * nz_dev *  szreal)
    call cuda_err_check(istat)



! divergence kernel
    istat = cudaMalloc(cDu1_dev, nx_div_coef1 * ny_div_coef1 * szreal)
    call cuda_err_check(istat)
    istat = cudaMalloc(cDv2_dev, nx_div_coef2 * ny_div_coef2 * szreal)
    call cuda_err_check(istat)
    istat = cudaMalloc(cDw3_dev, nx_div_coef3 * ny_div_coef3 * szreal)
    call cuda_err_check(istat)

    istat = cudaMalloc(phi_div_dev, nx_dev * ny_dev * nz_dev * szreal)
    call cuda_err_check(istat)
    istat = cudaMalloc(div_dev, nx_dev * ny_dev * nz_dev *  szreal)
    call cuda_err_check(istat)


! divergence 2 kernel

    istat = cudaMalloc(div2_dev, nx_dev * ny_dev * nz_dev *  szreal)
    call cuda_err_check(istat)
    istat = cudaMalloc(phi_div2_dev, nx_dev * ny_dev * nz_dev * 3 * szreal)
    call cuda_err_check(istat)



! interpolate_pre_vel kernel

    istat = cudaMalloc(cIpu_dev, nxIpu1 * nyIpu1 * szreal)
    call cuda_err_check(istat)
    istat = cudaMalloc(cIpv_dev, nxIpu2 * nyIpu2 * szreal)
    call cuda_err_check(istat)
    istat = cudaMalloc(cIpw_dev, nxIpu3 * nyIpu3 * szreal)
    call cuda_err_check(istat)


    istat = cudaMalloc(phi_inter_pv_dev, nx_dev * ny_dev * nz_dev * szreal)
    call cuda_err_check(istat)
    istat = cudaMalloc(inter_pv_dev, nx_dev * ny_dev * nz_dev *  szreal)
    call cuda_err_check(istat)



! interpolate_vel_pre kernel

    istat = cudaMalloc(cIup_dev, nxIup1 * nyIup1 * szreal)
    call cuda_err_check(istat)
    istat = cudaMalloc(cIvp_dev, nxIup2 * nyIup2 * szreal)
    call cuda_err_check(istat)
    istat = cudaMalloc(cIwp_dev, nxIup3 * nyIup3 * szreal)
    call cuda_err_check(istat)


    istat = cudaMalloc(phi_inter_vp_dev, nx_dev * ny_dev * nz_dev * szreal)
    call cuda_err_check(istat)
    ! phi and inter come in size size -- ref to mod_diff
    istat = cudaMalloc(inter_vp_dev, nx_dev * ny_dev * nz_dev *  szreal)
    call cuda_err_check(istat)


!  axpby kernel

! helmholtz kernel

    istat = cudaMalloc(cS1_helm_dev, nxs1 * nys1 * szreal)
    call cuda_err_check(istat)
    istat = cudaMalloc(cS2_helm_dev, nxs2 * nys2 * szreal)
    call cuda_err_check(istat)
    istat = cudaMalloc(cS3_helm_dev, nxs3 * nys3 * szreal)
    call cuda_err_check(istat)
    istat = cudaMalloc(cS4_helm_dev, nxs4 * nys4 * szreal)
    call cuda_err_check(istat)
    istat = cudaMalloc(cS5_helm_dev, nxs5 * nys5 * szreal)
    call cuda_err_check(istat)
    istat = cudaMalloc(cS6_helm_dev, nxs6 * nys6 * szreal)
    call cuda_err_check(istat)

    istat = cudaMalloc(helm_dev, nx_dev * ny_dev * nz_dev *  szreal)
    call cuda_err_check(istat)

    istat = cudaMalloc(vel1_dev, nx_dev * ny_dev * nz_dev *  szreal)
    call cuda_err_check(istat)
    istat = cudaMalloc(vel2_dev, nx_dev * ny_dev * nz_dev *  szreal)
    call cuda_err_check(istat)
    istat = cudaMalloc(vel3_dev, nx_dev * ny_dev * nz_dev *  szreal)
    call cuda_err_check(istat)



! nonlinear_upwind kernel

    istat = cudamalloc(cN1_dev, nxn1 * nyn1 * szreal)
    call cuda_err_check(istat)
    istat = cudamalloc(cN21_dev, nxn1 * nyn1 * szreal)
    call cuda_err_check(istat)
    istat = cudaMalloc(cN2_dev, nxn2 * nyn2 * szreal)
    call cuda_err_check(istat)
    istat = cudaMalloc(cN22_dev, nxn2 * nyn2 * szreal)
    call cuda_err_check(istat)
    istat = cudaMalloc(cN3_dev, nxn3 * nyn3 * szreal)
    call cuda_err_check(istat)
    istat = cudaMalloc(cN23_dev, nxn3 * nyn3 * szreal)
    call cuda_err_check(istat)
    istat = cudaMalloc(cN4_dev, nxn4 * nyn4 * szreal)
    call cuda_err_check(istat)
    istat = cudaMalloc(cN24_dev, nxn4 * nyn4 * szreal)
    call cuda_err_check(istat)
    istat = cudaMalloc(cN5_dev, nxn5 * nyn5 * szreal)
    call cuda_err_check(istat)
    istat = cudaMalloc(cN25_dev, nxn5 * nyn5 * szreal)
    call cuda_err_check(istat)
    istat = cudaMalloc(cN6_dev, nxn6 * nyn6 * szreal)
    call cuda_err_check(istat)
    istat = cudaMalloc(cN26_dev, nxn6 * nyn6 * szreal)
    call cuda_err_check(istat)


    istat = cudaMalloc(pp_dev, nx_dev * ny_dev * nz_dev * szreal)
    call cuda_err_check(istat)
    istat = cudaMalloc(nl_dev, nx_dev * ny_dev * nz_dev *  szreal)
    call cuda_err_check(istat)
    istat = cudaMalloc(vel_dev, nx_dev * ny_dev * nz_dev *  szreal)
    call cuda_err_check(istat)


! axpby kernel

    istat = cudaMalloc(d_x, nx_dev * ny_dev * nz_dev * szreal)
    call cuda_err_check(istat)
    istat = cudaMalloc(d_y, nx_dev * ny_dev * nz_dev * szreal)
    call cuda_err_check(istat)
    istat = cudaMalloc(d_z, nx_dev * ny_dev * nz_dev * szreal)
    call cuda_err_check(istat)

! CT geometry
    istat = cudaMalloc(ct_geometry_dev, nx_dev * ny_dev * nz_dev * szreal)
    call cuda_err_check(istat)




end subroutine GPU_malloc
!
!!> subroutine that does the memcopies for all the variables which are not time-dependent
subroutine GPU_memcopy

! Poisson solver
        istat = cudaMemcpy(cdg1_dev, c_loc(cdg1), nx_cdg1 * ny_cdg1 * nz_cdg1 * szreal, cudaMemcpyHostToDevice)
        call cuda_err_check(istat)
        istat = cudaMemcpy(cdg2_dev, c_loc(cdg2), nx_cdg2 * ny_cdg2 * nz_cdg2 * szreal, cudaMemcpyHostToDevice)
        call cuda_err_check(istat)
        istat = cudaMemcpy(cdg3_dev, c_loc(cdg3), nx_cdg3 * ny_cdg3 * nz_cdg3 * szreal, cudaMemcpyHostToDevice)
        call cuda_err_check(istat)
! gradient kernel
        istat = cudaMemcpy(cGp1_dev, c_loc(cGp1), nx_grad_coef1 * ny_grad_coef1 * szreal, cudaMemcpyHostToDevice)
        call cuda_err_check(istat)
        istat = cudaMemcpy(cGp2_dev, c_loc(cGp2), nx_grad_coef2 * ny_grad_coef2 * szreal, cudaMemcpyHostToDevice)
        call cuda_err_check(istat)
        istat = cudaMemcpy(cGp3_dev, c_loc(cGp3), nx_grad_coef3 * ny_grad_coef3 * szreal, cudaMemcpyHostToDevice)
        call cuda_err_check(istat)
! divergence kernel
        istat = cudaMemcpy(cDu1_dev, c_loc(cDu1), nx_div_coef1 * ny_div_coef1 * szreal, cudaMemcpyHostToDevice)
        call cuda_err_check(istat)
        istat = cudaMemcpy(cDv2_dev, c_loc(cDv2), nx_div_coef2 * ny_div_coef2 * szreal, cudaMemcpyHostToDevice)
        call cuda_err_check(istat)
        istat = cudaMemcpy(cDw3_dev, c_loc(cDw3), nx_div_coef3 * ny_div_coef3 * szreal, cudaMemcpyHostToDevice)
        call cuda_err_check(istat)

!  interpolate_pre_vel kernel
        istat = cudaMemcpy(cIpu_dev, c_loc(cIpu), nxIpu1 * nyIpu1 * szreal, cudaMemcpyHostToDevice)
        call cuda_err_check(istat)
        istat = cudaMemcpy(cIpv_dev, c_loc(cIpv), nxIpu2 * nyIpu2 * szreal, cudaMemcpyHostToDevice)
        call cuda_err_check(istat)
        istat = cudaMemcpy(cIpw_dev, c_loc(cIpw), nxIpu3 * nyIpu3 * szreal, cudaMemcpyHostToDevice)
        call cuda_err_check(istat)
!  interpolate_vel_pre kernel

        istat = cudaMemcpy(cIup_dev, c_loc(cIup), nxIup1 * nyIup1 * szreal, cudaMemcpyHostToDevice)
        call cuda_err_check(istat)
        istat = cudaMemcpy(cIvp_dev, c_loc(cIvp), nxIup2 * nyIup2 * szreal, cudaMemcpyHostToDevice)
        call cuda_err_check(istat)
        istat = cudaMemcpy(cIwp_dev, c_loc(cIwp), nxIup3 * nyIup3 * szreal, cudaMemcpyHostToDevice)
        call cuda_err_check(istat)
!   helmholtz kernel

        istat = cudaMemcpy(cS1_helm_dev, c_loc(cu11), nxs1 * nys1 * szreal, cudaMemcpyHostToDevice)
        call cuda_err_check(istat)
        istat = cudaMemcpy(cS2_helm_dev, c_loc(cp22), nxs2 * nys2 * szreal, cudaMemcpyHostToDevice)
        call cuda_err_check(istat)
        istat = cudaMemcpy(cS3_helm_dev, c_loc(cp33), nxs3 * nys3 * szreal, cudaMemcpyHostToDevice)
        call cuda_err_check(istat)
        istat = cudaMemcpy(cS4_helm_dev, c_loc(cv22), nxs4 * nys4 * szreal, cudaMemcpyHostToDevice)
        call cuda_err_check(istat)
        istat = cudaMemcpy(cS5_helm_dev, c_loc(cp11), nxs5 * nys5 * szreal, cudaMemcpyHostToDevice)
        call cuda_err_check(istat)
        istat = cudaMemcpy(cS6_helm_dev, c_loc(cw33), nxs6 * nys6 * szreal, cudaMemcpyHostToDevice)
        call cuda_err_check(istat)

! nonlinear_upwind kernel

        istat = cudaMemcpy(cN1_dev, c_loc(cNu1U), nxn1 * nyn1 * szreal, cudaMemcpyHostToDevice)
        call cuda_err_check(istat)
        istat = cudaMemcpy(cN21_dev, c_loc(cNu1D), nxn1 * nyn1 * szreal, cudaMemcpyHostToDevice)
        call cuda_err_check(istat)
        istat = cudaMemcpy(cN2_dev, c_loc(cNp2U), nxn2 * nyn2 * szreal, cudaMemcpyHostToDevice)
        call cuda_err_check(istat)
        istat = cudaMemcpy(cN22_dev, c_loc(cNp2D), nxn2 * nyn2 * szreal, cudaMemcpyHostToDevice)
        call cuda_err_check(istat)
        istat = cudaMemcpy(cN3_dev, c_loc(cNp3U), nxn3 * nyn3 * szreal, cudaMemcpyHostToDevice)
        call cuda_err_check(istat)
        istat = cudaMemcpy(cN23_dev, c_loc(cNp3D), nxn3 * nyn3 * szreal, cudaMemcpyHostToDevice)
        call cuda_err_check(istat)
        istat = cudaMemcpy(cN4_dev, c_loc(cNp1U), nxn4 * nyn4 * szreal, cudaMemcpyHostToDevice)
        call cuda_err_check(istat)
        istat = cudaMemcpy(cN24_dev, c_loc(cNp1D), nxn4 * nyn4 * szreal, cudaMemcpyHostToDevice)
        call cuda_err_check(istat)
        istat = cudaMemcpy(cN5_dev, c_loc(cNv2U), nxn5 * nyn5 * szreal, cudaMemcpyHostToDevice)
        call cuda_err_check(istat)
        istat = cudaMemcpy(cN25_dev, c_loc(cNv2D), nxn5 * nyn5 * szreal, cudaMemcpyHostToDevice)
        call cuda_err_check(istat)
        istat = cudaMemcpy(cN6_dev, c_loc(cNw3U), nxn6 * nyn6 * szreal, cudaMemcpyHostToDevice)
        call cuda_err_check(istat)
        istat = cudaMemcpy(cN26_dev, c_loc(cNw3D), nxn6 * nyn6 * szreal, cudaMemcpyHostToDevice)
        call cuda_err_check(istat)

 ! CT geometry 

        istat = cudaMemcpy(ct_geometry_dev, c_loc(ct_geometry), nx_dev * ny_dev * nz_dev * szreal, cudaMemcpyHostToDevice)
        call cuda_err_check(istat)

end subroutine GPU_memcopy
!
!!> subroutine that frees all the GPU memory at the end of the run 
subroutine GPU_memfree
! Poisson solver
        istat = cudaFree(cdg1_dev)
        call cuda_err_check(istat)
        istat = cudaFree(cdg2_dev)
        call cuda_err_check(istat)
        istat = cudaFree(cdg3_dev)
        call cuda_err_check(istat)

        istat = cudaFree(rel_dev)
        call cuda_err_check(istat)
        istat = cudaFree(bb_dev)
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


! gradient kernel
        istat = cudaFree(cGp1_dev)
        call cuda_err_check(istat)
        istat = cudaFree(cGp2_dev)
        call cuda_err_check(istat)
        istat = cudaFree(cGp3_dev)
        call cuda_err_check(istat)

        istat = cudaFree(grad_dev)
        call cuda_err_check(istat)
        istat = cudaFree(phi_grad_dev)
        call cuda_err_check(istat)




! divergence kernel
        istat = cudaFree(cDu1_dev)
        call cuda_err_check(istat)
        istat = cudaFree(cDv2_dev)
        call cuda_err_check(istat)
        istat = cudaFree(cDw3_dev)
        call cuda_err_check(istat)

        istat = cudaFree(div_dev)
        call cuda_err_check(istat)
        istat = cudaFree(phi_div_dev)
        call cuda_err_check(istat)

! divergence2 kernel

        istat = cudaFree(div2_dev)
        call cuda_err_check(istat)
        istat = cudaFree(phi_div2_dev)
        call cuda_err_check(istat)


! interpolate_pre_vel kernel

        istat = cudaFree(cIpu_dev)
        call cuda_err_check(istat)
        istat = cudaFree(cIpv_dev)
        call cuda_err_check(istat)
        istat = cudaFree(cIpw_dev)
        call cuda_err_check(istat)

       istat = cudaFree(inter_pv_dev)
       call cuda_err_check(istat)
       istat = cudaFree(phi_inter_pv_dev)
       call cuda_err_check(istat)



! interpolate_vel_pre kernel

        istat = cudaFree(cIup_dev)
        call cuda_err_check(istat)
        istat = cudaFree(cIvp_dev)
        call cuda_err_check(istat)
        istat = cudaFree(cIwp_dev)
        call cuda_err_check(istat)


        istat = cudaFree(inter_vp_dev)
        call cuda_err_check(istat)
        istat = cudaFree(phi_inter_vp_dev)
        call cuda_err_check(istat)


! helmholtz kernel

        istat = cudaFree(cS1_helm_dev)
        call cuda_err_check(istat)
        istat = cudaFree(cS2_helm_dev)
        call cuda_err_check(istat)
        istat = cudaFree(cS3_helm_dev)
        call cuda_err_check(istat)
        istat = cudaFree(cS4_helm_dev)
        call cuda_err_check(istat)
        istat = cudaFree(cS5_helm_dev)
        call cuda_err_check(istat)
        istat = cudaFree(cS6_helm_dev)
        call cuda_err_check(istat)

        istat = cudaFree(helm_dev)
        call cuda_err_check(istat)
        istat = cudaFree(vel1_dev)
        call cuda_err_check(istat)
        istat = cudaFree(vel2_dev)
        call cuda_err_check(istat)
        istat = cudaFree(vel3_dev)
        call cuda_err_check(istat)


! nonlinear upwind kernel

        istat = cudaFree(cN1_dev)
        call cuda_err_check(istat)
        istat = cudaFree(cN21_dev)
        call cuda_err_check(istat)
        istat = cudaFree(cN2_dev)
        call cuda_err_check(istat)
        istat = cudaFree(cN22_dev)
        call cuda_err_check(istat)
        istat = cudaFree(cN3_dev)
        call cuda_err_check(istat)
        istat = cudaFree(cN23_dev)
        call cuda_err_check(istat)
        istat = cudaFree(cN4_dev)
        call cuda_err_check(istat)
        istat = cudaFree(cN24_dev)
        call cuda_err_check(istat)
        istat = cudaFree(cN5_dev)
        call cuda_err_check(istat)
        istat = cudaFree(cN25_dev)
        call cuda_err_check(istat)
        istat = cudaFree(cN6_dev)
        call cuda_err_check(istat)
        istat = cudaFree(cN26_dev)
        call cuda_err_check(istat)


        istat = cudaFree(nl_dev)
        call cuda_err_check(istat)
        istat = cudaFree(vel_dev)
        call cuda_err_check(istat)
        istat = cudaFree(pp_dev)
        call cuda_err_check(istat)

! axpby kernel

       istat = cudaFree(d_x)
       call cuda_err_check(istat)
       istat = cudaFree(d_y)
       call cuda_err_check(istat)
       istat = cudaFree(d_z)
       call cuda_err_check(istat)

! CT geometry
       istat = cudaFree(ct_geometry_dev)
       call cuda_err_check(istat)


end subroutine GPU_memfree
  
END MODULE mod_mem_GPU
