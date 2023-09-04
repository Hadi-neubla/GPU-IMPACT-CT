!**************************************************************************************************                                           
!* by Hadi Zolfaghari, ARTORG Center, University of Bern, (zolfaghari.haadi@gmail.com)            *      
!* October 2015 - March 2020                                                                      *            
!**************************************************************************************************

MODULE mod_vars_GPU

  USE mod_dims
  USE mod_vars 
  USE mod_cudafor 
  USE ISO_C_BINDING
  
  IMPLICIT NONE
  
  
  !===========================================================================================================
  !===== GPU acceleration mode ===============================================================================
  !===========================================================================================================
  integer                ::  GPU_accelerated 
  integer                ::  GPU_Poisson_Solver 
  integer                ::  GPU_verbose
  integer                ::  gpu_step

  !===========================================================================================================
  !===== pointer setup for memory operations  ================================================================
  !===========================================================================================================
  integer(c_int) :: istat

  integer, parameter :: szblock = 1

  ! CUDA constants
  integer(c_int), parameter :: cudaMemcpyHostToDevice = 1
  integer(c_int), parameter :: cudaMemcpyDeviceToHost = 2

  integer(c_size_t) :: szreal, szint, szptr, offset = 0

  !===========================================================================================================
  !===== tensors 
  !=========================================================================================================== 

   integer             ::   phi_device_size(3)
   integer, target     ::   nx_dev, ny_dev, nz_dev

  !===========================================================================================================
  !===== stencil coefficients
  !=========================================================================================================== 


  !Poisson solver
   integer             ::   cdg1_size(3)
   integer             ::   cdg2_size(3)
   integer             ::   cdg3_size(3)
   integer, target     ::   nx_cdg1, ny_cdg1, nz_cdg1 ! the grad and phi size components
   integer, target     ::   nx_cdg2, ny_cdg2, nz_cdg2 ! the grad and phi size components
   integer, target     ::   nx_cdg3, ny_cdg3, nz_cdg3 ! the grad and phi size component


  !gradient kernel
   integer, target     ::   nx_grad_coef1, ny_grad_coef1 !the cGp1 size components   
   integer, target     ::   nx_grad_coef2, ny_grad_coef2 !the cGp1 size components   
   integer, target     ::   nx_grad_coef3, ny_grad_coef3 !the cGp1 size components   
   integer             ::   cGp1_size(2)
   integer             ::   cGp2_size(2)
   integer             ::   cGp3_size(2)


   !divergence and divergence2 kernels
   integer             ::   cDu1_size(2)
   integer             ::   cDv2_size(2)
   integer             ::   cDw3_size(2)
   integer, target     ::   nx_div_coef1, ny_div_coef1 !the cDu1 size components
   integer, target     ::   nx_div_coef2, ny_div_coef2 !the cDu1 size components
   integer, target     ::   nx_div_coef3, ny_div_coef3 !the cDu1 size components

   ! interpolate_pre_vel kernel

   integer             ::    cIpu_size(2)
   integer             ::    cIpv_size(2)
   integer             ::    cIpw_size(2)
   integer, target     ::    nxIpu1, nyIpu1 !the cIup size components
   integer, target     ::    nxIpu2, nyIpu2 !the cIup size components
   integer, target     ::    nxIpu3, nyIpu3 !the cIup size components

   ! interpolate_vel_pre kernel
   integer             ::   cIup_size(2)
   integer             ::   cIvp_size(2)
   integer             ::   cIwp_size(2)


   integer, target     ::    nxIup1, nyIup1 !the cIup size components
   integer, target     ::    nxIup2, nyIup2 !the cIup size components
   integer, target     ::    nxIup3, nyIup3 !the cIup size components
 
   ! helmholtz explicit kernel

   integer            ::   cS1_size(2)
   integer            ::   cS2_size(2)
   integer            ::   cS3_size(2)
   integer            ::   cS4_size(2)
   integer            ::   cS5_size(2)
   integer            ::   cS6_size(2)



   integer, target     ::   nxs1, nys1
   integer, target     ::   nxs2, nys2
   integer, target     ::   nxs3, nys3 !the cS1 size components
   integer, target     ::   nxs4, nys4
   integer, target     ::   nxs5, nys5
   integer, target     ::   nxs6, nys6 !the cS1 size components

  !  nonlinear_upwind kernel
   integer            ::   cN1_size(2)
   integer            ::   cN2_size(2)
   integer            ::   cN3_size(2)
   integer            ::   cN4_size(2)
   integer            ::   cN5_size(2)
   integer            ::   cN6_size(2)




   integer, target     ::   nxn1, nyn1
   integer, target     ::   nxn2, nyn2
   integer, target     ::   nxn3, nyn3 !the cS1 size components
   integer, target     ::   nxn4, nyn4
   integer, target     ::   nxn5, nyn5
   integer, target     ::   nxn6, nyn6
 
  !===========================================================================================================
  !===== Device pointers =====================================================================================
  !===========================================================================================================

  ! Poisson solver
  type(c_ptr) :: rel_dev, bb_dev, comp_dev, cdg1_dev,  cdg2_dev, cdg3_dev
  type(c_ptr) :: lower_ghost_dev, upper_ghost_dev, front_ghost_dev, rear_ghost_dev, west_ghost_dev, east_ghost_dev

  ! gradient kernel
  type(c_ptr) :: grad_dev, phi_grad_dev, cGp1_dev, cGp2_dev, cGp3_dev  

  ! divergence kernel
  type(c_ptr) :: div_dev, phi_div_dev, cDu1_dev, cDv2_dev, cDw3_dev, div_dev_copy 

  ! divergence 2 kernel
  type(c_ptr) :: div2_dev, phi_div2_dev, cS1_dev, cS2_dev, cS3_dev 

  ! axpby kernel
  type(c_ptr) :: d_x, d_y, d_z

  ! helmholtz kernel
  type(c_ptr) :: helm_dev, vel1_dev, vel2_dev,vel3_dev, cS1_helm_dev, cS2_helm_dev, cS3_helm_dev, cS4_helm_dev, cS5_helm_dev, cS6_helm_dev
  
  ! interpolate vel_pre kernel (device_interpolate2)
  type(c_ptr) :: inter_vp_dev, phi_inter_vp_dev, cIup_dev, cIvp_dev, cIwp_dev
  
  ! interpolate pre_vel kernel (device_interpolate)
  type(c_ptr) :: inter_pv_dev, phi_inter_pv_dev, cIpu_dev, cIpv_dev, cIpw_dev

  ! nonlinear kernel
  type(c_ptr) :: nl_dev, pp_dev, vel_dev, cN1_dev, cN21_dev, cN2_dev, cN22_dev, cN3_dev, cN23_dev, cN4_dev, cN24_dev, cN5_dev, cN25_dev, cN6_dev, cN26_dev

  ! CT grid 
  type(c_ptr) :: ct_geometry_dev

  !===========================================================================================================
  !===== host (driver) variables =============================================================================
  !===========================================================================================================
   ! cuJOR pressure solver
   REAL(8)   , target, allocatable ::  west_ghost (:,:)
   REAL(8)   , target, allocatable ::  east_ghost (:,:)
   REAL(8)   , target, allocatable ::  front_ghost(:,:)
   REAL(8)   , target, allocatable ::  rear_ghost (:,:)
   REAL(8)   , target, allocatable ::  lower_ghost(:,:)
   REAL(8)   , target, allocatable ::  upper_ghost(:,:)


   REAL(8)   , target, allocatable ::  west_chunk (:,:)
   REAL(8)   , target, allocatable ::  east_chunk (:,:)
   REAL(8)   , target, allocatable ::  front_chunk(:,:)
   REAL(8)   , target, allocatable ::  rear_chunk (:,:)
   REAL(8)   , target, allocatable ::  lower_chunk(:,:)
   REAL(8)   , target, allocatable ::  upper_chunk(:,:)



 
  
END MODULE mod_vars_GPU
