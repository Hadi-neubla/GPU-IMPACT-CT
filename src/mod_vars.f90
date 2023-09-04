!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!* modified by Barna Becsek, ARTORG CVE, University of Bern (barna.becsek@artorg.unibe.ch)                   *
!* October 2014                                                                                              *
!* GPU version by Hadi Zolfaghari , ARTORG CVE, then DAMTP, Cambridge University (hz382@cam.ac.uk)           *
!* Oct 2015 - Sep 2023                                                                                       *
!*************************************************************************************************************


!*************************************************************************************************************
!* Note: (bbecsek 110216)                                                                                    *
!*   New variables of PROTECTED type had to be added in order to make them accessible from C.                *
!*   PARAMETER type values cannot be accessed from C.                                                        *
!*   Other variables are not assigned ISO_C_BINDINGS because we are only using GCC compilers, and thus       *
!*   we do not need compiler independent access.                                                             *
!*   Arrays are accessed through pointers from the binding.                                                  *
!*************************************************************************************************************

!> @file mod_vars.f90
!! File containing variable declarations

!> All the variables are declared here. It uses the module mod_vars
MODULE mod_vars
  
  USE mod_dims
  USE ISO_C_BINDING !bbecsek
  
  IMPLICIT NONE
  
  
#ifdef ALLOC
  
  !===========================================================================================================
  !===== GPU acceleration mode ===============================================================================
  !===========================================================================================================
  ! Hadi Zolfaghari
!  integer                ::  GPU_accelerated 
!  integer                ::  GPU_verbose
!  integer                ::  gpu_step
  !===========================================================================================================
  !===== Device pointers =====================================================================================
  !===========================================================================================================

!  type(c_ptr) :: rel_dev, bb_dev, comp_dev, cdg1_dev,  cdg2_dev, cdg3_dev
!  type(c_ptr) :: lower_ghost_dev, upper_ghost_dev, front_ghost_dev, rear_ghost_dev, west_ghost_dev, east_ghost_dev


  !===========================================================================================================
  !===== host variables =====================================================================================
  !===========================================================================================================

!   REAL(8)   , target, allocatable ::  west_ghost (:,:)
!   REAL(8)   , target, allocatable ::  east_ghost (:,:)
!   REAL(8)   , target, allocatable ::  front_ghost(:,:)
!   REAL(8)   , target, allocatable ::  rear_ghost (:,:)
!   REAL(8)   , target, allocatable ::  lower_ghost(:,:)
!   REAL(8)   , target, allocatable ::  upper_ghost(:,:)
!
!
!   REAL(8)   , target, allocatable ::  west_chunk (:,:)
!   REAL(8)   , target, allocatable ::  east_chunk (:,:)
!   REAL(8)   , target, allocatable ::  front_chunk(:,:)
!   REAL(8)   , target, allocatable ::  rear_chunk (:,:)
!   REAL(8)   , target, allocatable ::  lower_chunk(:,:)
!   REAL(8)   , target, allocatable ::  upper_chunk(:,:)
!
!   REAL(8)   , target,  ALLOCATABLE   ::  cdg1 (:,:,:)
!   REAL(8)   , target,  ALLOCATABLE   ::  cdg2 (:,:,:)
!   REAL(8)   , target,  ALLOCATABLE   ::  cdg3 (:,:,:)
 


  !===========================================================================================================
  !=== Raemliche Dimensionen =================================================================================
  !===========================================================================================================
  ! 2D wird ueber M3==2 eingeschaltet
  INTEGER                ::  dimens
  
  
  !===========================================================================================================
  !=== Domain- und Blockspezifikationen ======================================================================
  !===========================================================================================================
  !--- zulaessige Blockgroessen ------------------------------------------------------------------------------
  ! n  x 1   2   4   8   16   32   64   128   256  512   ...
  !--------------------------------------------------------------
  ! n2 = 2,  4,  8, 16,  32,  64, 128,  256,  512, 1024, ...
  ! n3 = 3,  6, 12, 24,  48,  96, 192,  384,  768, 1536, ...
  ! n5 = 5, 10, 20, 40,  80, 160, 320,  640, 1280, 2560, ...
  ! n7 = 7, 14, 28, 56, 112, 224, 448,  896, 1792, 3584, ...
  ! n9 = 9, 18, 36, 72, 144, 288, 576, 1152, 2304, 4608, ...
  ! ...
  INTEGER                ::  N1 
  INTEGER                ::  N2
  INTEGER                ::  N3
  
  !--- Anzahl grobe Gitter (Multigrid) -----------------------------------------------------------------------
  INTEGER, PARAMETER     ::  n_grids_max = 15
  INTEGER                ::  n_grids, n_grids_limit
  !--- C-FORTRAN interface ---
  INTEGER, PROTECTED, bind(C, name='_n_grids_max_c') :: n_grids_max_c = n_grids_max

  REAL(8)                   ::  rho_fluid               !< bulk density of blood

  
  !--- Dimensionen -------------------------------------------------------------------------------------------
  INTEGER, PARAMETER     ::  dim_ncb1c = SIZE(ncb1c)
  INTEGER, PARAMETER     ::  dim_ncb1g = SIZE(ncb1g)
  INTEGER, PARAMETER     ::  dim_ncb1d = SIZE(ncb1d)
  
  INTEGER, PARAMETER     ::  dim_ncb2c = SIZE(ncb2c)
  INTEGER, PARAMETER     ::  dim_ncb2g = SIZE(ncb2g)
  INTEGER, PARAMETER     ::  dim_ncb2d = SIZE(ncb2d)
  
  INTEGER, PARAMETER     ::  dim_ncb3c = SIZE(ncb3c)
  INTEGER, PARAMETER     ::  dim_ncb3g = SIZE(ncb3g)
  INTEGER, PARAMETER     ::  dim_ncb3d = SIZE(ncb3d)
  !--- C-Fortran interface ---
  INTEGER, PROTECTED, bind(C, name='_dim_ncb1c_c') :: dim_ncb1c_c = dim_ncb1c
  INTEGER, PROTECTED, bind(C, name='_dim_ncb1g_c') :: dim_ncb1g_c = dim_ncb1g
  INTEGER, PROTECTED, bind(C, name='_dim_ncb1d_c') :: dim_ncb1d_c = dim_ncb1d
      
  INTEGER, PROTECTED, bind(C, name='_dim_ncb2c_c') :: dim_ncb2c_c = dim_ncb2c
  INTEGER, PROTECTED, bind(C, name='_dim_ncb2g_c') :: dim_ncb2g_c = dim_ncb2g
  INTEGER, PROTECTED, bind(C, name='_dim_ncb2d_c') :: dim_ncb2d_c = dim_ncb2d
         
  INTEGER, PROTECTED, bind(C, name='_dim_ncb3c_c') :: dim_ncb3c_c = dim_ncb3c
  INTEGER, PROTECTED, bind(C, name='_dim_ncb3g_c') :: dim_ncb3g_c = dim_ncb3d
  INTEGER, PROTECTED, bind(C, name='_dim_ncb3d_c') :: dim_ncb3d_c = dim_ncb3g
  

  !--- Anzahl Stencil-Koeffizienten (Feld) -------------------------------------------------------------------
  ! Anzahl der Koeffizienten im Feld (zentrale Differenzen angenommen):
  INTEGER, PARAMETER     ::  nc1c = ncb1c(dim_ncb1c)
  INTEGER, PARAMETER     ::  nc1s = ncb1g(dim_ncb1g)
  
  INTEGER, PARAMETER     ::  nc2c = ncb2c(dim_ncb2c)
  INTEGER, PARAMETER     ::  nc2s = ncb2g(dim_ncb2g)
  
  INTEGER, PARAMETER     ::  nc3c = ncb3c(dim_ncb3c)
  INTEGER, PARAMETER     ::  nc3s = ncb3g(dim_ncb3g)
  !--- C-Fortran interface ---
  INTEGER, PROTECTED, bind(C, name='_nc1c_c') :: nc1c_c = nc1c
  INTEGER, PROTECTED, bind(C, name='_nc1s_c') :: nc1s_c = nc1s
  
  INTEGER, PROTECTED, bind(C, name='_nc2c_c') :: nc2c_c = nc2c
  INTEGER, PROTECTED, bind(C, name='_nc2s_c') :: nc2s_c = nc2s
  
  INTEGER, PROTECTED, bind(C, name='_nc3c_c') :: nc3c_c = nc3c
  INTEGER, PROTECTED, bind(C, name='_nc3s_c') :: nc3s_c = nc3s
  
  !===========================================================================================================
  !=== Intervallgrenzen der Differenzen-Koeffizienten-Arrays =================================================
  !===========================================================================================================
  !--- zentral -----------------------------------------------------------------------------------------------
  INTEGER, PARAMETER     ::  b1U = nc1s/2
  INTEGER, PARAMETER     ::  b2U = nc2s/2
  INTEGER, PARAMETER     ::  b3U = nc3s/2
  
  INTEGER, PARAMETER     ::  b1L = -b1U
  INTEGER, PARAMETER     ::  b2L = -b2U
  INTEGER, PARAMETER     ::  b3L = -b3U
  !--- C-FORTRAN interface ---
  INTEGER, PROTECTED, bind(C, name='_b1u_c') :: b1u_c = b1U
  INTEGER, PROTECTED, bind(C, name='_b2u_c') :: b2u_c = b2U
  INTEGER, PROTECTED, bind(C, name='_b3u_c') :: b3u_c = b3U

  INTEGER, PROTECTED, bind(C, name='_b1l_c') :: b1l_c = b1L
  INTEGER, PROTECTED, bind(C, name='_b2l_c') :: b2l_c = b2L
  INTEGER, PROTECTED, bind(C, name='_b3l_c') :: b3l_c = b3L
  
  !--- upwind (nicht-linear) ---------------------------------------------------------------------------------
  ! (aktuell wird nicht zwischen auf- und abwärtsgerichteten Stencils unterschieden, um auch am Rand arbeiten
  !  zu können, wo KEINE upwind-Differenzierung verwendet wird)
  INTEGER, PARAMETER     ::  n1L = b1L
  INTEGER, PARAMETER     ::  n2L = b2L
  INTEGER, PARAMETER     ::  n3L = b3L
  
  INTEGER, PARAMETER     ::  n1U = b1U
  INTEGER, PARAMETER     ::  n2U = b2U
  INTEGER, PARAMETER     ::  n3U = b3U
  !--- C-FORTRAN interface ---
  INTEGER, PROTECTED, bind(C, name='_n1u_c') :: n1u_c = n1U
  INTEGER, PROTECTED, bind(C, name='_n2u_c') :: n2u_c = n2U
  INTEGER, PROTECTED, bind(C, name='_n3u_c') :: n3u_c = n3U

  INTEGER, PROTECTED, bind(C, name='_n1l_c') :: n1l_c = n1L
  INTEGER, PROTECTED, bind(C, name='_n2l_c') :: n2l_c = n2L
  INTEGER, PROTECTED, bind(C, name='_n3l_c') :: n3l_c = n3L
  
  !--- Divergenz ---------------------------------------------------------------------------------------------
  INTEGER, PARAMETER     ::  d1L = b1L
  INTEGER, PARAMETER     ::  d2L = b2L
  INTEGER, PARAMETER     ::  d3L = b3L
  
  INTEGER, PARAMETER     ::  d1U = b1U-1
  INTEGER, PARAMETER     ::  d2U = b2U-1
  INTEGER, PARAMETER     ::  d3U = b3U-1
  !--- C-FORTRAN interface ---
  INTEGER, PROTECTED, bind(C, name='_d1u_c') :: d1u_c = d1U
  INTEGER, PROTECTED, bind(C, name='_d2u_c') :: d2u_c = d2U
  INTEGER, PROTECTED, bind(C, name='_d3u_c') :: d3u_c = d3U

  INTEGER, PROTECTED, bind(C, name='_d1l_c') :: d1l_c = d1L
  INTEGER, PROTECTED, bind(C, name='_d2l_c') :: d2l_c = d2L
  INTEGER, PROTECTED, bind(C, name='_d3l_c') :: d3l_c = d3L
  
  !--- Gradient ----------------------------------------------------------------------------------------------
  INTEGER, PARAMETER     ::  g1L = b1L+1
  INTEGER, PARAMETER     ::  g2L = b2L+1
  INTEGER, PARAMETER     ::  g3L = b3L+1
  
  INTEGER, PARAMETER     ::  g1U = b1U
  INTEGER, PARAMETER     ::  g2U = b2U
  INTEGER, PARAMETER     ::  g3U = b3U
  !--- C-FORTRAN interface ---
  INTEGER, PROTECTED, bind(C, name='_g1u_c') :: g1u_c = g1U
  INTEGER, PROTECTED, bind(C, name='_g2u_c') :: g2u_c = g2U
  INTEGER, PROTECTED, bind(C, name='_g3u_c') :: g3u_c = g3U

  INTEGER, PROTECTED, bind(C, name='_g1l_c') :: g1l_c = g1L
  INTEGER, PROTECTED, bind(C, name='_g2l_c') :: g2l_c = g2L
  INTEGER, PROTECTED, bind(C, name='_g3l_c') :: g3l_c = g3L
  

  !--- CT-based aortic model 
       !--- geometry
  INTEGER                                    :: i_start_aorta_inlet, j_start_aorta_inlet, i_end_aorta_inlet, j_end_aorta_inlet

  REAL(8)                                       :: aortic_diameter, maximum_inflow_velocity
  REAL(8)                                       :: aortic_diameter_star, maximum_inflow_velocity_star
  REAL(8)                                       :: CT_pixel_width_star, CT_pixel_width, upsampling_factor
  REAL(8)                                       :: reference_length, reference_velocity, reference_time
  REAL(8), PARAMETER                            :: blood_kinematic_viscosity = 0.0000035
  REAL(8), PARAMETER                            :: blood_density             = 1060.

       !---pulse
  ! --- TAVI000:75, TAVI001:36, TAVI002:36, TAVI003:91
  INTEGER, PARAMETER  :: number_of_digitized_points=74
  INTEGER             :: pulse_number
  REAL(8)                :: pulse_length           ! dimensionless 
  REAL(8)                :: shifted_time_in_pulse  ! this is: time - (pulse_number-1) * pulse * length 

  TYPE :: time_Q_pair
     REAL(8) :: time_in_pulse
     REAL(8) :: Q_in_pulse
  END TYPE time_Q_pair

  !--- the pair for reading original CSV file contating flowrate-time data
  TYPE(time_Q_pair), DIMENSION(number_of_digitized_points) :: waveform_data
  !--- the pair for dimensionless velocity-time-data
  TYPE(time_Q_pair), DIMENSION(number_of_digitized_points) :: waveform_data_converted

  !--- in-flow constriction parameter 
  REAL(8)                 ::  flow_rate_correction_factor

  !--- parameters for the idealized model
  INTEGER              ::   inlet_voxels_count, outlet_voxels_count
  INTEGER              ::   idealized_aorta = 1, CT_based_aorta = 0 
   
  !--- parameters for buffer zones at inlet and outlet
  INTEGER                        ::   k_start_index, counter_inlet, counter_outlet
  INTEGER, ALLOCATABLE           ::   x_inlet_index(:), y_inlet_index(:)
  INTEGER, ALLOCATABLE           ::   x_outlet_index(:), y_outlet_index(:)
  ! --- cranial outlets (outlet2)
  INTEGER                        ::   k_start_index_2, counter_outlet_2
  INTEGER, ALLOCATABLE           ::   x_outlet_2_index(:), y_outlet_2_index(:)


  !--- parameters for the two-element windkessel model
  REAL(8)                           ::   p_end_systole_in, p_end_systole_out, p_end_systole_out_2 
  REAL(8)                           ::   pre_inlet, q_inlet, pre_outlet, q_outlet, pre_outlet_2, q_outlet_2
  LOGICAL, parameter             ::   aorta_rank_set = .FALSE. 
  !===========================================================================================================
  !=== Differenzen-Koeffizienten-Arrays ======================================================================
  !===========================================================================================================
  !--- 1. Ableitung (zentral) --------------------------------------------------------------------------------
  !--- hz : target attributes are specified for GPU c_loc memory operations
  REAL(8)   , target,  ALLOCATABLE   ::  cp1  (:,:)
  REAL(8)   , target,  ALLOCATABLE   ::  cp2  (:,:)
  REAL(8)   , target,  ALLOCATABLE   ::  cp3  (:,:)
  
  REAL(8)   , target, ALLOCATABLE   ::  cu1  (:,:)
  REAL(8)   , target, ALLOCATABLE   ::  cv2  (:,:)
  REAL(8)   , target, ALLOCATABLE   ::  cw3  (:,:)
  
  !--- 1. Ableitung (upwind) ---------------------------------------------------------------------------------
  REAL(8)   , target, ALLOCATABLE   ::  cNp1D(:,:)
  REAL(8)   , target, ALLOCATABLE   ::  cNp2D(:,:)
  REAL(8)   , target, ALLOCATABLE   ::  cNp3D(:,:)
  
  REAL(8)   , target, ALLOCATABLE   ::  cNp1U(:,:)
  REAL(8)   , target, ALLOCATABLE   ::  cNp2U(:,:)
  REAL(8)   , target, ALLOCATABLE   ::  cNp3U(:,:)
  
  REAL(8)   , target, ALLOCATABLE   ::  cNu1D(:,:)
  REAL(8)   , target, ALLOCATABLE   ::  cNv2D(:,:)
  REAL(8)   , target, ALLOCATABLE   ::  cNw3D(:,:)
  
  REAL(8)   , target, ALLOCATABLE   ::  cNu1U(:,:)
  REAL(8)   , target, ALLOCATABLE   ::  cNv2U(:,:)
  REAL(8)   , target, ALLOCATABLE   ::  cNw3U(:,:)
  
  !--- Divergenz ---------------------------------------------------------------------------------------------
  REAL(8)   , target, ALLOCATABLE   ::  cDu1 (:,:)
  REAL(8)   , target, ALLOCATABLE   ::  cDv2 (:,:)
  REAL(8)   , target, ALLOCATABLE   ::  cDw3 (:,:)
  
  !--- Divergenz (transponiert) ------------------------------------------------------------------------------
  REAL(8)   , ALLOCATABLE   ::  cDu1T(:,:)
  REAL(8)   , ALLOCATABLE   ::  cDv2T(:,:)
  REAL(8)   , ALLOCATABLE   ::  cDw3T(:,:)
  
  !--- Gradient ----------------------------------------------------------------------------------------------
  REAL(8)   ,target,  ALLOCATABLE   ::  cGp1 (:,:)
  REAL(8)   ,target,  ALLOCATABLE   ::  cGp2 (:,:)
  REAL(8)   ,target,  ALLOCATABLE   ::  cGp3 (:,:)
  
  !--- Gradient (transponiert) -------------------------------------------------------------------------------
  REAL(8)   , ALLOCATABLE   ::  cGp1T(:,:)
  REAL(8)   , ALLOCATABLE   ::  cGp2T(:,:)
  REAL(8)   , ALLOCATABLE   ::  cGp3T(:,:)
  
  !--- 2. Ableitung (zentral) --------------------------------------------------------------------------------
  REAL(8)   , target, ALLOCATABLE   ::  cp11 (:,:)
  REAL(8)   , target, ALLOCATABLE   ::  cp22 (:,:)
  REAL(8)   , target, ALLOCATABLE   ::  cp33 (:,:)
  
  REAL(8)   , target, ALLOCATABLE   ::  cu11 (:,:)
  REAL(8)   , target, ALLOCATABLE   ::  cv22 (:,:)
  REAL(8)   , target, ALLOCATABLE   ::  cw33 (:,:)
  
  REAL(8)   , target, ALLOCATABLE   ::  cS1 (:,:)
  !--- Interpolation ----------------------------------------------------------------------------------------- 
  REAL(8)   , target, ALLOCATABLE   ::  cIpu(:,:)
  REAL(8)   , target, ALLOCATABLE   ::  cIpv(:,:)
  REAL(8)   , target, ALLOCATABLE   ::  cIpw(:,:)
  
  REAL(8)   , target, ALLOCATABLE   ::  cIup(:,:)
  REAL(8)   , target, ALLOCATABLE   ::  cIvp(:,:)
  REAL(8)   , target, ALLOCATABLE   ::  cIwp(:,:)
  
  !--- Filter ------------------------------------------------------------------------------------------------
  REAL(8)   , ALLOCATABLE   ::  cFp1(:,:)
  REAL(8)   , ALLOCATABLE   ::  cFp2(:,:)
  REAL(8)   , ALLOCATABLE   ::  cFp3(:,:)
  
  REAL(8)   , ALLOCATABLE   ::  cFu1(:,:)
  REAL(8)   , ALLOCATABLE   ::  cFv2(:,:)
  REAL(8)   , ALLOCATABLE   ::  cFw3(:,:)
  
  !--- Integrator (nur für Druckgitter) ----------------------------------------------------------------------
  REAL(8)   , ALLOCATABLE   ::  cInt1(:,:)
  REAL(8)   , ALLOCATABLE   ::  cInt2(:,:)
  REAL(8)   , ALLOCATABLE   ::  cInt3(:,:)  
  
  !--- 2. Ableitung (Multigrid) ------------------------------------------------------------------------------ 
  ! Anmerkung: Die Koeffizientensätze unterscheiden sich z.T. lediglich durch die Randbedingungen.
  REAL(8)   , ALLOCATABLE   ::  cp11R(:,:,:)
  REAL(8)   , ALLOCATABLE   ::  cp22R(:,:,:)
  REAL(8)   , ALLOCATABLE   ::  cp33R(:,:,:)
  
  REAL(8)   , ALLOCATABLE   ::  cu11R(:,:,:)
  REAL(8)   , ALLOCATABLE   ::  cv22R(:,:,:)
  REAL(8)   , ALLOCATABLE   ::  cw33R(:,:,:)
  
  REAL(8)   , target,  ALLOCATABLE   ::  cdg1 (:,:,:)
  REAL(8)   , target,  ALLOCATABLE   ::  cdg2 (:,:,:)
  REAL(8)   , target,  ALLOCATABLE   ::  cdg3 (:,:,:)
  
  !--- Interpolation (Multigrid) ----------------------------------------------------------------------------- 
  REAL(8)   , ALLOCATABLE   ::  cI1(:,:,:)
  REAL(8)   , ALLOCATABLE   ::  cI2(:,:,:)
  REAL(8)   , ALLOCATABLE   ::  cI3(:,:,:)
  
  REAL(8)   , ALLOCATABLE   ::  cIH1(:,:)
  REAL(8)   , ALLOCATABLE   ::  cIH2(:,:)
  REAL(8)   , ALLOCATABLE   ::  cIH3(:,:)
  
  !--- Restriktion (Multigrid) ------------------------------------------------------------------------------- 
  REAL(8)   , ALLOCATABLE   ::  cR1 (:,:,:)
  REAL(8)   , ALLOCATABLE   ::  cR2 (:,:,:)
  REAL(8)   , ALLOCATABLE   ::  cR3 (:,:,:)
  
  REAL(8)   , ALLOCATABLE   ::  cRest1(:,:,:) ! TEST!!!
  REAL(8)   , ALLOCATABLE   ::  cRest2(:,:,:)
  REAL(8)   , ALLOCATABLE   ::  cRest3(:,:,:)
  
  REAL(8)   , ALLOCATABLE   ::  cRH1(:,:)
  REAL(8)   , ALLOCATABLE   ::  cRH2(:,:)
  REAL(8)   , ALLOCATABLE   ::  cRH3(:,:)
 

  !=================================
  !==== CT-based DNS
  !=================================

  REAL(8) , ALLOCATABLE, TARGET   ::  x_global(:), y_global(:), z_global(:)
  REAL(8)  , ALLOCATABLE, TARGET     :: ct_array(:)
  REAL(8)  , ALLOCATABLE, TARGET     :: ct_grid(:,:,:) 
  
  !===========================================================================================================
  !=== Gitterspezifikationen =================================================================================
  !===========================================================================================================
  !--- physiklische Koordinaten (global) ---------------------------------------------------------------------
  REAL(8)   , ALLOCATABLE, TARGET   ::  y1p(:), y1u(:)
  REAL(8)   , ALLOCATABLE, TARGET   ::  y2p(:), y2v(:)
  REAL(8)   , ALLOCATABLE, TARGET   ::  y3p(:), y3w(:)
  !--- C-FORTRAN interface ---
  TYPE(C_PTR), bind(C, name='_y1p_c') :: y1pptr
  TYPE(C_PTR), bind(C, name='_y2p_c') :: y2pptr
  TYPE(C_PTR), bind(C, name='_y3p_c') :: y3pptr

  TYPE(C_PTR), bind(C, name='_y1u_c') :: y1uptr
  TYPE(C_PTR), bind(C, name='_y2v_c') :: y2vptr
  TYPE(C_PTR), bind(C, name='_y3w_c') :: y3wptr

  !--- physiklische Koordinaten (Block) ----------------------------------------------------------------------
  REAL(8)   , ALLOCATABLE, TARGET   ::  x1p(:), x1u(:)
  REAL(8)   , ALLOCATABLE, TARGET   ::  x2p(:), x2v(:)
  REAL(8)   , ALLOCATABLE, TARGET   ::  x3p(:), x3w(:)
  !--- C-FORTRAN interface --- 
  TYPE(C_PTR), bind(C, name='_x1p_c') :: x1pptr
  TYPE(C_PTR), bind(C, name='_x2p_c') :: x2pptr
  TYPE(C_PTR), bind(C, name='_x3p_c') :: x3pptr

  TYPE(C_PTR), bind(C, name='_x1u_c') :: x1uptr
  TYPE(C_PTR), bind(C, name='_x2v_c') :: x2vptr
  TYPE(C_PTR), bind(C, name='_x3w_c') :: x3wptr

  !--- physiklische Koordinaten (Block, Multigrid) -----------------------------------------------------------
  REAL(8)   , ALLOCATABLE   ::  x1pR(:,:), x1uR(:,:)
  REAL(8)   , ALLOCATABLE   ::  x2pR(:,:), x2vR(:,:)
  REAL(8)   , ALLOCATABLE   ::  x3pR(:,:), x3wR(:,:)
  
  !--- Gitterweiten (global) ---------------------------------------------------------------------------------
  REAL(8)   , ALLOCATABLE, TARGET   ::  dy1p(:), dy1u(:)
  REAL(8)   , ALLOCATABLE, TARGET   ::  dy2p(:), dy2v(:)
  REAL(8)   , ALLOCATABLE, TARGET   ::  dy3p(:), dy3w(:)
  !--- C-FORTRAN interface --- 
  TYPE(C_PTR), bind(C, name='_dy1p_c') :: dy1pptr
  TYPE(C_PTR), bind(C, name='_dy2p_c') :: dy2pptr
  TYPE(C_PTR), bind(C, name='_dy3p_c') :: dy3pptr

  TYPE(C_PTR), bind(C, name='_dy1u_c') :: dy1uptr
  TYPE(C_PTR), bind(C, name='_dy2v_c') :: dy2vptr
  TYPE(C_PTR), bind(C, name='_dy3w_c') :: dy3wptr

  !--- Gitterweiten (Block) ----------------------------------------------------------------------------------
  REAL(8)   , ALLOCATABLE, TARGET   ::  dx1p(:), dx1u(:)
  REAL(8)   , ALLOCATABLE, TARGET   ::  dx2p(:), dx2v(:)
  REAL(8)   , ALLOCATABLE, TARGET   ::  dx3p(:), dx3w(:)
  
  REAL(8)   , ALLOCATABLE   ::  dx1DM(:), dx1pM(:), ddx1pM(:)
  REAL(8)   , ALLOCATABLE   ::  dx2DM(:), dx2pM(:), ddx2pM(:)
  REAL(8)   , ALLOCATABLE   ::  dx3DM(:), dx3pM(:), ddx3pM(:)
  
  REAL(8)   , ALLOCATABLE   ::  dx1GM(:), dx1uM(:), ddx1uM(:)
  REAL(8)   , ALLOCATABLE   ::  dx2GM(:), dx2vM(:), ddx2vM(:)
  REAL(8)   , ALLOCATABLE   ::  dx3GM(:), dx3wM(:), ddx3wM(:)  
  !--- C-FORTRAN interface --- 
  TYPE(C_PTR), bind(C, name='_dx1p_c') :: dx1pptr
  TYPE(C_PTR), bind(C, name='_dx2p_c') :: dx2pptr
  TYPE(C_PTR), bind(C, name='_dx3p_c') :: dx3pptr

  TYPE(C_PTR), bind(C, name='_dx1u_c') :: dx1uptr
  TYPE(C_PTR), bind(C, name='_dx2v_c') :: dx2vptr
  TYPE(C_PTR), bind(C, name='_dx3w_c') :: dx3wptr

  
  !===========================================================================================================
  !=== Arbeitsfelder =========================================================================================
  !===========================================================================================================
  !--- Geschwindigkeiten -------------------------------------------------------------------------------------
  REAL(8)   , ALLOCATABLE, TARGET   ::  vel(:,:,:,:)
  !--- C-FORTAN interface ---
  TYPE(C_PTR), bind(C, name='_vel_c') :: velptr

  !--- nicht-linearer Term -----------------------------------------------------------------------------------
  REAL(8)   , ALLOCATABLE, TARGET   ::  nl (:,:,:,:)
  !--- C-FORTRAN interface ---
  TYPE(C_PTR), bind(C, name='_nl_c') :: nlptr
  
  !--- Recht-Hand-Seite --------------------------------------------------------------------------------------
  REAL(8)   , ALLOCATABLE, TARGET   ::  rhs(:,:,:,:)
  !--- C-FORTRAN interface ---
  TYPE(C_PTR), bind(C, name='_rhs_c') :: rhsptr
  
  !--- Druck -------------------------------------------------------------------------------------------------
  REAL(8)   , ALLOCATABLE, TARGET   ::  pre(:,:,:)
  !--- C-FORTRAN interface ---
  TYPE(C_PTR), bind(C, name='_pre_c') :: preptr
 
  !--- The local (on-process) container for CT data 
  REAL(8)   , ALLOCATABLE, TARGET   ::  CT_geometry(:,:,:)

 
  !--- Ausfluss-RB (Geschwindigkeitsfeld) --------------------------------------------------------------------
  ! Da die RHS für die Konzentrationsfelder nicht über die Runge-Kutta-Zwischenschritte hinweg gespeichert 
  ! werden, müssen mindestens die zugehörigen Randbedingungen gespeichert werden.
  REAL(8)   , ALLOCATABLE   ::  bc11(:,:,:), nlbc11(:,:,:)
  REAL(8)   , ALLOCATABLE   ::  bc12(:,:,:), nlbc12(:,:,:)
  REAL(8)   , ALLOCATABLE   ::  bc13(:,:,:), nlbc13(:,:,:)
  
  REAL(8)   , ALLOCATABLE   ::  bc21(:,:,:), nlbc21(:,:,:)
  REAL(8)   , ALLOCATABLE   ::  bc22(:,:,:), nlbc22(:,:,:)
  REAL(8)   , ALLOCATABLE   ::  bc23(:,:,:), nlbc23(:,:,:)
  
  REAL(8)   , ALLOCATABLE   ::  bc31(:,:,:), nlbc31(:,:,:)
  REAL(8)   , ALLOCATABLE   ::  bc32(:,:,:), nlbc32(:,:,:)
  REAL(8)   , ALLOCATABLE   ::  bc33(:,:,:), nlbc33(:,:,:)
  
  REAL(8)   , ALLOCATABLE   ::  drift1(:,:,:)
  REAL(8)   , ALLOCATABLE   ::  drift2(:,:,:)
  REAL(8)   , ALLOCATABLE   ::  drift3(:,:,:)
  !--- Residuum ----------------------------------------------------------------------------------------------
  REAL(8)   , ALLOCATABLE   ::  res (:,:,:)
  
  !--- Druckgradient (eine Komponente) -----------------------------------------------------------------------
  
  REAL(8), target   , ALLOCATABLE   ::  gpre(:,:,:)
  !--- Gewichte für Divergenzfreiheit ------------------------------------------------------------------------
  REAL(8)   , ALLOCATABLE   ::  weight(:,:,:)
  
   !--- Null-Raum-Vektor --------------------------------------------------------------------------------------
  REAL(8)   , ALLOCATABLE   ::  psi    (:,:,:)
  REAL(8)   , ALLOCATABLE   ::  psi_vel(:,:,:,:)
  
  REAL(8)   , ALLOCATABLE   ::  psi_rel1 (:,:,:)
  REAL(8)   , ALLOCATABLE   ::  psi_rel2 (:,:,:)
  REAL(8)   , ALLOCATABLE   ::  psi_rel3 (:,:,:)
  REAL(8)   , ALLOCATABLE   ::  psi_rel4 (:,:,:)
  REAL(8)   , ALLOCATABLE   ::  psi_rel5 (:,:,:)
  REAL(8)   , ALLOCATABLE   ::  psi_rel6 (:,:,:)
  REAL(8)   , ALLOCATABLE   ::  psi_rel7 (:,:,:)
  REAL(8)   , ALLOCATABLE   ::  psi_rel8 (:,:,:)
  REAL(8)   , ALLOCATABLE   ::  psi_rel9 (:,:,:)
  REAL(8)   , ALLOCATABLE   ::  psi_rel10(:,:,:)
  REAL(8)   , ALLOCATABLE   ::  psi_rel11(:,:,:)
  REAL(8)   , ALLOCATABLE   ::  psi_rel12(:,:,:)
  REAL(8)   , ALLOCATABLE   ::  psi_rel13(:,:,:)
  REAL(8)   , ALLOCATABLE   ::  psi_rel14(:,:,:)
  REAL(8)   , ALLOCATABLE   ::  psi_rel15(:,:,:)
  
  REAL(8)   , ALLOCATABLE   ::  th11(:,:,:)
  REAL(8)   , ALLOCATABLE   ::  th12(:,:,:)
  REAL(8)   , ALLOCATABLE   ::  th13(:,:,:)
  
  REAL(8)   , ALLOCATABLE   ::  th21(:,:,:)
  REAL(8)   , ALLOCATABLE   ::  th22(:,:,:)
  REAL(8)   , ALLOCATABLE   ::  th23(:,:,:)
  
  REAL(8)   , ALLOCATABLE   ::  th31(:,:,:)
  REAL(8)   , ALLOCATABLE   ::  th32(:,:,:)
  REAL(8)   , ALLOCATABLE   ::  th33(:,:,:)
  
  !--- Multigrid ---------------------------------------------------------------------------------------------
  REAL(8)   , ALLOCATABLE   ::  vec1C (:,:,:)
  
  REAL(8)   , ALLOCATABLE   ::  vec2A (:,:,:)
  REAL(8)   , ALLOCATABLE   ::  vec2B (:,:,:)
  REAL(8)   , ALLOCATABLE   ::  vec2C (:,:,:)
  
  REAL(8)   , ALLOCATABLE   ::  vec3A (:,:,:)
  REAL(8)   , ALLOCATABLE   ::  vec3B (:,:,:)
  REAL(8)   , ALLOCATABLE   ::  vec3C (:,:,:)
  
  REAL(8)   , ALLOCATABLE   ::  vec4A (:,:,:)
  REAL(8)   , ALLOCATABLE   ::  vec4B (:,:,:)
  REAL(8)   , ALLOCATABLE   ::  vec4C (:,:,:)
  
  REAL(8)   , ALLOCATABLE   ::  vec5A (:,:,:)
  REAL(8)   , ALLOCATABLE   ::  vec5B (:,:,:)
  REAL(8)   , ALLOCATABLE   ::  vec5C (:,:,:)
  
  REAL(8)   , ALLOCATABLE   ::  vec6A (:,:,:)
  REAL(8)   , ALLOCATABLE   ::  vec6B (:,:,:)
  REAL(8)   , ALLOCATABLE   ::  vec6C (:,:,:)
  
  REAL(8)   , ALLOCATABLE   ::  vec7A (:,:,:)
  REAL(8)   , ALLOCATABLE   ::  vec7B (:,:,:)
  REAL(8)   , ALLOCATABLE   ::  vec7C (:,:,:)
  
  REAL(8)   , ALLOCATABLE   ::  vec8A (:,:,:)
  REAL(8)   , ALLOCATABLE   ::  vec8B (:,:,:)
  REAL(8)   , ALLOCATABLE   ::  vec8C (:,:,:)
  
  REAL(8)   , ALLOCATABLE   ::  vec9A (:,:,:)
  REAL(8)   , ALLOCATABLE   ::  vec9B (:,:,:)
  REAL(8)   , ALLOCATABLE   ::  vec9C (:,:,:)
  
  REAL(8)   , ALLOCATABLE   ::  vec10A(:,:,:)
  REAL(8)   , ALLOCATABLE   ::  vec10B(:,:,:)
  REAL(8)   , ALLOCATABLE   ::  vec10C(:,:,:)
  
  REAL(8)   , ALLOCATABLE   ::  vec11A(:,:,:)
  REAL(8)   , ALLOCATABLE   ::  vec11B(:,:,:)
  REAL(8)   , ALLOCATABLE   ::  vec11C(:,:,:)
  
  REAL(8)   , ALLOCATABLE   ::  vec12A(:,:,:)
  REAL(8)   , ALLOCATABLE   ::  vec12B(:,:,:)
  REAL(8)   , ALLOCATABLE   ::  vec12C(:,:,:)
  
  REAL(8)   , ALLOCATABLE   ::  vec13A(:,:,:)
  REAL(8)   , ALLOCATABLE   ::  vec13B(:,:,:)
  REAL(8)   , ALLOCATABLE   ::  vec13C(:,:,:)
  
  REAL(8)   , ALLOCATABLE   ::  vec14A(:,:,:)
  REAL(8)   , ALLOCATABLE   ::  vec14B(:,:,:)
  REAL(8)   , ALLOCATABLE   ::  vec14C(:,:,:)
  
  REAL(8)   , ALLOCATABLE   ::  vec15A(:,:,:)
  REAL(8)   , ALLOCATABLE   ::  vec15B(:,:,:)
  
  
  !--- BiCGstab / Richardson ---------------------------------------------------------------------------------
  REAL(8)   , target, ALLOCATABLE   ::  pp(:,:,:)
  REAL(8)   , ALLOCATABLE   ::  Ap(:,:,:)
  REAL(8)   , ALLOCATABLE   ::  rr(:,:,:)
  REAL(8)   , ALLOCATABLE   ::  rh(:,:,:)
  REAL(8)   , ALLOCATABLE   ::  Ar(:,:,:)
  REAL(8)   , ALLOCATABLE   ::  z1(:,:,:)
  REAL(8)   , ALLOCATABLE   ::  z2(:,:,:)
  
  !--- product_div_grad --------------------------------------------------------------------------------------
  REAL(8)   , target, ALLOCATABLE   ::  dig(:,:,:)
  
  !--- Hilfsfelder (Druckiteration) --------------------------------------------------------------------------
  ! - rhs wird auch in test_moment nochmals verwendet und kann daher in outer_iteration nicht belegt werden!
  ! - wird auch fuer interpolierte Geschwindigkeiten in rhs_NS und rhs_conc verwendet.
  REAL(8)   , ALLOCATABLE, TARGET   ::  work1(:,:,:)
  REAL(8)   , ALLOCATABLE, TARGET   ::  work2(:,:,:)
  REAL(8)   , ALLOCATABLE, TARGET   ::  work3(:,:,:)
  !--- C-FORTRAN interface ---
  TYPE(C_PTR), bind(C, name='_work1_c') :: work1ptr
  TYPE(C_PTR), bind(C, name='_work2_c') :: work2ptr
  TYPE(C_PTR), bind(C, name='_work3_c') :: work3ptr
  
  !--- Linienrelaxation --------------------------------------------------------------------------------------
  REAL(8)   , ALLOCATABLE   ::  vec1(:), dia1(:), SOR1(:), band1(:,:) ! TEST!!! siehe unten ...
  REAL(8)   , ALLOCATABLE   ::  vec2(:), dia2(:), SOR2(:), band2(:,:)
  REAL(8)   , ALLOCATABLE   ::  vec3(:), dia3(:), SOR3(:), band3(:,:)
  
  
  !===========================================================================================================
  !=== Indizierung (Intervallgrenzen, Verschiebungen) ========================================================
  !===========================================================================================================
  !--- Block-Index -------------------------------------------------------------------------------------------
  !INTEGER               ::  iB, jB, kB ! TEST!!! iB(1:3,1:n_grids_max) hierher verschieben ...
  
  !--- Indexverschiebung (Block --> global) ------------------------------------------------------------------
  INTEGER                ::  iShift, jShift, kShift
  
  !--- Domaingrösse (Periodizität-bereinigt) -----------------------------------------------------------------
  INTEGER                ::  dim1, dim2, dim33
  
  !--- Druck / Konzentrationen (inklusive Rand) --------------------------------------------------------------
  INTEGER                ::  S1p, S2p, S3p
  INTEGER                ::  N1p, N2p, N3p
  
  !--- Geschwindigkeiten (inklusive Rand) --------------------------------------------------------------------
  INTEGER                ::  S11B, S21B, S31B
  INTEGER                ::  S12B, S22B, S32B
  INTEGER                ::  S13B, S23B, S33B
  
  INTEGER                ::  N11B, N21B, N31B
  INTEGER                ::  N12B, N22B, N32B
  INTEGER                ::  N13B, N23B, N33B
  
  !--- Geschwindigkeiten (exklusive Rand) --------------------------------------------------------------------
  INTEGER                ::  S11, S21, S31
  INTEGER                ::  S12, S22, S32
  INTEGER                ::  S13, S23, S33
  
  INTEGER                ::  N11, N21, N31
  INTEGER                ::  N12, N22, N32
  INTEGER                ::  N13, N23, N33
  
  !--- grobe Gitter (Multigrid, INklusive Rand) --------------------------------------------------------------
  INTEGER                ::  S1R, S2R, S3R
  INTEGER                ::  d1R, d2R, d3R
  
  !--- grobe Gitter (Multigrid, EXklusive Rand) --------------------------------------------------------------
  INTEGER                ::  S11R, S22R, S33R
  INTEGER                ::  d11R, d22R, d33R
  
  !--- Ueberlappungskonvention der Blöcke (Multigrid, siehe mod_setup) ---------------------------------------
  INTEGER, PARAMETER     ::  ls1 = -1
  INTEGER, PARAMETER     ::  ls2 = -1
  INTEGER, PARAMETER     ::  ls3 = -1
  !--- C-FORTRAN interface ---
  INTEGER, PROTECTED, bind(C, name='_ls1_c') :: ls1_c = ls1
  INTEGER, PROTECTED, bind(C, name='_ls2_c') :: ls2_c = ls2
  INTEGER, PROTECTED, bind(C, name='_ls3_c') :: ls3_c = ls3


  !--- Austauschrichtung (Multigrid) -------------------------------------------------------------------------
  ! ex = -1: unten <--  oben
  ! ex =  0: unten <--> oben
  ! ex =  1: unten  --> oben
  INTEGER                ::  ex1, ex2, ex3
  
  
  !===========================================================================================================
  !=== Randbedingungen =======================================================================================
  !===========================================================================================================
  !                              _
  !    Symmetrie-RB:   BC = -2    |
  !    periodische RB: BC = -1    |- symmetrische, zentrale Stencils
  !    Nachbar-Block:  BC =  0   _|
  !    Dirichlet-RB:   BC =  1    |
  !    Neumann-RB:     BC =  2    |- schiefe, nicht-zentrale Stencils
  !    Robin-RB:       BC =  3   _|
  !
  !--- global ------------------------------------------------------------------------------------------------
  LOGICAL                ::  outlet(1:3,1:2,1:3)
  
  INTEGER                ::  BC_1L_global, BC_1U_global
  INTEGER                ::  BC_2L_global, BC_2U_global
  INTEGER                ::  BC_3L_global, BC_3U_global
  
  !--- lokal (Block) -----------------------------------------------------------------------------------------
  INTEGER                ::  BC_1L, BC_1U
  INTEGER                ::  BC_2L, BC_2U
  INTEGER                ::  BC_3L, BC_3U

  !--- field properties --------------------------------------------------------------------------------------
  INTEGER                ::  n_gather(1:3,1:n_grids_max)
  INTEGER                ::  NN (1:3,1:n_grids_max)
  INTEGER                ::  NB (1:3,1:n_grids_max)
  INTEGER, TARGET        ::  iB (1:3,1:n_grids_max)
  INTEGER                ::  SNF(1:2,1:3,1:n_grids_max)
  INTEGER                ::  SNB(1:2,1:3,1:n_grids_max)
  INTEGER                ::  BC (1:2,1:3,1:n_grids_max)
  INTEGER                ::  ngb(1:2,1:3,1:n_grids_max)
  INTEGER                ::  comm1(1:n_grids_max), comm2(1:n_grids_max)
  INTEGER                ::  rankc2(1:n_grids_max)
  LOGICAL                ::  participate_yes(1:n_grids_max)
  INTEGER, ALLOCATABLE   ::  recvR(  :,:), recvI(  :,:)
  INTEGER, ALLOCATABLE   ::  dispR(  :,:), dispI(  :,:)
  INTEGER, ALLOCATABLE   ::  offsR(:,:,:), offsI(:,:,:)
  INTEGER, ALLOCATABLE   ::  sizsR(:,:,:), sizsI(:,:,:)
  !--- C-FORTRAN interface ---
  TYPE(C_PTR), bind(C, name='_ib_c') :: ibptr

  
  !===========================================================================================================
  !=== physikalische Parameter ===============================================================================
  !===========================================================================================================
  REAL(8)                   ::  L1, L2, L3
  REAL(8)                   ::  Re
  
  
  !===========================================================================================================
  !=== numerische Parameter ==================================================================================
  !===========================================================================================================
  !--- allgemein ---------------------------------------------------------------------------------------------
  REAL(8)                   ::  CFL
  REAL(8)                   ::  time, dtime, subtime, time_start, time_end, dtime_max, dtime0, dtime_old
  INTEGER                ::  timestep, timestep_old, substep, n_timesteps
  LOGICAL                ::  mapping_yes, upwind_yes
  LOGICAL                ::  Euler_yes, Stokes_yes, twostep_yes
  LOGICAL, PARAMETER     ::  filter_BC_yes = .TRUE. ! TEST!!!
  INTEGER                ::  timeint_mode, forcing_mode
  INTEGER                ::  bulkflow_dir
  INTEGER                ::  n_lp_vel    , n_hp_vel
  REAL(8)                   ::  chi_vel
  REAL(8)                   ::  vel_bulk ! TEST!!!
  !--- C-FORTRAN interface ---
  LOGICAL, PROTECTED, bind(C, name='_filter_BC_yes_c') :: filter_bc_yes_c = filter_BC_yes
  
  !--- Runge-Kutta-Koeffizienten -----------------------------------------------------------------------------
  REAL(8)   , PARAMETER     ::  aRK(1:3) = (/8./15.,  5./12., 3./ 4./)
  REAL(8)   , PARAMETER     ::  bRK(1:3) = (/  0.  ,-17./60.,-5./12./)
  INTEGER, PARAMETER     ::  RK_steps = 3
  !--- C-FORTRAN interface ---
  REAL(8), PROTECTED, TARGET  ::  ark_c(1:3) = ark(1:3)
  REAL(8), PROTECTED, TARGET  ::  brk_c(1:3) = brk(1:3)
  TYPE(C_PTR), bind(C, name='_ark_c') :: arkptr
  TYPE(C_PTR), bind(C, name='_brk_c') :: brkptr
  INTEGER, PROTECTED, bind(C, name='_rk_steps_c') :: rk_steps_c = RK_steps

  !--- look-up table fuer Stabilitaetsgebiet der Zeitintegration (angle = pi/2,pi) ---------------------------
  REAL(8)   , PARAMETER     ::  stabilitylimit(0:40) = (/1.732050813, 1.943689093, 2.089210537, 2.201001743,  &
                                      &               2.290031261, 2.361554127, 2.418567407, 2.462989697,  &
                                      &               2.496169963, 2.519146008, 2.532795254, 2.537935070,  &
                                      &               2.535397854, 2.526091466, 2.511046932, 2.491448818,  &
                                      &               2.468639045, 2.444084180, 2.419302172, 2.395757241,  &
                                      &               2.374745783, 2.357302135, 2.344145473, 2.335672458,  &
                                      &               2.331985072, 2.332936948, 2.338183901, 2.347230689,  &
                                      &               2.359471631, 2.374225928, 2.390769340, 2.408363261,  &
                                      &               2.426281290, 2.443832601, 2.460381269, 2.475360992,  &
                                      &               2.488285197, 2.498753090, 2.506452564, 2.511161051,  &
                                      &               2.512745327 /)
  INTEGER, PARAMETER     ::  n_stab = SIZE(stabilitylimit)
  !--- C-FORTRAN interface ---
  INTEGER, PROTECTED, bind(C, name='_n_stab_c') :: n_stab_c = n_stab

  !--- Helmholtz-Vorfaktoren ---------------------------------------------------------------------------------
  REAL(8), target           ::  thetaL, multL
  
  !--- zeitliche Steuerung -----------------------------------------------------------------------------------
  INTEGER                ::  Int_dtime, Int_lev_pre
  
  INTEGER                ::  stride_large(1:3), stride_med(1:3), stride_small(1:3)
  LOGICAL                ::  write_large, write_med, write_small
  REAL(8)                   ::  time_out_scal, dtime_out_scal
  REAL(8)                   ::  time_out_vect, dtime_out_vect
  
  LOGICAL                ::  write_out_scal, write_out_vect
  LOGICAL                ::  new_dtime, finish_yes
    
  INTEGER                ::  write_count
  INTEGER                ::  restart
  CHARACTER(LEN=3)       ::  restart_char
  
  INTEGER                ::  n_conc_old

  
  !===========================================================================================================
  !=== weitere Steuerungsoptionen ============================================================================
  !===========================================================================================================
  INTEGER                ::  task
  LOGICAL                ::  read_nullspace_yes
  LOGICAL                ::  nullspace_yes, nullspace_coarse_yes
  LOGICAL                ::  nullspace_ortho_yes
  
  LOGICAL                ::  write_stout_yes
  LOGICAL                ::  log_iteration_yes
  LOGICAL                ::  write_restart_yes
  LOGICAL                ::  write_lambda2_yes
  LOGICAL                ::  write_test_yes
  
  !--- globale Laufindizes -----------------------------------------------------------------------------------
  INTEGER                ::  direction
  
  !--- explizite Behandlung von Ecken bei Dirichlet-Randbedingungen ------------------------------------------
  ! (Hintergrund: Der Druck ist an diesen Orten unbestimmt, so dass er dort künstlich zu Null gesetzt wird.)
  LOGICAL, PARAMETER     ::  corner_yes = .TRUE.
  !--- C-FORTRAN interface ---
  LOGICAL, PROTECTED, bind(C, name='_corner_yes_c') :: corner_yes_c = corner_yes

  !--- Systemzeit --------------------------------------------------------------------------------------------
  INTEGER                ::  elatime, ctime(1:8)
  INTEGER                ::  day, hour, minu, sec, msec
  
  
  !===========================================================================================================
  !=== Iterationsparameter ===================================================================================
  !===========================================================================================================
  !--- Abbruchkriterium / Absolute Genauigkeit der Geschwindigkeiten -----------------------------------------
  REAL(8)                   ::  epsU, epsU0
  
  !--- Glaetter ----------------------------------------------------------------------------------------------
  LOGICAL                ::  Jacobi_yes
  
  !--- max. Anzahl Iterationen -------------------------------------------------------------------------------
  INTEGER                ::  n_it_outer
  INTEGER                ::  n_it_Poisson
  INTEGER                ::  n_it_Helmh_vel
  
  !--- erwartete Konvergenzrate (äussere Iteration) ----------------------------------------------------------
  REAL(8)                   ::  precRatio0 (1:RK_steps)
  REAL(8)                   ::  precOffset0(1:RK_steps)
  REAL(8)   , ALLOCATABLE   ::  precRatio  (:,:)
  REAL(8)   , ALLOCATABLE   ::  precOffset (:,:)
  
  !--- Null-Initialisierung (äussere Iteration) --------------------------------------------------------------
  LOGICAL, TARGET                ::  init_pre(1:RK_steps), init_vel(1:RK_steps)!, init_conc(1:RK_steps)
  !--- C-Fortran interface ---
  TYPE(C_PTR), bind(C, name='_init_pre_c') :: init_preptr
  
  !--- Vorkonditionierung (Multigrid) ------------------------------------------------------------------------
  !--- Vorkonditionierung (Multigrid) ------------------------------------------------------------------------
  INTEGER                ::  precond_outer
  INTEGER                ::  precond_Poisson
  INTEGER                ::  precond_Helmh_vel
  !INTEGER                ::  precond_Helmh_conc
  
  !--- Anzahl Glättungen pro Gitterlevel (Multigrid) ---------------------------------------------------------
  INTEGER                ::  n_relax_down, n_relax_up, n_relax_bottom
  
  !--- implizite Richtungen bei Linienrelaxation (Multigrid) -------------------------------------------------
  INTEGER                ::  impl_dir(1:3)
  
  !--- Anzahl Glättungen pro Gitterlevel (Multigrid) ---------------------------------------------------------
  LOGICAL                ::  weighting_yes
  
  
  !===========================================================================================================
  !=== Iterationsstatistik ===================================================================================
  !===========================================================================================================
  REAL(8)                   ::  dtime_average
  REAL(8)                   ::  max_div_init(1:2)
  INTEGER                ::  number_poisson
  
  !--- Zähler ------------------------------------------------------------------------------------------------
  INTEGER                ::  countO(1:RK_steps)
  INTEGER                ::  countP(1:RK_steps,1:2)
  INTEGER                ::  countH(1:RK_steps,1:3)
  
  !--- Konvergenzrate ----------------------------------------------------------------------------------------
  REAL(8)                   ::  ratioO(1:RK_steps)
  REAL(8)                   ::  ratioH(1:RK_steps,1:3)
  REAL(8)                   ::  ratioP(1:RK_steps,1:2)
  
  
  !===========================================================================================================
  !=== MPI ===================================================================================================
  !===========================================================================================================
  !--- Kommunikatoren ----------------------------------------------------------------------------------------
  INTEGER                ::  COMM_CART
  
  INTEGER                ::  COMM_SLICE1, COMM_BAR1
  INTEGER                ::  COMM_SLICE2, COMM_BAR2
  INTEGER                ::  COMM_SLICE3, COMM_BAR3
  
  !--- Dimension und Position der Blöcke innerhalb der Kommunikatoren (Gitter-Indizes) -----------------------
  ! (für MPI_GATHERv, MPI_ALLGATHERv, vgl. iShift, jShift, kShift)
  INTEGER, ALLOCATABLE   ::  bar1_size(:), bar1_offset(:)
  INTEGER, ALLOCATABLE   ::  bar2_size(:), bar2_offset(:)
  INTEGER, ALLOCATABLE   ::  bar3_size(:), bar3_offset(:)
  
  !--- Ränge der Prozesse ------------------------------------------------------------------------------------
  INTEGER                ::  rank
  INTEGER                ::  rank_bar1, rank_slice1
  INTEGER                ::  rank_bar2, rank_slice2
  INTEGER                ::  rank_bar3, rank_slice3
  
  !--- Ränge der Nachbarprozesse (in kartesischem Gitter) ----------------------------------------------------
  INTEGER                ::  rank1L, rank1U
  INTEGER                ::  rank2L, rank2U
  INTEGER                ::  rank3L, rank3U
  
  !--- Error-Handle ------------------------------------------------------------------------------------------
  INTEGER                ::  merror
  
  !--- Request-Handles ---------------------------------------------------------------------------------------
  ! (müssen offenbar global angelegt werden) 
  INTEGER                ::  req1L, req1U
  INTEGER                ::  req2L, req2U
  INTEGER                ::  req3L, req3U
  
  !===========================================================================================================
  !=== HDF5 ==================================================================================================
  !===========================================================================================================
  INTEGER                ::  herror
  
  !*** bbecsek 2015: IMMERSED BOUNDARY PARAMETERS ************************************************************
  !REAL(8)                   ::  delta_s                 !< boundary mesh width
  REAL(8)                   ::  reach                   !< reach of interpolation reach
  INTEGER                ::  M_bound                 !< number of grid points on boundary
  INTEGER                ::  M_elems                 !< number of finite elements
  REAL(8)   , ALLOCATABLE   ::  yb(:,:)                 !< global Lagrangian boundary coordinates (allocated in init_ib)
  REAL(8)   , ALLOCATABLE   ::  xb(:,:)                 !< reference Lagrangian boundary coordinates (if needed) 
                                                     !! (allocated in init_ib)
  REAL(8)   , ALLOCATABLE   ::  db(:,:)                 !< Lagrangian displacement field
  REAL(8)   , ALLOCATABLE   ::  ub(:,:)                 !< boundary velocities
  REAL(8)   , ALLOCATABLE   ::  fb(:,:)                 !< Lagrangian force (allocated in init_ib)
  REAL(8)   , ALLOCATABLE, TARGET   ::  fd(:,:,:,:)             !< force density
  REAL(8)   , ALLOCATABLE   ::  container(:,:)          !< temporal container for the time-lagged RK vel-term (may be omitted 
                                                     !! later on by either performing the respective additions before and after 
                                                     !! the boundary update or by using the 'buffer' containers that already 
                                                     !! store the velocity fields
  !REAL(8)                   ::  kappa                   !< proportionality constant used for the Lagrangian force computation
  INTEGER                ::  ddf_type                !< flag for choosing the spreading and interpolation kernel function
                                                     !! (ddf: discrete delta function)
  LOGICAL                ::  IB_on                   !< specifies whether to immerse a boundary
  REAL(8)                   ::  mu_blood                !< dynamic viscosity of blood
  REAL(8)                   ::  rho_blood               !< bulk density of blood
  REAL(8)                   ::  E_mod                   !< Young's modulus of elasticity (linear) 
  REAL(8)                   ::  nu_poiss                !< Poisson's ratio  

  REAL(8)                   ::  rho_solid               !< density of solid
  REAL(8)   , ALLOCATABLE   ::  strain(:,:)             !< element strains
  REAL(8)   , ALLOCATABLE   ::  stress(:,:)             !< element stress
  REAL(8)   , ALLOCATABLE   ::  k_stiff(:,:)            !< global stiffness matrix
  REAL(8)   , ALLOCATABLE   ::  m_mass(:,:)             !< global mass matrix
  REAL(8)   , ALLOCATABLE   ::  c_damp(:,:)             !< global damping matrix (C=alpha*M + beta*K)
  REAL(8)                   ::  C(3,3)                  !< matrix linking strain to stress


  INTEGER, ALLOCATABLE   ::  elems(:,:)              !< Finite elements

  REAL(8)                   ::  L_ref                   !< reference length scale
  REAL(8)                   ::  U_ref                   !< reference velocity

  INTEGER                ::  req_bcast(3), req_red(3)      !< request handles

  INTEGER                ::  excl_p(3)               !< points to exclude 

  REAL(8)                   ::  delta_b1,delta_b2,delta_b3                 !< boundary discretizations

  INTEGER                ::  total_ranks             !< total number of ranks

  INTEGER                ::  M_ebcs                  !< number of prescribed essential boundary conditions

  REAL(8)   , ALLOCATABLE   ::  ebcs(:,:)               !< container for essential boundary conditions

  REAL(8)   , ALLOCATABLE   ::  node_vol(:)             !< volume (area) associated with node i

  INTEGER                :: COMM_LOCAL               !< communicator for local group
  INTEGER                :: COMM_INTER               !< intercommunicator solid-fluid

  REAL(8)                   ::  max_ev_freq             !< highest ev (frequency) of structure

  LOGICAL                ::  write_force_yes         !< whether to output force density or not

  LOGICAL                ::  write_xdmf_yes          !< whether to output .xmf files
  LOGICAL                ::  scale_output_yes        !< whether to scale all values with U_ref, L_ref (dimensionalize)

  !--- Oscillator Aorta Model ---
  LOGICAL                ::  fem_yes
  REAL(8)                   ::  k_x, k_y, k_z
  REAL(8)                   ::  aorta_start, aorta_end, aorta_radius

  !--- C-FORTRAN interface ---
  TYPE(C_PTR), bind(C, name='_fd_c') :: fdptr


#else
  
  
  
  
  !===========================================================================================================
  !=== Raemliche Dimensionen =================================================================================
  !===========================================================================================================
  ! 2D wird ueber M3==2 eingeschaltet
  INTEGER                ::  dimens
  
  
  !===========================================================================================================
  !=== Domain- und Blockspezifikationen ======================================================================
  !===========================================================================================================
  !--- zulaessige Blockgroessen ------------------------------------------------------------------------------
  ! n  x 1   2   4   8   16   32   64   128   256  512   ...
  !--------------------------------------------------------------
  ! n2 = 2,  4,  8, 16,  32,  64, 128,  256,  512, 1024, ...
  ! n3 = 3,  6, 12, 24,  48,  96, 192,  384,  768, 1536, ...
  ! n5 = 5, 10, 20, 40,  80, 160, 320,  640, 1280, 2560, ...
  ! n7 = 7, 14, 28, 56, 112, 224, 448,  896, 1792, 3584, ...
  ! n9 = 9, 18, 36, 72, 144, 288, 576, 1152, 2304, 4608, ...
  ! ...
  INTEGER, PARAMETER     ::  N1 = 1+(M1-1)/NB1
  INTEGER, PARAMETER     ::  N2 = 1+(M2-1)/NB2
  INTEGER, PARAMETER     ::  N3 = 1+(M3-1)/NB3
  
  !--- Anzahl grobe Gitter (Multigrid) -----------------------------------------------------------------------
  INTEGER, PARAMETER     ::  n_grids_max = 15
  INTEGER                ::  n_grids, n_grids_limit
  
  !--- Dimensionen -------------------------------------------------------------------------------------------
  INTEGER, PARAMETER     ::  dim_ncb1c = SIZE(ncb1c)
  INTEGER, PARAMETER     ::  dim_ncb1g = SIZE(ncb1g)
  INTEGER, PARAMETER     ::  dim_ncb1d = SIZE(ncb1d)
  
  INTEGER, PARAMETER     ::  dim_ncb2c = SIZE(ncb2c)
  INTEGER, PARAMETER     ::  dim_ncb2g = SIZE(ncb2g)
  INTEGER, PARAMETER     ::  dim_ncb2d = SIZE(ncb2d)
                         
  INTEGER, PARAMETER     ::  dim_ncb3c = SIZE(ncb3c)
  INTEGER, PARAMETER     ::  dim_ncb3g = SIZE(ncb3g)
  INTEGER, PARAMETER     ::  dim_ncb3d = SIZE(ncb3d)
  
  !--- Anzahl Stencil-Koeffizienten (Feld)------------- ------------------------------------------------------
  ! Anzahl der Koeffizienten im Feld (zentrale Differenzen angenommen):
  INTEGER, PARAMETER     ::  nc1c = ncb1c(dim_ncb1c)
  INTEGER, PARAMETER     ::  nc1s = ncb1g(dim_ncb1g)
                         
  INTEGER, PARAMETER     ::  nc2c = ncb2c(dim_ncb2c)
  INTEGER, PARAMETER     ::  nc2s = ncb2g(dim_ncb2g)
                         
  INTEGER, PARAMETER     ::  nc3c = ncb3c(dim_ncb3c)
  INTEGER, PARAMETER     ::  nc3s = ncb3g(dim_ncb3g)
  
  
  !===========================================================================================================
  !=== Intervallgrenzen der Differenzen-Koeffizienten-Arrays =================================================
  !===========================================================================================================
  !--- zentral -----------------------------------------------------------------------------------------------
  INTEGER, PARAMETER     ::  b1U = nc1s/2
  INTEGER, PARAMETER     ::  b2U = nc2s/2
  INTEGER, PARAMETER     ::  b3U = nc3s/2
  
  INTEGER, PARAMETER     ::  b1L = -b1U
  INTEGER, PARAMETER     ::  b2L = -b2U
  INTEGER, PARAMETER     ::  b3L = -b3U
  
  !--- upwind (nicht-linear) ---------------------------------------------------------------------------------
  ! (aktuell wird nicht zwischen auf- und abwärtsgerichteten Stencils unterschieden, um auch am Rand arbeiten
  !  zu können, wo KEINE upwind-Differenzierung verwendet wird)
  INTEGER, PARAMETER     ::  n1L = b1L
  INTEGER, PARAMETER     ::  n2L = b2L
  INTEGER, PARAMETER     ::  n3L = b3L
  
  INTEGER, PARAMETER     ::  n1U = b1U
  INTEGER, PARAMETER     ::  n2U = b2U
  INTEGER, PARAMETER     ::  n3U = b3U
  
  !--- Divergenz ---------------------------------------------------------------------------------------------
  INTEGER, PARAMETER     ::  d1L = b1L
  INTEGER, PARAMETER     ::  d2L = b2L
  INTEGER, PARAMETER     ::  d3L = b3L
  
  INTEGER, PARAMETER     ::  d1U = b1U-1
  INTEGER, PARAMETER     ::  d2U = b2U-1
  INTEGER, PARAMETER     ::  d3U = b3U-1
  
  !--- Gradient ----------------------------------------------------------------------------------------------
  INTEGER, PARAMETER     ::  g1L = b1L+1
  INTEGER, PARAMETER     ::  g2L = b2L+1
  INTEGER, PARAMETER     ::  g3L = b3L+1
  
  INTEGER, PARAMETER     ::  g1U = b1U
  INTEGER, PARAMETER     ::  g2U = b2U
  INTEGER, PARAMETER     ::  g3U = b3U
  
  
  !===========================================================================================================
  !=== Differenzen-Koeffizienten-Arrays ======================================================================
  !===========================================================================================================
  !--- 1. Ableitung (zentral) --------------------------------------------------------------------------------
  REAL(8)                   ::  cp1  (b1L:b1U,0:N1)
  REAL(8)                   ::  cp2  (b2L:b2U,0:N2)
  REAL(8)                   ::  cp3  (b3L:b3U,0:N3)
  
  REAL(8)                   ::  cu1  (b1L:b1U,0:N1)
  REAL(8)                   ::  cv2  (b2L:b2U,0:N2)
  REAL(8)                   ::  cw3  (b3L:b3U,0:N3)
  
  !--- 1. Ableitung (upwind) ---------------------------------------------------------------------------------
  REAL(8)                   ::  cNp1D(n1L:n1U,0:N1)
  REAL(8)                   ::  cNp2D(n2L:n2U,0:N2)
  REAL(8)                   ::  cNp3D(n3L:n3U,0:N3)
  
  REAL(8)                   ::  cNp1U(n1L:n1U,0:N1)
  REAL(8)                   ::  cNp2U(n2L:n2U,0:N2)
  REAL(8)                   ::  cNp3U(n3L:n3U,0:N3)
  
  REAL(8)                   ::  cNu1D(n1L:n1U,0:N1)
  REAL(8)                   ::  cNv2D(n2L:n2U,0:N2)
  REAL(8)                   ::  cNw3D(n3L:n3U,0:N3)
  
  REAL(8)                   ::  cNu1U(n1L:n1U,0:N1)
  REAL(8)                   ::  cNv2U(n2L:n2U,0:N2)
  REAL(8)                   ::  cNw3U(n3L:n3U,0:N3)

  !--- Divergenz ---------------------------------------------------------------------------------------------
  REAL(8)                   ::  cDu1 (d1L:d1U,0:N1)
  REAL(8)                   ::  cDv2 (d2L:d2U,0:N2)
  REAL(8)                   ::  cDw3 (d3L:d3U,0:N3)
  
  !--- Divergenz (transponiert) ------------------------------------------------------------------------------
  REAL(8)                   ::  cDu1T(g1L:g1U,0:N1)
  REAL(8)                   ::  cDv2T(g2L:g2U,0:N2)
  REAL(8)                   ::  cDw3T(g3L:g3U,0:N3)
  
  !--- Gradient ----------------------------------------------------------------------------------------------
  REAL(8)                   ::  cGp1 (g1L:g1U,0:N1)
  REAL(8)                   ::  cGp2 (g2L:g2U,0:N2)
  REAL(8)                   ::  cGp3 (g3L:g3U,0:N3)
  
  !--- Gradient (transponiert) -------------------------------------------------------------------------------
  REAL(8)                   ::  cGp1T(d1L:d1U,0:N1)
  REAL(8)                   ::  cGp2T(d2L:d2U,0:N2)
  REAL(8)                   ::  cGp3T(d3L:d3U,0:N3)
  
  !--- 2. Ableitung (zentral) --------------------------------------------------------------------------------
  REAL(8)                   ::  cp11 (b1L:b1U,0:N1)
  REAL(8)                   ::  cp22 (b2L:b2U,0:N2)
  REAL(8)                   ::  cp33 (b3L:b3U,0:N3)
  
  REAL(8)                   ::  cu11 (b1L:b1U,0:N1)
  REAL(8)                   ::  cv22 (b2L:b2U,0:N2)
  REAL(8)                   ::  cw33 (b3L:b3U,0:N3)

  !--- Interpolation ----------------------------------------------------------------------------------------- 
  REAL(8)                   ::  cIpu(g1L:g1U,0:N1)
  REAL(8)                   ::  cIpv(g2L:g2U,0:N2)
  REAL(8)                   ::  cIpw(g3L:g3U,0:N3)
  
  REAL(8)                   ::  cIup(d1L:d1U,0:N1)
  REAL(8)                   ::  cIvp(d2L:d2U,0:N2)
  REAL(8)                   ::  cIwp(d3L:d3U,0:N3)

  !--- Filter ------------------------------------------------------------------------------------------------
  REAL(8)                   ::  cFp1(b1L:b1U,0:N1)
  REAL(8)                   ::  cFp2(b2L:b2U,0:N2)
  REAL(8)                   ::  cFp3(b3L:b3U,0:N3)
  
  REAL(8)                   ::  cFu1(b1L:b1U,0:N1)
  REAL(8)                   ::  cFv2(b2L:b2U,0:N2)
  REAL(8)                   ::  cFw3(b3L:b3U,0:N3)

  !--- Integrator (nur für Druckgitter) ----------------------------------------------------------------------
  REAL(8)                   ::  cInt1(b1L:b1U,0:N1)
  REAL(8)                   ::  cInt2(b2L:b2U,0:N2)
  REAL(8)                   ::  cInt3(b3L:b3U,0:N3)
  
  !--- 2. Ableitung (Multigrid) ------------------------------------------------------------------------------ 
  ! Anmerkung: Die Koeffizientensätze unterscheiden sich z.T. lediglich durch die Randbedingungen.
  REAL(8)                   ::  cp11R(-1:1,0:N1,1:n_grids_max)
  REAL(8)                   ::  cp22R(-1:1,0:N2,1:n_grids_max)
  REAL(8)                   ::  cp33R(-1:1,0:N3,1:n_grids_max)
  
  REAL(8)                   ::  cu11R(-1:1,0:N1,1:n_grids_max)
  REAL(8)                   ::  cv22R(-1:1,0:N2,1:n_grids_max)
  REAL(8)                   ::  cw33R(-1:1,0:N3,1:n_grids_max)

  REAL(8)                   ::  cdg1 (-1:1,1:N1,1:n_grids_max)
  REAL(8)                   ::  cdg2 (-1:1,1:N2,1:n_grids_max)
  REAL(8)                   ::  cdg3 (-1:1,1:N3,1:n_grids_max)
  
  !--- Interpolation (Multigrid) ----------------------------------------------------------------------------- 
  REAL(8)                   ::  cI1(1:2,1:N1,1:n_grids_max)
  REAL(8)                   ::  cI2(1:2,1:N2,1:n_grids_max)
  REAL(8)                   ::  cI3(1:2,1:N3,1:n_grids_max)
  
  REAL(8)                   ::  cIH1(1:2,0:N1)
  REAL(8)                   ::  cIH2(1:2,0:N2)
  REAL(8)                   ::  cIH3(1:2,0:N3)
  
  !--- Restriktion (Multigrid) ------------------------------------------------------------------------------- 
  REAL(8)                   ::  cR1 (-1:1,1:N1,2:n_grids_max)
  REAL(8)                   ::  cR2 (-1:1,1:N2,2:n_grids_max)
  REAL(8)                   ::  cR3 (-1:1,1:N3,2:n_grids_max)
  
  REAL(8)                   ::  cRest1(b1L:b1U,0:N1,1:n_grids_max-1) ! TEST!!!
  REAL(8)                   ::  cRest2(b2L:b2U,0:N2,1:n_grids_max-1)
  REAL(8)                   ::  cRest3(b3L:b3U,0:N3,1:n_grids_max-1)
  
  REAL(8)                   ::  cRH1(1:2,1:N1)
  REAL(8)                   ::  cRH2(1:2,1:N2)
  REAL(8)                   ::  cRH3(1:2,1:N3)
  
  
  !===========================================================================================================
  !=== Gitterspezifikationen =================================================================================
  !===========================================================================================================
  !--- physiklische Koordinaten (global) ---------------------------------------------------------------------
  REAL(8)                   ::  y1p(1:M1), y1u(0:M1)
  REAL(8)                   ::  y2p(1:M2), y2v(0:M2)
  REAL(8)                   ::  y3p(1:M3), y3w(0:M3)
  
  !--- physiklische Koordinaten (Block) ----------------------------------------------------------------------
  REAL(8)                   ::  x1p(b1L:(N1+b1U)), x1u(b1L:(N1+b1U))
  REAL(8)                   ::  x2p(b2L:(N2+b2U)), x2v(b2L:(N2+b2U))
  REAL(8)                   ::  x3p(b3L:(N3+b3U)), x3w(b3L:(N3+b3U))
  
  !--- physiklische Koordinaten (Block, Multigrid) -----------------------------------------------------------
  REAL(8)                   ::  x1pR(b1L:(N1+b1U),1:n_grids_max), x1uR(b1L:(N1+b1U),1:n_grids_max)
  REAL(8)                   ::  x2pR(b2L:(N2+b2U),1:n_grids_max), x2vR(b2L:(N2+b2U),1:n_grids_max)
  REAL(8)                   ::  x3pR(b3L:(N3+b3U),1:n_grids_max), x3wR(b3L:(N3+b3U),1:n_grids_max)
  
  !--- Gitterweiten (global) ---------------------------------------------------------------------------------
  REAL(8)                   ::  dy1p(1:M1), dy1u(0:M1)
  REAL(8)                   ::  dy2p(1:M2), dy2v(0:M2)
  REAL(8)                   ::  dy3p(1:M3), dy3w(0:M3)
  
  !--- Gitterweiten (Block) ----------------------------------------------------------------------------------
  REAL(8)                   ::  dx1p(1:N1), dx1u(0:N1)
  REAL(8)                   ::  dx2p(1:N2), dx2v(0:N2)
  REAL(8)                   ::  dx3p(1:N3), dx3w(0:N3)
  
  REAL(8)                   ::  dx1DM(1:N1), dx1pM(1:N1), ddx1pM(1:N1)
  REAL(8)                   ::  dx2DM(1:N2), dx2pM(1:N2), ddx2pM(1:N2)
  REAL(8)                   ::  dx3DM(1:N3), dx3pM(1:N3), ddx3pM(1:N3)
  
  REAL(8)                   ::  dx1GM(0:N1), dx1uM(0:N1), ddx1uM(0:N1)
  REAL(8)                   ::  dx2GM(0:N2), dx2vM(0:N2), ddx2vM(0:N2)
  REAL(8)                   ::  dx3GM(0:N3), dx3wM(0:N3), ddx3wM(0:N3)

  
  !===========================================================================================================
  !=== Arbeitsfelder =========================================================================================
  !===========================================================================================================
  !--- Geschwindigkeiten -------------------------------------------------------------------------------------
  REAL(8)                   ::  vel(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3)
  
  !--- nicht-linearer Term -----------------------------------------------------------------------------------
  REAL(8)                   ::  nl (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3)
  
  !--- Recht-Hand-Seite --------------------------------------------------------------------------------------
  REAL(8)                   ::  rhs(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3)
  
  !--- Druck -------------------------------------------------------------------------------------------------
  REAL(8)                   ::  pre(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  !--- Ausfluss-RB (Geschwindigkeitsfeld) --------------------------------------------------------------------
  ! Da die RHS für die Konzentrationsfelder nicht über die Runge-Kutta-Zwischenschritte hinweg gespeichert 
  ! werden, müssen mindestens die zugehörigen Randbedingungen gespeichert werden.
  REAL(8)                   ::  bc11(1:N2,1:N3,1:2), nlbc11(1:N2,1:N3,1:2)
  REAL(8)                   ::  bc12(0:N1,1:N3,1:2), nlbc12(0:N1,1:N3,1:2)
  REAL(8)                   ::  bc13(0:N1,1:N2,1:2), nlbc13(0:N1,1:N2,1:2)
  
  REAL(8)                   ::  bc21(0:N2,1:N3,1:2), nlbc21(0:N2,1:N3,1:2)
  REAL(8)                   ::  bc22(1:N1,1:N3,1:2), nlbc22(1:N1,1:N3,1:2)
  REAL(8)                   ::  bc23(1:N1,0:N2,1:2), nlbc23(1:N1,0:N2,1:2)
  
  REAL(8)                   ::  bc31(1:N2,0:N3,1:2), nlbc31(1:N2,0:N3,1:2)
  REAL(8)                   ::  bc32(1:N1,0:N3,1:2), nlbc32(1:N1,0:N3,1:2)
  REAL(8)                   ::  bc33(1:N1,1:N2,1:2), nlbc33(1:N1,1:N2,1:2)
  
  REAL(8)                   ::  drift1(b2L:(N2+b2U),b3L:(N3+b3U),1:2)
  REAL(8)                   ::  drift2(b1L:(N1+b1U),b3L:(N3+b3U),1:2)
  REAL(8)                   ::  drift3(b1L:(N1+b1U),b2L:(N2+b2U),1:2)
  
  !--- Residuum ----------------------------------------------------------------------------------------------
  REAL(8)                   ::  res (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  !--- Druckgradient (eine Komponente) -----------------------------------------------------------------------
  
  REAL(8)                   ::  gpre(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  !--- Gewichte für Divergenzfreiheit ------------------------------------------------------------------------
  !--- Gewichte für Divergenzfreiheit ------------------------------------------------------------------------
  REAL(8)                   ::  weight(1:N1,1:N2,1:N3)
  
  !--- Null-Raum-Vektor --------------------------------------------------------------------------------------
  REAL(8)                   ::  psi      (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL(8)                   ::  psi_vel  (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3)
  
  REAL(8)                   ::  psi_rel1 (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL(8)   , ALLOCATABLE   ::  psi_rel2 (:,:,:)
  REAL(8)   , ALLOCATABLE   ::  psi_rel3 (:,:,:)
  REAL(8)   , ALLOCATABLE   ::  psi_rel4 (:,:,:)
  REAL(8)   , ALLOCATABLE   ::  psi_rel5 (:,:,:)
  REAL(8)   , ALLOCATABLE   ::  psi_rel6 (:,:,:)
  REAL(8)   , ALLOCATABLE   ::  psi_rel7 (:,:,:)
  REAL(8)   , ALLOCATABLE   ::  psi_rel8 (:,:,:)
  REAL(8)   , ALLOCATABLE   ::  psi_rel9 (:,:,:)
  REAL(8)   , ALLOCATABLE   ::  psi_rel10(:,:,:)
  REAL(8)   , ALLOCATABLE   ::  psi_rel11(:,:,:)
  REAL(8)   , ALLOCATABLE   ::  psi_rel12(:,:,:)
  REAL(8)   , ALLOCATABLE   ::  psi_rel13(:,:,:)
  REAL(8)   , ALLOCATABLE   ::  psi_rel14(:,:,:)
  REAL(8)   , ALLOCATABLE   ::  psi_rel15(:,:,:)
  
  REAL(8)                   ::  th11(1:N2,1:N3,1:2)
  REAL(8)                   ::  th12(0:N1,1:N3,1:2)
  REAL(8)                   ::  th13(0:N1,1:N2,1:2)
  
  REAL(8)                   ::  th21(0:N2,1:N3,1:2)
  REAL(8)                   ::  th22(1:N1,1:N3,1:2)
  REAL(8)                   ::  th23(1:N1,0:N2,1:2)
  
  REAL(8)                   ::  th31(1:N2,0:N3,1:2)
  REAL(8)                   ::  th32(1:N1,0:N3,1:2)
  REAL(8)                   ::  th33(1:N1,1:N2,1:2)
  
  !--- Multigrid ---------------------------------------------------------------------------------------------
  REAL(8)                   ::  vec1C (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  REAL(8)   , ALLOCATABLE   ::  vec2A (:,:,:)
  REAL(8)   , ALLOCATABLE   ::  vec2B (:,:,:)
  REAL(8)   , ALLOCATABLE   ::  vec2C (:,:,:)
  
  REAL(8)   , ALLOCATABLE   ::  vec3A (:,:,:)
  REAL(8)   , ALLOCATABLE   ::  vec3B (:,:,:)
  REAL(8)   , ALLOCATABLE   ::  vec3C (:,:,:)
  
  REAL(8)   , ALLOCATABLE   ::  vec4A (:,:,:)
  REAL(8)   , ALLOCATABLE   ::  vec4B (:,:,:)
  REAL(8)   , ALLOCATABLE   ::  vec4C (:,:,:)
  
  REAL(8)   , ALLOCATABLE   ::  vec5A (:,:,:)
  REAL(8)   , ALLOCATABLE   ::  vec5B (:,:,:)
  REAL(8)   , ALLOCATABLE   ::  vec5C (:,:,:)
  
  REAL(8)   , ALLOCATABLE   ::  vec6A (:,:,:)
  REAL(8)   , ALLOCATABLE   ::  vec6B (:,:,:)
  REAL(8)   , ALLOCATABLE   ::  vec6C (:,:,:)
  
  REAL(8)   , ALLOCATABLE   ::  vec7A (:,:,:)
  REAL(8)   , ALLOCATABLE   ::  vec7B (:,:,:)
  REAL(8)   , ALLOCATABLE   ::  vec7C (:,:,:)
  
  REAL(8)   , ALLOCATABLE   ::  vec8A (:,:,:)
  REAL(8)   , ALLOCATABLE   ::  vec8B (:,:,:)
  REAL(8)   , ALLOCATABLE   ::  vec8C (:,:,:)
  
  REAL(8)   , ALLOCATABLE   ::  vec9A (:,:,:)
  REAL(8)   , ALLOCATABLE   ::  vec9B (:,:,:)
  REAL(8)   , ALLOCATABLE   ::  vec9C (:,:,:)
  
  REAL(8)   , ALLOCATABLE   ::  vec10A(:,:,:)
  REAL(8)   , ALLOCATABLE   ::  vec10B(:,:,:)
  REAL(8)   , ALLOCATABLE   ::  vec10C(:,:,:)
  
  REAL(8)   , ALLOCATABLE   ::  vec11A(:,:,:)
  REAL(8)   , ALLOCATABLE   ::  vec11B(:,:,:)
  REAL(8)   , ALLOCATABLE   ::  vec11C(:,:,:)
  
  REAL(8)   , ALLOCATABLE   ::  vec12A(:,:,:)
  REAL(8)   , ALLOCATABLE   ::  vec12B(:,:,:)
  REAL(8)   , ALLOCATABLE   ::  vec12C(:,:,:)
  
  REAL(8)   , ALLOCATABLE   ::  vec13A(:,:,:)
  REAL(8)   , ALLOCATABLE   ::  vec13B(:,:,:)
  REAL(8)   , ALLOCATABLE   ::  vec13C(:,:,:)
  
  REAL(8)   , ALLOCATABLE   ::  vec14A(:,:,:)
  REAL(8)   , ALLOCATABLE   ::  vec14B(:,:,:)
  REAL(8)   , ALLOCATABLE   ::  vec14C(:,:,:)
  
  REAL(8)   , ALLOCATABLE   ::  vec15A(:,:,:)
  REAL(8)   , ALLOCATABLE   ::  vec15B(:,:,:)
  
  
  !--- BiCGstab / Richardson ---------------------------------------------------------------------------------
  REAL(8)                   ::  pp(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL(8)                   ::  Ap(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL(8)                   ::  rr(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL(8)                   ::  rh(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL(8)                   ::  Ar(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL(8)                   ::  z1(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL(8)                   ::  z2(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  !--- product_div_grad --------------------------------------------------------------------------------------
  REAL(8)                   ::  dig(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  !--- Hilfsfelder (Druckiteration) --------------------------------------------------------------------------
  ! - rhs wird auch in test_moment nochmals verwendet und kann daher in outer_iteration nicht belegt werden!
  ! - wird auch fuer interpolierte Geschwindigkeiten in rhs_NS und rhs_conc verwendet.
  REAL(8)                   ::  work1(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL(8)                   ::  work2(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL(8)                   ::  work3(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))

  !--- Linienrelaxation --------------------------------------------------------------------------------------
  REAL(8)                   ::  vec1(1:N1), dia1(1:N1), SOR1(1:N1), band1(1:2,1:N1) ! TEST!!! vec1, dia muessen auch in mod_Helmholtz ausgetauscht werden!
  REAL(8)                   ::  vec2(1:N2), dia2(1:N2), SOR2(1:N2), band2(1:2,1:N2) ! TEST!!! SOR muesste man idealerweise auch in band eingliedern!
  REAL(8)                   ::  vec3(1:N3), dia3(1:N3), SOR3(1:N3), band3(1:2,1:N3) ! dia3(:) --> band3(1,:), vec3(:) --> band3(2,:), etc. ...
  
  
  !===========================================================================================================
  !=== Indizierung (Intervallgrenzen, Verschiebungen) ========================================================
  !===========================================================================================================
  !--- Block-Index -------------------------------------------------------------------------------------------
  !INTEGER                ::  iB, jB, kB ! TEST!!!
  
  !--- Indexverschiebung (Block --> global) ------------------------------------------------------------------
  INTEGER                ::  iShift, jShift, kShift
  
  !--- Domaingrösse (Periodizität-bereinigt) -----------------------------------------------------------------
  INTEGER                ::  dim1, dim2, dim33
  
  !--- Druck / Konzentrationen (inklusive Rand) --------------------------------------------------------------
  INTEGER                ::  S1p, S2p, S3p
  INTEGER                ::  N1p, N2p, N3p
  
  !--- Geschwindigkeiten (inklusive Rand) --------------------------------------------------------------------
  INTEGER                ::  S11B, S21B, S31B
  INTEGER                ::  S12B, S22B, S32B
  INTEGER                ::  S13B, S23B, S33B
  
  INTEGER                ::  N11B, N21B, N31B
  INTEGER                ::  N12B, N22B, N32B
  INTEGER                ::  N13B, N23B, N33B
  
  !--- Geschwindigkeiten (exklusive Rand) --------------------------------------------------------------------
  INTEGER                ::  S11, S21, S31
  INTEGER                ::  S12, S22, S32
  INTEGER                ::  S13, S23, S33
  
  INTEGER                ::  N11, N21, N31
  INTEGER                ::  N12, N22, N32
  INTEGER                ::  N13, N23, N33
 
  !--- grobe Gitter (Multigrid, INklusive Rand) --------------------------------------------------------------
  INTEGER                ::  S1R, S2R, S3R
  INTEGER                ::  d1R, d2R, d3R
  
  !--- grobe Gitter (Multigrid, EXklusive Rand) --------------------------------------------------------------
  INTEGER                ::  S11R, S22R, S33R
  INTEGER                ::  d11R, d22R, d33R
  
  !--- Ueberlappungskonvention der Blöcke (Multigrid, siehe mod_setup) ---------------------------------------
  INTEGER, PARAMETER     ::  ls1 = -1
  INTEGER, PARAMETER     ::  ls2 = -1
  INTEGER, PARAMETER     ::  ls3 = -1
  
  !--- Austauschrichtung (Multigrid) -------------------------------------------------------------------------
  ! ex = -1: unten <--  oben
  ! ex =  0: unten <--> oben
  ! ex =  1: unten  --> oben
  INTEGER                ::  ex1, ex2, ex3
  
  
  !===========================================================================================================
  !=== Randbedingungen =======================================================================================
  !===========================================================================================================
  !                              _
  !    Symmetrie-RB:   BC = -2    |
  !    periodische RB: BC = -1    |- symmetrische, zentrale Stencils
  !    Nachbar-Block:  BC =  0   _|
  !    Dirichlet-RB:   BC =  1    |
  !    Neumann-RB:     BC =  2    |- schiefe, nicht-zentrale Stencils
  !    Robin-RB:       BC =  3   _|
  !
  !--- global ------------------------------------------------------------------------------------------------
  LOGICAL                ::  outlet   (1:3,1:2,1:3)
  
  INTEGER                ::  BC_1L_global, BC_1U_global
  INTEGER                ::  BC_2L_global, BC_2U_global
  INTEGER                ::  BC_3L_global, BC_3U_global
  
  !--- lokal (Block) -----------------------------------------------------------------------------------------
  INTEGER                ::  BC_1L, BC_1U
  INTEGER                ::  BC_2L, BC_2U
  INTEGER                ::  BC_3L, BC_3U

  !--- field properties --------------------------------------------------------------------------------------
  INTEGER                ::  n_gather(1:3,1:n_grids_max)
  INTEGER                ::  NN (1:3,1:n_grids_max)
  INTEGER                ::  NB (1:3,1:n_grids_max)
  INTEGER                ::  iB (1:3,1:n_grids_max)
  INTEGER                ::  SNF(1:2,1:3,1:n_grids_max)
  INTEGER                ::  SNB(1:2,1:3,1:n_grids_max)
  INTEGER                ::  BC (1:2,1:3,1:n_grids_max)
  INTEGER                ::  ngb(1:2,1:3,1:n_grids_max)
  INTEGER                ::  comm1(1:n_grids_max), comm2(1:n_grids_max)
  INTEGER                ::  rankc2(1:n_grids_max)
  LOGICAL                ::  participate_yes(1:n_grids_max)
  INTEGER                ::  recvR(    1:NB1*NB2*NB3,1:n_grids_max)
  INTEGER                ::  recvI(    1:NB1*NB2*NB3,1:n_grids_max)
  INTEGER                ::  dispR(    1:NB1*NB2*NB3,1:n_grids_max)
  INTEGER                ::  dispI(    1:NB1*NB2*NB3,1:n_grids_max)
  INTEGER                ::  offsR(1:3,1:NB1*NB2*NB3,1:n_grids_max)
  INTEGER                ::  offsI(1:3,1:NB1*NB2*NB3,1:n_grids_max)
  INTEGER                ::  sizsR(1:3,1:NB1*NB2*NB3,1:n_grids_max)
  INTEGER                ::  sizsI(1:3,1:NB1*NB2*NB3,1:n_grids_max)
  
  
  !===========================================================================================================
  !=== physikalische Parameter ===============================================================================
  !===========================================================================================================
  REAL(8)                   ::  L1, L2, L3
  REAL(8)                   ::  Re
  
  
  !===========================================================================================================
  !=== numerische Parameter ==================================================================================
  !===========================================================================================================
  !--- allgemein ---------------------------------------------------------------------------------------------
  REAL(8)                   ::  CFL
  REAL(8)                   ::  time, dtime, subtime, time_start, time_end, dtime_max, dtime0, dtime_old
  INTEGER                ::  timestep, timestep_old, substep, n_timesteps
  LOGICAL                ::  mapping_yes, upwind_yes
  LOGICAL                ::  Euler_yes, Stokes_yes, twostep_yes
  LOGICAL, PARAMETER     ::  filter_BC_yes = .TRUE. ! TEST!!!
  INTEGER                ::  timeint_mode, forcing_mode
  INTEGER                ::  bulkflow_dir
  INTEGER                ::  n_lp_vel           , n_hp_vel
  REAL(8)                   ::  chi_vel
  REAL(8)                   ::  vel_bulk ! TEST!!!
  
  
  !--- Runge-Kutta-Koeffizienten -----------------------------------------------------------------------------
  REAL(8)   , PARAMETER     ::  aRK(1:3) = (/8./15.,  5./12., 3./ 4./)
  REAL(8)   , PARAMETER     ::  bRK(1:3) = (/  0.  ,-17./60.,-5./12./)
  INTEGER, PARAMETER     ::  RK_steps = 3
  
  !--- look-up table fuer Stabilitaetsgebiet der Zeitintegration (angle = pi/2,pi) ---------------------------
  REAL(8)   , PARAMETER     ::  stabilitylimit(0:40) = (/1.732050813, 1.943689093, 2.089210537, 2.201001743,  &
                                      &               2.290031261, 2.361554127, 2.418567407, 2.462989697,  &
                                      &               2.496169963, 2.519146008, 2.532795254, 2.537935070,  &
                                      &               2.535397854, 2.526091466, 2.511046932, 2.491448818,  &
                                      &               2.468639045, 2.444084180, 2.419302172, 2.395757241,  &
                                      &               2.374745783, 2.357302135, 2.344145473, 2.335672458,  &
                                      &               2.331985072, 2.332936948, 2.338183901, 2.347230689,  &
                                      &               2.359471631, 2.374225928, 2.390769340, 2.408363261,  &
                                      &               2.426281290, 2.443832601, 2.460381269, 2.475360992,  &
                                      &               2.488285197, 2.498753090, 2.506452564, 2.511161051,  &
                                      &               2.512745327 /)
  INTEGER, PARAMETER     ::  n_stab = SIZE(stabilitylimit)
  
  !--- Helmholtz-Vorfaktoren ---------------------------------------------------------------------------------
  REAL(8),target                   ::  thetaL, multL
  
  !--- zeitliche Steuerung -----------------------------------------------------------------------------------
  INTEGER                ::  Int_dtime, Int_lev_pre
  
  INTEGER                ::  stride_large(1:3), stride_med(1:3), stride_small(1:3)
  LOGICAL                ::  write_large, write_med, write_small
  REAL(8)                   ::  time_out_scal, dtime_out_scal
  REAL(8)                   ::  time_out_vect, dtime_out_vect
  
  LOGICAL                ::  write_out_scal, write_out_vect  
  LOGICAL                ::  new_dtime, finish_yes
    
  INTEGER                ::  write_count
  INTEGER                ::  restart
  CHARACTER(LEN=3)       ::  restart_char
  
  INTEGER                ::  n_conc_old


  !===========================================================================================================
  !=== weitere Steuerungsoptionen ============================================================================
  !===========================================================================================================
  INTEGER                ::  task
  LOGICAL                ::  read_nullspace_yes
  LOGICAL                ::  nullspace_yes, nullspace_coarse_yes
  LOGICAL                ::  nullspace_ortho_yes
  
  LOGICAL                ::  write_stout_yes
  LOGICAL                ::  log_iteration_yes
  LOGICAL                ::  write_restart_yes
  LOGICAL                ::  write_lambda2_yes
  LOGICAL                ::  write_test_yes
  
  !--- globale Laufindizes -----------------------------------------------------------------------------------
  INTEGER                ::  direction
  
  !--- explizite Behandlung von Ecken bei Dirichlet-Randbedingungen ------------------------------------------
  ! (Hintergrund: Der Druck ist an diesen Orten unbestimmt, so dass er dort künstlich zu Null gesetzt wird.)
  LOGICAL, PARAMETER     ::  corner_yes = .TRUE.
  
  !--- Systemzeit --------------------------------------------------------------------------------------------
  INTEGER                ::  elatime, ctime(1:8)
  INTEGER                ::  day, hour, minu, sec, msec


  !===========================================================================================================
  !=== Iterationsparameter ===================================================================================
  !===========================================================================================================
  !--- Abbruchkriterium / Absolute Genauigkeit der Geschwindigkeiten -----------------------------------------
  REAL(8)                   ::  epsU, epsU0
  
  !--- Glaetter ----------------------------------------------------------------------------------------------
  LOGICAL                ::  Jacobi_yes
  
  !--- max. Anzahl Iterationen -------------------------------------------------------------------------------
  INTEGER                ::  n_it_outer
  INTEGER                ::  n_it_Poisson
  INTEGER                ::  n_it_Helmh_vel
  
  !--- erwartete Konvergenzrate (äussere Iteration) ----------------------------------------------------------
  REAL(8)                   ::  precRatio0 (1:RK_steps)
  REAL(8)                   ::  precOffset0(1:RK_steps)
  REAL(8), ALLOCATABLE      ::  precRatio  (:,:)
  REAL(8), ALLOCATABLE      ::  precOffset (:,:)
  
  !--- Null-Initialisierung (äussere Iteration) --------------------------------------------------------------
  LOGICAL                ::  init_pre(1:RK_steps), init_vel(1:RK_steps)!, init_conc(1:RK_steps)
  
  !--- Vorkonditionierung (Multigrid) ------------------------------------------------------------------------
  INTEGER                ::  precond_outer
  INTEGER                ::  precond_Poisson
  INTEGER                ::  precond_Helmh_vel
  
  !--- Anzahl Glättungen pro Gitterlevel (Multigrid) ---------------------------------------------------------
  INTEGER                ::  n_relax_down, n_relax_up, n_relax_bottom
  
  !--- implizite Richtungen bei Linienrelaxation (Multigrid) -------------------------------------------------
  INTEGER                ::  impl_dir(1:3)
  
  !--- Anzahl Glättungen pro Gitterlevel (Multigrid) ---------------------------------------------------------
  LOGICAL                ::  weighting_yes
  
  !===========================================================================================================
  !=== Iterationsstatistik ===================================================================================
  !===========================================================================================================
  REAL(8)                   ::  dtime_average
  REAL(8)                   ::  max_div_init(1:2)
  INTEGER                ::  number_poisson
  
  !--- Zähler ------------------------------------------------------------------------------------------------
  INTEGER                ::  countO(1:RK_steps)
  INTEGER                ::  countP(1:RK_steps,1:2)
  INTEGER                ::  countH(1:RK_steps,1:3)
  
  !--- Konvergenzrate ----------------------------------------------------------------------------------------
  REAL(8)                   ::  ratioO(1:RK_steps)
  REAL(8)                   ::  ratioH(1:RK_steps,1:3)
  REAL(8)                   ::  ratioP(1:RK_steps,1:2)
  
  
  !===========================================================================================================
  !=== MPI ===================================================================================================
  !===========================================================================================================
  !--- Kommunikatoren ----------------------------------------------------------------------------------------
  INTEGER                ::  COMM_CART
  
  INTEGER                ::  COMM_SLICE1, COMM_BAR1
  INTEGER                ::  COMM_SLICE2, COMM_BAR2
  INTEGER                ::  COMM_SLICE3, COMM_BAR3
  
  !--- Dimension und Position der Blöcke innerhalb der Kommunikatoren (Gitter-Indizes) -----------------------
  ! (für MPI_GATHERv, MPI_ALLGATHERv, vgl. iShift, jShift, kShift)
  INTEGER                ::  bar1_size(1:NB1), bar1_offset(1:NB1)
  INTEGER                ::  bar2_size(1:NB2), bar2_offset(1:NB2)
  INTEGER                ::  bar3_size(1:NB3), bar3_offset(1:NB3)
  
  !--- Ränge der Prozesse ------------------------------------------------------------------------------------
  INTEGER                ::  rank
  INTEGER                ::  rank_bar1, rank_slice1
  INTEGER                ::  rank_bar2, rank_slice2
  INTEGER                ::  rank_bar3, rank_slice3
  
  !--- Ränge der Nachbarprozesse (in kartesischem Gitter) ----------------------------------------------------
  INTEGER                ::  rank1L, rank1U
  INTEGER                ::  rank2L, rank2U
  INTEGER                ::  rank3L, rank3U
  
  !--- Error-Handle ------------------------------------------------------------------------------------------
  INTEGER                ::  merror
  
  !--- Request-Handles ---------------------------------------------------------------------------------------
  ! (müssen offenbar global angelegt werden) 
  INTEGER                ::  req1L, req1U
  INTEGER                ::  req2L, req2U
  INTEGER                ::  req3L, req3U
  
  !===========================================================================================================
  !=== HDF5 ==================================================================================================
  !===========================================================================================================
  INTEGER                ::  herror
  
  !*** bbecsek 2015: IMMERSED BOUNDARY PARAMETERS ************************************************************
  !REAL(8)                   ::  delta_s                 !< boundary mesh width
  REAL(8)                   ::  reach                   !< reach of interpolation reach
  INTEGER                ::  M_bound                 !< number of grid points on boundary
  INTEGER                ::  M_elems                 !< number of finite elements
  REAL(8)   , ALLOCATABLE   ::  yb(:,:)                 !< global Lagrangian boundary coordinates (allocated in init_ib)
  REAL(8)   , ALLOCATABLE   ::  xb(:,:)                 !< reference Lagrangian boundary coordinates (if needed) 
                                                     !! (allocated in init_ib)
  REAL(8)   , ALLOCATABLE   ::  db(:,:)                 !< Lagrangian displacement field
  REAL(8)   , ALLOCATABLE   ::  ub(:,:)                 !< boundary velocities
  REAL(8)   , ALLOCATABLE   ::  fb(:,:)                 !< Lagrangian force (allocated in init_ib)
  REAL(8)   , ALLOCATABLE   ::  fd(:,:,:,:)             !< force density
  REAL(8)   , ALLOCATABLE   ::  container(:,:)          !< temporal container for the time-lagged RK vel-term (may be omitted 
                                                     !! later on by either performing the respective additions before and after 
                                                     !! the boundary update or by using the 'buffer' containers that already 
                                                     !! store the velocity fields
  !REAL(8)                   ::  kappa                   !< proportionality constant used for the Lagrangian force computation
  INTEGER                ::  ddf_type                !< flag for choosing the spreading and interpolation kernel function
                                                     !! (ddf: discrete delta function)
  LOGICAL                ::  IB_on                   !< specifies whether to immerse a boundary 
  REAL(8)                   ::  mu_blood                !< dynamic viscosity of blood
  REAL(8)                   ::  rho_blood               !< bulk density of blood

  REAL(8)                   ::  E_mod                   !< Young's modulus of elasticity (linear) 
  REAL(8)                   ::  nu_poiss                !< Poisson's ratio  
 
  REAL(8)                   ::  rho_solid               !< density of solid
  REAL(8)   , ALLOCATABLE   ::  strain(:,:)             !< element strains
  REAL(8)   , ALLOCATABLE   ::  stress(:,:)             !< element stress
  REAL(8)   , ALLOCATABLE   ::  k_stiff(:,:)            !< global stiffness matrix
  REAL(8)   , ALLOCATABLE   ::  m_mass(:,:)             !< global mass matrix
  REAL(8)   , ALLOCATABLE   ::  c_damp(:,:)             !< global damping matrix (C=alpha*M + beta*K)
  REAL(8)                   ::  C(3,3)                  !< matrix linking strain to stress


  INTEGER, ALLOCATABLE   ::  elems(:,:)              !< Finite elements

  REAL(8)                   ::  L_ref                   !< reference length scale
  REAL(8)                   ::  U_ref                   !< reference velocity

  INTEGER                ::  req_bcast(3), req_red(3)      !< request handles

  INTEGER                ::  excl_p(3)               !< points to exclude

  REAL(8)                   ::  delta_b1,delta_b2,delta_b3                 !< boundary discretizations

  INTEGER                ::  total_ranks             !< total number of ranks

  INTEGER                ::  M_ebcs                  !< number of prescribed essential boundary conditions

  REAL(8)   , ALLOCATABLE   ::  ebcs(:,:)               !< container for essential boundary conditions

  REAL(8)   , ALLOCATABLE   ::  node_vol(:)             !< volume (area) associated with node i

  INTEGER                :: COMM_LOCAL               !< communicator for local group
  INTEGER                :: COMM_INTER               !< intercommunicator solid-fluid

  REAL(8)                   ::  max_ev_freq             !< highest ev (frequency) of structure

  LOGICAL                ::  write_force_yes         !< whether to output force density or not

  LOGICAL                ::  write_xdmf_yes          !< whether to output .xmf files
  LOGICAL                ::  scale_output_yes        !< whether to scale all values with U_ref, L_ref (dimensionalize)

  !--- Oscillator Aorta Model ---
  LOGICAL                ::  fem_yes
  REAL(8)                   ::  k_x, k_y, k_z
  REAL(8)                   ::  aorta_start, aorta_end, aorta_radius
#endif
  
  
END MODULE mod_vars
