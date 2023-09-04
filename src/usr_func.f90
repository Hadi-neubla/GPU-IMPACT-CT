!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!* modified by Barna Becsek, ARTORG CVE, University of Bern (barna.becsek@artorg.unibe.ch)                   *
!* October 2014                                                                                              *
!* GPU version by Hadi Zolfaghari , ARTORG CVE, then DAMTP, Cambridge University (hz382@cam.ac.uk)           *
!* Oct 2015 - Sep 2023                                                                                       *
!*************************************************************************************************************

!> module containing functions specified by the user
MODULE usr_func
  
  
  USE mod_dims
  USE mod_vars
  USE mod_lib
  USE mod_diff
  USE mod_exchange
  USE usr_vars
  !USE mpi  
  
  PRIVATE

  PUBLIC read_vtk, read_vtk_3D, get_process_CT_geometry, read_velocity_waveform, convert_waveform, get_pulse_number_and_length !hadi zolfaghari

  PUBLIC get_pressure_and_flow_rate_WKM, apply_CT_forcing, apply_pressure_BC_CT
  PUBLIC global_indicator, sphere_global, sphere_local    !hadi zolfaghari
  
  PUBLIC smooth_step !bbecsek
  PUBLIC linear_step !bbecsek
  PUBLIC fringe_coeff !bbecsek
  PUBLIC poiseuille_parabola !bbecsek
  PUBLIC interface, erf, atanh
  PUBLIC coord_tanh
  PUBLIC init_particles
  PUBLIC init_hdf5, init_mpi, finl_hdf5, finl_mpi, pass_mpi_comm, print_fcomm_size, mpi_bcast_fort !bbecsek
  PUBLIC windkessel_integration_step, pdot3EWK, pdot4EWK
  PUBLIC flow_and_rate, flow_and_rate2D, flow_and_rate3D, fitted_flow 
  PUBLIC apply_fringe_forcing, apply_windkessel_loading
  PUBLIC check_node_ids
  PUBLIC block_id, block_cart
  PUBLIC total_n_local_tet_elements, total_n_global_tet_elements
  PUBLIC total_n_local_tri_elements, total_n_global_tri_elements
  PUBLIC local_to_global_node_id_loc_con_allperiodic2
  PUBLIC global_to_local_node_id_loc_con_allperiodic2
  PUBLIC global_id2cart_loc_con, global_cart2id_loc_con
  PUBLIC n_local_nodes, n_global_nodes
  PUBLIC local_id2cart_loc_con, local_cart2id_loc_con
  PUBLIC local_tet_element_nodes
  PUBLIC local_tri_element_nodes
  PUBLIC interpolate_force_pre_vel
  PUBLIC local_to_global_tet_elem_id
  PUBLIC local_to_global_tri_elem_id
  PUBLIC residual2volume_force

  INCLUDE 'mpif.h'
  
  CONTAINS
  
!pgi$g unroll = n:8
!!pgi$r unroll = n:8
!!pgi$l unroll = n:8
  
  !=========================================================================================================!
  !--- user-specific subroutines and functions ---                                                          !
  !=========================================================================================================!
 
  SUBROUTINE read_vtk(x_offset, y_offset, z_offset, x, grid_set)
  IMPLICIT NONE
  CHARACTER(len=100) :: dummy
  CHARACTER(len=6) :: dataformat
  CHARACTER(len=20) :: dataset
  INTEGER           :: i, j
  INTEGER :: status, np, k, ln
  INTEGER, DIMENSION(3) :: dimensions
  REAL(8), DIMENSION(2) :: origin, dx, spacing
  integer, intent(in):: x_offset, y_offset, z_offset
  REAL(8), intent(out), allocatable :: x(:)
  REAL(8), intent(out) :: grid_set(1:M1,1:M2,1:M3) 
  REAL(8)               :: d1, d2, d3, d4, d5 


  OPEN(UNIT=8, FILE='CT1.vtk', IOSTAT=status, STATUS='old')
  
!  First two lines contain the version number and a title
  READ(8, '(A)') dummy
  READ(8, '(A)') dummy
  
!  The third line contains ASCII or BINARY
  READ(8, '(A)') dataformat
  
!  Fourth line contains the type of dataset
  READ(8, '(A7,1X,A)') dummy, dataset
  print*, dataset, dataformat
!  Dimensions
  READ(8, *) dummy, dimensions
  print*, dimensions
!  SPACING
  READ(8, *) dummy, dx
  print*, dummy, dx

!  Origin
  READ(8, *) dummy, origin
  print*, dummy, origin


! Point data
  READ(8, *) dummy, np
  print*, dummy, np
  
! Scalars
  READ(8, *) dummy,dataset
  print*, dummy, dataset
! LOOKUP Table
  READ(8, *) dummy
  print*, dummy


  allocate(x(np))

  READ(8, *) x
!  print*, x(200*dimensions(1)+1: 201*dimensions(1))


!  Finally read data
! DO ln = 1, 1
 !   READ(8, *) x
! END DO
  
  CLOSE(UNIT=8, IOSTAT=status)

 !-- fill up the ct_grid with CT data 
 !-- beyond offset is the fringe region indicated by 2.0

 grid_set  = 2.0
 DO i=1, dimensions(1)
    Do j=1, dimensions(2)
       grid_set(j+y_offset,i+x_offset,1) = x((j-1)*dimensions(1)+i) 
    ENDDO
 ENDDO

 grid_set(1:M1,:,:) = grid_set(M1:1:-1,:,:)

 !--- get the index bounds of aortic inlet

outer: DO j=1, 200
    Do i=10, 130
       IF (grid_set(i,j,1) .eq. 1.0) then
           j_start_aorta_inlet = j
           i_start_aorta_inlet = i
           if (rank == 0) print*, i, j, 'check'
           EXIT outer
       ENDIF
    ENDDO
 ENDDO outer


 DO j=j_start_aorta_inlet, 200
    IF (grid_set(i_start_aorta_inlet,j,1) .eq. 0.0) then
        j_end_aorta_inlet = j-1
        exit
    ENDIF
 ENDDO

!if (rank == 0) print*, 'length, start, end', (y2p(j_end_aorta_inlet) - y2p(j_start_aorta_inlet)), j_start_aorta_inlet, j_end_aorta_inlet 

 !--- update the Reynolds number according to inlet geometry
      !--- star values are set by user either directly (velocity) or indirectly (length): latter by setting the box dimensionless size

 aortic_diameter_star = (y2p(j_end_aorta_inlet) - y2p(j_start_aorta_inlet))
 maximum_inflow_velocity_star = 1.0


!--- mask the inlet fringe zone

 DO j=1, 200
    Do i=10, 130
       IF (grid_set(i,j,1) .eq. 1.0) then
           grid_set(i,j,1) = 5.0
       ENDIF
    ENDDO
 ENDDO

 END SUBROUTINE read_vtk

  SUBROUTINE read_vtk_3D(x_offset, y_offset, z_offset, x, grid_set)
  IMPLICIT NONE
  CHARACTER(len=100)             :: dummy
  CHARACTER(len=6)               :: dataformat
  CHARACTER(len=20)              :: dataset
  INTEGER                        :: i, j
  INTEGER                        :: status, np, k, ln
  INTEGER, DIMENSION(3)          :: dimensions
  REAL(8), DIMENSION(2)             :: origin, dx, spacing
  INTEGER, INTENT(IN)            :: x_offset, y_offset, z_offset
  REAL(8), INTENT(OUT), ALLOCATABLE :: x(:)
  REAL(8), INTENT(OUT)              :: grid_set(1:M1,1:M2,1:M3) 
  REAL(8)                           :: d1, d2, d3, d4, d5 
  REAL(8)                           :: x_mean_centre_inlet, y_mean_centre_inlet, R_inlet
!  INTEGER, ALLOCATABLE           :: x_inlet_index(:), y_inlet_index(:)
!  INTEGER                        :: counter_inlet, counter_outlet
  REAL(8)                           :: x_mean_centre_outlet, y_mean_centre_outlet
!  INTEGER, ALLOCATABLE           :: x_outlet_index(:), y_outlet_index(:)



! --- ideal case resolution 256
!  OPEN(UNIT=8, FILE='CT_3D_ideal.vtk', IOSTAT=status, STATUS='old')
! --- ideal case resolution 512
!  OPEN(UNIT=8, FILE='CT_3D_ideal_highRes.vtk', IOSTAT=status, STATUS='old')
! --- ideal case with upstream brachiocephalic branch - resolution 512
  OPEN(UNIT=8, FILE='CT_3D_ideal_highRes_revision.vtk', IOSTAT=status, STATUS='old')

! --- TAVI 000
!  OPEN(UNIT=8, FILE='CT_3D_TAVI000.vtk', IOSTAT=status, STATUS='old')
! --- low resolution 
!  OPEN(UNIT=8, FILE='CT_3D_lowRes.vtk', IOSTAT=status, STATUS='old')

! --- TAVI 001
!  OPEN(UNIT=8, FILE='CT_3D_TAVI001_3.vtk', IOSTAT=status, STATUS='old')
! --- TAVI 002
 ! OPEN(UNIT=8, FILE='CT_3D_TAVI002.vtk', IOSTAT=status, STATUS='old')


!  First two lines contain the version number and a title
  READ(8, '(A)') dummy
  READ(8, '(A)') dummy
  
!  The third line contains ASCII or BINARY
  READ(8, '(A)') dataformat
  
!  Fourth line contains the type of dataset
  READ(8, '(A7,1X,A)') dummy, dataset
  print*, dataset, dataformat
!  Dimensions
  READ(8, *) dummy, dimensions
  print*, dimensions
!  SPACING
  READ(8, *) dummy, dx
  print*, dummy, dx

!  Origin
  READ(8, *) dummy, origin
  print*, dummy, origin


! Point data
  READ(8, *) dummy, np
  print*, dummy, np
  
! Scalars
  READ(8, *) dummy,dataset
  print*, dummy, dataset
! LOOKUP Table
  READ(8, *) dummy
  print*, dummy


  allocate(x(np))

  READ(8, *) x
!  print*, x(200*dimensions(1)+1: 201*dimensions(1))


!  Finally read data
! DO ln = 1, 1
 !   READ(8, *) x
! END DO
  
  CLOSE(UNIT=8, IOSTAT=status)

 !-- fill up the ct_grid with CT data 
 !-- beyond offset is the fringe region indicated by 2.0

 !-- levelset key:
 !---- 2.0 : outside the CT box
 !---- 1.0 : the lumen (except the fringe)
 !---- 5.0 : the jet part of the fringe for inlet
 !---- 7.0 : the peripheral part of the fringe inlet which does damping 
 !---- 9.0 : the flow part of the fringe for outlet 

 grid_set  = 2.0
 DO i=1, dimensions(1)
    Do j=1, dimensions(2)
       Do k=1, dimensions(3)
          grid_set(j+y_offset,i+x_offset,k+z_offset) = x((k-1)*dimensions(2)*dimensions(1)+(j-1)*dimensions(1)+i) 
       ENDDO
    ENDDO
 ENDDO

 
! grid_set(1:M1,:,:) = grid_set(M1:1:-1,:,:)


 !--- mask aortic inlet

!--- no inlet profile specification (e.g. TAVI or given profile, thus flat profile)
! outer: DO j=1, floor(M2/2.)
!          Do i=1, M1
!             DO k =1, floor(M3/12.)
!                IF (grid_set(i,j,k) .eq. 1.0) grid_set(i,j,k) = 5.0
!               ENDDO
!          ENDDO
!       ENDDO outer


!-- INLET LABELING
!--- find the centre of aortic jet
        x_mean_centre_inlet = 0.
        y_mean_centre_inlet = 0.
        counter_inlet = 0
        k_start_index = z_offset + 1
outer1: DO j=1, M2 ! was floor(M2/2.)
          Do i=1, floor(M1/2.)
                IF ((grid_set(i,j,k_start_index) .ge. 0.5) .and. (grid_set(i,j,k_start_index) .lt. 2.0)) THEN
                   x_mean_centre_inlet = x_mean_centre_inlet + y1p(i)
                   y_mean_centre_inlet = y_mean_centre_inlet + y2p(j)
                   counter_inlet = counter_inlet + 1
                ENDIF       
          ENDDO
       ENDDO outer1

       allocate(x_inlet_index(counter_inlet))
       allocate(y_inlet_index(counter_inlet)) 

       !inlet_voxels_count = counter_inlet ! saving to a global variable


       counter_inlet = 0

       DO j=1, M2 ! was floor(M2/2.)
          Do i=1, floor(M1/2.)
                IF ((grid_set(i,j,k_start_index) .ge. 0.5) .and. (grid_set(i,j,k_start_index) .lt. 2.0)) THEN
                   counter_inlet = counter_inlet + 1
                   x_inlet_index(counter_inlet) = i
                   y_inlet_index(counter_inlet) = j
                ENDIF       
          ENDDO
       ENDDO 



       x_mean_centre_inlet = x_mean_centre_inlet/counter_inlet
       y_mean_centre_inlet = y_mean_centre_inlet/counter_inlet


       IF (RANK == 0) print*, x_mean_centre_inlet, y_mean_centre_inlet, counter_inlet, 'centre coords inlet'




!-- OUTLET LABELING
!--- find the centre of aortic jet
        x_mean_centre_outlet = 0.
        y_mean_centre_outlet = 0.
        counter_outlet = 0
        k_start_index = z_offset + 1
        DO j= 1,M2 ! 1 was floor(M2/2.)
          Do i=floor(M1/2.),M1
                IF ((grid_set(i,j,k_start_index) .ge. 0.5) .and. (grid_set(i,j,k_start_index) .lt. 2.0)) THEN
                   x_mean_centre_outlet = x_mean_centre_outlet + y1p(i)
                   y_mean_centre_outlet = y_mean_centre_outlet + y2p(j)
                   counter_outlet = counter_outlet + 1
                ENDIF       
          ENDDO
       ENDDO 

       allocate(x_outlet_index(counter_outlet))
       allocate(y_outlet_index(counter_outlet)) 

       counter_outlet = 0

       DO j=1, M2 ! 1 was floor(M2/2.)
          Do i=floor(M1/2.),M1
                IF ((grid_set(i,j,k_start_index) .ge. 0.5) .and. (grid_set(i,j,k_start_index) .lt. 2.0)) THEN
                   counter_outlet = counter_outlet + 1
                   x_outlet_index(counter_outlet) = i
                   y_outlet_index(counter_outlet) = j               
                ENDIF       
          ENDDO
       ENDDO 


       x_mean_centre_outlet = x_mean_centre_outlet/counter_outlet
       y_mean_centre_outlet = y_mean_centre_outlet/counter_outlet


       IF (RANK == 0) print*, x_mean_centre_outlet, y_mean_centre_outlet, counter_outlet, 'centre coords outlet'

!-- Cranial OUTLETS LABELING: outlet2
!--- find the centre of aortic jet
        counter_outlet_2 = 0
        k_start_index_2 = z_offset + dimensions(3) ! this is where the upper domain ends (buffer for outlets starts)
        DO j= 1,M2 ! 1 was floor(M2/2.)
           Do i=1,M1
                IF ((grid_set(i,j,k_start_index_2) .ge. 0.5) .and. (grid_set(i,j,k_start_index_2) .lt. 2.0)) THEN
                   counter_outlet_2 = counter_outlet_2 + 1
                ENDIF       
          ENDDO
       ENDDO 

       allocate(x_outlet_2_index(counter_outlet_2))
       allocate(y_outlet_2_index(counter_outlet_2)) 

       counter_outlet_2 = 0

       DO j=1, M2 ! 1 was floor(M2/2.)
          Do i=1,M1
                IF ((grid_set(i,j,k_start_index_2) .ge. 0.5) .and. (grid_set(i,j,k_start_index_2) .lt. 2.0)) THEN
                   counter_outlet_2 = counter_outlet_2 + 1
                   x_outlet_2_index(counter_outlet_2) = i
                   y_outlet_2_index(counter_outlet_2) = j               
                ENDIF       
          ENDDO
       ENDDO 



       IF (RANK == 0) print*, counter_outlet_2, 'number of cells within outflow 2 boundary'




!--- consitricted inflow  (centered jets)
!     - TAVI 00 xc=0.32,0.35 (TAVI_000_4)  yc=0.23, rc=0.07, 05 (TAVI_000_2, TAVI_000_5), 0.025(TAVI_000_3, TAVI_000_4)
!     - based on the reference values defined on Jan 11 2022
!     - TAVI 00 xc=0.32,0.35 (TAVI_000_4)  yc=0.23, rc=0.07, 05 (TAVI_000_2, TAVI_000_5), 0.025(TAVI_000_3, TAVI_000_4)
!     - TAVI 02    0.15     0.26     0.06
! outer: DO j=1, floor(M2/2.)
!          Do i=1, floor(M1/2.)
!             DO k =1, floor(M3/12.)
!                IF ((grid_set(i,j,k) .ge. 0.5) .and. (grid_set(i,j,k) .lt. 2.)) THEN
!                      grid_set(i,j,k) = 7.0  ! TAVI_000_1
!                   IF (sqrt((y1p(i)-0.32)**2+(y2p(j)-0.23)**2) .le. 0.025) THEN
!                      grid_set(i,j,k) = 5.0
!                   ENDIF
!                ENDIF       
!             ENDDO
!          ENDDO
!       ENDDO outer
!



 !-- marking up the inlet fringe zone
     !-- getting the radius R_inlet for the incoming jet when Q_max = 0.0004 m3/s
     R_inlet = sqrt(0.0004/(4. * atan(1.) * reference_velocity))/reference_length

     IF (rank == 0) print*, R_inlet, 'the incoming jet radius'
     inlet_voxels_count = 0
     IF (idealized_aorta == 0) THEN
        DO i=1, counter_inlet
           DO k =-10, -1
              grid_set(x_inlet_index(i), y_inlet_index(i), k_start_index+k) = 7.0
              IF (sqrt((y1p(x_inlet_index(i))-x_mean_centre_inlet)**2+(y2p(y_inlet_index(i))-y_mean_centre_inlet)**2) .le. R_inlet) THEN
                 grid_set(x_inlet_index(i),y_inlet_index(i),k_start_index+k) = 5.0
                 if (k == -1) inlet_voxels_count = inlet_voxels_count + 1
              ENDIF
          ENDDO
        ENDDO
     ELSE
         DO i=1, counter_inlet
           DO k =-10, -1
              grid_set(x_inlet_index(i), y_inlet_index(i), k_start_index+k) = 7.0
              IF (sqrt((y1p(x_inlet_index(i))-x_mean_centre_inlet)**2+(y2p(y_inlet_index(i))-y_mean_centre_inlet)**2) .le. R_inlet) THEN
                 grid_set(x_inlet_index(i),y_inlet_index(i),k_start_index+k) = 5.0
                 if (k == -1) inlet_voxels_count = inlet_voxels_count + 1
              ENDIF
          ENDDO
        ENDDO

     ENDIF

     IF (rank == 0) print*, 'number of cells within the incoming jet ',  inlet_voxels_count

 !-- marking up the outlet2 fringe zone

     DO i=1, counter_outlet
        DO k =-10, -1
           grid_set(x_outlet_index(i), y_outlet_index(i), k_start_index+k) = 9.0
!           IF (sqrt((x_outlet_index(i)-x_mean_centre_outlet)**2+(y_outlet_index(i)-y_mean_centre_outlet)**2) .le. 0.025) THEN
!              grid_set(i,j,k_start_index+k) = 9.0
!           ENDIF
       ENDDO
     ENDDO

 !-- marking up the outlet 2 fringe zone

     DO i=1, counter_outlet_2
        DO k =1, 10
           grid_set(x_outlet_2_index(i), y_outlet_2_index(i), k_start_index_2+k) = 11.0
!           IF (sqrt((x_outlet_index(i)-x_mean_centre_outlet)**2+(y_outlet_index(i)-y_mean_centre_outlet)**2) .le. 0.025) THEN
!              grid_set(i,j,k_start_index+k) = 9.0
!           ENDIF
       ENDDO
     ENDDO







 !--- get the index bounds of aortic inlet

!**outer: DO j=1, 200
!**          Do i=10, 130
!**             DO k =1, 100
!**                IF (grid_set(i,j,1) .eq. 1.0) then
!**                   j_start_aorta_inlet = j
!**                   i_start_aorta_inlet = i
!**                   if (rank == 0) print*, i, j, 'check'
!**                   EXIT outer
!**               ENDDO
!**             ENDIF
!**          ENDDO
!**       ENDDO outer
!**
!**
!** DO j=j_start_aorta_inlet, 200
!**    IF (grid_set(i_start_aorta_inlet,j,1) .eq. 0.0) then
!**        j_end_aorta_inlet = j-1
!**        exit
!**    ENDIF
!** ENDDO

!if (rank == 0) print*, 'length, start, end', (y2p(j_end_aorta_inlet) - y2p(j_start_aorta_inlet)), j_start_aorta_inlet, j_end_aorta_inlet 

 !--- update the Reynolds number according to inlet geometry
      !--- star values are set by user either directly (velocity) or indirectly (length): latter by setting the box dimensionless size

!** aortic_diameter_star = (y2p(j_end_aorta_inlet) - y2p(j_start_aorta_inlet))
 ! aortic_diameter_star = 0.15  ! APPROXIMATE VALUE, should be made model-based 

! upsampling_factor = 2. ! TAVI003
! upsampling_factor = 2. ! TAVI002
! upsampling_factor = 2. ! TAVI001
!upsampling_factor = 1. ! TAVI000
   
! CT_pixel_width_star = (1./M1) * upsampling_factor  
! maximum_inflow_velocity_star = 1.0


!--- mask the inlet fringe zone

!** DO j=1, 200
!**    Do i=10, 130
!**       IF (grid_set(i,j,1) .eq. 1.0) then
!**           grid_set(i,j,1) = 5.0
!**       ENDIF
!**    ENDDO
!** ENDDO

 END SUBROUTINE read_vtk_3D

 SUBROUTINE get_process_CT_geometry
 IMPLICIT NONE
 
 INTEGER :: i, j, k

 DO i=S1p,N1p
    DO j=S2p,N2p
       DO k=S3p,N3p
          CT_geometry(i,j,k) = ct_grid(iShift+i, jShift+j, kShift+k)
       ENDDO
    ENDDO
 ENDDO
  
 END SUBROUTINE get_process_CT_geometry



  SUBROUTINE read_velocity_waveform
  IMPLICIT NONE

  INTEGER              :: i
!  INTEGER, INTENT(IN)  :: number_of_digitized_points
!
!  TYPE :: time_Q_pair
!     REAL(8) :: time_in_pulse
!     REAL(8) :: Q_in_pulse
!  END TYPE time_Q_pair
!
!  TYPE(time_Q_pair), DIMENSION(number_of_digitized_points) :: waveform_data

  OPEN(UNIT=480, FILE='waveform.csv' , STATUS='old', &
     &  ACCESS ='sequential', FORM='formatted')!,recl=71781*10)


  DO i=1,number_of_digitized_points
      read(480,*) waveform_data(i)
  ENDDO
   

  END SUBROUTINE read_velocity_waveform


  SUBROUTINE convert_waveform
  IMPLICIT NONE
  INTEGER     ::    i


!  aortic_diameter = 0.024 !unit: [m] 
!  CT_pixel_width = 0.000789   ! TAVI 03
  !CT_pixel_width = 0.000779  ! TAVI 02
!  CT_pixel_width = 0.000746   ! TAVI 01
  CT_pixel_width = 0.0004   ! TAVI 00


! upsampling_factor = 2. ! TAVI003
! upsampling_factor = 2. ! TAVI002
! upsampling_factor = 2. ! TAVI001
 upsampling_factor = 1. ! TAVI000
   
 CT_pixel_width_star = (1./M1) * upsampling_factor  
 maximum_inflow_velocity_star = 1.0


 ! -- for the healthy TAVI 000 cases
   !waveform_data%Q_in_pulse = 0.33333333 * waveform_data%Q_in_pulse

!  maximum_inflow_velocity = (4. * maxval(waveform_data%Q_in_pulse) /60./1000.)/(2.*asin(1.)*aortic_diameter**2)  ! unit: [l/min]
  maximum_inflow_velocity = maxval(waveform_data%Q_in_pulse) ! Q_in_pulse from Echo holds velocity in [m/s]


  reference_velocity = maximum_inflow_velocity/maximum_inflow_velocity_star
!  reference_length   = aortic_diameter/aortic_diameter_star 
  reference_length  = CT_pixel_width/CT_pixel_width_star  

  !-- overwrite the reference values for idealized case
  IF (idealized_aorta == 1) THEN
     reference_length = 0.2
     reference_velocity = 2.7
  ENDIF


  Re = reference_velocity * reference_length / blood_kinematic_viscosity

  reference_time = reference_length/reference_velocity
  if (rank == 0) print*, reference_time,reference_velocity, reference_length, 'this is reference time, velocity, length' 
  DO i=1,number_of_digitized_points
       waveform_data_converted(i)%time_in_pulse = waveform_data(i)%time_in_pulse/(1000.*reference_time)
       ! --- convert Q to velocity in the data type 
       waveform_data_converted(i)%Q_in_pulse = waveform_data(i)%Q_in_pulse/reference_velocity 
  ENDDO

 


  END SUBROUTINE convert_waveform

  !> subroutine that calculates a pressure waveform based on a given echo-based velocity waveform and a WKM
  SUBROUTINE get_pressure_waveform_echo_WKM
  IMPLICIT NONE
  INTEGER         ::    i 
  REAL(8)            ::    q_in, time_in, q_max, vel_max, inlet_area, cardiac_output_factor
  REAL(8)            ::    C, R
  REAL(8)            ::    dt
  REAL(8)            ::    shifted_time_in_pulse


  !--- default values
  cardiac_output_factor = 1.0

  !--- inlet
  !-- TAVI 000
  R =  0.78 * 1000000.
  C =  2.166/1000000.
  !--- TAVI 000 excersie
  !R = 0.58 * 1000000. 
  !C = 1.07777/1000000.
  cardiac_output_factor = 1.
  !-- TAVI 001
  !R =  1.032 * 1000000.
  !C =  1.5597/1000000.
  !-- TAVI 002
  !R =  0.9755(0.938) * 1000000. !in paranthesis value for exactly periodic.. 
  !C =  1.9867/1000000.
  !-- TAVI 003
  !R =  1.2769 * 1000000.
  !C =  1.45/1000000.





  q_max = 0.0004 * cardiac_output_factor ! in m3/s
  vel_max = maxval(waveform_data(1:number_of_digitized_points)%Q_in_pulse)
  IF (rank == 0) print*, 'vel max in pulse is',vel_max 
  inlet_area = q_max/vel_max
  pre_inlet = 80. ! in mmHg
  DO i=1,number_of_digitized_points-1   
       call get_shifted_time_in_pulse
       if ((shifted_time_in_pulse .ge. waveform_data_converted(i)%time_in_pulse)  &
                             & .and. (shifted_time_in_pulse .le. waveform_data_converted(i+1)%time_in_pulse)) then    
       dt = waveform_data(i+1)%time_in_pulse/1000. - waveform_data(i)%time_in_pulse/1000.
       q_in = waveform_data(i)%Q_in_pulse * inlet_area
       pre_inlet = pre_inlet * (1. - dt/(C*R)) - dt / C * q_in
       endif
  ENDDO

  ! --- convert from mmHg to nondimensional
  pre_inlet = pre_inlet / (blood_density * reference_velocity**2.) * 133.322       



  !--- outlet
  !-- TAVI 000
  R =  0.898 * 1000000.
  C =  2.066/1000000.
  !-- TAVI 000 exercise
  !R = 0.6737*1000000. 
  !C = 0.95/1000000.
  !-- TAVI 001
  !R =  1.21 * 1000000.
  !C =  1.351/1000000.
  !-- TAVI 002
 ! R =  1.1225(1.0825) * 1000000. ! in paranthesis value for exactly periodic
 ! C =  1.9/1000000.
  !-- TAVI 003
 ! R =  1.498 * 1000000.
 ! C =  1.25/1000000.





  q_max = 0.0004 * 0.85 * cardiac_output_factor ! in m3/s
  vel_max = maxval(waveform_data(1:number_of_digitized_points)%Q_in_pulse)
  IF (rank == 0) print*, 'vel max in pulse is',vel_max 
  inlet_area = q_max/vel_max
  pre_outlet = 80. ! in mmHg
  DO i=1,number_of_digitized_points-1   
       call get_shifted_time_in_pulse
       if ((shifted_time_in_pulse .ge. waveform_data_converted(i)%time_in_pulse)  &
                             & .and. (shifted_time_in_pulse .le. waveform_data_converted(i+1)%time_in_pulse)) then    
       dt = waveform_data(i+1)%time_in_pulse/1000. - waveform_data(i)%time_in_pulse/1000.
       q_in = waveform_data(i)%Q_in_pulse * inlet_area
       pre_outlet = pre_outlet * (1. - dt/(C*R)) - dt / C * q_in
       endif
  ENDDO

  ! --- convert from mmHg to nondimensional
  pre_outlet = pre_outlet / (blood_density * reference_velocity**2.) * 133.322       


  if (rank == 0) print*, 'WK pressure based on echo data', pre_inlet, pre_outlet

  END SUBROUTINE get_pressure_waveform_echo_WKM


  SUBROUTINE get_pulse_number_and_length
  IMPLICIT NONE

  pulse_length = waveform_data_converted(number_of_digitized_points)%time_in_pulse 

  pulse_number = floor(time/pulse_length) + 1

  END SUBROUTINE get_pulse_number_and_length

  SUBROUTINE get_pressure_and_flow_rate_WKM(branch_number, p,q)
  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: branch_number
  REAL(8), INTENT(OUT) :: p,q
  REAL(8)              :: time_phys                   ! [in sec] the physiological time in pulse
  REAL(8)              :: T_beat, T_systole           ! [in sec] the beat length and the systolic phase length
  REAL(8)              :: Q_peak_systole              ! [in m3/s] the peak systolic flow rate
  REAL(8)              :: RWK, CWK                    ! resistance and capacitance of 2-element WK model
  REAL(8)              :: p_diastole   ! systolic pressure and diastolic pressure
  REAL(8)              :: pi

  pi = 4. *atan(1.)

  T_beat    =  1.  ! was 60./72.
  T_systole = 0.4 * T_beat

  !p_end_systole  = 120.  !* 133.322 ! 1 mmHg = 133.322 Pa 
  p_diastole     =  80.  !* 133.322


  IF (branch_number == 1) THEN  ! branch_number == 1 means inflow

     
     CWK = 1.5496/1000000.       ! was 1.1960 / 1000000. for T_beat = 60./72. 
     RWK = 0.9758*1000000.       ! was 0.9814 * 1000000.
     Q_peak_systole = 0.0004  ! [m3/s] 
     time_phys = time * reference_time
     IF (time_phys .le. T_systole) THEN
        p = (p_diastole + Q_peak_systole * T_systole * RWK * CWK * pi * RWK/(T_systole**2. +  CWK**2.*pi**2.*RWK**2.)) * &
                & exp(-1.*time_phys/(RWK*CWK)) &
                & - (T_systole * Q_peak_systole * RWK * (CWK * pi * RWK * cos(pi*time_phys/T_systole)) - &
                &   T_systole * sin(pi*time_phys/T_systole))/(T_systole**2. + CWK**2.*pi**2. * RWK**2.) 
        q = Q_peak_systole * sin(pi*(time_phys)/T_systole)
        p_end_systole_in = p
     !CALL  MPI_BCAST(p_end_systole_in,1,MPI_REAL8,0,COMM_CART,merror)
     ELSE
        p = p_end_systole_in * exp(-1.*(time_phys-T_systole)/(RWK*CWK))
        q = 0.
     ENDIF 

     !-- make non-dimenstional

     p = p / (blood_density * reference_velocity**2.) * 133.322
     q = q / (reference_length ** 3. / reference_time)

  ELSEIF (branch_number == 2) THEN  ! branch_number == 2 means outflow
     CWK = 1.5 /1000000.         ! was 1.153 / 1000000. for T_beat = 60.72 
     RWK = 1.12 * 1000000.       ! was 1.127 * 1000000. for T_beat=60/72  
     Q_peak_systole = 0.0004 * 0.85  ! [m3/s] 
     time_phys = time * reference_time
     IF (time_phys .le. T_systole) THEN
        p = (p_diastole + Q_peak_systole * T_systole * RWK * CWK * pi * RWK/(T_systole**2. +  CWK**2.*pi**2.*RWK**2.)) * &
                & exp(-1.*time_phys/(RWK*CWK)) &
                & - (T_systole * Q_peak_systole * RWK * (CWK * pi * RWK * cos(pi*time_phys/T_systole)) - &
                &   T_systole * sin(pi*time_phys/T_systole))/(T_systole**2. + CWK**2.*pi**2. * RWK**2.) 
        q = Q_peak_systole * sin(pi*time_phys/T_systole)
        p_end_systole_out = p
     ELSE
        p = p_end_systole_out * exp(-1.*(time_phys - T_systole)/(RWK*CWK))
        q = 0.
     ENDIF 

     !-- mak ep non-dimenstional

     p = p / (blood_density * reference_velocity**2.) * 133.322
     q = q / (reference_length ** 3. / reference_time)


  ELSEIF (branch_number == 3) THEN  ! branch_number == 2 means outflow
     CWK = 0.059 / 1000000. 
     RWK = 19.68 * 1000000.
     Q_peak_systole = 0.0004 * 0.05  ! [m3/s] 
     time_phys = time * reference_time
     IF (time_phys .le. T_systole) THEN
        p = (p_diastole + Q_peak_systole * T_systole * RWK * CWK * pi * RWK/(T_systole**2. +  CWK**2.*pi**2.*RWK**2.)) * &
                & exp(-1.*time_phys/(RWK*CWK)) &
                & - (T_systole * Q_peak_systole * RWK * (CWK * pi * RWK * cos(pi*time_phys/T_systole)) - &
                &   T_systole * sin(pi*time_phys/T_systole))/(T_systole**2. + CWK**2.*pi**2. * RWK**2.) 
        q = Q_peak_systole * sin(pi*time_phys/T_systole)
        p_end_systole_out_2 = p
     ELSE
        p = p_end_systole_out_2 * exp(-1.*(time_phys - T_systole)/(RWK*CWK))
        q = 0.
     ENDIF 

     !-- mak ep non-dimenstional

     p = p / (blood_density * reference_velocity**2.) * 133.322
     q = q / (reference_length ** 3. / reference_time)

  ENDIF

  END SUBROUTINE get_pressure_and_flow_rate_WKM

  SUBROUTINE apply_pressure_BC_CT(rel)
    implicit none

  integer :: ii,jj,kk
  REAL(8), INTENT(INOUT)    ::   rel(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)) 

   if ((ct_based_aorta == 1) .or. (idealized_aorta == 1)) then

      do ii=s1p,n1p
         do jj =s2p, n2p
             do kk =s3p, n3p
                if  (ct_grid(ii+ishift,jj+jshift,kk + kshift) == 5.0) then        
                     !nl(ii,jj,kk,3) =  nl(ii,jj,kk,3) - 1. *  pre_inlet/(y3p(k_start_index-1) - y3p(k_start_index-10)) 
                     !pre(ii,jj,kk) =   - (y3p(k_start_index-10) - x3p(kk))  *  pre_inlet/(y3p(k_start_index) - y3p(k_start_index-10)) 
                     if ((kk+kShift) == k_start_index-1) then
                    ! if (((kk+kShift) .le. k_start_index-6) .and. ((kk+kShift) .ge. k_start_index-6)) then
                             rel(ii,jj,kk) = pre_inlet *  (aRK(substep)+bRK(substep)) * dtime
                             !print*, 'pressure applied inlet', pre_inlet
                     endif 
                else if  (ct_grid(ii+ishift,jj+jshift,kk + kshift) == 9.0) then 
                  !   nl(ii,jj,kk,3) =  nl(ii,jj,kk,3) - 0.8 *  pre_outlet/(y3p(k_start_index-1) - y3p(k_start_index-10)) 
                     !pre(ii,jj,kk) =   -(y3p(k_start_index-10) - x3p(kk)) *  pre_outlet/(y3p(k_start_index) - y3p(k_start_index-10)) 
                     if ((kk+kShift) == k_start_index-1) then
              !               rel(ii,jj,kk) = pre_outlet * (aRK(substep)+bRK(substep)) * dtime
                             !print*, 'pressure applied outlet',  pre_outlet
                     endif
                else if  (ct_grid(ii+ishift,jj+jshift,kk + kshift) == 11.0) then 
                  !   nl(ii,jj,kk,3) =  nl(ii,jj,kk,3) - 0.8 *  pre_outlet_2/(y3p(k_start_index_2+1) - y3p(k_start_index_2+10)) 
                     !pre(ii,jj,kk) =   -(y3p(k_start_index_2+10) - x3p(kk)) *  pre_outlet_2/(y3p(k_start_index_2) - y3p(k_start_index_2+10)) 
                     !pre(ii,jj,kk) = pre_outlet_2
                endif
             end do
          end do        
      end do

   endif




  END SUBROUTINE apply_pressure_BC_CT 

  SUBROUTINE get_shifted_time_in_pulse
  IMPLICIT NONE
   
      CALL get_pulse_number_and_length
      shifted_time_in_pulse = time - (pulse_number - 1) * pulse_length
      !if (rank == 0) print*, pulse_number, time, shifted_time_in_pulse

 END SUBROUTINE get_shifted_time_in_pulse 


  SUBROUTINE apply_CT_forcing
    implicit none

  integer                ::          rank_inlet, rank_inlet_global, rank_outlet, rank_outlet_global, blocks_count, req1_in, req2_out
  INTEGER                ::          rank_outlet_2, rank_outlet_2_global
  LOGICAL                ::          inlet_checked, outlet_checked
  INTEGER                ::          status(MPI_STATUS_SIZE)
  INTEGER                ::          ii,jj,kk, iii





!--- flow forcing for CT-based DNS case
     IF (CT_based_aorta == 1) THEN
     !--- get pressure based on 1EM WKM and the ECHO waveform
     CALL get_pressure_waveform_echo_WKM

       flow_rate_correction_factor = 1.0  ! --- should be one by defualt - greater than one means higher output due to,
       !e.g.,exercise
        CALL get_shifted_time_in_pulse
        do ii = S11B, N11B
            do jj = S21B,N21B
               do kk = S31B, N31B
                  if ((ct_grid(ii+iShift, jj+jShift, kk+kShift) .lt. 0.5) .or. (ct_grid(ii+iShift, jj+jShift, kk+kShift) ==  7.0)) then
                      vel(ii,jj,kk,1) = .0
                      vel(ii,jj,kk,2) = .0
                      vel(ii,jj,kk,3) = .0
                  elseif (ct_grid(ii+iShift, jj+jShift, kk+kShift) == 5.0)  then
                      do iii = 1, number_of_digitized_points
                         if ((shifted_time_in_pulse .ge. waveform_data_converted(iii)%time_in_pulse)  &
                             & .and. (shifted_time_in_pulse .le. waveform_data_converted(iii+1)%time_in_pulse)) then
                            vel(ii,jj,kk,3) = waveform_data_converted(iii)%Q_in_pulse * flow_rate_correction_factor
                            vel(ii,jj,kk,2) = 0.
                            vel(ii,jj,kk,1) = 0.
                         endif
                      enddo
                  endif
               enddo
            enddo
         enddo
      ENDIF


     IF (idealized_aorta == 1) THEN
     rank_inlet = 0
     rank_inlet_global = 0

     rank_outlet = 0
     rank_outlet_global = 0

     rank_outlet_2 = 0
     rank_outlet_2_global = 0

     inlet_checked = .FALSE.
        !flow_rate_correction_factor = 0.33  ! --- should be one for no constriction
        !call get_pulse_number_and_length
        !shifted_time_in_pulse = time - (pulse_number - 1) * pulse_length
        !if (rank == 0) print*, pulse_number, time, shifted_time_in_pulse
        DO ii = S11B, N11B
             DO jj = S21B,N21B
                 DO kk = S31B, N31B
                     IF ((ct_grid(ii+iShift, jj+jShift, kk+kShift) .lt. 0.5) .or. (ct_grid(ii+iShift, jj+jShift, kk+kShift) ==  7.0)) then
                            vel(ii,jj,kk,1) = .0
                            vel(ii,jj,kk,2) = .0
                            vel(ii,jj,kk,3) = .0
                     ELSEIF (ct_grid(ii+iShift, jj+jShift, kk+kShift) == 5.0)  then                     
                            CALL get_pressure_and_flow_rate_WKM(1,pre_inlet, q_inlet)
                            rank_inlet = rank
                            !if (rank == 20) print*, 'we have inlet', rank, q_inlet
                            blocks_count = NB1 * NB2 * NB3
                            !pre_slope_inlet = pre_inlet/(y3p(k_start_index-1) - y3p(k_start_index-10))

                            !IF (inlet_checked == .FALSE.) THEN 
                            !   IF (rank .ne. 0) THEN
                            !      CALL MPI_SEND(p_end_systole_in,1,MPI_REAL8,0,1,COMM_CART,req1_in,merror)
                            !   ENDIF
                            !ENDIF
                            !inlet_checked = .TRUE.    
                            !blocks_count = NB1 * NB2 * NB3
                            !DO iii = 1,  blocks_count
                            !   IF iii == rank_init 

                            !IF (rank == 0) print*, 'passed'
                            !if (rank == 0) print*, p_end_systole_in
                            vel(ii,jj,kk,3) = q_inlet/(inlet_voxels_count * (L1/(M1-1))**2) 
                            vel(ii,jj,kk,2) = 0.
                            vel(ii,jj,kk,1) = 0.
                            !pre(ii,jj,kk) = pre_inlet
                     ELSEIF (ct_grid(ii+iShift, jj+jShift, kk+kShift) == 9.0)  then                     
                            CALL get_pressure_and_flow_rate_WKM(2,pre_outlet, q_outlet)
                            !print*, 'we have outlet', rank
                            rank_outlet = rank
                            !pre_slope_outlet = pre_outlet/(y3p(k_start_index-1) - y3p(k_start_index-10))
                            !CALL MPI_BCAST(p_end_systole_out,1,MPI_REAL8,rank,COMM_CART,merror)
                            !IF (rank == 4) print*, 'passed 2'
                            !vel(ii,jj,kk,3) = -q_outlet/(inlet_voxels_count * (L1/(M1-1))**2) !TODO change for outlet_voxel_count 
                            !vel(ii,jj,kk,2) = 0.
                            !vel(ii,jj,kk,1) = 0.
                            !pre(ii,jj,kk) = pre_outlet
                     ELSEIF (ct_grid(ii+iShift, jj+jShift, kk+kShift) == 11.0)  then                     
                            CALL get_pressure_and_flow_rate_WKM(3,pre_outlet_2, q_outlet_2)
                            rank_outlet_2 = rank
                            !print*, 'we have outlet_2', rank
                            !pre_slope_outlet_2 = pre_outlet_2/(y3p(k_start_index_2+1) - y3p(k_start_index_2+10))
                           !print*, rank
                     ELSEIF (ct_grid(ii+iShift, jj+jShift, kk+kShift) == 2.0) THEN
                            !pre(ii,jj,kk) = 0.
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO




       CALL MPI_ALLREDUCE(rank_inlet,rank_inlet_global,1,MPI_INT,MPI_MAX,COMM_CART,merror)    
       CALL MPI_ALLREDUCE(rank_outlet,rank_outlet_global,1,MPI_INT,MPI_MAX,COMM_CART,merror)    
       CALL MPI_ALLREDUCE(rank_outlet_2,rank_outlet_2_global,1,MPI_INT,MPI_MAX,COMM_CART,merror)    
       
       rank_outlet   = rank_outlet_global
       rank_outlet_2 = rank_outlet_2_global
       rank_inlet    = rank_inlet_global
       
       IF (rank == 0) print*, rank_inlet, rank_outlet, rank_outlet_2
       !CALL MPI_IRECV(p_end_systole_in,1,MPI_REAL8,rank,1,COMM_CART,req1_in,merror)
       !CALL MPI_WAIT(req1_in(rank), status, merror)

       CALL MPI_BCAST(p_end_systole_in,1,MPI_REAL8,rank_inlet,COMM_CART,merror)
       CALL MPI_BCAST(p_end_systole_out,1,MPI_REAL8,rank_outlet,COMM_CART,merror)
       CALL MPI_BCAST(p_end_systole_out_2,1,MPI_REAL8,rank_outlet_2,COMM_CART,merror)

       CALL MPI_BCAST(pre_inlet,1,MPI_REAL8,rank_inlet,COMM_CART,merror)
       CALL MPI_BCAST(pre_outlet,1,MPI_REAL8,rank_outlet,COMM_CART,merror)
       CALL MPI_BCAST(pre_outlet_2,1,MPI_REAL8,rank_outlet_2,COMM_CART,merror)

       if (rank == 0) print*,  pre_inlet, pre_outlet, pre_outlet_2

        ENDIF

  END SUBROUTINE apply_CT_forcing





 subroutine sphere_local
   implicit none

   integer :: i, j, k
   REAL(8)    :: xo, yo, r, dist

   xo = 0.5
   yo = 0.5
   r  = 0.05


   do i =S11,N11
     do  j =S21,N21
         do k =S31,N31
             dist = sqrt((x1u(i)-xo)**2 + (x2p(j)-yo)**2)
            if (dist .le. r)  then
              vel(i,j,:,1)=0.
            endif
            dist = sqrt((x1p(i)-xo)**2 + (x2v(j)-yo)**2)
            if (dist .le. r)  then
              vel(i,j,:,2)=0.
            endif
         enddo
     enddo
   enddo 


 end subroutine sphere_local


 subroutine sphere_global
   implicit none

   integer :: i, j, k
   REAL(8)    :: xo, yo, r, dist

   xo = 0.5
   yo = 0.5
   r  = 0.1


   do i =1,M1
     do  j =1,M2
         do k =1,M3
             dist = sqrt(r**2 - y1u(i)**2 - y2p(j)**2)
            if (dist .le. r)  then
              vel(i-iShift,j-jShift,:,1)=0.
            endif
            dist = sqrt(r**2 - y1p(i)**2 - y2v(j)**2)
            if (dist .le. r)  then
              vel(i-iShift,j-jShift,:,2)=0.
            endif
         enddo
     enddo
   enddo 


 end subroutine sphere_global


 subroutine global_indicator
 
 implicit none
 integer    ::    i, j, k


 ! -- first define 

 do i=1,M1
    x_global(i) = (i-1)/(M1-1)*L1
 enddo

 do j=1,M2
    y_global(j) = (j-1)/(M2-1)*L2
 enddo

 do k=1,M3
    z_global(k) = (k-1)/(M3-1)*L3
 enddo


 end subroutine global_indicator

 
  
  FUNCTION smooth_step(step_yes , xx) RESULT(fn_val)

  IMPLICIT NONE

  LOGICAL                :: step_yes
  REAL(8)   , INTENT(IN)    :: xx
  REAL(8)                   :: fn_val
  
  IF (step_yes) THEN
  IF       (xx .LE. 0.) THEN
    fn_val = 0.
  ELSE IF  ( (xx .GT. 0.) .AND. (xx .LT. 1.) ) THEN
    fn_val = 1./( 1. + EXP( (1./(xx-1.)) + (1./xx)) )
  ELSE IF  (xx .GE. 1.) THEN
    fn_val = 1.
  ELSE
    write(*,*)'     WARNING: invalid input to function smooth_step.'
  END IF
  ELSE
    fn_val = 1.
  END IF

  RETURN

  END FUNCTION smooth_step


  
  FUNCTION linear_step(xx) RESULT(fn_val)

  IMPLICIT NONE

  REAL(8)   , INTENT(IN)     :: xx
  REAL(8)                    :: fn_val

  IF ( xx .LE. 0 ) THEN
    fn_val = 0.
  ELSEIF ( (xx .GT. 0. ) .AND. (xx .LT. 1.) ) THEN
    fn_val = xx
  ELSEIF (xx .GE. 1.) THEN
    fn_val = 1.
  ELSE
    write(*,*)'     WARNING: invalid input to function linear_step.'
  END IF

  END FUNCTION linear_step



  !> function that computes the fringe_coeff lambda(x)
  !! @note: 
  SUBROUTINE fringe_coeff(lam_max, x_start, x_end, d_rise, d_fall, xx, lam_fringe)

  IMPLICIT NONE

  REAL(8)   , INTENT(IN)    :: lam_max
  REAL(8)   , INTENT(IN)    :: x_start, x_end
  REAL(8)   , INTENT(IN)    :: d_rise, d_fall
  REAL(8)   , INTENT(IN)    :: xx

  REAL(8)   , INTENT(OUT)   :: lam_fringe

  IF ((d_rise .NE. d_fall) .AND. (x_start .LT. L1) .AND. (x_start .LT. L2))  THEN
    !write(*,*)'      WARNING: make sure rise and fall intervals equal.'
  END IF

  lam_fringe = lam_max*( smooth_step(.TRUE.,(xx-x_start)/d_rise) - smooth_step(.TRUE.,1. + ( (xx-x_end)/d_fall )) )

  IF ( (x_end .LE. (x_start + d_rise)) .OR. (x_start .GE. (x_end - d_fall)) ) THEN
    lam_fringe = lam_fringe + lam_max
  END IF

  END SUBROUTINE fringe_coeff





  SUBROUTINE poiseuille_parabola(x_start , x_end , xx, parab_val)
 
  IMPLICIT NONE

  REAL(8)   , INTENT(IN)    :: x_start
  REAL(8)   , INTENT(IN)    :: x_end
  REAL(8)   , INTENT(IN)    :: xx
  REAL(8)                   :: a, b, c

  REAL(8)   , INTENT(OUT)   :: parab_val


  IF ((xx .GE. x_start) .AND. (xx .LE. x_end)) THEN
    a = -4.                 /((x_end - x_start)**2.)
    b = (4.*(x_end+x_start))/((x_end - x_start)**2.)
    c = -(4.*x_end*x_start) /((x_end - x_start)**2.)

    parab_val = a*(xx**2.) + (b*xx) + c
  ELSE
    parab_val = 0.
    
  END IF


  END SUBROUTINE poiseuille_parabola



  SUBROUTINE poiseuille_paraboloid(xx, yy, x_center, y_center, radius, parab_val)

  IMPLICIT NONE

  REAL(8)   , INTENT(IN)    :: xx, yy
  REAL(8)   , INTENT(IN)    :: x_center, y_center
  REAL(8)   , INTENT(IN)    :: radius

  REAL(8)   , INTENT(OUT)   :: parab_val


  parab_val = 0.

  IF ( ((xx - x_center)**2. + (yy - y_center)**2.) .LE. radius**2. ) THEN

    parab_val = -1. * ( ( (xx - x_center)**2. + (yy - y_center)**2. ) &
                         /(radius**2.) - 1. )

  END IF


  END SUBROUTINE poiseuille_paraboloid



  SUBROUTINE windkessel_integration_step

  IMPLICIT NONE

  REAL(8)                  :: WK_predot_old(1:3)

  WK_predot_old = WK_predot

  IF (WK_type .EQ. 3) THEN
    CALL pdot3EWK( WK_pre(2) , WK_predot(2) )
  ELSE IF (WK_type .EQ. 4) THEN
    CALL pdot4EWK( WK_pre , WK_predot )
  ELSE
    WRITE(*,*)'WARNING: Invalid Windkessel model specified. Aborting ...'
    CALL MPI_FINALIZE(merror)
    STOP
  END IF

  WK_pre = WK_pre + dtime*aRK(substep)*WK_predot + dtime*bRK(substep)*WK_predot_old
  
  !--- Algebraic equation for 3EWK, Q_AV is global and already computed ---
  IF (WK_type .EQ. 3) WK_pre(1) = WK_pre(2) + R_c*Q_AV

!*** debugging
!print*, subtime, WK_pre(1), WK_pre(2), WK_pre(3), Q_AV, dQ_AV_dt
  END SUBROUTINE windkessel_integration_step






  !> subroutine that returns the RHS of the 3EWK model
  SUBROUTINE pdot3EWK( pb , dpb_dt )

  IMPLICIT NONE

  REAL(8)   , INTENT(IN)    :: pb
  REAL(8)   , INTENT(OUT)   :: dpb_dt

  Q_AV = flow_and_rate( WK_flow_dir , WK_flow_pos , WK_flow_center , WK_flow_radius )
  Q_AV = Q_AV*U_ref*L_ref**2. ! make dimensional
  !Q_AV = fitted_flow(MOD(subtime,1.))

  dpb_dt = -1./(R_p * C_art) * pb + Q_AV/C_art

  END SUBROUTINE pdot3EWK




  !> subroutine that returns the RHS of the 4EWK model
  SUBROUTINE pdot4EWK( p , dp_dt )

  IMPLICIT NONE

  REAL(8)   , INTENT(IN)    :: p(1:3)
  REAL(8)   , INTENT(OUT)   :: dp_dt(1:3)
  REAL(8)                   :: Q_AV_old

  !Q_AV_old = Q_AV
  !Q_AV     = fitted_flow(MOD(subtime,1.))

!=== This computes the flow rate by integrating the RHS ===
  dQ_AV_dt = flow_and_rate( WK_flow_dir , WK_flow_pos , WK_flow_center , WK_flow_radius )
  dQ_AV_dt = dQ_AV_dt*L_ref*U_ref**2. ! make dimensional  
!=== This computes the flow rate from the flow by FD differentiation ===
!--- forward Euler ---
!  dQ_AV_dt = (Q_AV - Q_AV_old) / (dtime*(aRK(substep)+bRK(substep)))

!--- inverted RK3 scheme ---
!  IF (substep == 1) THEN
!    Q_AV_0   = Q_AV
!    !Q_AV    = flow( WK_flow_dir , WK_flow_pos , WK_flow_center , WK_flow_radius )
!    Q_AV_1   = fitted_flow(MOD(subtime,1.))
!    dQ_AV_dt = (Q_AV_1 - Q_AV_0)/(dtime*aRK(1))
!  ELSE IF (substep == 2) THEN
!    Q_AV_2   = fitted_flow(MOD(subtime,1.))
!    dQ_AV_dt = (-(aRK(1)+bRK(2))*(Q_AV_1-Q_AV_0)/(aRK(1)*aRK(2)) + (Q_AV_2 - Q_AV_0)/aRK(2) )/dtime
!  ELSE IF (substep == 3) THEN
!    Q_AV_3   = fitted_flow(MOD(subtime,1.))
!    dQ_AV_dt = ( bRK(3)*(aRK(1) + bRK(2))*(Q_AV_1 - Q_AV_0)/(aRK(1)*aRK(2)*aRK(3)) &
!                       -(aRK(2) + bRK(3))*(Q_AV_2 - Q_AV_0)/(aRK(2)*aRK(3)) &
!                       +                  (Q_AV_3 - Q_AV_0)/aRK(3))/dtime
!    Q_AV   = Q_AV_3
!  END IF
!========================================================================

  dp_dt(1) = -R_c/L_blood*p(1) + R_c/L_blood*p(2) +               p(3) + R_c * dQ_AV_dt 
  dp_dt(2) =                                                      p(3)
  dp_dt(3) =                                      -1./(R_p*C_art)*p(3) + dQ_AV_dt / C_art

  END SUBROUTINE pdot4EWK



  !> function returns the correct flow or flow rate depending on the WK type at a desired location
  !! WK4 -> flow rate; WK3 -> flow
  FUNCTION flow_and_rate(dir , pos , center , radius) RESULT(fn_val)

  IMPLICIT NONE
 
  INTEGER, INTENT(IN)    :: dir
  REAL(8)   , INTENT(IN)    :: pos
  REAL(8)   , INTENT(IN)    :: center(1:3) !< center(3)==1 in 2D
  REAL(8)   , INTENT(IN)    :: radius
  REAL(8)                   :: fn_val
  REAL(8)                   :: fn_val_global

  IF (dimens .EQ. 2) THEN
    fn_val = flow_and_rate2D(dir, pos, center, radius)
  ELSE IF (dimens .EQ. 3) THEN
    fn_val = flow_and_rate3D(dir, pos, center, radius)
  END IF

  !--- Sum up the flow contribution of all MPI processes
  fn_val_global = 0.
  CALL MPI_ALLREDUCE(fn_val, fn_val_global, 1, MPI_REAL8, MPI_SUM, COMM_CART, merror)
  fn_val = fn_val_global

  END FUNCTION flow_and_rate


  !> function returns the 2D flow along an axis
  FUNCTION flow_and_rate2D(dir, pos , center, radius) RESULT(fn_val)

  IMPLICIT NONE

  INTEGER, INTENT(IN)    :: dir
  REAL(8)   , INTENT(IN)    :: pos
  REAL(8)   , INTENT(IN)    :: center(1:3) !< center(3)==1 in 2D
  REAL(8)   , INTENT(IN)    :: radius
  REAL(8)                   :: x1start, x1end, x2start, x2end
  INTEGER                :: x1loc, x2loc
  REAL(8)                   :: fn_val
  INTEGER                :: i,j,k


  fn_val = 0.
 
  !=== in x-direction =======================================================================================
  IF (dir .EQ. 1) THEN
    x2start = center(2) - radius
    x2end   = center(2) + radius
    !--- 3EWK -> return flow --------------------------------------------------------------------------------
    IF (WK_type .EQ. 3) THEN
      IF ((pos .GE. x1p(S1p)) .AND. (pos .LE. x1p(N1p))) THEN ! check if desired position is in process
        x1loc   = MINLOC(ABS(x1p - pos), 1) !< finds the index closest to the specified location

        DO k = S3p, N3p ! this should be 1 anyway in 2D
          DO j = S2p, N2p
            DO i = x1loc, x1loc
              IF ( (x2p(j) .GE. x2start) .AND. (x2p(j) .LE. x2end) ) THEN 
                fn_val     = fn_val + work1(i,j,k)*dx1p(i)*dx2p(j)*dx3p(k)
              END IF
            END DO 
          END DO
        END DO
      ELSE
       fn_val = 0.
      END IF
    !--- 4EWK -> return flow rate ---------------------------------------------------------------------------
    ELSE IF (WK_type .EQ. 4) THEN
      IF ((pos .GE. x1u(S11B)) .AND. (pos .LE. x1u(N11B))) THEN ! check if desired position is in process
        x1loc   = MINLOC(ABS(x1u - pos), 1) !< finds the index closest to the specified location
     
        DO k = S31B, N31B
          DO j = S21B, N21B
            DO i = x1loc, x1loc
              IF ( (x2p(j) .GE. x2start) .AND. (x2p(j) .LE. x2end) ) THEN
                fn_val     = fn_val + rhs(i,j,k,1)*dx1u(i)*dx1p(j)*dx3p(k)
              END IF
            END DO
          END DO
        END DO
      ELSE
        fn_val = 0.
      END IF
    END IF
  !=== in y-direction ======================================================================================= 
  ELSE IF (dir .EQ. 2) THEN
    x1start = center(1) - radius
    x1end   = center(1) + radius
    !--- 3EWK -> return flow --------------------------------------------------------------------------------
    IF (WK_type .EQ. 3) THEN
      IF ((pos .GE. x2p(S2p)) .AND. (pos .LE. x2p(N2p))) THEN ! check if desired position is in process
        x2loc   = MINLOC(ABS(x2p - pos), 1) !< finds the index closest to the specified location 

        DO k = S3p, N3p ! this should be 1 anyway in 2D
          DO j = x2loc, x2loc
            DO i = S1p, N1p
              IF ( (x1p(i) .GE. x1start) .AND. (x1p(i) .LE. x1end) ) THEN
                fn_val     = fn_val + work2(i,j,k)*dx1p(i)*dx2p(j)*dx3p(k)
              END IF
            END DO
          END DO
        END DO
      ELSE
        fn_val = 0.
      END IF
    !--- 4EWK -> return flow rate ---------------------------------------------------------------------------
    ELSE IF (WK_type .EQ. 4) THEN
      IF ((pos .GE. x2v(S22B)) .AND. (pos .LE. x2v(N22B))) THEN ! check if desired position is in process
        x2loc   = MINLOC(ABS(x2v - pos), 1) !< finds the index closest to the specified location 

        DO k = S32B, N32B ! this should be 1 anyway in 2D
          DO j = x2loc, x2loc
            DO i = S12B, N12B
              IF ( (x1p(i) .GE. x1start) .AND. (x1p(i) .LE. x1end) ) THEN
                fn_val     = fn_val + rhs(i,j,k,2)*dx1p(i)*dx2v(j)*dx3p(k)
              END IF
            END DO
          END DO
        END DO
      ELSE
        fn_val = 0.
      END IF
    END IF
  !=== in z-direction: NOT POSSIBLE =========================================================================
  ELSE IF (dir .EQ. 3 .AND. dimens .NE. 3) THEN
     WRITE(*,*)'WARNING: cannot compute the flux in 3-direction for a 2D flow. Aborting...'
     CALL MPI_FINALIZE(merror)
     STOP
  !==========================================================================================================
  ELSE
     WRITE(*,*)'WARNING: something went wrong with calculating the flow. Aborting...'
     CALL MPI_FINALIZE(merror)
     STOP
  END IF

  END FUNCTION flow_and_rate2D



  !> returns the flow in 3D through a circular opening in the direction of an axis
  FUNCTION flow_and_rate3D(dir , pos , center , radius) RESULT(fn_val)

  IMPLICIT NONE

  INTEGER, INTENT(IN)    :: dir
  REAL(8)   , INTENT(IN)    :: pos
  REAL(8)   , INTENT(IN)    :: center(1:3) !< center(3)==1 in 2D
  REAL(8)   , INTENT(IN)    :: radius
  REAL(8)                   :: x1start, x1end, x2start, x2end, x3start, x3end
  INTEGER                :: x1loc, x2loc, x3loc
  REAL(8)                   :: fn_val
  INTEGER                :: i,j,k

  !==========================================================================================================
  ! Depending on wether we use 4EWK or 3EWK we need either the flow or the flow rate to be returned.
  ! We make use of the fact that the integral over du/dt = rhs equals the flow rate.
  ! Because the components u,v,w are stored in containers work1,2,3 on the pressure grid, we can use the 
  ! pressure indices to loop. However, the RHS is only given on the velocity grids, i.e. we need to loop 
  ! over different indices. It may at a certain point be nicer to reduce the code and include the WK type 
  ! checking into the triple loops and use vel(:,:,:,i) insted of the worki(:,:,:)
  !==========================================================================================================


  fn_val = 0.

  !=== in x-direction =======================================================================================
  IF (dir .EQ. 1) THEN
    !--- 3EWK -> return flow --------------------------------------------------------------------------------
    IF (WK_type .EQ. 3) THEN
      IF ((pos .GE. x1p(S1p)) .AND. (pos .LE. x1p(N1p))) THEN ! check if desired position is in process
        x1loc = MINLOC(ABS(x1p - pos), 1) !< finds the index closest to the specified location

        DO k = S3p, N3p
          DO j = S2p, N2p
            DO i = x1loc, x1loc
              IF ( (x2p(j) - center(2))**2 + (x3p(k) - center(3))**2 .LE. radius**2) THEN
                fn_val     = fn_val + work1(i,j,k)*dx1p(i)*dx2p(j)*dx3p(k)
              END IF
            END DO
          END DO
        END DO
      ELSE
        fn_val = 0.
      END IF
    !--- 4EWK -> return flow rate ---------------------------------------------------------------------------
    ELSE IF (WK_type .EQ. 4) THEN
      IF ((pos .GE. x1u(S11B)) .AND. (pos .LE. x1u(N11B))) THEN ! check if desired position is in process
        x1loc = MINLOC(ABS(x1u - pos), 1) !< finds the index closest to the specified location

        DO k = S31B, N31B
          DO j = S21B, N21B
            DO i = x1loc, x1loc
              IF ( (x2p(j) - center(2))**2 + (x3p(k) - center(3))**2 .LE. radius**2) THEN
                fn_val     = fn_val + rhs(i,j,k,1)*dx1u(i)*dx2p(j)*dx3p(k)
              END IF
            END DO
          END DO
        END DO
      ELSE
        fn_val = 0.
      END IF
    END IF
  
  !=== in y-direction =======================================================================================
  ELSE IF (dir .EQ. 2) THEN
    !--- 3EWK -> return flow --------------------------------------------------------------------------------
    IF (WK_type .EQ. 3) THEN
      IF ((pos .GE. x2p(S2p)) .AND. (pos .LE. x2p(N2p))) THEN ! check if desired position is in process
        x2loc = MINLOC(ABS(x2p - pos), 1) !< finds the index closest to the specified location

        DO k = S3p, N3p
          DO j = x2loc, x2loc
            DO i = S1p, N1p
              IF ( (x1p(i) - center(1))**2 + (x3p(k) - center(3))**2 .LE. radius**2) THEN
                fn_val     = fn_val + work2(i,j,k)*dx1p(i)*dx2p(j)*dx3p(k)
              END IF
            END DO
          END DO
        END DO
      ELSE
        fn_val = 0.
      END IF
    !--- 4EWK -> return flow rate ---------------------------------------------------------------------------
    ELSE IF (WK_type .EQ. 4) THEN
      IF ((pos .GE. x2v(S22B)) .AND. (pos .LE. x2v(N22B))) THEN ! check if desired position is in process
        x2loc = MINLOC(ABS(x2v - pos), 1) !< finds the index closest to the specified location

        DO k = S32B, N32B
          DO j = x2loc, x2loc
            DO i = S12B, N12B
              IF ( (x1p(i) - center(1))**2 + (x3p(k) - center(3))**2 .LE. radius**2) THEN
                fn_val     = fn_val + rhs(i,j,k,2)*dx1p(i)*dx2v(j)*dx3p(k)
              END IF
            END DO
          END DO
        END DO
      ELSE
        fn_val = 0.
      END IF
    END IF

  !=== in z-direction =======================================================================================
  ELSE IF (dir .EQ. 3) THEN
    !--- 3EWK -> return flow --------------------------------------------------------------------------------
    IF (WK_type .EQ. 3) THEN
      IF ((pos .GE. x3p(S3p)) .AND. (pos .LE. x3p(N3p))) THEN ! check if desired position is in process
        x3loc = MINLOC(ABS(x3p - pos), 1) !< finds the index closest to the specified location

        DO k = x3loc, x3loc
          DO j = S2p, N2p
            DO i = S1p, N1p
              IF ( (x1p(i) - center(1))**2 + (x2p(j) - center(2))**2 .LE. radius**2) THEN
                fn_val     = fn_val + work3(i,j,k)*dx1p(i)*dx2p(j)*dx3p(k)
              END IF
            END DO
          END DO
        END DO
      ELSE
        fn_val = 0.
      END IF
    !--- 4EWK -> return flow rate ---------------------------------------------------------------------------
    ELSE IF (WK_type .EQ. 4) THEN
      IF ((pos .GE. x3w(S33B)) .AND. (pos .LE. x3w(N33B))) THEN ! check if desired position is in process
        x3loc = MINLOC(ABS(x3w - pos), 1) !< finds the index closest to the specified location

        DO k = x3loc, x3loc
          DO j = S23B, N23B
            DO i = S13B, N13B
              IF ( (x1p(i) - center(1))**2 + (x2p(j) - center(2))**2 .LE. radius**2) THEN
                fn_val     = fn_val + rhs(i,j,k,3)*dx1p(i)*dx2p(j)*dx3w(k)
              END IF
            END DO
          END DO
        END DO
      ELSE
        fn_val = 0.
      END IF
    END IF

  !===========================================================================================================
  ELSE
     WRITE(*,*)'WARNING: something went wrong with calculating the flow. Aborting...'
     CALL MPI_FINALIZE(merror)
     STOP
  END IF

  END FUNCTION flow_and_rate3D



  FUNCTION fitted_flow(xx) RESULT(fn_val)

  IMPLICIT NONE

  REAL(8)                   ::  alpha1, beta1, gamma1
  REAL(8)                   ::  alpha2, beta2, gamma2
  REAL(8)                   ::  fn_val
  REAL(8), INTENT(IN)       ::  xx

  alpha1 = 0.0003117
  beta1  = 0.3399
  gamma1 = 0.04485
  alpha2 = 0.0004466
  beta2  = 0.2415
  gamma2 = 0.07956
  
  fn_val = alpha1*EXP(-((xx-beta1)/gamma1)**2) + alpha2*EXP(-((xx-beta2)/gamma2)**2)

  END FUNCTION fitted_flow

  
  
  
  FUNCTION fitted_pressure(xx) RESULT(fn_val)
 
  IMPLICIT NONE
 
  REAL(8)                   :: alpha1, beta1, gamma1
  REAL(8)                   :: alpha2, beta2, gamma2
  REAL(8)                   :: fn_val
  REAL(8), INTENT(IN)       :: xx

  alpha1 = 1.647e04
  beta1  = 0.2898
  gamma1 = 0.2067
  alpha2 = 2.637e16
  beta2  = 8.607
  gamma2 = 1.387

  fn_val = alpha1*EXP(-((xx-beta1)/gamma1)**2) + alpha2*EXP(-((xx-beta2)/gamma2)**2)
  
  END FUNCTION fitted_pressure



  SUBROUTINE apply_fringe_forcing

  IMPLICIT NONE

  INTEGER            :: i,j,k
  REAL(8)               :: lamb_fringe, parab
  REAL(8), PARAMETER    :: pi = 2.*ABS(ACOS(0.))


  !--- x-component ------------------------------------------------------------------------------------------
  DO k = S31, N31
     DO j = S21, N21
        DO i = S11, N11
          !--- forcing in x-direction ---
          IF (fringe_dir .EQ. 1) THEN
            !--- only force within desired bounds ---
            IF (dimens .NE. 3) fringe_center(3) = x3p(k) ! transform circle equation to interval
            IF (((x2p(j) - fringe_center(2))**2 + (x3p(k) - fringe_center(3))**2) .LE. fringe_radius**2) THEN
              CALL fringe_coeff(fringe_amp , fringe_start , fringe_end , fringe_rise , &
                    fringe_fall , x1u(i) , lamb_fringe)
              CALL poiseuille_paraboloid(x2p(j), x3p(k), fringe_center(2), &
                    fringe_center(3), fringe_radius, parab)
              
              !nl(i,j,k,1) = nl(i,j,k,1) - lamb_fringe*(fitted_pressure(MOD(subtime,1.)) - pre(i,j,k)) &
              !                / (rho_blood*U_ref**2.) &
              !                / (fringe_end - fringe_start) &
              !                * (2./(pi*fringe_radius**4.)) * parab
              nl(i,j,k,1) = nl(i,j,k,1) - lamb_fringe*(1. - vel(i,j,k,1))


              !nl(i,j,k,1) = nl(i,j,k,1) - lamb_fringe*(  (2./(pi*fringe_radius**4.))*parab &
              !              * fitted_flow(MOD(subtime,1.))/(U_ref*L_ref**2) - vel(i,j,k,1) )  !&
!                            *smooth_step(step_yes, subtime/t_rampup)! Poiseuille profile 
              
              !nl(i,j,k,1) = nl(i,j,k,1) - lamb_fringe*( (L1/ (( fringe_end-fringe_start )))*(2./Re)  )! poiseuille channel 
            END IF
          END IF
          !--- forcing in y-direction ---
          IF (fringe_dir .EQ. 2) THEN
            !--- only force within desired bounds ---
            IF (dimens .NE. 3) fringe_center(3) = x3p(k) ! transforms circle equation to interval
            IF (((x1p(i) - fringe_center(1))**2 + (x3p(k) - fringe_center(3))**2) .LE. fringe_radius**2) THEN
              CALL fringe_coeff(fringe_amp , fringe_start , fringe_end , fringe_rise , &
                    fringe_fall , x2p(j) , lamb_fringe)
              nl(i,j,k,1) = nl(i,j,k,1) - 100.*lamb_fringe*( 0. - vel(i,j,k,1) )
            END IF
          END IF
          !--- forcing in z-direction ---
          IF (fringe_dir .EQ. 3 .AND. dimens .EQ. 3) THEN
            !--- only force within desired bounds ---
            IF (((x1p(i) - fringe_center(1))**2 + (x2p(j) - fringe_center(2))**2) .LE. fringe_radius**2) THEN
              CALL fringe_coeff(fringe_amp , fringe_start , fringe_end , fringe_rise , &
                    fringe_fall , x3p(k) , lamb_fringe)
              nl(i,j,k,1) = nl(i,j,k,1) - 100.*lamb_fringe*( 0. - vel(i,j,k,1) )
            END IF
          END IF
        END DO
     END DO
  END DO
  !--- y-component ------------------------------------------------------------------------------------------
  DO k = S32, N32
     DO j = S22, N22
        DO i = S12, N12
          !--- forcing in x-direction ---
          IF (fringe_dir .EQ. 1) THEN
            !--- only force within desired bounds ---
            IF (dimens .NE. 3) fringe_center(3) = x3p(k) ! transform circle equation to interval
            IF (((x2p(j) - fringe_center(2))**2 + (x3p(k) - fringe_center(3))**2) .LE. fringe_radius**2) THEN
              CALL fringe_coeff(fringe_amp , fringe_start , fringe_end , fringe_rise , &
                    fringe_fall , x1p(i) , lamb_fringe)
              nl(i,j,k,2) = nl(i,j,k,2) - 100.*lamb_fringe*( 0. - vel(i,j,k,2) )
            END IF
          END IF
          !--- forcing in y-direction ---
          IF (fringe_dir .EQ. 2) THEN
            !--- only force within desired bounds ---
            IF (dimens .NE. 3) fringe_center(3) = x3p(k) ! transforms circle equation to interval
            IF (((x1p(i) - fringe_center(1))**2 + (x3p(k) - fringe_center(3))**2) .LE. fringe_radius**2) THEN
              CALL fringe_coeff(fringe_amp , fringe_start , fringe_end , fringe_rise , &
                    fringe_fall , x2v(j) , lamb_fringe)
              CALL poiseuille_parabola(x1p(S1p),x1p(N1p),x1p(i),parab)


              nl(i,j,k,2) = nl(i,j,k,2) - lamb_fringe*( parab  - vel(i,j,k,2) ) & 
                             *smooth_step(step_yes, subtime/t_rampup)! Poiseuille profile 
            
              nl(i,j,k,2) = nl(i,j,k,2) - lamb_fringe*( (L2/ (( fringe_end-fringe_start )))*(2./Re)  )! poiseuille channel 
            END IF
          END IF
          !--- forcing in z-direction ---
          IF (fringe_dir .EQ. 3 .AND. dimens .EQ. 3) THEN
            !--- only force within desired bounds ---
            IF (((x1p(i) - fringe_center(1))**2 + (x2p(j) - fringe_center(2))**2) .LE. fringe_radius**2) THEN
              CALL fringe_coeff(fringe_amp , fringe_start , fringe_end , fringe_rise , &
                    fringe_fall , x3p(k) , lamb_fringe)
              nl(i,j,k,2) = nl(i,j,k,2) - 100.*lamb_fringe*( 0. - vel(i,j,k,2) ) 
            END IF
          END IF
        END DO
     END DO
  END DO
  !--- z-component ------------------------------------------------------------------------------------------
  IF (dimens .EQ. 3) THEN
    DO k = S33, N33
      DO j = S23, N23
        DO i = S13, N13
          !--- forcing in x-direction ---
          IF (fringe_dir .EQ. 1) THEN
            !--- only force within desired bounds ---
            IF (((x2p(j) - fringe_center(2))**2 + (x3p(k) - fringe_center(3))**2) .LE. fringe_radius**2) THEN
              CALL fringe_coeff(fringe_amp , fringe_start , fringe_end , fringe_rise , &
                    fringe_fall , x1p(i) , lamb_fringe)
              nl(i,j,k,3) = nl(i,j,k,3) - 100.*lamb_fringe*( 0. - vel(i,j,k,3) )
            END IF
          END IF
          !--- forcing in y-direction ---
          IF (fringe_dir .EQ. 2) THEN
            !--- only force within desired bounds ---
            IF (((x1p(i) - fringe_center(1))**2 + (x3p(k) - fringe_center(3))**2) .LE. fringe_radius**2) THEN
                CALL fringe_coeff(fringe_amp , fringe_start , fringe_end , fringe_rise , &
                      fringe_fall , x2p(j) , lamb_fringe)
                nl(i,j,k,3) = nl(i,j,k,3) - 100.*lamb_fringe*( 0. - vel(i,j,k,3) )
            END IF
          END IF
          !--- forcing in z-direction ---
          IF (fringe_dir .EQ. 3) THEN
            !--- only force within desired bounds ---
            IF (((x1p(i) - fringe_center(1))**2 + (x2p(j) - fringe_center(2))**2) .LE. fringe_radius**2) THEN
                CALL fringe_coeff(fringe_amp , fringe_start , fringe_end , fringe_rise , &
                      fringe_fall , x3w(k) , lamb_fringe)
                ! something for nl(i,j,k,3) ...
                nl(i,j,k,3) = nl(i,j,k,3) - lamb_fringe*( 1. - vel(i,j,k,3))
                nl(i,j,k,3) = nl(i,j,k,3) - lamb_fringe*( (L3/ (( fringe_end-fringe_start )))*(2./Re)  )! poiseuille channel 
            END IF
          END IF
        END DO
      END DO
    END DO
  END IF
  !----------------------------------------------------------------------------------------------------------

  END SUBROUTINE apply_fringe_forcing




  SUBROUTINE apply_windkessel_loading

  IMPLICIT NONE

  INTEGER            :: i,j,k
  REAL(8)               :: lamb_fringe
  INTEGER            :: x1loc, x2loc, x3loc

  !--- x-direction ------------------------------------------------------------------------------------------
  IF (WK_flow_dir .EQ. 1) THEN
    DO k = S31, N31
      DO j = S21, N21
        DO i = S11, N11
          CALL fringe_coeff(WK_frge_amp , WK_frge_start , WK_frge_end , WK_frge_rise , &
                WK_frge_fall , x1u(i) , lamb_fringe)

          IF ( dimens .NE. 3) WK_flow_center(3) = x3p(k) ! transforms circle equation into linear interval
          IF ( ((x2p(j) - WK_flow_center(2))**2 + (x3p(k) - WK_flow_center(3))**2 .LE. WK_flow_radius**2 ) ) THEN
                 nl(i,j,k,1) = nl(i,j,k,1) + lamb_fringe*( WK_pre(1) ) &
                             / ( WK_frge_end - WK_frge_start ) &
                             / ( rho_blood * U_ref**2. )
!                             / ((Re**2.)*(mu_blood**2.))*(rho_blood*(L_ref**2.))
          END IF
        END DO
      END DO
    END DO
  END IF
  !--- y-direction ------------------------------------------------------------------------------------------
  IF (WK_flow_dir .EQ. 2) THEN
    DO k = S32, N32
      DO j = S22, N22
        DO i = S12, N12
          CALL fringe_coeff(WK_frge_amp , WK_frge_start , WK_frge_end , WK_frge_rise , &
                WK_frge_fall , x2v(j) , lamb_fringe)

          IF ( dimens .NE. 3) WK_flow_center(3) = x3p(k) ! transforms circle equation into linear interval
          IF ( ((x1p(i) - WK_flow_center(1))**2 + (x3p(k) - WK_flow_center(3))**2 .LE. WK_flow_radius**2 ) ) THEN
                 nl(i,j,k,2) = nl(i,j,k,2) + lamb_fringe*( WK_pre(1) ) &
                             / ( WK_frge_end - WK_frge_start ) &
                             / ( rho_blood * U_ref**2. )
!                             / ((Re**2.)*(mu_blood**2.))*(rho_blood*(L_ref**2.))
          END IF
        END DO
      END DO
    END DO
  END IF
  !--- z-direction ------------------------------------------------------------------------------------------
  IF (WK_flow_dir .EQ. 3 .AND. dimens .EQ. 3) THEN
    DO k = S33, N33
      DO j = S23, N23
        DO i = S13, N13
          CALL fringe_coeff(WK_frge_amp , WK_frge_start , WK_frge_end , WK_frge_rise , &
                WK_frge_fall , x3w(k) , lamb_fringe)

          IF ( ((x1p(i) - WK_flow_center(1))**2 + (x2p(j) - WK_flow_center(2))**2 .LE. WK_flow_radius**2 ) ) THEN
                 nl(i,j,k,3) = nl(i,j,k,3) + lamb_fringe*( WK_pre(1) ) & 
                             /( WK_frge_end - WK_frge_start ) &
                             / ( rho_blood * U_ref**2. )
!                             / ((Re**2.)*(mu_blood**2.))*(rho_blood*(L_ref**2.))
          END IF
        END DO
      END DO
    END DO
  END IF
  !----------------------------------------------------------------------------------------------------------


  END SUBROUTINE apply_windkessel_loading





  !> function that provides a smoothed step function (erf) from 0 (if input<-cutoff=-3) to 1 (if input>cutoff=3)
  !! @param[in] xx coordinate
  FUNCTION interface(xx) RESULT(fn_val)
  ! (sample function)
  
  IMPLICIT NONE
  
  REAL(8)   , INTENT(in   ) ::  xx
  REAL(8)                   ::  fn_val
  REAL(8)                   ::  cutoff
  
  
  ! underflow ...
  cutoff = 3.
  
  
  IF      (xx .GT.  cutoff) THEN
     fn_val = 1.
  ELSE IF (xx .LT. -cutoff) THEN
     fn_val = 0.
  ELSE
     fn_val = 0.5*(1.+erf(SQRT(pi)*xx))
  END IF
  
  
  END FUNCTION interface
  
  
  
  
  
  
  
  
  
  
  !> subroutine that transforms and maps coordinates to user-specified distributions
  !! @param[in] Lmax physical domain extent
  !! @param[in] iimax Number of gridpoints in the spacial dimension
  !! @param[in] ii0L lower fix bound of coordinate transform
  !! @param[in] ii0U upper fix bound of coordinate transform
  !! @param[in] ii current discrete coordinate point
  !! @param[out] xx transformed and mapped coordinate
  !! @param[out] dx immediate transformed and mapped coordinate step
  SUBROUTINE coord_tanh(Lmax,iimax,ii0L,ii0U,ii,xx,dx)
  ! (sample subroutine)
  
  IMPLICIT NONE
  
  REAL(8)   , INTENT(in)    ::  Lmax
  REAL(8)   , INTENT(in)    ::  iimax
  
  REAL(8)   , INTENT(in)    ::  ii0L
  REAL(8)   , INTENT(in)    ::  ii0U
  
  REAL(8)   , INTENT(in)    ::  ii
  REAL(8)   , INTENT(out)   ::  xx
  REAL(8)   , INTENT(out)   ::  dx
  
  REAL(8)                   ::  yy, cmin, cmax
  
  
  IF (ii0U == ii0L) THEN
     ! equidistant grid:
     xx = (ii-1.)*Lmax/(iimax-1.)
     dx =         Lmax/(iimax-1.)
  ELSE
     cmax =  TANH(ii0U)
     cmin = -TANH(ii0L)
     
     ! coordinate transformation (index i=1. is the origin):
     ! y = (i-1.)/(imax-1.)*(cmax-cmin) + cmin
     yy = (ii-1.)/(iimax-1.)*(cmax-cmin) + cmin
     
     ! mapping funktion f(yy)
     ! x = L * (f(y)-f(cmin)) / (f(cmax)-f(cmin))
     xx = Lmax * (atanh(yy)-atanh(cmin)) / (ii0U+ii0L)
     
     ! dx/di = L / (f(cmax)-f(cmin)) * dy/di                 * df(y(i))/dy
     !       = L / (f(cmax)-f(cmin)) * (cmax-cmin)/(imax-1.) * df(y(i))/dy
     dx = Lmax / (atanh(cmax)-atanh(cmin)) * (cmax-cmin)/(iimax-1.) * 1./(1.-yy**2)
  END IF
  
  
  END SUBROUTINE coord_tanh
  
  
  
  
  
  
  
  
  
  
  !> function that calculates the inverse of the hyperbolic tangent. This function has not
  !! been made an intrinsic fortran function until the 2008 standard. The way it is calculated
  !! here, only the REAL(8) parts are corresponding.
  !! @warning beware of the imaginary part if ever needed!
  FUNCTION atanh(x) RESULT(fn_val)
  ! (sample function)
  
  IMPLICIT NONE
  
  REAL(8)   , INTENT(IN)  :: x
  REAL(8)                 :: fn_val
  
  
  IF (x == 0.) THEN
     fn_val = 0.
     RETURN
  END IF
  
  fn_val = 0.5* LOG((1.+x)/(1.-x))
  
  
  END FUNCTION atanh
  
  
  
  
  
  
  
  
  
  
  !> Evaluation of the REAL(8) error function. Based upon a Fortran 66 routine in the Naval Surface Warfare
  !! Center's Mathematics Library (1993 version). Adapten by Alan.Miller@vic.cmis.csiro.au.
  FUNCTION erf(x) RESULT(fn_val)
  ! (sample function)
  
  !-----------------------------------------------------------------------
  !             EVALUATION OF THE REAL(8) ERROR FUNCTION
  ! Based upon a Fortran 66 routine in the Naval Surface Warfare Center's
  ! Mathematics Library (1993 version).
  ! Adapted by Alan.Miller @ vic.cmis.csiro.au
  !-----------------------------------------------------------------------
  
  IMPLICIT NONE
  
  REAL(8)   , INTENT(IN)  ::  x
  REAL(8)                 ::  fn_val
  
  REAL(8)   , PARAMETER   ::  c = .564189583547756, one = 1.0, half = 0.5, zero = 0.0
  REAL(8)   , PARAMETER   ::  &
           a(5) = (/  7.71058495001320D-05, -1.33733772997339D-03, 3.23076579225834D-02, 4.79137145607681D-02, 1.28379167095513D-01 /),  &
           b(3) = (/  3.01048631703895D-03,  5.38971687740286D-02, 3.75795757275549D-01 /),  &
           p(8) = (/ -1.36864857382717D-07,  5.64195517478974D-01, 7.21175825088309D+00, 4.31622272220567D+01, 1.52989285046940D+02, 3.39320816734344D+02, 4.51918953711873D+02, 3.00459261020162D+02 /), &
           q(8) = (/  1.00000000000000D+00,  1.27827273196294D+01, 7.70001529352295D+01, 2.77585444743988D+02, 6.38980264465631D+02, 9.31354094850610D+02, 7.90950925327898D+02, 3.00459260956983D+02 /), &
           r(5) = (/  2.10144126479064D+00,  2.62370141675169D+01, 2.13688200555087D+01, 4.65807828718470D+00, 2.82094791773523D-01 /),  &
           s(4) = (/  9.41537750555460D+01,  1.87114811799590D+02, 9.90191814623914D+01, 1.80124575948747D+01 /)
  REAL(8)                 ::  ax, bot, t, top, x2
  
  
  ax = ABS(x)
  
  IF (ax <= half) THEN
     t = x*x
     top = ((((a(1)*t + a(2))*t + a(3))*t + a(4))*t + a(5)) + one
     bot = ((b(1)*t + b(2))*t + b(3))*t + one
     fn_val = x*(top/bot)
     RETURN
  END IF
  
  IF (ax <= 4.0) THEN
     top = ((((((p(1)*ax + p(2))*ax + p(3))*ax + p(4))*ax + p(5))*ax  &
           + p(6))*ax + p(7))*ax + p(8)
     bot = ((((((q(1)*ax + q(2))*ax + q(3))*ax + q(4))*ax + q(5))*ax  &
           + q(6))*ax + q(7))*ax + q(8)
     fn_val = half + (half - EXP(-x*x)*top/bot)
     IF (x < zero) fn_val = -fn_val
     RETURN
  END IF
  
  IF (ax < 5.8) THEN
     x2 = x*x
     t = one / x2
     top = (((r(1)*t + r(2))*t + r(3))*t + r(4))*t + r(5)
     bot = (((s(1)*t + s(2))*t + s(3))*t + s(4))*t + one
     fn_val = (c - top/(x2*bot)) / ax
     fn_val = half + (half - EXP(-x2)*fn_val)
     IF (x < zero) fn_val = -fn_val
     RETURN
  END IF
  
  fn_val = SIGN(one, x)
  
  
  END FUNCTION erf
 
  
  SUBROUTINE init_hdf5()
  
  USE HDF5

  IMPLICIT NONE

  CALL h5open_f(herror)

  END SUBROUTINE init_hdf5


  SUBROUTINE finl_hdf5()
  
  USE HDF5

  IMPLICIT NONE

  CALL h5close_f(herror)

  END SUBROUTINE finl_hdf5



  SUBROUTINE init_mpi()
  
  IMPLICIT NONE

  CALL MPI_INIT(merror)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,merror)

  END SUBROUTINE init_mpi


  SUBROUTINE finl_mpi()
  
  IMPLICIT NONE

  CALL MPI_FINALIZE(merror)

  END SUBROUTINE finl_mpi

  SUBROUTINE print_fcomm_size()

  IMPLICIT NONE

  INTEGER :: fcomm_size
  INTEGER :: fcomm_rank
  INTEGER :: name_len
  CHARACTER(len=MPI_MAX_OBJECT_NAME) :: comm_name

  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,fcomm_size,merror)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,fcomm_rank,merror)
  CALL MPI_COMM_GET_NAME(MPI_COMM_WORLD,comm_name,name_len,merror)

  print*, 'MPI_COMM_WORLD with name ',comm_name,' in Fortran has size: ', fcomm_size, ' in process ',fcomm_rank,'.'

  END SUBROUTINE print_fcomm_size

  SUBROUTINE mpi_bcast_fort()

  IMPLICIT NONE

  CALL MPI_BCAST(finish_yes,1,MPI_LOGICAL,0,COMM_CART,merror)


  END SUBROUTINE mpi_bcast_fort


  SUBROUTINE check_node_ids

  IMPLICIT NONE

  INTEGER                ::  msize
  INTEGER                ::  i,j,k
  CHARACTER(LEN=2)       ::  schar1,schar2,schar3
  INTEGER                ::  ii,jj,kk,nid,gnid
  INTEGER                ::  nn1l, nn1u, nn2l, nn2u, nn3l, nn3u
  INTEGER                ::  dummy1, dummy2
  INTEGER                ::  nn1, nn2, nn3
  INTEGER                ::  dir
  LOGICAL                ::  boundary_yes, vel_grid_yes
  INTEGER                ::  n_local, n_ghost
  INTEGER                ::  n_elem_local

  ii = n1p
  jj = 1
  kk = 1
  dir = 1
  boundary_yes = .TRUE.
  vel_grid_yes = .TRUE.

  CALL get_block_dims(nn1, nn2, nn3, .FALSE., .FALSE., -1)

  write(*,*) 'process',rank,'nn1:',nn1,'N1p-S1p+1',N1p-S1p+1
  write(*,*) 'process',rank,'nn2:',nn2,'N2p-S2p+1',N2p-S2p+1
  write(*,*) 'process',rank,'nn3:',nn3,'N3p-S3p+1',N3p-S3p+1

  CALL get_block_dims(nn1, nn2, nn3, .TRUE., .FALSE., 1)

  write(*,*) 'process',rank,'nn1:',nn1,'N11-S11+1',N11-S11+1
  write(*,*) 'process',rank,'nn2:',nn2,'N21-S21+1',N21-S21+1
  write(*,*) 'process',rank,'nn3:',nn3,'N31-S31+1',N31-S31+1

  CALL get_block_dims(nn1, nn2, nn3, .TRUE., .TRUE., 1)

  write(*,*) 'process',rank,'nn1:',nn1,'N11B-S11B+1',N11B-S11B+1
  write(*,*) 'process',rank,'nn2:',nn2,'N21B-S21B+1',N21B-S21B+1
  write(*,*) 'process',rank,'nn3:',nn3,'N31B-S31B+1',N31B-S31B+1

  !--- Get the global number of nodes on which the values are stored
  CALL get_boundary_block_dims(nn1l, nn1u, nn2l, nn2u, nn3l, nn3u, vel_grid_yes, boundary_yes, dir)
  nn1 = nn1l + (NB1 - 2)*(N1 - 1) + nn1u
  nn2 = nn2l + (NB2 - 2)*(N2 - 1) + nn2u
  nn3 = nn3l + (NB3 - 2)*(N3 - 1) + nn3u

  !--- If in any direction only one block exists override dimension with local one
  IF (NB1 .EQ. 1) THEN
    CALL get_block_dims(nn1, dummy1, dummy2, vel_grid_yes, boundary_yes, dir)
  ENDIF
  IF (NB2 .EQ. 1) THEN
    CALL get_block_dims(dummy1, nn2, dummy2, vel_grid_yes, boundary_yes, dir)
  ENDIF
  IF (NB3 .EQ. 1) THEN
    CALL get_block_dims(dummy1, dummy2, nn3, vel_grid_yes, boundary_yes, dir)
  END IF

  

  CALL num_to_string(2,iB(1,1),schar1)
  CALL num_to_string(2,iB(2,1),schar2)
  CALL num_to_string(2,iB(3,1),schar3)
  OPEN(rank, FILE = "node_ids_proc_"//schar1//"_"//schar2//"_"//schar3//".txt",STATUS='UNKNOWN')
  OPEN(rank+100, FILE = "node_ids_local_cont_proc_"//schar1//"_"//schar2//"_"//schar3//".txt",STATUS='UNKNOWN')
  OPEN(rank+200, FILE = "local_node_ids_local_cont_proc_"//schar1//"_"//schar2//"_"//schar3//".txt",STATUS='UNKNOWN')


  IF (dimens .EQ. 3) THEN
    DO k = s3p-1, n3p
      DO j = s2p-1, n2p
        DO i = s1p-1, n1p

        !=== local coordinates with ghost cells ===
        CALL local_cart2id_loc_con(nid,i,j,k,n_local,n_ghost)
        write(rank+200,'(3i4,1a,1i8)') i,j,k, '-->', nid
        CALL local_id2cart_loc_con(nid,ii,jj,kk,n_local,n_ghost)
        write(rank+200,'(1i8,1a,3i4)') nid, '-->', ii,jj,kk
        write(rank+200,'(a)') '============================='
        FLUSH(rank+200)
        END DO
      END DO
    END DO
  ELSE
    DO k = s3p, n3p
      DO j = s2p-1, n2p
        DO i = s1p-1, n1p

        !=== local coordinates with ghost cells ===
        CALL local_cart2id_loc_con(nid,i,j,k,n_local,n_ghost)
        write(rank+200,'(3i4,1a,1i8)') i,j,k, '-->', nid
        CALL local_id2cart_loc_con(nid,ii,jj,kk,n_local,n_ghost)
        write(rank+200,'(1i8,1a,3i4)') nid, '-->', ii,jj,kk
        write(rank+200,'(a)') '============================='
        FLUSH(rank+200)
        END DO
      END DO
    END DO
  END IF



  DO k = s3p, n3p
    DO j = s2p, n2p
      DO i = s1p, n1p

      !=== global coordinates locally contiguous ===
      CALL global_cart2id_loc_con(nid,i,j,k)
      write(rank+100,'(3i4,1a,1i8)') i,j,k, '-->', nid
      CALL global_id2cart_loc_con(nid,ii,jj,kk)
      write(rank+100,'(1i8,1a,3i4)') nid, '-->', ii,jj,kk
      write(rank+100,'(a)') '============================='
      !=== output all nodes and their global ids ===
      CALL global_cart2id_glob_con(nid,i,j,k,.FALSE.,.FALSE.,dir)
      write(rank,'(3i4,1a,1i8)') i,j,k, '-->', nid
      !write(*,*) '=== Pressure Grid ==='
      !write(*,*) 'Process ',rank,': (',iB(1,1),',',iB(2,1),',',iB(3,1),'), Node (',i,',',j,',',k,') --> NodeID: ',nid
      CALL global_id2cart_glob_con(nid,ii,jj,kk,.FALSE.,.FALSE.,dir)
      write(rank,'(1i8,1a,3i4)') nid, '-->', ii,jj,kk
      write(rank,'(a)') '============================='
      !write(*,*) 'Process ',rank,': (',iB(1,1),',',iB(2,1),',',iB(3,1),'), NodeID: ',nid,' --> Node (',i,',',j,',',k,')' 
      !write(*,*) 'Node has ID: ', nid,' of a total ', nn1, 'x', nn2, 'x',nn3 , '=',nn1*nn2*nn3
      FLUSH(rank)
      FLUSH(rank+100)

      END DO
    END DO
  END DO

  CLOSE(rank)
  CLOSE(rank+100)
  CLOSE(rank+200)

  !=== Do Element node checking ============================
  IF (dimens .EQ. 3) THEN
    n_elem_local = total_n_local_tet_elements()
  ELSE
    n_elem_local = total_n_local_tri_elements()
  END IF
  OPEN(rank+300, FILE = "local_elem_node_ids_local_cont_proc_"//schar1//"_"//schar2//"_"//schar3//".txt",STATUS='UNKNOWN')
  DO i = 0, n_elem_local-1
    write(rank+300,'(1a,1i8)') "element",i
    write(rank+300,'(a)') "=========="
    IF (dimens .EQ. 3) THEN
      DO j = 0, 3
        CALL local_tet_element_nodes(i,j,nid)
        CALL local_id2cart_loc_con(nid,ii,jj,kk,dummy1,dummy2)
        write(rank+300,'(1a,1i4,1a,1i8,1a,3i4)') 'node#',j,', node_id',nid,'  -->',ii,jj,kk
        CALL local_to_global_node_id_loc_con_allperiodic2(nid, gnid)
        write(rank+300,'(1a,1i4,1a,1i8,1a,3i4)') 'node#',j,', node_id_loc',nid,' node_id_glob', gnid
      END DO
    ELSE
      DO j = 0, 2
        CALL local_tri_element_nodes(i,j,nid)
        CALL local_id2cart_loc_con(nid,ii,jj,kk,dummy1,dummy2)
        write(rank+300,'(1a,1i4,1a,1i8,1a,3i4)') 'node#',j,', node_id',nid,'  -->',ii,jj,kk
        CALL local_to_global_node_id_loc_con_allperiodic2(nid, gnid)
        write(rank+300,'(1a,1i4,1a,1i8,1a,3i4)') 'node#',j,', node_id_loc',nid,' node_id_glob', gnid
      END DO
    END IF
    write(rank+300,'(a)') "=========="
    FLUSH(rank+300)
  END DO

  CLOSE(rank+300)

 !=== Do local -> global checking =========================
 OPEN(rank+400, FILE = "local_to_global_node_ids_local_cont_proc_"//schar1//"_"//schar2//"_"//schar3//".txt",STATUS='UNKNOWN')
 IF (dimens .EQ. 3) THEN
   DO kk = s3p - 1, n3p
     DO jj = s2p - 1, n2p
       DO ii = s1p - 1, n1p
          CALL local_cart2id_loc_con(nid,ii,jj,kk,dummy1,dummy2)
          CALL local_to_global_node_id_loc_con_allperiodic2(nid,gnid)
          write(rank+400,'(1a,1i8,1a,3i4,1a,1i8)') 'local:',nid,'; i,j,k:',ii,jj,kk,'; global conversion:',gnid
          CALL global_to_local_node_id_loc_con_allperiodic2(nid,gnid)
          write(rank+400,'(1a,1i8,1a,3i4,1a,1i8)') 'global:',gnid,'; i,j,k:',ii,jj,kk,'; local conversion:',nid
          write(rank+400,'(a)') "==============================="
          FLUSH(rank+400)
        END DO
      END DO
    END DO
  ELSE
     DO kk = s3p, n3p
       DO jj = s2p - 1, n2p
         DO ii = s1p - 1, n1p
            CALL local_cart2id_loc_con(nid,ii,jj,kk,dummy1,dummy2)
            CALL local_to_global_node_id_loc_con_allperiodic2(nid,gnid)
            write(rank+400,'(1a,1i8,1a,3i4,1a,1i8)') 'local:',nid,'; i,j,k:',ii,jj,kk,'; global conversion:',gnid
            CALL global_to_local_node_id_loc_con_allperiodic2(nid,gnid)
            write(rank+400,'(1a,1i8,1a,3i4,1a,1i8)') 'global:',gnid,'; i,j,k:',ii,jj,kk,'; local conversion:',nid
            write(rank+400,'(a)') "==============================="
            FLUSH(rank+400)
          END DO
        END DO
      END DO
  END IF


  CLOSE(rank+400)

  print*, 'x1p(N1p)', x1p(N1p), 'y1p(N1p)',y1p(2*N1p)

  print* , mod(-43, 40), modulo(-43,40)

  END SUBROUTINE check_node_ids



  SUBROUTINE local_to_global_node_id_loc_con(loc_node_id, glob_node_id)

  IMPLICIT NONE

  INTEGER, INTENT(IN)   ::  loc_node_id
  INTEGER, INTENT(OUT)  ::  glob_node_id
  INTEGER               ::  ii,jj,kk
  INTEGER               ::  n_local, n_ghost
  INTEGER               ::  max_current_global_id, min_current_global_id
  INTEGER               ::  nn1l, nn1u, nn2l, nn2u, nn3l, nn3u
  INTEGER, DIMENSION(3) ::  nn_bulk_block
  INTEGER, DIMENSION(3) ::  nn_face_1l_block, nn_face_2l_block, nn_face_3l_block
  INTEGER, DIMENSION(3) ::  nn_face_1u_block, nn_face_2u_block, nn_face_3u_block
  INTEGER, DIMENSION(3) ::  nn_edge_1l2l_block, nn_edge_1u2l_block, nn_edge_1l2u_block, nn_edge_1u2u_block
  INTEGER, DIMENSION(3) ::  nn_edge_1l3l_block, nn_edge_1u3l_block, nn_edge_1l3u_block, nn_edge_1u3u_block
  INTEGER, DIMENSION(3) ::  nn_edge_2l3l_block, nn_edge_2u3l_block, nn_edge_2l3u_block, nn_edge_2u3u_block
  INTEGER, DIMENSION(3) ::  nn_corner_1l2l3l_block, nn_corner_1l2u3l_block, nn_corner_1l2l3u_block, nn_corner_1l2u3u_block
  INTEGER, DIMENSION(3) ::  nn_corner_1u2l3l_block, nn_corner_1u2u3l_block, nn_corner_1u2l3u_block, nn_corner_1u2u3u_block



  !--- Get the current process' maximum id ---
  CALL global_cart2id_loc_con(max_current_global_id,N1p,N2p,N3p)
  !--- Add 1 to have starting global node on next process in 1 dir ---
  max_current_global_id = max_current_global_id + 1

  !--- Get the current process' minumum id ---
  CALL global_cart2id_loc_con(min_current_global_id,S1p,S2p,S3p)
  !--- Substract one to get end id of previous process in 1dir ---
  min_current_global_id = min_current_global_id - 1


  !--- Get the number of nodes on the blocks that lie at the boundary ---
  CALL get_boundary_block_dims(nn1l, nn1u, nn2l, nn2u, nn3l, nn3u, .FALSE., .FALSE., -1)

  !--- compute and store the number of nodes in each direction for different
  !--- types of blocks ---
  !--- bulk blocks ---
  nn_bulk_block(1:3) = (/(N1 - 1), (N2 - 1), (N3 - 1)/)
  !--- face blocks ---
  nn_face_1l_block(1:3) = (/ nn1l   , (N2 - 1), (N3 - 1)/)
  nn_face_1u_block(1:3) = (/ nn1u   , (N2 - 1), (N3 - 1)/)
  nn_face_2l_block(1:3) = (/(N1 - 1),  nn2l   , (N3 - 1)/)
  nn_face_2u_block(1:3) = (/(N1 - 1),  nn2u   , (N3 - 1)/)
  nn_face_3l_block(1:3) = (/(N1 - 1), (N2 - 1),  nn3l   /)
  nn_face_3u_block(1:3) = (/(N1 - 1), (N2 - 1),  nn3u   /)
  !--- edge blocks ---
  nn_edge_1l2l_block(1:3) = (/ nn1l   ,  nn2l   , (N3 - 1)/)
  nn_edge_1u2l_block(1:3) = (/ nn1u   ,  nn2l   , (N3 - 1)/)
  nn_edge_1l2u_block(1:3) = (/ nn1l   ,  nn2u   , (N3 - 1)/)
  nn_edge_1u2u_block(1:3) = (/ nn1u   ,  nn2u   , (N3 - 1)/)
  nn_edge_1l3l_block(1:3) = (/ nn1l   , (N2 - 1),  nn3l   /)
  nn_edge_1u3l_block(1:3) = (/ nn1u   , (N2 - 1),  nn3l   /)
  nn_edge_1l3u_block(1:3) = (/ nn1l   , (N2 - 1),  nn3u   /)
  nn_edge_1u3u_block(1:3) = (/ nn1u   , (N2 - 1),  nn3u   /)
  nn_edge_2l3l_block(1:3) = (/(N1 - 1),  nn2l   ,  nn3l   /)
  nn_edge_2u3l_block(1:3) = (/(N1 - 1),  nn2u   ,  nn3l   /)
  nn_edge_2l3u_block(1:3) = (/(N1 - 1),  nn2l   ,  nn3u   /)
  nn_edge_2u3u_block(1:3) = (/(N1 - 1),  nn2u   ,  nn3u   /)
  !--- corner blocks ---
  nn_corner_1l2l3l_block(1:3) = (/nn1l , nn2l , nn3l/)
  nn_corner_1u2l3l_block(1:3) = (/nn1u , nn2l , nn3l/)
  nn_corner_1l2u3l_block(1:3) = (/nn1l , nn2u , nn3l/)
  nn_corner_1l2l3u_block(1:3) = (/nn1l , nn2l , nn3u/)
  nn_corner_1u2u3l_block(1:3) = (/nn1u , nn2u , nn3l/)
  nn_corner_1u2l3u_block(1:3) = (/nn1u , nn2l , nn3u/)
  nn_corner_1l2u3u_block(1:3) = (/nn1l , nn2u , nn3u/)
  nn_corner_1u2u3u_block(1:3) = (/nn1u , nn2u , nn3u/)

  !--- first convert locally to i,j,k notation ---
  CALL local_id2cart_loc_con(loc_node_id, ii, jj, kk, n_local, n_ghost)
  !--- Inner vertices (ghost vertices excluded) ---
  IF (ii .GE. S1p .AND. ii .LE. N1p) THEN
    IF (jj .GE. S2p .AND. jj .LE. N2p) THEN
      IF (kk .GE. S3p .AND. kk .LE. N3p) THEN
        CALL global_cart2id_loc_con(glob_node_id, ii, jj, kk)
      END IF
    END IF
  END IF
  !--- ghost vertices perpendicularly in front of a block face ---
  ! for each direction we need to check whether in that direction we have a
  ! neighbor that is in bulk, on face or periodic (inverted face)
  !--- face 1
  !--- top
  IF (ii .GT. N1p) THEN
    IF (jj .GE. S2p .AND. jj .LE. N2p) THEN
      IF (kk .GE. S3p .AND. kk .LE. N3p) THEN
      !========================================================================
        IF (iB(3,1) .GT. 1) THEN
          IF (iB(3,1) .LT. NB3) THEN ! 3-dir top bulk or face neighbor
            IF (iB(2,1) .GT. 1) THEN
              IF (iB(2,1) .LT. NB2) THEN ! 2-dir top bulk or face neighbor
                IF (iB(1,1) .LT. NB1-1) THEN ! 1-dir top bulk neighbor
                  ! These blocks have bulk neighbors in all directions
                  glob_node_id = (jj-1)*nn_bulk_block(1) &
                               + (kk-1)*nn_bulk_block(1) * nn_bulk_block(2) &
                               + max_current_global_id
                ELSE IF (iB(1,1) .LT. NB1) THEN ! 1-dir top face neighbor
                  ! These blocks have a face neighbor in 1 dir otherwise bulk
                  glob_node_id = (jj-1)*nn_face_1u_block(1) &
                               + (kk-1)*nn_face_1u_block(1) * nn_face_1u_block(2) &
                               + max_current_global_id
                ELSE ! 1-dir top periodic neighbor
                  ! These blocks have a periodic neighbor in 1 dir otherwise bulk
                  glob_node_id = (jj-1)*nn_face_1l_block(1) &
                               + (kk-1)*nn_face_1l_block(1) * nn_face_1l_block(2) &
                               + max_current_global_id &
                               - PRODUCT(nn_face_1u_block) &
                               - (NB1 - 2)*PRODUCT(nn_bulk_block) &
                               - PRODUCT(nn_face_1l_block)
                END IF
              ELSE ! 2-dir top periodic neighbor
                IF (iB(1,1) .LT. NB1-1) THEN ! 1-dir top bulk neighbor
                  ! These blocks are adjacent to 2 dir face and will have an
                  ! upper 2-face block neighbor in 1 dir
                  ! otherwise
                  glob_node_id = (jj-1)*nn_face_2u_block(1) &
                               + (kk-1)*nn_face_2u_block(1) * nn_face_2u_block(2) &
                               + max_current_global_id
                ELSE IF (iB(1,1) .LT. NB1) THEN ! 1-dir top face neighbor
                  ! These blocks are adjacent to 2 dir face and will have an upper 2,
                  ! upper 1 edge block neighbor
                  glob_node_id = (jj-1)*nn_edge_1u2u_block(1) &
                               + (kk-1)*nn_edge_1u2u_block(1) * nn_edge_1u2u_block(2) &
                               + max_current_global_id
                ELSE ! 1-dir top periodic neighbor
                  ! These blocks are adjacent to upper 1 upper 2 edge
                  glob_node_id = (jj-1)*nn_edge_1l2u_block(1) &
                               + (kk-1)*nn_edge_1l2u_block(1) * nn_edge_1l2u_block(2) &
                               + max_current_global_id &
                               - PRODUCT(nn_edge_1u2u_block) &
                               - (NB1 - 2)*PRODUCT(nn_face_2u_block) &
                               - PRODUCT(nn_edge_1l2u_block)
                END IF
              END IF
            ELSE ! iB(2,1) .EQ. 1
              IF (iB(1,1) .LT. NB1-1) THEN ! 1-dir top bulk neighbor
                ! These blocks are adjacent to the lower 2 face and will have
                ! a lower 2 face block as neighbor in positive 1 dir
                glob_node_id = (jj-1)*nn_face_2l_block(1) &
                             + (kk-1)*nn_face_2l_block(1) * nn_face_2l_block(2) &
                             + max_current_global_id
              ELSE IF (iB(1,1) .LT. NB1) THEN ! 1-dir top face neighbor
                ! These blocks are adjacent to the lower 2 face and will have
                ! a lower 2, upper 1 edge block neighbor
                glob_node_id = (jj-1)*nn_edge_1u2l_block(1) &
                             + (kk-1)*nn_edge_1u2l_block(1) * nn_edge_1u2l_block(2) &
                             + max_current_global_id
              ELSE ! 1-dir top periodic neighbor
                ! These blocks are adjacent to upper 1 lower 2 edge
                glob_node_id = (jj-1)*nn_edge_1l2l_block(1) &
                             + (kk-1)*nn_edge_1l2l_block(1) * nn_edge_1l2l_block(2) &
                             + max_current_global_id &
                             - PRODUCT(nn_edge_1u2l_block) &
                             - (NB1 - 2)*PRODUCT(nn_face_2l_block) &
                             - PRODUCT(nn_edge_1l2l_block)
              END IF
            END IF
          ELSE ! iB(3,1) .EQ. NB3
            IF (iB(2,1) .GT. 1) THEN
              IF (iB(2,1) .LT. NB2) THEN ! 2-dir top bulk or face neighbor
                IF (iB(1,1) .LT. NB1-1) THEN ! 1-dir top bulk neighbor
                  ! These blocks are adjacent to 3 face and will have an upper 3
                  ! face block neighbor in 1 direction
                  glob_node_id = (jj-1)*nn_face_3u_block(1) &
                               + (kk-1)*nn_face_3u_block(1) * nn_face_3u_block(2) &
                               + max_current_global_id
                ELSE IF (iB(1,1) .LT. NB1) THEN ! 1-dir top face neighbor
                  ! These blocks are adjacent to 3 face and will have an upper
                  ! 3, upper 1 edge neighbor
                  glob_node_id = (jj-1)*nn_edge_1u3u_block(1) &
                               + (kk-1)*nn_edge_1u3u_block(1) * nn_edge_1u3u_block(2) &
                               + max_current_global_id
                ELSE ! 1-dir top periodic neighbor
                  ! These blocks are adjacent to upper 3 upper 1 edge
                  glob_node_id = (jj-1)*nn_edge_1l3u_block(1) &
                               + (kk-1)*nn_edge_1l3u_block(1) * nn_edge_1l3u_block(2) &
                               + max_current_global_id &
                               - PRODUCT(nn_edge_1u3u_block) &
                               - (NB1 - 2)*PRODUCT(nn_face_3u_block) &
                               - PRODUCT(nn_edge_1l3u_block)
                END IF
              ELSE ! iB(2,1) .EQ. NB2
                IF (iB(1,1) .LT. NB1-1) THEN ! 1-dir top bulk neighbor
                  ! These blocks are adjacent to 2 and 3 face (upper) and will
                  ! have an upper 2, upper 3 edge neighbor
                  glob_node_id = (jj-1)*nn_edge_2u3u_block(1) &
                               + (kk-1)*nn_edge_2u3u_block(1) * nn_edge_2u3u_block(2) &
                               + max_current_global_id
                ELSE IF (iB(1,1) .LT. NB1) THEN ! 1-dir top face neighbor
                  ! These blocks are adjacent to 2 and 3 face (upper) and will
                  ! have an upper 2, upper 3, upper 1 corner neighbor
                  glob_node_id = (jj-1)*nn_corner_1u2u3u_block(1) &
                               + (kk-1)*nn_corner_1u2u3u_block(1) * nn_corner_1u2u3u_block(2) &
                               + max_current_global_id
                ELSE ! 1-dir top periodic neighbor
                  ! These blocks are adjacent to upper 1 upper 2 upper 3 corner
                  glob_node_id = (jj-1)*nn_corner_1l2u3u_block(1) &
                               + (kk-1)*nn_corner_1l2u3u_block(1) * nn_corner_1l2u3u_block(2) &
                               + max_current_global_id &
                               - PRODUCT(nn_corner_1l2u3u_block) &
                               - (NB1 - 2)*PRODUCT(nn_edge_2u3u_block) &
                               - PRODUCT(nn_corner_1u2u3u_block)
                END IF
              END IF
            ELSE ! iB(2,1) .EQ. 1
              IF (iB(1,1) .LT. NB1-1) THEN ! 1-dir top bulk neighbor
                ! These blocks are adjacent to lower 2 face and upper 3 face and
                ! will have an lower 2, upper 3 edge neighbor
                glob_node_id = (jj-1)*nn_edge_2l3u_block(1) &
                             + (kk-1)*nn_edge_2l3u_block(1) * nn_edge_2l3u_block(2) &
                             + max_current_global_id
              ELSE IF (iB(1,1) .LT. NB1) THEN ! 1-dir top face neighbor
                ! These blocks are adjacent to lower 2 face and upper 3 face and
                ! will have a lower 2, upper 3, upper 1 corner neighbor
                glob_node_id = (jj-1)*nn_corner_1u2l3u_block(1) &
                             + (kk-1)*nn_corner_1u2l3u_block(1) * nn_corner_1u2l3u_block(2) &
                             + max_current_global_id
              ELSE ! 1-dir top periodic neighbor
                ! These blocks are adjacent to upper 1 lower 2 upper 3 corner
                glob_node_id = (jj-1)*nn_corner_1l2l3u_block(1) &
                             + (kk-1)*nn_corner_1l2l3u_block(1) * nn_corner_1l2l3u_block(2) &
                             + max_current_global_id &
                             - PRODUCT(nn_corner_1u2l3u_block) &
                             - (NB1 - 2)*PRODUCT(nn_edge_2l3u_block) &
                             - PRODUCT(nn_corner_1l2l3u_block)
              END IF
            END IF
          END IF
        ELSE ! iB(3,1) .EQ. 1
          IF (iB(2,1) .GT. 1) THEN
            IF (iB(2,1) .LT. NB2) THEN ! 2-dir top bulk or face neighbor
              IF (iB(1,1) .LT. NB1-1) THEN ! 1-dir top bulk neighbor
                ! These blocks are adjacent to lower 3 face and will have a
                ! lower 3 face block as neighbor in 1 dir
                glob_node_id = (jj-1)*nn_face_3l_block(1) &
                             + (kk-1)*nn_face_3l_block(1) * nn_face_3l_block(2) &
                             + max_current_global_id
              ELSE IF (iB(1,1) .LT. NB1) THEN ! 1-dir top face neighbor
                ! These blocks are adjacent to lower 3 face and will have a
                ! lower 3, upper 1 edge neighor
                glob_node_id = (jj-1)*nn_edge_1u3l_block(1) &
                             + (kk-1)*nn_edge_1u3l_block(1) * nn_edge_1u3l_block(2) &
                             + max_current_global_id
              ELSE ! 1-dir top periodic neighbor
                ! These block are adjacent to lower 3 upper 1 edge
                glob_node_id = (jj-1)*nn_edge_1l3l_block(1) &
                             + (kk-1)*nn_edge_1l3l_block(1) * nn_edge_1l3l_block(2) &
                             + max_current_global_id &
                             - PRODUCT(nn_edge_1u3l_block) &
                             - (NB1 - 2)*PRODUCT(nn_face_3l_block) &
                             - PRODUCT(nn_edge_1l3l_block)
              END IF
            ELSE ! iB(2,1) .EQ. NB2
              IF (iB(1,1) .LT. NB1-1) THEN ! 1-dir top bulk neighbor
                ! These blocks are adjacent to upper 2 and lower 3 face and will
                ! have an upper 2, lower 3 edge neighbor
                glob_node_id = (jj-1)*nn_edge_2u3l_block(1) &
                             + (kk-1)*nn_edge_2u3l_block(1) * nn_edge_2u3l_block(2) &
                             + max_current_global_id
              ELSE IF (iB(1,1) .LT. NB1) THEN ! 1-dir top face neighbor
                ! These blocks are adjacent to upper 2 and lower 3 face and will
                ! have and upper 2, lower 3, upper 1 neighbor
                glob_node_id = (jj-1)*nn_corner_1u2u3l_block(1) &
                             + (kk-1)*nn_corner_1u2u3l_block(1) * nn_corner_1u2u3l_block(2) &
                             + max_current_global_id
              ELSE ! 1-dir top periodic neighbor
                ! These blocks are adjacent to upper 1 upper2 lower 3 corner
                glob_node_id = (jj-1)*nn_corner_1l2u3l_block(1) &
                             + (kk-1)*nn_corner_1l2u3l_block(1) * nn_corner_1l2u3l_block(2) &
                             + max_current_global_id &
                             - PRODUCT(nn_corner_1u2u3l_block) &
                             - (NB1 - 2)*PRODUCT(nn_edge_2u3l_block) &
                             - PRODUCT(nn_corner_1l2u3l_block)
              END IF
            END IF
          ELSE ! iB(2,1) .EQ. 1
            IF (iB(1,1) .LT. NB1-1) THEN ! 1-dir top bulk neighbor
              ! These blocks are adjacent to lower 2 andd lower 3 face and will
              ! have a lower2, lower 3 edge neighbor
              glob_node_id = (jj-1)*nn_edge_2l3l_block(1) &
                           + (kk-1)*nn_edge_2l3l_block(1) * nn_edge_2l3l_block(2) &
                           + max_current_global_id
            ELSE IF (iB(1,1) .LT. NB1) THEN ! 1-dir top face neighbor
              ! These blocks are adjacent to lower 2 and lower 3 face and will
              ! have a lower 2, lower 3 upper 1 corner neighbor
              glob_node_id = (jj-1)*nn_corner_1u2l3l_block(1) &
                           + (kk-1)*nn_corner_1u2l3l_block(1) * nn_corner_1u2l3l_block(2) &
                           + max_current_global_id
            ELSE ! 1-dir top periodic neighbor
              ! These blocks are adjacent to upper 1 lower 2 lower 3 corner
              glob_node_id = (jj-1)*nn_corner_1l2l3l_block(1) &
                           + (kk-1)*nn_corner_1l2l3l_block(1) * nn_corner_1l2l3l_block(2) &
                           + max_current_global_id &
                           - PRODUCT(nn_corner_1u2l3l_block) &
                           - (NB1 - 2)*PRODUCT(nn_edge_2l3l_block) &
                           - PRODUCT(nn_corner_1l2l3l_block)
            END IF
          END IF
        END IF
      !========================================================================
      END IF
    END IF
  END IF
  !--- bottom
  IF (ii .LT. S1p) THEN
    IF (jj .GE. S2p .AND. jj .LE. N2p) THEN
      IF (kk .GE. S3p .AND. kk .LE. N3p) THEN
        IF (iB(3,1) .GT. 1) THEN
          IF (iB(3,1) .LT. NB3) THEN
            IF (iB(2,1) .GT. 1) THEN
              IF (iB(2,1) .LT. NB2) THEN
                IF (iB(1,1) .GT. 2) THEN
                  ! These blocks are bulk blocks with a bulk block in neg. 1dir
                  glob_node_id = min_current_global_id &
                               - (nn_bulk_block(2) - jj) * (nn_bulk_block(1)) &
                               - (nn_bulk_block(3) - kk) * (nn_bulk_block(1)) &
                                                         * (nn_bulk_block(2))
                ELSE IF (iB(1,1) .GT. 1) THEN
                  ! These blocks are bulk blocks with a lower 1 face block in
                  ! neg. 1 dir
                  glob_node_id = min_current_global_id &
                               - (nn_face_1l_block(2) - jj) * (nn_face_1l_block(1)) &
                               - (nn_face_1l_block(3) - kk) * (nn_face_1l_block(1)) &
                                                            * (nn_face_1l_block(2))
                ELSE ! iB(1,1) .EQ. 1
                  glob_node_id = min_current_global_id &
                               - (nn_face_1u_block(2) - jj) * (nn_face_1u_block(1)) &
                               - (nn_face_1u_block(3) - kk) * (nn_face_1u_block(1)) &
                                                            * (nn_face_1u_block(2)) &
                               + PRODUCT(nn_face_1l_block) &
                               + (NB1 - 2)*PRODUCT(nn_bulk_block) &
                               + PRODUCT(nn_face_1u_block)
                END IF
              ELSE ! iB(2,1) .EQ. NB2
                IF (iB(1,1) .GT. 2) THEN
                  ! These blocks are adjacent to upper 2 face with upper 2 face
                  ! block in neg. 1 dir
                  glob_node_id = min_current_global_id &
                               - (nn_face_2u_block(2) - jj) * (nn_face_2u_block(1)) &
                               - (nn_face_2u_block(3) - kk) * (nn_face_2u_block(1)) &
                                                            * (nn_face_2u_block(2))
                ELSE IF (iB(1,1) .GT. 1) THEN
                  ! These blocks are adjacent to upper 2 face with upper 2 lower
                  ! 1 edge neighbor
                  glob_node_id = min_current_global_id &
                               - (nn_edge_1l2u_block(2) - jj) * (nn_edge_1l2u_block(1)) &
                               - (nn_edge_1l2u_block(3) - kk) * (nn_edge_1l2u_block(1)) &
                                                              * (nn_edge_1l2u_block(2))
                ELSE ! iB(1,1) .EQ. 1
                  glob_node_id = min_current_global_id &
                               - (nn_edge_1u2u_block(2) - jj) * (nn_edge_1u2u_block(1)) &
                               - (nn_edge_1u2u_block(3) - kk) * (nn_edge_1u2u_block(1)) &
                                                              * (nn_edge_1u2u_block(2)) &
                               + PRODUCT(nn_edge_1l2u_block) &
                               + (NB1 -2)*PRODUCT(nn_face_2u_block) &
                               + PRODUCT(nn_edge_1u2u_block)
                END IF
              END IF
            ELSE ! iB(2,1) .EQ. 1
              IF (iB(1,1) .GT. 2) THEN
                ! These blocks are adjacent to lower 2 face with lower 2 face
                ! neighbor
                glob_node_id = min_current_global_id &
                             - (nn_face_2l_block(2) - jj) * (nn_face_2l_block(1)) &
                             - (nn_face_2l_block(3) - kk) * (nn_face_2l_block(1)) &
                                                          * (nn_face_2l_block(2))
              ELSE IF (iB(1,1) .GT. 1) THEN
                ! These blocks are adjacent to lower 2 face with lower 2 lower 1
                ! edge neighbor
                glob_node_id = min_current_global_id &
                             - (nn_edge_1l2l_block(2) - jj) * (nn_edge_1l2l_block(1)) &
                             - (nn_edge_1l2l_block(3) - kk) * (nn_edge_1l2l_block(1)) &
                                                            * (nn_edge_1l2l_block(2))
              ELSE ! iB(1,1) .EQ. 1
                glob_node_id = min_current_global_id &
                             - (nn_edge_1u2l_block(2) - jj) * (nn_edge_1u2l_block(1)) &
                             - (nn_edge_1u2l_block(3) - kk) * (nn_edge_1u2l_block(1)) &
                                                            * (nn_edge_1u2l_block(2)) &
                             + PRODUCT(nn_edge_1l2l_block) &
                             + (NB1 -2)*PRODUCT(nn_face_2l_block) &
                             + PRODUCT(nn_edge_1u2l_block)
              END IF
            END IF
          ELSE ! iB(3,1) .EQ. NB3
            IF (iB(2,1) .GT. 1) THEN
              IF (iB(2,1) .LT. NB2) THEN
                IF (iB(1,1) .GT. 2) THEN
                  ! These blocks are adjacent to upper 3 face with upper 3
                  ! face neighbor
                  glob_node_id = min_current_global_id &
                               - (nn_face_3u_block(2) - jj) * (nn_face_3u_block(1)) &
                               - (nn_face_3u_block(3) - kk) * (nn_face_3u_block(1)) &
                                                            * (nn_face_3u_block(2))
                ELSE IF (iB(1,1) .GT. 1) THEN
                  ! These blocks are adjacent to upper 3 face with upper 3 lower
                  ! 1 edge neighbor
                  glob_node_id = min_current_global_id &
                               - (nn_edge_1l3u_block(2) - jj) * (nn_edge_1l3u_block(1)) &
                               - (nn_edge_1l3u_block(3) - kk) * (nn_edge_1l3u_block(1)) &
                                                              * (nn_edge_1l3u_block(2))
                ELSE ! iB(1,1) .EQ. 1
                  glob_node_id = min_current_global_id &
                               - (nn_edge_1u3u_block(2) - jj) * (nn_edge_1u3u_block(1)) &
                               - (nn_edge_1u3u_block(3) - kk) * (nn_edge_1u3u_block(1)) &
                                                              * (nn_edge_1u3u_block(2)) &
                               + PRODUCT(nn_edge_1l3u_block) &
                               + (NB1 -2)*PRODUCT(nn_face_3u_block) &
                               + PRODUCT(nn_edge_1u3u_block)

                END IF
              ELSE ! iB(2,1) .EQ. NB2
                IF (iB(1,1) .GT. 2) THEN
                  ! These blocks are adjacent to upper 2 upper 3 edge with
                  ! corresponding neighbor
                  glob_node_id = min_current_global_id &
                               - (nn_edge_2u3u_block(2) - jj) * (nn_edge_2u3u_block(1)) &
                               - (nn_edge_2u3u_block(3) - kk) * (nn_edge_2u3u_block(1)) &
                                                              * (nn_edge_2u3u_block(2))
                ELSE IF (iB(1,1) .GT. 1) THEN
                  ! These blocks are adjacent to upper2 upper 3 edge with upper
                  ! 2 upper 3 lower 1 corner neighbor
                  glob_node_id = min_current_global_id &
                               - (nn_corner_1l2u3u_block(2) - jj) * (nn_corner_1l2u3u_block(1)) &
                               - (nn_corner_1l2u3u_block(3) - kk) * (nn_corner_1l2u3u_block(1)) &
                                                                  * (nn_corner_1l2u3u_block(2))
                ELSE ! iB(1,1) .EQ. 1
                  glob_node_id = min_current_global_id &
                               - (nn_corner_1u2u3u_block(2) - jj) * (nn_corner_1u2u3u_block(1)) &
                               - (nn_corner_1u2u3u_block(3) - kk) * (nn_corner_1u2u3u_block(1)) &
                                                                  * (nn_corner_1u2u3u_block(2)) &
                               + PRODUCT(nn_corner_1l2u3u_block) &
                               + (NB1 - 2)*PRODUCT(nn_edge_2u3u_block) &
                               + PRODUCT(nn_corner_1u2u3u_block)
                END IF
              END IF
            ELSE ! iB(2,1) .EQ. 1
              IF (iB(1,1) .GT. 2) THEN
                ! These blocks are adjacent to lower 2 upper 3 edge with
                ! corresponding neighbor
                glob_node_id = min_current_global_id &
                             - (nn_edge_2l3u_block(2) - jj) * (nn_edge_2l3u_block(1)) &
                             - (nn_edge_2l3u_block(3) - kk) * (nn_edge_2l3u_block(1)) &
                                                            * (nn_edge_2l3u_block(2))
              ELSE IF (iB(1,1) .GT. 1) THEN
                ! These blocks are adjacent to lower 2 upper 3 edge with lower
                ! 2 upper 3 lower 1 corner neighbor
                glob_node_id = min_current_global_id &
                             - (nn_corner_1l2l3u_block(2) - jj) * (nn_corner_1l2l3u_block(1)) &
                             - (nn_corner_1l2l3u_block(3) - kk) * (nn_corner_1l2l3u_block(1)) &
                                                                * (nn_corner_1l2l3u_block(2))
              ELSE ! iB(1,1) .EQ. 1
                glob_node_id = min_current_global_id &
                             - (nn_corner_1u2l3u_block(2) - jj) * (nn_corner_1u2l3u_block(1)) &
                             - (nn_corner_1u2l3u_block(3) - kk) * (nn_corner_1u2l3u_block(1)) &
                                                                * (nn_corner_1u2l3u_block(2)) &
                             + PRODUCT(nn_corner_1l2l3u_block) &
                             + (NB1 - 2)*PRODUCT(nn_edge_2l3u_block) &
                             + PRODUCT(nn_corner_1u2l3u_block)
              END IF
            END IF
          END IF
        ELSE ! iB(3,1) .EQ. 1
          IF (iB(2,1) .GT. 1) THEN
            IF (iB(2,1) .LT. NB2) THEN
              IF (iB(1,1) .GT. 2) THEN
                ! These blocks are adjacent to lower 3 face with corresponding
                ! neighbor
                glob_node_id = min_current_global_id &
                             - (nn_face_3l_block(2) - jj) * (nn_face_3l_block(1)) &
                             - (nn_face_3l_block(3) - kk) * (nn_face_3l_block(1)) &
                                                          * (nn_face_3l_block(2))
              ELSE IF (iB(1,1) .GT. 1) THEN
                ! These blocks are adjacent to lower 3 face with lower 3 lower
                ! 1 edge neighbor
                glob_node_id = min_current_global_id &
                             - (nn_edge_1l3l_block(2) - jj) * (nn_edge_1l3l_block(1)) &
                             - (nn_edge_1l3l_block(3) - kk) * (nn_edge_1l3l_block(1)) &
                                                            * (nn_edge_1l3l_block(2))
              ELSE ! iB(1,1) .EQ. 1
                glob_node_id = min_current_global_id &
                             - (nn_edge_1u3l_block(2) - jj) * (nn_edge_1u3l_block(1)) &
                             - (nn_edge_1u3l_block(3) - kk) * (nn_edge_1u3l_block(1)) &
                                                            * (nn_edge_1u3l_block(2)) &
                             + PRODUCT(nn_edge_1l3l_block) &
                             + (NB1 -2)*PRODUCT(nn_face_3l_block) &
                             + PRODUCT(nn_edge_1u3l_block)
              END IF
            ELSE ! iB(2,1) .EQ. NB2
              IF (iB(1,1) .GT. 2) THEN
                ! These blocks are adjacent to lower 3 upper 2 edge with
                ! corresponding neighbor
                glob_node_id = min_current_global_id &
                             - (nn_edge_2u3l_block(2) - jj) * (nn_edge_2u3l_block(1)) &
                             - (nn_edge_2u3l_block(3) - kk) * (nn_edge_2u3l_block(1)) &
                                                            * (nn_edge_2u3l_block(2))
              ELSE IF (iB(1,1) .GT. 1) THEN
                ! These locks are adjacent to lower 3 upper 2 edge with lower 3
                ! upper 2 lower 1 corner neighbor
                glob_node_id = min_current_global_id &
                             - (nn_corner_1l2u3l_block(2) - jj) * (nn_corner_1l2u3l_block(1)) &
                             - (nn_corner_1l2u3l_block(3) - kk) * (nn_corner_1l2u3l_block(1)) &
                                                                * (nn_corner_1l2u3l_block(2))
              ELSE ! iB(1,1) .EQ. 1
                glob_node_id = min_current_global_id &
                             - (nn_corner_1u2u3l_block(2) - jj) * (nn_corner_1u2u3l_block(1)) &
                             - (nn_corner_1u2u3l_block(3) - kk) * (nn_corner_1u2u3l_block(1)) &
                                                                * (nn_corner_1u2u3l_block(2)) &
                             + PRODUCT(nn_corner_1l2u3l_block) &
                             + (NB1 - 2)*PRODUCT(nn_edge_2u3l_block) &
                             + PRODUCT(nn_corner_1u2u3l_block)
              END IF
            END IF
          ELSE ! iB(2,1) .EQ. 1
            IF (iB(1,1) .GT. 2) THEN
              ! These blocks are adjacent to lower 3 lower 2 edge with
              ! corresponding neighbor
              glob_node_id = min_current_global_id &
                           - (nn_edge_2l3l_block(2) - jj) * (nn_edge_2l3l_block(1)) &
                           - (nn_edge_2l3l_block(3) - kk) * (nn_edge_2l3l_block(1)) &
                                                          * (nn_edge_2l3l_block(2))
            ELSE IF (iB(1,1) .GT. 1) THEN
              ! These locks are adjacent to lower 3 lower 2 edge with lower 3
              ! lower 2 lower 1 corner neighbor
              glob_node_id = min_current_global_id &
                           - (nn_corner_1l2l3l_block(2) - jj) * (nn_corner_1l2l3l_block(1)) &
                           - (nn_corner_1l2l3l_block(3) - kk) * (nn_corner_1l2l3l_block(1)) &
                                                              * (nn_corner_1l2l3l_block(2))
            ELSE ! iB(1,1) .EQ. 1
              glob_node_id = min_current_global_id &
                           - (nn_corner_1u2l3l_block(2) - jj) * (nn_corner_1u2l3l_block(1)) &
                           - (nn_corner_1u2l3l_block(3) - kk) * (nn_corner_1u2l3l_block(1)) &
                                                              * (nn_corner_1u2l3l_block(2)) &
                           + PRODUCT(nn_corner_1l2l3l_block) &
                           + (NB1 - 2)*PRODUCT(nn_edge_2l3l_block) &
                           + PRODUCT(nn_corner_1u2l3l_block)
            END IF
          END IF
        END IF
      END IF
    END IF
  END IF
  !--- face 2
  !--- top
  IF (ii .GE. S1p .AND. ii .LE. N1p) THEN
    IF (jj .GT. N2p) THEN
      IF (kk .GE. S3p .AND. kk .LE. N3p) THEN

      END IF
    END IF
  END IF
  !--- bottom
  IF (ii .GE. S1p .AND. ii .LE. N1p) THEN
    IF (jj .LT. S2p) THEN
      IF (kk .GE. S3p .AND. kk .LE. N3p) THEN

      END IF
    END IF
  END IF
  !--- face 3
  !--- top
  IF (ii .GE. S1p .AND. ii .LE. N1p) THEN
    IF (jj .GE. S2p .AND. jj .LE. N2p) THEN
      IF (kk .GT. N3p) THEN

      END IF
    END IF
  END IF
  !--- bottom
  IF (ii .GE. S1p .AND. ii .LE. N1p) THEN
    IF (jj .LT. S2p) THEN
      IF (kk .GE. S3p .AND. kk .LE. N3p) THEN

      END IF
    END IF
  END IF
  !--- ghost vertices diagonally off to edges ---

  !--- ghost vertices in corners ---

  END SUBROUTINE local_to_global_node_id_loc_con


  SUBROUTINE local_to_global_node_id_loc_con_allperiodic(loc_node_id, glob_node_id)

  IMPLICIT NONE

  INTEGER, INTENT(IN)    ::  loc_node_id
  INTEGER, INTENT(OUT)   ::  glob_node_id
  INTEGER                ::  nn1, nn2, nn3
  INTEGER                ::  nn_block
  INTEGER                ::  ii, jj, kk
  INTEGER                ::  max_current_global_id, min_current_global_id
  INTEGER                ::  n_local, n_ghost


  nn1 = N1 - 1
  nn2 = N2 - 1
  nn3 = N3 - 1

  nn_block = nn1 * nn2 * nn3

  !--- Get the current process' maximum id ---
  CALL global_cart2id_loc_con(max_current_global_id,N1p,N2p,N3p)

  !--- Get the current process' minumum id ---
  CALL global_cart2id_loc_con(min_current_global_id,S1p,S2p,S3p)
  
  !--- first convert locally to i,j,k notation ---
  CALL local_id2cart_loc_con(loc_node_id, ii, jj, kk, n_local, n_ghost)
  !--- Inner vertices (ghost vertices excluded) ---
  IF (ii .GE. S1p .AND. ii .LE. N1p) THEN
    IF (jj .GE. S2p .AND. jj .LE. N2p) THEN
      IF (kk .GE. S3p .AND. kk .LE. N3p) THEN
        CALL global_cart2id_loc_con(glob_node_id, ii, jj, kk)
      END IF
    END IF
  END IF
  !--- ghost vertices perpendicularly in front of a block face ---
  ! for each direction we need to check whether in that direction we have a
  ! neighbor that is in bulk, on face or periodic (inverted face)
  !--- face 1
  !--- top
  IF (ii .GT. N1p) THEN
    IF (jj .GE. S2p .AND. jj .LE. N2p) THEN
      IF (kk .GE. S3p .AND. kk .LE. N3p) THEN
        IF (iB(1,1) .LT. NB1) THEN
          glob_node_id = (jj - 1) * nn1 &
                       + (kk - 1) * nn1 * nn2 &
                       + max_current_global_id + 1
        ELSE ! iB(1,1) .EQ. NB1
          glob_node_id = (jj - 1) * nn1 &
                       + (kk - 1) * nn1 * nn2 &
                       + max_current_global_id + 1&
                       - NB1 * nn_block
        END IF
      END IF
    END IF
  END IF
  !--- bottom
  IF (ii .LT. S1p) THEN
    IF (jj .GE. S2p .AND. jj .LE. N2p) THEN
      IF (kk .GE. S3p .AND. kk .LE. N3p) THEN
        IF (iB(1,1) .GT. 1) THEN
          glob_node_id = min_current_global_id - 1 &
                       - (nn2 - jj) * nn1 &
                       - (nn3 - kk) * nn1 * nn2
        ELSE ! iB(1,1) .EQ. 1
          glob_node_id = min_current_global_id - 1 &
                       - (nn2 - jj) * nn1 &
                       - (nn3 - kk) * nn1 * nn2 &
                       + NB1 * nn_block
        END IF
      END IF
    END IF
  END IF
  !--- face 2
  !--- top
  IF (ii .GE. S1p .AND. ii .LE. N1p) THEN
    IF (jj .GT. N2p) THEN
      IF (kk .GE. S3p .AND. kk .LE. N3p) THEN
        IF (iB(2,1) .LT. NB2) THEN
          glob_node_id = min_current_global_id &
                       + (ii - 1) &
                       + (kk - 1) * nn1 * nn2 &
                       + NB1 * nn_block
        ELSE ! iB(2,1) .EQ. NB2
          glob_node_id = min_current_global_id &
                       + (ii - 1) &
                       + (kk - 1) * nn1 * nn2 &
                       + NB1 * nn_block &
                       - NB1 * NB2 * nn_block
        END IF
      END IF
    END IF
  END IF
  !--- bottom
  IF (ii .GE. S1p .AND. ii .LE. N1p) THEN
    IF (jj .LT. S2p) THEN
      IF (kk .GE. S3p .AND. kk .LE. N3p) THEN
        IF (iB(2,1) .GT. 1) THEN
          glob_node_id = max_current_global_id &
                       - (nn1 - ii) &
                       - (nn3 - kk) * nn1 * nn2 &
                       - NB1 * nn_block
        ELSE ! iB(2,1) .EQ. 1
          glob_node_id = max_current_global_id &
                       - (nn1 - ii) &
                       - (nn3 - kk) * nn1 * nn2 &
                       - NB1 * nn_block &
                       + NB1 * NB2 * nn_block
        END IF
      END IF
    END IF
  END IF
  !--- face 3
  !--- top
  IF (ii .GE. S1p .AND. ii .LE. N1p) THEN
    IF (jj .GE. S2p .AND. jj .LE. N2p) THEN
      IF (kk .GT. N3p) THEN
        IF (iB(3,1) .LT. NB3) THEN
          glob_node_id = min_current_global_id &
                       + (ii - 1) &
                       + (jj - 1) * nn1 &
                       + NB1 * NB2 * nn_block
        ELSE ! iB(3,1) .EQ. NB3
          glob_node_id = min_current_global_id &
                       + (ii - 1) &
                       + (jj - 1) * nn1 &
                       + NB1 * NB2 * nn_block &
                       - NB1 * NB2 * NB3 * nn_block
        END IF
      END IF
    END IF
  END IF
  !--- bottom
  IF (ii .GE. S1p .AND. ii .LE. N1p) THEN
    IF (jj .GE. S2p .AND. jj .LE. N2p) THEN
      IF (kk .LT. S3p) THEN
        IF (iB(3,1) .GT. 1) THEN
          glob_node_id = max_current_global_id &
                       - (nn1 - ii) &
                       - (nn2 - jj) * nn1 &
                       - NB1 * NB2 * nn_block
        ELSE ! iB(3,1) .EQ. 1
          glob_node_id = max_current_global_id &
                       - (nn1 - ii) &
                       - (nn2 - jj) * nn1 &
                       - NB1 * NB2 * nn_block &
                       + NB1 * NB2 * NB3 * nn_block
        END IF
      END IF
    END IF
  END IF

  !--- ghost vertices diagonally off to edges ---
  !--- 1L, 2L
  IF (ii .LT. S1p) THEN
    IF (jj .LT. S2p) THEN
      IF (kk .GE. S3p .AND. kk .LE. N3p) THEN
        IF (iB(1,1) .EQ. 1 .AND. iB(2,1) .EQ. 1) THEN
          ! These nodes go dianonaaly to 1U, 2U
          glob_node_id = min_current_global_id - 1 &
                       - (nn3 - kk) * nn1 * nn2 &
                       + NB1 * NB2 * nn_block
        ELSE IF (iB(2,1) .EQ. 1) THEN
          ! These nodes go to top 2U face
          glob_node_id = min_current_global_id - 1 &
                       - (nn3 - kk) * nn1 * nn2 &
                       + NB1 * (NB2 - 1) * nn_block
        ELSE IF (iB(1,1) .EQ. 1) THEN
          ! These node go top 1U face
          glob_node_id = max_current_global_id &
                       - (nn3 - kk) * nn1 * nn2 &
                       - nn_block
        ELSE ! iB(1:2,1) .GT. 1
          glob_node_id = min_current_global_id - 1 &
                       - (nn3 - kk) * nn1 * nn2 &
                       - NB1 * nn_block
        END IF
      END IF
    END IF
  END IF
  !--- 1U, 2L
  IF (ii .GT. N1p) THEN
    IF (jj .LT. S2p) THEN
      IF (kk .GE. S3p .AND. kk .LE. N3p) THEN
        IF (iB(1,1) .EQ. NB1 .AND. iB(2,1) .EQ. 1) THEN
          ! These nodes go diagonally to 1L, 2U
          glob_node_id = max_current_global_id + 1&
                       - (nn3 - kk) * nn1 * nn2 &
                       - nn1 &
                       + nn_block
        ELSE IF (iB(2,1) .EQ. 1) THEN
          ! These nodes go to top 2U face
          glob_node_id = max_current_global_id + 1 &
                       - (nn3 - kk) * nn1 * nn2 &
                       - nn1 &
                       + ((NB2 - 1) *NB1 + 1) * nn_block
        ELSE IF (iB(1,1) .EQ. NB1) THEN
          ! These nodes go to bottom 1L face
          glob_node_id = max_current_global_id + 1&
                       - (nn3 - kk) * nn1 * nn2 &
                       - nn1 &
                       - (2 * NB1 - 1) * nn_block
        ELSE
          ! These nodes do not have periodicity
          glob_node_id = max_current_global_id + 1&
                       - (nn3 - kk) * nn1 * nn2 &
                       - nn1 &
                       - (NB1 - 1) * nn_block
        END IF
      END IF
    END IF
  END IF
  !--- 1L, 2U
  IF (ii .LT. S1p) THEN
    IF (jj .GT. N2p) THEN
      IF (kk .GE. S3p .AND. kk .LE. N3p) THEN
        IF (iB(1,1) .EQ. 1 .AND. iB(2,1) .EQ. NB2) THEN
         ! These nodes go diagonally to 1U, 2L
          glob_node_id = min_current_global_id - 1 &
                       + nn1 &
                       + (kk - 1) * nn1 * nn2 &
                       - nn_block &
                       - (NB2 - 2) * NB1 *nn_block
        ELSE IF (iB(2,1) .EQ. NB2) THEN
          ! These nodes go to bottom 2L face
          glob_node_id = min_current_global_id - 1&
                       + nn1 &
                       + (kk - 1) * nn1 * nn2 &
                       - ((NB2 - 1) * NB1 + 1) * nn_block 
        ELSE IF (iB(1,1) .EQ. 1) THEN
          ! These nodes go to top 1U face
          glob_node_id = min_current_global_id - 1 &
                       + nn1 &
                       + (kk - 1) * nn1 * nn2 &
                       + (2*NB1 - 1) * nn_block
        ELSE
          ! These nodes have no periodicity
          glob_node_id = min_current_global_id - 1 &
                       + nn1 &
                       + (kk - 1) * nn1 * nn2 &
                       + (NB1 - 1) * nn_block
        END IF
      END IF
    END IF
  END IF
  !--- 1U, 2U
  IF (ii .GT. N1p) THEN
    IF (jj .GT. N2p) THEN
      IF (kk .GE. S3p .AND. kk .LE. N3p) THEN
        IF (iB(1,1) .EQ. NB1 .AND. iB(2,1) .EQ. NB2) THEN
          ! These nodes go diagonaly to 1L, 2L
          glob_node_id = min_current_global_id &
                       + (kk - 1) * nn1 * nn2 &
                       - (NB2 * NB1 - 1) * nn_block
        ELSE IF (iB(2,1) .EQ. NB2) THEN
          ! These nodes go to bottom 2L face
          glob_node_id = min_current_global_id &
                       + (kk - 1) * nn1 * nn2 &
                       - ((NB2 - 1) * NB1 - 1) * nn_block
        ELSE IF (iB(1,1) .EQ. NB1) THEN
          ! These nodes go to bottom 1L face
          glob_node_id = min_current_global_id &
                       + (kk - 1) * nn1 * nn2 &
                       + nn_block
        ELSE
          ! These nodes have no periodicity
          glob_node_id = min_current_global_id &
                       + (kk - 1) * nn1 * nn2 &
                       + (NB1 + 1) * nn_block 
        END IF
      END IF
    END IF
  END IF
  !--- 1L, 3L
  IF (ii .LT. S1p) THEN
    IF (jj .GE. S2p .AND. jj .LE. N2p) THEN
      IF (kk .LT. S3p) THEN
        IF (iB(1,1) .EQ. 1 .AND. iB(3,1) .EQ. 1) THEN
          ! These nodes go diagonally to 1U 3U
          glob_node_id = max_current_global_id &
                       - (nn2 - jj) * nn1 &
                       + (NB3 - 1) * NB2 * NB1 * nn_block &
                       + (NB1 - 1) * nn_block
        ELSE IF (iB(3,1) .EQ. 1) THEN
          ! These nodes go top 3U face
          glob_node_id = max_current_global_id &
                       - (nn2 - jj) * nn1 &
                       + ((NB3 - 1) * NB2 * NB1 - 1) * nn_block
        ELSE IF (iB(1,1) .EQ. 1) THEN
          ! These nodes go to top 1U face
          glob_node_id = max_current_global_id &
                       - (nn2 - jj) * nn1 &
                       - ((NB2 - 1) * NB1 + 1) * nn_block
        ELSE
          ! These nodes have no periodicity
          glob_node_id = max_current_global_id &
                       - (nn2 - jj) * nn1 &
                       - (NB1 * NB2 + 1) * nn_block
        END IF
      END IF
    END IF
  END IF
  !--- 1U, 3L
  IF (ii .GT. N1p) THEN
    IF (jj .GE. S2p .AND. jj .LE. N2p) THEN
      IF (kk .LT. S3p) THEN
        IF (iB(1,1) .EQ. NB1 .AND. iB(3,1) .EQ. 1) THEN
          ! These nodes go diagonally to 1L 3U
          glob_node_id = min_current_global_id  &
                       + nn1 * nn2 * (nn3 - 1) &
                       + (jj - 1)*nn1 &
                       + (NB3 - 2) * NB1 * NB2 * nn_block &
                       + (NB2 - 1) * NB1 * nn_block &
                       + nn_block
        ELSE IF (iB(3,1) .EQ. 1) THEN
          ! These nodes go to top 3U face
          glob_node_id = min_current_global_id  &
                       + nn1 * nn2 * (nn3 - 1) &
                       + (jj - 1)*nn1 &
                       + (NB3 - 1) * NB2 * NB1 * nn_block &
                       + nn_block
        ELSE IF (iB(1,1) .EQ. NB1) THEN
          ! These nodes go to bottom 1L face
          glob_node_id = min_current_global_id &
                       + nn1 * nn2 * (nn3 - 1) &
                       + (jj - 1)*nn1 &
                       - NB2 * NB1 * nn_block &
                       - (NB1 - 1) * nn_block
        ELSE
          ! These nodes do not have periodicity
          glob_node_id = min_current_global_id &
                       + nn1 * nn2 * (nn3 - 1) &
                       + (jj - 1)*nn1 &
                       - (NB2 - 1) * NB1 * nn_block &
                       - (NB1 - 1) * nn_block
        END IF
      END IF
    END IF
  END IF
  !--- 1L, 3U
  IF (ii .LT. S1p) THEN
    IF (jj .GE. S2p .AND. jj .LE. N2p) THEN
      IF (kk .GT. N3p) THEN
        IF (iB(1,1) .EQ. 1 .AND. iB(3,1) .EQ. NB3) THEN
          ! These nodes go diagonally to 1U, 3L
          glob_node_id = min_current_global_id - 1 &
                       + (jj - 1) * nn1 &
                       + nn1 &
                       - (NB3 - 2) * NB2 * NB1 * nn_block &
                       - (NB2 - 1) * NB1 * nn_block &
                       - nn_block
        ELSE IF (iB(3,1) .EQ. NB3) THEN
          ! These nodes go to lower 3L face 
          glob_node_id = min_current_global_id - 1&
                       + (jj - 1) * nn1 &
                       + nn1 &
                       - (NB3 - 1) * NB2 * NB1 * nn_block &
                       - nn_block
        ELSE IF (iB(1,1) .EQ. 1) THEN
          ! These nodes go to upper 1U face
          glob_node_id = min_current_global_id - 1&
                       + (jj - 1) * nn1 &
                       + nn1 &
                       + NB2 * NB1 * nn_block &
                       + (NB1 - 1) * nn_block 
        ELSE
          ! these nodes have no periodicity
          glob_node_id = min_current_global_id - 1&
                       + (jj - 1) * nn1 &
                       + nn1 &
                       + (NB2 - 1) * NB1 * nn_block &
                       + (NB1 - 1) * nn_block 
        END IF
      END IF
    END IF
  END IF
  !--- 1U, 3U
  IF (ii .GT. N1p) THEN
    IF (jj .GE. S2p .AND. jj .LE. N2p) THEN
      IF (kk .GT. N3p) THEN
        IF (iB(1,1) .EQ. NB1 .AND. iB(3,1) .EQ. NB3) THEN
          ! These nodes go diagonally to 1L 3L
          glob_node_id = min_current_global_id &
                       + (jj - 1) * nn1 &
                       - (NB3 - 1) * NB1 * NB2 * nn_block &
                       - (NB1 - 1) * nn_block 
        ELSE IF (iB(3,1) .EQ. NB3) THEN
          ! These nodes to bottom 3L face
          glob_node_id = min_current_global_id &
                       + (jj - 1) * nn1 &
                       - ((NB3 - 1) * NB2 * NB1 - 1) *nn_block
        ELSE IF (iB(1,1) .EQ. NB1) THEN
          ! These nodes go to bottom 1L face
          glob_node_id = min_current_global_id &
                       + (jj - 1) * nn1 &
                       + ((NB2 - 1) * NB1 + 1) * nn_block
        ELSE
          ! These nodes have no periodicity
          glob_node_id = min_current_global_id &
                       + (jj - 1) * nn1 &
                       + (NB1*NB2 + 1) * nn_block
        END IF
      END IF
    END IF
  END IF
  !--- 2L, 3L
  IF (ii .GE. S1p .AND. ii .LE. N1p) THEN
    IF (jj .LT. S2p) THEN
      IF (kk .LT. S3p) THEN
        IF (iB(2,1) .EQ. 1 .AND. iB(3,1) .EQ. 1) THEN
          ! These nodes go diagonally to 2U 3U
          glob_node_id = max_current_global_id &
                       - (nn1 - ii) &
                       + (NB3 - 1) * NB2 * NB1 * nn_block &
                       + (NB2 - 1) * NB1 * nn_block 
        ELSE IF (iB(3,1) .EQ. 1) THEN
          ! These go to top 3U face
          glob_node_id = max_current_global_id &
                       - (nn1 - ii) &
                       + (NB3 - 2) * NB2 * NB1 * nn_block &
                       + (NB2 - 1) * NB1 * nn_block
        ELSE IF (iB(2,1) .EQ. 1) THEN
          ! These go to top 2U face
          glob_node_id = max_current_global_id &
                       - (nn1 - ii) &
                       - NB1 * nn_block
        ELSE
          ! These do not have periodicity
          glob_node_id = max_current_global_id &
                       - (nn1 - ii) &
                       - (NB2 + 1) * NB1 * nn_block
        END IF
      END IF
    END IF
  END IF
  !--- 2U, 3L
  IF (ii .GE. S1p .AND. ii .LE. N1p) THEN
    IF (jj .GT. N2p) THEN
      IF (kk .LT. S3p) THEN
        IF (iB(2,1) .EQ. NB2 .AND. iB(3,1) .EQ. 1) THEN
          ! These nodes go diagonally to 2L 3U
          glob_node_id = min_current_global_id &
                       + (nn3 - 1) * nn2 *nn1 &
                       + (ii - 1) &
                       + (NB3 - 2) * NB2 * NB1 * nn_block &
                       + NB1 * nn_block 
        ELSE IF (iB(3,1) .EQ. 1) THEN
          ! These go to to 3U face
          glob_node_id = min_current_global_id &
                       + (nn3 - 1) * nn2 *nn1 &
                       + (ii - 1) &
                       + (NB3 - 1) * NB2 * NB1 * nn_block &
                       +  NB1 * nn_block
        ELSE IF (iB(2,1) .EQ. NB2) THEN
          ! These go to bottom 2L face
          glob_node_id = min_current_global_id &
                       + (nn3 - 1) * nn2 *nn1 &
                       + (ii - 1) &
                       - NB1*(2*NB2 - 1) * nn_block
        ELSE 
          ! These do not have any periodicity
          glob_node_id = min_current_global_id &
                       + (nn3 - 1) * nn2 *nn1 &
                       + (ii - 1) &
                       - NB1*(NB2 - 1)* nn_block 
        END IF
      END IF
    END IF
  END IF
  !--- 2L, 3U
  IF (ii .GE. S1p .AND. ii .LE. N1p) THEN
    IF (jj .LT. S2p) THEN
      IF (kk .GT. N3p) THEN
        IF (iB(2,1) .EQ. 1 .AND. iB(3,1) .EQ. NB3) THEN
          ! these nodes go diagonally to 2U 3l
          glob_node_id = min_current_global_id &
                       + (nn2 - 1) * nn1 &
                       + (ii - 1) &
                       - (NB3 - 2) * NB2 * NB1 * nn_block &
                       - NB1 * nn_block 
        ELSE IF (iB(3,1) .EQ. NB3) THEN
          ! these go to bottom 3L face
          glob_node_id = min_current_global_id &
                       + (nn2 - 1) * nn1 &
                       + (ii - 1) &
                       - (NB3 - 1) * NB2 * NB1 * nn_block &
                       -  NB1 * nn_block
        ELSE IF (iB(2,1) .EQ. 1) THEN
          ! These go to top 2U face
          glob_node_id = min_current_global_id &
                       + (nn2 - 1) * nn1 &
                       + (ii - 1) &
                       + NB1*(2*NB2 - 1) * nn_block
        ELSE
          ! these have no periodicity
          glob_node_id = min_current_global_id &
                       + (nn2 - 1) * nn1 &
                       + (ii - 1) &
                       + NB1*(NB2 - 1)* nn_block 
        END IF
      END IF
    END IF
  END IF
  !--- 2U, 3U
  IF (ii .GE. S1p .AND. ii .LE. N1p) THEN
    IF (jj .GT. N2p) THEN
      IF (kk .GT. N3p) THEN
        IF (iB(2,1) .EQ. NB2 .AND. iB(3,1) .EQ. NB3) THEN
          ! These go diagonally to 2L 3L
          glob_node_id = min_current_global_id &
                       + (ii - 1) & 
                       - (NB3 - 1) * NB2 * NB1 * nn_block &
                       - (NB2 - 1) * NB1 * nn_block 
        ELSE IF (iB(3,1) .EQ. NB3) THEN
          ! These go to bottom 3L face
          glob_node_id = min_current_global_id &
                       + (ii - 1) & 
                       - (NB3 - 2) * NB2 * NB1 * nn_block &
                       - (NB2 - 1) * NB1 * nn_block
        ELSE IF (iB(2,1) .EQ. NB2) THEN
          ! these go to bottom 2L face
          glob_node_id = min_current_global_id &
                       + (ii - 1) & 
                       + NB1 * nn_block
        ELSE
          ! These have no periodicity
          glob_node_id = min_current_global_id &
                       + (ii - 1) & 
                       + (NB2 + 1) * NB1 * nn_block
        END IF
      END IF
    END IF
  END IF
  !--- ghost vertices in corners ---
  !--- 1L, 2L, 3L
  IF (ii .LT. S1p) THEN
    IF (jj .LT. S2p) THEN
      IF (kk .LT. S3p) THEN
        IF (iB(1,1) .EQ. 1 .AND. iB(2,1) .EQ. 1 .AND. iB(3,1)  .EQ. 1) THEN
          glob_node_id = max_current_global_id &
                       + (NB1 * NB2 * NB3 - 1) * nn_block
        ELSE IF (iB(1,1) .EQ. 1 .AND. iB(2,1) .EQ. 1) THEN
          glob_node_id = min_current_global_id &
                       - 1
        ELSE IF (iB(1,1) .EQ. 1 .AND. iB(3,1) .EQ. 1) THEN
          glob_node_id = max_current_global_id &
                       + ((NB3 - 1) *NB2 *NB1 -1 ) * nn_block
        ELSE IF (iB(2,1) .EQ. 1 .AND. iB(3,1) .EQ. 1) THEN
          glob_node_id = max_current_global_id &
                       + (NB3 - 2) * NB2 * NB1 * nn_block &
                       + 2*(NB2-1) * NB1 * nn_block &
                       + NB1 * nn_block
        ELSE IF (iB(1,1) .EQ. 1) THEN
          glob_node_id = max_current_global_id &
                       - (NB2 * NB1 + 1 ) * nn_block
        ELSE IF (iB(2,1) .EQ. 1) THEN
          glob_node_id = max_current_global_id &
                       - (NB1 + 1) * nn_block
        ELSE IF (iB(3,1) .EQ. 1) THEN
          glob_node_id = max_current_global_id &
                       + (NB3 - 2) *NB2 *NB1 *nn_block &
                       + ((NB2 - 1) * NB1 - 1) * nn_block
        ELSE
          glob_node_id = max_current_global_id &
                       - NB2 * NB1 * nn_block &
                       - (NB1 + 1) * nn_block

        END IF
      END IF
    END IF
  END IF
  !--- 1U, 2L, 3L
  IF (ii .GT. N1p) THEN
    IF (jj .LT. S2p) THEN
      IF (kk .LT. S3p) THEN

      END IF
    END IF
  END IF
  !--- 1L, 2U, 3L
  IF (ii .LT. S1p) THEN
    IF (jj .GT. N2p) THEN
      IF (kk .LT. S3p) THEN

      END IF
    END IF
  END IF
  !--- 1L, 2L, 3U
  IF (ii .LT. S1p) THEN
    IF (jj .LT. S2p) THEN
      IF (kk .GT. N3p) THEN

      END IF
    END IF
  END IF
  !--- 1U, 2U, 3L
  IF (ii .LT. S1p) THEN
    IF (jj .GT. N2p) THEN
      IF (kk .LT. S3p) THEN

      END IF
    END IF
  END IF
  !--- 1U, 2L, 3U
  IF (ii .GT. N1p) THEN
    IF (jj .LT. S2p) THEN
      IF (kk .GT. N3p) THEN

      END IF
    END IF
  END IF
  !--- 1L, 2U, 3U
  IF (ii .LT. S1p) THEN
    IF (jj .GT. N2p) THEN
      IF (kk .GT. N3p) THEN

      END IF
    END IF
  END IF
  !--- 1U, 2U, 3U
  IF (ii .GT. N1p) THEN
    IF (jj .GT. N2p) THEN
      IF (kk .GT. N3p) THEN

      END IF
    END IF
  END IF


  END SUBROUTINE local_to_global_node_id_loc_con_allperiodic


  SUBROUTINE local_to_global_node_id_loc_con_allperiodic3(loc_node_id, glob_node_id)

  IMPLICIT NONE

  INTEGER, INTENT(IN)    ::  loc_node_id
  INTEGER, INTENT(OUT)   ::  glob_node_id
  INTEGER                ::  nn1, nn2, nn3
  INTEGER                ::  nn_block
  INTEGER                ::  ii, jj, kk
  INTEGER                ::  max_current_global_id, min_current_global_id
  INTEGER                ::  n_local, n_ghost
  INTEGER                ::  from_block, to_block

  nn1 = N1 - 1
  nn2 = N2 - 1
  nn3 = N3 - 1

  nn_block = nn1 * nn2 * nn3

  !--- Get the current process' maximum id ---
  CALL global_cart2id_loc_con(max_current_global_id,N1p,N2p,N3p)

  !--- Get the current process' minumum id ---
  CALL global_cart2id_loc_con(min_current_global_id,S1p,S2p,S3p)

  !--- first convert locally to i,j,k notation ---
  CALL local_id2cart_loc_con(loc_node_id, ii, jj, kk, n_local, n_ghost)
  !--- Inner vertices (ghost vertices excluded) ---
  IF (ii .GE. S1p .AND. ii .LE. N1p) THEN
    IF (jj .GE. S2p .AND. jj .LE. N2p) THEN
      IF (kk .GE. S3p .AND. kk .LE. N3p) THEN
        CALL global_cart2id_loc_con(glob_node_id, ii, jj, kk)
      END IF
    END IF
  END IF

  !--- ghost vertices perpendicularly in front of a block face ---
  ! for each direction we need to check whether in that direction we have a
  ! neighbor that is in bulk, on face or periodic (inverted face)
  !--- face 1
  !--- top
  IF (ii .GT. N1p) THEN
    IF (jj .GE. S2p .AND. jj .LE. N2p) THEN
      IF (kk .GE. S3p .AND. kk .LE. N3p) THEN
        from_block   = block_id(iB(1,1)    , iB(2,1), iB(3,1))
        to_block     = block_id(iB(1,1) + 1, iB(2,1), iB(3,1))

        glob_node_id = min_current_global_id &
                     + (jj - 1) * nn1 &
                     + (kk - 1) * nn1 * nn2 &
                     + (to_block - from_block) * nn_block

      END IF
    END IF
  END IF
  !--- bottom
  IF (ii .LT. S1p) THEN
    IF (jj .GE. S2p .AND. jj .LE. N2p) THEN
      IF (kk .GE. S3p .AND. kk .LE. N3p) THEN
        from_block   = block_id(iB(1,1)    , iB(2,1), iB(3,1))
        to_block     = block_id(iB(1,1) - 1, iB(2,1), iB(3,1))

        glob_node_id = max_current_global_id &
                     - (nn2 - jj) * nn1 &
                     - (nn3 - kk) * nn1 * nn2 &
                     + (to_block - from_block) * nn_block

      END IF
    END IF
  END IF
  !--- face 2
  !--- top
  IF (ii .GE. S1p .AND. ii .LE. N1p) THEN
    IF (jj .GT. N2p) THEN
      IF (kk .GE. S3p .AND. kk .LE. N3p) THEN
        from_block   = block_id(iB(1,1), iB(2,1)     , iB(3,1))
        to_block     = block_id(iB(1,1), iB(2,1) +  1, iB(3,1))

        glob_node_id = min_current_global_id &
                     + (ii - 1) &
                     + (kk - 1) * nn1 * nn2 &
                     + (to_block - from_block) * nn_block

      END IF
    END IF
  END IF
  !--- bottom
  IF (ii .GE. S1p .AND. ii .LE. N1p) THEN
    IF (jj .LT. S2p) THEN
      IF (kk .GE. S3p .AND. kk .LE. N3p) THEN
        from_block   = block_id(iB(1,1), iB(2,1)     , iB(3,1))
        to_block     = block_id(iB(1,1), iB(2,1) -  1, iB(3,1))

        glob_node_id = max_current_global_id &
                     - (nn1 - ii) &
                     - (nn3 - kk) * nn1 * nn2 &
                     + (to_block - from_block) * nn_block

      END IF
    END IF
  END IF
  !--- face 3
  !--- top
  IF (ii .GE. S1p .AND. ii .LE. N1p) THEN
    IF (jj .GE. S2p .AND. jj .LE. N2p) THEN
      IF (kk .GT. N3p) THEN
        from_block   = block_id(iB(1,1), iB(2,1), iB(3,1)    )
        to_block     = block_id(iB(1,1), iB(2,1), iB(3,1) + 1)

        glob_node_id = min_current_global_id &
                     + (ii - 1) &
                     + (jj - 1) * nn1 &
                     + (to_block - from_block) * nn_block

      END IF
    END IF
  END IF
  !--- bottom
  IF (ii .GE. S1p .AND. ii .LE. N1p) THEN
    IF (jj .GE. S2p .AND. jj .LE. N2p) THEN
      IF (kk .LT. S3p) THEN
        from_block   = block_id(iB(1,1), iB(2,1), iB(3,1)    )
        to_block     = block_id(iB(1,1), iB(2,1), iB(3,1) - 1)

        glob_node_id = max_current_global_id &
                     - (nn1 - ii) &
                     - (nn2 - jj) * nn1 &
                     + (to_block - from_block) * nn_block

      END IF
    END IF
  END IF

  !--- ghost vertices diagonally off to edges ---
  !--- 1L, 2L
  IF (ii .LT. S1p) THEN
    IF (jj .LT. S2p) THEN
      IF (kk .GE. S3p .AND. kk .LE. N3p) THEN
        from_block   = block_id(iB(1,1)    , iB(2,1)    , iB(3,1))
        to_block     = block_id(iB(1,1) - 1, iB(2,1) - 1, iB(3,1))

        glob_node_id = max_current_global_id &
                     - (nn3 - kk) * nn1 * nn2 &
                     + (to_block - from_block) * nn_block

      END IF
    END IF
  END IF
  !--- 1U, 2L
  IF (ii .GT. N1p) THEN
    IF (jj .LT. S2p) THEN
      IF (kk .GE. S3p .AND. kk .LE. N3p) THEN
        from_block   = block_id(iB(1,1)    , iB(2,1)    , iB(3,1))
        to_block     = block_id(iB(1,1) + 1, iB(2,1) - 1, iB(3,1))

        glob_node_id = max_current_global_id &
                     - (nn1 - 1) &
                     - (nn3 - kk) * nn1 * nn2 &
                     + (to_block - from_block) * nn_block

      END IF
    END IF
  END IF
  !--- 1L, 2U
  IF (ii .LT. S1p) THEN
    IF (jj .GT. N2p) THEN
      IF (kk .GE. S3p .AND. kk .LE. N3p) THEN
        from_block   = block_id(iB(1,1)    , iB(2,1)    , iB(3,1))
        to_block     = block_id(iB(1,1) - 1, iB(2,1) + 1, iB(3,1))

        glob_node_id = min_current_global_id &
                     + (nn1 - 1) &
                     + (kk - 1) * nn1 * nn2 &
                     + (to_block - from_block) * nn_block

      END IF
    END IF
  END IF
  !--- 1U, 2U
  IF (ii .GT. N1p) THEN
    IF (jj .GT. N2p) THEN
      IF (kk .GE. S3p .AND. kk .LE. N3p) THEN
        from_block   = block_id(iB(1,1)    , iB(2,1)    , iB(3,1))
        to_block     = block_id(iB(1,1) + 1, iB(2,1) + 1, iB(3,1))

        glob_node_id = min_current_global_id &
                     + (kk - 1) * nn1 * nn2 &
                     + (to_block - from_block) * nn_block

      END IF
    END IF
  END IF
  !--- 1L, 3L
  IF (ii .LT. S1p) THEN
    IF (jj .GE. S2p .AND. jj .LE. N2p) THEN
      IF (kk .LT. S3p) THEN
        from_block   = block_id(iB(1,1)    , iB(2,1), iB(3,1)    )
        to_block     = block_id(iB(1,1) - 1, iB(2,1), iB(3,1) - 1)

        glob_node_id = max_current_global_id &
                     - (nn2 - jj) * nn1 &
                     + (to_block - from_block) * nn_block

      END IF
    END IF
  END IF
  !--- 1U, 3L
  IF (ii .GT. N1p) THEN
    IF (jj .GE. S2p .AND. jj .LE. N2p) THEN
      IF (kk .LT. S3p) THEN
        from_block   = block_id(iB(1,1)    , iB(2,1), iB(3,1)    )
        to_block     = block_id(iB(1,1) + 1, iB(2,1), iB(3,1) - 1)

        glob_node_id = min_current_global_id  &
                     + nn1 * nn2 * (nn3 - 1) &
                     + (jj - 1)*nn1 &
                     + (to_block - from_block) * nn_block

      END IF
    END IF
  END IF
  !--- 1L, 3U
  IF (ii .LT. S1p) THEN
    IF (jj .GE. S2p .AND. jj .LE. N2p) THEN
      IF (kk .GT. N3p) THEN
        from_block   = block_id(iB(1,1)    , iB(2,1), iB(3,1)    )
        to_block     = block_id(iB(1,1) - 1, iB(2,1), iB(3,1) + 1)

        glob_node_id = min_current_global_id &
                     + (nn1 - 1) &
                     + (jj - 1) * nn1 &
                     + (to_block - from_block) * nn_block

      END IF
    END IF
  END IF
  !--- 1U, 3U
  IF (ii .GT. N1p) THEN
    IF (jj .GE. S2p .AND. jj .LE. N2p) THEN
      IF (kk .GT. N3p) THEN
        from_block   = block_id(iB(1,1)    , iB(2,1), iB(3,1)    )
        to_block     = block_id(iB(1,1) + 1, iB(2,1), iB(3,1) + 1)

        glob_node_id = min_current_global_id &
                     + (jj - 1) * nn1 &
                     + (to_block - from_block) * nn_block

      END IF
    END IF
  END IF
  !--- 2L, 3L
  IF (ii .GE. S1p .AND. ii .LE. N1p) THEN
    IF (jj .LT. S2p) THEN
      IF (kk .LT. S3p) THEN
        from_block   = block_id(iB(1,1), iB(2,1)    , iB(3,1)    )
        to_block     = block_id(iB(1,1), iB(2,1) - 1, iB(3,1) - 1)

        glob_node_id = max_current_global_id &
                     - (nn1 - ii) &
                     + (to_block - from_block) * nn_block

      END IF
    END IF
  END IF
  !--- 2U, 3L
  IF (ii .GE. S1p .AND. ii .LE. N1p) THEN
    IF (jj .GT. N2p) THEN
      IF (kk .LT. S3p) THEN
        from_block   = block_id(iB(1,1), iB(2,1)    , iB(3,1)    )
        to_block     = block_id(iB(1,1), iB(2,1) + 1, iB(3,1) - 1)

        glob_node_id = min_current_global_id &
                     + (ii - 1) &
                     + (nn3 - 1) * nn2 *nn1 &
                     + (to_block - from_block) * nn_block

      END IF
    END IF
  END IF
  !--- 2L, 3U
  IF (ii .GE. S1p .AND. ii .LE. N1p) THEN
    IF (jj .LT. S2p) THEN
      IF (kk .GT. N3p) THEN
        from_block   = block_id(iB(1,1), iB(2,1)    , iB(3,1)    )
        to_block     = block_id(iB(1,1), iB(2,1) - 1, iB(3,1) + 1)

        glob_node_id = min_current_global_id &
                     + (ii - 1) &
                     + (nn2 - 1) * nn1 &
                     + (to_block - from_block) * nn_block

      END IF
    END IF
  END IF
  !--- 2U, 3U
  IF (ii .GE. S1p .AND. ii .LE. N1p) THEN
    IF (jj .GT. N2p) THEN
      IF (kk .GT. N3p) THEN
        from_block   = block_id(iB(1,1), iB(2,1)    , iB(3,1)    )
        to_block     = block_id(iB(1,1), iB(2,1) + 1, iB(3,1) + 1)

        glob_node_id = min_current_global_id &
                     + (ii - 1) & 
                     + (to_block - from_block) * nn_block

      END IF
    END IF
  END IF
  !--- ghost vertices in corners ---
  !--- 1L, 2L, 3L
  IF (ii .LT. S1p) THEN
    IF (jj .LT. S2p) THEN
      IF (kk .LT. S3p) THEN
        from_block   = block_id(iB(1,1)    , iB(2,1)    , iB(3,1)    )
        to_block     = block_id(iB(1,1) - 1, iB(2,1) - 1, iB(3,1) - 1)

        glob_node_id = max_current_global_id &
                     + (to_block - from_block) * nn_block

      END IF
    END IF
  END IF
  !--- 1U, 2L, 3L
  IF (ii .GT. N1p) THEN
    IF (jj .LT. S2p) THEN
      IF (kk .LT. S3p) THEN
        from_block   = block_id(iB(1,1)    , iB(2,1)    , iB(3,1)    )
        to_block     = block_id(iB(1,1) + 1, iB(2,1) - 1, iB(3,1) - 1)

        glob_node_id = max_current_global_id &
                     - (nn1 - 1) &
                     + (to_block - from_block) * nn_block

      END IF
    END IF
  END IF
  !--- 1L, 2U, 3L
  IF (ii .LT. S1p) THEN
    IF (jj .GT. N2p) THEN
      IF (kk .LT. S3p) THEN
        from_block   = block_id(iB(1,1)    , iB(2,1)    , iB(3,1)    )
        to_block     = block_id(iB(1,1) - 1, iB(2,1) + 1, iB(3,1) - 1)

        glob_node_id = min_current_global_id &
                     + (nn1 - 1) &
                     + (nn3 - 1) * nn2 * nn1 &
                     + (to_block - from_block) * nn_block

      END IF
    END IF
  END IF
  !--- 1L, 2L, 3U
  IF (ii .LT. S1p) THEN
    IF (jj .LT. S2p) THEN
      IF (kk .GT. N3p) THEN
        from_block   = block_id(iB(1,1)    , iB(2,1)    , iB(3,1)    )
        to_block     = block_id(iB(1,1) - 1, iB(2,1) - 1, iB(3,1) + 1)

        glob_node_id = min_current_global_id &
                     + (nn1 - 1) &
                     + (nn2 - 1) * nn1 &
                     + (to_block - from_block) * nn_block

      END IF
    END IF
  END IF
  !--- 1U, 2U, 3L
  IF (ii .GT. N1p) THEN
    IF (jj .GT. N2p) THEN
      IF (kk .LT. S3p) THEN
        from_block   = block_id(iB(1,1)    , iB(2,1)    , iB(3,1)    )
        to_block     = block_id(iB(1,1) + 1, iB(2,1) + 1, iB(3,1) - 1)

        glob_node_id = min_current_global_id &
                     + (nn3 - 1) * nn2 * nn1 &
                     + (to_block - from_block) * nn_block

      END IF
    END IF
  END IF
  !--- 1U, 2L, 3U
  IF (ii .GT. N1p) THEN
    IF (jj .LT. S2p) THEN
      IF (kk .GT. N3p) THEN
        from_block   = block_id(iB(1,1)    , iB(2,1)    , iB(3,1)    )
        to_block     = block_id(iB(1,1) + 1, iB(2,1) - 1, iB(3,1) + 1)

        glob_node_id = min_current_global_id &
                     + (nn2 - 1) * nn1 &
                     + (to_block - from_block) * nn_block

      END IF
    END IF
  END IF
  !--- 1L, 2U, 3U
  IF (ii .LT. S1p) THEN
    IF (jj .GT. N2p) THEN
      IF (kk .GT. N3p) THEN
        from_block   = block_id(iB(1,1)    , iB(2,1)    , iB(3,1)    )
        to_block     = block_id(iB(1,1) - 1, iB(2,1) + 1, iB(3,1) + 1)

        glob_node_id = min_current_global_id &
                     + (nn1 - 1) &
                     + (to_block - from_block) * nn_block

      END IF
    END IF
  END IF
  !--- 1U, 2U, 3U
  IF (ii .GT. N1p) THEN
    IF (jj .GT. N2p) THEN
      IF (kk .GT. N3p) THEN
        from_block   = block_id(iB(1,1)    , iB(2,1)    , iB(3,1)    )
        to_block     = block_id(iB(1,1) + 1, iB(2,1) + 1, iB(3,1) + 1)

        glob_node_id = min_current_global_id &
                     + (to_block - from_block) * nn_block

      END IF
    END IF
  END IF


  END SUBROUTINE local_to_global_node_id_loc_con_allperiodic3



  SUBROUTINE local_to_global_node_id_loc_con_allperiodic2(loc_node_id, glob_node_id)

  IMPLICIT NONE

  INTEGER, INTENT(IN)    ::  loc_node_id
  INTEGER, INTENT(OUT)   ::  glob_node_id
  INTEGER                ::  nn1, nn2, nn3
  INTEGER                ::  nn_block
  INTEGER                ::  ii, jj, kk
  INTEGER                ::  min_current_global_id
  INTEGER                ::  n_local, n_ghost
  INTEGER                ::  to_block, from_block
  INTEGER                ::  i_to_block, j_to_block, k_to_block

  nn1 = N1 - 1
  nn2 = N2 - 1
  nn3 = N3 - 1

  nn_block = nn1 * nn2 * nn3

  i_to_block = iB(1,1)
  j_to_block = iB(2,1)
  k_to_block = iB(3,1)

  !--- convert to i, j, k
  CALL local_id2cart_loc_con(loc_node_id, ii, jj, kk, n_local, n_ghost)

  !--- check if indices are greater or lesser than nip, sip and add/substract
  !--- nip.
  IF (ii .GT. N1p) THEN
    ii = ii - N1p
    i_to_block = i_to_block + 1
  ELSE IF (ii .LT. S1p) THEN
    ii = ii + N1p
    i_to_block = i_to_block - 1
  END IF
  IF (jj .GT. N2p) THEN
    jj = jj - N2p
    j_to_block = j_to_block + 1
  ELSE IF (jj .LT. S2p) THEN
    jj = jj + N2p
    j_to_block = j_to_block - 1
  END IF
  IF (kk .GT. N3p) THEN
    kk = kk - N3p
    k_to_block = k_to_block + 1
  ELSE IF (kk .LT. S3p) THEN
    kk = kk + N3p
    k_to_block = k_to_block - 1
  END IF

  !--- create local global variable ---
  ii = ii - 1
  jj = jj - 1
  kk = kk - 1
  glob_node_id = ii + jj * nn1 + kk * nn1 * nn2

  !--- augment global variable to be of global scope ---
  CALL global_cart2id_loc_con(min_current_global_id, S1p, S2p, S3p)
  from_block = block_id(iB(1,1), iB(2,1), iB(3,1))
  to_block   = block_id(i_to_block, j_to_block, k_to_block)
  glob_node_id = glob_node_id + min_current_global_id + (to_block - from_block) * nn_block

  END SUBROUTINE local_to_global_node_id_loc_con_allperiodic2



  SUBROUTINE global_to_local_node_id_loc_con_allperiodic(loc_node_id, glob_node_id, ii, jj, kk)

  IMPLICIT NONE

  INTEGER, INTENT(OUT)    :: loc_node_id
  INTEGER, INTENT(IN)     :: glob_node_id
  INTEGER                 :: nn1, nn2, nn3
  INTEGER                 :: nn_block
  INTEGER, INTENT(OUT)                 :: ii, jj, kk ! returnedonly for debugging
  INTEGER                 :: max_current_global_id, min_current_global_id
  INTEGER                 :: max_neighbor_global_id, min_neighbor_global_id
  INTEGER                 :: n_local, n_ghost
  INTEGER                 :: from_block, to_block
  INTEGER                 :: temp_id, ijn1

  nn1 = N1 - 1
  nn2 = N2 - 1
  nn3 = N3 - 1

  nn_block = nn1 * nn2 * nn3

  !--- Get the current process' maximum id ---
  CALL global_cart2id_loc_con(max_current_global_id,N1p,N2p,N3p)

  !--- Get the current process' minimum id ---
  CALL global_cart2id_loc_con(min_current_global_id,S1p,S2p,S3p)

  !--- Check whether our glob_node_id is within this process' global ids ---
  !--- If yes the conversion is straight forward ---
  IF (glob_node_id .GE. min_current_global_id .AND. glob_node_id .LE. max_current_global_id) THEN
    CALL global_id2cart_loc_con(glob_node_id, ii, jj, kk)
    CALL local_cart2id_loc_con(loc_node_id, ii, jj, kk, n_local, n_ghost)
    RETURN
  END IF

  !--- If it is not inside the current process' global ids then do the following ---
  !--- check whether it's within the global coordinates of face neighbors ---
  from_block = block_id(iB(1,1), iB(2,1), iB(3,1))
  !--- face 1
  !--- top
  to_block   = block_id(iB(1,1) + 1, iB(2,1), iB(3,1))
  min_neighbor_global_id = min_current_global_id + (to_block - from_block) * nn_block
  max_neighbor_global_id = max_current_global_id + (to_block - from_block) * nn_block
  IF (glob_node_id .GE. min_neighbor_global_id .AND. glob_node_id .LE. max_neighbor_global_id) THEN
    temp_id = glob_node_id - min_current_global_id &
                           - (to_block - from_block) * nn_block
    CALL id2cart(temp_id, ii, jj, kk, nn1, nn2, nn3)
    print*, ii, jj, kk
    !--- now account for correction
    IF (ii + 1 .EQ. S1p) THEN
      ii = ii + 1 + N1p
      jj = jj + 1
      kk = kk + 1
    !--- this catches the exception where the top and bottom neighbor is the
    !    same block
    ELSE IF (ii + 1 .EQ. N1p) THEN
      ii = ii + 1 - N1p
      jj = jj + 1
      kk = kk + 1
    END IF

    CALL local_cart2id_loc_con(loc_node_id, ii, jj, kk, n_local, n_ghost)
    print*, 'hello -1'
    RETURN
  END IF
  !--- bottom
  to_block = block_id(iB(1,1) - 1, iB(2,1), iB(3,1))
  min_neighbor_global_id = min_current_global_id + (to_block - from_block) * nn_block
  max_neighbor_global_id = max_current_global_id + (to_block - from_block) * nn_block
  IF (glob_node_id .GE. min_neighbor_global_id .AND. glob_node_id .LE. max_neighbor_global_id) THEN
    temp_id = - glob_node_id + max_current_global_id &
                             + (to_block - from_block) * nn_block
    CALL id2cart(temp_id, ii, jj, kk, nn1, nn2, nn3)
    !--- now account for correction
    ii = ii
    jj = nn2 - jj
    kk = nn3 - kk

    CALL local_cart2id_loc_con(loc_node_id, ii, jj, kk, n_local, n_ghost)
    print*, 'hello 0'
    RETURN
  END IF

  !--- face 2
  !--- top
  to_block = block_id(iB(1,1), iB(2,1) + 1, iB(3,1))
  min_neighbor_global_id = min_current_global_id + (to_block - from_block) * nn_block
  max_neighbor_global_id = max_current_global_id + (to_block - from_block) * nn_block
  IF (glob_node_id .GE. min_neighbor_global_id .AND. glob_node_id .LE. max_neighbor_global_id) THEN
    temp_id = glob_node_id - min_current_global_id &
                           - (to_block - from_block) * nn_block
    CALL id2cart(temp_id, ii, jj, kk, nn1, nn2, nn3)
    !--- now account for correction
    IF (jj + 1 .EQ. S2p) THEN
      ii = ii + 1
      jj = jj + 1 + N2p
      kk = kk + 1
    !--- this catches the exception where the top and bottom neighbor is the
    !    same block
    ELSE IF (jj + 1 .EQ. N2p) THEN
      ii = ii + 1
      jj = jj + 1 - N2p
      kk = kk + 1
    END IF

    CALL local_cart2id_loc_con(loc_node_id, ii, jj, kk, n_local, n_ghost)
    print*, 'hello 1'
    RETURN
  END IF
  !--- bottom
  to_block = block_id(iB(1,1), iB(2,1) - 1, iB(3,1))
  min_neighbor_global_id = min_current_global_id + (to_block - from_block) * nn_block
  max_neighbor_global_id = max_current_global_id + (to_block - from_block) * nn_block
  IF (glob_node_id .GE. min_neighbor_global_id .AND. glob_node_id .LE. max_neighbor_global_id) THEN
    temp_id = - glob_node_id + max_current_global_id &
                             + (to_block - from_block) * nn_block
    CALL id2cart(temp_id, ii, jj, kk, nn1, nn2, nn3)
    !--- now account for correction
    ii = nn1 - ii
    jj = jj
    kk = nn3 - kk

    CALL local_cart2id_loc_con(loc_node_id, ii, jj, kk, n_local, n_ghost)
    print*, 'hello 2'
    RETURN
  END IF

  !--- face 3
  !--- top
  to_block = block_id(iB(1,1), iB(2,1), iB(3,1) + 1)
  min_neighbor_global_id = min_current_global_id + (to_block - from_block) * nn_block
  max_neighbor_global_id = max_current_global_id + (to_block - from_block) * nn_block
  IF (glob_node_id .GE. min_neighbor_global_id .AND. glob_node_id .LE. max_neighbor_global_id) THEN
    temp_id = glob_node_id - min_current_global_id &
                           - (to_block - from_block) * nn_block
    CALL id2cart(temp_id, ii, jj, kk, nn1, nn2, nn3)
    !--- now account for correction
    IF (kk + 1 .EQ. S3p) THEN
      ii = ii + 1
      jj = jj + 1
      kk = kk + 1 + N3p
    !--- this catches the exception where the top and bottom neighbor is the
    !    same block
    ELSE IF (kk + 1 .EQ. N3p) THEN
      ii = ii + 1
      jj = jj + 1
      kk = kk + 1 - N3p
    END IF

    CALL local_cart2id_loc_con(loc_node_id, ii, jj, kk, n_local, n_ghost)
    print*, 'hello 3'
    RETURN
  END IF
  !--- bottom
  to_block = block_id(iB(1,1), iB(2,1), iB(3,1) - 1)
  min_neighbor_global_id = min_current_global_id + (to_block - from_block) * nn_block
  max_neighbor_global_id = max_current_global_id + (to_block - from_block) * nn_block
  IF (glob_node_id .GE. min_neighbor_global_id .AND. glob_node_id .LE. max_neighbor_global_id) THEN
    temp_id = - glob_node_id + max_current_global_id &
                             + (to_block - from_block) * nn_block
    CALL id2cart(temp_id, ii, jj, kk, nn1, nn2, nn3)
    !--- now account for correction
    ii = nn1 - ii
    jj = nn2 - jj
    kk = kk

    CALL local_cart2id_loc_con(loc_node_id, ii, jj, kk, n_local, n_ghost)
    print*, 'hello 4'
    RETURN
  END IF

  !--- check whether it's in the global ids of edge neighbors ---
  !--- 1U, 2U
  to_block = block_id(iB(1,1) + 1, iB(2,1) + 1, iB(3,1))
  min_neighbor_global_id = min_current_global_id + (to_block - from_block) * nn_block
  max_neighbor_global_id = max_current_global_id + (to_block - from_block) * nn_block
  IF (glob_node_id .GE. min_neighbor_global_id .AND. glob_node_id .LE. max_neighbor_global_id) THEN
    temp_id = glob_node_id - min_current_global_id &
                           - (to_block - from_block) * nn_block
    CALL id2cart(temp_id, ii, jj, kk, nn1, nn2, nn3)
    !--- now account for correction
    IF (ii + 1 .EQ. S1p .AND. jj + 1 .EQ. S2p) THEN
      ii = ii + 1 + N1p
      jj = jj + 1 + N2p
      kk = kk + 1
    !--- this catches the exception where the top and bottom neighbor is the
    !    same block
    ELSE IF (ii + 1 .EQ. N1p .AND. jj + 1 .EQ. N2p) THEN
      ii = ii + 1 - N1p
      jj = jj + 1 - N2p
      kk = kk + 1
    ELSE IF (ii + 1 .EQ. S1p .AND. jj + 1 .EQ. N2p) THEN
      ii = ii + 1 + N1p
      jj = jj + 1 - N1p
      kk = kk + 1
    ELSE IF (ii + 1 .EQ. N1p .AND. jj + 1 .EQ. S2p) THEN
      ii = ii + 1 - N1p
      jj = jj + 1 + N1p
      kk = kk + 1
    END IF

    CALL local_cart2id_loc_con(loc_node_id, ii, jj, kk, n_local, n_ghost)
    print*, 'hello 5'
    RETURN
  END IF
  !--- 1L, 2L
  to_block = block_id(iB(1,1) - 1, iB(2,1) - 1, iB(3,1))
  min_neighbor_global_id = min_current_global_id + (to_block - from_block) * nn_block
  max_neighbor_global_id = max_current_global_id + (to_block - from_block) * nn_block
  IF (glob_node_id .GE. min_neighbor_global_id .AND. glob_node_id .LE. max_neighbor_global_id) THEN
    temp_id = - glob_node_id + max_current_global_id &
                             + (to_block - from_block) * nn_block
    CALL id2cart(temp_id, ii, jj, kk, nn1, nn2, nn3)
    !--- now account for correction
    ii = ii
    jj = jj
    kk = nn3 - kk

    CALL local_cart2id_loc_con(loc_node_id, ii, jj, kk, n_local, n_ghost)
    print*, 'hello 66'
    RETURN
  END IF
  !--- 1L, 2U
  to_block = block_id(iB(1,1) - 1, iB(2,1) + 1, iB(3,1))
  min_neighbor_global_id = min_current_global_id + (to_block - from_block) * nn_block
  max_neighbor_global_id = max_current_global_id + (to_block - from_block) * nn_block
  IF (glob_node_id .GE. min_neighbor_global_id .AND. glob_node_id .LE. max_neighbor_global_id) THEN
    temp_id = glob_node_id - min_current_global_id &
                           - (to_block - from_block) * nn_block
    CALL id2cart(temp_id, ii, jj, kk, nn1, nn2, nn3)
    !--- now account for correction
    IF (ii + 1 .EQ. N1p .AND. jj + 1 .EQ. S2p) THEN
      ii = ii + 1 - N1p
      jj = jj + 1 + N2p
      kk = kk + 1
    !--- this catches the exception where the top and bottom neighbor is the
    !    same block
    ELSE IF (ii + 1 .EQ. S1p .AND. jj + 1 .EQ. N2p) THEN
      ii = ii + 1 + N1p
      jj = jj + 1 - S2p
      kk = kk + 1
    END IF

    CALL local_cart2id_loc_con(loc_node_id, ii, jj, kk, n_local, n_ghost)
    print*, 'hello 7'
    RETURN
  END IF
  !--- 1U, 2L
  to_block = block_id(iB(1,1) + 1, iB(2,1) - 1, iB(3,1))
  min_neighbor_global_id = min_current_global_id + (to_block - from_block) * nn_block
  max_neighbor_global_id = max_current_global_id + (to_block - from_block) * nn_block
  IF (glob_node_id .GE. min_neighbor_global_id .AND. glob_node_id .LE. max_neighbor_global_id) THEN
    temp_id = - glob_node_id + max_current_global_id &
                             + (to_block - from_block) * nn_block
    CALL id2cart(temp_id, ii, jj, kk, nn1, nn2, nn3)
    !--- now account for correction
    ii = ii - 2 + N1p
    jj = jj
    kk = nn3 - kk

    CALL local_cart2id_loc_con(loc_node_id, ii, jj, kk, n_local, n_ghost)
    print*, 'hello 6'
    RETURN
  END IF
  !--- 1U, 3U
  to_block = block_id(iB(1,1) + 1, iB(2,1), iB(3,1) + 1)
  min_neighbor_global_id = min_current_global_id + (to_block - from_block) * nn_block
  max_neighbor_global_id = max_current_global_id + (to_block - from_block) * nn_block
  IF (glob_node_id .GE. min_neighbor_global_id .AND. glob_node_id .LE. max_neighbor_global_id) THEN
    temp_id = glob_node_id - min_current_global_id &
                           - (to_block - from_block) * nn_block
    CALL id2cart(temp_id, ii, jj, kk, nn1, nn2, nn3)
    !--- now account for correction
    IF (ii + 1 .EQ. S1p .AND. kk + 1 .EQ. S3p) THEN
      ii = ii + 1 + N1p
      jj = jj + 1
      kk = kk + 1 + N3p
    ELSE IF (ii + 1 .EQ. N1p .AND. kk + 1 .EQ. N3p) THEN
      ii = ii + 1 - N1p
      jj = jj + 1
      kk = kk + 1 - N3p
    ELSE IF (ii + 1 .EQ. S1p .AND. kk + 1 .EQ. N3p) THEN
      ii = ii + 1 + N1p
      jj = jj + 1
      kk = kk + 1 - N3p
    ELSE IF (ii + 1 .EQ. N1p .AND. kk + 1 .EQ. S3p) THEN
      ii = ii + 1 - N1p
      jj = jj + 1
      kk = kk + 1 + N3p
    END IF


    CALL local_cart2id_loc_con(loc_node_id, ii, jj, kk, n_local, n_ghost)
    print*, 'hello 12'
    RETURN
  END IF
  !--- 1L, 3L
  to_block = block_id(iB(1,1) - 1, iB(2,1), iB(3,1) - 1)
  min_neighbor_global_id = min_current_global_id + (to_block - from_block) * nn_block
  max_neighbor_global_id = max_current_global_id + (to_block - from_block) * nn_block
  IF (glob_node_id .GE. min_neighbor_global_id .AND. glob_node_id .LE. max_neighbor_global_id) THEN
    temp_id = - glob_node_id + max_current_global_id &
                             + (to_block - from_block) * nn_block
    CALL id2cart(temp_id, ii, jj, kk, nn1, nn2, nn3)
    !--- now account for correction
    ii = ii
    jj = nn2 - jj
    kk = kk

    CALL local_cart2id_loc_con(loc_node_id, ii, jj, kk, n_local, n_ghost)
    print*, 'hello 9'
    RETURN
  END IF
  !--- 1L, 3U
  to_block = block_id(iB(1,1) - 1, iB(2,1), iB(3,1) + 1)
  min_neighbor_global_id = min_current_global_id + (to_block - from_block) * nn_block
  max_neighbor_global_id = max_current_global_id + (to_block - from_block) * nn_block
  IF (glob_node_id .GE. min_neighbor_global_id .AND. glob_node_id .LE. max_neighbor_global_id) THEN
    temp_id = glob_node_id - min_current_global_id &
                           - (to_block - from_block) * nn_block
    CALL id2cart(temp_id, ii, jj, kk, nn1, nn2, nn3)
    !--- now account for correction
    IF (ii + 1 .EQ. N1p .AND. kk + 1 .EQ. S3p) THEN
      ii = ii + 1 - N1p
      jj = jj + 1
      kk = kk + 1 + N3p
    ELSE IF (ii + 1 .EQ. S1p .AND. kk + 1 .EQ. N3p) THEN
      ii = ii + 1 + N1p
      jj = jj + 1
      kk = kk + 1 - N3p
    END IF

    CALL local_cart2id_loc_con(loc_node_id, ii, jj, kk, n_local, n_ghost)
    print*, 'hello 11'
    RETURN
  END IF
  !--- 1U, 3L
  to_block = block_id(iB(1,1) + 1, iB(2,1), iB(3,1) - 1)
  min_neighbor_global_id = min_current_global_id + (to_block - from_block) * nn_block
  max_neighbor_global_id = max_current_global_id + (to_block - from_block) * nn_block
  IF (glob_node_id .GE. min_neighbor_global_id .AND. glob_node_id .LE. max_neighbor_global_id) THEN
    temp_id = glob_node_id - min_current_global_id &
                           - (to_block - from_block) * nn_block
    CALL id2cart(temp_id, ii, jj, kk, nn1, nn2, nn3)
    !--- now account for correction
    ii = ii + 1 + N1p
    jj = jj + 1
    kk = kk + 1 - N3p
    CALL local_cart2id_loc_con(loc_node_id, ii, jj, kk, n_local, n_ghost)
    print*, 'hello 10'
    RETURN
  END IF
  !--- 2U, 3U
  to_block = block_id(iB(1,1), iB(2,1) + 1, iB(3,1) + 1)
  min_neighbor_global_id = min_current_global_id + (to_block - from_block) * nn_block
  max_neighbor_global_id = max_current_global_id + (to_block - from_block) * nn_block
  IF (glob_node_id .GE. min_neighbor_global_id .AND. glob_node_id .LE. max_neighbor_global_id) THEN
    temp_id = glob_node_id - min_current_global_id &
                           - (to_block - from_block) * nn_block
    CALL id2cart(temp_id, ii, jj, kk, nn1, nn2, nn3)
    !--- now account for correctioin
    IF (jj + 1 .EQ. S2p .AND. kk + 1 .EQ. S3p) THEN
      ii = ii + 1
      jj = jj + 1 + N2p
      kk = kk + 1 + N3p
    ELSE IF (jj + 1 .EQ. N2p .AND. kk + 1 .EQ. N3p) THEN
      ii = ii + 1
      jj = jj + 1 - N2p
      kk = kk + 1 - N3p
    ELSE IF (jj + 1 .EQ. S2p .AND. kk + 1 .EQ. N3p) THEN
      ii = ii + 1
      jj = jj + 1 + N2p
      kk = kk + 1 - N3p
    ELSE IF (jj + 1 .EQ. N2p .AND. kk + 1 .EQ. S3p) THEN
      ii = ii + 1
      jj = jj + 1 - N2p
      kk = kk + 1 + N3p
    END IF

    CALL local_cart2id_loc_con(loc_node_id, ii, jj, kk, n_local, n_ghost)
    print*, 'hello 16'
    RETURN
  END IF
  !--- 2L, 3L
  to_block = block_id(iB(1,1), iB(2,1) - 1, iB(3,1) - 1)
  min_neighbor_global_id = min_current_global_id + (to_block - from_block) * nn_block
  max_neighbor_global_id = max_current_global_id + (to_block - from_block) * nn_block
  IF (glob_node_id .GE. min_neighbor_global_id .AND. glob_node_id .LE. max_neighbor_global_id) THEN
    temp_id = - glob_node_id + max_current_global_id &
                             + (to_block - from_block) * nn_block
    CALL id2cart(temp_id, ii, jj, kk, nn1, nn2, nn3)
    !--- now account for correction
    ii = nn1 - ii
    jj = jj
    kk = kk

    CALL local_cart2id_loc_con(loc_node_id, ii, jj, kk, n_local, n_ghost)
    print*, 'hello 13'
    RETURN
  END IF
  !--- 2L, 3U
  to_block = block_id(iB(1,1), iB(2,1) - 1, iB(3,1) + 1)
  min_neighbor_global_id = min_current_global_id + (to_block - from_block) * nn_block
  max_neighbor_global_id = max_current_global_id + (to_block - from_block) * nn_block
  IF (glob_node_id .GE. min_neighbor_global_id .AND. glob_node_id .LE. max_neighbor_global_id) THEN
    temp_id = glob_node_id - min_current_global_id &
                           - (to_block - from_block) * nn_block
    CALL id2cart(temp_id, ii, jj, kk, nn1, nn2, nn3)
    !--- now account for correction
    IF (jj + 1 .EQ. N2p .AND. kk + 1 .EQ. S3p) THEN
      ii = ii + 1
      jj = jj + 1 - N2p
      kk = kk + 1 + N3p
    ELSE IF (jj + 1 .EQ. S2p .AND. kk + 1 .EQ. N3p) THEN
      ii = ii + 1
      jj = jj + 1 + N2p
      kk = kk + 1 - N3p
    END IF

    CALL local_cart2id_loc_con(loc_node_id, ii, jj, kk, n_local, n_ghost)
    print*, 'hello 15'
    RETURN
  END IF
  !--- 2U, 3L
  to_block = block_id(iB(1,1), iB(2,1) + 1, iB(3,1) - 1)
  min_neighbor_global_id = min_current_global_id + (to_block - from_block) * nn_block
  max_neighbor_global_id = max_current_global_id + (to_block - from_block) * nn_block
  IF (glob_node_id .GE. min_neighbor_global_id .AND. glob_node_id .LE. max_neighbor_global_id) THEN
    temp_id = glob_node_id - min_current_global_id &
                           - (to_block - from_block) * nn_block
    CALL id2cart(temp_id, ii, jj, kk, nn1, nn2, nn3)
    !--- now account for correction
    ii = ii + 1
    jj = jj + 1 + N2p
    kk = kk + 1 - N3p

    CALL local_cart2id_loc_con(loc_node_id, ii, jj, kk, n_local, n_ghost)
    print*, 'hello 14'
    RETURN
  END IF

  !--- check whether it's in global ids of corner neighbors ---
  !--- 1U, 2U, 3U
  to_block = block_id(iB(1,1) + 1, iB(2,1) + 1, iB(3,1) + 1)
  min_neighbor_global_id = min_current_global_id + (to_block - from_block) * nn_block
  max_neighbor_global_id = max_current_global_id + (to_block - from_block) * nn_block
  IF (glob_node_id .GE. min_neighbor_global_id .AND. glob_node_id .LE. max_neighbor_global_id) THEN
    temp_id = glob_node_id - min_current_global_id &
                           - (to_block - from_block) * nn_block
    CALL id2cart(temp_id, ii, jj, kk, nn1, nn2, nn3)
    !--- now account for correction
    IF (ii + 1 .EQ. S1p .AND. jj + 1 .EQ. S2p .AND. kk + 1 .EQ. S3p) THEN
      ii = ii + 1 + N1p
      jj = jj + 1 + N2p
      kk = kk + 1 + N3p
    ELSE IF (ii + 1 .EQ. N1p .AND. jj + 1 .EQ. S2p .AND. kk + 1 .EQ. S3p) THEN
      ii = ii + 1 - N1p
      jj = jj + 1 + N2p
      kk = kk + 1 + N3p
    ELSE IF (ii + 1 .EQ. S1p .AND. jj + 1 .EQ. N2p .AND. kk + 1 .EQ. S3p) THEN
      ii = ii + 1 + N1p
      jj = jj + 1 - N2p
      kk = kk + 1 + N3p
    ELSE IF (ii + 1 .EQ. S1p .AND. jj + 1 .EQ. S2p .AND. kk + 1 .EQ. N3p) THEN
      ii = ii + 1 + N1p
      jj = jj + 1 + N2p
      kk = kk + 1 - N3p
    ELSE IF (ii + 1 .EQ. N1p .AND. jj + 1 .EQ. N2p .AND. kk + 1 .EQ. S3p) THEN
      ii = ii + 1 - N1p
      jj = jj + 1 - N2p
      kk = kk + 1 + N3p
    ELSE IF (ii + 1 .EQ. N1p .AND. jj + 1 .EQ. S2p .AND. kk + 1 .EQ. N3p) THEN
      ii = ii + 1 - N1p
      jj = jj + 1 + N2p
      kk = kk + 1 - N3p
    ELSE IF (ii + 1 .EQ. S1p .AND. jj + 1 .EQ. N2p .AND. kk + 1 .EQ. N3p) THEN
      ii = ii + 1 + N1p
      jj = jj + 1 - N2p
      kk = kk + 1 - N3p
    ELSE IF (ii + 1 .EQ. N1p .AND. jj + 1 .EQ. N2p .AND. kk + 1 .EQ. N3p) THEN
      ii = ii + 1 - N1p
      jj = jj + 1 - N2p
      kk = kk + 1 - N3p
    END IF

    CALL local_cart2id_loc_con(loc_node_id, ii, jj, kk, n_local, n_ghost)
    print*, 'hello 24'
    RETURN
  END IF
  !--- 1L, 2L, 3L
  to_block = block_id(iB(1,1) - 1, iB(2,1) - 1, iB(3,1) - 1)
  min_neighbor_global_id = min_current_global_id + (to_block - from_block) * nn_block
  max_neighbor_global_id = max_current_global_id + (to_block - from_block) * nn_block
  IF (glob_node_id .GE. min_neighbor_global_id .AND. glob_node_id .LE. max_neighbor_global_id) THEN
    temp_id = glob_node_id - max_current_global_id &
                           - (to_block - from_block) * nn_block
    CALL id2cart(temp_id, ii, jj, kk, nn1, nn2, nn3)
    !--- now account for correction
    ii = ii
    jj = jj
    kk = kk

    CALL local_cart2id_loc_con(loc_node_id, ii, jj, kk, n_local, n_ghost)
    print*, 'hello 17'
    RETURN
  END IF
  !--- 1U, 2L, 3L
  to_block = block_id(iB(1,1) + 1, iB(2,1) - 1, iB(3,1) - 1)
  min_neighbor_global_id = min_current_global_id + (to_block - from_block) * nn_block
  max_neighbor_global_id = max_current_global_id + (to_block - from_block) * nn_block
  IF (glob_node_id .GE. min_neighbor_global_id .AND. glob_node_id .LE. max_neighbor_global_id) THEN
    temp_id = glob_node_id - max_current_global_id &
                           - (to_block - from_block) * nn_block
    CALL id2cart(temp_id, ii, jj, kk, nn1, nn2, nn3)
    !--- now account for correction
    ii = ii - 1 + nn1
    jj = jj
    kk = kk

    CALL local_cart2id_loc_con(loc_node_id, ii, jj, kk, n_local, n_ghost)
    print*, 'hello 18'
    RETURN
  END IF
  !--- 1L, 2U, 3L
  to_block = block_id(iB(1,1) - 1, iB(2,1) + 1, iB(3,1) - 1)
  min_neighbor_global_id = min_current_global_id + (to_block - from_block) * nn_block
  max_neighbor_global_id = max_current_global_id + (to_block - from_block) * nn_block
  IF (glob_node_id .GE. min_neighbor_global_id .AND. glob_node_id .LE. max_neighbor_global_id) THEN
    temp_id = glob_node_id - min_current_global_id &
                           - (to_block - from_block) * nn_block
    CALL id2cart(temp_id, ii, jj, kk, nn1, nn2, nn3)
    !--- now account for correction
    ii = ii + 1 - nn1
    jj = jj + N2p
    kk = kk + 1 - nn3

    CALL local_cart2id_loc_con(loc_node_id, ii, jj, kk, n_local, n_ghost)
    print*, 'hello 19'
    RETURN
  END IF
  !--- 1L, 2L, 3U
  to_block = block_id(iB(1,1) - 1, iB(2,1) - 1, iB(3,1) + 1)
  min_neighbor_global_id = min_current_global_id + (to_block - from_block) * nn_block
  max_neighbor_global_id = max_current_global_id + (to_block - from_block) * nn_block
  IF (glob_node_id .GE. min_neighbor_global_id .AND. glob_node_id .LE. max_neighbor_global_id) THEN
    temp_id = glob_node_id - min_current_global_id &
                           - (to_block - from_block) * nn_block
    CALL id2cart(temp_id, ii, jj, kk, nn1, nn2, nn3)
    !--- now account for correction
    ii = ii + 1 - nn1
    jj = jj + 1 - nn2
    kk = kk + N3p

    CALL local_cart2id_loc_con(loc_node_id, ii, jj, kk, n_local, n_ghost)
    print*, 'hello 20'
    RETURN
  END IF
  !--- 1U, 2U, 3L
  to_block = block_id(iB(1,1) + 1, iB(2,1) + 1, iB(3,1) - 1)
  min_neighbor_global_id = min_current_global_id + (to_block - from_block) * nn_block
  max_neighbor_global_id = max_current_global_id + (to_block - from_block) * nn_block
  IF (glob_node_id .GE. min_neighbor_global_id .AND. glob_node_id .LE. max_neighbor_global_id) THEN
    temp_id = glob_node_id - min_current_global_id &
                           - (to_block - from_block) * nn_block
    CALL id2cart(temp_id, ii, jj, kk, nn1, nn2, nn3)
    !--- now account for correction
    ii = ii + N1p
    jj = jj + N2p
    kk = kk + 1 - nn3

    CALL local_cart2id_loc_con(loc_node_id, ii, jj, kk, n_local, n_ghost)
    print*, 'hello 21'
    RETURN
  END IF
  !--- 1U, 2L, 3U
  to_block = block_id(iB(1,1) + 1, iB(2,1) - 1, iB(3,1) + 1)
  min_neighbor_global_id = min_current_global_id + (to_block - from_block) * nn_block
  max_neighbor_global_id = max_current_global_id + (to_block - from_block) * nn_block
  IF (glob_node_id .GE. min_neighbor_global_id .AND. glob_node_id .LE. max_neighbor_global_id) THEN
    temp_id = glob_node_id - min_current_global_id &
                           - (to_block - from_block) * nn_block
    CALL id2cart(temp_id, ii, jj, kk, nn1, nn2, nn3)
    !--- now account for correction
    ii = ii + N1p
    jj = jj + 1 - nn2
    kk = kk + N3p

    CALL local_cart2id_loc_con(loc_node_id, ii, jj, kk, n_local, n_ghost)
    print*, 'hello 22'
    RETURN
  END IF
  !--- 1L, 2U, 3U
  to_block = block_id(iB(1,1) - 1, iB(2,1) + 1, iB(3,1) + 1)
  min_neighbor_global_id = min_current_global_id + (to_block - from_block) * nn_block
  max_neighbor_global_id = max_current_global_id + (to_block - from_block) * nn_block
  IF (glob_node_id .GE. min_neighbor_global_id .AND. glob_node_id .LE. max_neighbor_global_id) THEN
    temp_id = glob_node_id - min_current_global_id &
                           - (to_block - from_block) * nn_block
    CALL id2cart(temp_id, ii, jj, kk, nn1, nn2, nn3)
    !--- now account for correction
    ii = ii + 1 - nn1
    jj = jj + 1 + N2p
    kk = kk + 1 + N3p

    CALL local_cart2id_loc_con(loc_node_id, ii, jj, kk, n_local, n_ghost)
    print*, 'hello 23'
    RETURN
  END IF


  END SUBROUTINE global_to_local_node_id_loc_con_allperiodic



  SUBROUTINE global_to_local_node_id_loc_con_allperiodic2(loc_node_id, glob_node_id)

  IMPLICIT NONE

  INTEGER, INTENT(OUT)    :: loc_node_id
  INTEGER, INTENT(IN)     :: glob_node_id
  INTEGER                 :: nn1, nn2, nn3
  INTEGER                 :: ii, jj, kk
  INTEGER                 :: nn_block
  INTEGER                 :: min_current_global_id, max_current_global_id
  INTEGER                 :: n_local, n_ghost
  INTEGER                 :: temp_id, rel_id
  INTEGER                 :: nn_block_mult
  LOGICAL                 :: i_block_shift, j_block_shift, k_block_shift
  INTEGER                 :: current_block_id, dest_block_id
  INTEGER                 :: i_block_dest, j_block_dest, k_block_dest

  nn1 = N1 - 1
  nn2 = N2 - 1
  nn3 = N3 - 1

  nn_block = nn1 * nn2 * nn3

  !--- Get the current process' minimum and maximum id ---
  CALL global_cart2id_loc_con(min_current_global_id, S1p, S2p, S3p)
  CALL global_cart2id_loc_con(max_current_global_id, N1p, N2p, N3p)

  !--- Check whether our glob_node_id is within this process' global ids ---
  !--- If yes the conversion is straight forward ---
  IF (glob_node_id .GE. min_current_global_id .AND. glob_node_id .LE. max_current_global_id) THEN
    CALL global_id2cart_loc_con(glob_node_id, ii, jj, kk)
    CALL local_cart2id_loc_con(loc_node_id, ii, jj, kk, n_local, n_ghost)
  ELSE
  !--- relate the global id to the current block's minimum id ---
  rel_id = glob_node_id - min_current_global_id

  !--- calculate the remainder (local relative global index) ---
  temp_id = MOD(rel_id, nn_block) ! >/< 0

  !--- if temp_id is negative, account for that ---
  temp_id = MODULO(temp_id, nn_block)

  !--- now convert to i,j,k ---
  CALL id2cart(temp_id, ii, jj, kk, nn1, nn2, nn3)

  !--- check which of the block-in-question's positions lie shifted
  !--- we have to do this because the mapping is not bijective by itself
  !--- floating point arithmetic
  i_block_shift = .FALSE.
  j_block_shift = .FALSE.
  k_block_shift = .FALSE.
  IF (rel_id .LT. 1) rel_id = rel_id - nn_block
  ! This tells us how many block we move forward(+)/backward(-)
  nn_block_mult = rel_id / nn_block ! this should be an integer!
  ! Get current block's id
  current_block_id = block_id(iB(1,1), iB(2,1), iB(3,1))
  ! Calculate the desination block_id
  dest_block_id = current_block_id + nn_block_mult
  ! Revert the destination block id to i,j,k notation
  CALL block_cart(dest_block_id, i_block_dest, j_block_dest, k_block_dest)

  ! set the shift flags depending on whether the desination block's 
  ! and the current block's indices differ
  IF (iB(1,1) .NE. i_block_dest) i_block_shift = .TRUE.
  IF (iB(2,1) .NE. j_block_dest) j_block_shift = .TRUE.
  IF (iB(3,1) .NE. k_block_dest) k_block_shift = .TRUE.


  !--- check if any of the indices match sip or nip, if they do, 
  !--- add or substract nip. In any case add 1 ---
  ii = ii + 1
  jj = jj + 1
  kk = kk + 1
  
  IF (ii .EQ. S1p) THEN
    IF (i_block_shift) ii = ii + N1p
  ELSE IF (ii .EQ. N1p) THEN
    IF (i_block_shift) ii = ii - N1p
  END IF
  IF (jj .EQ. S2p) THEN
    IF (j_block_shift) jj = jj + N2p
  ELSE IF (jj .EQ. N2p) THEN
    IF (j_block_shift) jj = jj - N2p
  END IF
  IF (kk .EQ. S3p) THEN
    IF (k_block_shift) kk = kk + N3p
  ELSE IF (kk .EQ. N3p) THEN
    IF (k_block_shift) kk = kk - N3p
  END IF

  !--- now convert back to local index ---
  CALL local_cart2id_loc_con(loc_node_id, ii, jj, kk, n_local, n_ghost)

  END IF

  END SUBROUTINE global_to_local_node_id_loc_con_allperiodic2






  FUNCTION block_id(ib1,ib2,ib3) RESULT(fn_val)

  IMPLICIT NONE 

  INTEGER                ::  fn_val
  INTEGER, INTENT(IN)    ::  ib1, ib2, ib3
  INTEGER                ::  ii, jj, kk
 
  !--- respect periodicity ---
  ii = MODULO(ib1-1, NB1) + 1
  jj = MODULO(ib2-1, NB2) + 1
  kk = MODULO(ib3-1, NB3) + 1

  !--- column major ---
  !fn_val = (ii - 1) + (jj - 1) * NB1 + (kk - 1) * NB1 * NB2

  !--- row major ---
  fn_val = (kk - 1) + (jj - 1) * NB3 + (ii - 1) * NB2 * NB3

  END FUNCTION block_id


  SUBROUTINE block_cart(block_id, ib, jb, kb)

  IMPLICIT NONE

  INTEGER, INTENT(IN)    :: block_id
  INTEGER, INTENT(OUT)   :: ib, jb, kb
  INTEGER                :: bid

  bid = block_id

  !--- column major ---
  !CALL id2cart(bid, ib, jb, kb, NB1, NB2, NB3)

  !--- row major ---
  CALL id2cart(bid, kb, jb, ib, NB3, NB2, NB1)

  ib = ib + 1
  jb = jb + 1
  kb = kb + 1

  END SUBROUTINE block_cart 


  SUBROUTINE global_cart2id_glob_con(node_id, i, j, k, vel_grid_yes, boundary_yes, dir)

  IMPLICIT NONE

  INTEGER, INTENT(IN)     ::  i,j,k
  INTEGER                 ::  ii,jj,kk
  INTEGER                 ::  nn1,nn2,nn3
  INTEGER                 ::  nn1l, nn1u, nn2l, nn2u, nn3l, nn3u
  INTEGER                 ::  dummy1, dummy2
  LOGICAL, INTENT(IN)     ::  vel_grid_yes, boundary_yes
  INTEGER, INTENT(IN)     ::  dir

  INTEGER, INTENT(OUT)    ::  node_id

  !--- Get the global number of nodes on which the values are stored
  CALL get_boundary_block_dims(nn1l, nn1u, nn2l, nn2u, nn3l, nn3u, vel_grid_yes, boundary_yes, dir)
  nn1 = nn1l + (NB1 - 2)*(N1 - 1) + nn1u
  nn2 = nn2l + (NB2 - 2)*(N2 - 1) + nn2u
  nn3 = nn3l + (NB3 - 2)*(N3 - 1) + nn3u

  !--- If in any direction only one block exists override dimension with local one
  IF (NB1 .EQ. 1) THEN
    CALL get_block_dims(nn1, dummy1, dummy2, vel_grid_yes, boundary_yes, dir)
  ENDIF
  IF (NB2 .EQ. 1) THEN
    CALL get_block_dims(dummy1, nn2, dummy2, vel_grid_yes, boundary_yes, dir)
  ENDIF
  IF (NB3 .EQ. 1) THEN
    CALL get_block_dims(dummy1, dummy2, nn3, vel_grid_yes, boundary_yes, dir)
  END IF

  !--- correct indices i,j,k locally
  IF (.NOT. vel_grid_yes) THEN
    !--- pressure grid indices ----------------------------------------------------------------------
    ! bbecsek: 
    ! if we want to convert the pressure grid indices into nodal IDs we have to
    ! consider that 
    ! if  BC_1U <= 0:  i = [1, N1-1]
    ! if  BC_1U >  0:  i = [1, N1  ]
    ! if  BC_2U <= 0:  j = [1, N2-1]
    ! if  BC_2U >  0:  j = [1, N2  ]
    ! if  BC_3U <= 0:  k = [1, N3-1]
    ! if  BC_3U >  0:  k = [1, N3  ]
    ! nodal IDs run from 0 to nn1xnn2xnn3-1
    ii  = i - 1
    jj  = j - 1
    kk  = k - 1

  ELSE
    !--- velocity grid indices ----------------------------------------------------------------------
    ! bbecsek: 
    ! if we want to convert the velocity grid indices into nodal IDs we have to
    ! consider that
    ! if we exclude boundary points, i.e. use sii and nii
    !   i = [1, N1-1]
    !   j = [1, N2-1]
    !   k = [1, N3-1]
    ! if we include boundary points, i.e. use siib and niib
    ! if  BC_1L <=0 and BC_1U <= 0:  i = [1, N1-1]
    ! if  BC_1L <=0 and BC_1U >  0:  i = [1, N1  ]
    ! if  BC_1L > 0 and BC_1U <= 0:  i = [0, N1-1]
    ! if  BC_1L > 0 and BC_1U >  0:  i = [0, N1  ]
    ! ... and analogously for the other spacial directions
    ! nodal IDs run from 0 to nn1xnn2xnn3-1
    IF ((BC_1L .LE. 0 .AND. dir .EQ. 1) .OR. (.NOT. boundary_yes .AND. dir .EQ. 1) .OR. (dir .NE. 1)) THEN
      ii = i - 1
    ELSE
      ii = i
    END IF

    IF ((BC_2L .LE. 0 .AND. dir .EQ. 2) .OR. (.NOT. boundary_yes .AND. dir .EQ. 2) .OR. (dir .NE. 2)) THEN
      jj = j - 1
    ELSE
      jj = j
    END IF

    IF ((BC_3L .LE. 0 .AND. dir .EQ. 3) .OR. (.NOT. boundary_yes .AND. dir .EQ. 3) .OR. (dir .NE. 3)) THEN
      kk = k - 1
    ELSE
      kk = k
    END IF

  END IF

  !--- This corrects the local indices to become global ones
  !--- 1-direction ---
  IF (iB(1,1) .GT. 1) THEN
    ii = ii + nn1l + (iB(1,1) - 2)*(N1 - 1)
  END IF
  !--- 2-direction ---
  IF (iB(2,1) .GT. 1) THEN
    jj = jj + nn2l + (iB(2,1) - 2)*(N2 - 1)
  END IF
  !--- 3-direction ---
  IF (iB(3,1) .GT. 1) THEN
    kk = kk + nn3l + (iB(3,1) - 2)*(N3 - 1)
  END IF

  node_id = ii + jj * nn1 + kk * nn1 * nn2

  END SUBROUTINE global_cart2id_glob_con


  SUBROUTINE global_id2cart_glob_con(node_id, i, j, k, vel_grid_yes, boundary_yes, dir)

  IMPLICIT NONE

  INTEGER, INTENT(IN)     ::  node_id
  INTEGER                 ::  ijn1
  INTEGER                 ::  nn1,nn2,nn3
  INTEGER                 ::  nn1l, nn1u, nn2l, nn2u, nn3l, nn3u
  INTEGER                 ::  dummy1, dummy2, dummy3
  LOGICAL, INTENT(IN)     ::  vel_grid_yes, boundary_yes
  INTEGER, INTENT(IN)     ::  dir

  INTEGER, INTENT(OUT)    ::  i,j,k

  !--- Get the global number of nodes on which the values are stored
  CALL get_boundary_block_dims(nn1l, nn1u, nn2l, nn2u, nn3l, nn3u, vel_grid_yes, boundary_yes, dir)
  nn1 = nn1l + (NB1 - 2)*(N1 - 1) + nn1u
  nn2 = nn2l + (NB2 - 2)*(N2 - 1) + nn2u
  nn3 = nn3l + (NB3 - 2)*(N3 - 1) + nn3u

  !--- If in any direction only one block exists override dimension with local one
  IF (NB1 .EQ. 1) THEN
    CALL get_block_dims(nn1, dummy1, dummy2, vel_grid_yes, boundary_yes, dir)
  ENDIF
  IF (NB2 .EQ. 1) THEN
    CALL get_block_dims(dummy1, nn2, dummy2, vel_grid_yes, boundary_yes, dir)
  ENDIF
  IF (NB3 .EQ. 1) THEN
    CALL get_block_dims(dummy1, dummy2, nn3, vel_grid_yes, boundary_yes, dir)
  END IF

  !IF (node_id .LT. nn1) THEN
  !  i = node_id
  !  j = 0
  !  k = 0
  !ELSE IF (node_id .GE. nn1 .AND. node_id .LT. nn1*nn2) THEN
  !  i = MOD(node_id, nn1)
  !  j = (node_id - i)/nn1
  !  k = 0
  !ELSE IF (node_id .GE. nn1*nn2) THEN
  !  ijn1 = MOD(node_id, nn1*nn2)
  !  i = MOD(ijn1, nn1)
  !  j = (ijn1 - i)/nn1
  !  k = (node_id - i - j*nn1)/(nn1*nn2)
  !END IF
  CALL id2cart(node_id, i, j, k, nn1, nn2, nn3)

  !--- This corrects the local indices to become global ones
  !--- 1-direction ---
  IF (iB(1,1) .GT. 1) THEN
    i = i - nn1l - (iB(1,1) - 2)*(N1 - 1)
  END IF
  !--- 2-direction ---
  IF (iB(2,1) .GT. 1) THEN
    j = j - nn2l - (iB(2,1) - 2)*(N2 - 1)
  END IF
  !--- 3-direction ---
  IF (iB(3,1) .GT. 1) THEN
    k = k - nn3l - (iB(3,1) - 2)*(N3 - 1)
  END IF

  !--- correct indices i,j,k locally
  IF (.NOT. vel_grid_yes) THEN
    i = i + 1
    j = j + 1
    k = k + 1
  ELSE
    IF ((BC_1L .LE. 0 .AND. dir .EQ. 1) .OR. (.NOT. boundary_yes .AND. dir .EQ. 1) .OR. (dir .NE. 1)) THEN
      i = i + 1
    END IF

    IF ((BC_2L .LE. 0 .AND. dir .EQ. 2) .OR. (.NOT. boundary_yes .AND. dir .EQ. 2) .OR. (dir .NE. 2)) THEN
      j = j + 1
    END IF

    IF ((BC_3L .LE. 0 .AND. dir .EQ. 3) .OR. (.NOT. boundary_yes .AND. dir .EQ. 3) .OR. (dir .NE. 3)) THEN
      k = k + 1
    END IF

  END IF

  END SUBROUTINE global_id2cart_glob_con


  ! Cartesian to index translation, only works for the pressure grid at the
  ! moment with periodic BCs
  SUBROUTINE global_cart2id_loc_con(node_id, i, j, k)

  IMPLICIT NONE

  INTEGER, INTENT(IN)     ::  i,j,k
  INTEGER                 ::  ii,jj,kk
  INTEGER                 ::  nn1,nn2,nn3
  INTEGER                 ::  nn_block, blockid
  INTEGER, INTENT(OUT)    ::  node_id

  nn1 = N1 - 1
  nn2 = N2 - 1
  nn3 = N3 - 1
  
  nn_block = nn1 * nn2 * nn3

  !--- Pressure grid indices
  ii = i - 1
  jj = j - 1
  kk = k - 1

  !--- global node id in this current block w/o consideration of other blocks
  node_id = ii + jj * nn1 + kk * nn1 * nn2

  !--- get the current block's id
  blockid = block_id(iB(1,1), iB(2,1), iB(3,1))

  !--- augment the node id to consider its location in the entire domain
  node_id = node_id + blockid * nn_block

  END SUBROUTINE global_cart2id_loc_con


  ! Cartesian to index translation, only works for the pressure grid at the
  ! moment
  ! buggy
  SUBROUTINE global_cart2id_loc_con2(node_id, i, j, k)

  IMPLICIT NONE

  INTEGER, INTENT(IN)     ::  i,j,k
  INTEGER                 ::  ii,jj,kk
  INTEGER                 ::  nn1,nn2,nn3
  INTEGER                 ::  nn1l, nn1u, nn2l, nn2u, nn3l, nn3u
  INTEGER                 ::  dummy1, dummy2
  INTEGER                 ::  nn_bulk_block
  INTEGER                 ::  nn_face_1l_block, nn_face_2l_block, nn_face_3l_block
  INTEGER                 ::  nn_face_1u_block, nn_face_2u_block, nn_face_3u_block
  INTEGER                 ::  nn_edge_1l2l_block, nn_edge_1u2l_block, nn_edge_1l2u_block, nn_edge_1u2u_block
  INTEGER                 ::  nn_edge_1l3l_block, nn_edge_1u3l_block, nn_edge_1l3u_block, nn_edge_1u3u_block
  INTEGER                 ::  nn_edge_2l3l_block, nn_edge_2u3l_block, nn_edge_2l3u_block, nn_edge_2u3u_block
  INTEGER                 ::  nn_corner_1l2l3l_block, nn_corner_1l2u3l_block, nn_corner_1l2l3u_block, nn_corner_1l2u3u_block
  INTEGER                 ::  nn_corner_1u2l3l_block, nn_corner_1u2u3l_block, nn_corner_1u2l3u_block, nn_corner_1u2u3u_block


  INTEGER, INTENT(OUT)    ::  node_id

  !--- Get the number of nodes on the blocks that lie at the boundary ---
  CALL get_boundary_block_dims(nn1l, nn1u, nn2l, nn2u, nn3l, nn3u, .FALSE., .FALSE., -1)

  !--- Compute block sizes for different types of process blocks ---
  !--- bulk blocks ---
  nn_bulk_block = (N1 - 1)*(N2 - 1)*(N3 - 1)
  !--- face blocks ---
  nn_face_1l_block =  nn1l   *(N2 - 1)*(N3 - 1)
  nn_face_1u_block =  nn1u   *(N2 - 1)*(N3 - 1)
  nn_face_2l_block = (N1 - 1)* nn2l   *(N3 - 1)
  nn_face_2u_block = (N1 - 1)* nn2u   *(N3 - 1)
  nn_face_3l_block = (N1 - 1)*(N2 - 1)* nn3l
  nn_face_3u_block = (N1 - 1)*(N2 - 1)* nn3u
  !--- edge blocks ---
  nn_edge_1l2l_block =  nn1l  * nn2l   *(N3 - 1)
  nn_edge_1u2l_block =  nn1u  * nn2l   *(N3 - 1)
  nn_edge_1l2u_block =  nn1l  * nn2u   *(N3 - 1)
  nn_edge_1u2u_block =  nn1u  * nn2u   *(N3 - 1)
  nn_edge_1l3l_block =  nn1l  *(N2 - 1)* nn3l
  nn_edge_1u3l_block =  nn1u  *(N2 - 1)* nn3l
  nn_edge_1l3u_block =  nn1l  *(N2 - 1)* nn3u
  nn_edge_1u3u_block =  nn1u  *(N2 - 1)* nn3u
  nn_edge_2l3l_block =(N1 - 1)* nn2l   * nn3l
  nn_edge_2u3l_block =(N1 - 1)* nn2u   * nn3l
  nn_edge_2l3u_block =(N1 - 1)* nn2l   * nn3u
  nn_edge_2u3u_block =(N1 - 1)* nn2u   * nn3u
  !--- corner blocks ---
  nn_corner_1l2l3l_block = nn1l * nn2l * nn3l
  nn_corner_1u2l3l_block = nn1u * nn2l * nn3l
  nn_corner_1l2u3l_block = nn1l * nn2u * nn3l
  nn_corner_1l2l3u_block = nn1l * nn2l * nn3u
  nn_corner_1u2u3l_block = nn1u * nn2u * nn3l
  nn_corner_1u2l3u_block = nn1u * nn2l * nn3u
  nn_corner_1l2u3u_block = nn1l * nn2u * nn3u
  nn_corner_1u2u3u_block = nn1u * nn2u * nn3u

  !--- pressure grid indices ---
  ii = i - 1
  jj = j - 1
  kk = k - 1

  !--- compute the local node ID ---
  CALL get_block_dims(nn1, nn2, nn3, .FALSE., .FALSE., -1)

  node_id = ii + jj * nn1 + kk * nn1 * nn2

  !--- now correct the local node_id to become global depending on in which block this is called ---
  IF (iB(1,1)*iB(2,1)*iB(3,1) .GT. 1) THEN
    node_id = node_id + nn_corner_1l2l3l_block
    IF (iB(3,1) .GT. 1) THEN
      !--- add corner blocks at lower 3 boundary (first corner was already added above) ---
      node_id = node_id + nn_corner_1u2l3l_block + nn_corner_1l2u3l_block + nn_corner_1u2u3l_block
      !--- add edge blocks at lower 3 boundary ---
      node_id = node_id + (NB1 - 2)*(nn_edge_2l3l_block + nn_edge_2u3l_block) &
                        + (NB2 - 2)*(nn_edge_1l3l_block + nn_edge_1u3l_block)
      !--- add face blocks at lower 3 boundary ---
      node_id = node_id + (NB1 - 2)*(NB2 - 2)*nn_face_3l_block
      !--- add slices in 3 plane that are not a boundary slice ---
      node_id = node_id + (iB(3,1) - 2)*(   nn_edge_1l2l_block + nn_edge_1u2l_block &
                                         +  nn_edge_1l2u_block + nn_edge_1u2u_block &
                                         + (NB1 - 2)*(NB2 - 2)*nn_bulk_block  )
      !--- add incomplete slice at current 3 direction position ---
      IF (iB(3,1) .LT. NB3) THEN
        IF (iB(2,1) .GT. 1) THEN
          !--- add lower 2 face along 1 direction ---
          node_id = node_id + nn_edge_1l2l_block + (NB1 - 2)*nn_face_2l_block + nn_edge_1u2l_block
          !--- add line in bulk along 1 direction ---
          node_id = node_id + (iB(2,1) - 2)*(nn_face_1l_block + (NB1 - 2)*nn_bulk_block + nn_face_1u_block)
          !--- add current line's blocks ---
          IF (iB(2,1) .LT. NB2) THEN
            IF (iB(1,1) .GT. 1) THEN
              node_id = node_id + nn_face_1l_block + (iB(1,1) - 2)*nn_bulk_block
            END IF
          END IF
          IF (iB(2,1) .EQ. NB2) THEN
            IF (iB(1,1) .GT. 1) THEN
              node_id = node_id + nn_edge_1l2u_block + (iB(1,1) - 2)*nn_face_2u_block
            END IF
          END IF
        ELSE ! iB(2,1) .EQ. 1
          IF (iB(1,1) .GT. 1) THEN
            node_id = node_id + nn_edge_1l2l_block + (iB(1,1) - 2)*nn_face_2l_block
          END IF
        END IF
      END IF
      !--- add upper 3 face ---
      IF (iB(3,1) .EQ. NB3) THEN
        IF (iB(2,1) .GT. 1) THEN
          !--- add lower 2 edge along 1 direction ---
          node_id = node_id + nn_corner_1l2l3u_block + (NB1 - 2)*nn_edge_2l3u_block + nn_corner_1u2l3u_block
          !--- add lines at upper 3 face along 1 direction ---
          node_id = node_id + (iB(2,1) - 2)*(nn_edge_1l3u_block + (NB1 - 2)*nn_face_3u_block + nn_edge_1u3u_block)
          !--- add current line's blocks ---
          IF (iB(2,1) .LT. NB2) THEN
            IF (iB(1,1) .GT. 1) THEN
              node_id = node_id + nn_edge_1l3u_block + (iB(1,1) - 2)*nn_face_3u_block
            END IF
          END IF
          IF (iB(2,1) .EQ. NB2) THEN
            IF (iB(1,1) .GT. 1) THEN
              node_id = node_id + nn_corner_1l2u3u_block + (iB(1,1) - 2)*nn_edge_2u3u_block
            END IF
          END IF
        ELSE ! iB(3,1) .EQ. NB3 .AND. iB(2,1) .EQ. 1
          IF (iB(1,1) .GT. 1) THEN
            node_id = node_id + nn_corner_1l2l3u_block + (iB(1,1) - 2)*nn_edge_2l3u_block
          END IF
        END IF
      END IF

    ELSE ! iB(3,1) .EQ. 1
      IF (iB(2,1) .GT. 1) THEN
        !--- add lower 2 edge along 1 direction ---
        node_id = node_id + (NB1 - 2)*nn_edge_2l3l_block + nn_corner_1u2l3l_block
        !--- add lines at lower 3 face along 1 direction ---
        node_id = node_id + (iB(2,1) - 2)*(nn_edge_1l3l_block + (NB1 - 2)*nn_face_3l_block + nn_edge_1u3l_block)
        !--- add current line's blocks ---
        IF (iB(2,1) .LT. NB2) THEN
          IF (iB(1,1) .GT. 1) THEN
            node_id = node_id + nn_edge_1l3l_block + (iB(1,1) - 2)*nn_face_3l_block
          END IF
        END IF
      !--- add upper 2 edge ---
        IF (iB(2,1) .EQ. NB2) THEN
          IF (iB(1,1) .GT. 1) THEN
            node_id = node_id + nn_corner_1l2u3l_block + (iB(1,1) - 2)*nn_edge_2u3l_block
          END IF
        END IF
      ELSE ! iB(3,1) .EQ. 1 .AND. iB(2,1) .EQ. 1
        IF (iB(1,1) .GT. 1) THEN
          node_id = node_id + (iB(1,1) - 2)*nn_edge_2l3l_block
        END IF
      END IF
    END IF

  END IF


  END SUBROUTINE global_cart2id_loc_con2


  SUBROUTINE global_id2cart_loc_con(node_id, i, j, k)

  IMPLICIT NONE

  INTEGER, INTENT(OUT)     ::  i,j,k
  INTEGER                  ::  ii,jj,kk
  INTEGER                  ::  nn1,nn2,nn3
  INTEGER                  ::  nn_block, blockid, nid
  INTEGER, INTENT(IN)      ::  node_id

  nn1 = N1 - 1
  nn2 = N2 - 1
  nn3 = N3 - 1

  nn_block = nn1 * nn2 * nn3

  !--- get the current block's id
  blockid = block_id(iB(1,1), iB(2,1), iB(3,1))

  !--- reduce the node id to be within this blocks ids
  nid = node_id - blockid * nn_block

  !--- regain indices
  CALL id2cart(nid, ii, jj, kk, nn1, nn2, nn3)

  !--- pressure grid indices
  i = ii + 1
  j = jj + 1
  k = kk + 1



  END SUBROUTINE global_id2cart_loc_con


  ! Index to cartesian translation, only works for the pressure grid at the
  ! moment
  SUBROUTINE global_id2cart_loc_con2(node_id, i, j, k)

  IMPLICIT NONE

  INTEGER, INTENT(OUT)     ::  i,j,k
  INTEGER                 ::  ii,jj,kk
  INTEGER                 ::  nn1,nn2,nn3
  INTEGER                 ::  nn1l, nn1u, nn2l, nn2u, nn3l, nn3u
  INTEGER                 ::  dummy1, dummy2
  INTEGER                 ::  nn_bulk_block
  INTEGER                 ::  nn_face_1l_block, nn_face_2l_block, nn_face_3l_block
  INTEGER                 ::  nn_face_1u_block, nn_face_2u_block, nn_face_3u_block
  INTEGER                 ::  nn_edge_1l2l_block, nn_edge_1u2l_block, nn_edge_1l2u_block, nn_edge_1u2u_block
  INTEGER                 ::  nn_edge_1l3l_block, nn_edge_1u3l_block, nn_edge_1l3u_block, nn_edge_1u3u_block
  INTEGER                 ::  nn_edge_2l3l_block, nn_edge_2u3l_block, nn_edge_2l3u_block, nn_edge_2u3u_block
  INTEGER                 ::  nn_corner_1l2l3l_block, nn_corner_1l2u3l_block, nn_corner_1l2l3u_block, nn_corner_1l2u3u_block
  INTEGER                 ::  nn_corner_1u2l3l_block, nn_corner_1u2u3l_block, nn_corner_1u2l3u_block, nn_corner_1u2u3u_block


  INTEGER, INTENT(IN)    ::  node_id
  INTEGER                ::  nid, ijn1

  nid = node_id

  !--- Get the number of nodes on the blocks that lie at the boundary ---
  CALL get_boundary_block_dims(nn1l, nn1u, nn2l, nn2u, nn3l, nn3u, .FALSE., .FALSE., -1)

  !--- Compute block sizes for different types of process blocks ---
  !--- bulk blocks ---
  nn_bulk_block = (N1 - 1)*(N2 - 1)*(N3 - 1)
  !--- face blocks ---
  nn_face_1l_block =  nn1l   *(N2 - 1)*(N3 - 1)
  nn_face_1u_block =  nn1u   *(N2 - 1)*(N3 - 1)
  nn_face_2l_block = (N1 - 1)* nn2l   *(N3 - 1)
  nn_face_2u_block = (N1 - 1)* nn2u   *(N3 - 1)
  nn_face_3l_block = (N1 - 1)*(N2 - 1)* nn3l
  nn_face_3u_block = (N1 - 1)*(N2 - 1)* nn3u
  !--- edge blocks ---
  nn_edge_1l2l_block =  nn1l  * nn2l   *(N3 - 1)
  nn_edge_1u2l_block =  nn1u  * nn2l   *(N3 - 1)
  nn_edge_1l2u_block =  nn1l  * nn2u   *(N3 - 1)
  nn_edge_1u2u_block =  nn1u  * nn2u   *(N3 - 1)
  nn_edge_1l3l_block =  nn1l  *(N2 - 1)* nn3l
  nn_edge_1u3l_block =  nn1u  *(N2 - 1)* nn3l
  nn_edge_1l3u_block =  nn1l  *(N2 - 1)* nn3u
  nn_edge_1u3u_block =  nn1u  *(N2 - 1)* nn3u
  nn_edge_2l3l_block =(N1 - 1)* nn2l   * nn3l
  nn_edge_2u3l_block =(N1 - 1)* nn2u   * nn3l
  nn_edge_2l3u_block =(N1 - 1)* nn2l   * nn3u
  nn_edge_2u3u_block =(N1 - 1)* nn2u   * nn3u
  !--- corner blocks ---
  nn_corner_1l2l3l_block = nn1l * nn2l * nn3l
  nn_corner_1u2l3l_block = nn1u * nn2l * nn3l
  nn_corner_1l2u3l_block = nn1l * nn2u * nn3l
  nn_corner_1l2l3u_block = nn1l * nn2l * nn3u
  nn_corner_1u2u3l_block = nn1u * nn2u * nn3l
  nn_corner_1u2l3u_block = nn1u * nn2l * nn3u
  nn_corner_1l2u3u_block = nn1l * nn2u * nn3u
  nn_corner_1u2u3u_block = nn1u * nn2u * nn3u

  !--- now correct the global node_id to become local depending on in which block this is called ---
  IF (iB(1,1)*iB(2,1)*iB(3,1) .GT. 1) THEN
    nid = nid - nn_corner_1l2l3l_block
    IF (iB(3,1) .GT. 1) THEN
      !--- remove corner blocks at lower 3 boundary (first corner was already removeed above) ---
      nid = nid - nn_corner_1u2l3l_block - nn_corner_1l2u3l_block - nn_corner_1u2u3l_block
      !--- remove edge blocks at lower 3 boundary ---
      nid = nid - (NB1 - 2)*(nn_edge_2l3l_block + nn_edge_2u3l_block) &
                        - (NB2 - 2)*(nn_edge_1l3l_block + nn_edge_1u3l_block)
      !--- remove face blocks at lower 3 boundary ---
      nid = nid - (NB1 - 2)*(NB2 - 2)*nn_face_3l_block
      !--- remove slices in 3 plane that are not a boundary slice ---
      nid = nid - (iB(3,1) - 2)*(   nn_edge_1l2l_block + nn_edge_1u2l_block &
                                         +  nn_edge_1l2u_block + nn_edge_1u2u_block &
                                         + (NB1 - 2)*(NB2 - 2)*nn_bulk_block  )
      !--- remove incomplete slice at current 3 direction position ---
      IF (iB(3,1) .LT. NB3) THEN
        IF (iB(2,1) .GT. 1) THEN
          !--- remove lower 2 face along 1 direction ---
          nid = nid - nn_edge_1l2l_block - (NB1 - 2)*nn_face_2l_block + nn_edge_1u2l_block
          !--- remove line in bulk along 1 direction ---
          nid = nid - (iB(2,1) - 2)*(nn_face_1l_block + (NB1 - 2)*nn_bulk_block + nn_face_1u_block)
          !--- remove current line's blocks ---
          IF (iB(2,1) .LT. NB2) THEN
            IF (iB(1,1) .GT. 1) THEN
              nid = nid - nn_face_1l_block - (iB(1,1) - 2)*nn_bulk_block
            END IF
          END IF
          IF (iB(2,1) .EQ. NB2) THEN
            IF (iB(1,1) .GT. 1) THEN
              nid = nid - nn_edge_1l2u_block - (iB(1,1) - 2)*nn_face_2u_block
            END IF
          END IF
        ELSE ! iB(2,1) .EQ. 1
          IF (iB(1,1) .GT. 1) THEN
            nid = nid - nn_edge_1l2l_block - (iB(1,1) - 2)*nn_face_2l_block
          END IF
        END IF
      END IF
      !--- remove upper 3 face ---
      IF (iB(3,1) .EQ. NB3) THEN
        IF (iB(2,1) .GT. 1) THEN
          !--- remove lower 2 edge along 1 direction ---
          nid = nid - nn_corner_1l2l3u_block - (NB1 - 2)*nn_edge_2l3u_block - nn_corner_1u2l3u_block
          !--- remove lines at upper 3 face along 1 direction ---
          nid = nid - (iB(2,1) - 2)*(nn_edge_1l3u_block + (NB1 - 2)*nn_face_3u_block + nn_edge_1u3u_block)
          !--- remove current line's blocks ---
          IF (iB(2,1) .LT. NB2) THEN
            IF (iB(1,1) .GT. 1) THEN
              nid = nid - nn_edge_1l3u_block - (iB(1,1) - 2)*nn_face_3u_block
            END IF
          END IF
          IF (iB(2,1) .EQ. NB2) THEN
            IF (iB(1,1) .GT. 1) THEN
              nid = nid - nn_corner_1l2u3u_block - (iB(1,1) - 2)*nn_edge_2u3u_block
            END IF
          END IF
        ELSE ! iB(3,1) .EQ. NB3 .AND. iB(2,1) .EQ. 1
          IF (iB(1,1) .GT. 1) THEN
            nid = nid - nn_corner_1l2l3u_block - (iB(1,1) - 2)*nn_edge_2l3u_block
          END IF
        END IF
      END IF

    ELSE ! iB(3,1) .EQ. 1
      IF (iB(2,1) .GT. 1) THEN
        !--- remove lower 2 edge along 1 direction ---
        nid = nid - (NB1 - 2)*nn_edge_2l3l_block - nn_corner_1u2l3l_block
        !--- remove lines at lower 3 face along 1 direction ---
        nid = nid - (iB(2,1) - 2)*(nn_edge_1l3l_block + (NB1 - 2)*nn_face_3l_block + nn_edge_1u3l_block)
        !--- remove current line's blocks ---
        IF (iB(2,1) .LT. NB2) THEN
          IF (iB(1,1) .GT. 1) THEN
            nid = nid - nn_edge_1l3l_block - (iB(1,1) - 2)*nn_face_3l_block
          END IF
        END IF
      !--- remove upper 2 edge ---
        IF (iB(2,1) .EQ. NB2) THEN
          IF (iB(1,1) .GT. 1) THEN
            nid = nid - nn_corner_1l2u3l_block - (iB(1,1) - 2)*nn_edge_2u3l_block
          END IF
        END IF
      ELSE ! iB(3,1) .EQ. 1 .AND. iB(2,1) .EQ. 1
        IF (iB(1,1) .GT. 1) THEN
          nid = nid - (iB(1,1) - 2)*nn_edge_2l3l_block
        END IF
      END IF
    END IF

  END IF
  
  !--- compute the local node indices ---
  CALL get_block_dims(nn1, nn2, nn3, .FALSE., .FALSE., -1)

  !IF (nid .LT. nn1) THEN
  !  i = nid
  !  j = 0
  !  k = 0
  !ELSE IF (nid .GE. nn1 .AND. nid .LT. nn1*nn2) THEN
  !  i = MOD(nid, nn1)
  !  j = (nid - i)/nn1
  !  k = 0
  !ELSE IF (nid .GE. nn1*nn2) THEN
  !  ijn1 = MOD(nid, nn1*nn2)
  !  i = MOD(ijn1, nn1)
  !  j = (ijn1 - i)/nn1
  !  k = (nid - i - j*nn1)/(nn1*nn2)
  !END IF
  CALL id2cart(nid, i, j, k, nn1, nn2, nn3)

  !--- pressure grid indices ---
  i = i + 1
  j = j + 1
  k = k + 1

  END SUBROUTINE global_id2cart_loc_con2


  FUNCTION n_local_nodes() RESULT (fn_val)

  IMPLICIT NONE
  INTEGER                :: nn1, nn2, nn3
  INTEGER                :: fn_val

  CALL get_block_dims(nn1, nn2, nn3, .FALSE., .FALSE., -1)

  fn_val =  nn1 * nn2 * nn3


  END FUNCTION n_local_nodes

  FUNCTION n_global_nodes() RESULT(fn_val)

  IMPLICIT NONE
  INTEGER                :: nn1, nn2, nn3
  INTEGER                :: fn_val

  CALL get_block_dims(nn1, nn2, nn3, .FALSE., .FALSE., -1)

  fn_val = nn1 * nn2 * nn3 * NB1 * NB2 * NB3
  
  END FUNCTION n_global_nodes


  !> local node IDs with one row of ghost cells included
  SUBROUTINE local_cart2id_loc_con(node_id,i,j,k,n_local,n_ghost)

  IMPLICIT NONE

  INTEGER, INTENT(IN)    :: i, j, k
  INTEGER, INTENT(OUT)   :: node_id
  INTEGER, INTENT(OUT)   :: n_local, n_ghost
  INTEGER                :: nn1, nn2, nn3
  INTEGER                :: ii, jj, kk

  CALL get_block_dims(nn1,nn2,nn3,.FALSE.,.FALSE.,-1)
  n_local = nn1*nn2*nn3

  !--- correct for rows of ghost vertices (one on each side) ---
  !--- 1-direction
  IF (BC_1L .LE. 0) THEN
    nn1 = nn1 + 1
    ii  = i
  ELSE
    ii  = i - 1
  END IF
  IF (BC_1U .LE. 0) THEN
!    nn1 = nn1 + 1
  END IF
  !--- 2-direction
  IF (BC_2L .LE. 0) THEN
    nn2 = nn2 + 1
    jj  = j
  ELSE
    jj  = j - 1
  END IF
  IF (BC_2U .LE. 0) THEN
!    nn2 = nn2 + 1
  END IF
  !--- 3-direction
  IF (dimens .EQ. 3) THEN
    IF (BC_3L .LE. 0) THEN
      nn3 = nn3 + 1
      kk  = k
    ELSE
      kk  = k - 1
    END IF
    IF (BC_3U .LE. 0) THEN
!      nn3 = nn3 + 1
    END IF
  ELSE
    nn3 = nn3
    kk  = 0
  END IF

  !--- calculate number of ghost nodes ---
  n_ghost = nn1*nn2*nn3 - n_local

  !--- compute node id ---
  node_id = ii + jj*nn1 + kk*nn1*nn2


  END SUBROUTINE local_cart2id_loc_con


  SUBROUTINE local_id2cart_loc_con(node_id,i,j,k,n_local,n_ghost)

  IMPLICIT NONE

  INTEGER, INTENT(OUT)   :: i, j, k
  INTEGER, INTENT(IN)    :: node_id
  INTEGER, INTENT(OUT)   :: n_local, n_ghost
  INTEGER                :: nn1, nn2, nn3
  INTEGER                :: ijn1

  CALL get_block_dims(nn1,nn2,nn3,.FALSE.,.FALSE.,-1)
  n_local = nn1*nn2*nn3

  !--- correct for rows of ghost vertices (one on each side) ---
  !--- 1-direction
  IF (BC_1L .LE. 0) THEN
    nn1 = nn1 + 1
  END IF
  IF (BC_1U .LE. 0) THEN
!    nn1 = nn1 + 1
  END IF
  !--- 2-direction
  IF (BC_2L .LE. 0) THEN
    nn2 = nn2 + 1
  END IF
  IF (BC_2U .LE. 0) THEN
!    nn2 = nn2 + 1
  END IF
  !--- 3-direction
  IF (dimens .EQ. 3) THEN
    IF (BC_3L .LE. 0) THEN
      nn3 = nn3 + 1
    END IF
    IF (BC_3U .LE. 0) THEN
!      nn3 = nn3 + 1
    END IF
  ELSE
    nn3 = nn3
  END IF

  !--- calculate number of ghost nodes ---
  n_ghost = nn1*nn2*nn3 - n_local

  !--- convert back ---
  !IF (node_id .LT. nn1) THEN
  !  i = node_id
  !  j = 0
  !  k = 0
  !ELSE IF (node_id .GE. nn1 .AND. node_id .LT. nn1*nn2) THEN
  !  i = MOD(node_id, nn1)
  !  j = (node_id - i)/nn1
  !  k = 0
  !ELSE IF (node_id .GE. nn1*nn2) THEN
  !  ijn1 = MOD(node_id, nn1*nn2)
  !  i = MOD(ijn1, nn1)
  !  j = (ijn1 - i)/nn1
  !  k = (node_id - i - j*nn1)/(nn1*nn2)
  !END IF
  CALL id2cart(node_id, i, j, k, nn1, nn2, nn3)

  !--- correct indices for rows of ghostcells ---
  IF (BC_1L .GT. 0) THEN
    i = i + 1
  END IF
  IF (BC_2L .GT. 0) THEN
    j = j + 1
  END IF
  IF (BC_3L .GT. 0 .AND. dimens .EQ.  3) THEN
    k = k + 1
  END IF

  IF (dimens .NE. 3) k = 1

  END SUBROUTINE local_id2cart_loc_con




  SUBROUTINE local_elem_base_loc_con(node_id,i,j,k)

  IMPLICIT NONE

  INTEGER, INTENT(OUT)   :: i, j, k
  INTEGER, INTENT(IN)    :: node_id
  INTEGER                :: nn1, nn2, nn3
  INTEGER                :: ijn1

  CALL get_block_dims(nn1,nn2,nn3,.FALSE.,.FALSE.,-1)

  !--- correct for 1 rows of ghost vertices (one on each side) ---
  ! per default we want one less node per direction as the base i
  ! (POTENTIALLY DEPRECATED)
  ! nn1 = nn1 - 1
  ! nn2 = nn2 - 1
  ! nn3 = nn3 - 1
  ! per default we want two less nodes per direction to avoid duplicate
  ! elements across process interfaces
  nn1 = nn1 - 2
  nn2 = nn2 - 2
  nn3 = nn3 - 2
  IF (BC_1L .LE. 0) THEN
    nn1 = nn1 + 1
  END IF
  IF (BC_1U .LE. 0) THEN
    nn1 = nn1 + 1
  END IF
  IF (BC_2L .LE. 0) THEN
    nn2 = nn2 + 1
  END IF
  IF (BC_2U .LE. 0) THEN
    nn2 = nn2 + 1
  END IF
  IF (BC_3L .LE. 0) THEN
    nn3 = nn3 + 1
  END IF
  IF (BC_3U .LE. 0) THEN
    nn3 = nn3 + 1
  END IF


  !--- convert back ---
  !IF (node_id .LT. nn1) THEN
  !  i = node_id
  !  j = 0
  !  k = 0
  !ELSE IF (node_id .GE. nn1 .AND. node_id .LT. nn1*nn2) THEN
  !  i = MOD(node_id, nn1)
  !  j = (node_id - i)/nn1
  !  k = 0
  !ELSE IF (node_id .GE. nn1*nn2) THEN
  !  ijn1 = MOD(node_id, nn1*nn2)
  !  i = MOD(ijn1, nn1)
  !  j = (ijn1 - i)/nn1
  !  k = (node_id - i - j*nn1)/(nn1*nn2)
  !END IF
  CALL id2cart(node_id, i, j, k, nn1, nn2, nn3)

  !--- correct indices for rows of ghostcells ---
  IF (BC_1L .GT. 0) THEN
    i = i + 1
  END IF
  IF (BC_2L .GT. 0) THEN
    j = j + 1
  END IF
  IF (BC_3L .GT. 0) THEN
    k = k + 1
  END IF

  IF (dimens .NE. 3) k = 1

  END SUBROUTINE local_elem_base_loc_con



  SUBROUTINE id2cart(id, ii, jj, kk, nn1, nn2, nn3)

  IMPLICIT NONE

  INTEGER, INTENT(IN)   :: id
  INTEGER, INTENT(OUT)  :: ii, jj, kk
  INTEGER, INTENT(IN)   :: nn1, nn2, nn3
  INTEGER               :: ijn1


  IF (id .LT. nn1) THEN
    ii = id
    jj = 0
    kk = 0
  ELSE IF (id .GE. nn1 .AND. id .LT. nn1*nn2) THEN
    ii = MOD(id, nn1)
    jj = (id - ii)/nn1
    kk = 0
  ELSE IF (id .GE. nn1*nn2) THEN
    ijn1 = MOD(id, nn1*nn2)
    ii = MOD(ijn1, nn1)
    jj = (ijn1 - ii)/nn1
    kk = (id - ii - jj*nn1)/(nn1*nn2)
  END IF

  END SUBROUTINE id2cart


  !> returns the local node_id to the node_no-th node of elem_no element
  SUBROUTINE local_tet_element_nodes(elem_no, node_no, node_id)

  IMPLICIT NONE

  INTEGER, INTENT(IN)    :: elem_no, node_no
  INTEGER, INTENT(OUT)   :: node_id
  INTEGER                :: base_node_id, tet_type
  INTEGER                :: ii,jj,kk
  INTEGER                :: n_local, n_ghost
  INTEGER                :: n_elem

  !--- check whether the entered elem_no is actually within this process ---
  n_elem = total_n_local_tet_elements()
  
  IF (elem_no .GE. n_elem) THEN
      WRITE(*,*)"WARNING: requested tetrahedral element out of range on process ",rank,". Aborting ..."
      CALL MPI_FINALIZE(merror)
      STOP
  END IF

  !--- compute the type of tetrahedron inside the hexahedron ---
  tet_type = MOD(elem_no,6)
  
  !--- now that we know the tet-type in [0,5], compute the base node ID ---
  base_node_id = (elem_no - tet_type)/6
  
  !--- convert the local base_node_id to local i,j,k ---
  CALL local_elem_base_loc_con(base_node_id,ii,jj,kk)

  !--- based on the tet_type compute the node_no node's node_id ---
  SELECT CASE (tet_type)
    CASE (0)
      SELECT CASE (node_no)
        CASE (0)
          ii = ii
          jj = jj
          kk = kk
        CASE (1)
          ii = ii + 1
          jj = jj + 1
          kk = kk + 1
        CASE (2)
          ii = ii
          jj = jj
          kk = kk + 1
        CASE (3)
          ii = ii
          jj = jj + 1
          kk = kk + 1
        CASE DEFAULT
          !--- check whether requested element node is existing ---
          WRITE(*,*)"WARNING: requested tetrahedron node out of range on process",rank,". Aborting ..."
          CALL MPI_FINALIZE(merror)
          STOP
      END SELECT
    CASE (1)
      SELECT CASE (node_no)
        CASE (0)
          ii = ii
          jj = jj
          kk = kk
        CASE (1)
          ii = ii + 1
          jj = jj + 1
          kk = kk + 1
        CASE (2)
          ii = ii
          jj = jj + 1
          kk = kk + 1
        CASE (3)
          ii = ii
          jj = jj + 1
          kk = kk
        CASE DEFAULT
          !--- check whether requested element node is existing ---
          WRITE(*,*)"WARNING: requested tetrahedron node out of range on process",rank,". Aborting ..."
          CALL MPI_FINALIZE(merror)
          STOP
      END SELECT
    CASE (2)
      SELECT CASE (node_no)
        CASE (0)
          ii = ii
          jj = jj
          kk = kk
        CASE (1)
          ii = ii + 1
          jj = jj
          kk = kk + 1
        CASE (2)
          ii = ii
          jj = jj
          kk = kk + 1
        CASE (3)
          ii = ii + 1
          jj = jj + 1
          kk = kk + 1
        CASE DEFAULT
          !--- check whether requested element node is existing ---
          WRITE(*,*)"WARNING: requested tetrahedron node out of range on process",rank,". Aborting ..."
          CALL MPI_FINALIZE(merror)
          STOP
      END SELECT
    CASE (3)
      SELECT CASE (node_no)
        CASE (0)
          ii = ii
          jj = jj
          kk = kk
        CASE (1)
          ii = ii + 1
          jj = jj + 1
          kk = kk
        CASE (2)
          ii = ii + 1
          jj = jj + 1
          kk = kk + 1
        CASE (3)
          ii = ii
          jj = jj + 1
          kk = kk
        CASE DEFAULT
          !--- check whether requested element node is existing ---
          WRITE(*,*)"WARNING: requested tetrahedron node out of range on process",rank,". Aborting ..."
          CALL MPI_FINALIZE(merror)
          STOP
      END SELECT
    CASE (4)
      SELECT CASE (node_no)
        CASE (0)
          ii = ii
          jj = jj
          kk = kk
        CASE (1)
          ii = ii + 1
          jj = jj
          kk = kk
        CASE (2)
          ii = ii + 1
          jj = jj + 1
          kk = kk + 1
        CASE (3)
          ii = ii + 1
          jj = jj + 1
          kk = kk
        CASE DEFAULT
          !--- check whether requested element node is existing ---
          WRITE(*,*)"WARNING: requested tetrahedron node out of range on process",rank,". Aborting ..."
          CALL MPI_FINALIZE(merror)
          STOP
      END SELECT
    CASE (5)
      SELECT CASE (node_no)
        CASE (0)
          ii = ii
          jj = jj
          kk = kk
        CASE (1)
          ii = ii + 1
          jj = jj
          kk = kk
        CASE (2)
          ii = ii + 1
          jj = jj
          kk = kk + 1
        CASE (3)
          ii = ii + 1
          jj = jj + 1
          kk = kk + 1
        CASE DEFAULT
          !--- check whether requested element node is existing ---
          WRITE(*,*)"WARNING: requested tetrahedron node out of range on process",rank,". Aborting ..."
          CALL MPI_FINALIZE(merror)
          STOP
      END SELECT
    CASE DEFAULT
      WRITE(*,*)"WARNING: invalid tetrahedron type encountered. Aborting..."
      CALL MPI_FINALIZE(merror)
      STOP
  END SELECT

  !--- convert new i,j,k back to node_id ---
  CALL local_cart2id_loc_con(node_id, ii, jj, kk, n_local, n_ghost)

  END SUBROUTINE local_tet_element_nodes


  FUNCTION local_to_global_tet_elem_id(loc_elem_id) RESULT(fn_val)

  IMPLICIT NONE

  INTEGER, INTENT(IN)    :: loc_elem_id
  INTEGER                :: fn_val

  fn_val = loc_elem_id + block_id(iB(1,1), iB(2,1), iB(3,1)) &
                       * total_n_local_tet_elements()

  END FUNCTION local_to_global_tet_elem_id


  FUNCTION local_to_global_tri_elem_id(loc_elem_id) RESULT(fn_val)

  IMPLICIT NONE

  INTEGER, INTENT(IN)    :: loc_elem_id
  INTEGER                :: fn_val

  fn_val = loc_elem_id + block_id(iB(1,1), iB(2,1), iB(3,1)) &
                       * total_n_local_tri_elements()

  END FUNCTION local_to_global_tri_elem_id



  SUBROUTINE local_tri_element_nodes(elem_no, node_no, node_id)

  IMPLICIT NONE

  INTEGER, INTENT(IN)    :: elem_no, node_no
  INTEGER, INTENT(OUT)   :: node_id
  INTEGER                :: base_node_id, tri_type
  INTEGER                :: ii,jj,kk
  INTEGER                :: n_elem
  INTEGER                :: n_local, n_ghost

  !--- check whether the entered elem_no is actually within this process ---
  n_elem = total_n_local_tri_elements()

  IF (elem_no .GE. n_elem) THEN
    WRITE(*,*)"WARNING: requested triangular element out of range on process ",rank,". Aborting ..."
    CALL MPI_FINALIZE(merror)
    STOP
  END IF

  !--- compute the type of triangle inside the quaadrilateral ---
  tri_type = MOD(elem_no,2)

  !--- now that we know the tri-type in [0,1], compute the base node ID ---
  base_node_id = (elem_no - tri_type)/2

  !--- convert the local base_node_id to local i,j,k ---
  CALL local_elem_base_loc_con(base_node_id,ii,jj,kk)

  !--- based on the tri_type compute the node_no node's node_id ---
  SELECT CASE (tri_type)
    CASE (0)
      SELECT CASE (node_no)
        CASE (0)
          ii = ii
          jj = jj
        CASE (1)
          ii = ii + 1
          jj = jj + 1
        CASE (2)
          ii = ii
          jj = jj + 1
        CASE DEFAULT
          !--- check whether requested element node is existing ---
          WRITE(*,*)"WARNING: requested triangle node out of range on process",rank,". Aborting ..."
          CALL MPI_FINALIZE(merror)
          STOP
      END SELECT
    CASE (1)
      SELECT CASE (node_no)
        CASE (0)
          ii = ii
          jj = jj
        CASE (1)
          ii = ii + 1
          jj = jj
        CASE (2)
          ii = ii + 1
          jj = jj + 1
        CASE DEFAULT
          !--- check whether requested element node is existing ---
          WRITE(*,*)"WARNING: requested triangle node out of range on process",rank,". Aborting ..."
          CALL MPI_FINALIZE(merror)
          STOP
      END SELECT
    CASE DEFAULT
      WRITE(*,*)"WARNING: invalid triangle type encountered. Aborting..."
      CALL MPI_FINALIZE(merror)
      STOP
  END SELECT

  !--- convert new i,j,k back to node_id ---
  CALL local_cart2id_loc_con(node_id, ii, jj, kk, n_local, n_ghost)

  END SUBROUTINE local_tri_element_nodes



  FUNCTION total_n_local_tri_elements() RESULT(fn_val)

  IMPLICIT NONE

  INTEGER                :: fn_val
  INTEGER                :: nn1, nn2, nn3

  IF (dimens .EQ. 3) THEN
    WRITE(*,*)"WARNING: Tri elements only supported for 2D simulations. Aborting ..."
    CALL MPI_FINALIZE(merror)
    STOP
  END IF



  CALL get_block_dims(nn1,nn2,nn3,.FALSE.,.FALSE.,-1)

  !--- correct for rows of ghost vertices (one on each side) ---
  !--- 1-direction ---
  IF (BC_1L .LE. 0) THEN
    nn1 = nn1 + 1
  END IF
  IF (BC_1U .LE. 0) THEN
    nn1 = nn1 + 1
  END IF
  !--- 2-direction ---
  IF (BC_2L .LE. 0) THEN
    nn2 = nn2 + 1
  END IF
  IF (BC_2U .LE. 0) THEN
    nn2 = nn2 + 1
  END IF

  !--- compute number of elements ---
  ! minus 2 to avoid duplicate elements across process boundaries
  fn_val = (nn1-2)*(nn2-2)*2

  END FUNCTION total_n_local_tri_elements



  FUNCTION total_n_local_tet_elements() RESULT(fn_val)

  IMPLICIT NONE

  INTEGER                :: fn_val
  INTEGER                :: nn1, nn2, nn3

  IF (dimens .EQ. 2) THEN
    WRITE(*,*)"WARNING: Tet elements only supported for 3D simulations. Aborting ..."
    CALL MPI_FINALIZE(merror)
    STOP
  END IF


  CALL get_block_dims(nn1,nn2,nn3,.FALSE.,.FALSE.,-1)

  !--- correct for rows of ghost vertices (one on each side) ---
  !--- 1-direction
  IF (BC_1L .LE. 0) THEN
    nn1 = nn1 + 1
  END IF
  IF (BC_1U .LE. 0) THEN
    nn1 = nn1 + 1
  END IF
  !--- 2-direction
  IF (BC_2L .LE. 0) THEN
    nn2 = nn2 + 1
  END IF
  IF (BC_2U .LE. 0) THEN
    nn2 = nn2 + 1
  END IF
  !--- 3-direction
  IF (BC_3L .LE. 0) THEN
    nn3 = nn3 + 1
  END IF
  IF (BC_3U .LE. 0) THEN
    nn3 = nn3 + 1
  END IF

  !--- compute number of elements ---
  ! minus 2 to avoid duplicate elements across process boundaries
  fn_val = (nn1-2)*(nn2-2)*(nn3-2)*6


  END FUNCTION total_n_local_tet_elements


  FUNCTION total_n_global_tri_elements() RESULT(fn_val)

  IMPLICIT NONE

  INTEGER                :: fn_val
  INTEGER                :: n_local_tri_elems

  fn_val = total_n_local_tri_elements() * NB1 * NB2

  END FUNCTION total_n_global_tri_elements


  FUNCTION total_n_global_tet_elements() RESULT(fn_val)

  IMPLICIT NONE

  INTEGER                :: fn_val
  INTEGER                :: n_local_tet_elems

  fn_val = total_n_local_tet_elements() * NB1 * NB2 * NB3

  END FUNCTION total_n_global_tet_elements



  !> The stuff in here should be replaced by comuting the dims
  !! through the start and end indices for the respective grids, duh.
  SUBROUTINE get_block_dims(nn1, nn2, nn3, vel_grid_yes, boundary_yes, dir)

  IMPLICIT NONE

  INTEGER, INTENT(OUT)   ::  nn1, nn2, nn3
  LOGICAL, INTENT(IN )   ::  vel_grid_yes, boundary_yes
  INTEGER, INTENT(IN )   ::  dir

  IF (.NOT. vel_grid_yes) THEN
    !--- pressure grid indices ----------------------------------------------------------------------
    ! bbecsek: 
    ! if we want to convert the pressure grid indices into nodal IDs we have to
    ! consider that 
    ! if  BC_1U <= 0:  i = [1, N1-1]
    ! if  BC_1U >  0:  i = [1, N1  ]
    ! if  BC_2U <= 0:  j = [1, N2-1]
    ! if  BC_2U >  0:  j = [1, N2  ]
    ! if  BC_3U <= 0:  k = [1, N3-1]
    ! if  BC_3U >  0:  k = [1, N3  ]
    ! nodal IDs run from 0 to nn1xnn2xnn3-1
    IF (BC_1U .GT. 0) THEN
      nn1 = N1
    ELSE
      nn1 = N1 - 1
    END IF

    IF (BC_2U .GT. 0) THEN
      nn2 = N2
    ELSE
      nn2 = N2 - 1
    END IF

    IF (BC_3U .GT. 0) THEN
      nn3 = N3
    ELSE
      nn3 = N3 - 1
    END IF

  ELSE
    !--- velocity grid indices ----------------------------------------------------------------------
    ! bbecsek:
    ! if we want to convert the velocity grid indices into nodal IDs we have to
    ! consider that
    ! if we exclude boundary points, i.e. use sii and nii
    !   i = [1, N1-1]
    !   j = [1, N2-1]
    !   k = [1, N3-1]
    ! if we include boundary points, i.e. use siib and niib
    ! if  BC_1L <=0 and BC_1U <= 0:  i = [1, N1-1]
    ! if  BC_1L <=0 and BC_1U >  0:  i = [1, N1  ]
    ! if  BC_1L > 0 and BC_1U <= 0:  i = [0, N1-1]
    ! if  BC_1L > 0 and BC_1U >  0:  i = [0, N1  ]
    ! ... and analogously for the other spacial directions
    ! nodal IDs run from 0 to nn1xnn2xnn3-1
    IF      ((BC_1L .LE. 0 .AND. BC_1U .LE. 0 .AND. dir .EQ. 1) .OR. (BC_1U .LE. 0 .AND. dir .NE. 1) .OR. (.NOT. boundary_yes .AND. dir .EQ. 1)) THEN
      nn1 = N1 - 1
    ELSE IF ( BC_1L .GT. 0 .AND. BC_1U .GT. 0 .AND. dir .EQ. 1) THEN
      nn1 = N1 + 1
    ELSE
      nn1 = N1
    END IF

    IF      ((BC_2L .LE. 0 .AND. BC_2U .LE. 0 .AND. dir .EQ. 2) .OR. (BC_2U .LE. 0 .AND. dir .NE. 2) .OR. (.NOT. boundary_yes .AND. dir .EQ. 2)) THEN
      nn2 = N2 - 1
    ELSE IF ( BC_2L .GT. 0 .AND. BC_2U .GT. 0 .AND. dir .EQ. 2) THEN
      nn2 = N2 + 1
    ELSE
      nn2 = N2
    END IF
    
    IF      ((BC_3L .LE. 0 .AND. BC_3U .LE. 0 .AND. dir .EQ. 3) .OR. (BC_3U .LE. 0 .AND. dir .NE. 3) .OR. (.NOT. boundary_yes .AND. dir .EQ. 3)) THEN
      nn3 = N3 - 1
    ELSE IF ( BC_3L .GT. 0 .AND. BC_3U .GT. 0 .AND. dir .EQ. 3) THEN
      nn3 = N3 + 1
    ELSE
      nn3 = N3
    END IF

  END IF

  END SUBROUTINE get_block_dims

  SUBROUTINE get_boundary_block_dims(nn1l, nn1u, nn2l, nn2u, nn3l, nn3u, vel_grid_yes, boundary_yes, dir)

  IMPLICIT NONE

  INTEGER, INTENT(IN )   ::  dir
  LOGICAL, INTENT(IN )   ::  vel_grid_yes, boundary_yes
  
  INTEGER, INTENT(OUT)   ::  nn1l, nn1u, nn2l, nn2u, nn3l, nn3u

  IF (.NOT. vel_grid_yes) THEN
    !--- 1-direction ---
    nn1l = N1 - 1
    IF (BC_1U_global .GT. 0) THEN
      nn1u = N1
    ELSE 
      nn1u = N1 - 1
    END IF
    !--- 2-direction ---
    nn2l = N2 - 1
    IF (BC_2U_global .GT. 0) THEN
      nn2u = N2
    ELSE 
      nn2u = N2 - 1
    END IF
    !--- 3-direction ---
    nn3l = N3 - 1
    IF (BC_3U_global .GT. 0) THEN
      nn3u = N3
    ELSE 
      nn3u = N3 - 1
    END IF
    !-------------------
  ELSE
    !--- 1-direction ---
    !--- lowest block ---
    IF (BC_1L_global .GT. 0 .AND. dir .EQ. 1) THEN
      nn1l = N1
    ELSE
      nn1l = N1 - 1
    END IF
    !--- topmost block ---
    IF ((.NOT. boundary_yes .AND. BC_1U_global .GT. 0 .AND. dir .EQ. 1) .OR. (BC_1U_global .GT. 0 .AND. dir .NE. 1)) THEN
      nn1u = N1
    ELSE
      nn1u = N1 - 1
    END IF
    !--- 2-direction ---
    !--- lowest block ---
    IF (BC_2L_global .GT. 0 .AND. dir .EQ. 2) THEN
      nn2l = N2
    ELSE
      nn2l = N2 - 1
    END IF
    !--- topmost block ---
    IF ((.NOT. boundary_yes .AND. BC_2U_global .GT. 0 .AND. dir .EQ. 2) .OR. (BC_2U_global .GT. 0 .AND. dir .NE. 2)) THEN
      nn2u = N2
    ELSE
      nn2u = N2 - 1
    END IF
    !--- 3-direction ---
    !--- lowest block ---
    IF (BC_3L_global .GT. 0 .AND. dir .EQ. 3) THEN
      nn3l = N3
    ELSE
      nn3l = N3 - 1
    END IF
    !--- topmost block ---
    IF ((.NOT. boundary_yes .AND. BC_3U_global .GT. 0 .AND. dir .EQ. 3) .OR. (BC_3U_global .GT. 0 .AND. dir .NE. 3)) THEN
      nn3u = N3
    ELSE
      nn3u = N3 - 1
    END IF
    !---------------------
  END IF

  END SUBROUTINE get_boundary_block_dims


  SUBROUTINE interpolate_force_pre_vel(dir)
 
  IMPLICIT NONE
  
  INTEGER, INTENT(IN)    :: dir
  REAL(8)                   :: inter(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))

  inter = 0.

  IF (dir .EQ. 1) THEN
    ! exchange across boundaries on pressure grid
    ! also updates ghost cells
    CALL exchange(1,0,fd(b1L,b2L,b3L,1))
    CALL exchange(2,0,fd(b1L,b2L,b3L,1))
    CALL exchange(3,0,fd(b1L,b2L,b3L,1))
    ! Force comp 1 onto vel 1 grid
    CALL interpolate2_pre_vel(.FALSE., 1, fd(b1L,b2L,b3L,1), inter)
    fd(:,:,:,1) = inter(:,:,:)
    ! exchange across boundaries on vel grid
    CALL exchange(1,1,fd(b1L,b2L,b3L,1))
    CALL exchange(2,1,fd(b1L,b2L,b3L,1))
    CALL exchange(3,1,fd(b1L,b2L,b3L,1))
  END IF

  IF (dir .EQ. 2) THEN
    ! exchange across boundaries on pressure grid
    ! also updates ghost cells
    CALL exchange(1,0,fd(b1L,b2L,b3L,2))
    CALL exchange(2,0,fd(b1L,b2L,b3L,2))
    CALL exchange(3,0,fd(b1L,b2L,b3L,2))
    ! Force comp 2 onto vel 2 grid
    CALL interpolate2_pre_vel(.FALSE., 2, fd(b1L,b2L,b3L,2), inter)
    fd(:,:,:,2) = inter(:,:,:)
    ! exchange across boundaries on vel grid
    CALL exchange(1,2,fd(b1L,b2L,b3L,2))
    CALL exchange(2,2,fd(b1L,b2L,b3L,2))
    CALL exchange(3,2,fd(b1L,b2L,b3L,2))
  END IF

  IF (dir .EQ. 3) THEN
    ! exchange across boundaries on pressure grid
    ! also updates ghost cells
    CALL exchange(1,0,fd(b1L,b2L,b3L,3))
    CALL exchange(2,0,fd(b1L,b2L,b3L,3))
    CALL exchange(3,0,fd(b1L,b2L,b3L,3))
    ! Force comp 3 onto vel 3 grid
    CALL interpolate2_pre_vel(.FALSE., 3, fd(b1L,b2L,b3L,3), inter)
    fd(:,:,:,3) = inter(:,:,:)
    ! exchange across boundaries on vel grid
    CALL exchange(1,3,fd(b1L,b2L,b3L,3))
    CALL exchange(2,3,fd(b1L,b2L,b3L,3))
    CALL exchange(3,3,fd(b1L,b2L,b3L,3))
  END IF

  END SUBROUTINE interpolate_force_pre_vel



  SUBROUTINE residual2volume_force(dir)

  IMPLICIT NONE

  INTEGER, INTENT(IN)    :: dir
  INTEGER                :: i,j,k
  REAL(8)                   :: vol

  ! At the moment it is all done on the pressure grid, so before interpolating
  ! the force onto the velocity grids!

  IF (dimens .EQ. 3) THEN
    DO k = S3p, N3p
      DO j = S2p, N2p
        DO i = S1p, N1p
          vol = (dx1p(i  ) * dx2p(j  ) * dx3p(k  ))/ 4. &
              + (dx1p(i+1) * dx2p(j  ) * dx3p(k  ))/12. &
              + (dx1p(i  ) * dx2p(j+1) * dx3p(k  ))/12. &
              + (dx1p(i  ) * dx2p(j  ) * dx3p(k+1))/12. &
              + (dx1p(i+1) * dx2p(j+1) * dx3p(k  ))/12. &
              + (dx1p(i+1) * dx2p(j  ) * dx3p(k+1))/12. &
              + (dx1p(i  ) * dx2p(j+1) * dx3p(k+1))/12. &
              + (dx1p(i+1) * dx2p(j+1) * dx3p(k+1))/4.
          IF (dir .EQ. 1) fd(i,j,k,1) = fd(i,j,k,1) / vol
          IF (dir .EQ. 2) fd(i,j,k,2) = fd(i,j,k,2) / vol
          IF (dir .EQ. 3) fd(i,j,k,3) = fd(i,j,k,3) / vol
        END DO
      END DO
    END DO
  ELSE
    k = 1
      DO j = S2p, N2p
        DO i = S1p, N1p
          vol = (dx1p(i  ) * dx2p(j  ))/3. &
              + (dx1p(i+1) * dx2p(j  ))/6. &
              + (dx1p(i  ) * dx2p(j+1))/6. &
              + (dx1p(i+1) * dx2p(j+1))/3.
          IF (dir .EQ. 1) fd(i,j,k,1) = fd(i,j,k,1) / vol
          IF (dir .EQ. 2) fd(i,j,k,2) = fd(i,j,k,2) / vol
          IF (dir .EQ. 3) WRITE(*,*) 'WARNING: There is no 3rd direction in a 2D force.'
        END DO
      END DO
    ! end k
  END IF


  END SUBROUTINE residual2volume_force





END MODULE usr_func
