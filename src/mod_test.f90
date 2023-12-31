!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!* modified by Barna Becsek, ARTORG CVE, University of Bern (barna.becsek@artorg.unibe.ch)	             *
!* October 2014
!*************************************************************************************************************

!> @file mod_test.f90
!! File containing mod_test with various subroutines for testing


!> Module containing public subroutines to test parameters and operators.
!! It uses mod_dims, mod_vars, mod_exchange, mod_diff, mod_helmholtz, mod_laplace.
MODULE mod_test
  
  USE mod_dims
  USE mod_vars
  USE mod_exchange
  USE mod_diff      ! CALL divergence (analyze_matrix), CALL gradient (analyze_matrix)
  USE mod_laplace   ! CALL product_div_grad (analyze_matrix)
  !USE mpi  
  
  PRIVATE
  
  PUBLIC test_parameter
  PUBLIC test_divergence, test_momentum
  PUBLIC analyze_matrix
  
  
  
  INCLUDE 'mpif.h'
  
  CONTAINS
  
!pgi$g unroll = n:8
!!pgi$r unroll = n:8
!!pgi$l unroll = n:8
  
  
  !> Subroutine for testing the compatibility of the global dimensions, the block division, FD stencil dimensions,
  !! periodic and advective boundaries and temporal parameters.
  !! @todo: remove parts related to particles and concentrations and compact differences
  !! @todo: write test routine for dtime_out_scal .LT. 0 & dtime_out_vect .LT. 0.
  !! @todo: write test for ndL <= n1U
  !! @todo: test dimXX > 1
  !! @todo: ndR =< b1U, ndR =< b2U, ndR =< b3U
  SUBROUTINE test_parameter
  
  ! revised: 26.05.08
  
  IMPLICIT NONE
  
  LOGICAL                ::  error(1:4)
  INTEGER                ::  m
  INTEGER                ::  comm_size
  
  
  !*******************************************************************************
  IF (M3 == 2) THEN ! TEST!!! wird auch in mod_setup::init_general nochmals ausgefuehrt!
     dimens = 2
  ELSE
     dimens = 3
  END IF
  
  IF (M3 == 2) THEN ! TEST!!! Notwendig? Pruefen, ggf. aendern und diese Fallunterscheidung loeschen ...
     BC_3L_global = -1
     BC_3U_global = -1
     
     impl_dir(3) = 0
     
     outlet   (3,:,:) = .FALSE.
  END IF
  !*******************************************************************************
  
  error(1:4) = .FALSE.
  
  
  ! TEST!!! - Testroutine schreiben f�r dtime_out_scal .LT. 0. & dtime_out_vect .LT. 0. ...
  !         - noch einen Test schreiben um ndL <= n1U zu testen ...
  !         - dimXX > 1 ebenfalls testen ...
  !         - ndR =< b1U, ndR =< b2U, ndR =< b3U
  IF (rank == 0) THEN
     
     !========================================================================================================
     !=== Dimensionen ========================================================================================
     !========================================================================================================
     ! TEST!!! unschoen ...
#ifdef ALLOC
     N1 = 1+(M1-1)/NB1
     N2 = 1+(M2-1)/NB2
     N3 = 1+(M3-1)/NB3
#endif
     !--------------------------------------------------------------------------------------------------------
     IF (dimens == 2 .AND. N3 /= 2) THEN
        WRITE(*,*) 'ERROR! Choose N3 = 2 for 2D!'
        error(1) = .TRUE.
     END IF
     !--------------------------------------------------------------------------------------------------------

     !--------------------------------------------------------------------------------------------------------
     IF (N1 .LT. 3) THEN
        WRITE(*,*) 'ERROR! N1 < 3!'
        error(1) = .TRUE.
     END IF
     IF (N2 .LT. 3) THEN
        WRITE(*,*) 'ERROR! N2 < 3!'
        error(1) = .TRUE.
     END IF
     IF (N3 .LT. 3 .AND. dimens == 3) THEN
        WRITE(*,*) 'ERROR! N3 < 3!'
        error(1) = .TRUE.
     END IF
     !--------------------------------------------------------------------------------------------------------
     IF (MOD(N1-1,2) /= 0) THEN
        WRITE(*,'(a)') 'ERROR! Dimension N1 cannot be used for multigrid!'
        error(1) = .TRUE.
     END IF
     IF (MOD(N2-1,2) /= 0) THEN
        WRITE(*,'(a)') 'ERROR! Dimension N2 cannot be used for multigrid!'
        error(1) = .TRUE.
     END IF
     IF (MOD(N3-1,2) /= 0 .AND. dimens == 3) THEN
        WRITE(*,'(a)') 'ERROR! Dimension N3 cannot be used for multigrid!'
        error(1) = .TRUE.
     END IF
     !========================================================================================================
     
     
     !========================================================================================================
     !=== Bockaufteilung =====================================================================================
     !========================================================================================================
     CALL MPI_COMM_SIZE(MPI_COMM_WORLD,comm_size,merror)
     IF (NB1*NB2*NB3 /= comm_size) THEN
        IF (rank == 0) WRITE(*,'(a,i4,a,i4,a)') 'ERROR! Number of blocks (=', NB1*NB2*NB3, ') differs from number of allocated processors (=', comm_size, ')!'
        error(1) = .TRUE.
     END IF
     !--------------------------------------------------------------------------------------------------------
     IF (dimens == 2 .AND. NB3 /= 1) THEN
        WRITE(*,*) 'ERROR! Choose NB3 = 1 for 2D!'
        error(1) = .TRUE.
     END IF
     !--------------------------------------------------------------------------------------------------------
     IF (NB1 .LT. 1) THEN
        WRITE(*,*) 'ERROR! NB1 < 1!'
        error(1) = .TRUE.
     END IF
     IF (NB2 .LT. 1) THEN
        WRITE(*,*) 'ERROR! NB2 < 1!'
        error(1) = .TRUE.
     END IF
     IF (NB3 .LT. 1) THEN
        WRITE(*,*) 'ERROR! NB3 < 1!'
        error(1) = .TRUE.
     END IF
     !--------------------------------------------------------------------------------------------------------
     IF (ls1 /= 0 .AND. ls1 /= -1) THEN
        WRITE(*,*) 'ERROR! ls1 /= 0,-1!'
        error(1) = .TRUE.
     END IF
     IF (ls2 /= 0 .AND. ls2 /= -1) THEN
        WRITE(*,*) 'ERROR! ls2 /= 0,-1!'
        error(1) = .TRUE.
     END IF
     IF (ls3 /= 0 .AND. ls3 /= -1) THEN
        WRITE(*,*) 'ERROR! ls3 /= 0,-1!'
        error(1) = .TRUE.
     END IF
     !========================================================================================================
     
     
     !========================================================================================================
     !=== FD-Stencil-Dimensionen =============================================================================
     !========================================================================================================
     IF (b1L .GT. -1) THEN
        WRITE(*,*) 'ERROR! b1L > -1!'
        error(1) = .TRUE.
     END IF
     IF (b2L .GT. -1) THEN
        WRITE(*,*) 'ERROR! b2L > -1!'
        error(1) = .TRUE.
     END IF
     IF (b3L .GT. -1) THEN
        WRITE(*,*) 'ERROR! b3L > -1!'
        error(1) = .TRUE.
     END IF
     
     IF (b1U .LT. 1) THEN
        WRITE(*,*) 'ERROR! b1U < 1!'
        error(1) = .TRUE.
     END IF
     IF (b2U .LT. 1) THEN
        WRITE(*,*) 'ERROR! b2U < 1!'
        error(1) = .TRUE.
     END IF
     IF (b3U .LT. 1) THEN
        WRITE(*,*) 'ERROR! b3U < 1!'
        error(1) = .TRUE.
     END IF
     !--------------------------------------------------------------------------------------------------------

     !========================================================================================================
     
     
     !========================================================================================================
     !=== Periodische Raender pr�fen =========================================================================
     !========================================================================================================
     ! 1L/1U-Rand:
     IF ((BC_1L_global == -1 .AND. BC_1U_global /= -1) .OR. (BC_1L_global /= -1 .AND. BC_1U_global == -1)) error(2) = .TRUE.
     IF ((BC_2L_global == -1 .AND. BC_2U_global /= -1) .OR. (BC_2L_global /= -1 .AND. BC_2U_global == -1)) error(2) = .TRUE.
     IF ((BC_3L_global == -1 .AND. BC_3U_global /= -1) .OR. (BC_3L_global /= -1 .AND. BC_3U_global == -1)) error(2) = .TRUE.
     
     IF (error(2)) WRITE(*,*) 'ERROR! A periodic boundary cannot be combined with a boundary condition!'
     !========================================================================================================
     
     
     !========================================================================================================
     !=== Ausflussrand pr�fen ================================================================================
     !========================================================================================================
     ! 1L/1U-Rand:
     IF (outlet(1,1,1) .AND. BC_1L_global /= 1 .AND. BC_1L_global /= 3) error(3) = .TRUE.
     IF (outlet(1,2,1) .AND. BC_1U_global /= 1 .AND. BC_1U_global /= 3) error(3) = .TRUE.
     IF (outlet(1,1,2) .AND. BC_1L_global /= 1 .AND. BC_1L_global /= 2) error(3) = .TRUE.
     IF (outlet(1,2,2) .AND. BC_1U_global /= 1 .AND. BC_1U_global /= 2) error(3) = .TRUE.
     IF (outlet(1,1,3) .AND. BC_1L_global /= 1 .AND. BC_1L_global /= 3) error(3) = .TRUE.
     IF (outlet(1,2,3) .AND. BC_1U_global /= 1 .AND. BC_1U_global /= 3) error(3) = .TRUE.
     
     ! 2L/2U-Rand:
     IF (outlet(2,1,1) .AND. BC_2L_global /= 1 .AND. BC_2L_global /= 2) error(3) = .TRUE.
     IF (outlet(2,2,1) .AND. BC_2U_global /= 1 .AND. BC_2U_global /= 2) error(3) = .TRUE.
     IF (outlet(2,1,2) .AND. BC_2L_global /= 1 .AND. BC_2L_global /= 3) error(3) = .TRUE.
     IF (outlet(2,2,2) .AND. BC_2U_global /= 1 .AND. BC_2U_global /= 3) error(3) = .TRUE.
     IF (outlet(2,1,3) .AND. BC_2L_global /= 1 .AND. BC_2L_global /= 3) error(3) = .TRUE.
     IF (outlet(2,2,3) .AND. BC_2U_global /= 1 .AND. BC_2U_global /= 3) error(3) = .TRUE.
     
     ! 3L/3U-Rand:
     IF (outlet(3,1,1) .AND. BC_3L_global /= 1 .AND. BC_3L_global /= 3) error(3) = .TRUE.
     IF (outlet(3,2,1) .AND. BC_3U_global /= 1 .AND. BC_3U_global /= 3) error(3) = .TRUE.
     IF (outlet(3,1,2) .AND. BC_3L_global /= 1 .AND. BC_3L_global /= 3) error(3) = .TRUE.
     IF (outlet(3,2,2) .AND. BC_3U_global /= 1 .AND. BC_3U_global /= 3) error(3) = .TRUE.
     IF (outlet(3,1,3) .AND. BC_3L_global /= 1 .AND. BC_3L_global /= 2) error(3) = .TRUE.
     IF (outlet(3,2,3) .AND. BC_3U_global /= 1 .AND. BC_3U_global /= 2) error(3) = .TRUE.
     
     IF (error(3)) WRITE(*,*) 'ERROR! Choice of outlet configuartion is not suitable!'
     !========================================================================================================
     
     
     !========================================================================================================

     !========================================================================================================
     
     
     !========================================================================================================
     !=== Sonstiges pr�fen ===================================================================================
     !========================================================================================================
     IF (time_start .GT. time_end) THEN 
        WRITE(*,*) 'ERROR! time_start > time_end!'
        error(1) = .TRUE.
     END IF
     
     !========================================================================================================
     
  END IF
  
  
  CALL MPI_BCAST(error(1:4),4,MPI_LOGICAL,0,MPI_COMM_WORLD,merror) ! MPI_COMM_CART ist hier noch nicht initialisiert ...
  
  IF (error(1) .OR. error(2) .OR. error(3) .OR. error(4)) THEN
     IF (rank == 0) WRITE(*,*) 'Exiting ...'
     CALL MPI_FINALIZE(merror)
     STOP
  END IF
  
  
  !===========================================================================================================
  !=== Zeitintegration =======================================================================================
  !===========================================================================================================
  ! Euler_yes ==> thetaL ==> timeint_mode ==> ...
  
  IF (thetaL == 0. .AND. timeint_mode == 0) THEN
     WRITE(*,*) 'WARNING! Setting timeint_mode = 1 for thetaL == 0. ...'
     timeint_mode = 1
  END IF
  
  IF (Euler_yes .AND. timeint_mode == 0) THEN
     WRITE(*,*) 'WARNING! Setting timeint_mode = 1 for Euler_yes == .TRUE. ...'
     timeint_mode = 1
  END IF
  !===========================================================================================================
  
  
  END SUBROUTINE test_parameter
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE test_divergence
  
  IMPLICIT NONE
  
  INTEGER                ::  i, j, k, m
  REAL(8)                   ::  max_div(1:2), max_div_global(1:2)
  
  
  ! TEST!!! handle_corner_rhs(res) evtl. zu Divergenz hineinziehen?
  CALL divergence2(vel,res)
  
  !IF (corner_yes) CALL handle_corner_rhs(res) ! TEST!!! mod_solvers wird z.Zt. nicht gelinked ...
  ! weight ist in den Ecken allerdings auch zu Null gesetzt ...
  
  max_div = 0.
  
  DO k = S3p, N3p
     DO j = S2p, N2p
        DO i = S1p, N1p
           IF (ABS(res(i,j,k)) .GE. max_div(1)) THEN
              max_div(1) = ABS(res(i,j,k))
           END IF
           IF (ABS(res(i,j,k)*weight(i,j,k)) .GE. max_div(2)) THEN
              max_div(2) = ABS(res(i,j,k)*weight(i,j,k))
           END IF
        END DO
     END DO
  END DO
  
  CALL MPI_REDUCE(max_div,max_div_global,2,MPI_REAL8,MPI_MAX,0,COMM_CART,merror)
  
  IF (rank == 0) THEN
     WRITE(*,'(a)')
     WRITE(*,'(a,E25.17)') 'MAX(div(u))        =', max_div_global(1)
     WRITE(*,'(a,E25.17)') 'MAX(div(u)*weight) =', max_div_global(2)
  END IF
  
  ! F�r Auswertung der Iterationsstatistiken:
  IF (timestep == 0) max_div_init = max_div_global
  
  
  END SUBROUTINE test_divergence
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE test_momentum
  
  IMPLICIT NONE
  
  INTEGER                ::  i, j, k, m
  REAL(8)                   ::  max_err, max_err_global
  REAL(8)                   ::  dd1, dd2
  
  
  dd1 = aRK(RK_steps)+bRK(RK_steps)
  
  !===========================================================================================================
  m         = 1
  direction = 1
  
  CALL gradient(m,pre,gpre)
  
  max_err = 0.
  
  DO k = S31B, N31B
     DO j = S21B, N21B
        DO i = S11B, N11B
           dd2 = ABS(rhs(i,j,k,m)-res(i,j,k)-dd1*gpre(i,j,k))
           IF (dd2 .GT. max_err) max_err = dd2
        END DO
     END DO
  END DO
  
  CALL MPI_REDUCE(max_err,max_err_global,1,MPI_REAL8,MPI_MAX,0,COMM_CART,merror)
  
  IF (rank == 0) WRITE(*,'(a      )')
  IF (rank == 0) WRITE(*,'(a,E13.5)') 'MAX(res(vel_1)) =', max_err_global
  !-----------------------------------------------------------------------------------------------------------
  m         = 2
  direction = 2
  
  CALL gradient(m,pre,gpre)
  
  max_err = 0.
  
  DO k = S32B, N32B
     DO j = S22B, N22B
        DO i = S12B, N12B
           dd2 = ABS(rhs(i,j,k,m)-res(i,j,k)-dd1*gpre(i,j,k))
           IF (dd2 .GT. max_err) max_err = dd2
        END DO
     END DO
  END DO
  
  CALL MPI_REDUCE(max_err,max_err_global,1,MPI_REAL8,MPI_MAX,0,COMM_CART,merror)
  
  IF (rank == 0) WRITE(*,'(a      )')
  IF (rank == 0) WRITE(*,'(a,E13.5)') 'MAX(res(vel_2)) =', max_err_global
  !-----------------------------------------------------------------------------------------------------------
  m         = 3
  direction = 3
  
  CALL gradient(m,pre,gpre)
  
  max_err = 0.
  
  DO k = S33B, N33B
     DO j = S23B, N23B
        DO i = S13B, N13B
           dd2 = ABS(rhs(i,j,k,m)-res(i,j,k)-dd1*gpre(i,j,k))
           IF (dd2 .GT. max_err) max_err = dd2
        END DO
     END DO
  END DO
  
  CALL MPI_REDUCE(max_err,max_err_global,1,MPI_REAL8,MPI_MAX,0,COMM_CART,merror)
  
  IF (rank == 0) WRITE(*,'(a      )')
  IF (rank == 0) WRITE(*,'(a,E13.5)') 'MAX(res(vel_3)) =', max_err_global
  !===========================================================================================================
  
  
  END SUBROUTINE test_momentum
  
  

  
 
  
  
END MODULE mod_test
