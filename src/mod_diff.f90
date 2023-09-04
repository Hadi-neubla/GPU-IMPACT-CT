!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!* modified by Barna Becsek, ARTORG CVE, University of Bern (barna.becsek@artorg.unibe.ch)                   *
!* October 2014                                                                                              *
!* GPU version by Hadi Zolfaghari , ARTORG CVE, then DAMTP, Cambridge University (hz382@cam.ac.uk)           *
!* Oct 2015 - Sep 2023                                                                                       *
!*************************************************************************************************************


!> module containing subroutines for constructing difference operators and for interpolating values. It uses
!! modules mod_dims, mod_vars and mod_exchange.
MODULE mod_diff
  
  
  USE mod_dims
  USE mod_vars
  USE mod_vars_GPU
  USE mod_exchange
  !USE mpi
  
  PRIVATE
  
  
  PUBLIC divergence, divergence2, divergence_transp
  PUBLIC gradient, gradient_transp
  PUBLIC Helmholtz, Helmholtz_explicit
  PUBLIC nonlinear, nonlinear_conc
  PUBLIC interpolate_vel
  PUBLIC outflow_bc
  PUBLIC bc_extrapolation, bc_extrapolation_transp
  PUBLIC interpolate_pre_vel , interpolate_vel_pre
  PUBLIC interpolate2_pre_vel, interpolate2_vel_pre
  PUBLIC first_pre_vel, first_vel_pre
  PUBLIC first_adv_pre, first_adv_vel
  PUBLIC Helmholtz_pre_explicit ! TEST!!!
  
  
  INCLUDE 'mpif.h'
  
  CONTAINS
  
!pgi$g unroll = n:8
!!pgi$r unroll = n:8
!!pgi$l unroll = n:8
  
  
  ! TEST!!! Generell bei Randbehandlung bei ii=0,1 beginnen, bzw. bis ii=N1 rechnen!  
 
  !> subroutine that computes the divergece operator (on the pressure grid) in one spacial direction.
  SUBROUTINE divergence(m,phi,div)
  
  IMPLICIT NONE
  
  INTEGER, INTENT(in)    ::  m                                           !< spacial dimension
  
  REAL(8)   , INTENT(inout) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)) !< ?
  REAL(8)   , INTENT(inout) ::  div(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)) !< divergence operator
  
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  
  
  CALL exchange(m,m,phi)

  !===========================================================================================================
  IF (m == 1) THEN

     !--------------------------------------------------------------------------------------------------------
        DO k = S3p, N3p
           DO j = S2p, N2p
              DO i = S1p, N1p
                 div(i,j,k) = cDu1(d1L,i)*phi(i+d1L,j,k)
!pgi$ unroll = n:8
                 DO ii = d1L+1, d1U
                    div(i,j,k) = div(i,j,k) + cDu1(ii,i)*phi(i+ii,j,k)
                 END DO
              END DO
           END DO
        END DO
     !--------------------------------------------------------------------------------------------------------
  END IF
  !===========================================================================================================
  IF (m == 2) THEN
     !--------------------------------------------------------------------------------------------------------
        DO k = S3p, N3p
           DO j = S2p, N2p
              DO i = S1p, N1p
!pgi$ unroll = n:8
                 DO jj = d2L, d2U
                    div(i,j,k) = div(i,j,k) + cDv2(jj,j)*phi(i,j+jj,k)
                 END DO
              END DO
           END DO
        END DO
     !--------------------------------------------------------------------------------------------------------
  END IF
  !===========================================================================================================
  IF (m == 3) THEN
     !--------------------------------------------------------------------------------------------------------
        DO k = S3p, N3p
           DO j = S2p, N2p
              DO i = S1p, N1p
!pgi$ unroll = n:8
                 DO kk = d3L, d3U
                    div(i,j,k) = div(i,j,k) + cDw3(kk,k)*phi(i,j,k+kk)
                 END DO
              END DO
           END DO
        END DO
     !--------------------------------------------------------------------------------------------------------
  END IF
  !=========================================================================================================
  
  
  END SUBROUTINE divergence
  
  
  
  
  
  
  
  
  
  
  !> subroutine that computes the divergence (on the pressure grid) in all spacial dimensions 
  SUBROUTINE divergence2(phi,div)
  
  IMPLICIT NONE
  
  REAL(8)   , INTENT(inout) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3)
  REAL(8)   , INTENT(  out) ::  div(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
!*  REAL(8)    ::  phi_copy(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3)
!*  REAL(8)    ::  div_copy(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  
  
  CALL exchange(1,1,phi(b1L,b2L,b3L,1))
  CALL exchange(2,2,phi(b1L,b2L,b3L,2))
  CALL exchange(3,3,phi(b1L,b2L,b3L,3))
 

!  print*, b1L, b1U, d1L, d1U, g1L, g1U, 'b1L, b1U, d1L, d1U, g1L, g1U' 
  
  !===========================================================================================================
     !--------------------------------------------------------------------------------------------------------
     IF (dimens == 3) THEN
!*     phi_copy =phi
!*     div_copy =div


     if (GPU_accelerated == 1) then

        call device_divergence2(phi, div)
     else

!*      call device_divergence2(phi_copy, div_copy) 
        DO k = S3p, N3p
           DO j = S2p, N2p
              DO i = S1p, N1p
                 div(i,j,k) = cDu1(d1L,i)*phi(i+d1L,j,k,1)
!pgi$ unroll = n:8
                 DO ii = d1L+1, d1U
                    div(i,j,k) = div(i,j,k) + cDu1(ii,i)*phi(i+ii,j,k,1)
                 END DO
!pgi$ unroll = n:8
                 DO jj = d2L, d2U
                    div(i,j,k) = div(i,j,k) + cDv2(jj,j)*phi(i,j+jj,k,2)
                 END DO
!pgi$ unroll = n:8
                 DO kk = d3L, d3U
                    div(i,j,k) = div(i,j,k) + cDw3(kk,k)*phi(i,j,k+kk,3)
                 END DO


!*        if (abs(div(i,j,k)-div_copy(i,j,k)) .ge. 0.000000000001) print*, 'mismatch', div_copy(i,j,k), div(i,j,k)
              END DO
           END DO
        END DO


      endif
        
     !--------------------------------------------------------------------------------------------------------
     ELSE
        
        DO k = S3p, N3p
           DO j = S2p, N2p
              DO i = S1p, N1p
                 div(i,j,k) = cDu1(d1L,i)*phi(i+d1L,j,k,1)
!pgi$ unroll = n:8
                 DO ii = d1L+1, d1U
                    div(i,j,k) = div(i,j,k) + cDu1(ii,i)*phi(i+ii,j,k,1)
                 END DO
!pgi$ unroll = n:8
                 DO jj = d2L, d2U
                    div(i,j,k) = div(i,j,k) + cDv2(jj,j)*phi(i,j+jj,k,2)
                 END DO
              END DO
           END DO
        END DO
        
     END IF
     !--------------------------------------------------------------------------------------------------------
  !===========================================================================================================
  
  
  END SUBROUTINE divergence2
  
  
  
  
  
  
  
  
  
  
  !> compute the transpose of the divergence operator 
  SUBROUTINE divergence_transp(m,phi,div)
  
  IMPLICIT NONE
  
  INTEGER, INTENT(in   ) ::  m                                          !< spacial direction
  
  REAL(8)   , INTENT(inout) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))!< ?
  REAL(8)   , INTENT(inout) ::  div(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))!< divergence (transpose) operator
  
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  
  !===========================================================================================================
  IF (m == 1) THEN
     !--------------------------------------------------------------------------------------------------------
        CALL exchange(m,0,phi)
        
        DO k = S31B, N31B ! "B"-Grenzen etwas konsistenter, aber nicht wirklich notwendig, da ohnehin mit bc_extrapolation nachmultipliziert wird ...
           DO j = S21B, N21B
              DO i = S11B, N11B
                 div(i,j,k) = cDu1T(g1L,i)*phi(i+g1L,j,k)
!pgi$ unroll = n:8
                 DO ii = g1L+1, g1U
                    div(i,j,k) = div(i,j,k) + cDu1T(ii,i)*phi(i+ii,j,k)
                 END DO
              END DO
           END DO
        END DO
     !--------------------------------------------------------------------------------------------------------
  END IF
  !===========================================================================================================
  IF (m == 2) THEN
     !--------------------------------------------------------------------------------------------------------
        CALL exchange(m,0,phi)
        
        DO k = S32B, N32B
           DO j = S22B, N22B
              DO i = S12B, N12B
                 div(i,j,k) = cDv2T(g2L,j)*phi(i,j+g2L,k)
!pgi$ unroll = n:8
                 DO jj = g2L+1, g2U
                    div(i,j,k) = div(i,j,k) + cDv2T(jj,j)*phi(i,j+jj,k)
                 END DO
              END DO
           END DO
        END DO
     !--------------------------------------------------------------------------------------------------------
  END IF
  !===========================================================================================================
  IF (m == 3) THEN
     !--------------------------------------------------------------------------------------------------------
        CALL exchange(m,0,phi)
        
        DO k = S33B, N33B
           DO j = S23B, N23B
              DO i = S13B, N13B
                 div(i,j,k) = cDw3T(g3L,k)*phi(i,j,k+g3L)
!pgi$ unroll = n:8
                 DO kk = g3L+1, g3U
                    div(i,j,k) = div(i,j,k) + cDw3T(kk,k)*phi(i,j,k+kk)
                 END DO
              END DO
           END DO
        END DO
     !--------------------------------------------------------------------------------------------------------
  END IF
  !===========================================================================================================
  
  
  END SUBROUTINE divergence_transp
  
  
  
  
  
  
  
  
  
  
  !> subroutine that computes the gradient in direction m
  SUBROUTINE gradient(m,phi,grad)
  
  IMPLICIT NONE
  
  INTEGER, INTENT(in   ) ::  m                                           !< spacial direction
  real(8)                   ::  start, finish
  REAL(8)   , INTENT(inout) ::  phi (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))!< ? 
  REAL(8)   , INTENT(  out) ::  grad(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))!< gradient operator
  
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Randbedingungen könnten nur zum Teil in die Stencils eingebaut werden, so dass sich das   !
  !                vermutlich nicht wirklich lohnt.                                                          !
  !----------------------------------------------------------------------------------------------------------!
  
  
  CALL exchange(m,0,phi)
 
 
 
  !===========================================================================================================
  IF (m == 1) THEN
     !--------------------------------------------------------------------------------------------------------
  !print*, S11, N11, S21, N21, s31, N31, 'S11,N11, S21, N21, s31, N31'
!!!$omp parallel do
        DO k = S31, N31
           DO j = S21, N21
              DO i = S11, N11
                 grad(i,j,k) = cGp1(g1L,i)*phi(i+g1L,j,k)
!pgi$ unroll = n:8
                 DO ii = g1L+1, g1U
                    grad(i,j,k) = grad(i,j,k) + phi(i+ii,j,k) * cGp1(ii,i) 
                 END DO
              END DO
           END DO
        END DO
     !--------------------------------------------------------------------------------------------------------
     
     !--- Randbedingungen ------------------------------------------------------------------------------------
     IF (BC_1L .GT. 0) grad(0 ,S21B:N21B,S31B:N31B) = 0.
     IF (BC_1U .GT. 0) grad(N1,S21B:N21B,S31B:N31B) = 0.
     IF (BC_2L .GT. 0) grad(S11B:N11B,1 ,S31B:N31B) = 0.
     IF (BC_2U .GT. 0) grad(S11B:N11B,N2,S31B:N31B) = 0.
     IF (BC_3L .GT. 0) grad(S11B:N11B,S21B:N21B,1 ) = 0.
     IF (BC_3U .GT. 0) grad(S11B:N11B,S21B:N21B,N3) = 0.
     
  END IF
  !===========================================================================================================
  IF (m == 2) THEN
     !--------------------------------------------------------------------------------------------------------
      !!!$omp parallel do 
      DO k = S32, N32
           DO j = S22, N22
              DO i = S12, N12
                 grad(i,j,k) = cGp2(g2L,j)*phi(i,j+g2L,k)
!pgi$ unroll = n:8
                 DO jj = g2L+1, g2U
                    grad(i,j,k) = grad(i,j,k) + cGp2(jj,j)*phi(i,j+jj,k)
                 END DO
              END DO
           END DO
        END DO
     !--------------------------------------------------------------------------------------------------------
     
     !--- Randbedingungen ------------------------------------------------------------------------------------
     IF (BC_1L .GT. 0) grad(1 ,S22B:N22B,S32B:N32B) = 0.
     IF (BC_1U .GT. 0) grad(N1,S22B:N22B,S32B:N32B) = 0.
     IF (BC_2L .GT. 0) grad(S12B:N12B,0 ,S32B:N32B) = 0.
     IF (BC_2U .GT. 0) grad(S12B:N12B,N2,S32B:N32B) = 0.
     IF (BC_3L .GT. 0) grad(S12B:N12B,S22B:N22B,1 ) = 0.
     IF (BC_3U .GT. 0) grad(S12B:N12B,S22B:N22B,N3) = 0.
     
  END IF
  !===========================================================================================================
  IF (m == 3) THEN
     !--------------------------------------------------------------------------------------------------------
       DO k = S33, N33
           DO j = S23, N23
              DO i = S13, N13
                 grad(i,j,k) = cGp3(g3L,k)*phi(i,j,k+g3L)
!pgi$ unroll = n:8
                 DO kk = g3L+1, g3U
                    grad(i,j,k) = grad(i,j,k) + cGp3(kk,k)*phi(i,j,k+kk)
                 END DO
              END DO
           END DO
        END DO
    !--------------------------------------------------------------------------------------------------------
     
     !--- Randbedingungen ------------------------------------------------------------------------------------
     IF (BC_1L .GT. 0) grad(1 ,S23B:N23B,S33B:N33B) = 0.
     IF (BC_1U .GT. 0) grad(N1,S23B:N23B,S33B:N33B) = 0.
     IF (BC_2L .GT. 0) grad(S13B:N13B,1 ,S33B:N33B) = 0.
     IF (BC_2U .GT. 0) grad(S13B:N13B,N2,S33B:N33B) = 0.
     IF (BC_3L .GT. 0) grad(S13B:N13B,S23B:N23B,0 ) = 0.
     IF (BC_3U .GT. 0) grad(S13B:N13B,S23B:N23B,N3) = 0.
     
  END IF
  !===========================================================================================================
  
  
  END SUBROUTINE gradient
  
  
  
  
  
  
  
  
  
  
  !> computes the transpose of the gradient operator 
  SUBROUTINE gradient_transp(m,phi,grad)
  
  IMPLICIT NONE
  
  INTEGER, INTENT(in   ) ::  m                                           !< spacial direction
  
  REAL(8)   , INTENT(inout) ::  phi (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))!< ?
  REAL(8)   , INTENT(  out) ::  grad(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))!< gradient (transpose) operator
  
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Umgekehrte Reihenfolge im Vergleich zu Subroutine "gradient".                             !
  !              - Randbedingungen könnten nur zum Teil in die Stencils eingebaut werden, so dass sich das   !
  !                vermutlich nicht wirklich lohnt.                                                          !
  !----------------------------------------------------------------------------------------------------------!
  
  
  !===========================================================================================================
  IF (m == 1) THEN
     
     !--- Randbedingungen ------------------------------------------------------------------------------------
     IF (BC_1L .GT. 0) phi(0 ,S21B:N21B,S31B:N31B) = 0.
     IF (BC_1U .GT. 0) phi(N1,S21B:N21B,S31B:N31B) = 0.
     IF (BC_2L .GT. 0) phi(S11B:N11B,1 ,S31B:N31B) = 0.
     IF (BC_2U .GT. 0) phi(S11B:N11B,N2,S31B:N31B) = 0.
     IF (BC_3L .GT. 0) phi(S11B:N11B,S21B:N21B,1 ) = 0.
     IF (BC_3U .GT. 0) phi(S11B:N11B,S21B:N21B,N3) = 0.
     
     !--------------------------------------------------------------------------------------------------------
        CALL exchange(m,m,phi) ! Muss nach Randbedingungen kommen!
        
        DO k = S3p, N3p
           DO j = S2p, N2p
              DO i = S1p, N1p
                 grad(i,j,k) = cGp1T(d1L,i)*phi(i+d1L,j,k)
!pgi$ unroll = n:8
                 DO ii = d1L+1, d1U
                    grad(i,j,k) = grad(i,j,k) + cGp1T(ii,i)*phi(i+ii,j,k)
                 END DO
              END DO
           END DO
        END DO
     !--------------------------------------------------------------------------------------------------------
     
  END IF
  !===========================================================================================================
  IF (m == 2) THEN
     
     !--- Randbedingungen ------------------------------------------------------------------------------------
     IF (BC_1L .GT. 0) phi(1 ,S22B:N22B,S32B:N32B) = 0.
     IF (BC_1U .GT. 0) phi(N1,S22B:N22B,S32B:N32B) = 0.
     IF (BC_2L .GT. 0) phi(S12B:N12B,0 ,S32B:N32B) = 0.
     IF (BC_2U .GT. 0) phi(S12B:N12B,N2,S32B:N32B) = 0.
     IF (BC_3L .GT. 0) phi(S12B:N12B,S22B:N22B,1 ) = 0.
     IF (BC_3U .GT. 0) phi(S12B:N12B,S22B:N22B,N3) = 0.
     
     !--------------------------------------------------------------------------------------------------------
        CALL exchange(m,m,phi) ! Muss nach Randbedingungen kommen!
        
        DO k = S3p, N3p
           DO j = S2p, N2p
              DO i = S1p, N1p
!pgi$ unroll = n:8
                 DO jj = d2L, d2U
                    grad(i,j,k) = grad(i,j,k) + cGp2T(jj,j)*phi(i,j+jj,k)
                 END DO
              END DO
           END DO
        END DO
     !--------------------------------------------------------------------------------------------------------
     
  END IF
  !===========================================================================================================
  IF (m == 3) THEN
     
     !--- Randbedingungen ------------------------------------------------------------------------------------
     IF (BC_1L .GT. 0) phi(1 ,S23B:N23B,S33B:N33B) = 0.
     IF (BC_1U .GT. 0) phi(N1,S23B:N23B,S33B:N33B) = 0.
     IF (BC_2L .GT. 0) phi(S13B:N13B,1 ,S33B:N33B) = 0.
     IF (BC_2U .GT. 0) phi(S13B:N13B,N2,S33B:N33B) = 0.
     IF (BC_3L .GT. 0) phi(S13B:N13B,S23B:N23B,0 ) = 0.
     IF (BC_3U .GT. 0) phi(S13B:N13B,S23B:N23B,N3) = 0.
     
     !--------------------------------------------------------------------------------------------------------
        CALL exchange(m,m,phi) ! Muss nach Randbedingungen kommen!
        
        DO k = S3p, N3p
           DO j = S2p, N2p
              DO i = S1p, N1p
!pgi$ unroll = n:8
                 DO kk = d3L, d3U
                    grad(i,j,k) = grad(i,j,k) + cGp3T(kk,k)*phi(i,j,k+kk)
                 END DO
              END DO
           END DO
        END DO
     !--------------------------------------------------------------------------------------------------------
     
  END IF
  !===========================================================================================================
  
  
  END SUBROUTINE gradient_transp
  
  
  
  
  
  
  
  
  
 
  
  
  
  
 
  
  SUBROUTINE Helmholtz_explicit(exch_yes)
  
  IMPLICIT NONE
  
  LOGICAL, INTENT(in   ) ::  exch_yes
  
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  
  REAL(8)                   ::  dd1

! for unit testing hz
  REAL(8)                   ::  nl_copy(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),3) !> Laplacian
  REAL(8)                   ::  vel_copy(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),3) !> Laplacian

  integer                :: counter1, counter2, counter3

  real(8)                   ::  t1, t2, t3

  include 'mpif.h' 
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Für den Aufbau der RHS bei expliziter Zeitintegration.                                    !
  !              - Randbedingungen müssen daher nicht berücksichtigt werden.                                 !
  !----------------------------------------------------------------------------------------------------------!

 


  
  ! for GPU unit testing uncomment !* 
  
  
  !===========================================================================================================
  !-----------------------------------------------------------------------------------------------------------
     IF (dimens == 3) THEN
         if (GPU_accelerated == 1) then
    


 !    call MPI_Barrier(MPI_COMM_WORLD, merror)
 !    t1 =   MPI_Wtime()


#ifdef CUDAVERSION10
         call device_helmholtz(1,vel(:,:,:,1), nl(:,:,:,1))
#endif 

!     call MPI_Barrier(MPI_COMM_WORLD, merror)
!     t2 =   MPI_Wtime()
   
        else
!*         nl_copy(:,:,:,1) = nl(:,:,:,1)
!*         vel_copy(:,:,:,1) = vel(:,:,:,1)
!*         call device_helmholtz(1,vel_copy(:,:,:,1), nl_copy(:,:,:,1))
!*         counter1=0
            DO k = S31, N31
               DO j = S21, N21
                  DO i = S11, N11
                     dd1 = cu11(b1L,i)*vel(i+b1L,j,k,1)
    !pgi$ unroll = n:8
                     DO ii = b1L+1, b1U
                        dd1 = dd1 + cu11(ii,i)*vel(i+ii,j,k,1)
                     END DO
    !pgi$ unroll = n:8
                     DO jj = b2L, b2U
                        dd1 = dd1 + cp22(jj,j)*vel(i,j+jj,k,1)
                     END DO
    !pgi$ unroll = n:8
                     DO kk = b3L, b3U
                        dd1 = dd1 + cp33(kk,k)*vel(i,j,k+kk,1)
                     END DO
                     nl(i,j,k,1) = nl(i,j,k,1) - multL*dd1

  !*          if (abs(nl_copy(i,j,k,1) - nl(i,j,k,1)) .le. 0.000000000001) counter1=counter1+1 
                  END DO
               END DO
            END DO
  !*          print*, counter1, 'dir1'
         endif
!     call MPI_Barrier(MPI_COMM_WORLD, merror)
!     t3 =   MPI_Wtime()

 !    print*, 'GPU time = ', t2-t1, 'CPU time = ', t3-t2

     ELSE
        DO k = S31, N31
           DO j = S21, N21
              DO i = S11, N11
                 dd1 = cu11(b1L,i)*vel(i+b1L,j,k,1)
!pgi$ unroll = n:8
                 DO ii = b1L+1, b1U
                    dd1 = dd1 + cu11(ii,i)*vel(i+ii,j,k,1)
                 END DO
!pgi$ unroll = n:8
                 DO jj = b2L, b2U
                    dd1 = dd1 + cp22(jj,j)*vel(i,j+jj,k,1)
                 END DO
                 nl(i,j,k,1) = nl(i,j,k,1) - multL*dd1
              END DO
           END DO
        END DO
     END IF
  !-----------------------------------------------------------------------------------------------------------
  !===========================================================================================================
  !-----------------------------------------------------------------------------------------------------------
     IF (dimens == 3) THEN

        if (GPU_accelerated == 1) then
        !*  nl_copy = nl
        !*  vel_copy = vel

#ifdef CUDAVERSION10
           call device_helmholtz(2, vel(:,:,:,2), nl(:,:,:,2)) 
#endif
   
        else

    !*     nl_copy(:,:,:,2) = nl(:,:,:,2)
    !*     vel_copy(:,:,:,2) = vel(:,:,:,2)
    !*     call device_helmholtz(2,vel_copy(:,:,:,2), nl_copy(:,:,:,2))
    !*     counter2=0

           DO k = S32, N32
              DO j = S22, N22
                 DO i = S12, N12
                    dd1 = cp11(b1L,i)*vel(i+b1L,j,k,2)
   !pgi$ unroll = n:8
                    DO ii = b1L+1, b1U
                       dd1 = dd1 + cp11(ii,i)*vel(i+ii,j,k,2)
                    END DO
   !pgi$ unroll = n:8
                    DO jj = b2L, b2U
                       dd1 = dd1 + cv22(jj,j)*vel(i,j+jj,k,2)
                    END DO
   !pgi$ unroll = n:8
                    DO kk = b3L, b3U
                       dd1 = dd1 + cp33(kk,k)*vel(i,j,k+kk,2)
                    END DO
                    nl(i,j,k,2) = nl(i,j,k,2) - multL*dd1
     !*         if (abs(nl_copy(i,j,k,2) - nl(i,j,k,2)) .le. 0.000000000001) then
     !*           counter2=counter2+1
     !*         endif
                 END DO
              END DO
           END DO

     !*    print*, counter2, 'dir 2'
        endif

     ELSE
        DO k = S32, N32
           DO j = S22, N22
              DO i = S12, N12
                 dd1 = cp11(b1L,i)*vel(i+b1L,j,k,2)
                 DO ii = b1L+1, b1U
                    dd1 = dd1 + cp11(ii,i)*vel(i+ii,j,k,2)
                 END DO
!pgi$ unroll = n:8
                 DO jj = b2L, b2U
                    dd1 = dd1 + cv22(jj,j)*vel(i,j+jj,k,2)
                 END DO
                 nl(i,j,k,2) = nl(i,j,k,2) - multL*dd1
              END DO
           END DO
        END DO
     END IF
  !-----------------------------------------------------------------------------------------------------------
  !===========================================================================================================
  !-----------------------------------------------------------------------------------------------------------
  IF (dimens == 3) THEN

      if (GPU_accelerated == 1) then
 
#ifdef CUDAVERSION10
          call device_helmholtz(3,vel(:,:,:,3), nl(:,:,:,3))
#endif
      else
       !*   nl_copy = nl
       !*   vel_copy = vel
       !*   call device_helmholtz(3,vel_copy(:,:,:,3), nl_copy(:,:,:,3))
       !* counter3=0 
            DO k = S33, N33
               DO j = S23, N23
                  DO i = S13, N13
                     dd1 = cp11(b1L,i)*vel(i+b1L,j,k,3)
    !pgi$ unroll = n:8
                     DO ii = b1L+1, b1U
                        dd1 = dd1 + cp11(ii,i)*vel(i+ii,j,k,3)
                     END DO
    !pgi$ unroll = n:8
                     DO jj = b2L, b2U
                        dd1 = dd1 + cp22(jj,j)*vel(i,j+jj,k,3)
                     END DO
    !pgi$ unroll = n:8
                     DO kk = b3L, b3U
                        dd1 = dd1 + cw33(kk,k)*vel(i,j,k+kk,3)
                     END DO
                     nl(i,j,k,3) = nl(i,j,k,3) - multL*dd1

        !*   if (abs(nl_copy(i,j,k,3) - nl(i,j,k,3)) .le. 0.000000000001) counter3=counter3+1  
                  END DO
               END DO
            END DO
        !*  print*, counter3, 'dir3'
      endif

  END IF
  !-----------------------------------------------------------------------------------------------------------
  !===========================================================================================================
  
  
  END SUBROUTINE Helmholtz_explicit
  
  
  
  
  
  
  
  
  
  
 
  
  
  
  
  
  
  ! TEST!!! Teile davon (nur zentrale Operationen!) koennten jeweils durch interpolate2_pre_vel/interpolate2_vel_pre, first_adv_pre/first_adv_vel
  !         ersetzt werden (beachte aber Addition von nl!)
  ! TEST!!! umbenennen in advect... (?)
  !> computes advective (nonlinear) terms.
  SUBROUTINE nonlinear(exch_yes)
  
  IMPLICIT NONE
  
  LOGICAL, INTENT(in   ) ::  exch_yes
  
  REAL(8)                   ::  dd1
  
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk

! unit testing nonlinear kernel - hz 
  REAL(8)                   ::  nl_copy(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),3) !> Laplacian
  REAL(8)                   ::  vel_copy(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),3) !> Laplacian
  REAL(8)                   ::  pp_copy(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)) !> Laplacian

 
  real(8)                   ::  t1, t2
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - [Sij,Nij] ist immer eine Untermenge von [Sip,Nip]                                         !
  !              - Feld "res" wird mehrfach ausgetauscht, um mit einem statt drei skalaren Feldern arbeiten  !
  !                zu können. Im Prinzip könnte auch rhs(:,:,:,1:3) für die Zwischenspeicherung verwendet    !
  !                werden, jedoch sind dazu einige Umbaumassnahmen notwendig bei geringem Effizienzgewinn.   !
  !----------------------------------------------------------------------------------------------------------!
  
  
  ! worki muss bereits ausgetauscht sein!
  IF (exch_yes) CALL exchange_all_all(.TRUE.,vel)
  
  !===========================================================================================================
  !=== u*du/dx ===============================================================================================
  !===========================================================================================================
  CALL interpolate2_pre_vel(.FALSE.,1,work1,pp)
  !-----------------------------------------------------------------------------------------------------------
 ! print* , n1L, n1U, n2L, n2U, n3L, n3U, 'bs: ', b1L, b1U, b2L, b2U, b3L, b3U
 ! they are all -3,3
  
  IF (upwind_yes) THEN
  
!*  pp_copy = pp
!*  nl_copy = nl
!*  vel_copy = vel


  if ((GPU_accelerated == 1) .and. (dimens == 3)) then 
              call device_Nonlinear_upwind(1, vel(:,:,:,1), nl(:,:,:,1))
  else

!*    call device_Nonlinear_upwind(1, vel_copy(:,:,:,1), nl_copy(:,:,:,1))
    DO k = S31, N31
        DO j = S21, N21
           DO i = S11, N11
              IF (pp(i,j,k) .GE. 0.) THEN
                 dd1 = cNu1U(n1L,i)*vel(i+n1L,j,k,1)
!pgi$ unroll = n:8
                 DO ii = n1L+1, n1U
                    dd1 = dd1 + cNu1U(ii,i)*vel(i+ii,j,k,1)
                 END DO
              ELSE
                 dd1 = cNu1D(n1L,i)*vel(i+n1L,j,k,1)
!pgi$ unroll = n:8
                 DO ii = n1L+1, n1U
                    dd1 = dd1 + cNu1D(ii,i)*vel(i+ii,j,k,1)
                 END DO
              END IF
              
              nl(i,j,k,1) = dd1*pp(i,j,k)
              
!*              if ((abs(nl(i,j,k,1) - nl_copy(i,j,k,1)) .ge. 0.000000000001)) print*, 'mismatch', nl(i,j,k,1) , nl_copy(i,j,k,1)
              

           END DO
        END DO
     END DO

   endif

  ELSE

      if ((GPU_accelerated == 1) .and.  (dimens == 3))  then

        ! call cpu_time(t1)
     !    call device_Nonlinear(1, vel(:,:,:,1), nl(:,:,:,1))
       !  call cpu_time(t2)
        ! print*, t2-t1, 'this is Nonlinear On GPU'
      else
  !       call cpu_time(t1)
            DO k = S31, N31
               DO j = S21, N21
                  DO i = S11, N11
                     dd1 = cu1(b1L,i)*vel(i+b1L,j,k,1)
    !pgi$ unroll = n:8
                     DO ii = b1L+1, b1U
                        dd1 = dd1 + cu1(ii,i)*vel(i+ii,j,k,1)
                     END DO
                     
                     nl(i,j,k,1) = dd1*pp(i,j,k)
                  END DO
               END DO
            END DO
   !      call cpu_time(t2)
        ! print*, t2-t1, 'this is Nonlinear on CPU'
        ! print*, N1, N2, N3 
      endif

  END IF
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== v*du/dy ===============================================================================================
  !===========================================================================================================
  CALL interpolate2_pre_vel(.FALSE.,1,work2,pp)
  !-----------------------------------------------------------------------------------------------------------
  IF (upwind_yes) THEN

!*  pp_copy = pp
!*  nl_copy = nl
!*  vel_copy = vel

  if ((GPU_accelerated == 1) .and. (dimens == 3)) then 
              call device_Nonlinear_upwind(2, vel(:,:,:,1), nl(:,:,:,1))
  else

!*              call device_Nonlinear_upwind(2, vel_copy(:,:,:,1), nl_copy(:,:,:,1))

     DO k = S31, N31
        DO j = S21, N21
           DO i = S11, N11
              IF (pp(i,j,k) .GE. 0.) THEN
                 dd1 = cNp2U(n2L,j)*vel(i,j+n2L,k,1)
!pgi$ unroll = n:8
                 DO jj = n2L+1, n2U
                    dd1 = dd1 + cNp2U(jj,j)*vel(i,j+jj,k,1)
                 END DO
              ELSE
                 dd1 = cNp2D(n2L,j)*vel(i,j+n2L,k,1)
!pgi$ unroll = n:8
                 DO jj = n2L+1, n2U
                    dd1 = dd1 + cNp2D(jj,j)*vel(i,j+jj,k,1)
                 END DO
              END IF
              
              nl(i,j,k,1) = nl(i,j,k,1) + dd1*pp(i,j,k)

!*              if ((abs(nl(i,j,k,1) - nl_copy(i,j,k,1)) .ge. 0.000000000001)) print*, 'mismatch', nl(i,j,k,1) , nl_copy(i,j,k,1)
           END DO
        END DO
     END DO
  endif

  ELSE

      if ((GPU_accelerated == 1) .and. (dimens == 3)) then
    
      !   call device_Nonlinear(2, vel(:,:,:,1), nl(:,:,:,1))
    
      else
    
            DO k = S31, N31
               DO j = S21, N21
                  DO i = S11, N11
                     dd1 = cp2(b2L,j)*vel(i,j+b2L,k,1)
    !pgi$ unroll = n:8
                     DO jj = b2L+1, b2U
                        dd1 = dd1 + cp2(jj,j)*vel(i,j+jj,k,1)
                     END DO
                     
                     nl(i,j,k,1) = nl(i,j,k,1) + dd1*pp(i,j,k)
                  END DO
               END DO
            END DO
    
      endif

  END IF



  !===========================================================================================================
  
  
  IF (dimens == 3) THEN
  !===========================================================================================================
  !=== w*du/dz ===============================================================================================
  !===========================================================================================================
  CALL interpolate2_pre_vel(.FALSE.,1,work3,pp)
  !-----------------------------------------------------------------------------------------------------------
  IF (upwind_yes) THEN

  if (GPU_accelerated == 1) then 
              call device_Nonlinear_upwind(3, vel(:,:,:,1), nl(:,:,:,1))
  else


     DO k = S31, N31
        DO j = S21, N21
           DO i = S11, N11
              IF (pp(i,j,k) .GE. 0.) THEN
                 dd1 = cNp3U(n3L,k)*vel(i,j,k+n3L,1)
!pgi$ unroll = n:8
                 DO kk = n3L+1, n3U
                    dd1 = dd1 + cNp3U(kk,k)*vel(i,j,k+kk,1)
                 END DO
              ELSE
                 dd1 = cNp3D(n3L,k)*vel(i,j,k+n3L,1)
!pgi$ unroll = n:8
                 DO kk = n3L+1, n3U
                    dd1 = dd1 + cNp3D(kk,k)*vel(i,j,k+kk,1)
                 END DO
              END IF
              
              nl(i,j,k,1) = nl(i,j,k,1) + dd1*pp(i,j,k)
           END DO
        END DO
     END DO
  endif

  ELSE

      if (GPU_accelerated == 1) then
    
       !  call device_Nonlinear(3, vel(:,:,:,1), nl(:,:,:,1))
    
      else
    
            DO k = S31, N31
               DO j = S21, N21
                  DO i = S11, N11
                     dd1 = cp3(b3L,k)*vel(i,j,k+b3L,1)
    !pgi$ unroll = n:8
                     DO kk = b3L+1, b3U
                        dd1 = dd1 + cp3(kk,k)*vel(i,j,k+kk,1)
                     END DO
                     
                     nl(i,j,k,1) = nl(i,j,k,1) + dd1*pp(i,j,k)
                  END DO
               END DO
            END DO
       endif
  END IF
  !===========================================================================================================
  END IF
  
  
  
  
  
  
  !===========================================================================================================
  !=== u*dv/dx ===============================================================================================
  !===========================================================================================================
  CALL interpolate2_pre_vel(.FALSE.,2,work1,pp)
  !-----------------------------------------------------------------------------------------------------------
  IF (upwind_yes) THEN

  if ((GPU_accelerated == 1) .and. (dimens == 3)) then 
              call device_Nonlinear_upwind(4, vel(:,:,:,2), nl(:,:,:,2))
  else

     DO k = S32, N32
        DO j = S22, N22
           DO i = S12, N12
              IF (pp(i,j,k) .GE. 0.) THEN
                 dd1 = cNp1U(n1L,i)*vel(i+n1L,j,k,2)
!pgi$ unroll = n:8
                 DO ii = n1L+1, n1U
                    dd1 = dd1 + cNp1U(ii,i)*vel(i+ii,j,k,2)
                 END DO
              ELSE
                 dd1 = cNp1D(n1L,i)*vel(i+n1L,j,k,2)
!pgi$ unroll = n:8
                 DO ii = n1L+1, n1U
                    dd1 = dd1 + cNp1D(ii,i)*vel(i+ii,j,k,2)
                 END DO
              END IF
              
              nl(i,j,k,2) = dd1*pp(i,j,k)
           END DO
        END DO
     END DO
   endif
  ELSE

      if ((GPU_accelerated == 1) .and. (dimens == 3)) then
    
        ! call device_Nonlinear(4, vel(:,:,:,2), nl(:,:,:,2))
    
      else
    
            DO k = S32, N32
               DO j = S22, N22
                  DO i = S12, N12
                     dd1 = cp1(b1L,i)*vel(i+b1L,j,k,2)
    !pgi$ unroll = n:8
                     DO ii = b1L+1, b1U
                        dd1 = dd1 + cp1(ii,i)*vel(i+ii,j,k,2)
                     END DO
                     
                     nl(i,j,k,2) = dd1*pp(i,j,k)
                  END DO
               END DO
            END DO
    
    
      endif
  END IF

  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== v*dv/dy ===============================================================================================
  !===========================================================================================================
  CALL interpolate2_pre_vel(.FALSE.,2,work2,pp)
  !-----------------------------------------------------------------------------------------------------------
  IF (upwind_yes) THEN

  if ((GPU_accelerated == 1) .and. (dimens == 3)) then 
              call device_Nonlinear_upwind(5, vel(:,:,:,2), nl(:,:,:,2))
 else


     DO k = S32, N32
        DO j = S22, N22
           DO i = S12, N12
              IF (pp(i,j,k) .GE. 0.) THEN
                 dd1 = cNv2U(n2L,j)*vel(i,j+n2L,k,2)
!pgi$ unroll = n:8
                 DO jj = n2L+1, n2U
                    dd1 = dd1 + cNv2U(jj,j)*vel(i,j+jj,k,2)
                 END DO
              ELSE
                 dd1 = cNv2D(n2L,j)*vel(i,j+n2L,k,2)
!pgi$ unroll = n:8
                 DO jj = n2L+1, n2U
                    dd1 = dd1 + cNv2D(jj,j)*vel(i,j+jj,k,2)
                 END DO
              END IF
              
              nl(i,j,k,2) = nl(i,j,k,2) + dd1*pp(i,j,k)
           END DO
        END DO
     END DO
   endif
  ELSE

      if ((GPU_accelerated == 1) .and. (dimens == 3)) then
    
        
         ! call device_Nonlinear(5, vel(:,:,:,2), nl(:,:,:,2))
    
      else
    
            DO k = S32, N32
               DO j = S22, N22
                  DO i = S12, N12
                     dd1 = cv2(b2L,j)*vel(i,j+b2L,k,2)
    !pgi$ unroll = n:8
                     DO jj = b2L+1, b2U
                        dd1 = dd1 + cv2(jj,j)*vel(i,j+jj,k,2)
                     END DO
                     
                     nl(i,j,k,2) = nl(i,j,k,2) + dd1*pp(i,j,k)
                  END DO
               END DO
            END DO
     
       endif
 
   END IF
  !===========================================================================================================
  
  
  IF (dimens == 3) THEN
  !===========================================================================================================
  !=== w*dv/dz ===============================================================================================
  !===========================================================================================================
  CALL interpolate2_pre_vel(.FALSE.,2,work3,pp)
  !-----------------------------------------------------------------------------------------------------------
  IF (upwind_yes) THEN

    if (GPU_accelerated == 1) then 
              call device_Nonlinear_upwind(6, vel(:,:,:,2), nl(:,:,:,2))
   else



     DO k = S32, N32
        DO j = S22, N22
           DO i = S12, N12
              IF (pp(i,j,k) .GE. 0.) THEN
                 dd1 = cNp3U(n3L,k)*vel(i,j,k+n3L,2)
!pgi$ unroll = n:8
                 DO kk = n3L+1, n3U
                    dd1 = dd1 + cNp3U(kk,k)*vel(i,j,k+kk,2)
                 END DO
              ELSE
                 dd1 = cNp3D(n3L,k)*vel(i,j,k+n3L,2)
!pgi$ unroll = n:8
                 DO kk = n3L+1, n3U
                    dd1 = dd1 + cNp3D(kk,k)*vel(i,j,k+kk,2)
                 END DO
              END IF
              
              nl(i,j,k,2) = nl(i,j,k,2) + dd1*pp(i,j,k)
           END DO
        END DO
     END DO
   
    endif

  ELSE

      if (GPU_accelerated == 1) then
    
        ! call device_Nonlinear(6, vel(:,:,:,2), nl(:,:,:,2))
    
      else
    
            DO k = S32, N32
               DO j = S22, N22
                  DO i = S12, N12
                     dd1 = cp3(b3L,k)*vel(i,j,k+b3L,2)
    !pgi$ unroll = n:8
                     DO kk = b3L+1, b3U
                        dd1 = dd1 + cp3(kk,k)*vel(i,j,k+kk,2)
                     END DO
                     
                     nl(i,j,k,2) = nl(i,j,k,2) + dd1*pp(i,j,k)
                  END DO
               END DO
            END DO
    
      endif

  END IF
  !===========================================================================================================
  END IF
  
  
  
  
  
  IF (dimens == 3) THEN
  !===========================================================================================================
  !=== u*dw/dx ===============================================================================================
  !===========================================================================================================
  CALL interpolate2_pre_vel(.FALSE.,3,work1,pp)
  !-----------------------------------------------------------------------------------------------------------
  IF (upwind_yes) THEN

  if (GPU_accelerated == 1) then 
              call device_Nonlinear_upwind(7, vel(:,:,:,3), nl(:,:,:,3))
  else



     DO k = S33, N33
        DO j = S23, N23
           DO i = S13, N13
              IF (pp(i,j,k) .GE. 0.) THEN
                 dd1 = cNp1U(n1L,i)*vel(i+n1L,j,k,3)
!pgi$ unroll = n:8
                 DO ii = n1L+1, n1U
                    dd1 = dd1 + cNp1U(ii,i)*vel(i+ii,j,k,3)
                 END DO
              ELSE
                 dd1 = cNp1D(n1L,i)*vel(i+n1L,j,k,3)
!pgi$ unroll = n:8
                 DO ii = n1L+1, n1U
                    dd1 = dd1 + cNp1D(ii,i)*vel(i+ii,j,k,3)
                 END DO
              END IF
              
              nl(i,j,k,3) = dd1*pp(i,j,k)
           END DO
        END DO
     END DO
    endif

  ELSE

      if (GPU_accelerated == 1) then
    
        ! call device_Nonlinear(7, vel(:,:,:,3), nl(:,:,:,3))
    
      else
    
            DO k = S33, N33
               DO j = S23, N23
                  DO i = S13, N13
                     dd1 = cp1(b1L,i)*vel(i+b1L,j,k,3)
    !pgi$ unroll = n:8
                     DO ii = b1L+1, b1U
                        dd1 = dd1 + cp1(ii,i)*vel(i+ii,j,k,3)
                     END DO
                     
                     nl(i,j,k,3) = dd1*pp(i,j,k)
                  END DO
               END DO
            END DO
    
      endif

  END IF
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== v*dw/dy ===============================================================================================
  !===========================================================================================================
  CALL interpolate2_pre_vel(.FALSE.,3,work2,pp)
  !-----------------------------------------------------------------------------------------------------------
  IF (upwind_yes) THEN

    if (GPU_accelerated == 1) then 
              call device_Nonlinear_upwind(8, vel(:,:,:,3), nl(:,:,:,3))
   else

     DO k = S33, N33
        DO j = S23, N23
           DO i = S13, N13
              IF (pp(i,j,k) .GE. 0.) THEN
                 dd1 = cNp2U(n2L,j)*vel(i,j+n2L,k,3)
!pgi$ unroll = n:8
                 DO jj = n2L+1, n2U
                    dd1 = dd1 + cNp2U(jj,j)*vel(i,j+jj,k,3)
                 END DO
              ELSE
                 dd1 = cNp2D(n2L,j)*vel(i,j+n2L,k,3)
!pgi$ unroll = n:8
                 DO jj = n2L+1, n2U
                    dd1 = dd1 + cNp2D(jj,j)*vel(i,j+jj,k,3)
                 END DO
              END IF
              
              nl(i,j,k,3) = nl(i,j,k,3) + dd1*pp(i,j,k)
           END DO
        END DO
     END DO
   endif

  ELSE

      if (GPU_accelerated == 1) then
    
        ! call device_Nonlinear(8, vel(:,:,:,3), nl(:,:,:,3))
    
      else
    
            DO k = S33, N33
               DO j = S23, N23
                  DO i = S13, N13
                     dd1 = cp2(b2L,j)*vel(i,j+b2L,k,3)
    !pgi$ unroll = n:8
                     DO jj = b2L+1, b2U
                        dd1 = dd1 + cp2(jj,j)*vel(i,j+jj,k,3)
                     END DO
                     
                     nl(i,j,k,3) = nl(i,j,k,3) + dd1*pp(i,j,k)
                  END DO
               END DO
            END DO
 
      endif

  END IF


  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== w*dw/dz ===============================================================================================
  !===========================================================================================================
  CALL interpolate2_pre_vel(.FALSE.,3,work3,pp)
  !-----------------------------------------------------------------------------------------------------------
  IF (upwind_yes) THEN

  if (GPU_accelerated == 1) then 
              call device_Nonlinear_upwind(9, vel(:,:,:,3), nl(:,:,:,3))
  else


     DO k = S33, N33
        DO j = S23, N23
           DO i = S13, N13
              IF (pp(i,j,k) .GE. 0.) THEN
                 dd1 = cNw3U(n3L,k)*vel(i,j,k+n3L,3)
!pgi$ unroll = n:8
                 DO kk = n3L+1, n3U
                    dd1 = dd1 + cNw3U(kk,k)*vel(i,j,k+kk,3)
                 END DO
              ELSE
                 dd1 = cNw3D(n3L,k)*vel(i,j,k+n3L,3)
!pgi$ unroll = n:8
                 DO kk = n3L+1, n3U
                    dd1 = dd1 + cNw3D(kk,k)*vel(i,j,k+kk,3)
                 END DO
              END IF
              
              nl(i,j,k,3) = nl(i,j,k,3) + dd1*pp(i,j,k)
           END DO
        END DO
     END DO
   endif

  ELSE

      if (GPU_accelerated == 1) then
    
        ! call device_Nonlinear(9, vel(:,:,:,3), nl(:,:,:,3))
    
      else
    
    
            DO k = S33, N33
               DO j = S23, N23
                  DO i = S13, N13
                     dd1 = cw3(b3L,k)*vel(i,j,k+b3L,3)
    !pgi$ unroll = n:8
                     DO kk = b3L+1, b3U
                        dd1 = dd1 + cw3(kk,k)*vel(i,j,k+kk,3)
                     END DO
                     
                     nl(i,j,k,3) = nl(i,j,k,3) + dd1*pp(i,j,k)
                  END DO
               END DO
            END DO
    
      endif


  END IF
  !===========================================================================================================
  END IF
  
  
  END SUBROUTINE nonlinear
  
  
  
  
  
  
  
  
  
  ! TEST!!! anderes Modul?
  ! TEST!!! basiert neu auf interpolate2_vel_pre ... ok??
  SUBROUTINE interpolate_vel(exch_yes)
  
  IMPLICIT NONE
  
  LOGICAL, INTENT(in   ) ::  exch_yes
  
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: -                                                                                           !
  !----------------------------------------------------------------------------------------------------------!
  
  
                   CALL interpolate2_vel_pre(exch_yes,1,vel(b1L,b2L,b3L,1),work1)
                   CALL interpolate2_vel_pre(exch_yes,2,vel(b1L,b2L,b3L,2),work2)
  IF (dimens == 3) CALL interpolate2_vel_pre(exch_yes,3,vel(b1L,b2L,b3L,3),work3)
  
  
  CALL exchange(1,0,work1)
  CALL exchange(2,0,work1)
  CALL exchange(3,0,work1)
  
  CALL exchange(1,0,work2)
  CALL exchange(2,0,work2)
  CALL exchange(3,0,work2)
  
  IF (dimens == 3) THEN
     CALL exchange(1,0,work3)
     CALL exchange(2,0,work3)
     CALL exchange(3,0,work3)
  END IF
  
  
  END SUBROUTINE interpolate_vel
  
  
  
  
  
  
  
  
  
  
  ! TEST!!! noch relativ ungetestet!
  SUBROUTINE interpolate_pre_vel(exch_yes,m,SS1,SS2,SS3,NN1,NN2,NN3,phi,inter)
  
  IMPLICIT NONE
  
  LOGICAL, INTENT(in   ) ::  exch_yes
  INTEGER, INTENT(in   ) ::  m
  INTEGER, INTENT(in   ) ::  SS1, SS2, SS3, NN1, NN2, NN3
  
  REAL(8)   , INTENT(inout) ::  phi  (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL(8)   , INTENT(inout) ::  inter(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - bei kompakter Differenzierung werden immer ganze Linien behandelt.                        !
  !              - die Punkte auf (bzw. hinter) dem Rand der Wand-normalen Komponente werden auch bei        !
  !                kompakter Differenzierung im Feld immer explizit gerechnet, um nur eine Variante          !
  !                kompakter Differenzen abspeichern zu muessen (Extrapolation!).                            !
  !              - bei kompakter Differenzierung werden immer ganze Linien behandelt.                        !
  !----------------------------------------------------------------------------------------------------------!
  
  
  IF (exch_yes) CALL exchange2(m,0,SS1,SS2,SS3,NN1,NN2,NN3,phi)
  
  
  !===========================================================================================================
  IF (m == 1) THEN
        DO k = SS3, NN3
           DO j = SS2, NN2
              DO i = S11B, N11B ! TEST!!! hier koennte man auch SS1, NN1 nehmen! gilt auch fuer andere Routinen!!
                 inter(i,j,k) = cIpu(g1L,i)*phi(i+g1L,j,k)
                 DO ii = g1L+1, g1U
                    inter(i,j,k) = inter(i,j,k) + cIpu(ii,i)*phi(i+ii,j,k)
                 END DO
              END DO
           END DO
        END DO
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (m == 2) THEN
        DO k = SS3, NN3
           DO j = S22B, N22B ! TEST!!! hier koennte man auch SS2, NN2 nehmen!
              DO i = SS1, NN1
                 inter(i,j,k) = cIpv(g2L,j)*phi(i,j+g2L,k)
                 DO jj = g2L+1, g2U
                    inter(i,j,k) = inter(i,j,k) + cIpv(jj,j)*phi(i,j+jj,k)
                 END DO
              END DO
           END DO
        END DO
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (m == 3) THEN
        DO k = S33B, N33B ! TEST!!! hier koennte man auch SS3, NN3 nehmen!
           DO j = SS2, NN2
              DO i = SS1, NN1
                 inter(i,j,k) = cIpw(g3L,k)*phi(i,j,k+g3L)
                 DO kk = g3L+1, g3U
                    inter(i,j,k) = inter(i,j,k) + cIpw(kk,k)*phi(i,j,k+kk)
                 END DO
              END DO
           END DO
        END DO
  END IF
  !===========================================================================================================
  
  
  END SUBROUTINE interpolate_pre_vel
  
  
  
  
  
  
  
  
  
  
  
  ! Wie interpolate_pre_vel, allerdings mit fixen Index-Limiten (ohne Rand)
  SUBROUTINE interpolate2_pre_vel(exch_yes,m,phi,inter)
  
  IMPLICIT NONE
  
  LOGICAL, INTENT(in   ) ::  exch_yes
  INTEGER, INTENT(in   ) ::  m
  
  REAL(8)   , INTENT(inout) ::  phi  (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL(8)   , INTENT(inout) ::  inter(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  
  
  IF (exch_yes) CALL exchange(m,0,phi)
  
 
 
  !===========================================================================================================
  IF (m == 1) THEN

      if ((GPU_accelerated == 1) .and. (dimens == 3)) then
         call device_interpolate(m, phi, inter) 
      else
    
           DO k = S31, N31
               DO j = S21, N21
                  DO i = S11, N11
                     inter(i,j,k) = cIpu(g1L,i)*phi(i+g1L,j,k)
                     DO ii = g1L+1, g1U
                        inter(i,j,k) = inter(i,j,k) + cIpu(ii,i)*phi(i+ii,j,k)
                     END DO
                  END DO
               END DO
            END DO       
      endif

  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (m == 2) THEN

      if ((GPU_accelerated == 1) .and. (dimens == 3)) then
         call device_interpolate(m, phi, inter)   
      else
    
    
           DO k = S32, N32
               DO j = S22, N22
                  DO i = S12, N12
                     inter(i,j,k) = cIpv(g2L,j)*phi(i,j+g2L,k)
                     DO jj = g2L+1, g2U
                        inter(i,j,k) = inter(i,j,k) + cIpv(jj,j)*phi(i,j+jj,k)
                     END DO
                  END DO
               END DO
            END DO
      endif
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (m == 3) THEN

      if ((GPU_accelerated == 1) .and. (dimens == 3)) then
         call device_interpolate(m, phi, inter)   
      else
    
            DO k = S33, N33
               DO j = S23, N23
                  DO i = S13, N13
                     inter(i,j,k) = cIpw(g3L,k)*phi(i,j,k+g3L)
                     DO kk = g3L+1, g3U
                        inter(i,j,k) = inter(i,j,k) + cIpw(kk,k)*phi(i,j,k+kk)
                     END DO
                  END DO
               END DO
            END DO
      endif
  END IF
  !===========================================================================================================
  
  
  END SUBROUTINE interpolate2_pre_vel
  
  
  
  
  
  
  
  
  
  
  
  ! TEST!!! noch relativ ungetestet!
  SUBROUTINE interpolate_vel_pre(exch_yes,m,SS1,SS2,SS3,NN1,NN2,NN3,phi,inter)
  
  IMPLICIT NONE
  
  LOGICAL, INTENT(in   ) ::  exch_yes
  INTEGER, INTENT(in)    ::  m
  INTEGER, INTENT(in   ) ::  SS1, SS2, SS3, NN1, NN2, NN3
  
  REAL(8)   , INTENT(inout) ::  phi  (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL(8)   , INTENT(inout) ::  inter(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - bei kompakter Differenzierung werden immer ganze Linien behandelt.                        !
  !----------------------------------------------------------------------------------------------------------!
  
  
  IF (exch_yes) CALL exchange2(m,m,SS1,SS2,SS3,NN1,NN2,NN3,phi)
  
  
  !===========================================================================================================
  IF (m == 1) THEN
        DO k = SS3, NN3
           DO j = SS2, NN2
              DO i = S1p, N1p
                 inter(i,j,k) = cIup(d1L,i)*phi(i+d1L,j,k)
                 DO ii = d1L+1, d1U
                    inter(i,j,k) = inter(i,j,k) + cIup(ii,i)*phi(i+ii,j,k)
                 END DO
              END DO
           END DO
        END DO
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (m == 2) THEN
        DO k = SS3, NN3
           DO j = S2p, N2p
              DO i = SS1, NN1
                 inter(i,j,k) = cIvp(d2L,j)*phi(i,j+d2L,k)
                 DO jj = d2L+1, d2U
                    inter(i,j,k) = inter(i,j,k) + cIvp(jj,j)*phi(i,j+jj,k)
                 END DO
              END DO
           END DO
        END DO
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (m == 3) THEN
        DO k = S3p, N3p
           DO j = SS2, NN2
              DO i = SS1, NN1
                 inter(i,j,k) = cIwp(d3L,k)*phi(i,j,k+d3L)
                 DO kk = d3L+1, d3U
                    inter(i,j,k) = inter(i,j,k) + cIwp(kk,k)*phi(i,j,k+kk)
                 END DO
              END DO
           END DO
        END DO
  END IF
  !===========================================================================================================
  
  
  END SUBROUTINE interpolate_vel_pre
  
  
  
  
  
  
  
  
  
  
  
  ! Wie interpolate_vel_pre, allerdings mit fixen Index-Limiten (ohne Rand)
  SUBROUTINE interpolate2_vel_pre(exch_yes,m,phi,inter)
  
  IMPLICIT NONE
  
  LOGICAL, INTENT(in   ) ::  exch_yes
  INTEGER, INTENT(in)    ::  m
  
  REAL(8)   , INTENT(inout) ::  phi  (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL(8)   , INTENT(inout) ::  inter(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  !for unit testing
!  REAL(8)    ::  phi_test  (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
!  REAL(8)    ::  inter_test(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
 

  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  
  integer :: counter
  
  IF (exch_yes) CALL exchange(m,m,phi)
  
  
  !===========================================================================================================
  IF (m == 1) THEN

!  phi_test=phi
!  inter_test=inter
!  counter = 0
  if (GPU_accelerated == 1) then
     call device_interpolate2(m,phi,inter)
  else
        DO k = S3p, N3p
           DO j = S2p, N2p
              DO i = S1p, N1p
                 inter(i,j,k) = cIup(d1L,i)*phi(i+d1L,j,k)
                 DO ii = d1L+1, d1U
                    inter(i,j,k) = inter(i,j,k) + cIup(ii,i)*phi(i+ii,j,k)
                 END DO
                 
!                    if (abs(inter_test(i,j,k) - inter(i,j,k)) .ge. 0.000000000000001) print*, 'mismatch', inter_test(i,j,k), inter(i,j,k) 
              END DO
           END DO
        END DO
  endif
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (m == 2) THEN

  if (GPU_accelerated == 1) then
     call device_interpolate2(m,phi,inter)
  else
       DO k = S3p, N3p
           DO j = S2p, N2p
              DO i = S1p, N1p
                 inter(i,j,k) = cIvp(d2L,j)*phi(i,j+d2L,k)
                 DO jj = d2L+1, d2U
                    inter(i,j,k) = inter(i,j,k) + cIvp(jj,j)*phi(i,j+jj,k)
                 END DO
              END DO
           END DO
        END DO
  endif
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (m == 3) THEN
 
  if (GPU_accelerated == 1) then
     call device_interpolate2(m,phi,inter)
  else
        DO k = S3p, N3p
           DO j = S2p, N2p
              DO i = S1p, N1p
                 inter(i,j,k) = cIwp(d3L,k)*phi(i,j,k+d3L)
                 DO kk = d3L+1, d3U
                    inter(i,j,k) = inter(i,j,k) + cIwp(kk,k)*phi(i,j,k+kk)
                 END DO
              END DO
           END DO
        END DO
  endif
  END IF
  !===========================================================================================================
  
  
  END SUBROUTINE interpolate2_vel_pre
  
  
  
  
  
  
  
  
  
  
  
 
  
  
  
  
  
 
  
  
  
  
  
  
  
 
  
  
  
  
  
 
  
  
END MODULE mod_diff
