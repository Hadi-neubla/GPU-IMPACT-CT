!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!* modified by Barna Becske, ARTORG CVE, University of Bern (barna.becsek@artorg.unibe.ch)                   *
!* October 2014                                                                                              *
!*************************************************************************************************************

!> module containing subroutines for computing the rhs. It uses modules mod_dims, mod_vars, mod_diff, 
!! mod_particles and mod_les
MODULE mod_rhs
  
  
  USE mod_dims
  USE mod_vars
  USE mod_diff
  !USE mpi  
  
  PRIVATE
  
  PUBLIC rhs_vel
  
  INCLUDE 'mpif.h'
  
  CONTAINS
  
!pgi$g unroll = n:8
!!pgi$r unroll = n:8
!!pgi$l unroll = n:8
  
  
  !> computes the rhs of the velocity
  SUBROUTINE rhs_vel
  
  IMPLICIT NONE
  
  INTEGER                ::  i, j, k, m
  REAL(8)                   ::  mult
  REAL(8)                   ::  flux, flux_global
 

  
  ! TEST!!! generell muesste eigentlich nur der "innere" Bereich (d.h. ohne Raender) von rhs" belegt werden!
  
  !===========================================================================================================
  !=== rhs = rhs + Helmholtz(vel) ============================================================================
  !===========================================================================================================
  !TODO: port to GPU_rhs
  IF (timeint_mode == 1 .OR. thetaL == 1.) THEN
     DO m = 1, dimens
        IF (m == 1) rhs(S11B:N11B,S21B:N21B,S31B:N31B,1) = vel(S11B:N11B,S21B:N21B,S31B:N31B,1)
        IF (m == 2) rhs(S12B:N12B,S22B:N22B,S32B:N32B,2) = vel(S12B:N12B,S22B:N22B,S32B:N32B,2)
        IF (m == 3) rhs(S13B:N13B,S23B:N23B,S33B:N33B,3) = vel(S13B:N13B,S23B:N23B,S33B:N33B,3)
        
     END DO
 END IF
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== rhs = rhs - nl_old ====================================================================================
  !===========================================================================================================
  !TODO: port to GPU_rhs
  IF (substep /= 1) THEN
     mult = dtime*bRK(substep)
     
     DO m = 1, dimens
        IF (m == 1) rhs(S11B:N11B,S21B:N21B,S31B:N31B,1) = rhs(S11B:N11B,S21B:N21B,S31B:N31B,1) - mult*nl(S11B:N11B,S21B:N21B,S31B:N31B,1)
        IF (m == 2) rhs(S12B:N12B,S22B:N22B,S32B:N32B,2) = rhs(S12B:N12B,S22B:N22B,S32B:N32B,2) - mult*nl(S12B:N12B,S22B:N22B,S32B:N32B,2)
        IF (m == 3) rhs(S13B:N13B,S23B:N23B,S33B:N33B,3) = rhs(S13B:N13B,S23B:N23B,S33B:N33B,3) - mult*nl(S13B:N13B,S23B:N23B,S33B:N33B,3)
     END DO
  END IF
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== rhs = rhs - nl ========================================================================================
  !===========================================================================================================
  !--- advective terms ---------------------------------------------------------------------------------------
   
  !--- viscose terms -----------------------------------------------------------------------------------------
  IF (timeint_mode == 1 .AND. (.NOT. Euler_yes)) THEN
     multL = 1./Re
     CALL Helmholtz_explicit(.FALSE.)
  END IF
  
  !--- forcing -----------------------------------------------------------------------------------------------
  !TODO port to GPU_rhs : This has to be open to customize
  CALL forcing_vel
  
  !--- rhs = rhs - nl ----------------------------------------------------------------------------------------
  mult = dtime*aRK(substep)
  !TODO port to GPU_rhs
  DO m = 1, dimens
     IF (m == 1) rhs(S11B:N11B,S21B:N21B,S31B:N31B,1) = rhs(S11B:N11B,S21B:N21B,S31B:N31B,1) - mult*nl(S11B:N11B,S21B:N21B,S31B:N31B,1)
     IF (m == 2) rhs(S12B:N12B,S22B:N22B,S32B:N32B,2) = rhs(S12B:N12B,S22B:N22B,S32B:N32B,2) - mult*nl(S12B:N12B,S22B:N22B,S32B:N32B,2)
     IF (m == 3) rhs(S13B:N13B,S23B:N23B,S33B:N33B,3) = rhs(S13B:N13B,S23B:N23B,S33B:N33B,3) - mult*nl(S13B:N13B,S23B:N23B,S33B:N33B,3)
  END DO
  !===========================================================================================================
  
  
  
    
  END SUBROUTINE rhs_vel
  
  
  
  
END MODULE mod_rhs
