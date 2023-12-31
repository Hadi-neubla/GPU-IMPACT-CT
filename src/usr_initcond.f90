!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!* modified by Barna Becsek, ARTORG CVE, University of Bern (barna.becsek@artorg.unibe.ch)                   *
!* October 2014                                                                                              *
!*************************************************************************************************************

!> @file usr_initcond.f90
!! file containing three subroutines for setting initial conditions:
!! one for velocities (and pressure).

!pgi$g unroll = n:8
!!pgi$r unroll = n:8
!!pgi$l unroll = n:8
  
  
  !> subroutine for setting initial conditions for the velocity and pressure fields.
  !! It uses mod_dims, mod_vars, usr_vars and usr_func.
  SUBROUTINE initial_conditions_vel
  ! (basic subroutine)
  
  USE mod_dims
  USE mod_vars
  USE usr_vars
  USE usr_func
  
  IMPLICIT NONE
  
  INTEGER                ::  i, j, k
  integer                ::  ii,jj,kk
  integer, parameter     ::  TGV_sim=0
  REAL(8)                   ::  parab
  
  
  !--- initial conditions for pressure ---
  ! note: - initialization is generally not necessary
  !       - specification provides only an initial guess for the first pressure iteration
  !       - cf. sketch in file "usr_geometry.f90"
  !
  !        grid points in the domain and on the boundary
  !        |       |       |
  ! pre(S1p:N1p,S2p:N2p,S3p:N3p)
  !
  pre = 0.
  
  
  !--- initial conditions for velocity ---
  ! note: - cf. sketch in file "usr_geometry.f90"
  !
  !         grid points in the domain and on the boundary
  !         |         |         |     velocity component
  !         |         |         |     |
  ! vel(S11B:N11B,S21B:N21B,S31B:N31B,1)
  ! vel(S12B:N12B,S22B:N22B,S32B:N32B,2)
  ! vel(S13B:N13B,S23B:N23B,S33B:N33B,3)
  !
  vel = 0.
  DO k = S31B, N31B
     DO j = S21B, N21B
        DO i = S11B, N11B
           !vel(i,j,k,2) = vel(i,j,k,2) + interface((x1p(i)-0.75*L1)/(0.01*L1))
           !vel(i,j,k,2) = vel(i,j,k,2) + SIN(2*pi*x1p(i)/L2)
           !CALL poiseuille_parabola(x2p(S2p),x2p(N2p),x2p(j),parab)

           !vel(i,j,k,1) = parab
        END DO
     END DO
  END DO
  
    pi = 3.14159265359

 if (TGV_sim == 1) then
    do ii=S1p, N1p
       do jj=S2p, N2p
          do kk=S3p, N3p
               vel(ii,jj,kk,1) = 1.*sin(x1p(ii)-pi)*cos(x2p(jj)-pi)*cos(x3p(kk)-pi)
               vel(ii,jj,kk,2) = -1.*cos(x1p(ii)-pi)*sin(x2p(jj)-pi)*cos(x3p(kk)-pi)
               vel(ii,jj,kk,3) = 0.0
          enddo
       enddo
    enddo
 endif



  !*** debugging
  !open(unit=58,file=('init_vel_field.txt'),status='unknown',recl=4002)
  !do i=S12B, N12B
  !write(58,*), (vel(i,j,1,2), j=S22B,N22B)
  !end do  


  !--- advection velocity ---
  ! note: - only boundary-normal velocity component is specified, i.e. velocity 
  !         components tangential to the boundary are not considered
  !       - must be specified for the pressure grid points (values for advection
  !         velocity of boundary-tangential velocity components are interpolated)
  !       - cf. sketch in file "usr_geometry.f90"
  !
  !      orientation of the boundary normal
  !      |    grid points on the boundary
  !      |    |       |     lower/upper boundary
  !      |    |       |     |
  ! drift1(S2p:N2p,S3p:N3p,1:2) = 1.
  ! drift2(S1p:N1p,S3p:N3p,1:2) = 1.
  ! drift3(S1p:N1p,S2p:N2p,1:2) = 1.
  !
  drift1 = 1.
  drift2 = 1.
  drift3 = 1.
  
  
  !--- specify vector for flux corrections ---
  ! note: - needs to be specified only if orthogonal projection is NOT used, i.e. for
  !         "nullspace_ortho_yes = .FALSE."
  !       - correction can be applied only at boundaries (not symmetry)
  !       - correction vector must not be orthogonal to the null space vector "psi_vel"
  !       - only the shape matters, not the amplitude or the sign
  !       - the correction vector will be normalized subsequently
  !       - the correction vector is applied only in sub-domains located at the corresponding boundaries
  !       - cf. sketch in file "usr_geometry.f90"
  !
  !   velocity component
  !   |orientation of the boundary normal
  !   ||     grid points on the boundary
  !   ||     |         |      lower/upper boundary
  !   ||     |         |      |
  ! th11(S21B:N21B,S31B:N31B,1:2)
  ! th12(S11B:N11B,S31B:N31B,1:2)
  ! th13(S11B:N11B,S21B:N21B,1:2)
  ! th21(S22B:N22B,S32B:N32B,1:2)
  ! th22(S12B:N12B,S32B:N32B,1:2)
  ! th23(S12B:N12B,S22B:N22B,1:2)
  ! th31(S23B:N23B,S33B:N33B,1:2)
  ! th32(S13B:N13B,S33B:N33B,1:2)
  ! th33(S13B:N13B,S23B:N23B,1:2)
  !
  th11 = 1.
  th12 = 0.
  th13 = 0.
  
  th21 = 0.
  th22 = 1.
  th23 = 0.
  
  th31 = 0.
  th32 = 0.
  th33 = 1.
  
  
  !--- additional, user-defined initializations ---
  !
  ! for statistics:
  energy_visc       = 0.
  diss_viscInt_old  = 0.
  
  
  END SUBROUTINE initial_conditions_vel
