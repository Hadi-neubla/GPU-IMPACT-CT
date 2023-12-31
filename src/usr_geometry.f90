!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!* modified by Barna Becsek, ARTORG CVE, University of Bern (barna.becsek@artorg.unibe.ch)                   *
!* October 2014                                                                                              *
!* GPU version by Hadi Zolfaghari , ARTORG CVE, then DAMTP, Cambridge University (hz382@cam.ac.uk)           *
!* Oct 2015 - Sep 2023                                                                                       *
!*************************************************************************************************************

 
!pgi$g unroll = n:8
!!pgi$r unroll = n:8
!!pgi$l unroll = n:8
  
  
  !--- grid, index and index-range convention---
  !
  ! coordinates and indices on staggered grid (2D):
  !   - boundary conditions are set on pressure grid points
  !   - assignments to grid points outside of the given index ranges are ignored or overwritten
  !   - x1p, x1u, x2p, x2v, x3p, x3w: local  coordinates (inside each processor-block)
  !   - y1p, y1u, y2p, y2v, y3p, y3w: global coordinates (include all blocks)
  !
  !
  !    O: pressure, concentrations
  !    >: vel(:,:,:1)
  !    v: vel(:,:,:2)
  !   //: boundary
  !
  !               y1p(1)  y1p(2)   ...    y1p(M1)
  !           y1u(0)  y1u(1)  y1u(2)   ...    y1u(M1)
  !
  ! x2v(N22B)       v       v       v       v       y2v(M2)
  !               //|///////|///////|///////|//
  ! x2p(N2p)    >---O--->---O--->---O--->---O--->   y2p(M2)
  !  .            //|       |       |       |//      .
  !  .            //v       v       v       v//      .
  !  .            //|       |       |       |//      .
  ! x2p(S1p +2) >---O--->---O--->---O--->---O--->   y2p(3)
  !               //|       |       |       |//
  ! x2v(S22B+2)   //v       v       v       v//     y2v(2)
  !               //|       |       |       |//
  ! x2p(S1p +1) >---O--->---O--->---O--->---O--->   y2p(2)
  !               //|       |       |       |//
  ! x2v(S22B+1)   //v       v       v       v//     y2v(1)
  !               //|       |       |       |//
  ! x2p(S1p)    >---O--->---O--->---O--->---O--->   y2p(1)
  !               //|///////|///////|///////|//
  ! x2v(S22B)       v       v       v       v       y2v(0)
  !
  !              x1p(S1p) x1p(S1p+1) ...  x1p(N1p)
  !         x1u(S11B) x1u(S11B+1) x1u(S11B+2)...x1u(N11B)
  !
  !
  !
  ! pressure / concentrations:
  !
  !  start/end index
  !  |spatial direction
  !  ||
  !  S1p
  !  N1p
  !
  !
  ! velocity:
  !
  !  start/end index
  !  |spatial direction
  !  ||velocity component
  !  |||includes/excludes the outermost grid points on the boundary or outside the domain
  !  ||||                             (periodicity/symmetry is not treated as a boundary)
  !  S11B
  !  N11B
  !  S11
  !  N11
  !
  !
  ! Note that
  !   S12B = S13B = S1p
  !   N12B = N13B = N1p
  !   S21B = S23B = S2p
  !   N21B = N23B = N2p
  !   S31B = S32B = S3p
  !   N31B = N32B = N3p
  !
  
  
  !> subroutine that sets up the global coordinates (2 coordinates in each direction, 1 for pressure, 1 for one 
  !! velocity component). Velocity coordinate components in the dimension that correspond to the velocity
  !! component in question are shifted by 1/2 (in between nodes). It uses modules mod_dims, mod_vars, usr_vars
  !! and usr_func.
  SUBROUTINE get_coords
  ! (basic subroutine)
  
  USE mod_dims
  USE mod_vars
  USE usr_vars
  USE usr_func
  
  IMPLICIT NONE
  
  INTEGER                ::  i, j, k
  REAL(8)                ::  M_dp, i_dp, j_dp, k_dp 
  
  !--- specify global coordinates and grid spacings ---
  ! note: - local processor-block coordinates and grid spacings are automatically derived from global grid
  !       - dy3p = dy3w = 1. for 2D (may simplify programming)
  !       - ensure that for all i
  !            y1p(i) < y1p(i+1)
  !            y1u(i) < y1u(i+1)
  !            y1p(i) < y1u(i) < y1p(i+1)
  !               etc.
  !       - code is tested only for
  !            y1p(1 ) = 0.
  !            y1p(M1) = L1
  !               etc.
  !
  ! example:
  

  DO i = 1, M1
     M_dp = dble(M1) 
     i_dp = dble(i) 
     CALL coord_tanh(L1,M_dp,dble(0.0),dble(0.0),i_dp,y1p(i),dy1p(i))
  END DO
  DO i = 0, M1
     CALL coord_tanh(L1,dble(M1*1.),dble(0.0),dble(0.0),dble(i*1.)+0.5,y1u(i),dy1u(i))
  END DO
  y1p = y1p - y1_origin
  y1u = y1u - y1_origin

  
  DO j = 1, M2
     CALL coord_tanh(L2,dble(M2*1.),dble(0.0),dble(0.0),dble(j*1.)    ,y2p(j),dy2p(j))
  END DO
  DO j = 0, M2
     CALL coord_tanh(L2,dble(M2*1.),dble(0.0),dble(0.0),dble(j*1.)+0.5,y2v(j),dy2v(j))
  END DO
  y2p = y2p - y2_origin
  y2v = y2v - y2_origin
  
  
  IF (M3 /= 2) THEN
  DO k = 1, M3
     CALL coord_tanh(L3,dble(M3*1.),dble(0.0),dble(0.0),dble(k*1.)    ,y3p(k),dy3p(k))
  END DO
  DO k = 0, M3
     CALL coord_tanh(L3,dble(M3*1.),dble(0.0),dble(0.0),dble(k*1.)+0.5,y3w(k),dy3w(k))
  END DO
  y3p = y3p - y3_origin
  y3w = y3w - y3_origin
  END IF
  
  
  END SUBROUTINE get_coords


   SUBROUTINE coord_stretch(direction_id,x_center,Lmax,lambda_stretch,N_points,xx,dx)
  USE mod_dims
  USE mod_vars
  USE usr_vars
  USE usr_func
 
  IMPLICIT NONE

  INTEGER, INTENT(in)    ::  direction_id
  REAL(8)   , INTENT(in)    ::  x_center
  REAL(8)   , INTENT(in)    ::  lambda_stretch

  INTEGER, INTENT(in)    ::  N_points

  REAL(8)   , INTENT(in)    ::  Lmax
  REAL(8)   , INTENT(out)   ::  xx(M1)
  REAL(8)   , INTENT(out)   ::  dx(M1)

  Real                   ::  psi_array(N_points)  !equidistant array on unit length
  integer                ::  i

  xx(1) = 0.
  psi_array(1) = 0.
  DO i=2,N_points
     psi_array(i) =  1./(N_points-1) * (i-1)
     xx(i) = Lmax * psi_array(i) + lambda_stretch * (x_center - Lmax*psi_array(i)) * (1.-psi_array(i)) * psi_array(i)
  !   print*, xx(i)
     dx(i-1) = xx(i) - xx(i-1)
  ENDDO
  dx(N_points) = dx(N_points-1)

  END SUBROUTINE coord_stretch
