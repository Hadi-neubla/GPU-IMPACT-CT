!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!* modified by Barna Becsek, ARTORG CVE, University of Bern (barna.becsek@arotrg.unibe.ch)                   *
!* October 2014                                                                                              *
!*************************************************************************************************************
 
!> @file usr_force.f90
!! file containing additional forcing terms for governing equaitons

!pgi$g unroll = n:8
!!pgi$r unroll = n:8
!!pgi$l unroll = n:8
  
  
  !> volume forcing of the momentum equation
  !! @note This could serve as the (most direct/straight-forward) interface between structural and fluid
  !!       solvers.
  SUBROUTINE forcing_vel
  ! (basic subroutine)
  
  USE mod_dims
  USE mod_vars
  USE usr_vars
  USE usr_func
!  USE mod_ibm !bbecsek
  
  IMPLICIT NONE
  
  INTEGER                ::  i, j, k, m
  INTEGER                ::  ii, jj, kk
  REAL(8)                   ::  lamb_fringe,parab,time_fact
  integer,parameter      ::  bileaflet=0, mixing=0, sphere_forcing=0
  REAL(8)                   ::  y1,y2,y3,y4,y5,y6,y7,y8
  real(8)                   ::  v1, v2, v3

!  real(8)                   ::  Qx_damp(1:(N11+1), 1:(N22+1), 1:(N33+1)), Qy_damp(1:(N11+1), 1:(N22+1), 1:(N33+1)), Qx_bdamp(1:(N11+1), 1:(N22+1), 1:(N33+1)), Qy_bdamp(1:(N11+1), 1:(N22+1), 1:(N33+1))        ! q and \bar{q} in the selective damping formulation

  real(8)                   ::  Qx_bdamp, Qy_bdamp
  real(8)                   ::  Delta_f, c_gain   ! filter width (Delta=1/omega_c) and control gain (x) in the formulation



  !--- additional volume forces in the momentum equation ---
  ! note: - du/dt = RHS-nl
  !       - cf. sketch in file "usr_geometry.f90"
  !
  !       grid points in the domain
  !       |       |       |    velocity component
  !       |       |       |    |
  ! nl(S11:N11,S21:N21,S31:N31,1)
  ! nl(S12:N12,S22:N22,S32:N32,2)
  ! nl(S13:N13,S23:N23,S33:N33,3)
  !

  !fd = 0.

 !  ---------------------------------------------------------------------------------------------
 !  --------------------------Selective Frequency Damping ---------------------------------------

   
 Delta_f = 0.03
 c_gain = 1.0
 
 ! at the first substep of RK we do Euler time marching for q
 ! first update Q_bar, then add the relaxation to rhs
 ! first dir
!  DO k = 1, N31+1
!     DO j = 1, N21+1
!        DO i = 1, N11+1
!        !    CALL poiseuille_parabola(0.0,L2,x2p(j),parab)
!           lamb_fringe =0.0
!           if (((x1p(i) .gt. 0.8 ) .and. (x2p(j) .le. 5.*L2/6. )) ) then
!              if (substep == 1) then
!                 if (abs(vel(i,j,k,1)) .gt. 0.0) then
!                    Qx_bdamp = (Delta_f/(Delta_f+dtime)) * (Qx_bdamp + dtime/Delta_f*vel(i, j, k,1))
!                 endif
!                    nl(i, j, k,1) = nl(i, j, k,1) -c_gain * ((vel(i, j, k,1))- Qx_bdamp) 
!               
!                    ! second dir
!                    Qy_bdamp = (Delta_f/(Delta_f+dtime)) * (Qy_bdamp + dtime/Delta_f*vel(i, j, k,2))
!                    nl(i, j, k,2) = nl(i, j, k,2) -c_gain * ((vel(i, j, k,2))- Qy_bdamp) 
!             endif
!          endif
!        enddo
!      enddo
!    enddo       
 ! endif









 

  !---- Fringe forcing for Rami 2D simulations -  in the main channel forcing
  !---- HZ

  DO k = 1, N31+1
     DO j = 1, N21+1
        DO i = 1, N11+1
!       !   CALL fringe_coeff(10. , .5 , 3.5 ,0.25  ,0.25 , x1u(i) , lamb_fringe)
          CALL poiseuille_parabola(dble(0.0),L2,x2p(j),parab)
!           lamb_fringe =0.0
!           if ((x1p(i) .le. 0.5 )  .and. (x1p(i) .ge. 0.0) ) lamb_fringe = 20.0
!           if ((x1p(i) .le. 1.0 )  .and. (x1p(i) .ge. .1) ) lamb_fringe = 10.0
!
!           nl(i,j,k,1) = nl(i,j,k,1) -lamb_fringe * (( parab  - vel(i,j,k,1))) ! Poiseuille profile 
!            vel(i,j,k,1) = parab

!--- regurgitation MS star
!!!           if ((x1p(i) .le. 21. )  .and. (x1p(i) .ge. 1.0) ) then
!!!               nl(i,j,k,1) =  nl(i,j,k,1) - .05  !*smooth_step(time/20.)!
!!!            endif

!!          nl(i,j,k,1) = nl(i,j,k,1) - lamb_fringe*alpha*( (L1/ (( 6. )))*(2./Re)  )!*smooth_step(time/20.)! 
!          nl(i,j,k,1) = nl(i,j,k,1) - 8./Re 

!
      END DO
     END DO
  END DO


  !--- Forcing inside the restriction

  if (bileaflet == 1) then


  DO k = 1, N31+1
     DO j = 1, N21+1
        DO i = 1, N11+1
        !    CALL poiseuille_parabola(0.0,L2,x2p(j),parab)
           lamb_fringe =0.0
           if ((x1p(i) .le. 0.6 )  .and. (x1p(i) .ge. 0.2) .and. (sqrt((x2p(j) -0.5)**2  + (x3p(k) -0.5)**2) .le. 0.33))   lamb_fringe = 20.0
           
               nl(i,j,k,1) = nl(i,j,k,1) -lamb_fringe * (( 2.0  - vel(i,j,k,1))) ! Poiseuille profile 
              ! nl(i,j,k,2) = nl(i,j,k,2) -lamb_fringe * (( 0.0  - vel(i,j,k,2))) ! Poiseuille profile 
              ! nl(i,j,k,3) = nl(i,j,k,3) -lamb_fringe * (( 0.0  - vel(i,j,k,3))) ! Poiseuille profile 
              
               ! nl(i,j,k,1) =  nl(i,j,k,1) - lamb_fringe/20.0 * (8./Re)

       END DO
     END DO
  END DO

  endif



!--- forcing for the sphere case
   if (sphere_forcing ==1) then

        do ii=S1p,N1p
           do jj =S2p, N2p
              do kk =S3p, N3p
                  if ((x3p(kk) .ge. 0.0) .and. (x3p(kk) .le. 0.6)) then
                      nl(ii,jj,kk,3) = nl(ii,jj,kk,3) - 10. *(1. - vel(ii,jj,kk,3))
                      nl(ii,jj,kk,2) = nl(ii,jj,kk,2) - 10. *(0. - vel(ii,jj,kk,2))
                      nl(ii,jj,kk,1) = nl(ii,jj,kk,1) - 10. *(0. - vel(ii,jj,kk,1))
                  endif
              enddo
           enddo
        enddo


   endif 


  







IF(IB_on) THEN
!  CALL spread_force_dens_to_vel_grid
  
  !nl = nl + ( fd/( ((mu_blood**2.) * (Re**2.))/rho_blood ) )!*(L_ref**2.) ! is the L_ref necessary?
  !write(*,*) 'in process',rank,'fd=',fd

  END IF

!  nl = nl + fd * L_ref/(rho_blood*U_ref**2.)

  !=== FRINGE FORCING =======================================================================================
  IF (fringe_yes) THEN
    CALL apply_fringe_forcing
  END IF
  !==========================================================================================================

  !=== WINDKESSEL LOADING ===================================================================================
  IF (WK_yes) THEN
    CALL windkessel_integration_step
    CALL apply_windkessel_loading
  END IF
  !==========================================================================================================

  !*** debugging
  !write(*,*)'x1p(65)=',x1p(65),' x1u(64)=',x1u(64)

  !DO k = S31, N31
  !   DO j = S21, N21
  !      DO i = S11, N11
  !         nl(i,j,k,1) = nl(i,j,k,1) - 1.*(0. - vel(i,j,k,1))*interface((x2p(j)-0.8*L2)/(0.5*L1))
  !      END DO
  !   END DO
  !END DO
  
  
  END SUBROUTINE forcing_vel
