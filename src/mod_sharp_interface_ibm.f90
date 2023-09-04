!Immersed Boundary Method for Rigid Body Handling in IMPACT
!by Hadi Zolfaghari

!***************************************************************************************************************

MODULE mod_sharp_interface_ibm




  USE mod_dims
  USE mod_vars
!  USE mod_geometry

  IMPLICIT NONE

        real(8), parameter       ::    level=.6
        real(8), parameter       ::    flat = 48.
        real(8), parameter       ::    flat1 =48.

!

        real(8)                  ::    x_cs,y_cs     ! coords of center of sinus arch
        real(8)                  ::    x_cc,y_cc     ! coords of the center of the comissure arch
        real(8)                  ::    r_sinus       ! radius of sinus arch
        real(8)                  ::    r_com         ! radius of comissure arch
        real(8)                  ::    x_ref, y_ref  ! refernce point coordinates
        real(8)                  ::    x2_s, y2_s
        real(8)                  ::    x3_s, y3_s
        real(8)                  ::    x1_c, y1_c
        real(8)                  ::    x2_c, y2_c
        real(8)                  ::    x3_c, y3_c
        real(8)                  ::    r_r
        real(8)                  ::    h1, h2
        real(8)                  ::    a1, b1, c1, a2, b2, c2
        real(8), parameter       ::    lf1 =0.15, lf2=0.1275 ,  theta1 = 1.047198, theta2 = 1.570796 


   contains

     subroutine parallel_shape_function(ycurve,yhcurve, ycurve_down, yhcurve_down)

        real(8)                                 ::    Lx, Ly, r, pi, h, d1,d2,d3,d4, amp, a, b, c, delta
        integer                              ::    n
        integer                              ::    i
        real(8)                                 ::    x1, x2, x3, x4, x5, x6
        integer                              ::    i1, i2, i3, i4, i5, i6
        real(8)                                 ::    m1, m2
        real(8)                                 ::    factor
        real(8)                                 ::    radi

        real(8), intent(out), dimension(N1p)     ::    ycurve, yhcurve
        real(8), intent(out), dimension(N1p)     ::    ycurve_down, yhcurve_down



!! The rami data



!-----This subroutine generates the surface function

        if (rank == 0)       open(unit=1,file='funcdata.dat')
        if (rank == 0)       open(unit=2,file='funcdata_down.dat')

!!!!!=========================Rami et. al 2D================================
!!!!!====================================================================!!!

!!!====Geometry Characteristics

      x_ref = L2/6.
!!      y_ref = 6.
      y_ref = 1.0

      r_r = L2/2. - x_ref
      h1 = 0.6*r_r
      h2 = 0.8*r_r


      x2_s = L2/2. - 1.4*r_r
      y2_s = y_ref + h1

      x3_s = L2/2. - 1.1*r_r
      y3_s = y_ref + h1 + h2

      x1_c = L2/2. + r_r
      y1_c = y_ref

      x2_c = L2/2. + 1.25*r_r
      y2_c = y_ref +h1

      x3_c = L2/2. + 1.1*r_r
      y3_c = y_ref + h1 + h2


!!!====Finding the archs's radi and centers

!!====For the sinus

      a1 = 2.*(x2_s - x_ref)
      b1 = 2.*(y2_s - y_ref)
      a2 = 2.*(x3_s - x_ref)
      b2 = 2.*(y3_s - y_ref)

      c1 = x2_s**2 + y2_s**2 - x_ref**2 - y_ref**2
      c2 = x3_s**2 + y3_s**2 - x_ref**2 - y_ref**2

      y_cs = (a1*c2 - a2*c1) / (a1*b2 - a2*b1)
      x_cs = (c1-b1*y_cs )/a1

      r_sinus = sqrt((x2_s-x_cs)**2 + (y2_s-y_cs)**2)

!!====For the comissures

      a1 = 2.*(x2_c - x1_c)
      b1 = 2.*(y2_c - y1_c)
      a2 = 2.*(x3_c - x1_c)
      b2 = 2.*(y3_c - y1_c)

      c1 = x2_c**2 + y2_c**2 - x1_c**2 - y1_c**2
      c2 = x3_c**2 + y3_c**2 - x1_c**2 - y1_c**2

      y_cc = (a1*c2 - a2*c1) / (a1*b2 - a2*b1)
      x_cc = (c1-b1*y_cc )/a1


      r_com = sqrt((x2_c-x_cc)**2 + (y2_c-y_cc)**2)


       m1 = (1.-flat/(N2p-1))*L2/(x1-x2)
       m2 = (1.-flat/(N2p-1))*L2/(x4-x3)


!!       x1 = 5.0
!       x2 = 5.0
!       x3 = 8.0
!       x4 = 8.0
!       x5 = y_ref
!       x6 = x5  + h1 + h2

       x1 = 0.0
       x2 = 0.0
       x3 = 5.0
       x4 = 5.0
       x5 = y_ref
       x6 = x5  + h1 + h2




       do i=1,N1p

          if (((x1p(i)) .ge. (x5)) .and. ((x1p(i)) .le. (x6))) then
   
              ycurve(i) = x_cc + sqrt(r_com**2-(x1p(i)-y_cc)**2)
              ycurve_down(i) = x_cs - sqrt(r_sinus**2-(x1p(i)-y_cs)**2)
   
          else if (((x1p(i)) .ge. x2) .and. ((x1p(i)) .lt. x5)) then

              ycurve(i) = x_ref + 2.*r_r
              ycurve_down(i) = x_ref


!---- The oval shape inlet

!             ycurve(i) = L2 - L2/3. * sqrt(1- ((x1p(i)-x5)/(x5-x1))**2)
!             ycurve_down(i) = L2/3. * sqrt(1- ((x1p(i)-x5)/(x5-x1))**2)


          else if (((x1p(i)) .gt. x6) .and. ((x1p(i)) .lt. x3)) then

              ycurve(i) = x3_c
              ycurve_down(i) =x3_s

          else

              ycurve(i)=L2
              ycurve_down(i)=0.

          endif
   
          if (((x1u(i)) .ge. (x5)) .and. ((x1u(i)) .le. (x6))) then

              yhcurve(i) = x_cc + sqrt(r_com**2-(x1u(i)-y_cc)**2)
              yhcurve_down(i) = x_cs - sqrt(r_sinus**2-(x1u(i)-y_cs)**2)

          else if (((x1u(i)) .ge. x2) .and. ((x1u(i)) .lt. x5)) then

!---- The straight inlet

               yhcurve(i) = x_ref + 2.*r_r
               yhcurve_down(i) = x_ref


!---- The oval shape inlet
!                yhcurve(i) = L2 - L2/3. * sqrt(1- ((x1u(i)-x5)/(x5-x1))**2)
!                yhcurve_down(i) = L2/3. * sqrt(1- ((x1u(i)-x5)/(x5-x1))**2)

          else if (((x1u(i)) .gt. x6) .and. ((x1u(i)) .lt. x3)) then

               yhcurve(i) = x3_c
               yhcurve_down(i) = x3_s

          else

               yhcurve(i)=L2
               yhcurve_down(i)=0.

          endif


       enddo

      end subroutine parallel_shape_function



      subroutine parallel_inner_body_point_detection(x,y,ll,flag, ycurve, yhcurve, ycurve_down, yhcurve_down)


        real(8), intent(in)                     ::    x , y
        integer, intent(in)                  ::    flag
        integer, intent(out)                 ::    ll
        integer                              ::    i


        real(8), intent(in), dimension(N1p)     ::    ycurve, yhcurve
        real(8), intent(in), dimension(N1p)     ::    ycurve_down, yhcurve_down


  !      print*, ycurve



        if (flag .eq. 2) then
            ll=0
            !print*, x
            do i=1,N1p
           !     print*, x1p(i), x
                if ((x .eq. x1p(i)))  then
                    ! print*, 'we are here'                   
                     if ((-ycurve(i)+(y)-dx1p(1)/2.0) .ge. 0.0) then
                          ll=3
                     endif
                endif
            enddo
            do i=1, N1p
                if ((x .eq. x1p(i)))  then
                     if ((-ycurve_down(i)+(y)-dx1p(1)/2.0) .le. 0.0) then
                          ll=3
                    endif
                endif
            enddo
        endif


        if (flag .eq. 3) then
            ll=0
            do i=1, N1p
                 if ((x .eq. x1u(i))) then
                       if ((-yhcurve(i)+(y)) .ge. 0.0) then
                             ll=3
                       endif
                 endif
            enddo
            do i=1, N1p
                 if ((x .eq. x1u(i))) then
                       if ((-yhcurve_down(i)+(y)) .le. 0.0) then
                             ll=3
                       endif
                 endif
            enddo
        endif


     end subroutine parallel_inner_body_point_detection




     subroutine parallel_save(InUsize, InVsize, Inner_V_i, Inner_V_j, Inner_U_i, Inner_U_j, ycurve, yhcurve, ycurve_down, yhcurve_down)
     !! TODO get the Inner_U,V_i,j to the timeint module variable header           

          integer, intent(out)                 ::   InVsize, InUsize
          integer                              ::   i, j
          integer                              ::   ll

          integer, intent(out)                 ::   Inner_V_i(200000), Inner_V_j(200000)
          integer, intent(out)                 ::   Inner_U_i(200000), Inner_U_j(200000)

          real(8), intent(in), dimension(N1p)     ::    ycurve, yhcurve
          real(8), intent(in), dimension(N1p)     ::    ycurve_down, yhcurve_down



          InVsize = 0
          InUsize = 0

          ll = 0
          do i = 1, N1p
               do j = 1, N2p
                   call parallel_inner_body_point_detection(x1p(i), x2v(j) , ll , 2, ycurve, yhcurve, ycurve_down, yhcurve_down)
                   if (ll .eq. 3) then
                       ll = 0
                       InVsize = InVsize+1
                       Inner_V_i(InVsize) = i
                       Inner_V_j(InVsize) = j
                   endif
                   call parallel_inner_body_point_detection(x1u(i), x2p(j) , ll , 3, ycurve, yhcurve, ycurve_down, yhcurve_down)
                   if (ll .eq. 3) then
                       ll = 0
                       InUsize = InUsize + 1
                       Inner_U_i(InUsize) = i
                       Inner_U_j(InUsize) = j
                   endif
             enddo
        enddo

 !       print*, InUsize, InVsize, rank

     end subroutine parallel_save


     subroutine parallel_apply_IB(InUsize, InVsize, Inner_V_i, Inner_V_j, Inner_U_i, Inner_U_j,u_in, v_in)


          integer, intent(in)                   ::   InVsize, InUsize
          integer, intent(in)                   ::   Inner_V_i(200000), Inner_V_j(200000)
          integer, intent(in)                   ::   Inner_U_i(200000), Inner_U_j(200000)


          real(8), intent(inout), dimension(N1p,N2p) ::   u_in,  v_in


          integer :: i
          !print*, rank, InUsize

          do i=1,InUsize
             u_in(Inner_U_i(i),Inner_U_j(i))=0.0
          enddo

          do i=1,InVsize
             v_in(Inner_V_i(i),Inner_V_j(i))=0.0
          enddo


     end subroutine parallel_apply_IB

      subroutine projection_vel_IBM_line_2D(u, v, xo , yo, xe, ye)


         integer             ::     ii,jj,kk



         integer             ::     x_index,  y_index

         integer             ::     n_ghost                           !   number of ghost cells


         real(8)                ::     ghost_cells(2, 10000)               !   ghost cells
         real(8)                ::     image_cells(2, 10000)               !   image cells
         real(8), intent(in)    ::     xo, yo, xe, ye                      !   origin point and end point of the line
         real(8)                ::     xt, yt, xg, yg, xin, yin, xip, yip  !   test, ghost,intercept and image point coordinates 
         real(8)                ::     m1, m2                              !   slope of line and the normal
         real(8)                ::     g_dist                              !   ghost cell distance from the line
         real(8)                ::     u(N1p,N2p),  v(N1p,N2p)
         real(8)                ::     x1 , y1, x2, y2, x3, y3, x4, y4
         integer             ::     Xindex, Yindex
         real(8)                ::     Bvec(4,1), C1, C2, C3, C4
         character*20        ::     file_name1

! calculating computational geometry which has to be done outside of the loop

         m1 = (ye - yo) / (xe -xo)
         m2 = -1. / m1


         n_ghost = 0
         do ii = S11B, N11B
            do jj = S21B, N21B
! on the u grid
                xt = x1u(ii)
                yt = x2p(jj)

                if ((xt .ge. xo)  .and.  (xt .le. xe)) then
                    g_dist = abs(yt - yo - m1*(xt - xo))/sqrt(m1**2 + 1.)
                    if (g_dist .le. 1.*dx1p(ii)) then
                        xg = xt
                        yg = yt
                        n_ghost = n_ghost + 1

                       !call Coordfinder(3,x1p,x2p,x1u,x2v,xg,yg,Xindex,Yindex,N1p,N2p)

                        u(ii,jj) = 0.
                    endif
                endif


! on the v grid

                xt = x1p(ii)
                yt = x2v(jj)
                if ((xt .ge. xo)  .and.  (xt .le. xe)) then
                    g_dist = abs(yt - yo - m1*(xt - xo))/sqrt(m1**2 + 1.)

                    if (g_dist .le. 1.*dx1p(ii)) then

                      !  n_ghost = n_ghost+1

                        xg = xt
                        yg = yt
                       ! call Coordfinder(2,x1p,x2p,x1u,x2v,xg,yg,Xindex,Yindex,N1p,N2p)

                        v(ii,jj) = 0.

                    endif
                endif


            enddo
         enddo

      end subroutine projection_vel_IBM_line_2D
      subroutine vel_IBM_vertical_line_2D(u, v, xo , yo, l, theta, angle_amplitude, freq, t)


         integer             ::     ii,jj,kk



         integer             ::     x_index,  y_index

         integer             ::     n_ghost                           !   number of ghost cells

         real(8), intent(in)    ::     theta, angle_amplitude, freq, t
         real(8)                ::     ghost_cells(2, 10000)               !   ghost cells
         real(8)                ::     image_cells(2, 10000)               !   image cells
         real(8), intent(in)    ::     xo, yo,l                            !   origin point and end point of the line
         real(8)                ::     xt, yt, xg, yg, xin, yin, xip, yip  !   test, ghost,intercept and image point coordinates 
         real(8)                ::     xe, ye
         real(8)                ::     m1, m2                              !   slope of line and the normal
         real(8)                ::     g_dist                              !   ghost cell distance from the line
         real(8)                ::     u(N1p,N2p),  v(N1p,N2p)
         real(8)                ::     x1 , y1, x2, y2, x3, y3, x4, y4
         integer             ::     Xindex, Yindex
         character*20        ::     file_name1

         real(8)                ::     radius

         real(8)                ::     phase, omega


         real(8),parameter      ::     pi = 3.14159265359 


! calculating computational geometry which has to be done outside of the loop

         phase = theta + angle_amplitude * sin (2. * pi * freq * t)

         xe = xo + l * cos(phase)
         ye = yo + l * sin(phase)


        ! m1 = (ye - yo) / (xe -xo)
        ! m2 = -1. / m1


         n_ghost = 0
         do ii = S11B, N11B
            do jj = S21B, N21B
! on the u grid
                xt = x1u(ii)
                yt = x2p(jj)

                if (((yt .ge. yo)  .and.  (yt .le. ye)) .or. ((yt .ge. ye)  .and.  (yt .le. yo)) ) then
                    g_dist = abs(xt-xo)
                    !g_dist = abs(yt - yo - m1*(xt - xo))/sqrt(m1**2 + 1.)
                    !if (((g_dist) .ge. -1.5*dx1p(ii)) .and.  (g_dist .le. 0.0)) then
                    if (g_dist .le. 0.5) then 
                        xg = xt
                        yg = yt
                        !n_ghost = n_ghost + 1

                        u(ii,jj) = 0.
                    endif

! the smooth ends
                     
                endif
!                    g_dist = sqrt((xt - xo)**2 + (yt-yo)**2) 
!                    if (g_dist .le. 0.005) then
!                        u(ii,jj) = 0.
!                    endif
 
!                    g_dist = sqrt((xt - xe)**2 + (yt-ye)**2)
!                    if (g_dist .le. 0.005) then
!                        u(ii,jj) = 0.
!                    endif


! on the v grid

                xt = x1p(ii)
                yt = x2v(jj)
                if (((yt .ge. yo)  .and.  (yt .le. ye)) .or.  ((yt .ge. ye)  .and.  (yt .le. yo)) )  then
                    g_dist = abs(xt-xo)
                    !g_dist = (yt - yo - m1*(xt - xo))/sqrt(m1**2 + 1.)

                    !if (((g_dist) .ge. -1.5*dx1p(ii)) .and. (g_dist .le. 0.)) then
                    if (g_dist .le. 0.5) then
                      !  n_ghost = n_ghost+1

                        xg = xt
                        yg = yt
                        radius = sqrt((xt-xo)**2 + (yt-yo)**2)

                        !call Coordfinder(3,x1p,x2p,x1u,x2v,xg,yg,Xindex,Yindex,N1p,N2p)

                     !   phase = theta + angle_amplitude * sin (2. * pi * freq * t)
                        omega =  2. * pi * freq * angle_amplitude * cos(2.*pi*freq*t)
                    
                        v(ii,jj) = 0.                         
 
                     !!   v(ii,jj) = sign(1., cos(2. * pi * freq * t)) * radius * omega  * cos(phase)
                        !v(ii,jj) = sign(1., cos(2. * pi * freq * t)) * radius * 2. * pi * freq * cos(phase)

!                        call Coordfinder(2,x1p,x2p,x1u,x2v,xg,yg,Xindex,Yindex,N1p,N2p)

!                        v(ii,jj) = 0.

                    endif
               endif
!                    g_dist = sqrt((xt - xo)**2 + (yt-yo)**2) 
!                    if (g_dist .le. 0.005) then
!                        v(ii,jj) = 0.
!                    endif
!                    g_dist = sqrt((xt - xe)**2 + (yt-ye)**2)
!                    if (g_dist .le. 0.005) then
!                        v(ii,jj) = 0.
!                    endif
             enddo
         enddo

      end subroutine vel_IBM_vertical_line_2D

      subroutine modify_Reynolds(i, j, xo , yo, l, theta)

         integer, intent(in) ::     i,j                          !   number of ghost cells
         
         real(8), intent(in)    ::     theta
         real(8), intent(in)    ::     xo, yo,l                            !   origin point and end point of the line
         real(8)                ::     xt, yt
         real(8)                ::     xe, ye
         real(8)                ::     g_dist                              !   ghost cell distance from the line
         character*20        ::     file_name1

         real(8)                ::     radius, r0

         real(8)                ::     phase, omega

         real(8),parameter      ::     pi = 3.14159265359 


         ! take only the horizontal profiles
         phase = 0. 
         xe = xo + l 
         ye = yo

         ! the gap radius exactly at the valve inlet
         r0 = 0.03 + (xo-0.839) * tan(theta)

                xt = x1p(i)
                yt = x2p(j)

                if ((xt .ge. xo)  .and.  (xt .le. xe)) then
                    radius = r0 + (xt - xo) * tan(theta)
                    g_dist = (yt - yo)
                    if ((g_dist .le. radius) .and. (g_dist .ge. 0.)) then 
                        Re = 5000;
                    endif
                endif
      end subroutine modify_Reynolds

      subroutine body_projection_vel_IBM_volume_3D(xc, yc, r)


         integer             ::     ii,jj,kk



         integer             ::     x_index,  y_index

         integer             ::     n_ghost                           !   number of ghost cells

         real(8), intent(in)    ::     xc,yc,r
         real(8)                ::     ghost_cells(2, 10000)                      !   ghost cells
         real(8)                ::     image_cells(2, 10000)                      !   image cells
!         real(8), intent(in)    ::     xo, yo,l                                   !   origin point and end point of the line
         real(8)                ::     xt, yt, zt, xg, yg, xin, yin, z_min, z_max  !   test, ghost,intercept and image point coordinates 
         real(8)                ::     arm,y_arm
         real(8)                ::     m1, m2                                     !   slope of line and the normal
         real(8)                ::     g_dist                                     !   ghost cell distance from the line
         real(8)                ::     u(N1p,N2p),  v(N1p,N2p)
         real(8)                ::     x1 , y1, x2, y2, x3, y3, x4, y4
         integer             ::     Xindex, Yindex
         character*20        ::     file_name1

         real(8)                ::     radius

         real(8)                ::     phase, omega


         real(8)                ::     p1(2), p2(2), p3(2)
         real(8)                ::     h1, h2, h3, p1p2, p1p3, p2p3

         real(8),parameter      ::     pi = 3.14159265359 
         real(8), parameter     ::     delta = 0.03 

! calculating computational geometry which has to be done outside of the loop

! on the u grid

      do kk =S31B, N31B
            do ii = S11B, N11B
                  do jj = S21B, N21B
                     if (.not.(sqrt((x2p(jj)-0.5)**2 + (x3p(kk)-0.5)**2) .le. 1./3.) .and. .not.(sqrt((x1u(ii) -xc)**2+(x2p(jj)-0.5)**2 + (x3p(kk)-0.5)**2) .le. r)) then 
                        vel(ii,jj,kk,1) = 0.
                     endif
                  enddo
            enddo
      enddo

! on the v grid
      do kk =S31B, N31B
            do ii = S11B, N11B
                  do jj = S21B, N21B
                     if (.not.(sqrt((x2v(jj)-0.5)**2 + (x3p(kk)-0.5)**2) .le. 1./3.) .and. .not.(sqrt((x1p(ii) -xc)**2+(x2v(jj)-0.5)**2 + (x3p(kk)-0.5)**2) .le. r)) then 
                        vel(ii,jj,kk,2) = 0.
                     endif
                  enddo
            enddo
      enddo


! on the w grid
      do kk =S31B, N31B
            do ii = S11B, N11B
                  do jj = S21B, N21B
                     if (.not.(sqrt((x2p(jj)-0.5)**2 + (x3w(kk)-0.5)**2) .le. 1./3.) .and. .not.(sqrt((x1p(ii) -xc)**2+(x2p(jj)-0.5)**2 + (x3w(kk)-0.5)**2) .le. r)) then 
                        vel(ii,jj,kk,3) = 0.
                     endif
                  enddo
            enddo
      enddo




     end subroutine body_projection_vel_IBM_volume_3D


      subroutine body_projection_vel_IBM_volume_3D_sinused(xc, yc, radi)


         integer             ::     ii,jj,kk



         integer             ::     x_index,  y_index

         integer             ::     n_ghost                           !   number of ghost cells

         real(8), intent(in)    ::     xc,yc,radi
         real(8)                ::     ghost_cells(2, 10000)                      !   ghost cells
         real(8)                ::     image_cells(2, 10000)                      !   image cells
!         real(8), intent(in)    ::     xo, yo,l                                   !   origin point and end point of the line
         real(8)                ::     xt, yt, zt, xg, yg, xin, yin, z_min, z_max  !   test, ghost,intercept and image point coordinates 
         real(8)                ::     arm,y_arm
         real(8)                ::     m1, m2                                     !   slope of line and the normal
         real(8)                ::     g_dist                                     !   ghost cell distance from the line
         real(8)                ::     u(N1p,N2p),  v(N1p,N2p)
         real(8)                ::     x1 , y1, x2, y2, x3, y3, x4, y4
         integer             ::     Xindex, Yindex
         character*20        ::     file_name1

         real(8)                ::     radius, r

         real(8)                ::     phase, omega


         real(8)                ::     p1(2), p2(2), p3(2)
         real(8)                ::     h1, h2, h3, p1p2, p1p3, p2p3

         real(8),parameter      ::     pi = 3.14159265359 
         real(8), parameter     ::     delta = 0.03 

! calculating computational geometry which has to be done outside of the loop

! on the u grid

      r=1./3.

      do kk =S31B, N31B
            do ii = S11B, N11B
                  do jj = S21B, N21B
                     if (.not.(sqrt((x2p(jj)-0.5)**2 + (x3p(kk)-0.5)**2) .le. 1./3.) .and. &
                         .not.(sqrt((x1u(ii) -xc)**2+(x3p(kk)-0.5-r/2.)**2 + (x2p(jj)-0.5)**2) .le. (r*sqrt(3.)/2.)) .and. &
                         .not.(sqrt((x1u(ii) -xc)**2+(x3p(kk)-0.5+r/4.)**2 + (x2p(jj)-0.5+r/2.*sqrt(3.)/2.)**2) .le. (r*sqrt(3.)/2.)) .and. &
                         .not.(sqrt((x1u(ii) -xc)**2+(x3p(kk)-0.5+r/4.)**2 + (x2p(jj)-0.5-r/2.*sqrt(3.)/2.)**2) .le. (r*sqrt(3.)/2.))  &
                           ) then 
                        vel(ii,jj,kk,1) = 0.
                     endif
                  enddo
            enddo
      enddo

! on the v grid
      do kk =S31B, N31B
            do ii = S11B, N11B
                  do jj = S21B, N21B
                      if (.not.(sqrt((x2v(jj)-0.5)**2 + (x3p(kk)-0.5)**2) .le. 1./3.) .and. &
                         .not.(sqrt((x1p(ii) -xc)**2+(x3p(kk)-0.5-r/2.)**2 + (x2v(jj)-0.5)**2) .le. (r*sqrt(3.)/2.)) .and. &
                         .not.(sqrt((x1p(ii) -xc)**2+(x3p(kk)-0.5+r/4.)**2 + (x2v(jj)-0.5+r/2.*sqrt(3.)/2.)**2) .le. (r*sqrt(3.)/2.)) .and. &
                         .not.(sqrt((x1p(ii) -xc)**2+(x3p(kk)-0.5+r/4.)**2 + (x2v(jj)-0.5-r/2.*sqrt(3.)/2.)**2) .le. (r*sqrt(3.)/2.))  &
                           ) then 
                      vel(ii,jj,kk,2) = 0.
                     endif
                  enddo
            enddo
      enddo


! on the w grid
      do kk =S31B, N31B
            do ii = S11B, N11B
                  do jj = S21B, N21B
                       if (.not.(sqrt((x2p(jj)-0.5)**2 + (x3w(kk)-0.5)**2) .le. 1./3.) .and. &
                         .not.(sqrt((x1p(ii) -xc)**2+(x3w(kk)-0.5-r/2.)**2 + (x2p(jj)-0.5)**2) .le. (r*sqrt(3.)/2.)) .and. &
                         .not.(sqrt((x1p(ii) -xc)**2+(x3w(kk)-0.5+r/4.)**2 + (x2p(jj)-0.5+r/2.*sqrt(3.)/2.)**2) .le. (r*sqrt(3.)/2.)) .and. &
                         .not.(sqrt((x1p(ii) -xc)**2+(x3w(kk)-0.5+r/4.)**2 + (x2p(jj)-0.5-r/2.*sqrt(3.)/2.)**2) .le. (r*sqrt(3.)/2.))   &
                           ) then 
                       vel(ii,jj,kk,3) = 0.
                     endif
                  enddo
            enddo
      enddo




     end subroutine body_projection_vel_IBM_volume_3D_sinused





      subroutine oscillator_projection_vel_IBM_volume_3D(xo , yo, l, theta)


         integer             ::     ii,jj,kk



         integer             ::     x_index,  y_index

         integer             ::     n_ghost                           !   number of ghost cells

         real(8), intent(in)    ::     theta
         real(8)                ::     ghost_cells(2, 10000)                      !   ghost cells
         real(8)                ::     image_cells(2, 10000)                      !   image cells
         real(8), intent(in)    ::     xo, yo,l                                   !   origin point and end point of the line
         real(8)                ::     xt, yt, zt, xg, yg, xin, yin, z_min, z_max  !   test, ghost,intercept and image point coordinates 
         real(8)                ::     arm,y_arm
         real(8)                ::     m1, m2                                     !   slope of line and the normal
         real(8)                ::     g_dist                                     !   ghost cell distance from the line
         real(8)                ::     u(N1p,N2p),  v(N1p,N2p)
         real(8)                ::     x1 , y1, x2, y2, x3, y3, x4, y4
         integer             ::     Xindex, Yindex
         character*20        ::     file_name1

         real(8)                ::     radius

         real(8)                ::     phase, omega


         real(8)                ::     p1(2), p2(2), p3(2)
         real(8)                ::     h1, h2, h3, p1p2, p1p3, p2p3

         real(8),parameter      ::     pi = 3.14159265359 
         real(8), parameter     ::     delta = 0.03 

! calculating computational geometry which has to be done outside of the loop

     z_min = 0.5 - 1./3. 
     z_max = 0.5 + 1./3.

! on the u grid
    !  if (rank ==0 ) print*, timestep
      do kk =S31B, N31B
         if ((x3p(kk) .le. z_max) .and. (x3p(kk) .ge. z_min)) then 
            do ii = S11B, N11B
                !get arm
                zt = abs(x3p(kk) - 0.5)  
                arm = sqrt(l**2 - zt**2)
                if ((x1u(ii) .ge. xo)  .and.  (x1u(ii) .le. (xo+arm))) then
                  !get  y_arm
                  y_arm = yo + tan(theta) * (x1u(ii) - xo)
                  do jj = S21B, N21B
                     if ((x2p(jj) .lt. (y_arm+delta/2.)) .and. (x2p(jj) .gt. (y_arm-delta/2.))) then 
                        vel(ii,jj,kk,1) = 0.
                      !  if (timestep .le. 2)   print*, 'ata'
                     endif
                  enddo
               endif
            enddo
         endif
      enddo

! on the v grid

      do kk =S31B, N31B
         if ((x3p(kk) .le. z_max) .and. (x3p(kk) .ge. z_min)) then
            do ii = S11B, N11B
                !get arm
                zt = abs(x3p(kk) - 0.5)  
                arm = sqrt(l**2 - zt**2)
                if ((x1p(ii) .ge. xo)  .and.  (x1p(ii) .le. (xo+arm))) then
                  !get  y_arm
                  y_arm = yo + tan(theta) * (x1p(ii) - xo)
                  do jj = S21B, N21B
                     if ((x2v(jj) .lt. (y_arm+delta/2.)) .and. (x2v(jj) .gt. (y_arm-delta/2.))) then 
                        vel(ii,jj,kk,2) = 0.
                     endif
                  enddo
               endif
            enddo
         endif
      enddo

! on the w grid

      do kk =S31B, N31B
         if ((x3w(kk) .le. z_max) .and. (x3w(kk) .ge. z_min)) then
            do ii = S11B, N11B
                !get arm
                zt = abs(x3w(kk) - 0.5)  
                arm = sqrt(l**2 - zt**2)
                if ((x1p(ii) .ge. xo)  .and.  (x1p(ii) .le. (xo+arm))) then
                  !get  y_arm
                  y_arm = yo + tan(theta) * (x1p(ii) - xo)
                  do jj = S21B, N21B
                     if ((x2p(jj) .lt. (y_arm+delta/2.)) .and. (x2p(jj) .gt. (y_arm-delta/2.))) then 
                        vel(ii,jj,kk,3) = 0.
                     endif
                  enddo
               endif
            enddo
         endif
      enddo




     end subroutine oscillator_projection_vel_IBM_volume_3D


      subroutine oscillator_projection_vel_IBM_volume_3D_edged(xo , yo, l, theta)


         integer             ::     ii,jj,kk



         integer             ::     x_index,  y_index

         integer             ::     n_ghost                           !   number of ghost cells

         real(8), intent(in)    ::     theta
         real(8)                ::     ghost_cells(2, 10000)                      !   ghost cells
         real(8)                ::     image_cells(2, 10000)                      !   image cells
         real(8), intent(in)    ::     xo, yo,l                                  !   origin point and end point of the line
         real(8)                ::     xt, yt, zt, xg, yg, xin, yin, z_min, z_max  !   test, ghost,intercept and image point coordinates 
         real(8)                ::     arm,y_arm, ye , xe
         real(8)                ::     m1, m2                                     !   slope of line and the normal
         real(8)                ::     g_dist                                     !   ghost cell distance from the line
         real(8)                ::     u(N1p,N2p),  v(N1p,N2p)
         real(8)                ::     x1 , y1, x2, y2, x3, y3, x4, y4
         integer             ::     Xindex, Yindex
         character*20        ::     file_name1

         real(8)                ::     radius

         real(8)                ::     phase, omega


         real(8)                ::     p1(2), p2(2), p3(2)
         real(8)                ::     h1, h2, h3, p1p2, p1p3, p2p3

         real(8),parameter      ::     pi = 3.14159265359 
         real(8), parameter     ::     delta = 0.03 

! calculating computational geometry which has to be done outside of the loop

     z_min = 0.5 - 1./3. 
     z_max = 0.5 + 1./3.

! on the u grid
    !  if (rank ==0 ) print*, timestep
      do kk =S31B, N31B
         if ((x3p(kk) .le. z_max) .and. (x3p(kk) .ge. z_min)) then 
            do ii = S11B, N11B
                !get arm
                zt = abs(x3p(kk) - 0.5)  
                arm = sqrt(l**2 - zt**2)
                if ((x1u(ii) .ge. xo)  .and.  (x1u(ii) .le. (xo+arm))) then
                  !get  y_arm
                  y_arm = yo + tan(theta) * (x1u(ii) - xo)
                  do jj = S21B, N21B
                     if ((x2p(jj) .lt. (y_arm+delta/2.)) .and. (x2p(jj) .gt. (y_arm-delta/2.))) then 
                        vel(ii,jj,kk,1) = 0.
                      !  if (timestep .le. 2)   print*, 'ata'
                     endif
                  enddo
               endif


     do jj = S21B,N21B
! the edges for the u grid

     if (theta .ge. 0.) then
!----triangle leading edge
                p1(1) = xo + delta/2. * sin(pi/36.)
                p1(2) = yo - delta/2. * cos(pi/36.)
                
                p2(1) = xo - delta/2. * sin(pi/36.)
                p2(2) = yo + delta/2. * cos(pi/36.)

                p3(1) = p2(1) - delta * 1./tan(pi/3.) * cos(pi/36.)
                p3(2) = p1(2) + delta / cos(pi/6.) * sin(pi/3. - pi/36.)


                xt = x1u(ii)
                yt = x2p(jj)  

                !-- p1 --p2
                h1 = abs(yt-p1(2)-tan(pi/2.+pi/36.)*(xt-p1(1)))/sqrt(1+(tan(pi/2.+pi/36.))**2 )
                !-- p1 --p3
                h2 = abs(yt-p1(2)-tan(pi + pi/36. - pi/3.)*(xt-p1(1)))/sqrt(1+(tan(pi + pi/36. - pi/3.))**2 )
                !-- p2 --p3
                h3 = abs(yt-p2(2)-tan(pi/36.)*(xt-p2(1)))/sqrt(1+(tan(pi/36.))**2)

                p1p2 = sqrt((p2(2)- p1(2))**2 + (p2(1)-p1(1))**2)
                p1p3 = sqrt((p3(2)- p1(2))**2 + (p3(1)-p1(1))**2)
                p2p3 = sqrt((p2(2)- p3(2))**2 + (p2(1)-p3(1))**2)

!                print*, h1, h2, h3


                if ((abs((h1*p1p2 + h2*p1p3 + h3*p2p3)- (delta*p2p3)) .le. 1e-12)) then
                  vel(ii,jj,kk,1) = 0.
                endif
!----triangle trailing edge
                xe = xo + arm * cos(theta)
                ye = yo + arm * sin(theta)




                p1(1) = xe - delta/2. * sin(pi/36.)
                p1(2) = ye + delta/2. * cos(pi/36.)

                p2(1) = xe + delta/2. * sin(pi/36.)
                p2(2) = ye - delta/2. * cos(pi/36.)

                p3(1) = p2(1) + delta * 1./tan(pi/3.) * cos(pi/36.)
                p3(2) = p1(2) - delta / cos(pi/6.) * sin(pi/3. - pi/36.)

                xt = x1u(ii)
                yt = x2p(jj)  

                !-- p1 --p2
                h1 = abs(yt-p1(2)-tan(pi/2.+pi/36.)*(xt-p1(1)))/sqrt(1+(tan(pi/2.+pi/36.))**2 )
                !-- p1 --p3
                h2 = abs(yt-p1(2)-tan(pi + pi/36. - pi/3.)*(xt-p1(1)))/sqrt(1+(tan(pi + pi/36. - pi/3.))**2 )
                !-- p2 --p3
                h3 = abs(yt-p2(2)-tan(pi/36.)*(xt-p2(1)))/sqrt(1+(tan(pi/36.))**2)

                p1p2 = sqrt((p2(2)- p1(2))**2 + (p2(1)-p1(1))**2)
                p1p3 = sqrt((p3(2)- p1(2))**2 + (p3(1)-p1(1))**2)
                p2p3 = sqrt((p2(2)- p3(2))**2 + (p2(1)-p3(1))**2)

!                print*, h1, h2, h3


                if ((abs((h1*p1p2 + h2*p1p3 + h3*p2p3)- (delta*p2p3)) .le. 1e-12)) then
                  vel(ii,jj,kk,1) = 0.
                endif
      endif



     if (theta .le. 0.) then
!----triangle leading edge
                p1(1) = xo + delta/2. * sin(pi/36.)
                p1(2) = yo + delta/2. * cos(pi/36.)

                p2(1) = xo - delta/2. * sin(pi/36.)
                p2(2) = yo - delta/2. * cos(pi/36.)

                p3(1) = p2(1) - delta * 1./tan(pi/3.) * cos(pi/36.)
                p3(2) = p1(2) - delta / cos(pi/6.) * sin(pi/3. - pi/36.)


                xt = x1u(ii)
                yt = x2p(jj)

                !-- p1 --p2
                h1 = abs(yt-p1(2)-tan(pi/2. - pi/36.)*(xt-p1(1)))/sqrt(1+(tan(pi/2. - pi/36.))**2 )
                !-- p1 --p3
                h2 = abs(yt-p1(2)-tan(pi/2. - pi/36. - pi/6.)*(xt-p1(1)))/sqrt(1+(tan(pi/2. - pi/36. - pi/6.))**2 )
                !-- p2 --p3
                h3 = abs(yt-p2(2)-tan(pi - pi/36.)*(xt-p2(1)))/sqrt(1+(tan(pi - pi/36.))**2)

                p1p2 = sqrt((p2(2)- p1(2))**2 + (p2(1)-p1(1))**2)
                p1p3 = sqrt((p3(2)- p1(2))**2 + (p3(1)-p1(1))**2)
                p2p3 = sqrt((p2(2)- p3(2))**2 + (p2(1)-p3(1))**2)

!                print*, h1, h2, h3


                if ((abs((h1*p1p2 + h2*p1p3 + h3*p2p3)- (delta*p2p3)) .le. 1e-12)) then
                  vel(ii,jj,kk,1) = 0.
                endif
!----triangle trailing edge
                xe = xo + arm * cos(theta)
                ye = yo + arm * sin(theta)


                p1(1) = xe - delta/2. * sin(pi/36.)
                p1(2) = ye - delta/2. * cos(pi/36.)

                p2(1) = xe + delta/2. * sin(pi/36.)
                p2(2) = ye + delta/2. * cos(pi/36.)

                p3(1) = p2(1) + delta * 1./tan(pi/3.) * cos(pi/36.)
                p3(2) = p1(2) + delta / cos(pi/6.) * sin(pi/3. - pi/36.)

!                xt = xo
!                yt = yo  

                !-- p1 --p2
                h1 = abs(yt-p1(2)-tan(pi/2. - pi/36.)*(xt-p1(1)))/sqrt(1+(tan(pi/2. - pi/36.))**2 )
                !-- p1 --p3
                h2 = abs(yt-p1(2)-tan(pi/2. - pi/36. - pi/6.)*(xt-p1(1)))/sqrt(1+(tan(pi/2. - pi/36. - pi/6.))**2 )
                !-- p2 --p3
                h3 = abs(yt-p2(2)-tan(pi - pi/36.)*(xt-p2(1)))/sqrt(1+(tan(pi - pi/36.))**2)

                p1p2 = sqrt((p2(2)- p1(2))**2 + (p2(1)-p1(1))**2)
                p1p3 = sqrt((p3(2)- p1(2))**2 + (p3(1)-p1(1))**2)
                p2p3 = sqrt((p2(2)- p3(2))**2 + (p2(1)-p3(1))**2)

!                print*, h1, h2, h3


                if ((abs((h1*p1p2 + h2*p1p3 + h3*p2p3)- (delta*p2p3)) .le. 1e-12)) then
                  vel(ii,jj,kk,1) = 0.
                endif
      endif


                 enddo
               enddo
         endif
     
  enddo


! on the v grid

      do kk =S31B, N31B
         if ((x3p(kk) .le. z_max) .and. (x3p(kk) .ge. z_min)) then
            do ii = S11B, N11B
                !get arm
                zt = abs(x3p(kk) - 0.5)  
                arm = sqrt(l**2 - zt**2)
                if ((x1p(ii) .ge. xo)  .and.  (x1p(ii) .le. (xo+arm))) then
                  !get  y_arm
                  y_arm = yo + tan(theta) * (x1p(ii) - xo)
                  do jj = S21B, N21B
                     if ((x2v(jj) .lt. (y_arm+delta/2.)) .and. (x2v(jj) .gt. (y_arm-delta/2.))) then 
                        vel(ii,jj,kk,2) = 0.
                     endif
                  enddo
               endif



  do jj = S21B, N21B


!  the edges for the v grid
     if (theta .ge. 0.) then
!----triangle leading edge
                p1(1) = xo + delta/2. * sin(pi/36.)
                p1(2) = yo - delta/2. * cos(pi/36.)
                
                p2(1) = xo - delta/2. * sin(pi/36.)
                p2(2) = yo + delta/2. * cos(pi/36.)

                p3(1) = p2(1) - delta * 1./tan(pi/3.) * cos(pi/36.)
                p3(2) = p1(2) + delta / cos(pi/6.) * sin(pi/3. - pi/36.)


                xt = x1p(ii)
                yt = x2v(jj)  

                !-- p1 --p2
                h1 = abs(yt-p1(2)-tan(pi/2.+pi/36.)*(xt-p1(1)))/sqrt(1+(tan(pi/2.+pi/36.))**2 )
                !-- p1 --p3
                h2 = abs(yt-p1(2)-tan(pi + pi/36. - pi/3.)*(xt-p1(1)))/sqrt(1+(tan(pi + pi/36. - pi/3.))**2 )
                !-- p2 --p3
                h3 = abs(yt-p2(2)-tan(pi/36.)*(xt-p2(1)))/sqrt(1+(tan(pi/36.))**2)

                p1p2 = sqrt((p2(2)- p1(2))**2 + (p2(1)-p1(1))**2)
                p1p3 = sqrt((p3(2)- p1(2))**2 + (p3(1)-p1(1))**2)
                p2p3 = sqrt((p2(2)- p3(2))**2 + (p2(1)-p3(1))**2)

!                print*, h1, h2, h3


                if ((abs((h1*p1p2 + h2*p1p3 + h3*p2p3)- (delta*p2p3)) .le. 1e-12)) then
                  vel(ii,jj,kk,2) = 0.
                endif
!----triangle trailing edge
                xe = xo + arm * cos(theta)
                ye = yo + arm * sin(theta)


                p1(1) = xe - delta/2. * sin(pi/36.)
                p1(2) = ye + delta/2. * cos(pi/36.)

                p2(1) = xe + delta/2. * sin(pi/36.)
                p2(2) = ye - delta/2. * cos(pi/36.)

                p3(1) = p2(1) + delta * 1./tan(pi/3.) * cos(pi/36.)
                p3(2) = p1(2) - delta / cos(pi/6.) * sin(pi/3. - pi/36.)

!                xt = xo
!                yt = yo  

                !-- p1 --p2
                h1 = abs(yt-p1(2)-tan(pi/2.+pi/36.)*(xt-p1(1)))/sqrt(1+(tan(pi/2.+pi/36.))**2 )
                !-- p1 --p3
                h2 = abs(yt-p1(2)-tan(pi + pi/36. - pi/3.)*(xt-p1(1)))/sqrt(1+(tan(pi + pi/36. - pi/3.))**2 )
                !-- p2 --p3
                h3 = abs(yt-p2(2)-tan(pi/36.)*(xt-p2(1)))/sqrt(1+(tan(pi/36.))**2)

                p1p2 = sqrt((p2(2)- p1(2))**2 + (p2(1)-p1(1))**2)
                p1p3 = sqrt((p3(2)- p1(2))**2 + (p3(1)-p1(1))**2)
                p2p3 = sqrt((p2(2)- p3(2))**2 + (p2(1)-p3(1))**2)

!                print*, h1, h2, h3


                if ((abs((h1*p1p2 + h2*p1p3 + h3*p2p3)- (delta*p2p3)) .le. 1e-12)) then
                  vel(ii,jj,kk,2) = 0.
                endif
      endif



     if (theta .le. 0.) then
!----triangle leading edge
                p1(1) = xo + delta/2. * sin(pi/36.)
                p1(2) = yo + delta/2. * cos(pi/36.)

                p2(1) = xo - delta/2. * sin(pi/36.)
                p2(2) = yo - delta/2. * cos(pi/36.)

                p3(1) = p2(1) - delta * 1./tan(pi/3.) * cos(pi/36.)
                p3(2) = p1(2) - delta / cos(pi/6.) * sin(pi/3. - pi/36.)


                xt = x1p(ii)
                yt = x2v(jj)  

                !-- p1 --p2
                h1 = abs(yt-p1(2)-tan(pi/2. - pi/36.)*(xt-p1(1)))/sqrt(1+(tan(pi/2. - pi/36.))**2 )
                !-- p1 --p3
                h2 = abs(yt-p1(2)-tan(pi/2. - pi/36. - pi/6.)*(xt-p1(1)))/sqrt(1+(tan(pi/2. - pi/36. - pi/6.))**2 )
                !-- p2 --p3
                h3 = abs(yt-p2(2)-tan(pi - pi/36.)*(xt-p2(1)))/sqrt(1+(tan(pi - pi/36.))**2)

                p1p2 = sqrt((p2(2)- p1(2))**2 + (p2(1)-p1(1))**2)
                p1p3 = sqrt((p3(2)- p1(2))**2 + (p3(1)-p1(1))**2)
                p2p3 = sqrt((p2(2)- p3(2))**2 + (p2(1)-p3(1))**2)

!                print*, h1, h2, h3


                if ((abs((h1*p1p2 + h2*p1p3 + h3*p2p3)- (delta*p2p3)) .le. 1e-12)) then
                  vel(ii,jj,kk,2) = 0.
                endif
!----triangle trailing edge

                xe = xo + arm * cos(theta)
                ye = yo + arm * sin(theta)


                p1(1) = xe - delta/2. * sin(pi/36.)
                p1(2) = ye - delta/2. * cos(pi/36.)

                p2(1) = xe + delta/2. * sin(pi/36.)
                p2(2) = ye + delta/2. * cos(pi/36.)

                p3(1) = p2(1) + delta * 1./tan(pi/3.) * cos(pi/36.)
                p3(2) = p1(2) + delta / cos(pi/6.) * sin(pi/3. - pi/36.)

                xt = x1p(ii)
                yt = x2v(jj)  

                !-- p1 --p2
                h1 = abs(yt-p1(2)-tan(pi/2. - pi/36.)*(xt-p1(1)))/sqrt(1+(tan(pi/2. - pi/36.))**2 )
                !-- p1 --p3
                h2 = abs(yt-p1(2)-tan(pi/2. - pi/36. - pi/6.)*(xt-p1(1)))/sqrt(1+(tan(pi/2. - pi/36. - pi/6.))**2 )
                !-- p2 --p3
                h3 = abs(yt-p2(2)-tan(pi - pi/36.)*(xt-p2(1)))/sqrt(1+(tan(pi - pi/36.))**2)

                p1p2 = sqrt((p2(2)- p1(2))**2 + (p2(1)-p1(1))**2)
                p1p3 = sqrt((p3(2)- p1(2))**2 + (p3(1)-p1(1))**2)
                p2p3 = sqrt((p2(2)- p3(2))**2 + (p2(1)-p3(1))**2)

!                print*, h1, h2, h3


                if ((abs((h1*p1p2 + h2*p1p3 + h3*p2p3)- (delta*p2p3)) .le. 1e-12)) then
                  vel(ii,jj,kk,2) = 0.
                endif
      endif

            enddo
           enddo
         endif
      enddo

! on the w grid

      do kk =S31B, N31B
         if ((x3w(kk) .le. z_max) .and. (x3w(kk) .ge. z_min)) then
            do ii = S11B, N11B
                !get arm
                zt = abs(x3w(kk) - 0.5)  
                arm = sqrt(l**2 - zt**2)
                if ((x1p(ii) .ge. xo)  .and.  (x1p(ii) .le. (xo+arm))) then
                  !get  y_arm
                  y_arm = yo + tan(theta) * (x1p(ii) - xo)
                  do jj = S21B, N21B
                     if ((x2p(jj) .lt. (y_arm+delta/2.)) .and. (x2p(jj) .gt. (y_arm-delta/2.))) then 
                        vel(ii,jj,kk,3) = 0.
                     endif
                  enddo
               endif

! the edges for the w grid

     do jj = S21B, N21B


!  the edges for the v grid
     if (theta .ge. 0.) then
!----triangle leading edge
                p1(1) = xo + delta/2. * sin(pi/36.)
                p1(2) = yo - delta/2. * cos(pi/36.)
                
                p2(1) = xo - delta/2. * sin(pi/36.)
                p2(2) = yo + delta/2. * cos(pi/36.)

                p3(1) = p2(1) - delta * 1./tan(pi/3.) * cos(pi/36.)
                p3(2) = p1(2) + delta / cos(pi/6.) * sin(pi/3. - pi/36.)


                xt = x1p(ii)
                yt = x2p(jj)  

                !-- p1 --p2
                h1 = abs(yt-p1(2)-tan(pi/2.+pi/36.)*(xt-p1(1)))/sqrt(1+(tan(pi/2.+pi/36.))**2 )
                !-- p1 --p3
                h2 = abs(yt-p1(2)-tan(pi + pi/36. - pi/3.)*(xt-p1(1)))/sqrt(1+(tan(pi + pi/36. - pi/3.))**2 )
                !-- p2 --p3
                h3 = abs(yt-p2(2)-tan(pi/36.)*(xt-p2(1)))/sqrt(1+(tan(pi/36.))**2)

                p1p2 = sqrt((p2(2)- p1(2))**2 + (p2(1)-p1(1))**2)
                p1p3 = sqrt((p3(2)- p1(2))**2 + (p3(1)-p1(1))**2)
                p2p3 = sqrt((p2(2)- p3(2))**2 + (p2(1)-p3(1))**2)

!                print*, h1, h2, h3


                if ((abs((h1*p1p2 + h2*p1p3 + h3*p2p3)- (delta*p2p3)) .le. 1e-12)) then
                  vel(ii,jj,kk,3) = 0.
                endif
!----triangle trailing edge
                xe = xo + arm * cos(theta)
                ye = yo + arm * sin(theta)


                p1(1) = xe - delta/2. * sin(pi/36.)
                p1(2) = ye + delta/2. * cos(pi/36.)

                p2(1) = xe + delta/2. * sin(pi/36.)
                p2(2) = ye - delta/2. * cos(pi/36.)

                p3(1) = p2(1) + delta * 1./tan(pi/3.) * cos(pi/36.)
                p3(2) = p1(2) - delta / cos(pi/6.) * sin(pi/3. - pi/36.)

!                xt = xo
!                yt = yo  

                !-- p1 --p2
                h1 = abs(yt-p1(2)-tan(pi/2.+pi/36.)*(xt-p1(1)))/sqrt(1+(tan(pi/2.+pi/36.))**2 )
                !-- p1 --p3
                h2 = abs(yt-p1(2)-tan(pi + pi/36. - pi/3.)*(xt-p1(1)))/sqrt(1+(tan(pi + pi/36. - pi/3.))**2 )
                !-- p2 --p3
                h3 = abs(yt-p2(2)-tan(pi/36.)*(xt-p2(1)))/sqrt(1+(tan(pi/36.))**2)

                p1p2 = sqrt((p2(2)- p1(2))**2 + (p2(1)-p1(1))**2)
                p1p3 = sqrt((p3(2)- p1(2))**2 + (p3(1)-p1(1))**2)
                p2p3 = sqrt((p2(2)- p3(2))**2 + (p2(1)-p3(1))**2)

!                print*, h1, h2, h3


                if ((abs((h1*p1p2 + h2*p1p3 + h3*p2p3)- (delta*p2p3)) .le. 1e-12)) then
                  vel(ii,jj,kk,3) = 0.
                endif
      endif



     if (theta .le. 0.) then
!----triangle leading edge
                p1(1) = xo + delta/2. * sin(pi/36.)
                p1(2) = yo + delta/2. * cos(pi/36.)

                p2(1) = xo - delta/2. * sin(pi/36.)
                p2(2) = yo - delta/2. * cos(pi/36.)

                p3(1) = p2(1) - delta * 1./tan(pi/3.) * cos(pi/36.)
                p3(2) = p1(2) - delta / cos(pi/6.) * sin(pi/3. - pi/36.)


                xt = x1p(ii)
                yt = x2p(jj)  

                !-- p1 --p2
                h1 = abs(yt-p1(2)-tan(pi/2. - pi/36.)*(xt-p1(1)))/sqrt(1+(tan(pi/2. - pi/36.))**2 )
                !-- p1 --p3
                h2 = abs(yt-p1(2)-tan(pi/2. - pi/36. - pi/6.)*(xt-p1(1)))/sqrt(1+(tan(pi/2. - pi/36. - pi/6.))**2 )
                !-- p2 --p3
                h3 = abs(yt-p2(2)-tan(pi - pi/36.)*(xt-p2(1)))/sqrt(1+(tan(pi - pi/36.))**2)

                p1p2 = sqrt((p2(2)- p1(2))**2 + (p2(1)-p1(1))**2)
                p1p3 = sqrt((p3(2)- p1(2))**2 + (p3(1)-p1(1))**2)
                p2p3 = sqrt((p2(2)- p3(2))**2 + (p2(1)-p3(1))**2)

!                print*, h1, h2, h3


                if ((abs((h1*p1p2 + h2*p1p3 + h3*p2p3)- (delta*p2p3)) .le. 1e-12)) then
                  vel(ii,jj,kk,3) = 0.
                endif
!----triangle trailing edge

                xe = xo + arm * cos(theta)
                ye = yo + arm * sin(theta)


                p1(1) = xe - delta/2. * sin(pi/36.)
                p1(2) = ye - delta/2. * cos(pi/36.)

                p2(1) = xe + delta/2. * sin(pi/36.)
                p2(2) = ye + delta/2. * cos(pi/36.)

                p3(1) = p2(1) + delta * 1./tan(pi/3.) * cos(pi/36.)
                p3(2) = p1(2) + delta / cos(pi/6.) * sin(pi/3. - pi/36.)

                xt = x1p(ii)
                yt = x2p(jj)  

                !-- p1 --p2
                h1 = abs(yt-p1(2)-tan(pi/2. - pi/36.)*(xt-p1(1)))/sqrt(1+(tan(pi/2. - pi/36.))**2 )
                !-- p1 --p3
                h2 = abs(yt-p1(2)-tan(pi/2. - pi/36. - pi/6.)*(xt-p1(1)))/sqrt(1+(tan(pi/2. - pi/36. - pi/6.))**2 )
                !-- p2 --p3
                h3 = abs(yt-p2(2)-tan(pi - pi/36.)*(xt-p2(1)))/sqrt(1+(tan(pi - pi/36.))**2)

                p1p2 = sqrt((p2(2)- p1(2))**2 + (p2(1)-p1(1))**2)
                p1p3 = sqrt((p3(2)- p1(2))**2 + (p3(1)-p1(1))**2)
                p2p3 = sqrt((p2(2)- p3(2))**2 + (p2(1)-p3(1))**2)

!                print*, h1, h2, h3


                if ((abs((h1*p1p2 + h2*p1p3 + h3*p2p3)- (delta*p2p3)) .le. 1e-12)) then
                  vel(ii,jj,kk,3) = 0.
                endif
      endif

            enddo
 



    enddo
         endif
      enddo




     end subroutine oscillator_projection_vel_IBM_volume_3D_edged






      subroutine oscillator_projection_vel_IBM_line_2D(u, v, xo , yo, l, theta, angle_amplitude, freq, t)


         integer             ::     ii,jj,kk



         integer             ::     x_index,  y_index

         integer             ::     n_ghost                           !   number of ghost cells

         real(8), intent(in)    ::     theta, angle_amplitude, freq, t
         real(8)                ::     ghost_cells(2, 10000)               !   ghost cells
         real(8)                ::     image_cells(2, 10000)               !   image cells
         real(8), intent(in)    ::     xo, yo,l                            !   origin point and end point of the line
         real(8)                ::     xt, yt, xg, yg, xin, yin, xip, yip  !   test, ghost,intercept and image point coordinates 
         real(8)                ::     xe, ye
         real(8)                ::     m1, m2                              !   slope of line and the normal
         real(8)                ::     g_dist                              !   ghost cell distance from the line
         real(8)                ::     u(N1p,N2p),  v(N1p,N2p)
         real(8)                ::     x1 , y1, x2, y2, x3, y3, x4, y4
         integer             ::     Xindex, Yindex
         character*20        ::     file_name1

         real(8)                ::     radius

         real(8)                ::     phase, omega


         real(8)                ::     p1(2), p2(2), p3(2)
         real(8)                ::     h1, h2, h3, p1p2, p1p3, p2p3

         real(8),parameter      ::     pi = 3.14159265359 
         real(8), parameter     ::     delta = 0.03 

! calculating computational geometry which has to be done outside of the loop

         phase = theta + angle_amplitude * sin (2. * pi * freq * t)

         xe = xo + l * cos(phase)
         ye = yo + l * sin(phase)


         m1 = (ye - yo) / (xe -xo)
         m2 = -1. / m1


         n_ghost = 0
         do ii = S11B, N11B
            do jj = S21B, N21B
! on the u grid
                xt = x1u(ii)
                yt = x2p(jj)

                if ((xt .ge. xo)  .and.  (xt .le. xe)) then
                    g_dist = abs(yt - yo - m1*(xt - xo))/sqrt(m1**2 + 1.)
                    !g_dist = abs(yt - yo - m1*(xt - xo))/sqrt(m1**2 + 1.)
                    !if (((g_dist) .ge. -1.5*dx1p(ii)) .and.  (g_dist .le. 0.0)) then
                    if (g_dist .le. delta/2) then 
                        xg = xt
                        yg = yt
                        !n_ghost = n_ghost + 1

                        radius = sqrt((xt-xo)**2 + (yt-yo)**2)

                        !call Coordfinder(3,x1p,x2p,x1u,x2v,xg,yg,Xindex,Yindex,N1p,N2p)
                        omega =  2. * pi * freq * angle_amplitude * cos(2.*pi*freq*t)

                      !  u(ii,jj) = -sign(1. , cos(2. * pi * freq * t)) * radius * omega *  sin(phase)  
                  
                        u(ii,jj) = 0.
                    endif

! the smooth ends
                     

                endif

!----round leading edge
!                    g_dist = sqrt((xt - xo)**2 + (yt-yo)**2) 
!                    if (g_dist .le. 0.0075) then
!                        u(ii,jj) = 0.
!                    endif
! 
!                    g_dist = sqrt((xt - xe)**2 + (yt-ye)**2)
!                    if (g_dist .le. 0.0075) then
!                        u(ii,jj) = 0.
!                    endif

     if (phase .ge. 0.) then
!----triangle leading edge
                p1(1) = xo + delta/2. * sin(pi/36.)
                p1(2) = yo - delta/2. * cos(pi/36.)
                
                p2(1) = xo - delta/2. * sin(pi/36.)
                p2(2) = yo + delta/2. * cos(pi/36.)

                p3(1) = p2(1) - delta * 1./tan(pi/3.) * cos(pi/36.)
                p3(2) = p1(2) + delta / cos(pi/6.) * sin(pi/3. - pi/36.)


!                xt = xo
!                yt = yo  

                !-- p1 --p2
                h1 = abs(yt-p1(2)-tan(pi/2.+pi/36.)*(xt-p1(1)))/sqrt(1+(tan(pi/2.+pi/36.))**2 )
                !-- p1 --p3
                h2 = abs(yt-p1(2)-tan(pi + pi/36. - pi/3.)*(xt-p1(1)))/sqrt(1+(tan(pi + pi/36. - pi/3.))**2 )
                !-- p2 --p3
                h3 = abs(yt-p2(2)-tan(pi/36.)*(xt-p2(1)))/sqrt(1+(tan(pi/36.))**2)

                p1p2 = sqrt((p2(2)- p1(2))**2 + (p2(1)-p1(1))**2)
                p1p3 = sqrt((p3(2)- p1(2))**2 + (p3(1)-p1(1))**2)
                p2p3 = sqrt((p2(2)- p3(2))**2 + (p2(1)-p3(1))**2)

!                print*, h1, h2, h3


                if ((abs((h1*p1p2 + h2*p1p3 + h3*p2p3)- (delta*p2p3)) .le. 1e-12)) then
                  u(ii,jj) = 0.
                endif
!----triangle trailing edge
                p1(1) = xe - delta/2. * sin(pi/36.)
                p1(2) = ye + delta/2. * cos(pi/36.)

                p2(1) = xe + delta/2. * sin(pi/36.)
                p2(2) = ye - delta/2. * cos(pi/36.)

                p3(1) = p2(1) + delta * 1./tan(pi/3.) * cos(pi/36.)
                p3(2) = p1(2) - delta / cos(pi/6.) * sin(pi/3. - pi/36.)

!                xt = xo
!                yt = yo  

                !-- p1 --p2
                h1 = abs(yt-p1(2)-tan(pi/2.+pi/36.)*(xt-p1(1)))/sqrt(1+(tan(pi/2.+pi/36.))**2 )
                !-- p1 --p3
                h2 = abs(yt-p1(2)-tan(pi + pi/36. - pi/3.)*(xt-p1(1)))/sqrt(1+(tan(pi + pi/36. - pi/3.))**2 )
                !-- p2 --p3
                h3 = abs(yt-p2(2)-tan(pi/36.)*(xt-p2(1)))/sqrt(1+(tan(pi/36.))**2)

                p1p2 = sqrt((p2(2)- p1(2))**2 + (p2(1)-p1(1))**2)
                p1p3 = sqrt((p3(2)- p1(2))**2 + (p3(1)-p1(1))**2)
                p2p3 = sqrt((p2(2)- p3(2))**2 + (p2(1)-p3(1))**2)

!                print*, h1, h2, h3


                if ((abs((h1*p1p2 + h2*p1p3 + h3*p2p3)- (delta*p2p3)) .le. 1e-12)) then
                  u(ii,jj) = 0.
                endif
      endif



     if (phase .le. 0.) then
!----triangle leading edge
                p1(1) = xo + delta/2. * sin(pi/36.)
                p1(2) = yo + delta/2. * cos(pi/36.)

                p2(1) = xo - delta/2. * sin(pi/36.)
                p2(2) = yo - delta/2. * cos(pi/36.)

                p3(1) = p2(1) - delta * 1./tan(pi/3.) * cos(pi/36.)
                p3(2) = p1(2) - delta / cos(pi/6.) * sin(pi/3. - pi/36.)


!                xt = xo
!                yt = yo  

                !-- p1 --p2
                h1 = abs(yt-p1(2)-tan(pi/2. - pi/36.)*(xt-p1(1)))/sqrt(1+(tan(pi/2. - pi/36.))**2 )
                !-- p1 --p3
                h2 = abs(yt-p1(2)-tan(pi/2. - pi/36. - pi/6.)*(xt-p1(1)))/sqrt(1+(tan(pi/2. - pi/36. - pi/6.))**2 )
                !-- p2 --p3
                h3 = abs(yt-p2(2)-tan(pi - pi/36.)*(xt-p2(1)))/sqrt(1+(tan(pi - pi/36.))**2)

                p1p2 = sqrt((p2(2)- p1(2))**2 + (p2(1)-p1(1))**2)
                p1p3 = sqrt((p3(2)- p1(2))**2 + (p3(1)-p1(1))**2)
                p2p3 = sqrt((p2(2)- p3(2))**2 + (p2(1)-p3(1))**2)

!                print*, h1, h2, h3


                if ((abs((h1*p1p2 + h2*p1p3 + h3*p2p3)- (delta*p2p3)) .le. 1e-12)) then
                  u(ii,jj) = 0.
                endif
!----triangle trailing edge
                p1(1) = xe - delta/2. * sin(pi/36.)
                p1(2) = ye - delta/2. * cos(pi/36.)

                p2(1) = xe + delta/2. * sin(pi/36.)
                p2(2) = ye + delta/2. * cos(pi/36.)

                p3(1) = p2(1) + delta * 1./tan(pi/3.) * cos(pi/36.)
                p3(2) = p1(2) + delta / cos(pi/6.) * sin(pi/3. - pi/36.)

!                xt = xo
!                yt = yo  

                !-- p1 --p2
                h1 = abs(yt-p1(2)-tan(pi/2. - pi/36.)*(xt-p1(1)))/sqrt(1+(tan(pi/2. - pi/36.))**2 )
                !-- p1 --p3
                h2 = abs(yt-p1(2)-tan(pi/2. - pi/36. - pi/6.)*(xt-p1(1)))/sqrt(1+(tan(pi/2. - pi/36. - pi/6.))**2 )
                !-- p2 --p3
                h3 = abs(yt-p2(2)-tan(pi - pi/36.)*(xt-p2(1)))/sqrt(1+(tan(pi - pi/36.))**2)

                p1p2 = sqrt((p2(2)- p1(2))**2 + (p2(1)-p1(1))**2)
                p1p3 = sqrt((p3(2)- p1(2))**2 + (p3(1)-p1(1))**2)
                p2p3 = sqrt((p2(2)- p3(2))**2 + (p2(1)-p3(1))**2)

!                print*, h1, h2, h3


                if ((abs((h1*p1p2 + h2*p1p3 + h3*p2p3)- (delta*p2p3)) .le. 1e-12)) then
                  u(ii,jj) = 0.
                endif
      endif



! on the v grid

                xt = x1p(ii)
                yt = x2v(jj)
                if ((xt .ge. xo)  .and.  (xt .le. xe)) then
                    g_dist = abs(yt - yo - m1*(xt - xo))/sqrt(m1**2 + 1.)
                    !g_dist = (yt - yo - m1*(xt - xo))/sqrt(m1**2 + 1.)

                    !if (((g_dist) .ge. -1.5*dx1p(ii)) .and. (g_dist .le. 0.)) then
                    if (g_dist .le. delta/2) then
                      !  n_ghost = n_ghost+1

                        xg = xt
                        yg = yt
                        radius = sqrt((xt-xo)**2 + (yt-yo)**2)

                        !call Coordfinder(3,x1p,x2p,x1u,x2v,xg,yg,Xindex,Yindex,N1p,N2p)

                     !   phase = theta + angle_amplitude * sin (2. * pi * freq * t)
                        omega =  2. * pi * freq * angle_amplitude * cos(2.*pi*freq*t)
                    
                        v(ii,jj) = 0.                         
 
                     !!   v(ii,jj) = sign(1., cos(2. * pi * freq * t)) * radius * omega  * cos(phase)
                        !v(ii,jj) = sign(1., cos(2. * pi * freq * t)) * radius * 2. * pi * freq * cos(phase)

!                        call Coordfinder(2,x1p,x2p,x1u,x2v,xg,yg,Xindex,Yindex,N1p,N2p)

!                        v(ii,jj) = 0.

                    endif
               endif
!----round leading edge
!                    g_dist = sqrt((xt - xo)**2 + (yt-yo)**2) 
!                    if (g_dist .le. 0.015) then
!                        v(ii,jj) = 0.
!                    endif
!                    g_dist = sqrt((xt - xe)**2 + (yt-ye)**2)
!                    if (g_dist .le. 0.015) then
!                        v(ii,jj) = 0.
!                    endif
!----triangle leading edge

     if (phase .ge. 0.) then
!----triangle leading edge
                p1(1) = xo + delta/2. * sin(pi/36.)
                p1(2) = yo - delta/2. * cos(pi/36.)
                
                p2(1) = xo - delta/2. * sin(pi/36.)
                p2(2) = yo + delta/2. * cos(pi/36.)

                p3(1) = p2(1) - delta * 1./tan(pi/3.) * cos(pi/36.)
                p3(2) = p1(2) + delta / cos(pi/6.) * sin(pi/3. - pi/36.)


!                xt = xo
!                yt = yo  

                !-- p1 --p2
                h1 = abs(yt-p1(2)-tan(pi/2.+pi/36.)*(xt-p1(1)))/sqrt(1+(tan(pi/2.+pi/36.))**2 )
                !-- p1 --p3
                h2 = abs(yt-p1(2)-tan(pi + pi/36. - pi/3.)*(xt-p1(1)))/sqrt(1+(tan(pi + pi/36. - pi/3.))**2 )
                !-- p2 --p3
                h3 = abs(yt-p2(2)-tan(pi/36.)*(xt-p2(1)))/sqrt(1+(tan(pi/36.))**2)

                p1p2 = sqrt((p2(2)- p1(2))**2 + (p2(1)-p1(1))**2)
                p1p3 = sqrt((p3(2)- p1(2))**2 + (p3(1)-p1(1))**2)
                p2p3 = sqrt((p2(2)- p3(2))**2 + (p2(1)-p3(1))**2)

!                print*, h1, h2, h3


                if ((abs((h1*p1p2 + h2*p1p3 + h3*p2p3)- (delta*p2p3)) .le. 1e-12)) then
                  v(ii,jj) = 0.
                endif
!----triangle trailing edge
                p1(1) = xe - delta/2. * sin(pi/36.)
                p1(2) = ye + delta/2. * cos(pi/36.)

                p2(1) = xe + delta/2. * sin(pi/36.)
                p2(2) = ye - delta/2. * cos(pi/36.)

                p3(1) = p2(1) + delta * 1./tan(pi/3.) * cos(pi/36.)
                p3(2) = p1(2) - delta / cos(pi/6.) * sin(pi/3. - pi/36.)

!                xt = xo
!                yt = yo  

                !-- p1 --p2
                h1 = abs(yt-p1(2)-tan(pi/2.+pi/36.)*(xt-p1(1)))/sqrt(1+(tan(pi/2.+pi/36.))**2 )
                !-- p1 --p3
                h2 = abs(yt-p1(2)-tan(pi + pi/36. - pi/3.)*(xt-p1(1)))/sqrt(1+(tan(pi + pi/36. - pi/3.))**2 )
                !-- p2 --p3
                h3 = abs(yt-p2(2)-tan(pi/36.)*(xt-p2(1)))/sqrt(1+(tan(pi/36.))**2)

                p1p2 = sqrt((p2(2)- p1(2))**2 + (p2(1)-p1(1))**2)
                p1p3 = sqrt((p3(2)- p1(2))**2 + (p3(1)-p1(1))**2)
                p2p3 = sqrt((p2(2)- p3(2))**2 + (p2(1)-p3(1))**2)

!                print*, h1, h2, h3


                if ((abs((h1*p1p2 + h2*p1p3 + h3*p2p3)- (delta*p2p3)) .le. 1e-12)) then
                  v(ii,jj) = 0.
                endif
      endif



     if (phase .le. 0.) then
!----triangle leading edge
                p1(1) = xo + delta/2. * sin(pi/36.)
                p1(2) = yo + delta/2. * cos(pi/36.)

                p2(1) = xo - delta/2. * sin(pi/36.)
                p2(2) = yo - delta/2. * cos(pi/36.)

                p3(1) = p2(1) - delta * 1./tan(pi/3.) * cos(pi/36.)
                p3(2) = p1(2) - delta / cos(pi/6.) * sin(pi/3. - pi/36.)


!                xt = xo
!                yt = yo  

                !-- p1 --p2
                h1 = abs(yt-p1(2)-tan(pi/2. - pi/36.)*(xt-p1(1)))/sqrt(1+(tan(pi/2. - pi/36.))**2 )
                !-- p1 --p3
                h2 = abs(yt-p1(2)-tan(pi/2. - pi/36. - pi/6.)*(xt-p1(1)))/sqrt(1+(tan(pi/2. - pi/36. - pi/6.))**2 )
                !-- p2 --p3
                h3 = abs(yt-p2(2)-tan(pi - pi/36.)*(xt-p2(1)))/sqrt(1+(tan(pi - pi/36.))**2)

                p1p2 = sqrt((p2(2)- p1(2))**2 + (p2(1)-p1(1))**2)
                p1p3 = sqrt((p3(2)- p1(2))**2 + (p3(1)-p1(1))**2)
                p2p3 = sqrt((p2(2)- p3(2))**2 + (p2(1)-p3(1))**2)

!                print*, h1, h2, h3


                if ((abs((h1*p1p2 + h2*p1p3 + h3*p2p3)- (delta*p2p3)) .le. 1e-12)) then
                  v(ii,jj) = 0.
                endif
!----triangle trailing edge
                p1(1) = xe - delta/2. * sin(pi/36.)
                p1(2) = ye - delta/2. * cos(pi/36.)

                p2(1) = xe + delta/2. * sin(pi/36.)
                p2(2) = ye + delta/2. * cos(pi/36.)

                p3(1) = p2(1) + delta * 1./tan(pi/3.) * cos(pi/36.)
                p3(2) = p1(2) + delta / cos(pi/6.) * sin(pi/3. - pi/36.)

!                xt = xo
!                yt = yo  

                !-- p1 --p2
                h1 = abs(yt-p1(2)-tan(pi/2. - pi/36.)*(xt-p1(1)))/sqrt(1+(tan(pi/2. - pi/36.))**2 )
                !-- p1 --p3
                h2 = abs(yt-p1(2)-tan(pi/2. - pi/36. - pi/6.)*(xt-p1(1)))/sqrt(1+(tan(pi/2. - pi/36. - pi/6.))**2 )
                !-- p2 --p3
                h3 = abs(yt-p2(2)-tan(pi - pi/36.)*(xt-p2(1)))/sqrt(1+(tan(pi - pi/36.))**2)

                p1p2 = sqrt((p2(2)- p1(2))**2 + (p2(1)-p1(1))**2)
                p1p3 = sqrt((p3(2)- p1(2))**2 + (p3(1)-p1(1))**2)
                p2p3 = sqrt((p2(2)- p3(2))**2 + (p2(1)-p3(1))**2)

!                print*, h1, h2, h3


                if ((abs((h1*p1p2 + h2*p1p3 + h3*p2p3)- (delta*p2p3)) .le. 1e-12)) then
                  v(ii,jj) = 0.
                endif
      endif





             enddo
         enddo

      end subroutine oscillator_projection_vel_IBM_line_2D


      subroutine oscillator_projection_vel_IBM_line_2D_arbit_incidence(u, v, xo , yo, l, theta, angle_amplitude, freq, t)


         integer             ::     ii,jj,kk



         integer             ::     x_index,  y_index

         integer             ::     n_ghost                           !   number of ghost cells

         real(8), intent(in)    ::     theta, angle_amplitude, freq, t
         real(8)                ::     ghost_cells(2, 10000)               !   ghost cells
         real(8)                ::     image_cells(2, 10000)               !   image cells
         real(8), intent(in)    ::     xo, yo,l                            !   origin point and end point of the line
         real(8)                ::     xt, yt, xg, yg, xin, yin, xip, yip  !   test, ghost,intercept and image point coordinates 
         real(8)                ::     xe, ye
         real(8)                ::     m1, m2                              !   slope of line and the normal
         real(8)                ::     g_dist                              !   ghost cell distance from the line
         real(8)                ::     u(N1p,N2p),  v(N1p,N2p)
         real(8)                ::     x1 , y1, x2, y2, x3, y3, x4, y4
         integer             ::     Xindex, Yindex
         character*20        ::     file_name1

         real(8)                ::     radius

         real(8)                ::     phase, omega


         real(8)                ::     p1(2), p2(2), p3(2)
         real(8)                ::     h1, h2, h3, p1p2, p1p3, p2p3

         real(8),parameter      ::     pi = 3.14159265359 
         real(8), parameter     ::     delta = 0.03 

! calculating computational geometry which has to be done outside of the loop

         phase = theta + angle_amplitude * sin (2. * pi * freq * t)

         xe = xo + l * cos(phase)
         ye = yo + l * sin(phase)


         m1 = (ye - yo) / (xe -xo)
         m2 = -1. / m1


         n_ghost = 0
         do ii = S11B, N11B
            do jj = S21B, N21B
! on the u grid
                xt = x1u(ii)
                yt = x2p(jj)

                if ((xt .ge. xo)  .and.  (xt .le. xe)) then
                    g_dist = abs(yt - yo - m1*(xt - xo))/sqrt(m1**2 + 1.)
                    !g_dist = abs(yt - yo - m1*(xt - xo))/sqrt(m1**2 + 1.)
                    !if (((g_dist) .ge. -1.5*dx1p(ii)) .and.  (g_dist .le. 0.0)) then
                    if (g_dist .le. delta/2) then 
                        xg = xt
                        yg = yt
                        !n_ghost = n_ghost + 1

                        radius = sqrt((xt-xo)**2 + (yt-yo)**2)

                        !call Coordfinder(3,x1p,x2p,x1u,x2v,xg,yg,Xindex,Yindex,N1p,N2p)
                        omega =  2. * pi * freq * angle_amplitude * cos(2.*pi*freq*t)

                      !  u(ii,jj) = -sign(1. , cos(2. * pi * freq * t)) * radius * omega *  sin(phase)  
                  
                        u(ii,jj) = 0.
                    endif

! the smooth ends
                     

                endif

!----round leading edge
!                    g_dist = sqrt((xt - xo)**2 + (yt-yo)**2) 
!                    if (g_dist .le. 0.0075) then
!                        u(ii,jj) = 0.
!                    endif
! 
!                    g_dist = sqrt((xt - xe)**2 + (yt-ye)**2)
!                    if (g_dist .le. 0.0075) then
!                        u(ii,jj) = 0.
!                    endif

     if (yo .ge. 0.5) then
!----triangle leading edge
                p1(1) = xo + delta/2. * sin(theta)
                p1(2) = yo - delta/2. * cos(theta)
                
                p2(1) = xo - delta/2. * sin(theta)
                p2(2) = yo + delta/2. * cos(theta)

                p3(1) = p2(1) - delta * 1./tan(pi/3.) * cos(theta)
                p3(2) = p1(2) + delta / cos(pi/6.) * sin(pi/3. - theta)


!                xt = xo
!                yt = yo  

                !-- p1 --p2
                h1 = abs(yt-p1(2)-tan(pi/2.+theta)*(xt-p1(1)))/sqrt(1+(tan(pi/2.+theta))**2 )
                !-- p1 --p3
                h2 = abs(yt-p1(2)-tan(pi + theta - pi/3.)*(xt-p1(1)))/sqrt(1+(tan(pi + theta - pi/3.))**2 )
                !-- p2 --p3
                h3 = abs(yt-p2(2)-tan(theta)*(xt-p2(1)))/sqrt(1+(tan(theta))**2)

                p1p2 = sqrt((p2(2)- p1(2))**2 + (p2(1)-p1(1))**2)
                p1p3 = sqrt((p3(2)- p1(2))**2 + (p3(1)-p1(1))**2)
                p2p3 = sqrt((p2(2)- p3(2))**2 + (p2(1)-p3(1))**2)

!                print*, h1, h2, h3


                if ((abs((h1*p1p2 + h2*p1p3 + h3*p2p3)- (delta*p2p3)) .le. 1e-12)) then
                  u(ii,jj) = 0.
                endif
!----triangle trailing edge
                p1(1) = xe - delta/2. * sin(theta)
                p1(2) = ye + delta/2. * cos(theta)

                p2(1) = xe + delta/2. * sin(theta)
                p2(2) = ye - delta/2. * cos(theta)

                p3(1) = p2(1) + delta * 1./tan(pi/3.) * cos(theta)
                p3(2) = p1(2) - delta / cos(pi/6.) * sin(pi/3. - theta)

!                xt = xo
!                yt = yo  

                !-- p1 --p2
                h1 = abs(yt-p1(2)-tan(pi/2.+theta)*(xt-p1(1)))/sqrt(1+(tan(pi/2.+theta))**2 )
                !-- p1 --p3
                h2 = abs(yt-p1(2)-tan(pi + theta - pi/3.)*(xt-p1(1)))/sqrt(1+(tan(pi + theta - pi/3.))**2 )
                !-- p2 --p3
                h3 = abs(yt-p2(2)-tan(theta)*(xt-p2(1)))/sqrt(1+(tan(theta))**2)

                p1p2 = sqrt((p2(2)- p1(2))**2 + (p2(1)-p1(1))**2)
                p1p3 = sqrt((p3(2)- p1(2))**2 + (p3(1)-p1(1))**2)
                p2p3 = sqrt((p2(2)- p3(2))**2 + (p2(1)-p3(1))**2)

!                print*, h1, h2, h3


                if ((abs((h1*p1p2 + h2*p1p3 + h3*p2p3)- (delta*p2p3)) .le. 1e-12)) then
                  u(ii,jj) = 0.
                endif
      endif



     if (yo .le. 0.5) then
!----triangle leading edge
                p1(1) = xo + delta/2. * sin(theta)
                p1(2) = yo + delta/2. * cos(theta)

                p2(1) = xo - delta/2. * sin(theta)
                p2(2) = yo - delta/2. * cos(theta)

                p3(1) = p2(1) - delta * 1./tan(pi/3.) * cos(theta)
                p3(2) = p1(2) - delta / cos(pi/6.) * sin(pi/3. - theta)


!                xt = xo
!                yt = yo  

                !-- p1 --p2
                h1 = abs(yt-p1(2)-tan(pi/2. - theta)*(xt-p1(1)))/sqrt(1+(tan(pi/2. - theta))**2 )
                !-- p1 --p3
                h2 = abs(yt-p1(2)-tan(pi/2. - theta - pi/6.)*(xt-p1(1)))/sqrt(1+(tan(pi/2. - theta - pi/6.))**2 )
                !-- p2 --p3
                h3 = abs(yt-p2(2)-tan(pi - theta)*(xt-p2(1)))/sqrt(1+(tan(pi - theta))**2)

                p1p2 = sqrt((p2(2)- p1(2))**2 + (p2(1)-p1(1))**2)
                p1p3 = sqrt((p3(2)- p1(2))**2 + (p3(1)-p1(1))**2)
                p2p3 = sqrt((p2(2)- p3(2))**2 + (p2(1)-p3(1))**2)

!                print*, h1, h2, h3


                if ((abs((h1*p1p2 + h2*p1p3 + h3*p2p3)- (delta*p2p3)) .le. 1e-12)) then
                  u(ii,jj) = 0.
                endif
!----triangle trailing edge
                p1(1) = xe - delta/2. * sin(theta)
                p1(2) = ye - delta/2. * cos(theta)

                p2(1) = xe + delta/2. * sin(theta)
                p2(2) = ye + delta/2. * cos(theta)

                p3(1) = p2(1) + delta * 1./tan(pi/3.) * cos(theta)
                p3(2) = p1(2) + delta / cos(pi/6.) * sin(pi/3. - theta)

!                xt = xo
!                yt = yo  

                !-- p1 --p2
                h1 = abs(yt-p1(2)-tan(pi/2. - theta)*(xt-p1(1)))/sqrt(1+(tan(pi/2. - theta))**2 )
                !-- p1 --p3
                h2 = abs(yt-p1(2)-tan(pi/2. - theta - pi/6.)*(xt-p1(1)))/sqrt(1+(tan(pi/2. - theta - pi/6.))**2 )
                !-- p2 --p3
                h3 = abs(yt-p2(2)-tan(pi - theta)*(xt-p2(1)))/sqrt(1+(tan(pi - theta))**2)

                p1p2 = sqrt((p2(2)- p1(2))**2 + (p2(1)-p1(1))**2)
                p1p3 = sqrt((p3(2)- p1(2))**2 + (p3(1)-p1(1))**2)
                p2p3 = sqrt((p2(2)- p3(2))**2 + (p2(1)-p3(1))**2)

!                print*, h1, h2, h3


                if ((abs((h1*p1p2 + h2*p1p3 + h3*p2p3)- (delta*p2p3)) .le. 1e-12)) then
                  u(ii,jj) = 0.
                endif
      endif



! on the v grid

                xt = x1p(ii)
                yt = x2v(jj)
                if ((xt .ge. xo)  .and.  (xt .le. xe)) then
                    g_dist = abs(yt - yo - m1*(xt - xo))/sqrt(m1**2 + 1.)
                    !g_dist = (yt - yo - m1*(xt - xo))/sqrt(m1**2 + 1.)

                    !if (((g_dist) .ge. -1.5*dx1p(ii)) .and. (g_dist .le. 0.)) then
                    if (g_dist .le. delta/2) then
                      !  n_ghost = n_ghost+1

                        xg = xt
                        yg = yt
                        radius = sqrt((xt-xo)**2 + (yt-yo)**2)

                        !call Coordfinder(3,x1p,x2p,x1u,x2v,xg,yg,Xindex,Yindex,N1p,N2p)

                     !   phase = theta + angle_amplitude * sin (2. * pi * freq * t)
                        omega =  2. * pi * freq * angle_amplitude * cos(2.*pi*freq*t)
                    
                        v(ii,jj) = 0.                         
 
                     !!   v(ii,jj) = sign(1., cos(2. * pi * freq * t)) * radius * omega  * cos(phase)
                        !v(ii,jj) = sign(1., cos(2. * pi * freq * t)) * radius * 2. * pi * freq * cos(phase)

!                        call Coordfinder(2,x1p,x2p,x1u,x2v,xg,yg,Xindex,Yindex,N1p,N2p)

!                        v(ii,jj) = 0.

                    endif
               endif
!----round leading edge
!                    g_dist = sqrt((xt - xo)**2 + (yt-yo)**2) 
!                    if (g_dist .le. 0.015) then
!                        v(ii,jj) = 0.
!                    endif
!                    g_dist = sqrt((xt - xe)**2 + (yt-ye)**2)
!                    if (g_dist .le. 0.015) then
!                        v(ii,jj) = 0.
!                    endif
!----triangle leading edge

     if (yo .ge. 0.5) then
!----triangle leading edge
                p1(1) = xo + delta/2. * sin(theta)
                p1(2) = yo - delta/2. * cos(theta)
                
                p2(1) = xo - delta/2. * sin(theta)
                p2(2) = yo + delta/2. * cos(theta)

                p3(1) = p2(1) - delta * 1./tan(pi/3.) * cos(theta)
                p3(2) = p1(2) + delta / cos(pi/6.) * sin(pi/3. - theta)


!                xt = xo
!                yt = yo  

                !-- p1 --p2
                h1 = abs(yt-p1(2)-tan(pi/2.+theta)*(xt-p1(1)))/sqrt(1+(tan(pi/2.+theta))**2 )
                !-- p1 --p3
                h2 = abs(yt-p1(2)-tan(pi + theta - pi/3.)*(xt-p1(1)))/sqrt(1+(tan(pi + theta - pi/3.))**2 )
                !-- p2 --p3
                h3 = abs(yt-p2(2)-tan(theta)*(xt-p2(1)))/sqrt(1+(tan(theta))**2)

                p1p2 = sqrt((p2(2)- p1(2))**2 + (p2(1)-p1(1))**2)
                p1p3 = sqrt((p3(2)- p1(2))**2 + (p3(1)-p1(1))**2)
                p2p3 = sqrt((p2(2)- p3(2))**2 + (p2(1)-p3(1))**2)

!                print*, h1, h2, h3


                if ((abs((h1*p1p2 + h2*p1p3 + h3*p2p3)- (delta*p2p3)) .le. 1e-12)) then
                  v(ii,jj) = 0.
                endif
!----triangle trailing edge
                p1(1) = xe - delta/2. * sin(theta)
                p1(2) = ye + delta/2. * cos(theta)

                p2(1) = xe + delta/2. * sin(theta)
                p2(2) = ye - delta/2. * cos(theta)

                p3(1) = p2(1) + delta * 1./tan(pi/3.) * cos(theta)
                p3(2) = p1(2) - delta / cos(pi/6.) * sin(pi/3. - theta)

!                xt = xo
!                yt = yo  

                !-- p1 --p2
                h1 = abs(yt-p1(2)-tan(pi/2.+theta)*(xt-p1(1)))/sqrt(1+(tan(pi/2.+theta))**2 )
                !-- p1 --p3
                h2 = abs(yt-p1(2)-tan(pi + theta - pi/3.)*(xt-p1(1)))/sqrt(1+(tan(pi + theta - pi/3.))**2 )
                !-- p2 --p3
                h3 = abs(yt-p2(2)-tan(theta)*(xt-p2(1)))/sqrt(1+(tan(theta))**2)

                p1p2 = sqrt((p2(2)- p1(2))**2 + (p2(1)-p1(1))**2)
                p1p3 = sqrt((p3(2)- p1(2))**2 + (p3(1)-p1(1))**2)
                p2p3 = sqrt((p2(2)- p3(2))**2 + (p2(1)-p3(1))**2)

!                print*, h1, h2, h3


                if ((abs((h1*p1p2 + h2*p1p3 + h3*p2p3)- (delta*p2p3)) .le. 1e-12)) then
                  v(ii,jj) = 0.
                endif
      endif



     if (yo .le. 0.5) then
!----triangle leading edge
                p1(1) = xo + delta/2. * sin(theta)
                p1(2) = yo + delta/2. * cos(theta)

                p2(1) = xo - delta/2. * sin(theta)
                p2(2) = yo - delta/2. * cos(theta)

                p3(1) = p2(1) - delta * 1./tan(pi/3.) * cos(theta)
                p3(2) = p1(2) - delta / cos(pi/6.) * sin(pi/3. - theta)


!                xt = xo
!                yt = yo  

                !-- p1 --p2
                h1 = abs(yt-p1(2)-tan(pi/2. - theta)*(xt-p1(1)))/sqrt(1+(tan(pi/2. - theta))**2 )
                !-- p1 --p3
                h2 = abs(yt-p1(2)-tan(pi/2. - theta - pi/6.)*(xt-p1(1)))/sqrt(1+(tan(pi/2. - theta - pi/6.))**2 )
                !-- p2 --p3
                h3 = abs(yt-p2(2)-tan(pi - theta)*(xt-p2(1)))/sqrt(1+(tan(pi - theta))**2)

                p1p2 = sqrt((p2(2)- p1(2))**2 + (p2(1)-p1(1))**2)
                p1p3 = sqrt((p3(2)- p1(2))**2 + (p3(1)-p1(1))**2)
                p2p3 = sqrt((p2(2)- p3(2))**2 + (p2(1)-p3(1))**2)

!                print*, h1, h2, h3


                if ((abs((h1*p1p2 + h2*p1p3 + h3*p2p3)- (delta*p2p3)) .le. 1e-12)) then
                  v(ii,jj) = 0.
                endif
!----triangle trailing edge
                p1(1) = xe - delta/2. * sin(theta)
                p1(2) = ye - delta/2. * cos(theta)

                p2(1) = xe + delta/2. * sin(theta)
                p2(2) = ye + delta/2. * cos(theta)

                p3(1) = p2(1) + delta * 1./tan(pi/3.) * cos(theta)
                p3(2) = p1(2) + delta / cos(pi/6.) * sin(pi/3. - theta)

!                xt = xo
!                yt = yo  

                !-- p1 --p2
                h1 = abs(yt-p1(2)-tan(pi/2. - theta)*(xt-p1(1)))/sqrt(1+(tan(pi/2. - theta))**2 )
                !-- p1 --p3
                h2 = abs(yt-p1(2)-tan(pi/2. - theta - pi/6.)*(xt-p1(1)))/sqrt(1+(tan(pi/2. - theta - pi/6.))**2 )
                !-- p2 --p3
                h3 = abs(yt-p2(2)-tan(pi - theta)*(xt-p2(1)))/sqrt(1+(tan(pi - theta))**2)

                p1p2 = sqrt((p2(2)- p1(2))**2 + (p2(1)-p1(1))**2)
                p1p3 = sqrt((p3(2)- p1(2))**2 + (p3(1)-p1(1))**2)
                p2p3 = sqrt((p2(2)- p3(2))**2 + (p2(1)-p3(1))**2)

!                print*, h1, h2, h3


                if ((abs((h1*p1p2 + h2*p1p3 + h3*p2p3)- (delta*p2p3)) .le. 1e-12)) then
                  v(ii,jj) = 0.
                endif
      endif





             enddo
         enddo

      end subroutine oscillator_projection_vel_IBM_line_2D_arbit_incidence



      subroutine oscillator_projection_vel_IBM_line_2D_onx(u, v, xo , yo, l, theta, angle_amplitude, freq, t)


         integer             ::     ii,jj,kk



         integer             ::     x_index,  y_index

         integer             ::     n_ghost                           !   number of ghost cells

         real(8), intent(in)    ::     theta, angle_amplitude, freq, t
         real(8)                ::     ghost_cells(2, 10000)               !   ghost cells
         real(8)                ::     image_cells(2, 10000)               !   image cells
         real(8), intent(in)    ::     xo, yo,l                            !   origin point and end point of the line
         real(8)                ::     xt, yt, xg, yg, xin, yin, xip, yip  !   test, ghost,intercept and image point coordinates 
         real(8)                ::     xe, ye
         real(8)                ::     m1, m2                              !   slope of line and the normal
         real(8)                ::     g_dist                              !   ghost cell distance from the line
         real(8)                ::     u(N1p,N2p),  v(N1p,N2p)
         real(8)                ::     x1 , y1, x2, y2, x3, y3, x4, y4
         integer             ::     Xindex, Yindex
         character*20        ::     file_name1

         real(8)                ::     radius

         real(8)                ::     phase, omega

         real(8)                ::     attack

         real(8)                ::     p1(2), p2(2), p3(2)
         real(8)                ::     h1, h2, h3, p1p2, p1p3, p2p3

         real(8),parameter      ::     pi = 3.14159265359 
         real(8), parameter     ::     delta = 0.03 

! calculating computational geometry which has to be done outside of the loop

         phase = theta + angle_amplitude * sin (2. * pi * freq * t)

         xe = xo + l * cos(phase)
         ye = yo + l * sin(phase)


         m1 = (ye - yo) / (xe -xo)
         m2 = -1. / m1


         n_ghost = 0
         do ii = S11B, N11B
            do jj = S21B, N21B
! on the u grid
                xt = x1u(ii)
                yt = x2p(jj)

                if ((xt .ge. xo)  .and.  (xt .le. xe)) then
                    g_dist = abs(yt - yo - m1*(xt - xo))/sqrt(m1**2 + 1.)
                    !g_dist = abs(yt - yo - m1*(xt - xo))/sqrt(m1**2 + 1.)
                    !if (((g_dist) .ge. -1.5*dx1p(ii)) .and.  (g_dist .le. 0.0)) then
                    if (g_dist .le. delta/2) then 
                        xg = xt
                        yg = yt
                        !n_ghost = n_ghost + 1

                        radius = sqrt((xt-xo)**2 + (yt-yo)**2)

                        !call Coordfinder(3,x1p,x2p,x1u,x2v,xg,yg,Xindex,Yindex,N1p,N2p)
                        omega =  2. * pi * freq * angle_amplitude * cos(2.*pi*freq*t)

                      !  u(ii,jj) = -sign(1. , cos(2. * pi * freq * t)) * radius * omega *  sin(phase)  
                  
                        u(ii,jj) = 0.
                    endif

! the smooth ends
                     

                endif

!----round leading edge
!                    g_dist = sqrt((xt - xo)**2 + (yt-yo)**2) 
!                    if (g_dist .le. 0.0075) then
!                        u(ii,jj) = 0.
!                    endif
! 
!                    g_dist = sqrt((xt - xe)**2 + (yt-ye)**2)
!                    if (g_dist .le. 0.0075) then
!                        u(ii,jj) = 0.
!                    endif

     attack = abs(phase)


     if (phase .ge. 0.) then
!----triangle leading edge
                p1(1) = xo + delta/2. * sin(attack)
                p1(2) = yo - delta/2. * cos(attack)
                
                p2(1) = xo - delta/2. * sin(attack)
                p2(2) = yo + delta/2. * cos(attack)

                p3(1) = p2(1) - delta * 1./tan(pi/3.) * cos(attack)
                p3(2) = p1(2) + delta / cos(pi/6.) * sin(pi/3. - attack)


!                xt = xo
!                yt = yo  

                !-- p1 --p2
                h1 = abs(yt-p1(2)-tan(pi/2.+attack)*(xt-p1(1)))/sqrt(1+(tan(pi/2.+attack))**2 )
                !-- p1 --p3
                h2 = abs(yt-p1(2)-tan(pi + attack - pi/3.)*(xt-p1(1)))/sqrt(1+(tan(pi + attack - pi/3.))**2 )
                !-- p2 --p3
                h3 = abs(yt-p2(2)-tan(attack)*(xt-p2(1)))/sqrt(1+(tan(attack))**2)

                p1p2 = sqrt((p2(2)- p1(2))**2 + (p2(1)-p1(1))**2)
                p1p3 = sqrt((p3(2)- p1(2))**2 + (p3(1)-p1(1))**2)
                p2p3 = sqrt((p2(2)- p3(2))**2 + (p2(1)-p3(1))**2)

!                print*, h1, h2, h3


                if ((abs((h1*p1p2 + h2*p1p3 + h3*p2p3)- (delta*p2p3)) .le. 1e-12)) then
                  u(ii,jj) = 0.
                endif
!----triangle trailing edge
                p1(1) = xe - delta/2. * sin(attack)
                p1(2) = ye + delta/2. * cos(attack)

                p2(1) = xe + delta/2. * sin(attack)
                p2(2) = ye - delta/2. * cos(attack)

                p3(1) = p2(1) + delta * 1./tan(pi/3.) * cos(attack)
                p3(2) = p1(2) - delta / cos(pi/6.) * sin(pi/3. - attack)

!                xt = xo
!                yt = yo  

                !-- p1 --p2
                h1 = abs(yt-p1(2)-tan(pi/2.+attack)*(xt-p1(1)))/sqrt(1+(tan(pi/2.+attack))**2 )
                !-- p1 --p3
                h2 = abs(yt-p1(2)-tan(pi + attack - pi/3.)*(xt-p1(1)))/sqrt(1+(tan(pi + attack - pi/3.))**2 )
                !-- p2 --p3
                h3 = abs(yt-p2(2)-tan(attack)*(xt-p2(1)))/sqrt(1+(tan(pi/attack))**2)

                p1p2 = sqrt((p2(2)- p1(2))**2 + (p2(1)-p1(1))**2)
                p1p3 = sqrt((p3(2)- p1(2))**2 + (p3(1)-p1(1))**2)
                p2p3 = sqrt((p2(2)- p3(2))**2 + (p2(1)-p3(1))**2)

!                print*, h1, h2, h3


                if ((abs((h1*p1p2 + h2*p1p3 + h3*p2p3)- (delta*p2p3)) .le. 1e-12)) then
                  u(ii,jj) = 0.
                endif
      endif



     if (phase .le. 0.) then
!----triangle leading edge
                p1(1) = xo + delta/2. * sin(attack)
                p1(2) = yo + delta/2. * cos(attack)

                p2(1) = xo - delta/2. * sin(attack)
                p2(2) = yo - delta/2. * cos(attack)

                p3(1) = p2(1) - delta * 1./tan(pi/3.) * cos(attack)
                p3(2) = p1(2) - delta / cos(pi/6.) * sin(pi/3. - attack)


!                xt = xo
!                yt = yo  

                !-- p1 --p2
                h1 = abs(yt-p1(2)-tan(pi/2. - attack)*(xt-p1(1)))/sqrt(1+(tan(pi/2. - attack))**2 )
                !-- p1 --p3
                h2 = abs(yt-p1(2)-tan(pi/2. - attack - pi/6.)*(xt-p1(1)))/sqrt(1+(tan(pi/2. - attack - pi/6.))**2 )
                !-- p2 --p3
                h3 = abs(yt-p2(2)-tan(pi - attack)*(xt-p2(1)))/sqrt(1+(tan(pi - attack))**2)

                p1p2 = sqrt((p2(2)- p1(2))**2 + (p2(1)-p1(1))**2)
                p1p3 = sqrt((p3(2)- p1(2))**2 + (p3(1)-p1(1))**2)
                p2p3 = sqrt((p2(2)- p3(2))**2 + (p2(1)-p3(1))**2)

!                print*, h1, h2, h3


                if ((abs((h1*p1p2 + h2*p1p3 + h3*p2p3)- (delta*p2p3)) .le. 1e-12)) then
                  u(ii,jj) = 0.
                endif
!----triangle trailing edge
                p1(1) = xe - delta/2. * sin(attack)
                p1(2) = ye - delta/2. * cos(attack)

                p2(1) = xe + delta/2. * sin(attack)
                p2(2) = ye + delta/2. * cos(attack)

                p3(1) = p2(1) + delta * 1./tan(pi/3.) * cos(attack)
                p3(2) = p1(2) + delta / cos(pi/6.) * sin(pi/3. - attack)

!                xt = xo
!                yt = yo  

                !-- p1 --p2
                h1 = abs(yt-p1(2)-tan(pi/2. - attack)*(xt-p1(1)))/sqrt(1+(tan(pi/2. - attack))**2 )
                !-- p1 --p3
                h2 = abs(yt-p1(2)-tan(pi/2. - attack - pi/6.)*(xt-p1(1)))/sqrt(1+(tan(pi/2. - attack - pi/6.))**2 )
                !-- p2 --p3
                h3 = abs(yt-p2(2)-tan(pi - attack)*(xt-p2(1)))/sqrt(1+(tan(pi - attack))**2)

                p1p2 = sqrt((p2(2)- p1(2))**2 + (p2(1)-p1(1))**2)
                p1p3 = sqrt((p3(2)- p1(2))**2 + (p3(1)-p1(1))**2)
                p2p3 = sqrt((p2(2)- p3(2))**2 + (p2(1)-p3(1))**2)

!                print*, h1, h2, h3


                if ((abs((h1*p1p2 + h2*p1p3 + h3*p2p3)- (delta*p2p3)) .le. 1e-12)) then
                  u(ii,jj) = 0.
                endif
      endif



! on the v grid

      ! phase = theta

                xt = x1p(ii)
                yt = x2v(jj)
                if ((xt .ge. xo)  .and.  (xt .le. xe)) then
                    g_dist = abs(yt - yo - m1*(xt - xo))/sqrt(m1**2 + 1.)
                    !g_dist = (yt - yo - m1*(xt - xo))/sqrt(m1**2 + 1.)

                    !if (((g_dist) .ge. -1.5*dx1p(ii)) .and. (g_dist .le. 0.)) then
                    if (g_dist .le. delta/2) then
                      !  n_ghost = n_ghost+1

                        xg = xt
                        yg = yt
                        radius = sqrt((xt-xo)**2 + (yt-yo)**2)

                        !call Coordfinder(3,x1p,x2p,x1u,x2v,xg,yg,Xindex,Yindex,N1p,N2p)

                     !   phase = theta + angle_amplitude * sin (2. * pi * freq * t)
                        omega =  2. * pi * freq * angle_amplitude * cos(2.*pi*freq*t)
                    
                        v(ii,jj) = 0.                         
 
                     !!   v(ii,jj) = sign(1., cos(2. * pi * freq * t)) * radius * omega  * cos(phase)
                        !v(ii,jj) = sign(1., cos(2. * pi * freq * t)) * radius * 2. * pi * freq * cos(phase)

!                        call Coordfinder(2,x1p,x2p,x1u,x2v,xg,yg,Xindex,Yindex,N1p,N2p)

!                        v(ii,jj) = 0.

                    endif
               endif
!----round leading edge
!                    g_dist = sqrt((xt - xo)**2 + (yt-yo)**2) 
!                    if (g_dist .le. 0.015) then
!                        v(ii,jj) = 0.
!                    endif
!                    g_dist = sqrt((xt - xe)**2 + (yt-ye)**2)
!                    if (g_dist .le. 0.015) then
!                        v(ii,jj) = 0.
!                    endif
!----triangle leading edge

     if (phase .ge. 0.) then

!----triangle leading edge
                p1(1) = xo + delta/2. * sin(attack)
                p1(2) = yo - delta/2. * cos(attack)
                
                p2(1) = xo - delta/2. * sin(attack)
                p2(2) = yo + delta/2. * cos(attack)

                p3(1) = p2(1) - delta * 1./tan(pi/3.) * cos(attack)
                p3(2) = p1(2) + delta / cos(pi/6.) * sin(pi/3. - attack)


!                xt = xo
!                yt = yo  

                !-- p1 --p2
                h1 = abs(yt-p1(2)-tan(pi/2.+attack)*(xt-p1(1)))/sqrt(1+(tan(pi/2.+attack))**2 )
                !-- p1 --p3
                h2 = abs(yt-p1(2)-tan(pi + attack - pi/3.)*(xt-p1(1)))/sqrt(1+(tan(pi + attack - pi/3.))**2 )
                !-- p2 --p3
                h3 = abs(yt-p2(2)-tan(attack)*(xt-p2(1)))/sqrt(1+(tan(attack))**2)

                p1p2 = sqrt((p2(2)- p1(2))**2 + (p2(1)-p1(1))**2)
                p1p3 = sqrt((p3(2)- p1(2))**2 + (p3(1)-p1(1))**2)
                p2p3 = sqrt((p2(2)- p3(2))**2 + (p2(1)-p3(1))**2)

!                print*, h1, h2, h3


                if ((abs((h1*p1p2 + h2*p1p3 + h3*p2p3)- (delta*p2p3)) .le. 1e-12)) then
                  v(ii,jj) = 0.
                endif
!----triangle trailing edge
                p1(1) = xe - delta/2. * sin(attack)
                p1(2) = ye + delta/2. * cos(attack)

                p2(1) = xe + delta/2. * sin(attack)
                p2(2) = ye - delta/2. * cos(attack)

                p3(1) = p2(1) + delta * 1./tan(pi/3.) * cos(attack)
                p3(2) = p1(2) - delta / cos(pi/6.) * sin(pi/3. - attack)

!                xt = xo
!                yt = yo  

                !-- p1 --p2
                h1 = abs(yt-p1(2)-tan(pi/2.+attack)*(xt-p1(1)))/sqrt(1+(tan(pi/2.+attack))**2 )
                !-- p1 --p3
                h2 = abs(yt-p1(2)-tan(pi + attack - pi/3.)*(xt-p1(1)))/sqrt(1+(tan(pi + attack - pi/3.))**2 )
                !-- p2 --p3
                h3 = abs(yt-p2(2)-tan(attack)*(xt-p2(1)))/sqrt(1+(tan(attack))**2)

                p1p2 = sqrt((p2(2)- p1(2))**2 + (p2(1)-p1(1))**2)
                p1p3 = sqrt((p3(2)- p1(2))**2 + (p3(1)-p1(1))**2)
                p2p3 = sqrt((p2(2)- p3(2))**2 + (p2(1)-p3(1))**2)

!                print*, h1, h2, h3


                if ((abs((h1*p1p2 + h2*p1p3 + h3*p2p3)- (delta*p2p3)) .le. 1e-12)) then
                  v(ii,jj) = 0.
                endif
      endif



     if (phase .le. 0.) then
!----triangle leading edge
                p1(1) = xo + delta/2. * sin(attack)
                p1(2) = yo + delta/2. * cos(attack)

                p2(1) = xo - delta/2. * sin(attack)
                p2(2) = yo - delta/2. * cos(attack)

                p3(1) = p2(1) - delta * 1./tan(pi/3.) * cos(attack)
                p3(2) = p1(2) - delta / cos(pi/6.) * sin(pi/3. - attack)


!                xt = xo
!                yt = yo  

                !-- p1 --p2
                h1 = abs(yt-p1(2)-tan(pi/2. - attack)*(xt-p1(1)))/sqrt(1+(tan(pi/2. - attack))**2 )
                !-- p1 --p3
                h2 = abs(yt-p1(2)-tan(pi/2. - attack - pi/6.)*(xt-p1(1)))/sqrt(1+(tan(pi/2. - attack - pi/6.))**2 )
                !-- p2 --p3
                h3 = abs(yt-p2(2)-tan(pi - attack)*(xt-p2(1)))/sqrt(1+(tan(pi - attack))**2)

                p1p2 = sqrt((p2(2)- p1(2))**2 + (p2(1)-p1(1))**2)
                p1p3 = sqrt((p3(2)- p1(2))**2 + (p3(1)-p1(1))**2)
                p2p3 = sqrt((p2(2)- p3(2))**2 + (p2(1)-p3(1))**2)

!                print*, h1, h2, h3


                if ((abs((h1*p1p2 + h2*p1p3 + h3*p2p3)- (delta*p2p3)) .le. 1e-12)) then
                  v(ii,jj) = 0.
                endif
!----triangle trailing edge
                p1(1) = xe - delta/2. * sin(attack)
                p1(2) = ye - delta/2. * cos(attack)

                p2(1) = xe + delta/2. * sin(attack)
                p2(2) = ye + delta/2. * cos(attack)

                p3(1) = p2(1) + delta * 1./tan(pi/3.) * cos(attack)
                p3(2) = p1(2) + delta / cos(pi/6.) * sin(pi/3. - attack)

!                xt = xo
!                yt = yo  

                !-- p1 --p2
                h1 = abs(yt-p1(2)-tan(pi/2. - attack)*(xt-p1(1)))/sqrt(1+(tan(pi/2. - attack))**2 )
                !-- p1 --p3
                h2 = abs(yt-p1(2)-tan(pi/2. - attack - pi/6.)*(xt-p1(1)))/sqrt(1+(tan(pi/2. - attack - pi/6.))**2 )
                !-- p2 --p3
                h3 = abs(yt-p2(2)-tan(pi - attack)*(xt-p2(1)))/sqrt(1+(tan(pi - attack))**2)

                p1p2 = sqrt((p2(2)- p1(2))**2 + (p2(1)-p1(1))**2)
                p1p3 = sqrt((p3(2)- p1(2))**2 + (p3(1)-p1(1))**2)
                p2p3 = sqrt((p2(2)- p3(2))**2 + (p2(1)-p3(1))**2)

!                print*, h1, h2, h3


                if ((abs((h1*p1p2 + h2*p1p3 + h3*p2p3)- (delta*p2p3)) .le. 1e-12)) then
                  v(ii,jj) = 0.
                endif
      endif





             enddo
         enddo

      end subroutine oscillator_projection_vel_IBM_line_2D_onx


      subroutine oscillator_projection_vel_IBM_line_2D_design(u, v, xo , yo, l, theta,alpha, angle_amplitude, freq, t)


         integer             ::     ii,jj,kk



         integer             ::     x_index,  y_index

         integer             ::     n_ghost                           !   number of ghost cells

         real(8), intent(in)    ::     theta, alpha, angle_amplitude, freq, t
         real(8)                ::     ghost_cells(2, 10000)               !   ghost cells
         real(8)                ::     image_cells(2, 10000)               !   image cells
         real(8), intent(in)    ::     xo, yo,l                            !   origin point and end point of the line
         real(8)                ::     xt, yt, xg, yg, xin, yin, xip, yip  !   test, ghost,intercept and image point coordinates 
         real(8)                ::     xe, ye
         real(8)                ::     m1, m2                              !   slope of line and the normal
         real(8)                ::     g_dist                              !   ghost cell distance from the line
         real(8)                ::     u(N1p,N2p),  v(N1p,N2p)
         real(8)                ::     x1 , y1, x2, y2, x3, y3, x4, y4
         integer             ::     Xindex, Yindex
         character*20        ::     file_name1

         real(8)                ::     radius

         real(8)                ::     phase, omega


         real(8)                ::     p1(2), p2(2), p3(2)
         real(8)                ::     h1, h2, h3, p1p2, p1p3, p2p3

         real(8),parameter      ::     pi = 3.14159265359 
         real(8), parameter     ::     delta = 0.03 

! calculating computational geometry which has to be done outside of the loop

         phase = theta + angle_amplitude * sin (2. * pi * freq * t)

         xe = xo + l * cos(phase)
         ye = yo + l * sin(phase)


         m1 = (ye - yo) / (xe -xo)
         m2 = -1. / m1


         n_ghost = 0
         do ii = S11B, N11B
            do jj = S21B, N21B
! on the u grid
                xt = x1u(ii)
                yt = x2p(jj)

                if ((xt .ge. xo)  .and.  (xt .le. xe)) then
                    g_dist = abs(yt - yo - m1*(xt - xo))/sqrt(m1**2 + 1.)
                    !g_dist = abs(yt - yo - m1*(xt - xo))/sqrt(m1**2 + 1.)
                    !if (((g_dist) .ge. -1.5*dx1p(ii)) .and.  (g_dist .le. 0.0)) then
                    if (g_dist .le. delta/2) then 
                        xg = xt
                        yg = yt
                        !n_ghost = n_ghost + 1

                        radius = sqrt((xt-xo)**2 + (yt-yo)**2)

                        !call Coordfinder(3,x1p,x2p,x1u,x2v,xg,yg,Xindex,Yindex,N1p,N2p)
                        omega =  2. * pi * freq * angle_amplitude * cos(2.*pi*freq*t)

                      !  u(ii,jj) = -sign(1. , cos(2. * pi * freq * t)) * radius * omega *  sin(phase)  
                  
                        u(ii,jj) = 0.
                    endif

! the smooth ends
                     

                endif

!----round leading edge
!                    g_dist = sqrt((xt - xo)**2 + (yt-yo)**2) 
!                    if (g_dist .le. 0.0075) then
!                        u(ii,jj) = 0.
!                    endif
! 
!                    g_dist = sqrt((xt - xe)**2 + (yt-ye)**2)
!                    if (g_dist .le. 0.0075) then
!                        u(ii,jj) = 0.
!                    endif

     if (phase .ge. 0.) then
!----triangle leading edge
                p1(1) = xo + delta/2. * sin(pi/36.)
                p1(2) = yo - delta/2. * cos(pi/36.)
                
                p2(1) = xo - delta/2. * sin(pi/36.)
                p2(2) = yo + delta/2. * cos(pi/36.)

                p3(1) = p2(1) - delta * 1./tan(alpha) * cos(pi/36.)
                p3(2) = p1(2) + delta / cos(pi/2.-alpha) * sin(alpha - pi/36.)


!                xt = xo
!                yt = yo  

                !-- p1 --p2
                h1 = abs(yt-p1(2)-tan(pi/2.+pi/36.)*(xt-p1(1)))/sqrt(1+(tan(pi/2.+pi/36.))**2 )
                !-- p1 --p3
                h2 = abs(yt-p1(2)-tan(pi + pi/36. - alpha)*(xt-p1(1)))/sqrt(1+(tan(pi + pi/36. - alpha))**2 )
                !-- p2 --p3
                h3 = abs(yt-p2(2)-tan(pi/36.)*(xt-p2(1)))/sqrt(1+(tan(pi/36.))**2)

                p1p2 = sqrt((p2(2)- p1(2))**2 + (p2(1)-p1(1))**2)
                p1p3 = sqrt((p3(2)- p1(2))**2 + (p3(1)-p1(1))**2)
                p2p3 = sqrt((p2(2)- p3(2))**2 + (p2(1)-p3(1))**2)

!                print*, h1, h2, h3


                if ((abs((h1*p1p2 + h2*p1p3 + h3*p2p3)- (delta*p2p3)) .le. 1e-12)) then
                  u(ii,jj) = 0.
                endif
!----triangle trailing edge
                p1(1) = xe - delta/2. * sin(pi/36.)
                p1(2) = ye + delta/2. * cos(pi/36.)

                p2(1) = xe + delta/2. * sin(pi/36.)
                p2(2) = ye - delta/2. * cos(pi/36.)

                p3(1) = p2(1) + delta * 1./tan(pi/3.) * cos(pi/36.)
                p3(2) = p1(2) - delta / cos(pi/6.) * sin(pi/3. - pi/36.)

!                xt = xo
!                yt = yo  

                !-- p1 --p2
                h1 = abs(yt-p1(2)-tan(pi/2.+pi/36.)*(xt-p1(1)))/sqrt(1+(tan(pi/2.+pi/36.))**2 )
                !-- p1 --p3
                h2 = abs(yt-p1(2)-tan(pi + pi/36. - pi/3.)*(xt-p1(1)))/sqrt(1+(tan(pi + pi/36. - pi/3.))**2 )
                !-- p2 --p3
                h3 = abs(yt-p2(2)-tan(pi/36.)*(xt-p2(1)))/sqrt(1+(tan(pi/36.))**2)

                p1p2 = sqrt((p2(2)- p1(2))**2 + (p2(1)-p1(1))**2)
                p1p3 = sqrt((p3(2)- p1(2))**2 + (p3(1)-p1(1))**2)
                p2p3 = sqrt((p2(2)- p3(2))**2 + (p2(1)-p3(1))**2)

!                print*, h1, h2, h3


                if ((abs((h1*p1p2 + h2*p1p3 + h3*p2p3)- (delta*p2p3)) .le. 1e-12)) then
                  u(ii,jj) = 0.
                endif
      endif



     if (phase .le. 0.) then
!----triangle leading edge
                p1(1) = xo + delta/2. * sin(pi/36.)
                p1(2) = yo + delta/2. * cos(pi/36.)

                p2(1) = xo - delta/2. * sin(pi/36.)
                p2(2) = yo - delta/2. * cos(pi/36.)

                p3(1) = p2(1) - delta * 1./tan(alpha) * cos(pi/36.)
                p3(2) = p1(2) - delta / cos(pi/2.-alpha) * sin(alpha - pi/36.)


!                xt = xo
!                yt = yo  

                !-- p1 --p2
                h1 = abs(yt-p1(2)-tan(pi/2. - pi/36.)*(xt-p1(1)))/sqrt(1+(tan(pi/2. - pi/36.))**2 )
                !-- p1 --p3
                h2 = abs(yt-p1(2)-tan(pi/2. - pi/36. - pi/2. + alpha)*(xt-p1(1)))/sqrt(1+(tan(pi/2. - pi/36. - pi/2. + alpha))**2 )
                !-- p2 --p3
                h3 = abs(yt-p2(2)-tan(pi - pi/36.)*(xt-p2(1)))/sqrt(1+(tan(pi - pi/36.))**2)

                p1p2 = sqrt((p2(2)- p1(2))**2 + (p2(1)-p1(1))**2)
                p1p3 = sqrt((p3(2)- p1(2))**2 + (p3(1)-p1(1))**2)
                p2p3 = sqrt((p2(2)- p3(2))**2 + (p2(1)-p3(1))**2)

!                print*, h1, h2, h3


                if ((abs((h1*p1p2 + h2*p1p3 + h3*p2p3)- (delta*p2p3)) .le. 1e-12)) then
                  u(ii,jj) = 0.
                endif
!----triangle trailing edge
                p1(1) = xe - delta/2. * sin(pi/36.)
                p1(2) = ye - delta/2. * cos(pi/36.)

                p2(1) = xe + delta/2. * sin(pi/36.)
                p2(2) = ye + delta/2. * cos(pi/36.)

                p3(1) = p2(1) + delta * 1./tan(pi/3.) * cos(pi/36.)
                p3(2) = p1(2) + delta / cos(pi/6.) * sin(pi/3. - pi/36.)

!                xt = xo
!                yt = yo  

                !-- p1 --p2
                h1 = abs(yt-p1(2)-tan(pi/2. - pi/36.)*(xt-p1(1)))/sqrt(1+(tan(pi/2. - pi/36.))**2 )
                !-- p1 --p3
                h2 = abs(yt-p1(2)-tan(pi/2. - pi/36. - pi/6.)*(xt-p1(1)))/sqrt(1+(tan(pi/2. - pi/36. - pi/6.))**2 )
                !-- p2 --p3
                h3 = abs(yt-p2(2)-tan(pi - pi/36.)*(xt-p2(1)))/sqrt(1+(tan(pi - pi/36.))**2)

                p1p2 = sqrt((p2(2)- p1(2))**2 + (p2(1)-p1(1))**2)
                p1p3 = sqrt((p3(2)- p1(2))**2 + (p3(1)-p1(1))**2)
                p2p3 = sqrt((p2(2)- p3(2))**2 + (p2(1)-p3(1))**2)

!                print*, h1, h2, h3


                if ((abs((h1*p1p2 + h2*p1p3 + h3*p2p3)- (delta*p2p3)) .le. 1e-12)) then
                  u(ii,jj) = 0.
                endif
      endif



! on the v grid

                xt = x1p(ii)
                yt = x2v(jj)
                if ((xt .ge. xo)  .and.  (xt .le. xe)) then
                    g_dist = abs(yt - yo - m1*(xt - xo))/sqrt(m1**2 + 1.)
                    !g_dist = (yt - yo - m1*(xt - xo))/sqrt(m1**2 + 1.)

                    !if (((g_dist) .ge. -1.5*dx1p(ii)) .and. (g_dist .le. 0.)) then
                    if (g_dist .le. delta/2) then
                      !  n_ghost = n_ghost+1

                        xg = xt
                        yg = yt
                        radius = sqrt((xt-xo)**2 + (yt-yo)**2)

                        !call Coordfinder(3,x1p,x2p,x1u,x2v,xg,yg,Xindex,Yindex,N1p,N2p)

                     !   phase = theta + angle_amplitude * sin (2. * pi * freq * t)
                        omega =  2. * pi * freq * angle_amplitude * cos(2.*pi*freq*t)
                    
                        v(ii,jj) = 0.                         
 
                     !!   v(ii,jj) = sign(1., cos(2. * pi * freq * t)) * radius * omega  * cos(phase)
                        !v(ii,jj) = sign(1., cos(2. * pi * freq * t)) * radius * 2. * pi * freq * cos(phase)

!                        call Coordfinder(2,x1p,x2p,x1u,x2v,xg,yg,Xindex,Yindex,N1p,N2p)

!                        v(ii,jj) = 0.

                    endif
               endif
!----round leading edge
!                    g_dist = sqrt((xt - xo)**2 + (yt-yo)**2) 
!                    if (g_dist .le. 0.015) then
!                        v(ii,jj) = 0.
!                    endif
!                    g_dist = sqrt((xt - xe)**2 + (yt-ye)**2)
!                    if (g_dist .le. 0.015) then
!                        v(ii,jj) = 0.
!                    endif
!----triangle leading edge

     if (phase .ge. 0.) then
!----triangle leading edge
                p1(1) = xo + delta/2. * sin(pi/36.)
                p1(2) = yo - delta/2. * cos(pi/36.)
                
                p2(1) = xo - delta/2. * sin(pi/36.)
                p2(2) = yo + delta/2. * cos(pi/36.)

                p3(1) = p2(1) - delta * 1./tan(alpha) * cos(pi/36.)
                p3(2) = p1(2) + delta / cos(pi/2. - alpha) * sin(alpha - pi/36.)


!                xt = xo
!                yt = yo  

                !-- p1 --p2
                h1 = abs(yt-p1(2)-tan(pi/2.+pi/36.)*(xt-p1(1)))/sqrt(1+(tan(pi/2.+pi/36.))**2 )
                !-- p1 --p3
                h2 = abs(yt-p1(2)-tan(pi + pi/36. - alpha)*(xt-p1(1)))/sqrt(1+(tan(pi + pi/36. - alpha))**2 )
                !-- p2 --p3
                h3 = abs(yt-p2(2)-tan(pi/36.)*(xt-p2(1)))/sqrt(1+(tan(pi/36.))**2)

                p1p2 = sqrt((p2(2)- p1(2))**2 + (p2(1)-p1(1))**2)
                p1p3 = sqrt((p3(2)- p1(2))**2 + (p3(1)-p1(1))**2)
                p2p3 = sqrt((p2(2)- p3(2))**2 + (p2(1)-p3(1))**2)

!                print*, h1, h2, h3


                if ((abs((h1*p1p2 + h2*p1p3 + h3*p2p3)- (delta*p2p3)) .le. 1e-12)) then
                  v(ii,jj) = 0.
                endif
!----triangle trailing edge
                p1(1) = xe - delta/2. * sin(pi/36.)
                p1(2) = ye + delta/2. * cos(pi/36.)

                p2(1) = xe + delta/2. * sin(pi/36.)
                p2(2) = ye - delta/2. * cos(pi/36.)

                p3(1) = p2(1) + delta * 1./tan(pi/3.) * cos(pi/36.)
                p3(2) = p1(2) - delta / cos(pi/6.) * sin(pi/3. - pi/36.)

!                xt = xo
!                yt = yo  

                !-- p1 --p2
                h1 = abs(yt-p1(2)-tan(pi/2.+pi/36.)*(xt-p1(1)))/sqrt(1+(tan(pi/2.+pi/36.))**2 )
                !-- p1 --p3
                h2 = abs(yt-p1(2)-tan(pi + pi/36. - pi/3.)*(xt-p1(1)))/sqrt(1+(tan(pi + pi/36. - pi/3.))**2 )
                !-- p2 --p3
                h3 = abs(yt-p2(2)-tan(pi/36.)*(xt-p2(1)))/sqrt(1+(tan(pi/36.))**2)

                p1p2 = sqrt((p2(2)- p1(2))**2 + (p2(1)-p1(1))**2)
                p1p3 = sqrt((p3(2)- p1(2))**2 + (p3(1)-p1(1))**2)
                p2p3 = sqrt((p2(2)- p3(2))**2 + (p2(1)-p3(1))**2)

!                print*, h1, h2, h3


                if ((abs((h1*p1p2 + h2*p1p3 + h3*p2p3)- (delta*p2p3)) .le. 1e-12)) then
                  v(ii,jj) = 0.
                endif
      endif



     if (phase .le. 0.) then
!----triangle leading edge
                p1(1) = xo + delta/2. * sin(pi/36.)
                p1(2) = yo + delta/2. * cos(pi/36.)

                p2(1) = xo - delta/2. * sin(pi/36.)
                p2(2) = yo - delta/2. * cos(pi/36.)

                p3(1) = p2(1) - delta * 1./tan(alpha) * cos(pi/36.)
                p3(2) = p1(2) - delta / cos(pi/2. - alpha) * sin(alpha - pi/36.)


!                xt = xo
!                yt = yo  

                !-- p1 --p2
                h1 = abs(yt-p1(2)-tan(pi/2. - pi/36.)*(xt-p1(1)))/sqrt(1+(tan(pi/2. - pi/36.))**2 )
                !-- p1 --p3
                h2 = abs(yt-p1(2)-tan(pi/2. - pi/36. - pi/2. + alpha)*(xt-p1(1)))/sqrt(1+(tan(pi/2. - pi/36. - pi/2. + alpha))**2 )
                !-- p2 --p3
                h3 = abs(yt-p2(2)-tan(pi - pi/36.)*(xt-p2(1)))/sqrt(1+(tan(pi - pi/36.))**2)

                p1p2 = sqrt((p2(2)- p1(2))**2 + (p2(1)-p1(1))**2)
                p1p3 = sqrt((p3(2)- p1(2))**2 + (p3(1)-p1(1))**2)
                p2p3 = sqrt((p2(2)- p3(2))**2 + (p2(1)-p3(1))**2)

!                print*, h1, h2, h3


                if ((abs((h1*p1p2 + h2*p1p3 + h3*p2p3)- (delta*p2p3)) .le. 1e-12)) then
                  v(ii,jj) = 0.
                endif
!----triangle trailing edge
                p1(1) = xe - delta/2. * sin(pi/36.)
                p1(2) = ye - delta/2. * cos(pi/36.)

                p2(1) = xe + delta/2. * sin(pi/36.)
                p2(2) = ye + delta/2. * cos(pi/36.)

                p3(1) = p2(1) + delta * 1./tan(pi/3.) * cos(pi/36.)
                p3(2) = p1(2) + delta / cos(pi/6.) * sin(pi/3. - pi/36.)

!                xt = xo
!                yt = yo  

                !-- p1 --p2
                h1 = abs(yt-p1(2)-tan(pi/2. - pi/36.)*(xt-p1(1)))/sqrt(1+(tan(pi/2. - pi/36.))**2 )
                !-- p1 --p3
                h2 = abs(yt-p1(2)-tan(pi/2. - pi/36. - pi/6.)*(xt-p1(1)))/sqrt(1+(tan(pi/2. - pi/36. - pi/6.))**2 )
                !-- p2 --p3
                h3 = abs(yt-p2(2)-tan(pi - pi/36.)*(xt-p2(1)))/sqrt(1+(tan(pi - pi/36.))**2)

                p1p2 = sqrt((p2(2)- p1(2))**2 + (p2(1)-p1(1))**2)
                p1p3 = sqrt((p3(2)- p1(2))**2 + (p3(1)-p1(1))**2)
                p2p3 = sqrt((p2(2)- p3(2))**2 + (p2(1)-p3(1))**2)

!                print*, h1, h2, h3


                if ((abs((h1*p1p2 + h2*p1p3 + h3*p2p3)- (delta*p2p3)) .le. 1e-12)) then
                  v(ii,jj) = 0.
                endif
      endif





             enddo
         enddo

      end subroutine oscillator_projection_vel_IBM_line_2D_design




      subroutine instability_filtering_setup(xo , yo, l,frac1,frac2, theta,u,v, u_average,v_average,n_time)


         integer             ::     ii,jj,kk



         integer             ::     x_index,  y_index

         integer             ::     n_filter                                  !   number of cells to be filtered
         integer, intent(in) ::     n_time

         real(8), intent(in)    ::     frac1, frac2                              !   start and end points of filter (in terms of leaflet length fraction)
         real(8), intent(in)    ::     theta
         real(8), intent(in)    ::     xo, yo,l                                  !   origin point and end point of the upper leaflet
         real(8), intent(inout) ::     u(N1p,N2p), v(N1p, N2p)                   !   velocity fields to be operated            
         real(8), intent(inout) ::     u_average(N1p,N2p), v_average(N1p, N2p)   !   Time averaged velocity fields
         real(8)                ::     xt, yt                                    !   test point
         real(8)                ::     xe, ye 
         real(8)                ::     m1, m2                                    !   slope of line and the normal
         real(8)                ::     g_dist                                    !   ghost cell distance from the line
         real(8)                ::     x1 , y1, x2, y2                           !   (x1,y1) is the location where the filter starts (x2,y2) is where it ends
         character*20        ::     file_name1

         real(8)                ::     r_s1, r_s2


         real(8),parameter      ::     pi = 3.14159265359 

! calculating computational geometry which has to be done outside of the loop


         xe = xo + l * cos(theta)
         ye = yo + l * sin(theta)


         x1 = xo + frac1 * l * cos(theta)
         x2 = xo + frac2 * l * cos(theta)
!for lower leaflet
         y1 = yo + frac1 * l * sin(theta)
         y2 = yo + frac2 * l * sin(theta)


         m1 = (ye - yo) / (xe -xo)

         m2 = -1. / m1

!             if (rank .eq. 1) then
!                print*, u(5,5), u_average(5,5)
!            endif


         
         r_s1 = abs(0.5 - y1)  ! y-centerline is considered at y=0.5 
!         r_s2 = abs(0.5 - yo - m1*(x2 - xo))/sqrt(m1**2 + 1.)  ! y-centerline is considered at y=0.5
 

         do ii = S11B, N11B
            do jj = S21B, N21B
! update u_average and v_average

!            u_average(ii,jj) = (n_time-1)*1.0/n_time * u_average(ii,jj) + u(ii,jj)/n_time

!            v_average(ii,jj) = (n_time-1)*1.0/n_time * v_average(ii,jj) + v(ii,jj)/n_time

! on the u grid
                xt = x1u(ii)
                yt = x2p(jj)
                if ((xt .ge. x1)  .and.  (xt .le. x2)) then
                    r_s2 = r_s1 - m1*(xt-x1)
                    !g_dist = abs(yt - yo - m1*(xt - xo))/sqrt(m1**2 + 1.)
                    !if (((g_dist) .ge. -1.5*dx1p(ii)) .and.  (g_dist .le. 0.0)) then
                    if (abs(yt-0.5) .le. r_s2) then 
                 !       print*, u_average(ii,jj)
                  !      u(ii,jj) = u_average(ii,jj)
                    endif
                endif
! on the v grid
                xt = x1p(ii)
                yt = x2v(jj)
                if ((xt .ge. x1)  .and.  (xt .le. x2)) then
                    r_s2 = r_s1 - m1*(xt-x1)
                    !g_dist = abs(yt - yo - m1*(xt - xo))/sqrt(m1**2 + 1.)
                    !if (((g_dist) .ge. -1.5*dx1p(ii)) .and.  (g_dist .le. 0.0)) then
                    if (abs(yt-0.5) .le. r_s2) then 
                    !    v(ii,jj) = v_average(ii,jj)
                         v(ii,jj) = 0 
                    endif
                endif

            enddo
        enddo

     end subroutine instability_filtering_setup

     subroutine thin_oscillator_projection_vel_IBM_line_2D(u, v, xo , yo, l, theta, angle_amplitude, freq, t)


         integer             ::     ii,jj,kk



         integer             ::     x_index,  y_index

         integer             ::     n_ghost                           !   number of ghost cells

         real(8), intent(in)    ::     theta, angle_amplitude, freq, t
         real(8)                ::     ghost_cells(2, 10000)               !   ghost cells
         real(8)                ::     image_cells(2, 10000)               !   image cells
         real(8), intent(in)    ::     xo, yo,l                            !   origin point and end point of the line
         real(8)                ::     xt, yt, xg, yg, xin, yin, xip, yip  !   test, ghost,intercept and image point coordinates 
         real(8)                ::     xe, ye
         real(8)                ::     m1, m2                              !   slope of line and the normal
         real(8)                ::     g_dist                              !   ghost cell distance from the line
         real(8)                ::     u(N1p,N2p),  v(N1p,N2p)
         real(8)                ::     x1 , y1, x2, y2, x3, y3, x4, y4
         integer             ::     Xindex, Yindex
         character*20        ::     file_name1

         real(8)                ::     radius

         real(8)                ::     phase, omega


         real(8),parameter      ::     pi = 3.14159265359 


! calculating computational geometry which has to be done outside of the loop

         phase = theta + angle_amplitude * sin (2. * pi * freq * t)

         xe = xo + l * cos(phase)
         ye = yo + l * sin(phase)


         m1 = (ye - yo) / (xe -xo)
         m2 = -1. / m1


         n_ghost = 0
         do ii = S11B, N11B
            do jj = S21B, N21B
! on the u grid
                xt = x1u(ii)
                yt = x2p(jj)

                if ((xt .ge. xo)  .and.  (xt .le. xe)) then
                    g_dist = (yt - yo - m1*(xt - xo))/sqrt(m1**2 + 1.)
                    if (abs(g_dist) .le. 1.*dx1p(ii)) then
                        xg = xt
                        yg = yt

                        xin = (yo - yg - m1*xo + m2*xg)/(m2 - m1)
                        yin = (yo + m1*(xin - xo))
                        
                        call FourPtFinder(3,xin,yin,x1p,x2p,x1u,x2v,x1,y1,x2,y2,x3,y3,x4,y4,dx1p,dx2p,dx1u,dx2v,N1p,N2p)
                        call Coordfinder(3,x1p,x2p,x1u,x2v,x2,y2,Xindex,Yindex,N1p,N2p)

                        !n_ghost = n_ghost + 1

                        radius = sqrt((xin-xo)**2 + (yin-yo)**2)

                        !call Coordfinder(3,x1p,x2p,x1u,x2v,xg,yg,Xindex,Yindex,N1p,N2p)
                        omega =  2. * pi * freq * angle_amplitude * cos(2. * pi * freq * t )


                        u(Xindex,Yindex) = -sign(1. , 1.*cos(2. * pi * freq * t)) * radius * omega *  sin(phase)  
                    endif
                endif

! on the v grid

                xt = x1p(ii)
                yt = x2v(jj)
                if ((xt .ge. xo)  .and.  (xt .le. xe)) then
                    g_dist = (yt - yo - m1*(xt - xo))/sqrt(m1**2 + 1.)

                    if (abs(g_dist) .le. 1.*dx1p(ii)) then

                      !  n_ghost = n_ghost+1

                        xg = xt
                        yg = yt

                        xin = (yo - yg - m1*xo + m2*xg)/(m2 - m1)
                        yin = (yo + m1*(xin - xo))
                       !wa 2 here
                        call FourPtFinder(2,xin,yin,x1p,x2p,x1u,x2v,x1,y1,x2,y2,x3,y3,x4,y4,dx1p,dx2p,dx1u,dx2v,N1p,N2p)
                        call Coordfinder(2,x1p,x2p,x1u,x2v,x2,y2,Xindex,Yindex,N1p,N2p)

                        radius = sqrt((xin-xo)**2 + (yin-yo)**2)
                       ! print*, radius
                        !call Coordfinder(3,x1p,x2p,x1u,x2v,xin,yin,Xindex,Yindex,N1p,N2p)

                     !   phase = theta + angle_amplitude * sin (2. * pi * freq * t)
                        omega =  2. * pi * freq * angle_amplitude * cos(2.*pi*freq*t)


                        v(Xindex,Yindex) = sign(1., 1.*cos(2. * pi * freq * t)) * radius * omega  * cos(phase)
                        !v(ii,jj) = sign(1., cos(2. * pi * freq * t)) * radius * 2. * pi * freq * cos(phase)

!                        call Coordfinder(2,x1p,x2p,x1u,x2v,xg,yg,Xindex,Yindex,N1p,N2p)

!                        v(ii,jj) = 0.

                    endif
                endif
            enddo
         enddo

      end subroutine thin_oscillator_projection_vel_IBM_line_2D




   subroutine thin_lid_fsi_find_pressure_sum(p, mass_rigid, xo , yo, l, theta, pressure_tot_of_process,rhs_thin)


         integer             ::     ii,jj,kk



         integer             ::     x_index,  y_index

         integer             ::     n_ghost                             !   number of ghost cells

         real(8), intent(in)    ::     theta
         real(8), intent(in)    ::     xo, yo,l, mass_rigid                            !   origin point and end point of the line
         real(8)                ::     xt, yt, xg, yg, xin, yin, xip, yip  !   test, ghost,intercept and image point coordinates 
         real(8)                ::     xe, ye
         real(8)                ::     m1, m2                              !   slope of line and the normal
         real(8)                ::     g_dist                              !   ghost cell distance from the line
         real(8),intent(in)     ::     p(N1p,N2p)
         integer             ::     Xindex, Yindex

         real(8)                ::     pip , pgp, deltaP                   !   pressure of ghost and image points 
         real(8)                ::     delP_tot

         real(8), intent(out)   ::     pressure_tot_of_process             !  

         integer             ::     rank_hinge                   
   
         integer             ::     point_exists 

         real(8)                ::     radius

         real(8),parameter      ::     pi = 3.14159265359 

         real(8)                ::     C1, C2, C3, C4                      !   bilinear coefficient interpolations
         real(8)                ::     x1, y1, x2, y2, x3, y3, x4, y4

         real(8)                ::     ksi, eta

         real(8)                ::     rhs_thin                            !   right hand side for the discrete rotation equation



! calculating computational geometry which has to be done outside of the loop

         call CheckInMPIProcess(2,xo,yo, point_exists,Xindex,Yindex)
         if (point_exists == 1) rank_hinge = rank

         xe = xo + l * cos(theta)
         ye = yo + l * sin(theta)

         m1 = (ye - yo) / (xe -xo)
         m2 = -1. / m1

         pressure_tot_of_process = 0.
         n_ghost = 0

         n_ghost = 0
         do ii = 1, N1p
            do jj = 1, N2p
! on the u grid
                xt = x1p(ii)
                yt = x2p(jj)

                if ((xt .ge. xo)  .and.  (xt .le. xe)) then
                    g_dist = (yt - yo - m1*(xt - xo))/sqrt(m1**2 + 1.)
                    if ((g_dist) .le. 1.*dx1p(ii) .and. (g_dist .ge. 0.) ) then
                        xg = xt
                        yg = yt

                        call CheckInMPIProcess(1,xg,yg, point_exists,Xindex,Yindex) 

                        if (point_exists == 1) then  

                            pgp = p(Xindex,Yindex)
                            
                            xin = (yo - yg - m1*xo + m2*xg)/(m2 - m1)
                            yin = (yo + m1*(xin - xo))

                            xip = 2.0 * xin - xt
                            yip = 2.0 * yin - yt
                        

                            n_ghost = n_ghost + 1
                            call CheckInMPIProcess(1,xip,yip, point_exists,Xindex,Yindex)
                            if (point_exists == 1) then
                            call find_local_coords(xip,yip, Xindex, Yindex, ksi, eta)
                            call bilinear_interpolate(p, xip, yip, Xindex, Yindex, ksi, eta, pip)
                               deltaP = pgp - pip
                            else
                               print*, 'WARNING : Image Point Out of Processor Domain'
                            endif

                        else

                            print*, 'WARNING : Ghost Point Out of Processor Domain'

                        endif

                   radius = sqrt((xin-xo)**2 + (yin-yo)**2)

                        !call Coordfinder(3,x1p,x2p,x1u,x2v,xg,yg,Xindex,Yindex,N1p,N2p)
                   pressure_tot_of_process = pressure_tot_of_process + radius * deltaP
                   rhs_thin = 2./(mass_rigid * n_ghost) * pressure_tot_of_process


                endif
              endif
           enddo
         enddo

      print*, n_ghost


      end subroutine thin_lid_fsi_find_pressure_sum


      subroutine thin_lid_velocity_mapping(u, v, xo , yo, l, theta, omega)


         integer             ::     ii,jj,kk



         integer             ::     x_index,  y_index

         integer             ::     n_ghost                           !   number of ghost cells

         real(8), intent(in)    ::     theta, omega
         real(8), intent(in)    ::     xo, yo,l                            !   origin point and end point of the line
         real(8)                ::     xt, yt, xg, yg, xin, yin, xip, yip  !   test, ghost,intercept and image point coordinates 
         real(8)                ::     xe, ye
         real(8)                ::     m1, m2                              !   slope of line and the normal
         real(8)                ::     g_dist                              !   ghost cell distance from the line
         real(8)                ::     u(N1p,N2p),  v(N1p,N2p)
         real(8)                ::     x1 , y1, x2, y2, x3, y3, x4, y4
         integer             ::     Xindex, Yindex
         character*20        ::     file_name1

         real(8)                ::     radius



         real(8),parameter      ::     pi = 3.14159265359 


! calculating computational geometry which has to be done outside of the loop


         xe = xo + l * cos(theta)
         ye = yo + l * sin(theta)


         m1 = (ye - yo) / (xe -xo)
         m2 = -1. / m1


         n_ghost = 0
         do ii = S11B, N11B
            do jj = S21B, N21B
! on the u grid
                xt = x1u(ii)
                yt = x2p(jj)

                if ((xt .ge. xo)  .and.  (xt .le. xe)) then
                    !g_dist = (yt - yo - m1*(xt - xo))/sqrt(m1**2 + 1.)
                    g_dist = abs(yt - yo - m1*(xt - xo))/sqrt(m1**2 + 1.)
                    !if (((g_dist) .ge. -1.5*dx1p(ii)) .and.  (g_dist .le. 0.0)) then
                    if (g_dist .le. 1.*dx1p(ii)) then 

                        radius = sqrt((xt-xo)**2 + (yt-yo)**2)

                      !  omega =  2. * pi * freq * angle_amplitude * cos(2.*pi*freq*t)

                        u(ii,jj) = -sign(1. , 1.*cos(omega)) * radius * omega *  sin(theta)  
                  
           !             u(ii,jj) = 0.
                     endif
                endif

! on the v grid

                xt = x1p(ii)
                yt = x2v(jj)
                if ((xt .ge. xo)  .and.  (xt .le. xe)) then
                    g_dist = abs(yt - yo - m1*(xt - xo))/sqrt(m1**2 + 1.)

                    if (g_dist .le. 1.*dx1p(ii)) then

                        radius = sqrt((xt-xo)**2 + (yt-yo)**2)


                   !     omega =  2. * pi * freq * angle_amplitude * cos(2.*pi*freq*t)
                    
                       ! v(ii,jj) = 0.                         
 
                        v(ii,jj) = sign(1., 1.*cos(omega)) * radius * omega  * cos(theta)
                    endif
                endif
            enddo
         enddo

      end subroutine thin_lid_velocity_mapping






 subroutine mesh_oscillator_projection_vel_IBM_line_2D(u, v, xo , yo, l, theta, angle_amplitude, freq, t)


         integer             ::     ii,jj,kk



         integer             ::     x_index,  y_index

         integer             ::     n_ghost                           !   number of ghost cells

         real(8), intent(in)    ::     theta, angle_amplitude, freq, t
         real(8), intent(in)    ::     xo, yo,l                            !   origin point and end point of the line
         real(8)                ::     xt, yt, xg, yg, xin, yin, xip, yip  !   test, ghost,intercept and image point coordinates 
         real(8)                ::     xe, ye
         real(8)                ::     m1, m2                              !   slope of line and the normal
         real(8)                ::     g_dist                              !   ghost cell distance from the line
         real(8)                ::     u(N1p,N2p),  v(N1p,N2p)
         real(8)                ::     x1 , y1, x2, y2, x3, y3, x4, y4
         integer             ::     Xindex, Yindex
         character*20        ::     file_name1

         real(8)                ::     radius  

         real(8)                ::     phase, omega
         integer             ::     point_exists

         real(8),parameter      ::     pi = 3.14159265359 
         integer, parameter  ::     mesh_size = 100
         real(8)                ::     mesh(2,mesh_size)

         real(8)                ::     x_start, y_start

         real(8)                ::     u_value
         real(8)                ::     v_value

         real(8)                ::     phi

         integer             ::     i , j

! calculating computational geometry which has to be done outside of the loop


         phase = theta + angle_amplitude * sin (2. * pi * freq * t)

         x_start = xo - 0.01*cos(phase)
         y_start = yo - 0.01*sin(phase)


         call  mesh_leaflet(x_start , y_start, l, phase, mesh_size, mesh)

!         print*, mesh
         xe = xo + l * cos(phase)
         ye = yo + l * sin(phase)


         m1 = (ye - yo) / (xe -xo)
         m2 = -1. / m1


         n_ghost = 0
         do ii =1,  mesh_size
! on the u grid
            xt = mesh(1,ii)
            yt = mesh(2,ii)
        
            call CheckInMPIProcess(2,xt,yt, point_exists,Xindex,Yindex)           
               if (point_exists == 1) then
                  radius = sqrt((xt-x_start)**2 + (yt-y_start)**2)
                  omega =  2. * pi * freq * angle_amplitude * cos(2. * pi * freq * t )
                  u_value = -sign(1. , 1.*cos(2. * pi * freq * t)) * radius * omega *  sin(phase) 
                  do i = -3, 3
                     do j = -3, 3
                        call peskin_phi_function(2, xt , yt , x1u(Xindex+i) , x1p(Yindex+j),  2.0,  phi)  
                        u(Xindex + i, Yindex + j) = phi * u_value
                     enddo
                  enddo
!                 print*, point_exists
               endif

            call CheckInMPIProcess(3,xt,yt, point_exists,Xindex,Yindex)

               if (point_exists == 1) then
                  radius = sqrt((xt-x_start)**2 + (yt-y_start)**2)
                  omega =  2. * pi * freq * angle_amplitude * cos(2. * pi * freq * t )
                  v_value = sign(1. , 1.*cos(2. * pi * freq * t)) * radius * omega *  cos(phase)

                  do i = -3, 3
                     do j = -3, 3
                        call peskin_phi_function(3, xt , yt , x1p(Xindex+i) , x2v(Yindex+j),  2.0,  phi)
                        u(Xindex + i, Yindex + j) = phi * u_value
                     enddo
                  enddo
               endif


         enddo

             
      end subroutine mesh_oscillator_projection_vel_IBM_line_2D



   subroutine mesh_leaflet(xo , yo, l, phase, mesh_size, mesh)


         integer             ::     ii,jj,kk


         integer             ::     n_ghost                           !   number of ghost cells
         integer, intent(in) ::     mesh_size                         !   mesh size

         real(8), intent(in)    ::     phase                             !   phasing of the leaflet
         real(8), intent(out)   ::     mesh(2, mesh_size)                !   mesh
         real(8), intent(in)    ::     xo, yo,l                          !   origin point and end point of the line

         real(8),parameter      ::     pi = 3.14159265359

         do ii=1,mesh_size
            mesh(1,ii) = xo + ii * (l/mesh_size) * cos(phase)
            mesh(2,ii) = yo + ii * (l/mesh_size) * sin(phase)
         enddo

   end subroutine mesh_leaflet

 subroutine spread_penalty_force_on_meshed_lid(u, v, xo , yo, l, theta, angle_amplitude, freq, t)


         integer             ::     ii,jj,kk



         integer             ::     x_index,  y_index

         integer             ::     n_ghost                           !   number of ghost cells

         real(8), intent(in)    ::     theta, angle_amplitude, freq, t
         real(8), intent(in)    ::     xo, yo,l                            !   origin point and end point of the line
         real(8)                ::     xt, yt, xg, yg, xin, yin, xip, yip  !   test, ghost,intercept and image point coordinates 
         real(8)                ::     xe, ye
         real(8)                ::     m1, m2                              !   slope of line and the normal
         real(8)                ::     g_dist                              !   ghost cell distance from the line
         real(8)                ::     u(N1p,N2p),  v(N1p,N2p)
         real(8)                ::     x1 , y1, x2, y2, x3, y3, x4, y4
         integer             ::     Xindex, Yindex
         character*20        ::     file_name1

         real(8)                ::     radius

         real(8)                ::     phase, omega
         integer             ::     point_exists

         real(8),parameter      ::     pi = 3.14159265359 
         integer, parameter  ::     mesh_size = 100
         real(8)                ::     mesh(2,mesh_size)

         real(8)                ::     x_start, y_start

! calculating computational geometry which has to be done outside of the loop


         phase = theta + angle_amplitude * sin (2. * pi * freq * t)

         x_start = xo - 0.01*cos(phase)
         y_start = yo - 0.01*sin(phase)


         call  mesh_leaflet(x_start , y_start, l, phase, mesh_size, mesh)

!         print*, mesh
         xe = xo + l * cos(phase)
         ye = yo + l * sin(phase)


         m1 = (ye - yo) / (xe -xo)
         m2 = -1. / m1


         n_ghost = 0
         do ii =1,  mesh_size
! on the u grid
            xt = mesh(1,ii)
            yt = mesh(2,ii)
        
            call CheckInMPIProcess(2,xt,yt, point_exists,Xindex,Yindex)           
               if (point_exists == 1) then
                  radius = sqrt((xt-x_start)**2 + (yt-y_start)**2)
                  omega =  2. * pi * freq * angle_amplitude * cos(2. * pi * freq * t )
                  u(Xindex+1,Yindex) = -sign(1. , 1.*cos(2. * pi * freq * t)) * radius * omega *  sin(phase) 
!                  print*, point_exists
               endif

            call CheckInMPIProcess(3,xt,yt, point_exists,Xindex,Yindex)

               if (point_exists == 1) then
                  radius = sqrt((xt-x_start)**2 + (yt-y_start)**2)
                  omega =  2. * pi * freq * angle_amplitude * cos(2. * pi * freq * t )
                  v(Xindex+1,Yindex) = sign(1. , 1.*cos(2. * pi * freq * t)) * radius * omega *  cos(phase)

               endif


         enddo

             
      end subroutine spread_penalty_force_on_meshed_lid





   subroutine FourPtFinder(flag,xip,yip,x,y,xh,yh,x1,y1,x2,y2,x3,y3,x4,y4,dx,dy,dxh,dyh,nxp2,nyp2)
!!     implicit double precision (a-h,o-z)
!!     integer i,j,ii,jj,nxp2,nyp2,nxp1,nyp1,flag
!!     dimension x(nxp2),y(nyp2),dx(nxp2),dy(nyp2)
!!     dimension xh(nxp2),dxh(nxp2),dyh(nyp2),yh(nyp2)
!!     common/grid/nxp1,nyp1,nxlast,nylast,nx,ny,xl,yl 


       integer,         intent(in)      ::     flag, nxp2, nyp2
       real(8),            intent(in)      ::     xip, yip
       real(8),            intent(in)      ::     x(nxp2), y(nyp2), xh(nxp2), yh(nyp2), dx(nxp2), dy(nyp2), dxh(nxp2), dyh(nyp2)
       real(8),            intent(out)     ::     x1, y1, x2, y2, x3, y3, x4, y4



       integer                          ::     i, j, ii, jj
       integer                          ::     nxp1, nyp1



!-  this is a subroutine for finding the four neighboring points of the image point
!- first we find the coordinates of the image point and then we comare it
!- by adding and deducing dx and dy we can find the coordinates of the points
!
!- flag=1 means the pressure grid
!- flag=2 means the u velocity grid
!- flag=3 means the v velocity grid
!- xi and yi are the coordinates of the four surounding points
!- xip and yip are the coordinates of the image point
!- x and y are the arrayes of the velocity grid locations 
!- xh and yh are the arrays of the pressure grid coordinates

!     print*, size(xh), size(yh), nxp2, nyp2, size(dx), size(dy), size(dxh), size(dyh) 



      nxp1=nxp2-1
      nyp1=nyp2-1


      if (flag .eq. 1) then
      do i=1,nxp1

          if ((xip .le. x(i+1)) .and. (xip .ge. x(i))) then

             x1=x(i)
             ii=i
          endif

      enddo
      do j=1,nyp1
          if ((yip .le. y(j+1)) .and. (yip .ge. y(j))) then
             y1=y(j)
             jj=j
          endif

      enddo


      x3=x1+dx(ii)
      y3=y1
      x2=x1
      y2=y1+dy(jj)
      x4=x1+dx(ii)
      y4=y1+dy(jj)
      endif
                     
      if (flag .eq. 2) then
      do i=1,nxp1
          if ((xip .le. xh(i+1)) .and. (xip .ge. xh(i))) then
             x1=xh(i)
             ii=i
          endif
      enddo
      do j=1,nyp1
          if ((yip .le. y(j+1)) .and. (yip .ge. y(j))) then
             y1=y(j)
             jj=j
          endif
      enddo
      x3=x1+dxh(ii)
      y3=y1
      x2=x1
      y2=y1+dy(jj)
      x4=x1+dxh(ii)
      y4=y1+dy(jj)
      endif



      if (flag .eq. 3) then
      do i=1,nxp1
          if ((xip .le. x(i+1)) .and. (xip .ge. x(i))) then
             x1=x(i)
             ii=i
          endif
      enddo
      do j=1,nyp1
          if ((yip .le. y(j+1)) .and. (yip .ge. y(j))) then
             y1=y(j)
             jj=j
          endif
      enddo
      x3=x1+dx(ii)
      y3=y1
      x2=x1
      y2=y1+dy(jj)
      x4=x1+dx(ii)
      y4=y1+dy(jj)
      endif


      end subroutine FourPtFinder
 
     subroutine CheckInMPIProcess(flag,x,y, point_exists,xindex,yindex)
!!     implicit double precision (a-h,o-z)
!!     integer i,j,ii,jj,nxp2,nyp2,nxp1,nyp1,flag
!!     dimension x(nxp2),y(nyp2),dx(nxp2),dy(nyp2)
!!     dimension xh(nxp2),dxh(nxp2),dyh(nyp2),yh(nyp2)
!!     common/grid/nxp1,nyp1,nxlast,nylast,nx,ny,xl,yl 


       integer,         intent(in)      ::     flag
       real(8),            intent(in)      ::     x, y
       integer,         intent(out)     ::     point_exists
       integer,         intent(out)     ::     xindex, yindex 

!
!  x. . . . .x
!  .         .
!  .         . 
!  .   .x,y  . 
!  .         . 
!  x . . . . x
!  |
!  --> this index is returned as xindex, yindex

       integer                          ::     i, j, ii, jj
       



!-  this is a subroutine for finding the four neighboring points of the image point
!- first we find the coordinates of the image point and then we comare it
!- by adding and deducing dx and dy we can find the coordinates of the points
!
!- flag=1 means the pressure grid
!- flag=2 means the u velocity grid
!- flag=3 means the v velocity grid
!- xi and yi are the coordinates of the four surounding points
!- xip and yip are the coordinates of the image point
!- x and y are the arrayes of the velocity grid locations 
!- xh and yh are the arrays of the pressure grid coordinates

!     print*, size(xh), size(yh), nxp2, nyp2, size(dx), size(dy), size(dxh), size(dyh) 


      point_exists = 0
      xindex = -5
      yindex = -5

      if (flag .eq. 1) then
      do i=1,N1p-1
        do j =1, N2p-1
          if ((x .le. x1p(i+1)) .and. (x .ge. x1p(i)) .and. ((y .le. x2p(j+1)) .and. (y .ge. x2p(j)))) then

             xindex = i
             yindex = j
             point_exists = 1 

          endif
        enddo
      enddo
      endif
                     
      if (flag .eq. 2) then
      do i=1,N1p-1
        do j =1, N2p-1
          if ((x .le. x1u(i+1)) .and. (x .ge. x1u(i)) .and. ((y .le. x2p(j+1)) .and. (y .ge. x2p(j)))) then

             xindex = i
             yindex = j
             point_exists = 1
      !       print*, 'what the fuck?'
          endif
        enddo
      enddo

      endif



      if (flag .eq. 3) then
      do i=1,N1p-1
        do j =1, N2p-1
          if ((x .le. x1p(i+1)) .and. (x .ge. x1p(i)) .and. ((y .le. x2v(j+1)) .and. (y .ge. x2v(j)))) then

             xindex = i
             yindex = j
             point_exists = 1

          endif
        enddo
      enddo
      endif


      end subroutine CheckInMPIProcess

      subroutine Coordfinder(flag,x,y,xh,yh,x1,y1,ii,jj,nxp2,nyp2)


       integer,   intent(in)    ::    nxp2, nyp2, flag
       real(8),      intent(in)    ::    x(nxp2), y(nyp2), xh(nxp2), yh(nyp2), x1, y1
       integer,   intent(out)   ::    ii, jj


       integer                  ::    i, j

!cc----Certified
!cc--  this subroutine take the coordinates and returns the indices
!cc--  flag=1 :pressure grid, flag=2: u velocity grid flag=3: v velocity grid

      if (flag .eq. 1) then
      do i=1,nxp2
         if (x1 .eq. x(i)) then
         ii=i
         endif
      enddo
      do j=1,nyp2
         if (y1 .eq. y(j)) then
         jj=j
         endif
      enddo
      endif

      if (flag .eq. 2) then
      do i=1,nxp2
         if (x1 .eq. xh(i)) then
         ii=i
         endif
      enddo
      do j=1,nyp2
         if (y1 .eq. y(j)) then
         jj=j
         endif
      enddo
      endif
      if (flag .eq. 3) then
      do i=1,nxp2
         if (x1 .eq. x(i)) then
         ii=i
         endif
      enddo
      do j=1,nyp2
         if (y1 .eq. yh(j)) then
         jj=j
         endif
      enddo
      endif
      end subroutine Coordfinder

       subroutine find_local_coords(x,y, Xindex, Yindex, ksi, eta)
!!     implicit double precision (a-h,o-z)
!!     integer i,j,ii,jj,nxp2,nyp2,nxp1,nyp1,flag
!!     dimension x(nxp2),y(nyp2),dx(nxp2),dy(nyp2)
!!     dimension xh(nxp2),dxh(nxp2),dyh(nyp2),yh(nyp2)
!!     common/grid/nxp1,nyp1,nxlast,nylast,nx,ny,xl,yl 


       real(8),            intent(in)      ::     x, y
       integer,         intent(in)      ::     Xindex, Yindex
       real(8),            intent(out)     ::     ksi, eta

!
!  x. . . . .x
!  .         .
!  .         . 
!  .   .x,y  . 
!  .         . 
!  x . . . . x
!  |
!  --> this index is returned as xindex, yindex




!-  this is a subroutine for finding the four neighboring points of the image point
!- first we find the coordinates of the image point and then we comare it
!- by adding and deducing dx and dy we can find the coordinates of the points
!
!- flag=1 means the pressure grid
!- flag=2 means the u velocity grid ! not here YET
!- flag=3 means the v velocity grid ! not here YET
!- xi and yi are the coordinates of the four surounding points
!- xip and yip are the coordinates of the image point
!- x and y are the arrayes of the velocity grid locations 
!- xh and yh are the arrays of the pressure grid coordinates

!     print*, size(xh), size(yh), nxp2, nyp2, size(dx), size(dy), size(dxh), size(dyh) 


!      if (flag .eq. 1) then
!      do i=1,N1p-1
!        do j =1, N2p-1
!          if ((x .le. x1p(i+1)) .and. (x .ge. x1p(i)) .and. ((y .le. x2p(j+1)) .and. (y .ge. x2p(j)))) then
!
!             xindex = i
!             yindex = j
!             point_exists = 1
!
!          endif
!        enddo
!      enddo
!      endif




    ksi = -1. + (x-x1p(Xindex))/dx1p(Xindex)
    eta  = -1. + (y-x2p(Yindex))/dx2p(Yindex) 


    end subroutine find_local_coords

     subroutine bilinear_interpolate(p,x,y, Xindex, Yindex, ksi, eta, evaluated)

!      interpolates with a bilinear shape function

!!     f77 header

!!     implicit double precision (a-h,o-z)
!!     integer i,j,ii,jj,nxp2,nyp2,nxp1,nyp1,flag
!!     dimension x(nxp2),y(nyp2),dx(nxp2),dy(nyp2)
!!     dimension xh(nxp2),dxh(nxp2),dyh(nyp2),yh(nyp2)
!!     common/grid/nxp1,nyp1,nxlast,nylast,nx,ny,xl,yl 


       real(8),            intent(in)      ::     x, y
       integer,         intent(in)      ::     Xindex, Yindex
       real(8),            intent(in)      ::     ksi, eta
       real(8),            intent(in)      ::     p(N1p,N2p)
       real(8)                             ::     phi1, phi2, phi3, phi4    ! bilinear shape function coefficients

       real(8),            intent(out)     ::     evaluated

!
!  x. . . . .x
!  .         .
!  .         . 
!  .   .x,y  . 
!  .         . 
!  x . . . . x
!  |
!  --> this index is returned as xindex, yindex




       !--- calculate the shape function (phi) components

       phi1 = 1./4. * (1. - ksi) * (1. - eta)
       phi2 = 1./4. * (1. + ksi) * (1. - eta)
       phi3 = 1./4. * (1. - ksi) * (1. + eta)
       phi4 = 1./4. * (1. + ksi) * (1. + eta)


       !--- calculate the dot product phi (transpose) . p

       evaluated = phi1 * p(Xindex,Yindex) + phi2 * p(Xindex+1, Yindex) + phi3 * p(Xindex, Yindex+1) + phi4 * p(Xindex+1, Yindex+1)


    end subroutine bilinear_interpolate





     subroutine advance_theta_in_time(theta1, theta2, delta_t, rhs, theta, omega)

!      does first order backward Euler to advance theta in time


       real(8),            intent(in)      ::    rhs, delta_t
       
       real(8),            intent(out)     ::    theta, omega

       real(8),            intent(inout)   ::    theta1, theta2

       theta =  rhs * delta_t**2 + 2. * theta1 - theta2 
       omega = (theta - theta2)/delta_t

       ! renew the  

       theta1 = theta2
       theta2 = theta
       

    end subroutine advance_theta_in_time






    subroutine peskin_phi_function(flag , xo , yo , x , y, radi,  phi)

!      does first order backward Euler to advance theta in time


       integer,         intent(in)    ::     flag

       real(8),            intent(in)    ::     xo , yo , x, y , radi

       real(8),            intent(out)   ::     phi 

       integer                        ::     i,j  

       integer                        ::     Xindex, Yindex
   
       integer                        ::     point_exists

       real(8)                           ::     dist, r
       
       real(8), parameter                ::     pi = 3.14159265359
!      --------------------------------------------------------

!flag = 1 : pressure , flag = 2 : u velocity , flag = 3 : v velocity 

       if  (flag == 1) then
            phi = 0.0
            call CheckInMPIProcess(flag,xo,yo, point_exists,Xindex,Yindex)
            if (point_exists == 1) then
                dist = sqrt((x-xo)**2 + (y-yo)**2) 
                r = dist/dx1p(Xindex)
                if (r .le. 2.) then
                    phi = 0.25 * (1 + cos (pi * r / 2.))
                  
                endif                       
            endif
       endif
 
       if  (flag == 2) then
            phi = 0.0
            call CheckInMPIProcess(flag,xo,yo, point_exists,Xindex,Yindex)
            if (point_exists == 1) then
                dist = sqrt((x-xo)**2 + (y-yo)**2)
                r = dist/dx1p(Xindex)
                if (r .le. 2.) then
                    phi = 0.25 * (1 + cos (pi * r / 2.))
                endif
            endif
       endif

       if  (flag == 3) then
            phi = 0.0
            call CheckInMPIProcess(flag,xo,yo, point_exists,Xindex,Yindex)
            if (point_exists == 1) then
                dist = sqrt((x-xo)**2 + (y-yo)**2)
                r = dist/dx1p(Xindex)
                if (r .le. 2.) then
                    phi = 0.25 * (1 + cos (pi * r / 2.))
                endif
            endif
       endif
  
    end subroutine peskin_phi_function


 end module mod_sharp_interface_ibm
