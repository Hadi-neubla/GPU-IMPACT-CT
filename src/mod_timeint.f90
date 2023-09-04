!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!* modified by Barna Becsek, ARTORG CVE, University of Bern (barna.becsek@artorg.unibe.ch)                   *
!* October 2014                                                                                              *
!* GPU version by Hadi Zolfaghari , ARTORG CVE, then DAMTP, Cambridge University (hz382@cam.ac.uk)           *
!* Oct 2015 - Sep 2023                                                                                       *
!*************************************************************************************************************


!> module containing the subroutine timeintegration, which is basically the central part of the DNS.
!! It uses modules mod_dims, mod_vars, mod_exchange, mod_diff, mod_inout, mod_coeffs, mod_lib, mod_test, mod_particles,
!! mod_solvers, mod_rhs
MODULE mod_timeint
  
  
  USE mod_dims
  USE mod_vars
  USE mod_exchange
  USE mod_diff
  USE mod_inout
  USE mod_coeffs
  USE mod_lib
  USE mod_test
  USE mod_solvers
  USE usr_func
  USE mod_rhs
  !USE mpi

  PRIVATE
  
  PUBLIC timeintegration
  
  INCLUDE 'mpif.h'
  
  CONTAINS
  
!pgi$g unroll = n:8
!!pgi$r unroll = n:8
!!pgi$l unroll = n:8
  
  
  !> public subroutine that forms the center of the DNS simulation.
  SUBROUTINE timeintegration
  
  IMPLICIT NONE
  
  INTEGER                ::  m

  !---Header variables for sharp interface immersed boundary method
  !---Hadi Zolfaghari 

  REAL(8)                   ::          ycurve(N1p),yhcurve(N1p)
  REAL(8)                   ::          ycurve_down(N1p),yhcurve_down(N1p)
  integer                ::          Inner_V_i(200000), Inner_V_j(200000)
  integer                ::          Inner_U_i(200000), Inner_U_j(200000)
  REAL(8)                   ::          u(N1p,N2p),  v(N1p,N2p), p(N1p,N2p)
  REAL(8)                   ::          u_average(N1p,N2p), v_average(N1p, N2p)
  integer                ::          ii, jj, kk, iii
  integer                ::          InVsize, InUsize
  integer, parameter     ::          sharp_interface_ibm_Rami2D_on = 0, MHV_3D_on=0, sphere_sim =0, TGV_sim=0

         !---leaflet  and commisure characteristics in 2D
  REAL(8)                   ::          freq, x1_oscil, y1_oscil, l_oscil, theta_axis_oscil, angular_amplitude,alpha
  REAL(8)                   ::          x1_commis, y1_commis, x2_commis, y2_commis
  REAL(8), parameter        ::          pi = 3.14159265359

  REAL(8)                   ::          pressure_tot_of_process
  REAL(8)                   ::          mass_rigid
  REAL(8)                   ::          rhs_thin , omega_thin, theta_thin 
  REAL(8)                   ::          theta1, theta2  
  REAL(8)                   ::          t1, t2
  REAL(8)                   ::          x_c,r_c,x_s, r_s   ! the hinge point on the lower side  
  REAL(8)                   ::          x_cont, y_cont, r_cont
         !--- TGV benchmark
  REAL(8)                   ::          E_new, E_new_global, E_old, E_old_global, KE_diss
!  REAL(8)                   ::          pre_inlet, q_inlet, pre_outlet, q_outlet
  integer                ::          rank_inlet, rank_inlet_global, rank_outlet, rank_outlet_global, blocks_count, req1_in, req2_out
  INTEGER                ::          rank_outlet_2, rank_outlet_2_global
  LOGICAL                ::          inlet_checked, outlet_checked
  INTEGER                ::          status(MPI_STATUS_SIZE)

  !--- Startvektor -------------------------------------------------------------------------------------------
  ! Hier bereits notwendig, um weitere Konzentrationen zuschalten zu koennen (BC werden aber weiter unten initialisiert!).
                         CALL initial_conditions_vel
  
  
  IF (restart == 0) THEN
     time          = time_start
     time_out_vect = time_start
     time_out_scal = time_start
     
     dtime         = 0.
     timestep      = 0
     
     new_dtime      = .TRUE.
     write_out_vect = .TRUE.
     write_out_scal = .TRUE.
     
     write_count = 0
     
     IF (dtime_out_vect == 0.) write_out_vect = .FALSE.
     IF (dtime_out_scal == 0.) write_out_scal = .FALSE.
  ELSE
                               CALL read_restart
!     IF (IB_on .AND. rank==0 ) CALL read_restart_ibm !bbecsek
     IF (dtime_out_scal /= 0.) CALL read_restart_stats
     IF (rank == 0 .AND. write_stout_yes) THEN
        WRITE(*,'(a,1E13.5)') '             time =', time
        WRITE(*,'(a,1i5   )') '         timestep =', timestep
        
        WRITE(*,'(a,1E13.5)') '            dtime =', dtime
        WRITE(*,'(a,1E13.5)') '    time_out_vect =', time_out_vect
        WRITE(*,'(a,1E13.5)') '    time_out_scal =', time_out_scal
        
        WRITE(*,'(a,1L2   )') '        new_dtime =', new_dtime
        WRITE(*,'(a,1L2   )') '   write_out_vect =', write_out_vect
        WRITE(*,'(a,1L2   )') '   write_out_scal =', write_out_scal
        
        WRITE(*,'(a,1i5   )') '      write_count =', write_count
        print*, p_end_systole_in, p_end_systole_out, 'p end systole in and out' 
        print*, 'this is done'
     END IF
  END IF
  
  timestep_old  = timestep
  dtime_old     = dtime
  dtime_average = 0.
  finish_yes    = .FALSE.
  
  !--- Null-Raeume bestimmen ---------------------------------------------------------------------------------
  ! Steht hier, weil Korrekturvektor "th" nach "configuration" erst alloziert und bestimmt werden muss
  ! ("initial_conditions_" werden danach als nächstes gerufen, s.o.)
  ! Alternative: eigene Subroutine für "th" kreieren ...
  IF (nullspace_yes) THEN
     CALL get_stencil_transp
!     CALL get_nullspace
     CALL get_stencil ! TEST!!! Unschoen! besser zu impact.f90 packen und Reihenfolge abaendern ...
  END IF
  
  !--- RB initialisieren -------------------------------------------------------------------------------------
  CALL init_BC
  
  !--- Divergenz-Freiheit testen -----------------------------------------------------------------------------
  CALL test_divergence
  
  !--- diverse Files öffnen ----------------------------------------------------------------------------------
  CALL open_stats
  
  !--- File fuer zwischenzeitlichen Abbruch der Zeitintegration neu erstellen --------------------------------
  IF (rank == 0) THEN
     OPEN (10,FILE='send_signal.txt',STATUS='UNKNOWN')
     WRITE(10,'(a)') '0'
     WRITE(10,*) n_timesteps
     WRITE(10,*) time_end
     CALL flush(10)
     CLOSE(10)
  END IF
  
  
  IF (rank == 0 .AND. write_stout_yes) THEN
     WRITE(*,'(a)')
     WRITE(*,'(a)') '================================================================================='
     WRITE(*,'(a)') '================================================================================='
     WRITE(*,'(a)') '================================================================================='
     WRITE(*,'(a)')
     WRITE(*,'(a)') '---------------------------- START TIME-INTEGRATION -----------------------------'
     WRITE(*,'(a)')
     WRITE(*,'(a,E12.5)')                 '                     Re =',Re
     WRITE(*,'(a,i4,a,i4,a,i4)')          '          box resolution:',M1,' x',M2,' x',M3
     WRITE(*,'(a,E12.5,a,E12.5,a,E12.5)') '          box dimension :',L1,' x',L2,' x',L3
     WRITE(*,'(a,E12.5)')                 '                   epsU =',epsU
     WRITE(*,'(a)') '================================================================================='
     WRITE(*,'(a)') '================================================================================='
     WRITE(*,'(a)') '================================================================================='
     WRITE(*,'(a)')
  END IF
  
  
  !--- Zeitmessung starten -----------------------------------------------------------------------------------
  IF (rank == 0) THEN
     CALL DATE_AND_TIME(values=ctime)
     day  = ctime(3)
     hour = ctime(5)
     minu = ctime(6)
     sec  = ctime(7)
     msec = ctime(8)
     
     elatime = msec+1000*(sec+60*(minu+60*hour))
     OPEN(99,FILE='test_wallclocktime_restart'//restart_char//'.txt',STATUS='UNKNOWN')
     WRITE(99,'(a,i2,a,i2,a,i4,a,i2,a,i2,a,i2,a,i3)') 'Begin time integration at ', ctime(3),'.',ctime(2),'.',ctime(1),    &
                                                           &      ', ',ctime(5),':',ctime(6),':',ctime(7),'.',ctime(8)
     CALL flush(99)
  END IF
  
  
  !--- Ausschreiben ------------------------------------------------------------------------------------------
  IF (write_xdmf_yes .AND. write_out_vect) CALL write_xdmf_xml ! bbecsek
  IF (write_out_scal) CALL compute_stats
  IF (write_out_vect) CALL write_fields

! -- for performance runs comment the line above and uncomment the two lines below
  time_out_vect  = time_out_vect + dtime_out_vect
  write_out_vect = .FALSE.

  !===========================================================================================================


  !--- initialization of the rigid geometry for the sharp interface ibm
  !--- Hadi Zolfaghari
     if ( (sharp_interface_ibm_Rami2D_on == 1) ) then

!      call parallel_shape_function(ycurve,yhcurve, ycurve_down, yhcurve_down)
!      call parallel_save(InUsize, InVsize, Inner_V_i, Inner_V_j, Inner_U_i, Inner_U_j,ycurve, yhcurve, ycurve_down, yhcurve_down)

     endif



  IF (TGV_sim == 1) THEN
  !--- intialize kinetic energy  (for TGV benchmarking)
      E_old=0.
      E_old_global=0.
      do ii=S1p, N1p
         do jj=S2p, N2p
            do kk=S3p, N3p
                E_old=E_old+(vel(ii,jj,kk,1)**2+vel(ii,jj,kk,2)**2+vel(ii,jj,kk,3)**2)/(2.*(M1-1)**3)
            enddo
         enddo
      enddo

   CALL MPI_REDUCE(E_old,E_old_global,1,MPI_REAL8,MPI_SUM,0,COMM_CART,merror)

   !--- end of initialize kinetic energy (for TGV benchmarking)
  ENDIF 
  

  !===========================================================================================================
  !=== Zeitintegration =======================================================================================
  !===========================================================================================================
  timeint: DO

    ! CALL test_divergence
     
     CALL get_dtime

     

     IF (rank == 0 .AND. write_stout_yes) THEN
        IF (timeint_mode == 0) THEN
           WRITE(*,'(a)')
           WRITE(*,'(a)') '================================================================================='
           WRITE(*,'(a)') '================================================================================='
        END IF
        WRITE(*,'(a,i8,a,E25.17,a,E25.17)') 'time step = ',timestep,' ; time =',time,' ; dtime =',dtime
        IF (timeint_mode == 0) THEN
           WRITE(*,'(a)') '================================================================================='
           WRITE(*,'(a)') '================================================================================='
        END IF
     END IF
     
     IF (rank == 0 .AND. log_iteration_yes) THEN
        OPEN(10,FILE='log_iterations.txt', STATUS='UNKNOWN')
        WRITE(10,'(a)')
        WRITE(10,'(a)') '================================================================================='
        WRITE(10,'(a)') '================================================================================='
        WRITE(10,'(a,i8,a,E25.17,a,E25.17)') 'time step = ',timestep,'; time =',time,'; dtime =',dtime
        WRITE(10,'(a)') '================================================================================='
        WRITE(10,'(a)') '================================================================================='
     END IF
     
     !========================================================================================================
     
     DO substep = 1, RK_steps
        
        IF (rank == 0 .AND. write_stout_yes .AND. timeint_mode == 0) THEN
           WRITE(*,'(a)') 
           WRITE(*,'(a)') '================================================================================='
           WRITE(*,'(a,i2,a)') 'Runge-Kutta sub-step',substep,':'
           WRITE(*,'(a)') '================================================================================='
        END IF
        IF (rank == 0 .AND. log_iteration_yes) THEN
           WRITE(10,'(a)') 
           WRITE(10,'(a)') '================================================================================='
           WRITE(10,'(a,i2,a)') 'Runge-Kutta sub-step',substep,':'
           WRITE(10,'(a)') '================================================================================='
        END IF
        
        
        !--- Zeit --------------------------------------------------------------------------------------------
        IF (substep == 1) subtime = time + dtime* aRK(1)
        IF (substep == 2) subtime = time + dtime*(aRK(1)+aRK(2)+bRK(2))
        IF (substep == 3) subtime = time + dtime
        
        
        !--- ghost cell update (fuer RHS) --------------------------------------------------------------------
        CALL exchange_all_all(.TRUE.,vel)

        
        
        !--- interpolate advection velocity + update ghost cells ---------------------------------------------
        ! vel(:,:,:,i) --> worki(:,:,:)
        CALL interpolate_vel(.FALSE.) ! TEST!!! Wurde teilweise schon bei Zeitschritt-Bestimmung erledigt!
        
        !--- rhs (ggf. Neumann-RB überschreiben) -------------------------------------------------------------
        ! Muss vor Konzentrationen kommen, weil
        !  - bcii für die no-flux-RB verwendet werden,
        !  - die Konzentrationen die Eddy-Viscosity des Geschwindigkeitsfeldes benoetigen.
        CALL rhs_vel
        
        
        !--- Helmholtz-Multiplikator -------------------------------------------------------------------------
        multL = thetaL*(aRK(substep)+bRK(substep))*dtime / Re
        
        
        !--- Umskalieren (Effizienz, initial guess) ----------------------------------------------------------
        IF (.NOT. init_pre(substep)) pre(S1p:N1p,S2p:N2p,S3p:N3p) = pre(S1p:N1p,S2p:N2p,S3p:N3p) * (aRK(substep)+bRK(substep)) * dtime
        
        !CALL apply_pressure_BC_CT


!   if ((ct_based_aorta == 1) .or. (idealized_aorta == 1)) then
!
!      do ii=s1p,n1p
!         do jj =s2p, n2p
!             do kk =s3p, n3p
!                if  (ct_grid(ii+ishift,jj+jshift,kk + kshift) == 5.0) then        
!                     !nl(ii,jj,kk,3) =  nl(ii,jj,kk,3) - 1. *  pre_inlet/(y3p(k_start_index-1) - y3p(k_start_index-10)) 
!                     !pre(ii,jj,kk) =   - (y3p(k_start_index-10) - x3p(kk))  *  pre_inlet/(y3p(k_start_index) - y3p(k_start_index-10)) 
!                     pre(ii,jj,kk) = pre_inlet
!                else if  (ct_grid(ii+ishift,jj+jshift,kk + kshift) == 9.0) then 
!                  !   nl(ii,jj,kk,3) =  nl(ii,jj,kk,3) - 0.8 *  pre_outlet/(y3p(k_start_index-1) - y3p(k_start_index-10)) 
!                     !pre(ii,jj,kk) =   -(y3p(k_start_index-10) - x3p(kk)) *  pre_outlet/(y3p(k_start_index) - y3p(k_start_index-10)) 
!                     pre(ii,jj,kk) = pre_outlet
!                else if  (ct_grid(ii+ishift,jj+jshift,kk + kshift) == 11.0) then 
!                  !   nl(ii,jj,kk,3) =  nl(ii,jj,kk,3) - 0.8 *  pre_outlet_2/(y3p(k_start_index_2+1) - y3p(k_start_index_2+10)) 
!                     !pre(ii,jj,kk) =   -(y3p(k_start_index_2+10) - x3p(kk)) *  pre_outlet_2/(y3p(k_start_index_2) - y3p(k_start_index_2+10)) 
!                     pre(ii,jj,kk) = pre_outlet_2
!                endif
!             end do
!          end do        
!      end do
!
!   endif






        !--- Löser -------------------------------------------------------------------------------------------
        IF (timeint_mode == 1 .OR. thetaL == 1.) THEN
!     call MPI_Barrier(MPI_COMM_WORLD, merror)
!     t1= MPI_Wtime()
     CALL explicit
!    call MPI_Barrier(MPI_COMM_WORLD, merror)
!     t2= MPI_Wtime()


!     print*, 'explicit substep', t2-t1

    ELSE
           IF (twostep_yes) THEN
        !      CALL twostep
           ELSE
         !     CALL outer_iteration
           END IF
        END IF
        
        !--- physikalischer Druck ----------------------------------------------------------------------------
        pre(S1p:N1p,S2p:N2p,S3p:N3p) = pre(S1p:N1p,S2p:N2p,S3p:N3p) / (aRK(substep)+bRK(substep)) / dtime
        
        !--- Undefinierte Ecken / Kanten auffüllen -----------------------------------------------------------
        CALL fill_corners(pre)
        
     END DO


    IF (TGV_sim == 1) THEN
    !--- Energy statistics for TGV benchmark    
        E_new=0.
        IF (rank==0)   OPEN(UNIT=139,FILE="KE_decay.txt",FORM="FORMATTED",position="append",STATUS="OLD",ACTION="READWRITE")
        IF (rank==0)   OPEN(UNIT=141,FILE="KE.txt",FORM="FORMATTED",position="append",STATUS="OLD",ACTION="READWRITE")
        DO ii=S1p, N1p
           DO jj=S2p, N2p
              DO kk=S3p, N3p
                  E_new=E_new+(vel(ii,jj,kk,1)**2+vel(ii,jj,kk,2)**2+vel(ii,jj,kk,3)**2)/(2.*(M1-1)**3)
              ENDDO
           ENDDO
        ENDDO
    
       CALL MPI_REDUCE(E_new,E_new_global,1,MPI_REAL8,MPI_SUM,0,COMM_CART,merror)
    
       IF ((rank==0)) THEN
          KE_diss = (E_old_global-E_new_global)/dtime
          E_old_global=E_new_global
          WRITE(UNIT=139, FMT=*) KE_diss
          WRITE(UNIT=141, FMT=*) E_new_global
       ENDIF

   !--- End of energy statistics for TGV benchmark
   ENDIF
     !========================================================================================================
     timestep = timestep + 1
     time     = time + dtime
     

     !--- Applying the sharp-interface method
     !--- Hadi Zolfaghari
     !--- TODO port to GPU


!     pre(N1, :, 1) = 0.0

      call apply_CT_forcing

!!--- flow forcing for CT-based DNS case
!     IF (CT_based_aorta == 1) THEN
!        flow_rate_correction_factor = 0.33  ! --- should be one for no constriction
!        call get_pulse_number_and_length
!        shifted_time_in_pulse = time - (pulse_number - 1) * pulse_length
!        if (rank == 0) print*, pulse_number, time, shifted_time_in_pulse
!        do ii = S11B, N11B
!            do jj = S21B,N21B
!               do kk = S31B, N31B
!                  if ((ct_grid(ii+iShift, jj+jShift, kk+kShift) .lt. 0.5) .or. (ct_grid(ii+iShift, jj+jShift, kk+kShift) ==  7.0)) then
!                      vel(ii,jj,kk,1) = .0
!                      vel(ii,jj,kk,2) = .0
!                      vel(ii,jj,kk,3) = .0
!                  elseif (ct_grid(ii+iShift, jj+jShift, kk+kShift) == 5.0)  then
!                      do iii = 1, number_of_digitized_points
!                         if ((shifted_time_in_pulse .ge. waveform_data_converted(iii)%time_in_pulse)  &
!                             & .and. (shifted_time_in_pulse .le. waveform_data_converted(iii+1)%time_in_pulse)) then
!                            vel(ii,jj,kk,3) = waveform_data_converted(iii)%Q_in_pulse * flow_rate_correction_factor
!                            vel(ii,jj,kk,2) = 0.
!                            vel(ii,jj,kk,1) = 0.
!                         endif
!                      enddo
!                  endif
!               enddo
!            enddo
!         enddo
!      ENDIF
!
!
!     inlet_checked = .FALSE.
!     IF (idealized_aorta == 1) THEN
!        !flow_rate_correction_factor = 0.33  ! --- should be one for no constriction
!        !call get_pulse_number_and_length
!        !shifted_time_in_pulse = time - (pulse_number - 1) * pulse_length
!        !if (rank == 0) print*, pulse_number, time, shifted_time_in_pulse
!        DO ii = S11B, N11B
!             DO jj = S21B,N21B
!                 DO kk = S31B, N31B
!                     IF ((ct_grid(ii+iShift, jj+jShift, kk+kShift) .lt. 0.5) .or. (ct_grid(ii+iShift, jj+jShift, kk+kShift) ==  7.0)) then
!                            vel(ii,jj,kk,1) = .0
!                            vel(ii,jj,kk,2) = .0
!                            vel(ii,jj,kk,3) = .0
!                     ELSEIF (ct_grid(ii+iShift, jj+jShift, kk+kShift) == 5.0)  then                     
!                            CALL get_pressure_and_flow_rate_WKM(1,pre_inlet, q_inlet)
!                            rank_inlet = rank
!                            blocks_count = NB1 * NB2 * NB3
!                            !pre_slope_inlet = pre_inlet/(y3p(k_start_index-1) - y3p(k_start_index-10))
!
!                            !IF (inlet_checked == .FALSE.) THEN 
!                            !   IF (rank .ne. 0) THEN
!                            !      CALL MPI_SEND(p_end_systole_in,1,MPI_REAL8,0,1,COMM_CART,req1_in,merror)
!                            !   ENDIF
!                            !ENDIF
!                            !inlet_checked = .TRUE.    
!                            !blocks_count = NB1 * NB2 * NB3
!                            !DO iii = 1,  blocks_count
!                            !   IF iii == rank_init 
!
!                            !IF (rank == 0) print*, 'passed'
!                            !if (rank == 0) print*, p_end_systole_in
!                            vel(ii,jj,kk,3) = q_inlet/(inlet_voxels_count * (L1/(M1-1))**2)
!                            vel(ii,jj,kk,2) = 0.
!                            vel(ii,jj,kk,1) = 0.
!                            !pre(ii,jj,kk) = pre_inlet
!                     ELSEIF (ct_grid(ii+iShift, jj+jShift, kk+kShift) == 9.0)  then                     
!                            CALL get_pressure_and_flow_rate_WKM(2,pre_outlet, q_outlet)
!                            rank_outlet = rank
!                            !pre_slope_outlet = pre_outlet/(y3p(k_start_index-1) - y3p(k_start_index-10))
!                            !CALL MPI_BCAST(p_end_systole_out,1,MPI_REAL8,rank,COMM_CART,merror)
!                            !IF (rank == 4) print*, 'passed 2'
!                            !vel(ii,jj,kk,3) = -q_outlet/(inlet_voxels_count * (L1/(M1-1))**2) !TODO change for outlet_voxel_count 
!                            !vel(ii,jj,kk,2) = 0.
!                            !vel(ii,jj,kk,1) = 0.
!                            !pre(ii,jj,kk) = pre_outlet
!                     ELSEIF (ct_grid(ii+iShift, jj+jShift, kk+kShift) == 11.0)  then                     
!                            CALL get_pressure_and_flow_rate_WKM(3,pre_outlet_2, q_outlet_2)
!                            rank_outlet_2 = rank
!                            !pre_slope_outlet_2 = pre_outlet_2/(y3p(k_start_index_2+1) - y3p(k_start_index_2+10))
!                           !print*, rank
!                     ELSEIF (ct_grid(ii+iShift, jj+jShift, kk+kShift) == 2.0) THEN
!                            pre(ii,jj,kk) = 0.
!                     ENDIF
!                  ENDDO
!               ENDDO
!            ENDDO
!        ENDIF
!
!       CALL MPI_ALLREDUCE(rank_inlet,rank_inlet_global,1,MPI_INT,MPI_MAX,COMM_CART,merror)    
!       CALL MPI_ALLREDUCE(rank_outlet,rank_outlet_global,1,MPI_INT,MPI_MAX,COMM_CART,merror)    
!       CALL MPI_ALLREDUCE(rank_outlet_2,rank_outlet_2_global,1,MPI_INT,MPI_MAX,COMM_CART,merror)    
!       
!       rank_outlet   = rank_outlet_global
!       rank_outlet_2 = rank_outlet_2_global
!       rank_inlet    = rank_inlet_global
!       
!       !CALL MPI_IRECV(p_end_systole_in,1,MPI_REAL8,rank,1,COMM_CART,req1_in,merror)
!       !CALL MPI_WAIT(req1_in(rank), status, merror)
!
!       CALL MPI_BCAST(p_end_systole_in,1,MPI_REAL8,rank_inlet,COMM_CART,merror)
!       CALL MPI_BCAST(p_end_systole_out,1,MPI_REAL8,rank_outlet,COMM_CART,merror)
!       CALL MPI_BCAST(p_end_systole_out_2,1,MPI_REAL8,rank_outlet_2,COMM_CART,merror)
!
!       CALL MPI_BCAST(pre_inlet,1,MPI_REAL8,rank_inlet,COMM_CART,merror)
!       CALL MPI_BCAST(pre_outlet,1,MPI_REAL8,rank_outlet,COMM_CART,merror)
!       CALL MPI_BCAST(pre_outlet_2,1,MPI_REAL8,rank_outlet_2,COMM_CART,merror)
!
!       if (rank == 0) print*, pre_inlet, pre_outlet, pre_outlet_2
!

     if (sphere_sim == 1) then
         do ii=S1p, N1p
            do jj=S2p, N2p
               do kk=S3p, N3p
                  if (sqrt((x1u(ii)-0.495)**2 + (x2p(jj)-0.5)**2 +(x3p(kk)- 1.)**2) .le. 0.0714)   vel(ii,jj,kk,1) = 0.0
                  if (sqrt((x1p(ii)-0.495)**2 + (x2v(jj)-0.5)**2 +(x3p(kk)- 1.)**2) .le. 0.0714)   vel(ii,jj,kk,2) = 0.0
                  if (sqrt((x1p(ii)-0.495)**2 + (x2p(jj)-0.5)**2 +(x3w(kk)- 1.)**2) .le. 0.0714)   vel(ii,jj,kk,3) = 0.0
               enddo
            enddo
         enddo
     endif


      

     !=======================================================================================================




     !--- send_signal.txt lesen ------------------------------------------------------------------------------
     CALL check_signal
     
     
     !--- Druck-Niveau festhalten ----------------------------------------------------------------------------
     CALL level_pressure
     !CALL level_pressure_aorta     
     
     !--- Ausschreiben ---------------------------------------------------------------------------------------
    
!     call MPI_Barrier(MPI_COMM_WORLD, merror)
!     t1= MPI_Wtime()

     IF (write_xdmf_yes .AND. write_out_vect) CALL write_xdmf_xml ! bbecsek
     IF (write_out_scal) CALL compute_stats
     IF (write_out_vect) CALL write_fields
 
!     call MPI_Barrier(MPI_COMM_WORLD, merror)
!     t2= MPI_Wtime()

 !    print*, 'printOutTime:', t2-t1 
     
     !--------------------------------------------------------------------------------------------------------
     IF (rank == 0 .AND. log_iteration_yes) CLOSE(10)
     
     IF (rank == 0 .AND. write_stout_yes .AND. timeint_mode == 0) WRITE(*,*)
     IF (rank == 0 .AND. write_stout_yes .AND. timeint_mode == 0) WRITE(*,*)
     
     CALL MPI_BCAST(finish_yes,1,MPI_LOGICAL,0,COMM_CART,merror) ! notwendig fuer "check_alarm"
     IF (finish_yes) EXIT
     
  END DO timeint
  !===========================================================================================================
  
  
  !--- Zeitmessung beenden -----------------------------------------------------------------------------------
  IF (rank == 0) THEN
     CALL DATE_AND_TIME(values=ctime)
     hour = ctime(5)
     minu = ctime(6)
     sec  = ctime(7)
     msec = ctime(8)
     
     IF (ctime(3) /= day) THEN
        ! Anmerkung: Gilt nur für Jobs <= 24h
        elatime = msec+1000*(sec+60*(minu+60*hour)) - elatime + 24*60*60*1000
     ELSE
        elatime = msec+1000*(sec+60*(minu+60*hour)) - elatime
     END IF
     
     WRITE(99,'(a,i2,a,i2,a,i4,a,i2,a,i2,a,i2,a,i3)') 'Finish time integration at ', ctime(3),'.',ctime(2),'.',ctime(1),    &
                                                           &       ', ',ctime(5),':',ctime(6),':',ctime(7),'.',ctime(8)
     WRITE(99,'(a,E13.5)') 'elapsed time [sec]', DBLE(elatime)/1000.
     CLOSE(99)
  END IF
  
  !--- Restart schreiben -------------------------------------------------------------------------------------
  restart = restart + 1
  CALL write_restart
!  IF (IB_on .AND. rank==0) CALL write_restart_ibm ! bbecsek
  CALL write_restart_stats
  
  
  !--- Iterationsstatistiken auswerten -----------------------------------------------------------------------
  CALL iteration_stats
  
  !--- diverse Files schliessen ------------------------------------------------------------------------------
  CALL close_stats
 
  !--- link all XDMF files together --------------------------------------------------------------------------
  IF (write_xdmf_yes) CALL write_xdmf_timecollection ! bbecsek

  END SUBROUTINE timeintegration
  
  
  
END MODULE mod_timeint
