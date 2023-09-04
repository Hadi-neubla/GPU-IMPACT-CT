!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!* modified by Barna Becsek, ARTORG CVE, University of Bern (barna.becsek@artorg.unibe.ch)                   *
!* October 2014                                                                                              *
!* GPU version by Hadi Zolfaghari , ARTORG CVE, then DAMTP, Cambridge University (hz382@cam.ac.uk)           *
!* Oct 2015 - Sep 2023                                                                                       *
!*************************************************************************************************************


!> module containing solver routines and preconditioners for solving the NSE. It uses modules mod_dim,
!! mod_vars, mod_exchange, mod_diff, mod_laplace, mod_helmholtz and mod_inout.
!! @note: what is the difference between BiCGstab and BiCGstab2? (one has an extra option for problem_type, 
!!        anything else?)
MODULE mod_solvers
  

!  use magma
  
  USE mod_dims
  USE mod_vars
  USE mod_vars_GPU
  USE mod_exchange
  USE mod_diff
  USE mod_laplace
  USE mod_inout
  !USE mpi

  PRIVATE
  
  PUBLIC outer_iteration, explicit, twostep
  PUBLIC force_massflow, apply_nullspace, get_nullspace, solve_nullspace
  PUBLIC solve_Helmholtz
  PUBLIC multigridV, multigridF, restrict, interpolate
  PUBLIC restrict_Helmholtz, interpolate_Helmholtz
  PUBLIC BiCGstab, Richardson
  PUBLIC get_norms, product_scalar, multadd1, multadd2
  PUBLIC status_iteration
  PUBLIC handle_corner_rhs
  
  PUBLIC BiCGstab2, get_norms2, product_scalar2, apply_nullspace2 ! TEST!!!
  PUBLIC relax_restrict, interpolate_relax, plain_restrict, interpolate_mg, relax_bottom
  
  
  INCLUDE 'mpif.h'
  
  CONTAINS
  
!pgi$g unroll = n:8
!!pgi$r unroll = n:8
!!pgi$l unroll = n:8
  
  
  
  
  

  
  
  
  
  
  SUBROUTINE explicit
  
  IMPLICIT NONE
  
  INTEGER                ::  i, j, k, m
  real(8), volatile         ::  t1, t2, t3  
  real(8)                   ::  val1, val2 
  integer                ::  n_match

  REAL(8)                   ::  gpre_test1(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL(8)                   ::  GPUres(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL(8)                   ::  gpre_test2(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3)
  REAL(8)                   ::  axpby_test(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))

  real(8)                   ::  norm_lap, norm_lap_global, norm_res_inf, norm_res_inf_global

!  CHARACTER(LEN=*)       :: filename


  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Beim zweiten Poisson-Problem im Vorkonditionierer hat es sich bislang bewährt, die Lösung !
  !                aus dem ersten Poisson-Problem als Startfeld zu benutzen.                                 !
  !              - Richardson-Iteration scheint langsamer als BiCGstab zu sein (~30%).                       !
  !----------------------------------------------------------------------------------------------------------!
  
  
  !===========================================================================================================
  !=== RHS fuer Druck-Gleichungssystem =======================================================================
  !===========================================================================================================
  DO m = 1, dimens
     
     ! TEST!!! weglassen und stattdessen rhs verwenden? Testroutinen lassen sich dann nicht mehr verwenden ...
     IF (m == 1) vel(S11B:N11B,S21B:N21B,S31B:N31B,1) = rhs(S11B:N11B,S21B:N21B,S31B:N31B,1)

     IF (m == 2)  vel(S12B:N12B,S22B:N22B,S32B:N32B,2) = rhs(S12B:N12B,S22B:N22B,S32B:N32B,2)

     IF (m == 3)  vel(S13B:N13B,S23B:N23B,S33B:N33B,3) = rhs(S13B:N13B,S23B:N23B,S33B:N33B,3)
     !--- Randbedingungen extrapolieren ----------------------------------------------------------------------
     
  END DO
  
  !--- Divergenz bilden --------------------------------------------------------------------------------------
  CALL divergence2(vel,res)
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== Poisson-Lösung ========================================================================================
  !===========================================================================================================
  ! Kanten/Ecken sind durch die Randbedingungen immer Divergenz-frei:
  
  IF (rank == 0 .AND. log_iteration_yes) WRITE(10,'(a)') 'Poisson problem:'
  
  number_poisson = 1

!  if ((substep ==1) .and. (rank==0)) OPEN(UNIT=139,FILE="IMPACTPressureSolverError.txt",FORM="FORMATTED",position="append",STATUS="OLD",ACTION="READWRITE")
  
  call MPI_Barrier(MPI_COMM_WORLD, merror)
  t1 =   MPI_Wtime()

 if ((timestep .le. 3) .or. (GPU_Poisson_Solver == 0)) then  
     CALL BiCGstab_orig  (epsU,n_it_Poisson,init_pre(substep),S1p,S2p,S3p,N1p,N2p,N3p,res,pre,2,.TRUE.,.TRUE.,precond_Poisson)
!   call device_jacobi_solver(1,0.8,S1p,S2p,S3p,N1p,N2p,N3p,pre,res)
 else

   call device_jacobi_solver(1,dble(0.8),S1p,S2p,S3p,N1p,N2p,N3p,pre,res)
  
 endif

  call MPI_Barrier(MPI_COMM_WORLD, merror)
  t2 =   MPI_Wtime()


  if (rank==0) print*, 'pressure solver time', t2-t1

!! uncomment for error output
      CALL exchange_relax(1,0,0,0,0,.FALSE.,pre)

      norm_lap = 0.
     ! norm_res_inf =0.
      DO k = S3p, N3p
         DO j = S2p, N2p
            DO i = S1p, N1p
               GPUres(i,j,k) =0.
               GPUres(i,j,k) = res(i,j,k)                                     &
                           &      - cdg1(-1,i,1)*pre(i-1,j,k) - cdg1(1,i,1)*pre(i+1,j,k)     &
                           &      - cdg2(-1,j,1)*pre(i,j-1,k) - cdg2(1,j,1)*pre(i,j+1,k)     &
                           &      - cdg3(-1,k,1)*pre(i,j,k-1) - cdg3(1,k,1)*pre(i,j,k+1)    &
                           &      - (cdg1(0,i,1) + cdg2(0,j,1) + cdg3(0,k,1))*pre(i,j,k)
               norm_lap = MAX(ABS(GPUres(i,j,k)),norm_lap)
               !norm_res_inf = MAX(ABS(res(i,j,k)),norm_res_inf)
           !     if (abs(GPUres(i,j,k)) .le. 1) print*, i,j,k, GPUres(i,j,k), substep
           END DO
         END DO
      END DO
     CALL MPI_ALLREDUCE(norm_lap,norm_lap_global,1,MPI_REAL8,MPI_MAX,COMM_CART,merror) ! hz
    ! CALL MPI_ALLREDUCE(norm_res_inf,norm_res_inf_global,1,MPI_REAL8,MPI_MAX,COMM_CART,merror) ! hz
     norm_lap = norm_lap_global
     !norm_res_inf = norm_res_inf_global
     !========================================================================================================
        if (rank ==0 ) print*, 'final residual', norm_lap



!  if ((substep ==1) .and. (rank==0)) WRITE(UNIT=139, FMT=*) (t2-t1)
!  if ((substep ==1) .and. (rank==0)) WRITE(UNIT=139, FMT=*) norm_lap

!  if ((substep ==1) .and. (rank==0)) CLOSE(UNIT=139)

!  if (rank == 0) print*, 'the BiCGstab converged in the time :', t2-t1

  !CALL Richardson(epsU,n_it_Poisson,init_pre(substep),S1p,S2p,S3p,N1p,N2p,N3p,res,pre,2,.TRUE.,.TRUE.,precond_Poisson)
  
  IF (rank == 0 .AND. log_iteration_yes) CALL flush(10)
  !===========================================================================================================
 
  
  !===========================================================================================================
  !=== Korrektur der Lösung ==================================================================================
  !===========================================================================================================
  DO m = 1, dimens
         
     !gpre_test1 = gpre
!#ifdef cuIMPACT
  if (GPU_accelerated == 1) then

 !    call cpu_time(t1)
!  call MPI_Barrier(MPI_COMM_WORLD, merror)

! t1 =   MPI_Wtime()

     call device_gradient(m,pre,gpre)

!  call MPI_Barrier(MPI_COMM_WORLD, merror)
!  t2 =   MPI_Wtime()

! print*, t2-t1, 'GPU time'
  !   call cpu_time(t2) 
!     print*, t2-t1, 'GPU time',m

!     print*, '---------- ---------- -----------on GPU gradient'


!     print*, 'on GPU'
!     print*, 'max value: ', maxval(gpre)
!#else
  else

 ! call cpu_time(t1)

!  call MPI_Barrier(MPI_COMM_WORLD, merror)
!  t1 =   MPI_Wtime()

     CALL gradient(m,pre,gpre)
    ! pre=1.0
 !   call device_gradient(m,pre,gpre_test1)
!print*,   'the pressure value from CPU', pre(:, :, :)

!  call MPI_Barrier(MPI_COMM_WORLD, merror)
!  t2 =   MPI_Wtime()


!  print*,  t2- t1, 'CPU time'
!     print*, 'on CPU'
!     print*, 'max value: ', maxval(gpre)
!  call cpu_time(t2)
!  print*, t2-t1, 'CPU_time'

  
  !endif

!#endif

!   n_match =0
!
!    do i=1, N11
!        do j=1, N21
!           do k=1, N31
!
!              if (abs(gpre_test1(i,j,k) - gpre(i,j,k)) .ge. 0.00000000000001 ) then
!                 n_match = n_match +1
!                 if (n_match .eq. 356) then 
!                    print*, gpre(i,j,k), 'CPU value'
!                    print*, gpre_test1(i,j,k), 'GPU value'
!                    print*, abs(gpre_test1(i,j,k) - gpre(i,j,k))
!                 endif
!               endif
!           enddo
!        enddo
!     enddo 
!
!     print*, m, n_match

 endif
!!! 
!!! 
!!!      print*, sum(gpre(:,:,:)), 'max of gpu result'
!!!      print*, sum(gpre_test1(:,:,:)), 'max of cpu result'
!!!
!!!
!!!  
!!!      print*, t3-t2, t2-t1, 'gpu time', 'cpu time'



     !--- Randbedingungen extrapolieren ----------------------------------------------------------------------
    !---HZ: Null for periodic boundary conditions 
    
    
     IF (m == 1) THEN
!*     gpre_test1 = gpre
!*     axpby_test = vel(:,:,:,1)


!  call MPI_Barrier(MPI_COMM_WORLD, merror)
!  t1 =   MPI_Wtime()
    if (GPU_accelerated == 1) then
        call device_axpby(dble(1.) ,dble(-1.) , vel(:,:,:,m), gpre, vel(:,:,:,m))
    else
!  call MPI_Barrier(MPI_COMM_WORLD, merror)
!  t2 =   MPI_Wtime()
!*      call device_axpby(1.,-1.,axpby_test, gpre_test1, axpby_test) 
         DO k = S31B, N31B
           DO j = S21B, N21B
!pgi$ unroll = n:8
              DO i = S11B, N11B
                 vel(i,j,k,m) = vel(i,j,k,m) - gpre(i,j,k)
!*                 if (abs(vel(i,j,k,m)-axpby_test(i,j,k)) .ge. 0.0000000000000001) print*, 'mismatch', vel(i,j,k,m), axpby_test(i,j,k) 
              END DO
           END DO
        END DO

    endif
 !    call MPI_Barrier(MPI_COMM_WORLD, merror)
 ! t3 =   MPI_Wtime()

 !   print*, 'GPU time == ', t3-t2, 'CPU time == ', t2-t1      


    END IF
     IF (m == 2) THEN
     if (GPU_accelerated == 1) then
     call device_axpby(dble(1.) ,dble(-1.) , vel(:,:,:,m), gpre, vel(:,:,:,m))
     else
        DO k = S32B, N32B
           DO j = S22B, N22B
!pgi$ unroll = n:8
              DO i = S12B, N12B
                 vel(i,j,k,m) = vel(i,j,k,m) - gpre(i,j,k)
              END DO
           END DO
        END DO
     endif
!        call device_axpby(1. ,-1. , vel(:,:,:,2), gpre(:,:,:), vel(:,:,:,2))

     END IF
     IF (m == 3) THEN
     if (GPU_accelerated == 1) then
     call device_axpby(dble(1.) ,dble(-1.) , vel(:,:,:,m), gpre, vel(:,:,:,m))
     else 
        DO k = S33B, N33B
           DO j = S23B, N23B
!pgi$ unroll = n:8
              DO i = S13B, N13B
                 vel(i,j,k,m) = vel(i,j,k,m) - gpre(i,j,k)
              END DO
           END DO
        END DO

     endif

!        call device_axpby(1. ,-1. , vel(:,:,:,3), gpre(:,:,:), vel(:,:,:,3))
     END IF

    
     
  END DO
  !===========================================================================================================
  
  
 
  
  END SUBROUTINE explicit
  
  
 
  
  
  !> performs a restriction step with a preliminary relaxation step
  SUBROUTINE relax_restrict(init_yes,nullspace_yes,g,psi,bb,phi,work,coarse,problem_type)
  
  IMPLICIT NONE
  
  LOGICAL, INTENT(in   ) ::  init_yes
  LOGICAL, INTENT(in   ) ::  nullspace_yes
  INTEGER, INTENT(in   ) ::  g
  
  REAL(8)   , INTENT(in   ) ::  psi   (b1L:(NN(1,g  )+b1U),b2L:(NN(2,g  )+b2U),b3L:(NN(3,g  )+b3U))
  REAL(8)   , INTENT(inout) ::  bb    (b1L:(NN(1,g  )+b1U),b2L:(NN(2,g  )+b2U),b3L:(NN(3,g  )+b3U))
  REAL(8)   , INTENT(inout) ::  phi   (b1L:(NN(1,g  )+b1U),b2L:(NN(2,g  )+b2U),b3L:(NN(3,g  )+b3U))
  REAL(8)   , INTENT(inout) ::  work  (b1L:(NN(1,g  )+b1U),b2L:(NN(2,g  )+b2U),b3L:(NN(3,g  )+b3U))
  REAL(8)   , INTENT(out  ) ::  coarse(b1L:(NN(1,g+1)+b1U),b2L:(NN(2,g+1)+b2U),b3L:(NN(3,g+1)+b3U))
  
  INTEGER, INTENT(in   ) ::  problem_type
  
  !-----------------------------------------------------------------------------------------------------------
     CALL relaxation_div_grad(init_yes,n_relax_down,g,bb,phi)
     CALL product_div_grad_relax(g,phi,work)
     CALL restrict(.TRUE.,g,coarse,bb,work)
  !-----------------------------------------------------------------------------------------------------------
  
  
  
  END SUBROUTINE relax_restrict
  
  
  
  
  
  
  
  
  
  
  !> interpolation with a relaxation step at each interpolation 
  SUBROUTINE interpolate_relax(g,bb,phi,work,coarse,problem_type)
  
  IMPLICIT NONE
  
  INTEGER, INTENT(in   ) ::  g
  
  REAL(8)   , INTENT(in   ) ::  bb    (b1L:(NN(1,g  )+b1U),b2L:(NN(2,g  )+b2U),b3L:(NN(3,g  )+b3U))
  REAL(8)   , INTENT(inout) ::  phi   (b1L:(NN(1,g  )+b1U),b2L:(NN(2,g  )+b2U),b3L:(NN(3,g  )+b3U))
  REAL(8)   , INTENT(out  ) ::  work  (b1L:(NN(1,g  )+b1U),b2L:(NN(2,g  )+b2U),b3L:(NN(3,g  )+b3U))
  REAL(8)   , INTENT(inout) ::  coarse(b1L:(NN(1,g+1)+b1U),b2L:(NN(2,g+1)+b2U),b3L:(NN(3,g+1)+b3U))
  
  INTEGER, INTENT(in   ) ::  problem_type
  
  
  !-----------------------------------------------------------------------------------------------------------
     CALL interpolate(.TRUE.,g,coarse,phi,work)
     CALL relaxation_div_grad_inv(.FALSE.,n_relax_up,g,bb,phi)
  !-----------------------------------------------------------------------------------------------------------
  
  END SUBROUTINE interpolate_relax
  
  
  
  
  
  
  
  
  
  
  !> performs a relaxation step at the coarsest level
  SUBROUTINE relax_bottom(init_yes,nullspace_yes,g,psi,bb,phi,problem_type)
  
  IMPLICIT NONE
  
  LOGICAL, INTENT(in   ) ::  init_yes
  LOGICAL, INTENT(in   ) ::  nullspace_yes
  INTEGER, INTENT(in   ) ::  g
  
  REAL(8)   , INTENT(in   ) ::  psi(b1L:(NN(1,g)+b1U),b2L:(NN(2,g)+b2U),b3L:(NN(3,g)+b3U))
  REAL(8)   , INTENT(inout) ::  bb (b1L:(NN(1,g)+b1U),b2L:(NN(2,g)+b2U),b3L:(NN(3,g)+b3U))
  REAL(8)   , INTENT(inout) ::  phi(b1L:(NN(1,g)+b1U),b2L:(NN(2,g)+b2U),b3L:(NN(3,g)+b3U))
  
  INTEGER, INTENT(in   ) ::  problem_type
  
  
  !-----------------------------------------------------------------------------------------------------------
     CALL relaxation_div_grad    (init_yes,n_relax_bottom,g,bb,phi)
     CALL relaxation_div_grad_inv(.FALSE. ,n_relax_bottom,g,bb,phi)
  !-----------------------------------------------------------------------------------------------------------
  
  
  !-----------------------------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------------------------
  
  
  !-----------------------------------------------------------------------------------------------------------
  ! bbecsek 150106: problem_type == 3 has been removed (associated with concentrations)
  !-----------------------------------------------------------------------------------------------------------
  
  
  END SUBROUTINE relax_bottom
  
  
  
  
  !> computes a V-cycle MG smoothing (preconditioning) 
  SUBROUTINE multigridV(init_yes,gstart,bb,phi,problem_type)
  
  IMPLICIT NONE
  
  LOGICAL, INTENT(in   ) ::  init_yes
  INTEGER, INTENT(in   ) ::  gstart
  INTEGER, INTENT(in   ) ::  problem_type
  
  REAL(8)   , INTENT(inout) ::  bb (b1L:(NN(1,gstart)+b1U),b2L:(NN(2,gstart)+b2U),b3L:(NN(3,gstart)+b3U))
  REAL(8)   , INTENT(inout) ::  phi(b1L:(NN(1,gstart)+b1U),b2L:(NN(2,gstart)+b2U),b3L:(NN(3,gstart)+b3U))
  
  INTEGER                ::  g
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Übergebene Felder sind generell zu gross (b1L:(N1+b1U) vs. 0:(N1+1)), andererseits stört  !
  !                der Umstand i.A. auch nicht!                                                              !
  !              - vec1C bzw. vecC können eingespart werden, wenn Produkt und Restriktion zusammengelegt     !
  !                werden. Darauf wird hier der Übersicht halber verzichtet, zumal vecC nicht auf dem        !
  !                feinsten Gitter existiert und somit der Speicher nicht überstrapaziert wird. Diese Felder !
  !                werden zudem an die Interpolationsroutinen als Arbeitsfelder übergeben, wären aber auch   !
  !                dort mit entsprechendem Programmieraufwand eliminierbar.                                  !
  !              - INTENT(inout) für bb ist wegen restrict_Helmholtz notwendig.                              !
  !              - Es wird vorausgesetzt, dass das feinstes Gitterlevel keine Nullraum-Korrektur benötigt    !
  !                (sollte normalerweise auch erfüllt sein).                                                 !
  !----------------------------------------------------------------------------------------------------------!
  
! print*, participate_yes 
  !===========================================================================================================
  DO g = gstart, n_grids-1
     
     IF (participate_yes(g)) THEN
        IF (g == 1 )    CALL relax_restrict(init_yes,.FALSE.             ,g,psi_rel1 ,bb    ,phi   ,vec1C ,vec2A ,problem_type)
        IF (g == gstart) THEN
           IF (g == 2 ) CALL relax_restrict(init_yes,.FALSE.             ,g,psi_rel2 ,bb    ,phi   ,vec2C ,vec3A ,problem_type)
           IF (g == 3 ) CALL relax_restrict(init_yes,.FALSE.             ,g,psi_rel3 ,bb    ,phi   ,vec3C ,vec4A ,problem_type)
           IF (g == 4 ) CALL relax_restrict(init_yes,.FALSE.             ,g,psi_rel4 ,bb    ,phi   ,vec4C ,vec5A ,problem_type)
           IF (g == 5 ) CALL relax_restrict(init_yes,.FALSE.             ,g,psi_rel5 ,bb    ,phi   ,vec5C ,vec6A ,problem_type)
           IF (g == 6 ) CALL relax_restrict(init_yes,.FALSE.             ,g,psi_rel6 ,bb    ,phi   ,vec6C ,vec7A ,problem_type)
           IF (g == 7 ) CALL relax_restrict(init_yes,.FALSE.             ,g,psi_rel7 ,bb    ,phi   ,vec7C ,vec8A ,problem_type)
           IF (g == 8 ) CALL relax_restrict(init_yes,.FALSE.             ,g,psi_rel8 ,bb    ,phi   ,vec8C ,vec9A ,problem_type)
           IF (g == 9 ) CALL relax_restrict(init_yes,.FALSE.             ,g,psi_rel9 ,bb    ,phi   ,vec9C ,vec10A,problem_type)
           IF (g == 10) CALL relax_restrict(init_yes,.FALSE.             ,g,psi_rel10,bb    ,phi   ,vec10C,vec11A,problem_type)
           IF (g == 11) CALL relax_restrict(init_yes,.FALSE.             ,g,psi_rel11,bb    ,phi   ,vec11C,vec12A,problem_type)
           IF (g == 12) CALL relax_restrict(init_yes,.FALSE.             ,g,psi_rel12,bb    ,phi   ,vec12C,vec13A,problem_type)
           IF (g == 13) CALL relax_restrict(init_yes,.FALSE.             ,g,psi_rel13,bb    ,phi   ,vec13C,vec14A,problem_type)
           IF (g == 14) CALL relax_restrict(init_yes,.FALSE.             ,g,psi_rel14,bb    ,phi   ,vec14C,vec15A,problem_type)
        ELSE
           IF (g == 2 ) CALL relax_restrict(.TRUE.  ,nullspace_coarse_yes,g,psi_rel2 ,vec2A ,vec2B ,vec2C ,vec3A ,problem_type)
           IF (g == 3 ) CALL relax_restrict(.TRUE.  ,nullspace_coarse_yes,g,psi_rel3 ,vec3A ,vec3B ,vec3C ,vec4A ,problem_type)
           IF (g == 4 ) CALL relax_restrict(.TRUE.  ,nullspace_coarse_yes,g,psi_rel4 ,vec4A ,vec4B ,vec4C ,vec5A ,problem_type)
           IF (g == 5 ) CALL relax_restrict(.TRUE.  ,nullspace_coarse_yes,g,psi_rel5 ,vec5A ,vec5B ,vec5C ,vec6A ,problem_type)
           IF (g == 6 ) CALL relax_restrict(.TRUE.  ,nullspace_coarse_yes,g,psi_rel6 ,vec6A ,vec6B ,vec6C ,vec7A ,problem_type)
           IF (g == 7 ) CALL relax_restrict(.TRUE.  ,nullspace_coarse_yes,g,psi_rel7 ,vec7A ,vec7B ,vec7C ,vec8A ,problem_type)
           IF (g == 8 ) CALL relax_restrict(.TRUE.  ,nullspace_coarse_yes,g,psi_rel8 ,vec8A ,vec8B ,vec8C ,vec9A ,problem_type)
           IF (g == 9 ) CALL relax_restrict(.TRUE.  ,nullspace_coarse_yes,g,psi_rel9 ,vec9A ,vec9B ,vec9C ,vec10A,problem_type)
           IF (g == 10) CALL relax_restrict(.TRUE.  ,nullspace_coarse_yes,g,psi_rel10,vec10A,vec10B,vec10C,vec11A,problem_type)
           IF (g == 11) CALL relax_restrict(.TRUE.  ,nullspace_coarse_yes,g,psi_rel11,vec11A,vec11B,vec11C,vec12A,problem_type)
           IF (g == 12) CALL relax_restrict(.TRUE.  ,nullspace_coarse_yes,g,psi_rel12,vec12A,vec12B,vec12C,vec13A,problem_type)
           IF (g == 13) CALL relax_restrict(.TRUE.  ,nullspace_coarse_yes,g,psi_rel13,vec13A,vec13B,vec13C,vec14A,problem_type)
           IF (g == 14) CALL relax_restrict(.TRUE.  ,nullspace_coarse_yes,g,psi_rel14,vec14A,vec14B,vec14C,vec15A,problem_type)
        END IF
     END IF
     
  END DO
  !===========================================================================================================
  
  !--- Grob-Gitter Lösung -------------------------------------
  IF (participate_yes(n_grids)) THEN
     IF (gstart == n_grids) THEN
                           CALL relax_bottom(init_yes,.FALSE.             ,n_grids,psi_rel1 ,bb    ,phi   ,problem_type) ! Achtung: psi_rel1 ist i.A. zu gross! (nullspace_coarse == .FALSE. <-- unschön!)
     ELSE
        IF (n_grids == 2 ) CALL relax_bottom(.TRUE.  ,nullspace_coarse_yes,n_grids,psi_rel2 ,vec2A ,vec2B ,problem_type)
        IF (n_grids == 3 ) CALL relax_bottom(.TRUE.  ,nullspace_coarse_yes,n_grids,psi_rel3 ,vec3A ,vec3B ,problem_type)
        IF (n_grids == 4 ) CALL relax_bottom(.TRUE.  ,nullspace_coarse_yes,n_grids,psi_rel4 ,vec4A ,vec4B ,problem_type)
        IF (n_grids == 5 ) CALL relax_bottom(.TRUE.  ,nullspace_coarse_yes,n_grids,psi_rel5 ,vec5A ,vec5B ,problem_type)
        IF (n_grids == 6 ) CALL relax_bottom(.TRUE.  ,nullspace_coarse_yes,n_grids,psi_rel6 ,vec6A ,vec6B ,problem_type)
        IF (n_grids == 7 ) CALL relax_bottom(.TRUE.  ,nullspace_coarse_yes,n_grids,psi_rel7 ,vec7A ,vec7B ,problem_type)
        IF (n_grids == 8 ) CALL relax_bottom(.TRUE.  ,nullspace_coarse_yes,n_grids,psi_rel8 ,vec8A ,vec8B ,problem_type)
        IF (n_grids == 9 ) CALL relax_bottom(.TRUE.  ,nullspace_coarse_yes,n_grids,psi_rel9 ,vec9A ,vec9B ,problem_type)
        IF (n_grids == 10) CALL relax_bottom(.TRUE.  ,nullspace_coarse_yes,n_grids,psi_rel10,vec10A,vec10B,problem_type)
        IF (n_grids == 11) CALL relax_bottom(.TRUE.  ,nullspace_coarse_yes,n_grids,psi_rel11,vec11A,vec11B,problem_type)
        IF (n_grids == 12) CALL relax_bottom(.TRUE.  ,nullspace_coarse_yes,n_grids,psi_rel12,vec12A,vec12B,problem_type)
        IF (n_grids == 13) CALL relax_bottom(.TRUE.  ,nullspace_coarse_yes,n_grids,psi_rel13,vec13A,vec13B,problem_type)
        IF (n_grids == 14) CALL relax_bottom(.TRUE.  ,nullspace_coarse_yes,n_grids,psi_rel14,vec14A,vec14B,problem_type)
        IF (n_grids == 15) CALL relax_bottom(.TRUE.  ,nullspace_coarse_yes,n_grids,psi_rel15,vec15A,vec15B,problem_type)
     END IF
  END IF
  
  !===========================================================================================================
  DO g = n_grids-1, gstart, -1
     
     IF (participate_yes(g)) THEN
        IF (g == gstart) THEN
           IF (g == 14) CALL interpolate_relax(g,bb    ,phi   ,vec14C,vec15B,problem_type)
           IF (g == 13) CALL interpolate_relax(g,bb    ,phi   ,vec13C,vec14B,problem_type)
           IF (g == 12) CALL interpolate_relax(g,bb    ,phi   ,vec12C,vec13B,problem_type)
           IF (g == 11) CALL interpolate_relax(g,bb    ,phi   ,vec11C,vec12B,problem_type)
           IF (g == 10) CALL interpolate_relax(g,bb    ,phi   ,vec10C,vec11B,problem_type)
           IF (g == 9 ) CALL interpolate_relax(g,bb    ,phi   ,vec9C ,vec10B,problem_type)
           IF (g == 8 ) CALL interpolate_relax(g,bb    ,phi   ,vec8C ,vec9B ,problem_type)
           IF (g == 7 ) CALL interpolate_relax(g,bb    ,phi   ,vec7C ,vec8B ,problem_type)
           IF (g == 6 ) CALL interpolate_relax(g,bb    ,phi   ,vec6C ,vec7B ,problem_type)
           IF (g == 5 ) CALL interpolate_relax(g,bb    ,phi   ,vec5C ,vec6B ,problem_type)
           IF (g == 4 ) CALL interpolate_relax(g,bb    ,phi   ,vec4C ,vec5B ,problem_type)
           IF (g == 3 ) CALL interpolate_relax(g,bb    ,phi   ,vec3C ,vec4B ,problem_type)
           IF (g == 2 ) CALL interpolate_relax(g,bb    ,phi   ,vec2C ,vec3B ,problem_type)
        ELSE
           IF (g == 14) CALL interpolate_relax(g,vec14A,vec14B,vec14C,vec15B,problem_type)
           IF (g == 13) CALL interpolate_relax(g,vec13A,vec13B,vec13C,vec14B,problem_type)
           IF (g == 12) CALL interpolate_relax(g,vec12A,vec12B,vec12C,vec13B,problem_type)
           IF (g == 11) CALL interpolate_relax(g,vec11A,vec11B,vec11C,vec12B,problem_type)
           IF (g == 10) CALL interpolate_relax(g,vec10A,vec10B,vec10C,vec11B,problem_type)
           IF (g == 9 ) CALL interpolate_relax(g,vec9A ,vec9B ,vec9C ,vec10B,problem_type)
           IF (g == 8 ) CALL interpolate_relax(g,vec8A ,vec8B ,vec8C ,vec9B ,problem_type)
           IF (g == 7 ) CALL interpolate_relax(g,vec7A ,vec7B ,vec7C ,vec8B ,problem_type)
           IF (g == 6 ) CALL interpolate_relax(g,vec6A ,vec6B ,vec6C ,vec7B ,problem_type)
           IF (g == 5 ) CALL interpolate_relax(g,vec5A ,vec5B ,vec5C ,vec6B ,problem_type)
           IF (g == 4 ) CALL interpolate_relax(g,vec4A ,vec4B ,vec4C ,vec5B ,problem_type)
           IF (g == 3 ) CALL interpolate_relax(g,vec3A ,vec3B ,vec3C ,vec4B ,problem_type)
           IF (g == 2 ) CALL interpolate_relax(g,vec2A ,vec2B ,vec2C ,vec3B ,problem_type)
        END IF
        IF (g == 1 )    CALL interpolate_relax(g,bb    ,phi   ,vec1C ,vec2B ,problem_type)
     END IF
     
  END DO
  !===========================================================================================================
  
  
  END SUBROUTINE multigridV
  
  
    SUBROUTINE restrict(add_yes,g,coarse,fine1,fine2) ! TEST!!! aufräumen ...
  
  IMPLICIT NONE
  
  LOGICAL, INTENT(in   ) ::  add_yes
  INTEGER, INTENT(in   ) ::  g
  
  REAL(8)   , INTENT(out  ) ::  coarse(b1L:(NN(1,g+1)+b1U),b2L:(NN(2,g+1)+b2U),b3L:(NN(3,g+1)+b3U))
  REAL(8)   , INTENT(inout) ::  fine1 (b1L:(NN(1,g  )+b1U),b2L:(NN(2,g  )+b2U),b3L:(NN(3,g  )+b3U))
  REAL(8)   , INTENT(inout) ::  fine2 (b1L:(NN(1,g  )+b1U),b2L:(NN(2,g  )+b2U),b3L:(NN(3,g  )+b3U))
  
  INTEGER                ::  i, ii, di, imax, iimax
  INTEGER                ::  j, jj, dj, jmax, jjmax
  INTEGER                ::  k, kk, dk, kmax, kkmax
  
  INTEGER                ::  sizsg(1:3), offsg(1:3), dispg    
  REAL(8)   , ALLOCATABLE   ::  sendbuf(:,:,:)
  REAL(8)                   ::  recvbuf(1:NN(1,g+1)*NN(2,g+1)*NN(3,g+1))
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - für allgemeine di, dj, dk geeignet!                                                       !
  !              - überlappende Schicht in Blöcken wird (der Einfachheit halber) ebenfalls ausgetauscht, ist !
  !                aber im Prinzip redundant (genauer: coarse(S1R:N1R,S2R:N2R,S3R:N3R) = ...).               !
  !              - Motivation für diese kurze Routine ist die Möglichkeit, auch Varianten wie Full-Weighting !
  !                etc. ggf. einzubauen, ansonsten könnte sie auch eingespaart werden.                       !
  !              - Die Block-überlappenden Stirnflächen werden ebenfalls mitverarbeitet, aber eigentlich     !
  !                nicht gebraucht (erleichtert die Programmierung), so dass eine Initialisierung notwendig  !
  !                ist. Dies wiederum bedingt die INTENT(inout)-Deklaration.                                 !
  !----------------------------------------------------------------------------------------------------------!
  
  
  imax    = NN(1,g)
  jmax    = NN(2,g)
  kmax    = NN(3,g)
  
  iimax   = (NN(1,g+1)-1)/n_gather(1,g+1)+1
  jjmax   = (NN(2,g+1)-1)/n_gather(2,g+1)+1
  kkmax   = (NN(3,g+1)-1)/n_gather(3,g+1)+1
  
  di      = (NN(1,g)-1)/(iimax-1)
  dj      = (NN(2,g)-1)/(jjmax-1)
  dk      = (NN(3,g)-1)/(kkmax-1)
  
  
  IF (BC(1,1,g) .GT. 0 .AND. BC(1,2,g) .GT. 0) fine1(  1      ,  1      ,1:NN(3,g)) = (fine1(2        ,1      ,1:NN(3,g)) + fine1(1      ,2        ,1:NN(3,g)))/2.
  IF (BC(1,1,g) .GT. 0 .AND. BC(2,2,g) .GT. 0) fine1(  1      ,  NN(2,g),1:NN(3,g)) = (fine1(2        ,NN(2,g),1:NN(3,g)) + fine1(1      ,NN(2,g)-1,1:NN(3,g)))/2.
  IF (BC(2,1,g) .GT. 0 .AND. BC(1,2,g) .GT. 0) fine1(  NN(1,g),  1      ,1:NN(3,g)) = (fine1(NN(1,g)-1,1      ,1:NN(3,g)) + fine1(NN(1,g),2        ,1:NN(3,g)))/2.
  IF (BC(2,1,g) .GT. 0 .AND. BC(2,2,g) .GT. 0) fine1(  NN(1,g),  NN(2,g),1:NN(3,g)) = (fine1(NN(1,g)-1,NN(2,g),1:NN(3,g)) + fine1(NN(1,g),NN(2,g)-1,1:NN(3,g)))/2.
  
  IF (BC(1,1,g) .GT. 0 .AND. BC(1,3,g) .GT. 0) fine1(  1      ,1:NN(2,g),  1      ) = (fine1(2        ,1:NN(2,g),1      ) + fine1(1      ,1:NN(2,g),2        ))/2.
  IF (BC(1,1,g) .GT. 0 .AND. BC(2,3,g) .GT. 0) fine1(  1      ,1:NN(2,g),  NN(3,g)) = (fine1(2        ,1:NN(2,g),NN(3,g)) + fine1(1      ,1:NN(2,g),NN(3,g)-1))/2.
  IF (BC(2,1,g) .GT. 0 .AND. BC(1,3,g) .GT. 0) fine1(  NN(1,g),1:NN(2,g),  1      ) = (fine1(NN(1,g)-1,1:NN(2,g),1      ) + fine1(NN(1,g),1:NN(2,g),2        ))/2.
  IF (BC(2,1,g) .GT. 0 .AND. BC(2,3,g) .GT. 0) fine1(  NN(1,g),1:NN(2,g),  NN(3,g)) = (fine1(NN(1,g)-1,1:NN(2,g),NN(3,g)) + fine1(NN(1,g),1:NN(2,g),NN(3,g)-1))/2.

  IF (BC(1,2,g) .GT. 0 .AND. BC(1,3,g) .GT. 0) fine1(1:NN(1,g),  1      ,  1      ) = (fine1(1:NN(1,g),2        ,1      ) + fine1(1:NN(1,g),1      ,2        ))/2.
  IF (BC(1,2,g) .GT. 0 .AND. BC(2,3,g) .GT. 0) fine1(1:NN(1,g),  1      ,  NN(3,g)) = (fine1(1:NN(1,g),2        ,NN(3,g)) + fine1(1:NN(1,g),1      ,NN(3,g)-1))/2.
  IF (BC(2,2,g) .GT. 0 .AND. BC(1,3,g) .GT. 0) fine1(1:NN(1,g),  NN(2,g),  1      ) = (fine1(1:NN(1,g),NN(2,g)-1,1      ) + fine1(1:NN(1,g),NN(2,g),2        ))/2.
  IF (BC(2,2,g) .GT. 0 .AND. BC(2,3,g) .GT. 0) fine1(1:NN(1,g),  NN(2,g),  NN(3,g)) = (fine1(1:NN(1,g),NN(2,g)-1,NN(3,g)) + fine1(1:NN(1,g),NN(2,g),NN(3,g)-1))/2.
  
  
  IF (BC(1,1,g) .GT. 0 .AND. BC(1,2,g) .GT. 0) fine2(  1      ,  1      ,1:NN(3,g)) = (fine2(2        ,1      ,1:NN(3,g)) + fine2(1      ,2        ,1:NN(3,g)))/2.
  IF (BC(1,1,g) .GT. 0 .AND. BC(2,2,g) .GT. 0) fine2(  1      ,  NN(2,g),1:NN(3,g)) = (fine2(2        ,NN(2,g),1:NN(3,g)) + fine2(1      ,NN(2,g)-1,1:NN(3,g)))/2.
  IF (BC(2,1,g) .GT. 0 .AND. BC(1,2,g) .GT. 0) fine2(  NN(1,g),  1      ,1:NN(3,g)) = (fine2(NN(1,g)-1,1      ,1:NN(3,g)) + fine2(NN(1,g),2        ,1:NN(3,g)))/2.
  IF (BC(2,1,g) .GT. 0 .AND. BC(2,2,g) .GT. 0) fine2(  NN(1,g),  NN(2,g),1:NN(3,g)) = (fine2(NN(1,g)-1,NN(2,g),1:NN(3,g)) + fine2(NN(1,g),NN(2,g)-1,1:NN(3,g)))/2.
  
  IF (BC(1,1,g) .GT. 0 .AND. BC(1,3,g) .GT. 0) fine2(  1      ,1:NN(2,g),  1      ) = (fine2(2        ,1:NN(2,g),1      ) + fine2(1      ,1:NN(2,g),2        ))/2.
  IF (BC(1,1,g) .GT. 0 .AND. BC(2,3,g) .GT. 0) fine2(  1      ,1:NN(2,g),  NN(3,g)) = (fine2(2        ,1:NN(2,g),NN(3,g)) + fine2(1      ,1:NN(2,g),NN(3,g)-1))/2.
  IF (BC(2,1,g) .GT. 0 .AND. BC(1,3,g) .GT. 0) fine2(  NN(1,g),1:NN(2,g),  1      ) = (fine2(NN(1,g)-1,1:NN(2,g),1      ) + fine2(NN(1,g),1:NN(2,g),2        ))/2.
  IF (BC(2,1,g) .GT. 0 .AND. BC(2,3,g) .GT. 0) fine2(  NN(1,g),1:NN(2,g),  NN(3,g)) = (fine2(NN(1,g)-1,1:NN(2,g),NN(3,g)) + fine2(NN(1,g),1:NN(2,g),NN(3,g)-1))/2.

  IF (BC(1,2,g) .GT. 0 .AND. BC(1,3,g) .GT. 0) fine2(1:NN(1,g),  1      ,  1      ) = (fine2(1:NN(1,g),2        ,1      ) + fine2(1:NN(1,g),1      ,2        ))/2.
  IF (BC(1,2,g) .GT. 0 .AND. BC(2,3,g) .GT. 0) fine2(1:NN(1,g),  1      ,  NN(3,g)) = (fine2(1:NN(1,g),2        ,NN(3,g)) + fine2(1:NN(1,g),1      ,NN(3,g)-1))/2.
  IF (BC(2,2,g) .GT. 0 .AND. BC(1,3,g) .GT. 0) fine2(1:NN(1,g),  NN(2,g),  1      ) = (fine2(1:NN(1,g),NN(2,g)-1,1      ) + fine2(1:NN(1,g),NN(2,g),2        ))/2.
  IF (BC(2,2,g) .GT. 0 .AND. BC(2,3,g) .GT. 0) fine2(1:NN(1,g),  NN(2,g),  NN(3,g)) = (fine2(1:NN(1,g),NN(2,g)-1,NN(3,g)) + fine2(1:NN(1,g),NN(2,g),NN(3,g)-1))/2.
  
  
  IF (ls1 ==  0 .AND. (BC(1,1,g) == 0 .OR. BC(1,1,g) == -1)) iimax = iimax-1
  IF (ls1 == -1 .AND. (BC(2,1,g) == 0 .OR. BC(2,1,g) == -1)) iimax = iimax-1
  
  IF (ls2 ==  0 .AND. (BC(1,2,g) == 0 .OR. BC(1,2,g) == -1)) jjmax = jjmax-1
  IF (ls2 == -1 .AND. (BC(2,2,g) == 0 .OR. BC(2,2,g) == -1)) jjmax = jjmax-1
  
  IF (ls3 ==  0 .AND. (BC(1,3,g) == 0 .OR. BC(1,3,g) == -1)) kkmax = kkmax-1
  IF (ls3 == -1 .AND. (BC(2,3,g) == 0 .OR. BC(2,3,g) == -1)) kkmax = kkmax-1
  
  
  IF (1 == 2) THEN ! TEST!!!
     IF (add_yes) THEN
!pgi$ unroll = n:8
        coarse(1:iimax,1:jjmax,1:kkmax) = fine1(1:imax:di,1:jmax:dj,1:kmax:dk) - fine2(1:imax:di,1:jmax:dj,1:kmax:dk)
     ELSE
        coarse(1:iimax,1:jjmax,1:kkmax) = fine1(1:imax:di,1:jmax:dj,1:kmax:dk)
     END IF
     
  ELSE
     
     IF (add_yes) THEN ! TEST!!! etwas seriöser einbauen ...
     
     CALL exchange_relax(g,0,0,0,0,.TRUE.,fine1)
     CALL exchange_relax(g,0,0,0,0,.TRUE.,fine2)
     
     IF (dimens == 3) THEN
        DO kk = 1, kkmax
           k = dk*(kk-1)+1
           DO jj = 1, jjmax
              j = dj*(jj-1)+1
              DO ii = 1, iimax
                 i = di*(ii-1)+1
                 coarse(ii,jj,kk) = ((cR1(0,ii,g+1)+cR2(0,jj,g+1)+cR3(0,kk,g+1))*(fine1(i,j,k)-fine2(i,j,k)) +   &
                         &              cR1(-1,ii,g+1)*(fine1(i-1,j,k)-fine2(i-1,j,k)) +  &
                         &              cR1( 1,ii,g+1)*(fine1(i+1,j,k)-fine2(i+1,j,k)) +  &
                         &              cR2(-1,jj,g+1)*(fine1(i,j-1,k)-fine2(i,j-1,k)) +  &
                         &              cR2( 1,jj,g+1)*(fine1(i,j+1,k)-fine2(i,j+1,k)) +  &
                         &              cR3(-1,kk,g+1)*(fine1(i,j,k-1)-fine2(i,j,k-1)) +  &
                         &              cR3( 1,kk,g+1)*(fine1(i,j,k+1)-fine2(i,j,k+1))) / 3.
              END DO
           END DO
        END DO
     ELSE
        k  = 1
        kk = 1
        DO jj = 1, jjmax
           j = dj*(jj-1)+1
           DO ii = 1, iimax
              i = di*(ii-1)+1
              coarse(ii,jj,kk) = ((cR1(0,ii,g+1)+cR2(0,jj,g+1))*(fine1(i,j,k)-fine2(i,j,k)) +   &
                      &              cR1(-1,ii,g+1)*(fine1(i-1,j,k)-fine2(i-1,j,k)) +  &
                      &              cR1( 1,ii,g+1)*(fine1(i+1,j,k)-fine2(i+1,j,k)) +  &
                      &              cR2(-1,jj,g+1)*(fine1(i,j-1,k)-fine2(i,j-1,k)) +  &
                      &              cR2( 1,jj,g+1)*(fine1(i,j+1,k)-fine2(i,j+1,k))) / 2.
           END DO
        END DO
     END IF
     
     ELSE
     
     CALL exchange_relax(g,0,0,0,0,.TRUE.,fine1)
     
     IF (dimens == 3) THEN
        DO kk = 1, kkmax
           k = dk*(kk-1)+1
           DO jj = 1, jjmax
              j = dj*(jj-1)+1
              DO ii = 1, iimax
                 i = di*(ii-1)+1
                 coarse(ii,jj,kk) = ((cR1(0,ii,g+1)+cR2(0,jj,g+1)+cR3(0,kk,g+1))*fine1(i,j,k) +   &
                         &              cR1(-1,ii,g+1)*fine1(i-1,j,k) +  &
                         &              cR1( 1,ii,g+1)*fine1(i+1,j,k) +  &
                         &              cR2(-1,jj,g+1)*fine1(i,j-1,k) +  &
                         &              cR2( 1,jj,g+1)*fine1(i,j+1,k) +  &
                         &              cR3(-1,kk,g+1)*fine1(i,j,k-1) +  &
                         &              cR3( 1,kk,g+1)*fine1(i,j,k+1)) / 3.
              END DO
           END DO
        END DO
     ELSE
        k  = 1
        kk = 1
        DO jj = 1, jjmax
           j = dj*(jj-1)+1
           DO ii = 1, iimax
              i = di*(ii-1)+1
              coarse(ii,jj,kk) = ((cR1(0,ii,g+1)+cR2(0,jj,g+1))*fine1(i,j,k) +   &
                      &              cR1(-1,ii,g+1)*fine1(i-1,j,k) +  &
                      &              cR1( 1,ii,g+1)*fine1(i+1,j,k) +  &
                      &              cR2(-1,jj,g+1)*fine1(i,j-1,k) +  &
                      &              cR2( 1,jj,g+1)*fine1(i,j+1,k)) / 2.
           END DO
        END DO
     END IF
     
     END IF
     
  END IF
  
  
  
  IF (n_gather(1,g+1)*n_gather(2,g+1)*n_gather(3,g+1) .GT. 1) THEN
     
     ALLOCATE(sendbuf(1:iimax,1:jjmax,1:kkmax)) ! Anmerkung: Besser nicht fest allocieren um Speicherplatz zu sparen, ODER gleich "coarse" verwenden!
     
     DO kk = 1, kkmax
        DO jj = 1, jjmax
           DO ii = 1, iimax
              sendbuf(ii,jj,kk) = coarse(ii,jj,kk)
           END DO
        END DO
     END DO
     
     CALL MPI_GATHERv(sendbuf,iimax*jjmax*kkmax,MPI_REAL8,recvbuf,recvR(1,g+1),dispR(1,g+1),MPI_REAL8,rankc2(g+1),comm2(g+1),merror)
     
     DEALLOCATE(sendbuf)
     
     
     IF (participate_yes(g+1)) THEN
        DO k = 1, n_gather(3,g+1)
           DO j = 1, n_gather(2,g+1)
              DO i = 1, n_gather(1,g+1)
                 
                 sizsg(1:3) = sizsR(1:3,i+(j-1)*n_gather(1,g+1)+(k-1)*n_gather(1,g+1)*n_gather(2,g+1),g+1)
                 offsg(1:3) = offsR(1:3,i+(j-1)*n_gather(1,g+1)+(k-1)*n_gather(1,g+1)*n_gather(2,g+1),g+1)
                 dispg      = dispR(    i+(j-1)*n_gather(1,g+1)+(k-1)*n_gather(1,g+1)*n_gather(2,g+1),g+1)
                 
                 DO kk = 1, sizsg(3)
                    DO jj = 1, sizsg(2)
                       DO ii = 1, sizsg(1)
                          coarse(ii+offsg(1),jj+offsg(2),kk+offsg(3)) = recvbuf(dispg+ii+(jj-1)*sizsg(1)+(kk-1)*sizsg(1)*sizsg(2))
                       END DO
                    END DO
                 END DO
                 
              END DO
           END DO
        END DO
     END IF
     
  END IF
  
  
  END SUBROUTINE restrict
 
  
  
  
  
  
  
  
  
  ! TEST!!! Hier lassen sich evtl. Operationen einsparen! (fine1-fine2) nur einmal berechnen, ist aber fraglich, ob das auch schneller ist ...
  SUBROUTINE restrict_shreded(add_yes,g,coarse,fine1,fine2) ! TEST!!! aufräumen ...
  
  IMPLICIT NONE
  
  LOGICAL, INTENT(in   ) ::  add_yes
  INTEGER, INTENT(in   ) ::  g
  
  REAL(8)   , INTENT(out  ) ::  coarse(b1L:(NN(1,g+1)+b1U),b2L:(NN(2,g+1)+b2U),b3L:(NN(3,g+1)+b3U))
  REAL(8)   , INTENT(inout) ::  fine1 (b1L:(NN(1,g  )+b1U),b2L:(NN(2,g  )+b2U),b3L:(NN(3,g  )+b3U))
  REAL(8)   , INTENT(inout) ::  fine2 (b1L:(NN(1,g  )+b1U),b2L:(NN(2,g  )+b2U),b3L:(NN(3,g  )+b3U))
  
  INTEGER                ::  i, ii, di, imax, iimax
  INTEGER                ::  j, jj, dj, jmax, jjmax
  INTEGER                ::  k, kk, dk, kmax, kkmax
  
  INTEGER                ::  sizsg(1:3), offsg(1:3), dispg    
  REAL(8)   , ALLOCATABLE   ::  sendbuf(:,:,:)
  REAL(8)                   ::  recvbuf(1:NN(1,g+1)*NN(2,g+1)*NN(3,g+1))
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - für allgemeine di, dj, dk geeignet!                                                       !
  !              - überlappende Schicht in Blöcken wird (der Einfachheit halber) ebenfalls ausgetauscht, ist !
  !                aber im Prinzip redundant (genauer: coarse(S1R:N1R,S2R:N2R,S3R:N3R) = ...).               !
  !              - Motivation für diese kurze Routine ist die Möglichkeit, auch Varianten wie Full-Weighting !
  !                etc. ggf. einzubauen, ansonsten könnte sie auch eingespaart werden.                       !
  !              - Die Block-überlappenden Stirnflächen werden ebenfalls mitverarbeitet, aber eigentlich     !
  !                nicht gebraucht (erleichtert die Programmierung), so dass eine Initialisierung notwendig  !
  !                ist. Dies wiederum bedingt die INTENT(inout)-Deklaration.                                 !
  !----------------------------------------------------------------------------------------------------------!
  
  
  imax    = NN(1,g)
  jmax    = NN(2,g)
  kmax    = NN(3,g)
  
  iimax   = (NN(1,g+1)-1)/n_gather(1,g+1)+1
  jjmax   = (NN(2,g+1)-1)/n_gather(2,g+1)+1
  kkmax   = (NN(3,g+1)-1)/n_gather(3,g+1)+1
  
  di      = (NN(1,g)-1)/(iimax-1)
  dj      = (NN(2,g)-1)/(jjmax-1)
  dk      = (NN(3,g)-1)/(kkmax-1)
  
  
 
  
  IF (ls1 ==  0 .AND. (BC(1,1,g) == 0 .OR. BC(1,1,g) == -1)) iimax = iimax-1
  IF (ls1 == -1 .AND. (BC(2,1,g) == 0 .OR. BC(2,1,g) == -1)) iimax = iimax-1
  
  IF (ls2 ==  0 .AND. (BC(1,2,g) == 0 .OR. BC(1,2,g) == -1)) jjmax = jjmax-1
  IF (ls2 == -1 .AND. (BC(2,2,g) == 0 .OR. BC(2,2,g) == -1)) jjmax = jjmax-1
  
  IF (ls3 ==  0 .AND. (BC(1,3,g) == 0 .OR. BC(1,3,g) == -1)) kkmax = kkmax-1
  IF (ls3 == -1 .AND. (BC(2,3,g) == 0 .OR. BC(2,3,g) == -1)) kkmax = kkmax-1
  
  
    
  IF (add_yes) THEN ! TEST!!! etwas seriöser einbauen ...
     
     CALL exchange_relax(g,0,0,0,0,.TRUE.,fine1)
     CALL exchange_relax(g,0,0,0,0,.TRUE.,fine2)
     
     IF (dimens == 3) THEN
        DO kk = 1, kkmax
           k = dk*(kk-1)+1
           DO jj = 1, jjmax
              j = dj*(jj-1)+1
              DO ii = 1, iimax
                 i = di*(ii-1)+1
                 coarse(ii,jj,kk) = ((cR1(0,ii,g+1)+cR2(0,jj,g+1)+cR3(0,kk,g+1))*(fine1(i,j,k)-fine2(i,j,k)) +   &
                         &              cR1(-1,ii,g+1)*(fine1(i-1,j,k)-fine2(i-1,j,k)) +  &
                         &              cR1( 1,ii,g+1)*(fine1(i+1,j,k)-fine2(i+1,j,k)) +  &
                         &              cR2(-1,jj,g+1)*(fine1(i,j-1,k)-fine2(i,j-1,k)) +  &
                         &              cR2( 1,jj,g+1)*(fine1(i,j+1,k)-fine2(i,j+1,k)) +  &
                         &              cR3(-1,kk,g+1)*(fine1(i,j,k-1)-fine2(i,j,k-1)) +  &
                         &              cR3( 1,kk,g+1)*(fine1(i,j,k+1)-fine2(i,j,k+1))) / 3.
              END DO
           END DO
        END DO
      END IF
     
    
  END IF
  
  
  IF (n_gather(1,g+1)*n_gather(2,g+1)*n_gather(3,g+1) .GT. 1) THEN
     ALLOCATE(sendbuf(1:iimax,1:jjmax,1:kkmax)) ! Anmerkung: Besser nicht fest allocieren um Speicherplatz zu sparen, ODER gleich "coarse" verwenden!
     
     DO kk = 1, kkmax
        DO jj = 1, jjmax
           DO ii = 1, iimax
              sendbuf(ii,jj,kk) = coarse(ii,jj,kk)
           END DO
        END DO
     END DO
     
     CALL MPI_GATHERv(sendbuf,iimax*jjmax*kkmax,MPI_REAL8,recvbuf,recvR(1,g+1),dispR(1,g+1),MPI_REAL8,rankc2(g+1),comm2(g+1),merror)
     
     DEALLOCATE(sendbuf)
     
     
     IF (participate_yes(g+1)) THEN
        DO k = 1, n_gather(3,g+1)
           DO j = 1, n_gather(2,g+1)
              DO i = 1, n_gather(1,g+1)
                 
                 sizsg(1:3) = sizsR(1:3,i+(j-1)*n_gather(1,g+1)+(k-1)*n_gather(1,g+1)*n_gather(2,g+1),g+1)
                 offsg(1:3) = offsR(1:3,i+(j-1)*n_gather(1,g+1)+(k-1)*n_gather(1,g+1)*n_gather(2,g+1),g+1)
                 dispg      = dispR(    i+(j-1)*n_gather(1,g+1)+(k-1)*n_gather(1,g+1)*n_gather(2,g+1),g+1)
                 
                 DO kk = 1, sizsg(3)
                    DO jj = 1, sizsg(2)
                       DO ii = 1, sizsg(1)
                          coarse(ii+offsg(1),jj+offsg(2),kk+offsg(3)) = recvbuf(dispg+ii+(jj-1)*sizsg(1)+(kk-1)*sizsg(1)*sizsg(2))
                       END DO
                    END DO
                 END DO
                 
              END DO
           END DO
        END DO
     END IF
     
  END IF
  
  
  END SUBROUTINE restrict_shreded
  
  
  
  
  
  
  
  
    SUBROUTINE interpolate(add_yes,g,coarse,fine,work)
  
  IMPLICIT NONE
  
  LOGICAL, INTENT(in   ) ::  add_yes
  INTEGER, INTENT(in   ) ::  g
  REAL(8)   , INTENT(inout) ::  coarse(b1L:(NN(1,g+1)+b1U),b2L:(NN(2,g+1)+b2U),b3L:(NN(3,g+1)+b3U))
  REAL(8)   , INTENT(inout) ::  fine  (b1L:(NN(1,g  )+b1U),b2L:(NN(2,g  )+b2U),b3L:(NN(3,g  )+b3U))
  REAL(8)   , INTENT(out  ) ::  work  (b1L:(NN(1,g  )+b1U),b2L:(NN(2,g  )+b2U),b3L:(NN(3,g  )+b3U))
  
  INTEGER                ::  i, di, imax, iimax, iiShift
  INTEGER                ::  j, dj, jmax, jjmax, jjShift
  INTEGER                ::  k, dk, kmax, kkmax, kkShift
  
  !******************************************************************
  INTEGER                ::  ii, jj, kk
  INTEGER                ::  dispg, offsg(1:3)
  REAL(8)   , ALLOCATABLE   ::  sendbuf(:,:,:)
  REAL(8)                   ::  recvbuf(1:(NN(1,g+1)+NB(1,g)-1)*(NN(2,g+1)+NB(2,g)-1)*(NN(3,g+1)+NB(3,g)-1)) ! TEST!!! Ist das richtig so?
  !******************************************************************
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - für allgemeine di, dj, dk geeignet                                                        !
  !              - di /= 1 <===> N1 /= 1                                                                     !
  !              - es wird nur in eine Richung ausgetauscht                                                  !
  !              - Null-Setzen am Rand nicht notwendig                                                       !
  !              - Es wird sequentiell über alle Raumrichtungen interpoliert, um keinen individuellen        !
  !                Interpolationsstencil für jeden Punkt im Raum speichern zu müssen.                        !
  !              - Durch das sequentielle Interpolieren kann der Interpolationsstencil klein und damit der   !
  !                Gesamtaufwand minimiert werden (Alternative: 8- bzw. 26-Punkt Stencil (!)). Nachteilig    !
  !                ist dabei das zusätzliche Arbeitsfeld auf dem feineren Gitterniveau (wird der Multigrid-  !
  !                Routine entliehen).                                               .                       !
  !              - Interpolationskoeffizienten werden auf dem jeweils feineren Gitter gespeichert, um nicht  !
  !                auf die entsprechenden Indizes des gröberen Gitters umrechnen zu müssen.                  !
  !              - Die Block-überlappenden Stirnflächen werden ebenfalls mitverarbeitet, aber eigentlich     !
  !                nicht gebraucht (erleichtert die Programmierung), so dass eigentlich eine Initialisierung !
  !                notwendig wäre. Dies wird jedoch zuvor schon in der korrespondierenden Restriktions-      !
  !                Routine erledigt, so dass dies hier nicht mehr notwendig ist.                             !
  !----------------------------------------------------------------------------------------------------------!
  
  
  ! ACHTUNG!!! Verwendung von parametrischen Strides di,dj,dk in Schleifen verhindert die Vektorisierung bzw.
  !            das Prefetching! Andererseits ist der Geschwindigkeitsgewinn nur sehr gering (schon getestet).
  !
  ! - Geschwindigkeit ist trotz Prefetching / Vektorisierung stark durch die Speicherzugriffszeit limitiert
  ! - Das wird z.B. deutlich bei Single- vs. Dualcorebetrieb
  
  
  imax    = NN(1,g)
  jmax    = NN(2,g)
  kmax    = NN(3,g)
  
  iimax   = (NN(1,g+1)-1)/n_gather(1,g+1)+1
  jjmax   = (NN(2,g+1)-1)/n_gather(2,g+1)+1
  kkmax   = (NN(3,g+1)-1)/n_gather(3,g+1)+1
  
  di      = (NN(1,g)-1)/(iimax-1)
  dj      = (NN(2,g)-1)/(jjmax-1)
  dk      = (NN(3,g)-1)/(kkmax-1)
  
  iiShift = (iimax-1)*MOD(iB(1,g)-1,n_gather(1,g+1))
  jjShift = (jjmax-1)*MOD(iB(2,g)-1,n_gather(2,g+1))
  kkShift = (kkmax-1)*MOD(iB(3,g)-1,n_gather(3,g+1))
  
  
  IF (BC(1,1,g+1) .GT. 0 .AND. BC(1,2,g+1) .GT. 0) coarse(1        ,1        ,1:NN(3,g+1)) = (coarse(2          ,1        ,1:NN(3,g+1)) + coarse(1        ,2          ,1:NN(3,g+1)) + coarse(2          ,2          ,1:NN(3,g+1)))/3.
  IF (BC(1,1,g+1) .GT. 0 .AND. BC(2,2,g+1) .GT. 0) coarse(1        ,NN(2,g+1),1:NN(3,g+1)) = (coarse(2          ,NN(2,g+1),1:NN(3,g+1)) + coarse(1        ,NN(2,g+1)-1,1:NN(3,g+1)) + coarse(2          ,NN(2,g+1)-1,1:NN(3,g+1)))/3.
  IF (BC(2,1,g+1) .GT. 0 .AND. BC(1,2,g+1) .GT. 0) coarse(NN(1,g+1),1        ,1:NN(3,g+1)) = (coarse(NN(1,g+1)-1,1        ,1:NN(3,g+1)) + coarse(NN(1,g+1),2          ,1:NN(3,g+1)) + coarse(NN(1,g+1)-1,2          ,1:NN(3,g+1)))/3.
  IF (BC(2,1,g+1) .GT. 0 .AND. BC(2,2,g+1) .GT. 0) coarse(NN(1,g+1),NN(2,g+1),1:NN(3,g+1)) = (coarse(NN(1,g+1)-1,NN(2,g+1),1:NN(3,g+1)) + coarse(NN(1,g+1),NN(2,g+1)-1,1:NN(3,g+1)) + coarse(NN(1,g+1)-1,NN(2,g+1)-1,1:NN(3,g+1)))/3.
  
  IF (BC(1,1,g+1) .GT. 0 .AND. BC(1,3,g+1) .GT. 0) coarse(1        ,1:NN(2,g+1),1        ) = (coarse(2          ,1:NN(2,g+1),1        ) + coarse(1        ,1:NN(2,g+1),2          ) + coarse(2          ,1:NN(2,g+1),2          ))/3.
  IF (BC(1,1,g+1) .GT. 0 .AND. BC(2,3,g+1) .GT. 0) coarse(1        ,1:NN(2,g+1),NN(3,g+1)) = (coarse(2          ,1:NN(2,g+1),NN(3,g+1)) + coarse(1        ,1:NN(2,g+1),NN(3,g+1)-1) + coarse(2          ,1:NN(2,g+1),NN(3,g+1)-1))/3.
  IF (BC(2,1,g+1) .GT. 0 .AND. BC(1,3,g+1) .GT. 0) coarse(NN(1,g+1),1:NN(2,g+1),1        ) = (coarse(NN(1,g+1)-1,1:NN(2,g+1),1        ) + coarse(NN(1,g+1),1:NN(2,g+1),2          ) + coarse(NN(1,g+1)-1,1:NN(2,g+1),2          ))/3.
  IF (BC(2,1,g+1) .GT. 0 .AND. BC(2,3,g+1) .GT. 0) coarse(NN(1,g+1),1:NN(2,g+1),NN(3,g+1)) = (coarse(NN(1,g+1)-1,1:NN(2,g+1),NN(3,g+1)) + coarse(NN(1,g+1),1:NN(2,g+1),NN(3,g+1)-1) + coarse(NN(1,g+1)-1,1:NN(2,g+1),NN(3,g+1)-1))/3.

  IF (BC(1,2,g+1) .GT. 0 .AND. BC(1,3,g+1) .GT. 0) coarse(1:NN(1,g+1),1        ,1        ) = (coarse(1:NN(1,g+1),2          ,1        ) + coarse(1:NN(1,g+1),1        ,2          ) + coarse(1:NN(1,g+1),2          ,2          ))/3.
  IF (BC(1,2,g+1) .GT. 0 .AND. BC(2,3,g+1) .GT. 0) coarse(1:NN(1,g+1),1        ,NN(3,g+1)) = (coarse(1:NN(1,g+1),2          ,NN(3,g+1)) + coarse(1:NN(1,g+1),1        ,NN(3,g+1)-1) + coarse(1:NN(1,g+1),2          ,NN(3,g+1)-1))/3.
  IF (BC(2,2,g+1) .GT. 0 .AND. BC(1,3,g+1) .GT. 0) coarse(1:NN(1,g+1),NN(2,g+1),1        ) = (coarse(1:NN(1,g+1),NN(2,g+1)-1,1        ) + coarse(1:NN(1,g+1),NN(2,g+1),2          ) + coarse(1:NN(1,g+1),NN(2,g+1)-1,2          ))/3.
  IF (BC(2,2,g+1) .GT. 0 .AND. BC(2,3,g+1) .GT. 0) coarse(1:NN(1,g+1),NN(2,g+1),NN(3,g+1)) = (coarse(1:NN(1,g+1),NN(2,g+1)-1,NN(3,g+1)) + coarse(1:NN(1,g+1),NN(2,g+1),NN(3,g+1)-1) + coarse(1:NN(1,g+1),NN(2,g+1)-1,NN(3,g+1)-1))/3.
  
  
  CALL exchange_relax(g+1,ex1,ex2,ex3,0,.FALSE.,coarse) ! Anmerkung: .FALSE. ist ok ...
  
  
  !***********************************************************************************************************
  IF (n_gather(1,g+1)*n_gather(2,g+1)*n_gather(3,g+1) .GT. 1) THEN
     
     IF (participate_yes(g+1)) THEN
        DO k = 1, n_gather(3,g+1)
           DO j = 1, n_gather(2,g+1)
              DO i = 1, n_gather(1,g+1)
                 
                 dispg      = dispI(    i+(j-1)*n_gather(1,g+1)+(k-1)*n_gather(1,g+1)*n_gather(2,g+1),g+1)
                 offsg(1:3) = offsI(1:3,i+(j-1)*n_gather(1,g+1)+(k-1)*n_gather(1,g+1)*n_gather(2,g+1),g+1)
                 
                 DO kk = 1, kkmax
                    DO jj = 1, jjmax
                       DO ii = 1, iimax
                          recvbuf(dispg+ii+(jj-1)*iimax+(kk-1)*iimax*jjmax) = coarse(ii+offsg(1),jj+offsg(2),kk+offsg(3))
                       END DO
                    END DO
                 END DO
                 
              END DO
           END DO
        END DO
     END IF
     
     
     ALLOCATE(sendbuf(1:iimax,1:jjmax,1:kkmax)) ! Anmerkung: Besser nicht fest allocieren um Speicherplatz zu sparen, ODER gleich "coarse" verwenden!
     
     CALL MPI_SCATTER(recvbuf,iimax*jjmax*kkmax,MPI_REAL8,sendbuf,iimax*jjmax*kkmax,MPI_REAL8,rankc2(g+1),comm2(g+1),merror)
     
     DO kk = 1, kkmax
        k = dk*(kk-1)+1
        DO jj = 1, jjmax
           j = dj*(jj-1)+1
!pgi$ unroll = n:8
           DO ii = 1, iimax
              i = di*(ii-1)+1
              work(i,j,k) = sendbuf(ii,jj,kk)
           END DO
        END DO
     END DO
     
     DEALLOCATE(sendbuf)
     
  ELSE
     
!pgi$ unroll = n:8
     work(1:imax:di,1:jmax:dj,1:kmax:dk) = coarse((1+iiShift):(iimax+iiShift),(1+jjShift):(jjmax+jjShift),(1+kkShift):(kkmax+kkShift))
     
  END IF
  !***********************************************************************************************************
  
  
  
  !===========================================================================================================
  IF (dk /= 1) THEN ! (dimens == 2) <==> (dk == 1) automatisch erfüllt!
     
     DO k = 2, kmax-1, dk
        DO j = 1, jmax, dj
!pgi$ unroll = n:8
           DO i = 1, imax, di
              work(i,j,k) = cI3(1,k,g)*work(i,j,k-1) + cI3(2,k,g)*work(i,j,k+1)
           END DO
        END DO
     END DO
     
  END IF
  !===========================================================================================================
  IF (dj /= 1) THEN ! TEST!!! in 2D wird hier doppelte Arbeit geleistet! (kmax == 2??)
     
     DO k = 1, kmax
        DO j = 2, jmax-1, dj
!pgi$ unroll = n:8
           DO i = 1, imax, di
              work(i,j,k) = cI2(1,j,g)*work(i,j-1,k) + cI2(2,j,g)*work(i,j+1,k)
           END DO
        END DO
     END DO
     
  END IF
  !===========================================================================================================
  IF (di /= 1) THEN
     
     DO k = 1, kmax
        DO j = 1, jmax
!pgi$ unroll = n:8
           DO i = 2, imax-1, di
              work(i,j,k) = cI1(1,i,g)*work(i-1,j,k) + cI1(2,i,g)*work(i+1,j,k)
           END DO
        END DO
     END DO
     
  END IF
  !===========================================================================================================
  
  
  IF (add_yes) THEN
!pgi$ unroll = n:8
     fine(1:imax,1:jmax,1:kmax) = fine(1:imax,1:jmax,1:kmax) + work(1:imax,1:jmax,1:kmax)
  END IF
  
  
  END SUBROUTINE interpolate
 
  
  
 
  
  
  
  
  
  
  
  
  SUBROUTINE interpolate_shreded(add_yes,g,coarse,fine,work)
  
  IMPLICIT NONE
  
  LOGICAL, INTENT(in   ) ::  add_yes
  INTEGER, INTENT(in   ) ::  g
  REAL(8)   , INTENT(inout) ::  coarse(b1L:(NN(1,g+1)+b1U),b2L:(NN(2,g+1)+b2U),b3L:(NN(3,g+1)+b3U))
  REAL(8)   , INTENT(inout) ::  fine  (b1L:(NN(1,g  )+b1U),b2L:(NN(2,g  )+b2U),b3L:(NN(3,g  )+b3U))
  REAL(8)   , INTENT(out  ) ::  work  (b1L:(NN(1,g  )+b1U),b2L:(NN(2,g  )+b2U),b3L:(NN(3,g  )+b3U))
  
  INTEGER                ::  i, di, imax, iimax, iiShift
  INTEGER                ::  j, dj, jmax, jjmax, jjShift
  INTEGER                ::  k, dk, kmax, kkmax, kkShift
  
  !******************************************************************
  INTEGER                ::  ii, jj, kk
  INTEGER                ::  dispg, offsg(1:3)
  REAL(8)   , ALLOCATABLE   ::  sendbuf(:,:,:)
  REAL(8)                   ::  recvbuf(1:(NN(1,g+1)+NB(1,g)-1)*(NN(2,g+1)+NB(2,g)-1)*(NN(3,g+1)+NB(3,g)-1)) ! TEST!!! Ist das richtig so?
  !******************************************************************
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - für allgemeine di, dj, dk geeignet                                                        !
  !              - di /= 1 <===> N1 /= 1                                                                     !
  !              - es wird nur in eine Richung ausgetauscht                                                  !
  !              - Null-Setzen am Rand nicht notwendig                                                       !
  !              - Es wird sequentiell über alle Raumrichtungen interpoliert, um keinen individuellen        !
  !                Interpolationsstencil für jeden Punkt im Raum speichern zu müssen.                        !
  !              - Durch das sequentielle Interpolieren kann der Interpolationsstencil klein und damit der   !
  !                Gesamtaufwand minimiert werden (Alternative: 8- bzw. 26-Punkt Stencil (!)). Nachteilig    !
  !                ist dabei das zusätzliche Arbeitsfeld auf dem feineren Gitterniveau (wird der Multigrid-  !
  !                Routine entliehen).                                               .                       !
  !              - Interpolationskoeffizienten werden auf dem jeweils feineren Gitter gespeichert, um nicht  !
  !                auf die entsprechenden Indizes des gröberen Gitters umrechnen zu müssen.                  !
  !              - Die Block-überlappenden Stirnflächen werden ebenfalls mitverarbeitet, aber eigentlich     !
  !                nicht gebraucht (erleichtert die Programmierung), so dass eigentlich eine Initialisierung !
  !                notwendig wäre. Dies wird jedoch zuvor schon in der korrespondierenden Restriktions-      !
  !                Routine erledigt, so dass dies hier nicht mehr notwendig ist.                             !
  !----------------------------------------------------------------------------------------------------------!
  
  
  ! ACHTUNG!!! Verwendung von parametrischen Strides di,dj,dk in Schleifen verhindert die Vektorisierung bzw.
  !            das Prefetching! Andererseits ist der Geschwindigkeitsgewinn nur sehr gering (schon getestet).
  !
  ! - Geschwindigkeit ist trotz Prefetching / Vektorisierung stark durch die Speicherzugriffszeit limitiert
  ! - Das wird z.B. deutlich bei Single- vs. Dualcorebetrieb
 
  imax    = NN(1,g)
  jmax    = NN(2,g)
  kmax    = NN(3,g)

  iimax   = (NN(1,g+1)-1)/n_gather(1,g+1)+1
  jjmax   = (NN(2,g+1)-1)/n_gather(2,g+1)+1
  kkmax   = (NN(3,g+1)-1)/n_gather(3,g+1)+1

  di      = (NN(1,g)-1)/(iimax-1)
  dj      = (NN(2,g)-1)/(jjmax-1)
  dk      = (NN(3,g)-1)/(kkmax-1)

  iiShift = (iimax-1)*MOD(iB(1,g)-1,n_gather(1,g+1))
  jjShift = (jjmax-1)*MOD(iB(2,g)-1,n_gather(2,g+1))
  kkShift = (kkmax-1)*MOD(iB(3,g)-1,n_gather(3,g+1))  
  !print*, g, BC(1,1,g),BC(1,2,g),BC(1,3,g), 'BC(1,1,g), BC(1,2,g), BC(1,3,g)' 
 
  CALL exchange_relax(g+1,ex1,ex2,ex3,0,.FALSE.,coarse) ! Anmerkung: .FALSE. ist ok ...
 
  !***********************************************************************************************************
  IF (n_gather(1,g+1)*n_gather(2,g+1)*n_gather(3,g+1) .GT. 1) THEN
     IF (participate_yes(g+1)) THEN
        DO k = 1, n_gather(3,g+1)
           DO j = 1, n_gather(2,g+1)
              DO i = 1, n_gather(1,g+1)
                 
                 dispg      = dispI(    i+(j-1)*n_gather(1,g+1)+(k-1)*n_gather(1,g+1)*n_gather(2,g+1),g+1)
                 offsg(1:3) = offsI(1:3,i+(j-1)*n_gather(1,g+1)+(k-1)*n_gather(1,g+1)*n_gather(2,g+1),g+1)
                 
                 DO kk = 1, kkmax
                    DO jj = 1, jjmax
                       DO ii = 1, iimax
                          recvbuf(dispg+ii+(jj-1)*iimax+(kk-1)*iimax*jjmax) = coarse(ii+offsg(1),jj+offsg(2),kk+offsg(3))
                       END DO
                    END DO
                 END DO
                 
              END DO
           END DO
        END DO
     END IF
    
    
     ALLOCATE(sendbuf(1:iimax,1:jjmax,1:kkmax)) ! Anmerkung: Besser nicht fest allocieren um Speicherplatz zu sparen, ODER gleich "coarse" verwenden!
     
     CALL MPI_SCATTER(recvbuf,iimax*jjmax*kkmax,MPI_REAL8,sendbuf,iimax*jjmax*kkmax,MPI_REAL8,rankc2(g+1),comm2(g+1),merror)
     
     DO kk = 1, kkmax
        k = dk*(kk-1)+1
        DO jj = 1, jjmax
           j = dj*(jj-1)+1
           DO ii = 1, iimax
              i = di*(ii-1)+1
              work(i,j,k) = sendbuf(ii,jj,kk)
           END DO
        END DO
     END DO
     
     DEALLOCATE(sendbuf)
     
  ELSE
    
    work(1:imax:di,1:jmax:dj,1:kmax:dk) = coarse((1+iiShift):(iimax+iiShift),(1+jjShift):(jjmax+jjShift),(1+kkShift):(kkmax+kkShift))
    
  END IF
  !***********************************************************************************************************
  
  
  !===========================================================================================================
     
     DO k = 2, kmax-1, dk
        DO j = 1, jmax, dj
!pgi$ unroll = n:8
           DO i = 1, imax, di
              work(i,j,k) = cI3(1,k,g)*work(i,j,k-1) + cI3(2,k,g)*work(i,j,k+1)
           END DO
        END DO
     END DO
     
  !===========================================================================================================
     
     DO k = 1, kmax
        DO j = 2, jmax-1, dj
!pgi$ unroll = n:8
           DO i = 1, imax, di
              work(i,j,k) = cI2(1,j,g)*work(i,j-1,k) + cI2(2,j,g)*work(i,j+1,k)
           END DO
        END DO
     END DO
     
  !===========================================================================================================
     
     DO k = 1, kmax
        DO j = 1, jmax
!pgi$ unroll = n:8
           DO i = 2, imax-1, di
              work(i,j,k) = cI1(1,i,g)*work(i-1,j,k) + cI1(2,i,g)*work(i+1,j,k)
           END DO
        END DO
     END DO
     
  !===========================================================================================================
  
  
  IF (add_yes) THEN
!pgi$ unroll = n:8
     fine(1:imax,1:jmax,1:kmax) = fine(1:imax,1:jmax,1:kmax) + work(1:imax,1:jmax,1:kmax)
  END IF
  
  
  END SUBROUTINE interpolate_shreded
  
  
  
  
  
  
  
    SUBROUTINE GPU_solver(eps,n_it_max,init_yes,SS1,SS2,SS3,NN1,NN2,NN3,bb,phi,problem_type,quiet_yes1,quiet_yes2,preconditioner)
  
  IMPLICIT NONE
  
  REAL(8)   , INTENT(in)    ::  eps      !< termination criterion 
  INTEGER, INTENT(in)    ::  n_it_max !< max iteration steps
  LOGICAL, INTENT(in)    ::  init_yes !< whether to initialize (what?)
  
  INTEGER, INTENT(in)    ::  SS1      !< lower index bound dimension 1
  INTEGER, INTENT(in)    ::  SS2      !< lower index bound dimension 2
  INTEGER, INTENT(in)    ::  SS3      !< lower index bound dimension 3
  
  INTEGER, INTENT(in)    ::  NN1      !< upper index bound dimension 1
  INTEGER, INTENT(in)    ::  NN2      !< upper index bound dimension 2
  INTEGER, INTENT(in)    ::  NN3      !< upper index bound dimension 3
  
  INTEGER, INTENT(in)    ::  preconditioner ! see above 
  INTEGER, INTENT(in)    ::  problem_type   !< problem type to be solved.
  
  REAL(8)   , INTENT(in)    ::  bb (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)) !< residuum 
  REAL(8)   , INTENT(inout) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)) !< solution vector
  REAL(8)                   ::  phi_old(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)) !< solution vector 
  REAL(8)                   ::  GPUres(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)) !< solution vector
  
  
  LOGICAL, INTENT(in)    ::  quiet_yes1 !< output suppression? 
  LOGICAL, INTENT(in)    ::  quiet_yes2 !< output suppression?
  
  INTEGER                ::  counter
  INTEGER                ::  i, j, k
  
  REAL(8)                   ::  norm2   !, norm2_global
  REAL(8)                   ::  norm_inf, norm_inf_global, norm_inf_prev
  REAL(8)                   ::  rhr, rhr_prev, ArAr, rAr, rhAp
  REAL(8)                   ::  scalar_global(1:2)
  REAL(8)                   ::  alpha, beta, omega
  
  LOGICAL                ::  exit_yes
  
  
  real(8)                   ::  t1, t2
 
  include 'mpif.h'

  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Gitterstreckung ist ein Problem bei der Gewichtung der Residuen, was beim Helmholtz-      !
  !                Problem aufgrund der Diagonaldominanz weit weniger ausgeprägt sein sollte als beim reinen !
  !                Poisson-Problem. Da der Multigrid-Glätter (Jacobi, Gauss-Seidel) davon generell unbeein-  !
  !                flusst sind, kann dieses Problem hier mehr oder weniger vernachlässigt werden. Alternativ !
  !                könnten die Stencils analog zu "weight" auch für das gesamte Feld gespeichert werden, was !
  !                allerdings einen sehr grossen Speicheraufwand und lange Speicherzugriffszeiten erwarten   !
  !                lassen, so dass hier eine evtl. geringere Konvergenzrate mit gutem Gewissen in Kauf       !
  !                genommen wird.                                                                            !
  !              - Es wird hier generell rh = rr angenommen, da somit immer <rh,rr> /= 0 gilt!               !
  !              - dAXPY ist nicht wirklich schneller als die direkt programmierten Operationen.             !
  !              - rhr muss nicht initialisiert werden!                                                      !
  !              - Um pp nur bei Bedarf berechnen zu müssen, wird dazu eine Fallunterscheidung nach erst     !
  !                hinter "CALL status_iteration" eingeführt.                                                !
  !              - z2 könnte zu Ungunsten der Geschwindigkeit eingespaart werden.                            !
  !----------------------------------------------------------------------------------------------------------!
  
  
  ! - Geschwindigkeit ist trotz Prefetching / Vektorisierung stark durch die Speicherzugriffszeit limitiert
  ! - Das wird z.B. deutlich bei Single- vs. Dualcorebetrieb
  
  
  !===========================================================================================================
  IF (.NOT. init_yes) THEN
     IF (problem_type == 2) CALL product_div_grad (phi,rr)
     !-----------------------------------------------------------------------------------------------------------
     ! bbecsek 150106: problem_type == 3 has been removed (associated with concentrations)
     !-----------------------------------------------------------------------------------------------------------
     
     rr(SS1:NN1,SS2:NN2,SS3:NN3) = bb(SS1:NN1,SS2:NN2,SS3:NN3) - rr(SS1:NN1,SS2:NN2,SS3:NN3)
     rh(SS1:NN1,SS2:NN2,SS3:NN3) = rr(SS1:NN1,SS2:NN2,SS3:NN3)
  END IF
  !===========================================================================================================
  
  !===========================================================================================================
  !=== Residuum ==============================================================================================
  !===========================================================================================================
  IF (init_yes) THEN
     CALL get_norms(SS1,SS2,SS3,NN1,NN2,NN3,bb,problem_type,.TRUE.,.TRUE. ,normInf=norm_inf,normTwo=norm2)
  ELSE
     CALL get_norms(SS1,SS2,SS3,NN1,NN2,NN3,rr,problem_type,.TRUE.,.FALSE.,normInf=norm_inf)
  END IF
  norm_inf_prev = norm_inf
  !===========================================================================================================
  
  counter = 0
  
  
  ITERATE: DO
     
     !========================================================================================================
     !=== Überprüfen des Konvergenzkriteriums ================================================================
     !========================================================================================================
!     CALL status_iteration(eps,norm_inf,counter,n_it_max,exit_yes,quiet_yes1,quiet_yes2)
     
   
     IF (exit_yes .AND. counter == 0 .AND. init_yes) phi(SS1:NN1,SS2:NN2,SS3:NN3) = 0.
     !========================================================================================================
     
     
     !========================================================================================================
     !=== nächster Durchlauf =================================================================================
     !========================================================================================================
     counter = counter + 1
     !========================================================================================================
     
     
    
    CALL exchange(1,1,phi(:,:,:))
    CALL exchange(2,2,phi(:,:,:))
    CALL exchange(3,3,phi(:,:,:))
   ! call exchange_all_all(.TRUE.,bb)

    call device_jacobi_mg_relaxation(1,.7,phi, bb, phi) 

     ! ACHTUNG!!! Zu teuer im Vergleich zur obigen Variante:
     !CALL product_scalar(SS1,SS2,SS3,NN1,NN2,NN3,Ar,rr,rAr)
     !CALL product_scalar(SS1,SS2,SS3,NN1,NN2,NN3,Ar,Ar,ArAr)
     !========================================================================================================
    !========================================================================================================
     norm_inf = 0.
     !--------------------------------------------------------------------------------------------------------
!    print*, counter, 'counter'


 !       DO k = SS3, NN3
 !          DO j = SS2, NN2
 !             DO i = SS1, NN1
 !                GPUres(i,j,k) = bb(i,j,k)                                     &
 !                            &      - cdg1(-1,i,1)*phi(i-1,j,k) - cdg1(1,i,1)*phi(i+1,j,k)     &
 !                            &      - cdg2(-1,j,1)*phi(i,j-1,k) - cdg2(1,j,1)*phi(i,j+1,k)     &
 !                            &      - cdg3(-1,k,1)*phi(i,j,k-1) - cdg3(1,k,1)*phi(i,j,k+1)    &
 !                            &      - (cdg1(0,i,1) + cdg2(0,j,1) + cdg3(0,k,1))*phi(i,j,k)
 !                norm_inf = MAX(ABS(GPUres(i,j,k)),norm_inf)
 !             END DO
 !          END DO
 !       END DO



              DO k = SS3, NN3
                 DO j = SS2, NN2
                    DO i = SS1, NN1
                       norm_inf = MAX(ABS((phi(i,j,k)-phi_old(i,j,k))),norm_inf)
                       phi_old(i,j,k) = phi(i,j,k)
                    END DO
                 END DO
              END DO
    !--------------------------------------------------------------------------------------------------------
     CALL MPI_ALLREDUCE(norm_inf,norm_inf_global,1,MPI_REAL8,MPI_MAX,COMM_CART,merror) ! MPI_REDUCE bringt nichts, weil exit_yes dann mit MPI_BCAST verteilt werden müsste ...
     norm_inf = norm_inf_global
     !========================================================================================================
     if (norm_inf .le. 1.e-6) EXIT ITERATE  
     if ((rank ==0) .and. (mod(counter,50)==0))  print*,counter, norm_inf
!     if (rank ==0)  print*, counter
     
     !========================================================================================================
     !=== Konvergenzstatistik ================================================================================
     !========================================================================================================
     ratioP(substep,number_poisson) = ratioP(substep,number_poisson) + LOG10(norm_inf/norm_inf_prev)
     
     countP(substep,number_poisson) = countP(substep,number_poisson) + 1
     
     norm_inf_prev = norm_inf
     !========================================================================================================
     
  END DO ITERATE
  
  
  END SUBROUTINE GPU_solver
  
  
  
  
 
  
  
  
   !> subroutine for solving LDEs with the iterative bi-conjugate gradient stabilized method.
  !! The subroutine is used both for the preconditioned pressure iteration (Poisson problem) as well as for the 
  !! velocity iteration (Helmholtz problems).
  !! The routine is only called with problem_type={1,2,3,4} unlike BiCGstab2 which is called only with problem_type=5.
  !! @param[in] preconditioner referring to the preconditioner for the Krylov-subspace methods (biCGstab), not the preconditioner
  !!        of the LDE matrices. 0=:no preconditioning, 1:=V-cycle MG smoothing, 2:=F-cycle MG smoothing.
  !! @note: why aren't BiCGstab and BiCGstab2 combined into a single routine?
   SUBROUTINE BiCGstab_orig(eps,n_it_max,init_yes,SS1,SS2,SS3,NN1,NN2,NN3,bb,phi,problem_type,quiet_yes1,quiet_yes2,preconditioner)
  
  IMPLICIT NONE
  
  REAL(8)   , INTENT(in)    ::  eps      !< termination criterion 
  INTEGER, INTENT(in)    ::  n_it_max !< max iteration steps
  LOGICAL, INTENT(in)    ::  init_yes !< whether to initialize (what?)
  
  INTEGER, INTENT(in)    ::  SS1      !< lower index bound dimension 1
  INTEGER, INTENT(in)    ::  SS2      !< lower index bound dimension 2
  INTEGER, INTENT(in)    ::  SS3      !< lower index bound dimension 3
  
  INTEGER, INTENT(in)    ::  NN1      !< upper index bound dimension 1
  INTEGER, INTENT(in)    ::  NN2      !< upper index bound dimension 2
  INTEGER, INTENT(in)    ::  NN3      !< upper index bound dimension 3
  
  INTEGER, INTENT(in)    ::  preconditioner ! see above 
  INTEGER, INTENT(in)    ::  problem_type   !< problem type to be solved.
  
  REAL(8)   , INTENT(in)    ::  bb (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)) !< residuum 
  REAL(8)   , INTENT(inout) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)) !< solution vector
  REAL(8)                   ::  phi_old(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)) !< solution vector
  
  

  LOGICAL, INTENT(in)    ::  quiet_yes1 !< output suppression? 
  LOGICAL, INTENT(in)    ::  quiet_yes2 !< output suppression?
  
  INTEGER                ::  counter
  INTEGER                ::  i, j, k
  
  REAL(8)                   ::  norm2   !, norm2_global
  REAL(8)                   ::  norm_inf, norm_inf_global, norm_lap, norm_lap_global,norm_inf2, norm_inf_global2, norm_inf_prev
  REAL(8)                   ::  rhr, rhr_prev, ArAr, rAr, rhAp
  REAL(8)                   ::  scalar_global(1:2)
  REAL(8)                   ::  alpha, beta, omega
  
  LOGICAL                ::  exit_yes
  
  
  real(8)                   ::  t1, t2
 
  include 'mpif.h'

  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Gitterstreckung ist ein Problem bei der Gewichtung der Residuen, was beim Helmholtz-      !
  !                Problem aufgrund der Diagonaldominanz weit weniger ausgeprägt sein sollte als beim reinen !
  !                Poisson-Problem. Da der Multigrid-Glätter (Jacobi, Gauss-Seidel) davon generell unbeein-  !
  !                flusst sind, kann dieses Problem hier mehr oder weniger vernachlässigt werden. Alternativ !
  !                könnten die Stencils analog zu "weight" auch für das gesamte Feld gespeichert werden, was !
  !                allerdings einen sehr grossen Speicheraufwand und lange Speicherzugriffszeiten erwarten   !
  !                lassen, so dass hier eine evtl. geringere Konvergenzrate mit gutem Gewissen in Kauf       !
  !                genommen wird.                                                                            !
  !              - Es wird hier generell rh = rr angenommen, da somit immer <rh,rr> /= 0 gilt!               !
  !              - dAXPY ist nicht wirklich schneller als die direkt programmierten Operationen.             !
  !              - rhr muss nicht initialisiert werden!                                                      !
  !              - Um pp nur bei Bedarf berechnen zu müssen, wird dazu eine Fallunterscheidung nach erst     !
  !                hinter "CALL status_iteration" eingeführt.                                                !
  !              - z2 könnte zu Ungunsten der Geschwindigkeit eingespaart werden.                            !
  !----------------------------------------------------------------------------------------------------------!
  
  
  ! - Geschwindigkeit ist trotz Prefetching / Vektorisierung stark durch die Speicherzugriffszeit limitiert
  ! - Das wird z.B. deutlich bei Single- vs. Dualcorebetrieb
  
  
  !===========================================================================================================
  IF (.NOT. init_yes) THEN
     IF (problem_type == 2) CALL product_div_grad (phi,rr)
     !-----------------------------------------------------------------------------------------------------------
     ! bbecsek 150106: problem_type == 3 has been removed (associated with concentrations)
     !-----------------------------------------------------------------------------------------------------------
     
     rr(SS1:NN1,SS2:NN2,SS3:NN3) = bb(SS1:NN1,SS2:NN2,SS3:NN3) - rr(SS1:NN1,SS2:NN2,SS3:NN3)
     rh(SS1:NN1,SS2:NN2,SS3:NN3) = rr(SS1:NN1,SS2:NN2,SS3:NN3)
  END IF
  !===========================================================================================================
  
  !===========================================================================================================
  !=== Residuum ==============================================================================================
  !===========================================================================================================
  IF (init_yes) THEN
     CALL get_norms(SS1,SS2,SS3,NN1,NN2,NN3,bb,problem_type,.TRUE.,.TRUE. ,normInf=norm_inf,normTwo=norm2)
  ELSE
     CALL get_norms(SS1,SS2,SS3,NN1,NN2,NN3,rr,problem_type,.TRUE.,.FALSE.,normInf=norm_inf)
  END IF
  norm_inf_prev = norm_inf
  !===========================================================================================================
  
  counter = 0
  
  
  ITERATE: DO
     
     !========================================================================================================
     !=== Überprüfen des Konvergenzkriteriums ================================================================
     !========================================================================================================
     CALL status_iteration(eps,norm_inf,counter,n_it_max,exit_yes,quiet_yes1,quiet_yes2)
     
     ! soll Konvergenz sicherstellen (oder offsetprec runtersetzen ...).
     !    norm_inf == 0. << eps ist ein Sonderfall, bei dem es keinen Sinn macht, eine weitere Iteration
     !    zu rechnen. 
     IF (problem_type == 1 .AND.                           counter .LT. 1 .AND. norm_inf /= 0.) exit_yes = .FALSE.
     IF (problem_type == 2 .AND. number_poisson == 1 .AND. counter .LT. 1 .AND. norm_inf /= 0.) exit_yes = .FALSE.
     IF (problem_type == 2 .AND. number_poisson == 2 .AND. counter .LT. 0 .AND. norm_inf /= 0.) exit_yes = .FALSE.
     
     IF (exit_yes .AND. counter == 0 .AND. init_yes) phi(SS1:NN1,SS2:NN2,SS3:NN3) = 0.
     IF (exit_yes) EXIT ITERATE
     !========================================================================================================
     
     
     !========================================================================================================
     !=== nächster Durchlauf =================================================================================
     !========================================================================================================
     counter = counter + 1
     !========================================================================================================
     
     
     !========================================================================================================
     rhr_prev = rhr
     IF (init_yes) THEN
        IF (counter == 1) THEN
           rhr = norm2
        ELSE
           CALL product_scalar(SS1,SS2,SS3,NN1,NN2,NN3,rr,bb,rhr)
        END IF
     ELSE
        CALL product_scalar(SS1,SS2,SS3,NN1,NN2,NN3,rr,rh,rhr)
     END IF
     !========================================================================================================
     IF (ABS(rhr) == 0.) THEN
        IF (rank == 0) WRITE(* ,'(a,E13.5)') 'rhr =', rhr
        IF (rank == 0) WRITE(10,'(a,E13.5)') 'rhr =', rhr
        
        rh(SS1:NN1,SS2:NN2,SS3:NN3) = rr(SS1:NN1,SS2:NN2,SS3:NN3) ! Neuer Referenzvektor ...
        CALL product_scalar(SS1,SS2,SS3,NN1,NN2,NN3,rr,rh,rhr)
     END IF
     !========================================================================================================
     IF (counter .GE. 2) THEN
        IF (omega == 0.) THEN
           IF (rank == 0) WRITE(* ,'(a,E13.5)') 'omega =', omega
           IF (rank == 0) WRITE(10,'(a,E13.5)') 'omega =', omega
           IF (problem_type == 2 .OR. problem_type == 4) THEN
              EXIT ITERATE
           ELSE
              CALL MPI_FINALIZE(merror)
              STOP
           END IF
        END IF
        IF (rhr_prev == 0.) THEN
           IF (rank == 0) WRITE(* ,'(a,E13.5)') 'rhr_prev =', rhr_prev
           IF (rank == 0) WRITE(10,'(a,E13.5)') 'rhr_prev =', rhr_prev
           IF (problem_type == 2 .OR. problem_type == 4) THEN
              EXIT ITERATE
           ELSE
              CALL MPI_FINALIZE(merror)
              STOP
           END IF
        END IF
        beta = (alpha/omega)*(rhr/rhr_prev)
        omega = -beta*omega
        pp(SS1:NN1,SS2:NN2,SS3:NN3) = rr(SS1:NN1,SS2:NN2,SS3:NN3) + beta*pp(SS1:NN1,SS2:NN2,SS3:NN3) + omega*Ap(SS1:NN1,SS2:NN2,SS3:NN3)
     ELSE
        IF (init_yes .AND. counter == 1) THEN
           pp(SS1:NN1,SS2:NN2,SS3:NN3) = bb(SS1:NN1,SS2:NN2,SS3:NN3)
        ELSE
           pp(SS1:NN1,SS2:NN2,SS3:NN3) = rr(SS1:NN1,SS2:NN2,SS3:NN3)
        END IF
     END IF
     !========================================================================================================
     IF (preconditioner == 0) THEN
        IF (problem_type == 2) CALL product_div_grad (pp,Ap)
        !-----------------------------------------------------------------------------------------------------------
        ! bbecsek 150106: problem_type == 3 has been removed (associated with concentrations)
        !-----------------------------------------------------------------------------------------------------------
     ELSE IF (preconditioner == 1 .OR. preconditioner == 2) THEN
        IF (preconditioner == 1) then
              !call MPI_Barrier(MPI_COMM_WORLD, merror)
              !t1 =   MPI_Wtime()                    
                    CALL multigridV(.TRUE.,1,pp,z1,problem_type)
              !call MPI_Barrier(MPI_COMM_WORLD, merror)
              !t2 =   MPI_Wtime()
              !print*, 'multigrid V cycle time', t1, t2, t2-t1          
        endif
        IF (problem_type == 2) CALL product_div_grad (z1,Ap)
        !-----------------------------------------------------------------------------------------------------------
        ! bbecsek 150106: problem_type == 3 has been removed (associated with concentrations)
        !-----------------------------------------------------------------------------------------------------------
     ELSE
        IF (rank == 0) WRITE(*,'(a)') 'ERROR! Specify valid preconditioner!'
        CALL MPI_FINALIZE(merror)
        STOP
     END IF
     !========================================================================================================
     IF (init_yes) THEN
        CALL product_scalar(SS1,SS2,SS3,NN1,NN2,NN3,Ap,bb,rhAp)
     ELSE
        CALL product_scalar(SS1,SS2,SS3,NN1,NN2,NN3,Ap,rh,rhAp)
     END IF
     !========================================================================================================
     IF (ABS(rhAp) == 0.) THEN
        IF (rank == 0) WRITE(* ,'(a,E13.5)') 'rhAp =', rhAp
        IF (rank == 0) WRITE(10,'(a,E13.5)') 'rhAp =', rhAp
        IF (problem_type == 2 .OR. problem_type == 4) THEN
           EXIT ITERATE
        ELSE
           IF (rhr /= 0.) THEN
              CALL MPI_FINALIZE(merror)
              STOP
           END IF
        END IF
        alpha = 0.
     ELSE
        alpha = rhr / rhAp
     END IF
     !========================================================================================================
     IF (init_yes .AND. counter == 1) THEN
        DO k = SS3, NN3
           DO j = SS2, NN2
!pgi$ unroll = n:8
              DO i = SS1, NN1
                 rr(i,j,k) = bb(i,j,k) - alpha*Ap(i,j,k)
              END DO
           END DO
        END DO
     ELSE
        DO k = SS3, NN3
           DO j = SS2, NN2
!pgi$ unroll = n:8
              DO i = SS1, NN1
                 rr(i,j,k) = rr(i,j,k) - alpha*Ap(i,j,k)
              END DO
           END DO
        END DO
     END IF
     !========================================================================================================
     IF (preconditioner == 0) THEN
     IF (problem_type == 2) CALL product_div_grad (rr,Ar)
        !-----------------------------------------------------------------------------------------------------------
        ! bbecsek 150106: problem_type == 3 has been removed (associated with concentrations)
        !-----------------------------------------------------------------------------------------------------------
     ELSE IF (preconditioner == 1 .OR. preconditioner == 2) THEN
        IF (preconditioner == 1) CALL multigridV(.TRUE.,1,rr,z2,problem_type)
        IF (problem_type == 2) CALL product_div_grad (z2,Ar)
        !-----------------------------------------------------------------------------------------------------------
        ! bbecsek 150106: problem_type == 3 has been removed (associated with concentrations)
        !-----------------------------------------------------------------------------------------------------------
     END IF
     !========================================================================================================
     rAr  = 0.
     ArAr = 0.
     DO k = SS3, NN3
        DO j = SS2, NN2
!pgi$ unroll = n:8
           DO i = SS1, NN1
              rAr  = rAr  + rr(i,j,k)*Ar(i,j,k)
              ArAr = ArAr + Ar(i,j,k)*Ar(i,j,k)
           END DO
        END DO
     END DO
     
     CALL MPI_ALLREDUCE((/rAr,ArAr/),scalar_global,2,MPI_REAL8,MPI_SUM,COMM_CART,merror)
     rAr  = scalar_global(1)
     ArAr = scalar_global(2)
     
     ! ACHTUNG!!! Zu teuer im Vergleich zur obigen Variante:
     !CALL product_scalar(SS1,SS2,SS3,NN1,NN2,NN3,Ar,rr,rAr)
     !CALL product_scalar(SS1,SS2,SS3,NN1,NN2,NN3,Ar,Ar,ArAr)
     !========================================================================================================
     IF (ABS(rAr) == 0.) THEN
        IF (rank == 0) WRITE(* ,'(a,E13.5)') 'rAr =', rAr
        IF (rank == 0) WRITE(10,'(a,E13.5)') 'rAr =', rAr
     END IF
     IF (ABS(ArAr) == 0.) THEN
        IF (rank == 0) WRITE(* ,'(a,E13.5)') 'ArAr =', ArAr
        IF (rank == 0) WRITE(10,'(a,E13.5)') 'ArAr =', ArAr
        IF (problem_type == 2 .OR. problem_type == 4) THEN
           EXIT ITERATE
        ELSE
           IF (rAr /= 0.) THEN
              CALL MPI_FINALIZE(merror)
              STOP
           END IF
        END IF
        omega = 0.
     ELSE
        omega = rAr / ArAr
     END IF
     !========================================================================================================
     norm_inf = 0.
     norm_inf2= 0.
     !--------------------------------------------------------------------------------------------------------
!    print*, counter, 'counter'
     IF (counter == 1 .AND. init_yes) THEN
        IF (preconditioner == 0) THEN
           IF (problem_type == 2 .AND. weighting_yes) THEN
              DO k = SS3, NN3
                 DO j = SS2, NN2
!pgi$ unroll = n:8
                    DO i = SS1, NN1
                       phi(i,j,k) = alpha*pp(i,j,k) + omega*rr(i,j,k)
                       rr (i,j,k) =       rr(i,j,k) - omega*Ar(i,j,k)
                       norm_inf = MAX(ABS(rr(i,j,k)*weight(i,j,k)),norm_inf)
                    END DO
                 END DO
              END DO
           ELSE
              DO k = SS3, NN3
                 DO j = SS2, NN2
!pgi$ unroll = n:8
                    DO i = SS1, NN1
                       phi(i,j,k) = alpha*pp(i,j,k) + omega*rr(i,j,k)
                       rr (i,j,k) =       rr(i,j,k) - omega*Ar(i,j,k)
                       norm_inf = MAX(ABS(rr(i,j,k)),norm_inf)
                    END DO
                 END DO
              END DO
           END IF
        !-----------------------------------------------------------------------------------------------------
        ELSE
           IF (problem_type == 2 .AND. weighting_yes) THEN
              DO k = SS3, NN3
                 DO j = SS2, NN2
!pgi$ unroll = n:8
                    DO i = SS1, NN1
                       phi(i,j,k) = alpha*z1(i,j,k) + omega*z2(i,j,k)
                       rr (i,j,k) =       rr(i,j,k) - omega*Ar(i,j,k)
                       norm_inf = MAX(ABS(rr(i,j,k)*weight(i,j,k)),norm_inf)
                    END DO
                 END DO
              END DO
           ELSE
              print*, 'first loop'
              DO k = SS3, NN3
                 DO j = SS2, NN2
!pgi$ unroll = n:8
                    DO i = SS1, NN1
                       phi(i,j,k) = alpha*z1(i,j,k) + omega*z2(i,j,k)
                       rr (i,j,k) =       rr(i,j,k) - omega*Ar(i,j,k)
                       norm_inf = MAX(ABS(rr(i,j,k)),norm_inf)
                    END DO
                 END DO
              END DO
           END IF
        END IF
     !--------------------------------------------------------------------------------------------------------
     ELSE
        IF (preconditioner == 0) THEN
           IF (problem_type == 2 .AND. weighting_yes) THEN
              DO k = SS3, NN3
                 DO j = SS2, NN2
!pgi$ unroll = n:8
                    DO i = SS1, NN1
                       phi(i,j,k) = phi(i,j,k) + omega*rr(i,j,k) + alpha*pp(i,j,k)
                       rr (i,j,k) = rr (i,j,k) - omega*Ar(i,j,k)
                       norm_inf = MAX(ABS(rr(i,j,k)*weight(i,j,k)),norm_inf)
                    END DO
                 END DO
              END DO
           ELSE
              DO k = SS3, NN3
                 DO j = SS2, NN2
!pgi$ unroll = n:8
                    DO i = SS1, NN1
                       phi(i,j,k) = phi(i,j,k) + omega*rr(i,j,k) + alpha*pp(i,j,k)
                       rr (i,j,k) = rr (i,j,k) - omega*Ar(i,j,k)
                       norm_inf = MAX(ABS(rr(i,j,k)),norm_inf)
                    END DO
                 END DO
              END DO
           END IF
        !-----------------------------------------------------------------------------------------------------
        ELSE
           IF (problem_type == 2 .AND. weighting_yes) THEN
              DO k = SS3, NN3
                 DO j = SS2, NN2
!pgi$ unroll = n:8
                    DO i = SS1, NN1
                       phi(i,j,k) = phi(i,j,k) + omega*z2(i,j,k) + alpha*z1(i,j,k)
                       rr (i,j,k) = rr (i,j,k) - omega*Ar(i,j,k)
                       norm_inf = MAX(ABS(rr(i,j,k)*weight(i,j,k)),norm_inf)
                    END DO
                 END DO
              END DO
           ELSE
              DO k = SS3, NN3
                 DO j = SS2, NN2
!pgi$ unroll = n:8
                    DO i = SS1, NN1
                       phi(i,j,k) = phi(i,j,k) + omega*z2(i,j,k) + alpha*z1(i,j,k)
                       rr (i,j,k) = rr (i,j,k) - omega*Ar(i,j,k)
                       norm_inf = MAX(ABS(rr(i,j,k)),norm_inf)
 !                      norm_inf2 = MAX(ABS(phi(i,j,k)-phi_old(i,j,k)),norm_inf2)
 !                      phi_old(i,j,k) = phi(i,j,k)
                    END DO
                 END DO
              END DO
           END IF
        END IF
     END IF


     !--------------------------------------------------------------------------------------------------------
     CALL MPI_ALLREDUCE(norm_inf,norm_inf_global,1,MPI_REAL8,MPI_MAX,COMM_CART,merror) ! MPI_REDUCE bringt nichts, weil exit_yes dann mit MPI_BCAST verteilt werden müsste ...
     norm_inf = norm_inf_global
!     CALL MPI_ALLREDUCE(norm_inf2,norm_inf_global2,1,MPI_REAL8,MPI_MAX,COMM_CART,merror) ! MPI_REDUCE bringt nichts, weil exit_yes dann mit MPI_BCAST verteilt werden müsste ...
!     norm_inf2 = norm_inf_global2
     !========================================================================================================
!     if (rank==1) print*, norm_inf2

     
     !========================================================================================================
     !=== Konvergenzstatistik ================================================================================
     !========================================================================================================
     IF (problem_type == 1) ratioH(substep,direction     ) = ratioH(substep,direction     ) + LOG10(norm_inf/norm_inf_prev)
     IF (problem_type == 2) ratioP(substep,number_poisson) = ratioP(substep,number_poisson) + LOG10(norm_inf/norm_inf_prev)
     
     IF (problem_type == 1) countH(substep,direction     ) = countH(substep,direction     ) + 1
     IF (problem_type == 2) countP(substep,number_poisson) = countP(substep,number_poisson) + 1
     
     norm_inf_prev = norm_inf
     !========================================================================================================
     
  END DO ITERATE
  
  
  END SUBROUTINE BiCGstab_orig
  
  
  
  
  
  
  
  !> subroutine for solving LDEs with the iterative bi-conjugate gradient stabilized method.
  !! The subroutine is used both for the preconditioned pressure iteration (Poisson problem) as well as for the 
  !! velocity iteration (Helmholtz problems).
  !! The routine is only called with problem_type={1,2,3,4} unlike BiCGstab2 which is called only with problem_type=5.
  !! @param[in] preconditioner referring to the preconditioner for the Krylov-subspace methods (biCGstab), not the preconditioner
  !!        of the LDE matrices. 0=:no preconditioning, 1:=V-cycle MG smoothing, 2:=F-cycle MG smoothing.
  !! @note: why aren't BiCGstab and BiCGstab2 combined into a single routine?
   SUBROUTINE BiCGstab(eps,n_it_max,init_yes,SS1,SS2,SS3,NN1,NN2,NN3,bb,phi,problem_type,quiet_yes1,quiet_yes2,preconditioner)
  
  IMPLICIT NONE
  
  REAL(8)   , INTENT(in)    ::  eps      !< termination criterion 
  INTEGER, INTENT(in)    ::  n_it_max !< max iteration steps
  LOGICAL, INTENT(in)    ::  init_yes !< whether to initialize (what?)
  
  INTEGER, INTENT(in)    ::  SS1      !< lower index bound dimension 1
  INTEGER, INTENT(in)    ::  SS2      !< lower index bound dimension 2
  INTEGER, INTENT(in)    ::  SS3      !< lower index bound dimension 3
  
  INTEGER, INTENT(in)    ::  NN1      !< upper index bound dimension 1
  INTEGER, INTENT(in)    ::  NN2      !< upper index bound dimension 2
  INTEGER, INTENT(in)    ::  NN3      !< upper index bound dimension 3
  
  INTEGER, INTENT(in)    ::  preconditioner ! see above 
  INTEGER, INTENT(in)    ::  problem_type   !< problem type to be solved.
  
  REAL(8)   , INTENT(in)    ::  bb (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)) !< residuum 
  REAL(8)   , INTENT(inout) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)) !< solution vector
  
  REAL(8)                   ::  GPUres(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)) !< solution vector
  
  LOGICAL, INTENT(in)    ::  quiet_yes1 !< output suppression? 
  LOGICAL, INTENT(in)    ::  quiet_yes2 !< output suppression?
  
  INTEGER                ::  counter
  INTEGER                ::  i, j, k
  
  REAL(8)                   ::  norm2   !, norm2_global
  REAL(8)                   ::  norm_inf, norm_inf_global,norm_lap, norm_lap_global, norm_inf2, norm_inf_global2, norm_inf_prev, norm_zero
  REAL(8)                   ::  rhr, rhr_prev, ArAr, rAr, rhAp
  REAL(8)                   ::  scalar_global(1:2)
  REAL(8)                   ::  alpha, beta, omega
  
  LOGICAL                ::  exit_yes
  
  
  real(8)                   ::  t1, t2
 
  include 'mpif.h'

  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Gitterstreckung ist ein Problem bei der Gewichtung der Residuen, was beim Helmholtz-      !
  !                Problem aufgrund der Diagonaldominanz weit weniger ausgeprägt sein sollte als beim reinen !
  !                Poisson-Problem. Da der Multigrid-Glätter (Jacobi, Gauss-Seidel) davon generell unbeein-  !
  !                flusst sind, kann dieses Problem hier mehr oder weniger vernachlässigt werden. Alternativ !
  !                könnten die Stencils analog zu "weight" auch für das gesamte Feld gespeichert werden, was !
  !                allerdings einen sehr grossen Speicheraufwand und lange Speicherzugriffszeiten erwarten   !
  !                lassen, so dass hier eine evtl. geringere Konvergenzrate mit gutem Gewissen in Kauf       !
  !                genommen wird.                                                                            !
  !              - Es wird hier generell rh = rr angenommen, da somit immer <rh,rr> /= 0 gilt!               !
  !              - dAXPY ist nicht wirklich schneller als die direkt programmierten Operationen.             !
  !              - rhr muss nicht initialisiert werden!                                                      !
  !              - Um pp nur bei Bedarf berechnen zu müssen, wird dazu eine Fallunterscheidung nach erst     !
  !                hinter "CALL status_iteration" eingeführt.                                                !
  !              - z2 könnte zu Ungunsten der Geschwindigkeit eingespaart werden.                            !
  !----------------------------------------------------------------------------------------------------------!
  
  
  ! - Geschwindigkeit ist trotz Prefetching / Vektorisierung stark durch die Speicherzugriffszeit limitiert
  ! - Das wird z.B. deutlich bei Single- vs. Dualcorebetrieb
  
  
  !===========================================================================================================
  IF (.NOT. init_yes) THEN
     IF (problem_type == 2) CALL product_div_grad (phi,rr)
     !-----------------------------------------------------------------------------------------------------------
     ! bbecsek 150106: problem_type == 3 has been removed (associated with concentrations)
     !-----------------------------------------------------------------------------------------------------------
     
     rr(SS1:NN1,SS2:NN2,SS3:NN3) = bb(SS1:NN1,SS2:NN2,SS3:NN3) - rr(SS1:NN1,SS2:NN2,SS3:NN3)
     rh(SS1:NN1,SS2:NN2,SS3:NN3) = rr(SS1:NN1,SS2:NN2,SS3:NN3)
  END IF
  !===========================================================================================================
  
  !===========================================================================================================
  !=== Residuum ==============================================================================================
  !===========================================================================================================
  IF (init_yes) THEN
     CALL get_norms(SS1,SS2,SS3,NN1,NN2,NN3,bb,problem_type,.TRUE.,.TRUE. ,normInf=norm_inf,normTwo=norm2)
  ELSE
     CALL get_norms(SS1,SS2,SS3,NN1,NN2,NN3,rr,problem_type,.TRUE.,.FALSE.,normInf=norm_inf)
  END IF
  norm_inf_prev = norm_inf
  !===========================================================================================================
  
  counter = 0
  norm_zero =0.  
  
  ITERATE: DO
  

   
     !========================================================================================================
     !=== Überprüfen des Konvergenzkriteriums ================================================================
     !========================================================================================================
     CALL status_iteration(eps,norm_inf,counter,n_it_max,exit_yes,quiet_yes1,quiet_yes2)
     
     ! soll Konvergenz sicherstellen (oder offsetprec runtersetzen ...).
     !    norm_inf == 0. << eps ist ein Sonderfall, bei dem es keinen Sinn macht, eine weitere Iteration
     !    zu rechnen. 
     IF (problem_type == 1 .AND.                           counter .LT. 1 .AND. norm_inf /= 0.) exit_yes = .FALSE.
     IF (problem_type == 2 .AND. number_poisson == 1 .AND. counter .LT. 1 .AND. norm_inf /= 0.) exit_yes = .FALSE.
     IF (problem_type == 2 .AND. number_poisson == 2 .AND. counter .LT. 0 .AND. norm_inf /= 0.) exit_yes = .FALSE.
     
     IF (exit_yes .AND. counter == 0 .AND. init_yes) phi(SS1:NN1,SS2:NN2,SS3:NN3) = 0.
     IF (exit_yes) EXIT ITERATE
     !========================================================================================================
     
     
     !========================================================================================================
     !=== nächster Durchlauf =================================================================================
     !========================================================================================================
     counter = counter + 1
     !========================================================================================================
     
     
     !========================================================================================================
     rhr_prev = rhr
     IF (init_yes) THEN
        IF (counter == 1) THEN
           rhr = norm2
        ELSE
           CALL product_scalar(SS1,SS2,SS3,NN1,NN2,NN3,rr,bb,rhr)
        END IF
     ELSE
        CALL product_scalar(SS1,SS2,SS3,NN1,NN2,NN3,rr,rh,rhr)
     END IF
     !========================================================================================================
     IF (ABS(rhr) == 0.) THEN
        IF (rank == 0) WRITE(* ,'(a,E13.5)') 'rhr =', rhr
        IF (rank == 0) WRITE(10,'(a,E13.5)') 'rhr =', rhr
        
        rh(SS1:NN1,SS2:NN2,SS3:NN3) = rr(SS1:NN1,SS2:NN2,SS3:NN3) ! Neuer Referenzvektor ...
        CALL product_scalar(SS1,SS2,SS3,NN1,NN2,NN3,rr,rh,rhr)
     END IF
     !========================================================================================================
     IF (counter .GE. 2) THEN
        IF (omega == 0.) THEN
           IF (rank == 0) WRITE(* ,'(a,E13.5)') 'omega =', omega
           IF (rank == 0) WRITE(10,'(a,E13.5)') 'omega =', omega
           IF (problem_type == 2 .OR. problem_type == 4) THEN
              EXIT ITERATE
           ELSE
              CALL MPI_FINALIZE(merror)
              STOP
           END IF
        END IF
        IF (rhr_prev == 0.) THEN
           IF (rank == 0) WRITE(* ,'(a,E13.5)') 'rhr_prev =', rhr_prev
           IF (rank == 0) WRITE(10,'(a,E13.5)') 'rhr_prev =', rhr_prev
           IF (problem_type == 2 .OR. problem_type == 4) THEN
              EXIT ITERATE
           ELSE
              CALL MPI_FINALIZE(merror)
              STOP
           END IF
        END IF
        beta = (alpha/omega)*(rhr/rhr_prev)
        omega = -beta*omega
        pp(SS1:NN1,SS2:NN2,SS3:NN3) = rr(SS1:NN1,SS2:NN2,SS3:NN3) + beta*pp(SS1:NN1,SS2:NN2,SS3:NN3) + omega*Ap(SS1:NN1,SS2:NN2,SS3:NN3)
     ELSE
        IF (init_yes .AND. counter == 1) THEN
           pp(SS1:NN1,SS2:NN2,SS3:NN3) = bb(SS1:NN1,SS2:NN2,SS3:NN3)
        ELSE
           pp(SS1:NN1,SS2:NN2,SS3:NN3) = rr(SS1:NN1,SS2:NN2,SS3:NN3)
        END IF
     END IF
     !========================================================================================================
!              call MPI_Barrier(MPI_COMM_WORLD, merror)
!              t1 =   MPI_Wtime()                    
                    CALL multigridV(.TRUE.,1,pp,z1,problem_type)
!              call MPI_Barrier(MPI_COMM_WORLD, merror)
!              t2 =   MPI_Wtime()
!              print*, 'multigrid V cycle time', t2-t1          
                    CALL product_div_grad (z1,Ap)
     !========================================================================================================
     IF (init_yes) THEN
        CALL product_scalar(SS1,SS2,SS3,NN1,NN2,NN3,Ap,bb,rhAp)
     ELSE
        CALL product_scalar(SS1,SS2,SS3,NN1,NN2,NN3,Ap,rh,rhAp)
     END IF
     !========================================================================================================
     IF (ABS(rhAp) == 0.) THEN
        IF (rank == 0) WRITE(* ,'(a,E13.5)') 'rhAp =', rhAp
        IF (rank == 0) WRITE(10,'(a,E13.5)') 'rhAp =', rhAp
        IF (problem_type == 2 .OR. problem_type == 4) THEN
           EXIT ITERATE
        ELSE
           IF (rhr /= 0.) THEN
              CALL MPI_FINALIZE(merror)
              STOP
           END IF
        END IF
        alpha = 0.
     ELSE
        alpha = rhr / rhAp
     END IF
     !========================================================================================================
     IF (init_yes .AND. counter == 1) THEN
        DO k = SS3, NN3
           DO j = SS2, NN2
!pgi$ unroll = n:8
              DO i = SS1, NN1
                 rr(i,j,k) = bb(i,j,k) - alpha*Ap(i,j,k)
              END DO
           END DO
        END DO
     ELSE
        DO k = SS3, NN3
           DO j = SS2, NN2
!pgi$ unroll = n:8
              DO i = SS1, NN1
                 rr(i,j,k) = rr(i,j,k) - alpha*Ap(i,j,k)
              END DO
           END DO
        END DO
     END IF
     !========================================================================================================
   
     CALL multigridV(.TRUE.,1,rr,z2,problem_type)
      CALL product_div_grad (z2,Ar)
     
     !========================================================================================================
     rAr  = 0.
     ArAr = 0.
     DO k = SS3, NN3
        DO j = SS2, NN2
!pgi$ unroll = n:8
           DO i = SS1, NN1
              rAr  = rAr  + rr(i,j,k)*Ar(i,j,k)
              ArAr = ArAr + Ar(i,j,k)*Ar(i,j,k)
           END DO
        END DO
     END DO
     
     CALL MPI_ALLREDUCE((/rAr,ArAr/),scalar_global,2,MPI_REAL8,MPI_SUM,COMM_CART,merror)
     rAr  = scalar_global(1)
     ArAr = scalar_global(2)
     
     ! ACHTUNG!!! Zu teuer im Vergleich zur obigen Variante:
     !CALL product_scalar(SS1,SS2,SS3,NN1,NN2,NN3,Ar,rr,rAr)
     !CALL product_scalar(SS1,SS2,SS3,NN1,NN2,NN3,Ar,Ar,ArAr)
     !========================================================================================================
     IF (ABS(rAr) == 0.) THEN
        IF (rank == 0) WRITE(* ,'(a,E13.5)') 'rAr =', rAr
        IF (rank == 0) WRITE(10,'(a,E13.5)') 'rAr =', rAr
     END IF
     IF (ABS(ArAr) == 0.) THEN
        IF (rank == 0) WRITE(* ,'(a,E13.5)') 'ArAr =', ArAr
        IF (rank == 0) WRITE(10,'(a,E13.5)') 'ArAr =', ArAr
        IF (problem_type == 2 .OR. problem_type == 4) THEN
           EXIT ITERATE
        ELSE
           IF (rAr /= 0.) THEN
              CALL MPI_FINALIZE(merror)
              STOP
           END IF
        END IF
        omega = 0.
     ELSE
        omega = rAr / ArAr
     END IF
     !========================================================================================================
     norm_inf = 0.
  !   norm_inf2 = 0.
    ! GPUres(:,:,:)=0.

     !--------------------------------------------------------------------------------------------------------
              DO k = SS3, NN3
                 DO j = SS2, NN2
!pgi$ unroll = n:8
                    DO i = SS1, NN1
                       phi(i,j,k) = phi(i,j,k) + omega*z2(i,j,k) + alpha*z1(i,j,k)
                       rr (i,j,k) = rr (i,j,k) - omega*Ar(i,j,k)
                       norm_inf = MAX(ABS(rr(i,j,k)),norm_inf)
                    END DO
                 END DO
              END DO
       CALL MPI_ALLREDUCE(norm_inf,norm_inf_global,1,MPI_REAL8,MPI_MAX,COMM_CART,merror) ! MPI_REDUCE bringt nichts, weil exit_yes dann mit MPI_BCAST verteilt werden müsste ...
       norm_inf = norm_inf_global

!   CALL exchange(1,1,phi(:,:,:))
!    CALL exchange(2,2,phi(:,:,:))
!    CALL exchange(3,3,phi(:,:,:))
!
!
!

! CALL exchange_relax(1,0,0,0,0,.FALSE.,phi)
!
!       norm_lap = 0.
!       DO k = SS3, NN3
!          DO j = SS2, NN2
!             DO i = SS1, NN1
!                GPUres(i,j,k) =0 
!                GPUres(i,j,k) = bb(i,j,k)                                     &
!                            &      - cdg1(-1,i,1)*phi(i-1,j,k) - cdg1(1,i,1)*phi(i+1,j,k)     &
!                            &      - cdg2(-1,j,1)*phi(i,j-1,k) - cdg2(1,j,1)*phi(i,j+1,k)     &
!                            &      - cdg3(-1,k,1)*phi(i,j,k-1) - cdg3(1,k,1)*phi(i,j,k+1)    &
!                            &      - (cdg1(0,i,1) + cdg2(0,j,1) + cdg3(0,k,1))*phi(i,j,k)
!                norm_lap = MAX(ABS(GPUres(i,j,k)),norm_lap)
!               ! if (norm_lap .ge. 1.) print*, i,j,k, norm_lap 
!            END DO
!          END DO
!       END DO
!      ! if ((counter ==3) .and. (timestep ==2)) print*, GPUres
!       CALL MPI_ALLREDUCE(norm_lap,norm_lap_global,1,MPI_REAL8,MPI_MAX,COMM_CART,merror) ! MPI_REDUCE bringt nichts, weil exit_yes dann mit MPI_BCAST verteilt werden müsste ...
!       norm_lap = norm_lap_global
!!
!!     !========================================================================================================
!        if (rank ==0 ) print*, counter, norm_lap, norm_inf



!        DO k = SS3, NN3
!           DO j = SS2, NN2
!              DO i = SS1, NN1
!                 GPUres(i,j,k) = bb(i,j,k)                                     &
!                             &      - cdg1(-1,i,1)*phi(i-1,j,k) - cdg1(1,i,1)*phi(i+1,j,k)     &
!                             &      - cdg2(-1,j,1)*phi(i,j-1,k) - cdg2(1,j,1)*phi(i,j+1,k)     &
!                             &      - cdg3(-1,k,1)*phi(i,j,k-1) - cdg3(1,k,1)*phi(i,j,k+1)    &
!                             &      - (cdg1(0,i,1) + cdg2(0,j,1) + cdg3(0,k,1))*phi(i,j,k)
!                 norm_inf2 = MAX(ABS(GPUres(i,j,k)),norm_inf)
!              END DO
!           END DO
!        END DO





    !--------------------------------------------------------------------------------------------------------
 !    CALL MPI_ALLREDUCE(norm_inf,norm_inf_global,1,MPI_REAL8,MPI_MAX,COMM_CART,merror) ! MPI_REDUCE bringt nichts, weil exit_yes dann mit MPI_BCAST verteilt werden müsste ...
 !    norm_inf = norm_inf_global
!      CALL MPI_ALLREDUCE(norm_inf2,norm_inf_global2,1,MPI_REAL8,MPI_MAX,COMM_CART,merror) ! MPI_REDUCE bringt nichts, weil exit_yes dann mit MPI_BCAST verteilt werden müsste ...
!     norm_inf2 = norm_inf_global2
     !========================================================================================================
 !    if (rank==1) print*, counter, norm_inf, norm_inf2 
     
     !========================================================================================================
     !=== Konvergenzstatistik ================================================================================
     !========================================================================================================
     ratioP(substep,number_poisson) = ratioP(substep,number_poisson) + LOG10(norm_inf/norm_inf_prev)
     
     countP(substep,number_poisson) = countP(substep,number_poisson) + 1
     
     norm_inf_prev = norm_inf
     !========================================================================================================
     
  END DO ITERATE
  
  
  END SUBROUTINE BiCGstab
  
  
  
  
  
  
  
  
  
 
  
 
  !> computed the element-wise product of two fields (arrays) 
  SUBROUTINE product_scalar(SS1,SS2,SS3,NN1,NN2,NN3,phi1,phi2,scalar)
  
  IMPLICIT NONE
  
  INTEGER, INTENT(in)    ::  SS1 !< lower bound dimension one
  INTEGER, INTENT(in)    ::  SS2 !< lower bound dimension two
  INTEGER, INTENT(in)    ::  SS3 !< lower bound dimension three
  
  INTEGER, INTENT(in)    ::  NN1 !< upper bound dimension one
  INTEGER, INTENT(in)    ::  NN2 !< upper bound dimension two
  INTEGER, INTENT(in)    ::  NN3 !< upper bound dimension three
  
  REAL(8)   , INTENT(in)    ::  phi1(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)) !< array 1
  REAL(8)   , INTENT(in)    ::  phi2(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)) !< array 2
  
  REAL(8)   , INTENT(out)   ::  scalar !< element-wise product
  REAL(8)                   ::  scalar_global
  INTEGER                ::  i, j, k
  
  
  scalar = 0.
  
  DO k = SS3, NN3
     DO j = SS2, NN2
!pgi$ unroll = n:8
        DO i = SS1, NN1
           scalar = scalar + phi1(i,j,k)*phi2(i,j,k)
        END DO
     END DO
  END DO
  
  CALL MPI_ALLREDUCE(scalar,scalar_global,1,MPI_REAL8,MPI_SUM,COMM_CART,merror)
  scalar = scalar_global
  
  
  END SUBROUTINE product_scalar
  
  
  !> subroutine that compute infinity and/or two norms.  
  SUBROUTINE get_norms(SS1,SS2,SS3,NN1,NN2,NN3,phi,problem_type,inf_yes,two_yes,normInf,normTwo)
  
  IMPLICIT NONE
  
  INTEGER, INTENT(in)    ::  SS1 !< lower bound first dimension
  INTEGER, INTENT(in)    ::  SS2 !< lower bound second dimension
  INTEGER, INTENT(in)    ::  SS3 !< lower bound third dimension
  
  INTEGER, INTENT(in)    ::  NN1 !< upper bound first dimension
  INTEGER, INTENT(in)    ::  NN2 !< upper bound second dimension
  INTEGER, INTENT(in)    ::  NN3 !< upper bound third dimension
  
  REAL(8)   , INTENT(in)    ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)) !< field to compute the norm of
  
  INTEGER, INTENT(in)    ::  problem_type !< type of problem
  LOGICAL, INTENT(in)    ::  inf_yes      !< whether to compute infinity norm
  LOGICAL, INTENT(in)    ::  two_yes      !< whether to compute two-norm
  
  REAL(8)   , OPTIONAL, INTENT(out) ::  normInf
  REAL(8)   , OPTIONAL, INTENT(out) ::  normTwo
  
  REAL(8)                   ::  normInf_global, normTwo_global
  INTEGER                ::  i, j, k
  
  
  IF (inf_yes .AND. two_yes) THEN
     
     normInf = 0.
     normTwo = 0.
     
     IF (problem_type == 2 .AND. weighting_yes) THEN
        
        DO k = SS3, NN3
           DO j = SS2, NN2
!pgi$ unroll = n:8
              DO i = SS1, NN1
                 normInf = MAX(ABS(phi(i,j,k)*weight(i,j,k)),normInf)
                 normTwo = normTwo + phi(i,j,k)**2
              END DO
           END DO
        END DO
        
     ELSE
        
        DO k = SS3, NN3
           DO j = SS2, NN2
!pgi$ unroll = n:8
              DO i = SS1, NN1
                 normInf = MAX(ABS(phi(i,j,k)),normInf)
                 normTwo = normTwo + phi(i,j,k)**2
              END DO
           END DO
        END DO
        
     END IF
     
     ! Lassen sich wegen MPI_SUM / MPI_MAX nicht zusammenlegen:
     CALL MPI_ALLREDUCE(normInf,normInf_global,1,MPI_REAL8,MPI_MAX,COMM_CART,merror)
     CALL MPI_ALLREDUCE(normTwo,normTwo_global,1,MPI_REAL8,MPI_SUM,COMM_CART,merror)
     normInf = normInf_global
     normTwo = normTwo_global
     
  ELSE IF (inf_yes) THEN
     
     normInf = 0.
     
     IF (problem_type == 2 .AND. weighting_yes) THEN
        
        DO k = SS3, NN3
           DO j = SS2, NN2
!pgi$ unroll = n:8
              DO i = SS1, NN1
                 normInf = MAX(ABS(phi(i,j,k)*weight(i,j,k)),normInf)
              END DO
           END DO
        END DO
        
     ELSE
        
        DO k = SS3, NN3
           DO j = SS2, NN2
!pgi$ unroll = n:8
              DO i = SS1, NN1
                 normInf = MAX(ABS(phi(i,j,k)),normInf)
              END DO
           END DO
        END DO
        
     END IF
     
     CALL MPI_ALLREDUCE(normInf,normInf_global,1,MPI_REAL8,MPI_MAX,COMM_CART,merror) ! MPI_REDUCE bringt nichts, weil exit_yes dann mit MPI_BCAST verteilt werden müsste ...
     normInf = normInf_global
     
  ELSE IF (two_yes) THEN
     
     normTwo = 0.
     
     DO k = SS3, NN3
        DO j = SS2, NN2
!pgi$ unroll = n:8
           DO i = SS1, NN1
              normTwo = normTwo + phi(i,j,k)**2
           END DO
        END DO
     END DO
     
     CALL MPI_ALLREDUCE(normTwo,normTwo_global,1,MPI_REAL8,MPI_SUM,COMM_CART,merror)
     normTwo = normTwo_global
     
  END IF
  
  
  END SUBROUTINE get_norms
  
  
  
  
  
  
  
  
   
  
  !> checks whether iteration criterion was met or whether maximum numbers of iterations have been reached.
  SUBROUTINE status_iteration(eps,norm,counter,N_restarts,exit_yes,quiet_yes1,quiet_yes2)
  
  IMPLICIT NONE
  
  REAL(8)   , INTENT(in)    ::  eps       !< convergence criterion
  REAL(8)   , INTENT(in)    ::  norm      !< (two, infinity) norm
  INTEGER, INTENT(in)    ::  counter   !< iteration counter (current)
  INTEGER, INTENT(in)    ::  N_restarts!< max iteration count
  LOGICAL, INTENT(out)   ::  exit_yes  !< whether to stop the iteration
  LOGICAL, INTENT(in)    ::  quiet_yes1!< whether to suppress output during iteration
  LOGICAL, INTENT(in)    ::  quiet_yes2!< whether to suppress output for termination
  
  
  exit_yes = .FALSE.
  
  IF (norm .LT. eps) THEN
     IF (rank == 0 .AND. write_stout_yes .AND. .NOT. quiet_yes2) WRITE(* ,'(a,i5,a,E24.17,a)') '  Iteration',counter,'; ||res|| =',norm,'  (Termination criterion satisfied)'
     IF (rank == 0 .AND. log_iteration_yes                     ) WRITE(10,'(a,i5,a,E24.17,a)') '  Iteration',counter,'; ||res|| =',norm,'  (Termination criterion satisfied)'
     exit_yes = .TRUE.
  END IF
  
  IF ((.NOT. exit_yes) .AND. counter == N_restarts) THEN
     IF (rank == 0 .AND. write_stout_yes .AND. .NOT. quiet_yes2) WRITE(* ,'(a,i5,a,E24.17,a)') '  Iteration',counter,'; ||res|| =',norm,'  WARNING! Too many iterations!'
     IF (rank == 0 .AND. log_iteration_yes                     ) WRITE(10,'(a,i5,a,E24.17,a)') '  Iteration',counter,'; ||res|| =',norm,'  WARNING! Too many iterations!'
     exit_yes = .TRUE.
  END IF
  
  IF (.NOT. exit_yes) THEN
     IF (rank == 0 .AND. write_stout_yes .AND. .NOT. quiet_yes1) WRITE(* ,'(a,i5,a,E24.17  )') '  Iteration',counter,'; ||res|| =',norm
     IF (rank == 0 .AND. log_iteration_yes                     ) WRITE(10,'(a,i5,a,E24.17  )') '  Iteration',counter,'; ||res|| =',norm
  END IF
  
  
  END SUBROUTINE status_iteration
  
  
  
  
  
  
 
  
  
END MODULE mod_solvers
