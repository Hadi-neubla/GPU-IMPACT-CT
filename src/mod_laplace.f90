!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!* modified by Barna Becsek, ARTORG CVE, University of Bern (barna.becsek@artorg.unibe.ch)                   *
!* October 2014                                                                                              *
!* GPU version by Hadi Zolfaghari , ARTORG CVE, then DAMTP, Cambridge University (hz382@cam.ac.uk)           *
!* Oct 2015 - Sep 2023                                                                                       *
!*************************************************************************************************************


!> module containing subroutines for computing the Laplacian operator. It uses modules mod_dims, mod-vars,
!! mod_diff and mod_exchange.
MODULE mod_laplace


  USE mod_dims
  USE mod_vars
  USE mod_vars_GPU
  USE mod_diff
  USE mod_exchange
  
  PRIVATE
  
  PUBLIC product_div_grad, product_div_grad_transp
  PUBLIC product_div_grad_relax, relaxation_div_grad, relaxation_div_grad_inv
  PUBLIC handle_corner_Lap  
  
  
  CONTAINS
  
!pgi$g unroll = n:8
!!pgi$r unroll = n:8
!!pgi$l unroll = n:8
  
  
  !> subroutine that computes the Laplacian operator 
  SUBROUTINE product_div_grad(phi,Lap)
  
  IMPLICIT NONE
 
 
  REAL(8)   , target, INTENT(inout) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)) !> ?
  
  REAL(8)   , target,  INTENT(out  ) ::  Lap(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)) !> Laplacian

  REAL(8)                   ::  save_phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)) !> Laplacian

  REAL(8)                   ::  test_dig(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)) !> Laplacian

  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  
  INTEGER                ::  m
  
  integer                ::  n_match

  real(8)                   ::  t1, t2, t3, t4

  include 'mpif.h'
  !----------------------------------------------------------------------------------------------------------!
  ! Achtung: - Ein direktes Aufstellen des Laplace-Operators funktioniert i.A. nur sehr schlecht als Vor-    !
  !            konditionierer.                                                                               !
  !          - Das Feld "grad" kann nur eliminiert werden, wenn der div-grad-Stencil komplett ausmultipli-   !
  !            ziert wird (aktuell zu aufwendig).                                                            !
  !----------------------------------------------------------------------------------------------------------!
  
  
  DO m = 1, dimens
      
     save_phi = phi

!#ifdef cuIMPACT
  if (GPU_accelerated == 9) then

!for performance measurement uncomment two lines below
!     call MPI_Barrier(MPI_COMM_WORLD, merror)
!     t1 =   MPI_Wtime() 
     call device_gradient(m,phi,dig)
!for performance measuremnt uncomment five lines below
!     call MPI_Barrier(MPI_COMM_WORLD, merror)
!     t2 =   MPI_Wtime()
!     call gradient(m,phi,dig)

!     call MPI_Barrier(MPI_COMM_WORLD, merror)
!     t3 =   MPI_Wtime()

     if (GPU_verbose == 1)   print*, '---------- ---------- -----------on GPU gradient'
     
!     print*, 'GPU time: ',  t2- t1, 'CPU time: ', t3-t2 


!     print*, 'max value: ', maxval(dig), 'from GPU'

  else
!#else
     CALL gradient(m, phi, dig)
!     print*, 'on CPU'

     ! print*, t4 -t3
!     print*, 'CPU time: ', t4 - t3,' GPU time: ', t2- t1
!     print*, 'max value: ', maxval(dig), ' from CPU'
  endif 

!#endif

!     n_match =0
!
!     do i=1, N11
!         do j=1, N21
!            do k=1, N31
! 
!               if (abs(test_dig(i,j,k) - dig(i,j,k)) .le. 0.00000000000001 ) then
!                  n_match = n_match +1
!                  if (n_match .eq. 356) then
!                     print*, test_dig(i,j,k), 'GPU value'
!                     print*, dig(i,j,k), 'out of gpu pre value'
!                  endif
!                endif
!            enddo
!         enddo
!      enddo
!      print*, n_match
 
 
 

     if (GPU_accelerated == 9) then

     if (GPU_verbose == 1)     print*, '---------- ---------- -----------on GPU divergence'
!-----uncomment for unit testing or speedup measuremeant
!     save_phi = Lap
!     test_dig = dig

      
     
    ! call cpu_time(t1)


!----uncomment for performance measurement
!    call MPI_Barrier(MPI_COMM_WORLD, merror)
!    t1 =   MPI_Wtime()

    !call  divergence    (m,test_dig, save_phi)
!    call MPI_Barrier(MPI_COMM_WORLD, merror)
!    t2 =   MPI_Wtime()
    call device_divergence(m, dig, Lap)     
    

!    call MPI_Barrier(MPI_COMM_WORLD, merror)
!    t3 =   MPI_Wtime()

!     print*, 'GPU time: ',  t3- t2, 'CPU time: ', t2-t1 

!      call divergence(m, dig, Lap)


!      print*, 'these are the limits', b1L, N1 + b1U, b2L, N2 + b2U, b3L, N3 + b3U 

      else ! if GPU_accelerated == 0 
 
!-----Uncomment for unit testing
!    save_phi = Lap
!    test_dig = dig

    call divergence(m,dig,Lap)

!---- Unit testing for the divergence kernel

    
!    call device_divergence(m,test_dig,save_phi)
!    n_match =0
!
!    do i = b1L, N1+b1U
!        do j = b2L, N2+b2U
!           do k = b3L, N3+b3U
!              if (abs(Lap(i,j,k) - save_phi(i,j,k)) .le. 0.00000000000001 ) then
!                 n_match = n_match +1
!             else
!               print*, 'here is not matching', Lap(i,j,k), save_phi(i,j,k), i,j,k
!               endif
!           enddo
!        enddo
!     enddo
!   
!     print*, n_match, m, (N1+b1U-b1L + 1) * (N2+b2U-b2L + 1) * (N3+b3U-b3L + 1)

!--- end of unit testing

      endif
     
  END DO
  
  
  
  END SUBROUTINE product_div_grad
  
  
  
  
  
  
  
  
  
  
 
  
  
  
  
  
  
  
  SUBROUTINE product_div_grad_relax(g,phi,Lap) ! TEST!!! aufraeumen und Variablen substituieren ...
  
  IMPLICIT NONE
  
  INTEGER, INTENT(in   ) ::  g
  
  REAL(8)   , INTENT(inout) ::  phi(b1L:(NN(1,g)+b1U),b2L:(NN(2,g)+b2U),b3L:(NN(3,g)+b3U))
  REAL(8)   , INTENT(out  ) ::  Lap(b1L:(NN(1,g)+b1U),b2L:(NN(2,g)+b2U),b3L:(NN(3,g)+b3U))

  ! for unit testing
!*  REAL(8)    ::  phi_test(b1L:(NN(1,g)+b1U),b2L:(NN(2,g)+b2U),b3L:(NN(3,g)+b3U))
!*  REAL(8)    ::  Lap_test(b1L:(NN(1,g)+b1U),b2L:(NN(2,g)+b2U),b3L:(NN(3,g)+b3U)) 
  
  INTEGER                ::  i, N1, S1R, N1R, S11R, N11R
  INTEGER                ::  j, N2, S2R, N2R, S22R, N22R
  INTEGER                ::  k, N3, S3R, N3R, S33R, N33R
  
  INTEGER                ::  BC_1L, BC_1U, BC_2L, BC_2U, BC_3L, BC_3U

  
  !*****************************************************************************************
  N1 = NN(1,g)
  N2 = NN(2,g)
  N3 = NN(3,g)
  
  BC_1L = BC(1,1,g)
  BC_1U = BC(2,1,g)
  BC_2L = BC(1,2,g)
  BC_2U = BC(2,2,g)
  BC_3L = BC(1,3,g)
  BC_3U = BC(2,3,g)
 

 
  S1R  = SNB(1,1,g)
  S2R  = SNB(1,2,g)
  S3R  = SNB(1,3,g)
  
  N1R  = SNB(2,1,g)
  N2R  = SNB(2,2,g)
  N3R  = SNB(2,3,g)
  
  S11R = SNF(1,1,g)
  S22R = SNF(1,2,g)
  S33R = SNF(1,3,g)
  
  N11R = SNF(2,1,g)
  N22R = SNF(2,2,g)
  N33R = SNF(2,3,g)
  !*****************************************************************************************
  
  
  !-----------------------------------------------------------------------------------------------------!
  !              der Interpolation verwendet (ebenfalls zur Vereinfachung der Programmierung) und ist   !
  !              somit ohnehin vorhanden. GeschwindigkeitsmÃ¤ssig ist vermutlich auch nicht mehr viel zu !
  !              holen.                                                                                 !        
  !-----------------------------------------------------------------------------------------------------!
  
  
  CALL exchange_relax(g,0,0,0,0,.TRUE.,phi)
  
  
  !===========================================================================================================
  IF (dimens == 3) THEN

  if ((GPU_accelerated == 0)  .and. ((N33R .ge. 32) .and. (N22R .ge. 32) .and. (N11R .ge. 32)))  then
     call device_product_dg_relax(g,phi,Lap)
  else 

!*   if  ((N33R .ge. 32) .and. (N22R .ge. 32) .and. (N11R .ge. 32)) then
!*        Lap_test=Lap
!*        phi_test =phi
!*        call device_product_dg_relax(g,phi_test,Lap_test)
!*   endif

     DO k = S33R, N33R
        DO j = S22R, N22R
!pgi$ unroll = n:8
           DO i = S11R, N11R
              Lap(i,j,k) =  cdg1(-1,i,g)*phi(i-1,j,k) + cdg1(1,i,g)*phi(i+1,j,k)   &
                      &  +  cdg2(-1,j,g)*phi(i,j-1,k) + cdg2(1,j,g)*phi(i,j+1,k)   &
                      &  +  cdg3(-1,k,g)*phi(i,j,k-1) + cdg3(1,k,g)*phi(i,j,k+1)   &
                      &  + (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g))*phi(i,j,k)
!*           if  ((N33R .ge. 32) .and. (N22R .ge. 32) .and. (N11R .ge. 32)) then
!*                      if (abs(Lap_test(i,j,k) - Lap(i,j,k)) .ge. 0.000000000001) print*, 'mismatch', Lap_test(i,j,k), Lap(i,j,k)
!*           endif
           END DO
        END DO
     END DO
  endif
  !===========================================================================================================
    
  END IF
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== Randbedingungen =======================================================================================
  !===========================================================================================================
  
 
  !===========================================================================================================
  
 !===========================================================================================================
  
 
  !===========================================================================================================
  
  
  
  
  END SUBROUTINE product_div_grad_relax
  
  
  
  
  
  
  
  
  
  ! Habe RB-Linienrelaxation herausgeschmissen, weil die Gewinne in der Geschwindigkeit zu gering waren.
  ! Gleichzeitig war der Code viel zu lang und unuebersichtlich.
  SUBROUTINE relaxation_div_grad(init_yes,n_relax,g,bb,rel) ! TEST!!! reine 2D-Variante fehlt noch ... ! TEST!!! aufraeumen und Variablen substituieren ...
  
  IMPLICIT NONE
  
  !*************************************************************************************************
  INTEGER                ::  BC_1L, BC_1U, BC_2L, BC_2U, BC_3L, BC_3U
  !*************************************************************************************************
  
  LOGICAL, INTENT(in   ) ::  init_yes
  INTEGER, INTENT(in   ) ::  n_relax
  
  INTEGER, INTENT(in   ) ::  g
  
  REAL(8)   , INTENT(in   ) ::  bb  (b1L:(NN(1,g)+b1U),b2L:(NN(2,g)+b2U),b3L:(NN(3,g)+b3U))
  REAL(8)   , INTENT(inout) ::  rel (b1L:(NN(1,g)+b1U),b2L:(NN(2,g)+b2U),b3L:(NN(3,g)+b3U))
  REAL(8)                   ::  comp(b1L:(NN(1,g)+b1U),b2L:(NN(2,g)+b2U),b3L:(NN(3,g)+b3U)) ! FELD!!!
!*  REAL(8)                   ::  comp_gpu(b1L:(NN(1,g)+b1U),b2L:(NN(2,g)+b2U),b3L:(NN(3,g)+b3U)) ! FELD!!!
!*  REAL(8)                   ::  bb_copy(b1L:(NN(1,g)+b1U),b2L:(NN(2,g)+b2U),b3L:(NN(3,g)+b3U)) ! FELD!!!
!*  REAL(8)                   ::  rel_copy(b1L:(NN(1,g)+b1U),b2L:(NN(2,g)+b2U),b3L:(NN(3,g)+b3U)) ! FELD!!!
  
  INTEGER                ::  i, N1, S1R, N1R, S11R, N11R
  INTEGER                ::  j, N2, S2R, N2R, S22R, N22R
  INTEGER                ::  k, N3, S3R, N3R, S33R, N33R
  
  INTEGER                ::  r, ss
  REAL(8)                   ::  mult
  
  ! TEST!!! herausziehen?
  ! RBGS_mode =< 0 : naturally ordered Gauss-Seidel (slow on CPU, cannot be vectorized)
  ! RBGS_mode == 1 : 3D-Red-Black Gauss-Seidel (faster on CPU due to vectorization, parallelization-independent)
  ! RBGS_mode >= 2 : 2D-Red-Black Gauss-Seidel (normally the fastest variant on CPU due to better cache utilization)
  INTEGER, PARAMETER     ::  RBGS_mode = 2
  LOGICAL, PARAMETER     ::  SOR_yes = .FALSE.
  !REAL(8)   , PARAMETER     ::  omega   = 0.8! 1.2 ! 1.27 !1.27  ! omega = 0.9
  REAL(8), target                   ::  omega


  integer                :: jacobi_iterant
  integer,parameter      :: upper_limit_jacobi = 1

  
  integer                ::  nmatch  

  real(8)                   ::  t1, t2, t3, t4

  integer                ::  gpu_limit

  include 'mpif.h'
  !----------------------------------------------------------------------------------------------------------!
  ! Optimization annotations:                                                                                !
  !    - line relaxation runs on Rosa (Cray XT5, Istanbul hex-core processors) in single core mode with      !
  !      about 5.6-5.8% peak performance in direction 1, 4.2-4.6% in direction 2, 3.9-4.0% in direction 3.   !
  !    - line relaxation in direction 3 is about 50% more expensive than in direction 1; direction 2 is      !
  !      somewhere in between.                                                                               !
  !    - line relaxation in direction 1 is typically 30% of the total execution time(!).                     !
  !    - RB-GS converges slightly worse than standard GS.                                                    !
  !    - RB-GS allows full vectorization (execpt some operations in direction 1), provided that the loops    !
  !      are reordered to let direction 1 be the innermost (speedup was not as good as hoped (<5%) for       !
  !      direction 1 and 2, direction 3 was not tested).                                                     !
  !    - unrolling the non-vectorizable loops 8 times gives roughly the best results.                        !
  !    - IF statements are automatically pulled out of the loops by the compiler.                            !
  !----------------------------------------------------------------------------------------------------------!
  
  ! Allgemein, pgf90: "Invariant if transformation" <==> IF-Statements werden aus dem Loop herauskopiert, muss man nicht von Hand machen!
  
  
  !*****************************************************************************************

!  print*, b1L, NN(1,g) +b1U, b2L, NN(2,g)+ b2U, b3L, NN(3,g) + b3U, 'bb shape'
 


  N1 = NN(1,g)
  N2 = NN(2,g)
  N3 = NN(3,g)
  
  BC_1L = BC(1,1,g)
  BC_1U = BC(2,1,g)
  BC_2L = BC(1,2,g)
  BC_2U = BC(2,2,g)
  BC_3L = BC(1,3,g)
  BC_3U = BC(2,3,g)


  
  S1R  = SNB(1,1,g)
  S2R  = SNB(1,2,g)
  S3R  = SNB(1,3,g)
  
  N1R  = SNB(2,1,g)
  N2R  = SNB(2,2,g)
  N3R  = SNB(2,3,g)
  
  S11R = SNF(1,1,g)
  S22R = SNF(1,2,g)
  S33R = SNF(1,3,g)
  
  N11R = SNF(2,1,g)
  N22R = SNF(2,2,g)
  N33R = SNF(2,3,g)
 

!print*, S11R, S22R, S33R, 'S11R, S22R, S33R'
!print*, N11R, N22R, N33R, 'N11R, N22R, N33R'
 
!print*, SNF(2,1, :)


 !*****************************************************************************************
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Null-Setzen am Rand nicht notwendig, da Startindizes entsprechend gewÃ¤hlt sind!           !
  !              - Bei der Initialisierung mÃ¼ssen die Intervallgrenzen SiR:NiR anstelle von SiiR:NiiR        !
  !                gewÃ¤hlt werden, da beim  Aufbau der RHS ("vec") innerhalb der Linienrelaxation IMMER auch !
  !                die Randbereiche SiR und NiR aufgebaut aber ggf. mit einer Korrektur wieder Ã¼berschrieben !
  !                werden. Ansonsten wÃ¤re eine Initialisierung allein im Feldbereich SiiR:NiiR ausreichend,  !
  !                kann aber z.B. zu Floating Point Exceptions fÃ¼hren (die fÃ¼r die eigentliche Rechnung      !
  !                allerdings irrelevant wÃ¤ren)!                                                             !
  !              - obere StirnflÃ¤chen werden bei der Initialisierung ebenfalls berÃ¼cksichtigt, da ggf.       !
  !                verschiedene Richtungen nacheinander bearbeitet werden und dies beim Aufbauen der rechten !
  !                Seite sonst berÃ¼cksichtigt werden mÃ¼sste. ==> geringerer Mehraufwand beim Rechnen,        !
  !                weniger Programmierdetails.                                                               !
  !              - LU-Zerlegung (d.h. band, mult) kann gespeichert werden, wenn mindestens eine Richtung     !
  !                Ã¤quidistant ist. Der LÃ¶sungsaufwand wÃ¼rde sich etwa halbieren!!                           !
  !              - "r == 1 .AND. init_yes" sollte idealerweise aus den Schleifen herausgezogen werden, was   !
  !                hier aber aus GrÃ¼nden der Ãbersicht bisher nicht ausgefÃ¼hrt wurde.                        !
  !----------------------------------------------------------------------------------------------------------!
  
  
  omega = 0.98
  
  
  IF (Jacobi_yes) THEN
  DO r = 1, n_relax
     
     IF (r == 1 .AND. init_yes) THEN
        !=====================================================================================================
       !=====================================================================================================
        DO k = S33R, N33R
           DO j = S22R, N22R
!pgi$ unroll = n:8
              DO i = S11R, N11R
                 rel(i,j,k) = omega*bb(i,j,k) / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g))
              END DO
           END DO
        END DO
        !=====================================================================================================
       !=====================================================================================================
        
     ELSE
        
        CALL exchange_relax(g,0,0,0,0,.TRUE.,rel)
        
        !=====================================================================================================
       !====================================================================================================
  
  if ( ((N11R .lt. 64) .and. (N22R .lt. 64) .and. (N33R .lt. 64) ) .or. (GPU_accelerated == 0) ) then

!*     comp_gpu =comp
!*     bb_copy = bb
!*     rel_copy = rel

 !    call MPI_Barrier(MPI_COMM_WORLD, merror)
 !    t1 =   MPI_Wtime()

!*   if ( ((N11R .ge. 32) .and. (N22R .ge. 32) .and. (N33R .ge. 32) ) ) then
!*       call device_jacobi_mg_relaxation(g,omega, rel_copy, bb_copy, comp_gpu)
!*  endif
        DO k = S33R, N33R
           DO j = S22R, N22R
!pgi$ unroll = n:8
              DO i = S11R, N11R
                 comp(i,j,k) = omega*(bb(i,j,k)                                                &
                             &      - cdg1(-1,i,g)*rel(i-1,j,k) - cdg1(1,i,g)*rel(i+1,j,k)     &
                             &      - cdg2(-1,j,g)*rel(i,j-1,k) - cdg2(1,j,g)*rel(i,j+1,k)     &
                             &      - cdg3(-1,k,g)*rel(i,j,k-1) - cdg3(1,k,g)*rel(i,j,k+1))    &
                             &      / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g)) + (1.-omega)*rel(i,j,k)
!*         if ( ((N11R .ge. 32) .and. (N22R .ge. 32) .and. (N33R .ge. 32) ) ) then
!*                  if (abs(comp_gpu(i,j,k) - comp(i,j,k)) .ge. 0.000000000000001) print*, 'mismatch'
!*        endif
              END DO
           END DO
        END DO
!     call MPI_Barrier(MPI_COMM_WORLD, merror)
!     t2 =   MPI_Wtime()


 ! if ( ((N11R .ge. 32) .and. (N22R .ge. 32) .and. (N33R .ge. 32) ) ) then

    else
          if (GPU_verbose == 1)    print*, '---------- ---------- -----------On GPU Jacobi is Relaxing'
             call device_jacobi_mg_relaxation(g,omega, rel, bb, comp)
!     call MPI_Barrier(MPI_COMM_WORLD, merror)
!     t3 =   MPI_Wtime()

!     print*, N11R, N22R, N33R, 'GPU time =', t3-t2, 'CPU time =', t2-t1
! endif 
    endif



!     print*, 'For g = ', g, ' and size: ', N11R, N22R, N33R, ' : CPU time = ',t2-t1, ' and GPU time = ',t3-t2  


!!  endif
!   
!
!         nmatch = 0
!         DO k = S33R, N33R
!            DO j = S22R, N22R
!!pgi$ unroll = n:8
!               DO i = S11R, N11R
!              !    print*, comp_gpu(i,j,k) 
!                  if (abs(comp_gpu(i,j,k) - comp(i,j,k)) .le. 0.000000000000001) then
!                                              nmatch = nmatch + 1
!
!                  else  
!
!               !       print*, comp(i,j,k), comp_gpu(i,j,k), 'comp and comp_gpu'
!               !       print*, i, j, k
!                       if (nmatch ==300)  then
!                       !      print*, comp_gpu(i,j,k) , 'gpu result'
!                       !      print*, comp(i,j,k) , 'cpu result'
!
!
!                       endif
!
!
!
!                  endif
!               END DO
!            END DO
!         END DO
!
!
!
!
!
!         print*, nmatch,  'grid: ', g
!         print*, 'the limits of computationas are: ', N33R, N22R, N11R
!         print*, maxval (comp) , maxval (comp_gpu), 'max of cpu and gpu'
!         print*, minval (comp) , minval (comp_gpu), 'min of cpu and gpu'
!         print*, sum(comp) , sum(comp_gpu), 'sum of cpu and gpu'
!
!    endif
!     print*, 'I am here anyway'


!    endif 



        !=====================================================================================================
       !=====================================================================================================
        
        rel(1:N1,1:N2,1:N3) = comp(1:N1,1:N2,1:N3)
        
     END IF
     
  END DO
  
  
  !***********************************************************************************************************
  !***********************************************************************************************************
  !***********************************************************************************************************
 !***********************************************************************************************************
  !***********************************************************************************************************
  !***********************************************************************************************************
  END IF
  
  
  
  
  END SUBROUTINE relaxation_div_grad
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE relaxation_div_grad_inv(init_yes,n_relax,g,bb,rel) ! TEST!!! 2D-Variante fehlt noch ... ! TEST!!! aufraeumen und Variablen substituieren ...
  
  IMPLICIT NONE
  
  !*************************************************************************************************
  INTEGER                ::  BC_1L, BC_1U, BC_2L, BC_2U, BC_3L, BC_3U
  !*************************************************************************************************
  
  LOGICAL, INTENT(in   ) ::  init_yes
  INTEGER, INTENT(in   ) ::  n_relax
  
  INTEGER, INTENT(in   ) ::  g
  
  REAL(8)   , INTENT(in   ) ::  bb  (b1L:(NN(1,g)+b1U),b2L:(NN(2,g)+b2U),b3L:(NN(3,g)+b3U))
  REAL(8)   , INTENT(inout) ::  rel (b1L:(NN(1,g)+b1U),b2L:(NN(2,g)+b2U),b3L:(NN(3,g)+b3U))
  REAL(8)                   ::  comp(b1L:(NN(1,g)+b1U),b2L:(NN(2,g)+b2U),b3L:(NN(3,g)+b3U)) ! FELD!!!
  
  INTEGER                ::  i, N1, S1R, N1R, S11R, N11R
  INTEGER                ::  j, N2, S2R, N2R, S22R, N22R
  INTEGER                ::  k, N3, S3R, N3R, S33R, N33R
  
  INTEGER                ::  r, ss
  REAL(8)                   ::  mult
  
  ! TEST!!! herausziehen?
  ! RBGS_mode =< 0 : naturally ordered Gauss-Seidel (slow on CPU, cannot be vectorized)
  ! RBGS_mode == 1 : 3D-Red-Black Gauss-Seidel (faster on CPU due to vectorization, parallelization-independent)
  ! RBGS_mode >= 2 : 2D-Red-Black Gauss-Seidel (normally the fastest variant on CPU due to better cache utilization)
  INTEGER, PARAMETER     ::  RBGS_mode = 2
  LOGICAL, PARAMETER     ::  SOR_yes = .FALSE.
  !REAL(8)   , PARAMETER     ::  omega   = 0.8! 1.2 ! 1.27 !1.27  ! omega = 0.9
  REAL(8), target                   ::  omega
 
  integer                :: jacobi_iterant
  integer,parameter      :: upper_limit_jacobi  = 1
  ! Allgemein, pgf90: "Invariant if transformation" <==> IF-Statements werden aus dem Loop herauskopiert, muss man nicht von Hand machen!
  
  
  !*****************************************************************************************
  N1 = NN(1,g)
  N2 = NN(2,g)
  N3 = NN(3,g)
  
  BC_1L = BC(1,1,g)
  BC_1U = BC(2,1,g)
  BC_2L = BC(1,2,g)
  BC_2U = BC(2,2,g)
  BC_3L = BC(1,3,g)
  BC_3U = BC(2,3,g)
  
  S1R  = SNB(1,1,g)
  S2R  = SNB(1,2,g)
  S3R  = SNB(1,3,g)
  
  N1R  = SNB(2,1,g)
  N2R  = SNB(2,2,g)
  N3R  = SNB(2,3,g)
  
  S11R = SNF(1,1,g)
  S22R = SNF(1,2,g)
  S33R = SNF(1,3,g)
  
  N11R = SNF(2,1,g)
  N22R = SNF(2,2,g)
  N33R = SNF(2,3,g)
  !*****************************************************************************************
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Null-Setzen am Rand nicht notwendig, da Startindizes entsprechend gewÃ¤hlt sind!           !
  !              - Bei der Initialisierung mÃ¼ssen die Intervallgrenzen SiR:NiR anstelle von SiiR:NiiR        !
  !                gewÃ¤hlt werden, da beim  Aufbau der RHS ("vec") innerhalb der Linienrelaxation IMMER auch !
  !                die Randbereiche SiR und NiR aufgebaut aber ggf. mit einer Korrektur wieder Ã¼berschrieben !
  !                werden. Ansonsten wÃ¤re eine Initialisierung allein im Feldbereich SiiR:NiiR ausreichend,  !
  !                kann aber z.B. zu Floating Point Exceptions fÃ¼hren (die fÃ¼r die eigentliche Rechnung      !
  !                allerdings irrelevant wÃ¤ren)!                                                             !
  !              - obere StirnflÃ¤chen werden bei der Initialisierung ebenfalls berÃ¼cksichtigt, da ggf.       !
  !                verschiedene Richtungen nacheinander bearbeitet werden und dies beim Aufbauen der rechten !
  !                Seite sonst berÃ¼cksichtigt werden mÃ¼sste. ==> geringerer Mehraufwand beim Rechnen,        !
  !                weniger Programmierdetails.                                                               !
  !              - LU-Zerlegung (d.h. band, mult) kann gespeichert werden, wenn mindestens eine Richtung     !
  !                Ã¤quidistant ist. Der LÃ¶sungsaufwand wÃ¼rde sich etwa halbieren!!                           !
  !              - "r == 1 .AND. init_yes" sollte idealerweise aus den Schleifen herausgezogen werden, was   !
  !                hier aber aus GrÃ¼nden der Ãbersicht bisher nicht ausgefÃ¼hrt wurde.                        !
  !----------------------------------------------------------------------------------------------------------!
  
  
  omega = 0.98
  
  
  IF (Jacobi_yes) THEN
  
  DO r = 1, n_relax
     
     IF (r == 1 .AND. init_yes) THEN
        !=====================================================================================================
       !=====================================================================================================
        DO k = S33R, N33R
           DO j = S22R, N22R
!pgi$ unroll = n:8
              DO i = S11R, N11R
                 rel(i,j,k) = omega*bb(i,j,k) / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g))
              END DO
           END DO
        END DO
        !=====================================================================================================
       !=====================================================================================================
        
     ELSE
        
        CALL exchange_relax(g,0,0,0,0,.TRUE.,rel)
        
        !=====================================================================================================
       !=====================================================================================================
   if (((N11R .lt. 64) .and. (N22R .lt. 64) .and. (N33R .lt. 64)) .or. (GPU_accelerated == 0)) then


    do jacobi_iterant =1, upper_limit_jacobi

        DO k = S33R, N33R
           DO j = S22R, N22R
!pgi$ unroll = n:8
              DO i = S11R, N11R
                 comp(i,j,k) = omega*(bb(i,j,k)                                                &
                             &      - cdg1(-1,i,g)*rel(i-1,j,k) - cdg1(1,i,g)*rel(i+1,j,k)     &
                             &      - cdg2(-1,j,g)*rel(i,j-1,k) - cdg2(1,j,g)*rel(i,j+1,k)     &
                             &      - cdg3(-1,k,g)*rel(i,j,k-1) - cdg3(1,k,g)*rel(i,j,k+1))    &
                             &      / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g)) + (1.-omega)*rel(i,j,k)
              END DO
           END DO
        END DO

    end do

!    endif


    else
    if (GPU_verbose == 1)    print*, '---------- ---------- -----------on GPU Jacobi is Inverse Relaxing ' 
          call device_jacobi_mg_relaxation(g,omega, rel , bb, comp)
    endif



        !=====================================================================================================
       !=====================================================================================================
        
        rel(1:N1,1:N2,1:N3) = comp(1:N1,1:N2,1:N3)
        
     END IF
  end do   
  !***********************************************************************************************************
  !***********************************************************************************************************
  !***********************************************************************************************************
  
  END IF
  
  
  
  
  END SUBROUTINE relaxation_div_grad_inv
  
  
  
  
  
  
  
  
  
  
 
  
  
 END MODULE mod_laplace 
