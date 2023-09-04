!**************************************************************************************************                                           
!* by Hadi Zolfaghari, ARTORG Center, University of Bern, (zolfaghari.haadi@gmail.com)            *      
!* October 2015 - March 2020                                                                      *            
!**************************************************************************************************
subroutine sgemm(transa,transb,m,n,k,alpha,&

        a,lda,b,ldb,beta,c,ldc)

      

      use iso_c_binding

      implicit none

      

!     .. Scalar Arguments ..

      real(8) alpha,beta

      integer k,lda,ldb,ldc,m,n

      character(len=*) transa,transb

      type(c_ptr):: devA, devB,devC

      integer retVal,retValA,retValB,retValC



!     .. Array Arguments ..

      real(8) a(lda,*),b(ldb,*),c(ldc,*)



! Define the INTERFACE to the NVIDIA C code cublasSgemm and to other helper 

! functions (cublasAlloc,cublasFree, cublasSetMatrix,cublasGetMatrix)



      character(1,c_char) cta, ctb

      interface

!

!      void cublasSgemm (char transa, char transb, int m, int n,

!                        int k, float alpha, const float *A, int lda,

!                  const float *B, int ldb, float beta, float *C, int ldc)

!

       subroutine c_sgemm(cta, ctb, m, n, k,&

         alpha, A, lda, B, ldb, beta, c, ldc) bind(C,name='cublasSgemm')

             use iso_c_binding

             character(1,c_char),value :: cta, ctb

             integer(c_int),value :: m,n,k,lda,ldb,ldc

             real(8)(c_float),value :: alpha,beta

             type(c_ptr),value :: A,B,C

       end subroutine c_sgemm



!

!     cublasStatus CUBLASAPI cublasAlloc (int n, int elemSize, void **devicePtr);

!

       integer  (c_int) function cublasAlloc(n,size,ptr) bind(C,name='cublasAlloc')

         use iso_c_binding

         integer(c_int), value:: n, size

         integer(c_int) :: retVal

         type(c_ptr) :: ptr

       end function cublasAlloc



!

!     cublasStatus CUBLASAPI cublasFree (const void *devicePtr);

!

       integer  (c_int) function cublasFree(ptr) bind(C,name='cublasFree')

         use iso_c_binding

         integer(c_int) :: retVal

         type(c_ptr),value :: ptr

       end function cublasFree

         

!

!    cublasStatus CUBLASAPI cublasSetMatrix (int rows, int cols, int elemSize,

!                                            const void *A, int lda, void *B,

!                                            int ldb);



       integer  (c_int) function cublasSetMatrix(rows,cols,elemSize,A,lda,devA,lda2) bind(C,name='cublasSetMatrix')

         use iso_c_binding

         integer(c_int), value :: rows,cols, elemSize,lda,lda2

         integer(c_int) :: retVal

         real(8)(c_float) :: A(lda,*)

         type(c_ptr), value :: devA

       end function cublasSetMatrix



!

!cublasStatus CUBLASAPI cublasGetMatrix (int rows, int cols, int elemSize,

!                                        const void *A, int lda, void *B,

!                                        int ldb);

!



        integer  (c_int) function cublasGetMatrix(rows,cols,elemSize,devA,lda,B,ldb) bind(C,name='cublasGetMatrix')

         use iso_c_binding

         integer(c_int), value:: rows,cols, elemSize,lda,ldb

         integer(c_int) :: retVal

         type(c_ptr), value:: devA

         real(8)(c_float):: B(ldb,*)

       end function cublasGetMatrix



      end interface 



!

! The calculation, including memory initialization and finalization,



      cta=transa(1:1); ctb=transb(1:1) ! Pass only first character.     



!

!   Allocate the memory on the GPU

!

      retValA=cublasAlloc(m*k,4,devA)

      retValB=cublasAlloc(k*n,4,devB)

      retValC=cublasAlloc(m*n,4,devC)

      if ( (retValA .ne. 0) .or. (retValB .ne. 0) .or. (retValC .ne. 0) ) print *,"Error in memory allocation"



!

!   Move the data from CPU memory to GPU memory

!

      retValA=cublasSetMatrix(m,k,4,A,lda,devA,m)

      retValB=cublasSetMatrix(k,n,4,B,ldb,devB,k)

      retValC=cublasSetMatrix(m,n,4,C,ldc,devC,m)

      if ( (retValA .ne. 0) .or. (retValB .ne. 0) .or. (retValC .ne. 0) ) print *,"Error in memory copy from CPU to GPU"



      call c_sgemm(cta, ctb, m, n, k, alpha, devA, m, devB, k, beta, devC, m)



!

!   Move the result from GPU memory to CPU memory

!

      retVal=cublasGetMatrix(m,n,4,devC,m,C,ldc)

      if ( (retVal .ne. 0) ) print *,"Error in memory copy from GPU to CPU"



!

!   Free the memory on the GPU

!

      retVal=cublasFree(devA)

      retVal=cublasFree(devB)

      retVal=cublasFree(devC)



      return     

      end
