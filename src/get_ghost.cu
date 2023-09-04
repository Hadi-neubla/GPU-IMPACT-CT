extern "C" __global__ void get_ghost(int nx, int ny, int nz, double* rel, double* rear_ghost, double* front_ghost, double* upper_ghost, double* lower_ghost, double* west_ghost, double* east_ghost)



{
for (int line_iter=1; line_iter<=1;line_iter++){

        int i = threadIdx.x + blockIdx.x * blockDim.x;
        int j = threadIdx.y + blockIdx.y * blockDim.y;
        int k = threadIdx.z + blockIdx.z * blockDim.z;



//!!//!!// The serial loop for testing the copy




//          for (int ii= 4; ii <= nx-5; ii++){  
//              for (int jj= 4; jj <= ny-5; jj++){
//                  for (int kk= 4; kk <= nz-5; kk++){
//
//                      comp[ii + jj * nx + kk * nx * ny] = 1.0  * ( bb[ii + jj * nx + kk * nx * ny]
//                                                               -  cdg1[(ii - 4) * nx_cdg + (grid_id - 1) * ny_cdg1 * nx_cdg] * rel[ii - 1 + jj * nx + kk * nx * ny]
//                                                               -  cdg1[(ii - 4) * nx_cdg + 2 + (grid_id - 1) * ny_cdg1 * nx_cdg] * rel[ii + 1 + jj * nx + kk * nx * ny]
//                                                               -  cdg2[(jj - 4) * nx_cdg + (grid_id - 1) * ny_cdg2 * nx_cdg] * rel[ii  + (jj - 1) * nx + kk * nx * ny]
//                                                               -  cdg2[(jj - 4) * nx_cdg + 2 + (grid_id - 1) * ny_cdg2 * nx_cdg] * rel[ii + (jj + 1) * nx + kk * nx * ny]
//                                                               -  cdg3[(kk - 4) * nx_cdg + (grid_id - 1) * ny_cdg3 * nx_cdg] * rel[ii  + jj * nx + (kk - 1) * nx * ny]
//                                                               -  cdg3[(kk - 4) * nx_cdg + 2 + (grid_id - 1) * ny_cdg3 * nx_cdg] * rel[ii + jj * nx + (kk + 1) * nx * ny])
//                                                               / (cdg1[(ii - 4) * nx_cdg + (grid_id - 1) * ny_cdg1 * nx_cdg + 1]
//                                                                + cdg2[(jj - 4) * nx_cdg + (grid_id - 1) * ny_cdg2 * nx_cdg + 1]
//                                                                + cdg3[(kk - 4) * nx_cdg + (grid_id - 1) * ny_cdg3 * nx_cdg + 1])
//                                                                + (1.0 - 1.0) * rel[ii + jj * nx + kk * nx * ny];
//                                                              
//                  }
//              }
//          }



//        DO k = S33R, N33R
//           DO j = S22R, N22R
//!pgi$ unroll = n:8
//              DO i = S11R, N11R
//                 comp(i,j,k) = omega*(bb(i,j,k)                                                &
//                             &      - cdg1(-1,i,g)*rel(i-1,j,k) - cdg1(1,i,g)*rel(i+1,j,k)     &
//                             &      - cdg2(-1,j,g)*rel(i,j-1,k) - cdg2(1,j,g)*rel(i,j+1,k)     &
//                             &      - cdg3(-1,k,g)*rel(i,j,k-1) - cdg3(1,k,g)*rel(i,j,k+1))    &
//                             &      / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g)) + (1.-omega)*rel(i,j,k)
//              END DO
//           END DO
//        END DO






//!//!//! this is the shifted operation -- to be converted to power 2 operation for coalesed memory call on GPU

///    if ((i < 4) || (i > (nx - 4)) || (j < 4) || (j > (ny - 4)) || (k < 4) || (k > (nz - 4))) return;
///
///
///        if ((i >= 4) && (i <= (nx-5))) {
///            if ((j >= 4) && (j <= (ny-5))) {     
///                if ((k >= 4) && (k <= (nz-5))) {
///
///
///
///                  comp[i + j * nx + k * nx * ny] =          ( bb[i + j * nx + k * nx * ny]
///                                                           -  cdg1[(i - 4) * nx_cdg + (grid_id - 1) * ny_cdg1 * nx_cdg] * rel[i - 1 + j * nx + k * nx * ny]
///                                                           -  cdg1[(i - 4) * nx_cdg + 2 + (grid_id - 1) * ny_cdg1 * nx_cdg] * rel[i + 1 + j * nx + k * nx * ny]
///                                                           -  cdg2[(j - 4) * nx_cdg + (grid_id - 1) * ny_cdg2 * nx_cdg] * rel[i  + (j - 1) * nx + k * nx * ny]
///                                                           -  cdg2[(j - 4) * nx_cdg + 2 + (grid_id - 1) * ny_cdg2 * nx_cdg] * rel[i + (j + 1) * nx + k * nx * ny]
///                                                           -  cdg3[(k - 4) * nx_cdg + (grid_id - 1) * ny_cdg3 * nx_cdg] * rel[i  + j * nx + (k - 1) * nx * ny]
///                                                           -  cdg3[(k - 4) * nx_cdg + 2 + (grid_id - 1) * ny_cdg3 * nx_cdg] * rel[i + j * nx + (k + 1) * nx * ny])
///                                                           / (cdg1[(i - 4) * nx_cdg + (grid_id - 1) * ny_cdg1 * nx_cdg + 1]
///                                                            + cdg2[(j - 4) * nx_cdg + (grid_id - 1) * ny_cdg2 * nx_cdg + 1]
///                                                            + cdg3[(k - 4) * nx_cdg + (grid_id - 1) * ny_cdg3 * nx_cdg + 1]) ;
/////                                                            + (1.0 - 1.0) * rel[i + j * nx + k * nx * ny];
///
///
///                 }
///               }
///            }







//!//!//! this is the power 2 implementation
        int indg_rf = i + j * (nx - 8);
        int indg_ew = j + k * (ny - 8);
        int indg_lu = i + k * (nx - 8);
        int ind = (i + 4) + (j + 4) * nx + (k + 4) * nx * ny;
        if ((i >= 1) && (i <= (nx-10))) {
            if ((j >= 1) && (j <= (ny-10))) {
               if (k == 0) {
//               rel[ind] =     omega *    ( bb[ind]
//                                                           -  cdg1[(i) * nx_cdg] * rel[ind-1]
//                                                           -  cdg1[(i) * nx_cdg + 2 ] * rel[ind+1]
//                                                           -  cdg2[(j) * nx_cdg ] * rel[ind-nx]
//                                                           -  cdg2[(j) * nx_cdg + 2 ] * rel[ind+nx]
//                                                           -  cdg3[(k) * nx_cdg  ] *  rear_ghost[indg_rf]                     //rel[ind -   nx * ny]
//                                                           -  cdg3[(k) * nx_cdg + 2 ] * rel[ind +  nx * ny])
//                                                           / (cdg1[(i) * nx_cdg + 1]
//                                                            + cdg2[(j) * nx_cdg + 1]
//                                                            + cdg3[(k) * nx_cdg + 1])
//                                                                 + (1.0 - omega) * rel[ind];
                  rear_ghost[indg_rf] = rel[ind];
                  }

               // front ghosts
                 if (k == (nz-9)) {
                 // int ind = (i + 4) + (j + 4) * nx + (k + 4) * nx * ny;
                 // int indg = i + j * (nx - 7);
//                  rel[ind] =     omega *    ( bb[ind]
//                                                           -  cdg1[(i) * nx_cdg] * rel[ind-1]
//                                                           -  cdg1[(i) * nx_cdg + 2 ] * rel[ind+1]
//                                                           -  cdg2[(j) * nx_cdg ] * rel[ind-nx]
//                                                           -  cdg2[(j) * nx_cdg + 2 ] * rel[ind+nx]
//                                                           -  cdg3[(k) * nx_cdg  ] * rel[ind -   nx * ny]
//                                                           -  cdg3[(k) * nx_cdg + 2 ] * front_ghost[indg_rf])                  // rel[ind +  nx * ny])
//                                                           / (cdg1[(i) * nx_cdg + 1]
//                                                            + cdg2[(j) * nx_cdg + 1]
//                                                            + cdg3[(k) * nx_cdg + 1])
//                                                                 + (1.0 - omega) * rel[ind];
                  front_ghost[indg_rf] = rel[ind];
                  }
             }
       }

      if ((i >= 1) && (i <= (nx-10))) {
         if ((k >= 1) && (k <= (nz-10))) {
             // lower ghost
              if (j==0) {
                 // int ind = (i + 4) + (j + 4) * nx + (k + 4) * nx * ny;
                 // int indg = i  + k * (nx - 7);
//                  rel[ind] =     omega *    ( bb[ind]
//                                                           -  cdg1[(i) * nx_cdg] * rel[ind-1]
//                                                           -  cdg1[(i) * nx_cdg + 2 ] * rel[ind+1]
//                                                           -  cdg2[(j) * nx_cdg ] *  lower_ghost[indg_lu]                     // rel[ind-nx]
//                                                           -  cdg2[(j) * nx_cdg + 2 ] * rel[ind+nx]
//                                                           -  cdg3[(k) * nx_cdg  ] * rel[ind -   nx * ny]
//                                                           -  cdg3[(k) * nx_cdg + 2 ] * rel[ind +  nx * ny])
//                                                           / (cdg1[(i) * nx_cdg + 1]
//                                                            + cdg2[(j) * nx_cdg + 1]
//                                                            + cdg3[(k) * nx_cdg + 1])
//                                                                 + (1.0 - omega) * rel[ind];
             lower_ghost[indg_lu] = rel[ind];
             }
             // upper ghost
             if (j==(ny-9)){
                 // int ind = (i + 4) + (j + 4) * nx + (k + 4) * nx * ny;
                 // int indg = i + k * (nx - 7) ;
//                  rel[ind] =     omega *    ( bb[ind]
//                                                           -  cdg1[(i) * nx_cdg] * rel[ind-1]
//                                                           -  cdg1[(i) * nx_cdg + 2 ] * rel[ind+1]
//                                                           -  cdg2[(j) * nx_cdg ] * rel[ind-nx]
//                                                           -  cdg2[(j) * nx_cdg + 2 ] * upper_ghost[indg_lu]                   // rel[ind+nx]
//                                                           -  cdg3[(k) * nx_cdg  ] * rel[ind -   nx * ny]
//                                                           -  cdg3[(k) * nx_cdg + 2 ] * rel[ind +  nx * ny])
//                                                           / (cdg1[(i) * nx_cdg + 1]
//                                                            + cdg2[(j) * nx_cdg + 1]
//                                                            + cdg3[(k) * nx_cdg + 1])
//                                                                 + (1.0 - omega) * rel[ind];
                 upper_ghost[indg_lu] = rel[ind];
             }
       }
     }

     if ((j >= 1) && (j <= (ny-10))) {
       if ((k >= 1) && (k <= (nz-10))) {

       // west ghost 
       if (i == 0){
                 // int ind = (i + 4) + (j + 4) * nx + (k + 4) * nx * ny;
                 // int indg = j  + k * (ny - 7);
//                  rel[ind] =     omega *    ( bb[ind]
//                                                           -  cdg1[(i) * nx_cdg] *  west_ghost[indg_ew]           //rel[ind-1]
//                                                           -  cdg1[(i) * nx_cdg + 2 ] * rel[ind+1]
//                                                           -  cdg2[(j) * nx_cdg ] * rel[ind-nx]
//                                                           -  cdg2[(j) * nx_cdg + 2 ] * rel[ind+nx]
//                                                           -  cdg3[(k) * nx_cdg  ] * rel[ind -   nx * ny]
//                                                           -  cdg3[(k) * nx_cdg + 2 ] * rel[ind +  nx * ny])
//                                                           / (cdg1[(i) * nx_cdg + 1]
//                                                            + cdg2[(j) * nx_cdg + 1]
//                                                            + cdg3[(k) * nx_cdg + 1])
//                                                                 + (1.0 - omega) * rel[ind];      
                 west_ghost[indg_ew] = rel[ind]; 
      }
       // east ghost
         if (i==(nx-9)){
                 // int ind = (i + 4) + (j + 4) * nx + (k + 4) * nx * ny;
                 // int indg = j + k * (ny - 7);
//                  rel[ind] =     omega *    ( bb[ind]
//                                                           -  cdg1[(i) * nx_cdg] * rel[ind-1]
//                                                           -  cdg1[(i) * nx_cdg + 2 ] * east_ghost[indg_ew]                      //rel[ind+1]
//                                                           -  cdg2[(j) * nx_cdg ] * rel[ind-nx]
//                                                           -  cdg2[(j) * nx_cdg + 2 ] * rel[ind+nx]
//                                                           -  cdg3[(k) * nx_cdg  ] * rel[ind -   nx * ny]
//                                                           -  cdg3[(k) * nx_cdg + 2 ] * rel[ind +  nx * ny])
//                                                           / (cdg1[(i) * nx_cdg + 1]
//                                                            + cdg2[(j) * nx_cdg + 1]
//                                                            + cdg3[(k) * nx_cdg + 1])
//                                                                 + (1.0 - omega) * rel[ind];
                 east_ghost[indg_ew] = rel[ind];
       }
   }
 }
   
/////////////////// handling borders: 2nd part: corners
      ///////////// i=0, j=0
       if ((i == 0) && (j==0)) {
          if ((k == 0)){
        // int ind = (i + 4) + (j + 4) * nx + (k + 4) * nx * ny;
//          rel[ind] =     omega *    ( bb[ind]
//                                                  -  cdg1[(i) * nx_cdg] * west_ghost[indg_ew]
//                                                  -  cdg1[(i) * nx_cdg + 2 ] * rel[ind+1]
//                                                  -  cdg2[(j) * nx_cdg ] * lower_ghost[indg_lu]
//                                                  -  cdg2[(j) * nx_cdg + 2 ] * rel[ind+nx]
//                                                  -  cdg3[(k) * nx_cdg  ] * rear_ghost[indg_rf]
//                                                  -  cdg3[(k) * nx_cdg + 2 ] * rel[ind +  nx * ny])
//                                                  / (cdg1[(i) * nx_cdg + 1]
//                                                   + cdg2[(j) * nx_cdg + 1]
//                                                   + cdg3[(k) * nx_cdg + 1])
//                                                        + (1.0 - omega) * rel[ind];
           west_ghost[indg_ew] = rel[ind];
           lower_ghost[indg_lu] = rel[ind];
           rear_ghost[indg_rf] = rel[ind];
         }
          if (k==(nz-9)) {
        // int ind = (i + 4) + (j + 4) * nx + (k + 4) * nx * ny;
//          rel[ind] =     omega *    ( bb[ind]
//                                                  -  cdg1[(i) * nx_cdg] * west_ghost[indg_ew]
//                                                  -  cdg1[(i) * nx_cdg + 2 ] * rel[ind+1]
//                                                  -  cdg2[(j) * nx_cdg ] * lower_ghost[indg_lu]
//                                                  -  cdg2[(j) * nx_cdg + 2 ] * rel[ind+nx]
//                                                  -  cdg3[(k) * nx_cdg  ] * rel[ind- nx*ny]
//                                                  -  cdg3[(k) * nx_cdg + 2 ] * front_ghost[indg_rf])
//                                                  / (cdg1[(i) * nx_cdg + 1]
//                                                   + cdg2[(j) * nx_cdg + 1]
//                                                   + cdg3[(k) * nx_cdg + 1])
//                                                        + (1.0 - omega) * rel[ind];
           west_ghost[indg_ew] = rel[ind];
           lower_ghost[indg_lu] = rel[ind];
           front_ghost[indg_rf] = rel[ind];
           
         }
//         else  {
//
//        // int ind = (i + 4) + (j + 4) * nx + (k + 4) * nx * ny;
//          rel[ind] =     omega *    ( bb[ind]
//                                                  -  cdg1[(i) * nx_cdg] * west_ghost[indg_ew]
//                                                  -  cdg1[(i) * nx_cdg + 2 ] * rel[ind+1]
//                                                  -  cdg2[(j) * nx_cdg ] * lower_ghost[indg_lu]
//                                                  -  cdg2[(j) * nx_cdg + 2 ] * rel[ind+nx]
//                                                  -  cdg3[(k) * nx_cdg  ] * rel[ind - nx*ny]
//                                                  -  cdg3[(k) * nx_cdg + 2 ] * rel[ind + nx*ny])
//                                                  / (cdg1[(i) * nx_cdg + 1]
//                                                   + cdg2[(j) * nx_cdg + 1]
//                                                   + cdg3[(k) * nx_cdg + 1])
//                                                        + (1.0 - omega) * rel[ind];
//         }
       
      }


      /////////// i=nx-9, j=0
       if ((i == (nx-9)) && (j==0)) {
          if ((k == 0)){
        // int ind = (i + 4) + (j + 4) * nx + (k + 4) * nx * ny;
//          rel[ind] =     omega *    ( bb[ind]
//                                                  -  cdg1[(i) * nx_cdg] * rel[ind-1]
//                                                  -  cdg1[(i) * nx_cdg + 2 ] * east_ghost[indg_ew]
//                                                  -  cdg2[(j) * nx_cdg ] * lower_ghost[indg_lu]
//                                                  -  cdg2[(j) * nx_cdg + 2 ] * rel[ind+nx]
//                                                  -  cdg3[(k) * nx_cdg  ] * rear_ghost[indg_rf]
//                                                  -  cdg3[(k) * nx_cdg + 2 ] * rel[ind +  nx * ny])
//                                                  / (cdg1[(i) * nx_cdg + 1]
//                                                   + cdg2[(j) * nx_cdg + 1]
//                                                   + cdg3[(k) * nx_cdg + 1])
//                                                        + (1.0 - omega) * rel[ind];
           east_ghost[indg_ew] = rel[ind];
           lower_ghost[indg_lu] = rel[ind];
           rear_ghost[indg_rf] = rel[ind];
 
         }
         if (k==(nz-9)) {
        // int ind = (i + 4) + (j + 4) * nx + (k + 4) * nx * ny;
//          rel[ind] =     omega *    ( bb[ind]
//                                                  -  cdg1[(i) * nx_cdg] * rel[ind-1]
//                                                  -  cdg1[(i) * nx_cdg + 2 ] * east_ghost[indg_ew]
//                                                  -  cdg2[(j) * nx_cdg ] * lower_ghost[indg_lu]
//                                                  -  cdg2[(j) * nx_cdg + 2 ] * rel[ind+nx]
//                                                  -  cdg3[(k) * nx_cdg  ] * rel[ind -nx*ny]
//                                                  -  cdg3[(k) * nx_cdg + 2 ] * front_ghost[indg_rf])
//                                                  / (cdg1[(i) * nx_cdg + 1]
//                                                   + cdg2[(j) * nx_cdg + 1]
//                                                   + cdg3[(k) * nx_cdg + 1])
//                                                        + (1.0 - omega) * rel[ind];
           east_ghost[indg_ew] = rel[ind];
           lower_ghost[indg_lu] = rel[ind];
           front_ghost[indg_rf] = rel[ind];
 
         }
//         else  {
//
//        // int ind = (i + 4) + (j + 4) * nx + (k + 4) * nx * ny;
//          rel[ind] =     omega *    ( bb[ind]
//                                                  -  cdg1[(i) * nx_cdg] * rel[ind-1]
//                                                  -  cdg1[(i) * nx_cdg + 2 ] * east_ghost[indg_ew]
//                                                  -  cdg2[(j) * nx_cdg ] * lower_ghost[indg_lu]
//                                                  -  cdg2[(j) * nx_cdg + 2 ] * rel[ind+nx]
//                                                  -  cdg3[(k) * nx_cdg  ] * rel[ind - nx*ny]
//                                                  -  cdg3[(k) * nx_cdg + 2 ] * rel[ind + nx*ny])
//                                                  / (cdg1[(i) * nx_cdg + 1]
//                                                   + cdg2[(j) * nx_cdg + 1]
//                                                   + cdg3[(k) * nx_cdg + 1])
//                                                        + (1.0 - omega) * rel[ind];
//         }
       }
      /////////// i=0, j=ny-9
       if ((i == 0) && (j==(ny-9))) {
          if ((k == 0)){
        // int ind = (i + 4) + (j + 4) * nx + (k + 4) * nx * ny;
//          rel[ind] =     omega *    ( bb[ind]
//                                                  -  cdg1[(i) * nx_cdg] * west_ghost[indg_ew]
//                                                  -  cdg1[(i) * nx_cdg + 2 ] * rel[ind+1]
//                                                  -  cdg2[(j) * nx_cdg ] * rel[ind-nx]
//                                                  -  cdg2[(j) * nx_cdg + 2 ] * upper_ghost[indg_lu]
//                                                  -  cdg3[(k) * nx_cdg  ] * rear_ghost[indg_rf]
//                                                  -  cdg3[(k) * nx_cdg + 2 ] * rel[ind +  nx * ny])
//                                                  / (cdg1[(i) * nx_cdg + 1]
//                                                   + cdg2[(j) * nx_cdg + 1]
//                                                   + cdg3[(k) * nx_cdg + 1])
//                                                        + (1.0 - omega) * rel[ind];
           west_ghost[indg_ew] = rel[ind];
           upper_ghost[indg_lu] = rel[ind];
           rear_ghost[indg_rf] = rel[ind];
 
         }
          if (k==(nz-9)) {

        // int ind = (i + 4) + (j + 4) * nx + (k + 4) * nx * ny;
//          rel[ind] =     omega *    ( bb[ind]
//                                                  -  cdg1[(i) * nx_cdg] * west_ghost[indg_ew]
//                                                  -  cdg1[(i) * nx_cdg + 2 ] * rel[ind+1]
//                                                  -  cdg2[(j) * nx_cdg ] * rel[ind-nx]
//                                                  -  cdg2[(j) * nx_cdg + 2 ] * upper_ghost[indg_lu]
//                                                  -  cdg3[(k) * nx_cdg  ] * rel[ind -nx*ny]
//                                                  -  cdg3[(k) * nx_cdg + 2 ] * front_ghost[indg_rf])
//                                                  / (cdg1[(i) * nx_cdg + 1]
//                                                   + cdg2[(j) * nx_cdg + 1]
//                                                   + cdg3[(k) * nx_cdg + 1])
//                                                        + (1.0 - omega) * rel[ind];
           west_ghost[indg_ew] = rel[ind];
           upper_ghost[indg_lu] = rel[ind];
           front_ghost[indg_rf] = rel[ind];
 
         }
//         else  {
//
//        // int ind = (i + 4) + (j + 4) * nx + (k + 4) * nx * ny;
//          rel[ind] =     omega *    ( bb[ind]
//                                                  -  cdg1[(i) * nx_cdg] * west_ghost[indg_ew]
//                                                  -  cdg1[(i) * nx_cdg + 2 ] * rel[ind+1]
//                                                  -  cdg2[(j) * nx_cdg ] * lower_ghost[indg_lu]
//                                                  -  cdg2[(j) * nx_cdg + 2 ] * rel[ind+nx]
//                                                  -  cdg3[(k) * nx_cdg  ] * rel[ind - nx*ny]
//                                                  -  cdg3[(k) * nx_cdg + 2 ] * rel[ind + nx*ny])
//                                                  / (cdg1[(i) * nx_cdg + 1]
//                                                   + cdg2[(j) * nx_cdg + 1]
//                                                   + cdg3[(k) * nx_cdg + 1])
//                                                        + (1.0 - omega) * rel[ind];
//         }
       }


      /////////// i=nx-9, j=ny-9
       if ((i == (nx-9)) && (j==(ny-9))) {
          if ((k == 0)){
        // int ind = (i + 4) + (j + 4) * nx + (k + 4) * nx * ny;
//          rel[ind] =     omega *    ( bb[ind]
//                                                  -  cdg1[(i) * nx_cdg] * rel[ind-1]
//                                                  -  cdg1[(i) * nx_cdg + 2 ] * east_ghost[indg_ew]
//                                                  -  cdg2[(j) * nx_cdg ] * rel[ind-nx]
//                                                  -  cdg2[(j) * nx_cdg + 2 ] * upper_ghost[indg_lu]
//                                                  -  cdg3[(k) * nx_cdg  ] * rear_ghost[indg_rf]
//                                                  -  cdg3[(k) * nx_cdg + 2 ] * rel[ind +  nx * ny])
//                                                  / (cdg1[(i) * nx_cdg + 1]
//                                                   + cdg2[(j) * nx_cdg + 1]
//                                                   + cdg3[(k) * nx_cdg + 1])
//                                                        + (1.0 - omega) * rel[ind];
           east_ghost[indg_ew] = rel[ind];
           upper_ghost[indg_lu] = rel[ind];
           rear_ghost[indg_rf] = rel[ind];
 
         }
         if (k==(nz-9)) {

        // int ind = (i + 4) + (j + 4) * nx + (k + 4) * nx * ny;
//          rel[ind] =     omega *    ( bb[ind]
//                                                  -  cdg1[(i) * nx_cdg] * rel[ind-1]
//                                                  -  cdg1[(i) * nx_cdg + 2 ] * east_ghost[indg_ew]
//                                                  -  cdg2[(j) * nx_cdg ] * rel[ind-nx]
//                                                  -  cdg2[(j) * nx_cdg + 2 ] * upper_ghost[indg_lu]
//                                                  -  cdg3[(k) * nx_cdg  ] * rel[ind -nx*ny]
//                                                  -  cdg3[(k) * nx_cdg + 2 ] * front_ghost[indg_rf])
//                                                  / (cdg1[(i) * nx_cdg + 1]
//                                                   + cdg2[(j) * nx_cdg + 1]
//                                                   + cdg3[(k) * nx_cdg + 1])
//                                                        + (1.0 - omega) * rel[ind];
           east_ghost[indg_ew] = rel[ind];
           upper_ghost[indg_lu] = rel[ind];
           front_ghost[indg_rf] = rel[ind];
 
                        }
//         else  {
//
//        // int ind = (i + 4) + (j + 4) * nx + (k + 4) * nx * ny;
//          rel[ind] =     omega *    ( bb[ind]
//                                                  -  cdg1[(i) * nx_cdg] * west_ghost[indg_ew]
//                                                  -  cdg1[(i) * nx_cdg + 2 ] * rel[ind+1]
//                                                  -  cdg2[(j) * nx_cdg ] * rel[ind-nx]
//                                                  -  cdg2[(j) * nx_cdg + 2 ] * upper_ghost[indg_lu]
//                                                  -  cdg3[(k) * nx_cdg  ] * rel[ind - nx*ny]
//                                                  -  cdg3[(k) * nx_cdg + 2 ] * rel[ind + nx*ny])
//                                                  / (cdg1[(i) * nx_cdg + 1]
//                                                   + cdg2[(j) * nx_cdg + 1]
//                                                   + cdg3[(k) * nx_cdg + 1])
//                                                        + (1.0 - omega) * rel[ind];
//         }
       }



//////////////////////////// handling borders, third part: bars
////////////////////// xz
      /////////// i=0, k=nz-9
       if ((i == 0) && (k==(nz-9))) {
          if ((j > 0) && (j < (ny-9))){ 
        // int ind = (i + 4) + (j + 4) * nx + (k + 4) * nx * ny;
//          rel[ind] =     omega *    ( bb[ind]
//                                                  -  cdg1[(i) * nx_cdg] * west_ghost[indg_ew]
//                                                  -  cdg1[(i) * nx_cdg + 2 ] * rel[ind+1]
//                                                  -  cdg2[(j) * nx_cdg ] * rel[ind-nx]
//                                                  -  cdg2[(j) * nx_cdg + 2 ] * rel[ind+nx]
//                                                  -  cdg3[(k) * nx_cdg  ] * rel[ind - nx*ny]
//                                                  -  cdg3[(k) * nx_cdg + 2 ] * front_ghost[indg_rf])
//                                                  / (cdg1[(i) * nx_cdg + 1]
//                                                   + cdg2[(j) * nx_cdg + 1]
//                                                   + cdg3[(k) * nx_cdg + 1])
//                                                        + (1.0 - omega) * rel[ind];
           west_ghost[indg_ew] = rel[ind];
           front_ghost[indg_rf] = rel[ind];
 
         }
       }


      /////////// i=0, k=0
       if ((i == 0) && (k==0)) {
          if ((j > 0) && (j < (ny-9))){
        // int ind = (i + 4) + (j + 4) * nx + (k + 4) * nx * ny;
//          rel[ind] =     omega *    ( bb[ind]
//                                                  -  cdg1[(i) * nx_cdg] * west_ghost[indg_ew]
//                                                  -  cdg1[(i) * nx_cdg + 2 ] * rel[ind+1]
//                                                  -  cdg2[(j) * nx_cdg ] * rel[ind-nx]
//                                                  -  cdg2[(j) * nx_cdg + 2 ] * rel[ind+nx]
//                                                  -  cdg3[(k) * nx_cdg  ] * rear_ghost[indg_rf]
//                                                  -  cdg3[(k) * nx_cdg + 2 ] * rel[ind + nx*ny])
//                                                  / (cdg1[(i) * nx_cdg + 1]
//                                                   + cdg2[(j) * nx_cdg + 1]
//                                                   + cdg3[(k) * nx_cdg + 1])
//                                                        + (1.0 - omega) * rel[ind];
           west_ghost[indg_ew] = rel[ind];
           rear_ghost[indg_rf] = rel[ind];
 
         }
       }

       if ((i == (nx-9)) && (k==0)) {
          if ((j >  0) && (j < (ny-9))){
        // int ind = (i + 4) + (j + 4) * nx + (k + 4) * nx * ny;
//          rel[ind] =     omega *    ( bb[ind]
//                                                  -  cdg1[(i) * nx_cdg] * rel[ind-1]
//                                                  -  cdg1[(i) * nx_cdg + 2 ] * east_ghost[indg_ew]
//                                                  -  cdg2[(j) * nx_cdg ] * rel[ind-nx]
//                                                  -  cdg2[(j) * nx_cdg + 2 ] * rel[ind+nx]
//                                                  -  cdg3[(k) * nx_cdg  ] * rear_ghost[indg_rf]
//                                                  -  cdg3[(k) * nx_cdg + 2 ] * rel[ind + nx*ny])
//                                                  / (cdg1[(i) * nx_cdg + 1]
//                                                   + cdg2[(j) * nx_cdg + 1]
//                                                   + cdg3[(k) * nx_cdg + 1])
//                                                        + (1.0 - omega) * rel[ind];
           east_ghost[indg_ew] = rel[ind];
           rear_ghost[indg_rf] = rel[ind];
 

         }
       }
       if ((i == (nx-9)) && (k==(nz-9))) {
          if ((j > 0) && (j < (ny-9))){
        // int ind = (i + 4) + (j + 4) * nx + (k + 4) * nx * ny;
//          rel[ind] =     omega *    ( bb[ind]
//                                                  -  cdg1[(i) * nx_cdg] * rel[ind-1]
//                                                  -  cdg1[(i) * nx_cdg + 2 ] * east_ghost[indg_ew]
//                                                  -  cdg2[(j) * nx_cdg ] * rel[ind-nx]
//                                                  -  cdg2[(j) * nx_cdg + 2 ] * rel[ind+nx]
//                                                  -  cdg3[(k) * nx_cdg  ] * rel[ind- nx*ny]
//                                                  -  cdg3[(k) * nx_cdg + 2 ] * front_ghost[indg_rf])
//                                                  / (cdg1[(i) * nx_cdg + 1]
//                                                   + cdg2[(j) * nx_cdg + 1]
//                                                   + cdg3[(k) * nx_cdg + 1])
//                                                        + (1.0 - omega) * rel[ind];
//
           east_ghost[indg_ew] = rel[ind];
           front_ghost[indg_rf] = rel[ind];
         }
       }

////////////////////// yz
      /////////// j=0, k=nz-9
       if ((j == 0) && (k==(nz-9))) {
          if ((i > 0) && (i < (nx-9))){
        // int ind = (i + 4) + (j + 4) * nx + (k + 4) * nx * ny;
//          rel[ind] =     omega *    ( bb[ind]
//                                                  -  cdg1[(i) * nx_cdg] * rel[ind-1]
//                                                  -  cdg1[(i) * nx_cdg + 2 ] * rel[ind+1]
//                                                  -  cdg2[(j) * nx_cdg ] * lower_ghost[indg_lu]
//                                                  -  cdg2[(j) * nx_cdg + 2 ] * rel[ind+nx]
//                                                  -  cdg3[(k) * nx_cdg  ] * rel[ind - nx*ny]
//                                                  -  cdg3[(k) * nx_cdg + 2 ] * front_ghost[indg_rf])
//                                                  / (cdg1[(i) * nx_cdg + 1]
//                                                   + cdg2[(j) * nx_cdg + 1]
//                                                   + cdg3[(k) * nx_cdg + 1])
//                                                        + (1.0 - omega) * rel[ind];
//
           lower_ghost[indg_lu] = rel[ind];
           front_ghost[indg_rf] = rel[ind];
         }
       }
       if ((j == 0) && (k == 0)) {
          if ((i > 0) && (i < (nx-9))){
        // int ind = (i + 4) + (j + 4) * nx + (k + 4) * nx * ny;
//          rel[ind] =     omega *    ( bb[ind]
//                                                  -  cdg1[(i) * nx_cdg] * rel[ind-1]
//                                                  -  cdg1[(i) * nx_cdg + 2 ] * rel[ind+1]
//                                                  -  cdg2[(j) * nx_cdg ] * lower_ghost[indg_lu]
//                                                  -  cdg2[(j) * nx_cdg + 2 ] * rel[ind+nx]
//                                                  -  cdg3[(k) * nx_cdg  ] * rear_ghost[indg_rf]
//                                                  -  cdg3[(k) * nx_cdg + 2 ] * rel[ind + nx * ny])
//                                                  / (cdg1[(i) * nx_cdg + 1]
//                                                   + cdg2[(j) * nx_cdg + 1]
//                                                   + cdg3[(k) * nx_cdg + 1])
//                                                        + (1.0 - omega) * rel[ind];
//
           lower_ghost[indg_lu] = rel[ind];
           rear_ghost[indg_rf] = rel[ind];
         }
       }
       if ((j == (ny-9)) && (k == 0)) {
          if ((i > 0) && (i < (nx-9))){
        // int ind = (i + 4) + (j + 4) * nx + (k + 4) * nx * ny;
//          rel[ind] =     omega *    ( bb[ind]
//                                                  -  cdg1[(i) * nx_cdg] * rel[ind-1]
//                                                  -  cdg1[(i) * nx_cdg + 2 ] * rel[ind+1]
//                                                  -  cdg2[(j) * nx_cdg ] * rel[ind-nx]
//                                                  -  cdg2[(j) * nx_cdg + 2 ] * upper_ghost[indg_lu]
//                                                  -  cdg3[(k) * nx_cdg  ] * rear_ghost[indg_rf]
//                                                  -  cdg3[(k) * nx_cdg + 2 ] * rel[ind + nx * ny])
//                                                  / (cdg1[(i) * nx_cdg + 1]
//                                                   + cdg2[(j) * nx_cdg + 1]
//                                                   + cdg3[(k) * nx_cdg + 1])
//                                                        + (1.0 - omega) * rel[ind];
           upper_ghost[indg_lu] = rel[ind];
           rear_ghost[indg_rf] = rel[ind];
         }
       }

       if ((j == (ny-9)) && (k == (nz-9))) {
          if ((i >0) && (i <(nx-9))){
        // int ind = (i + 4) + (j + 4) * nx + (k + 4) * nx * ny;
//          rel[ind] =     omega *    ( bb[ind]
//                                                  -  cdg1[(i) * nx_cdg] * rel[ind-1]
//                                                  -  cdg1[(i) * nx_cdg + 2 ] * rel[ind+1]
//                                                  -  cdg2[(j) * nx_cdg ] * rel[ind-nx]
//                                                  -  cdg2[(j) * nx_cdg + 2 ] * upper_ghost[indg_lu]
//                                                  -  cdg3[(k) * nx_cdg  ] * rel[ind - nx * ny]
//                                                  -  cdg3[(k) * nx_cdg + 2 ] * front_ghost[indg_rf])
//                                                  / (cdg1[(i) * nx_cdg + 1]
//                                                   + cdg2[(j) * nx_cdg + 1]
//                                                   + cdg3[(k) * nx_cdg + 1])
//                                                        + (1.0 - omega) * rel[ind];
           upper_ghost[indg_lu] = rel[ind];
           front_ghost[indg_rf] = rel[ind];
         }
       }
////////////////////// xy
      /////////// i=0, j=ny-9
       if ((i == 0) && (j==(ny-9))) {
          if ((k > 0) && (k < (nz-9))){ 
        // int ind = (i + 4) + (j + 4) * nx + (k + 4) * nx * ny;
//          rel[ind] =     omega *    ( bb[ind]
//                                                  -  cdg1[(i) * nx_cdg] * west_ghost[indg_ew]
//                                                  -  cdg1[(i) * nx_cdg + 2 ] * rel[ind+1]
//                                                  -  cdg2[(j) * nx_cdg ] * rel[ind-nx]
//                                                  -  cdg2[(j) * nx_cdg + 2 ] * upper_ghost[indg_lu]
//                                                  -  cdg3[(k) * nx_cdg  ] * rel[ind - nx*ny]
//                                                  -  cdg3[(k) * nx_cdg + 2 ] * rel[ind + nx*ny])
//                                                  / (cdg1[(i) * nx_cdg + 1]
//                                                   + cdg2[(j) * nx_cdg + 1]
//                                                   + cdg3[(k) * nx_cdg + 1])
//                                                        + (1.0 - omega) * rel[ind];
           upper_ghost[indg_lu] = rel[ind];
           west_ghost[indg_ew] = rel[ind];
         }
       }


      /////////// i=0, j=0
       if ((i == 0) && (j==0)) {
          if ((k > 0) && (k < (nz-9))){
        // int ind = (i + 4) + (j + 4) * nx + (k + 4) * nx * ny;
//          rel[ind] =     omega *    ( bb[ind]
//                                                  -  cdg1[(i) * nx_cdg] * west_ghost[indg_ew]
//                                                  -  cdg1[(i) * nx_cdg + 2 ] * rel[ind+1]
//                                                  -  cdg2[(j) * nx_cdg ] * lower_ghost[indg_lu]
//                                                  -  cdg2[(j) * nx_cdg + 2 ] * rel[ind+nx]
//                                                  -  cdg3[(k) * nx_cdg  ] * rel[ind - nx*ny]
//                                                  -  cdg3[(k) * nx_cdg + 2 ] * rel[ind + nx*ny])
//                                                  / (cdg1[(i) * nx_cdg + 1]
//                                                   + cdg2[(j) * nx_cdg + 1]
//                                                   + cdg3[(k) * nx_cdg + 1])
//                                                        + (1.0 - omega) * rel[ind];
           west_ghost[indg_ew] = rel[ind];
           lower_ghost[indg_lu] = rel[ind];
         }
       }
      /////////// i=nx-9 , j=0
       if ((i == (nx-9)) && (j==0)) {
          if ((k >  0) && (k < (nz-9))){
        // int ind = (i + 4) + (j + 4) * nx + (k + 4) * nx * ny;
//          rel[ind] =     omega *    ( bb[ind]
//                                                  -  cdg1[(i) * nx_cdg] * rel[ind-1]
//                                                  -  cdg1[(i) * nx_cdg + 2 ] * east_ghost[indg_ew]
//                                                  -  cdg2[(j) * nx_cdg ] * lower_ghost[indg_lu]
//                                                  -  cdg2[(j) * nx_cdg + 2 ] * rel[ind+nx]
//                                                  -  cdg3[(k) * nx_cdg  ] * rel[ind - nx*ny]
//                                                  -  cdg3[(k) * nx_cdg + 2 ] * rel[ind + nx*ny])
//                                                  / (cdg1[(i) * nx_cdg + 1]
//                                                   + cdg2[(j) * nx_cdg + 1]
//                                                   + cdg3[(k) * nx_cdg + 1])
//                                                        + (1.0 - omega) * rel[ind];
           lower_ghost[indg_lu] = rel[ind];
           east_ghost[indg_ew] = rel[ind];
         }
       }
       /////////// i =nx-9, j=ny-9
       if ((i == (nx-9)) && (j==(ny-9))) {
          if ((k > 0) && (k < (nz-9))){
        // int ind = (i + 4) + (j + 4) * nx + (k + 4) * nx * ny;
//          rel[ind] =     omega *    ( bb[ind]
//                                                  -  cdg1[(i) * nx_cdg] * rel[ind-1]
//                                                  -  cdg1[(i) * nx_cdg + 2 ] * east_ghost[indg_ew]
//                                                  -  cdg2[(j) * nx_cdg ] * rel[ind-nx]
//                                                  -  cdg2[(j) * nx_cdg + 2 ] * upper_ghost[indg_lu]
//                                                  -  cdg3[(k) * nx_cdg  ] * rel[ind- nx*ny]
//                                                  -  cdg3[(k) * nx_cdg + 2 ] * rel[ind + nx*ny])
//                                                  / (cdg1[(i) * nx_cdg + 1]
//                                                   + cdg2[(j) * nx_cdg + 1]
//                                                   + cdg3[(k) * nx_cdg + 1])
//                                                        + (1.0 - omega) * rel[ind];
//
           upper_ghost[indg_lu] = rel[ind];
           east_ghost[indg_ew] = rel[ind];
         }
       }



//cudaThreadSynchronize(); 
}
}



