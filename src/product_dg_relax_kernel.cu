// !**************************************************************************************************                                           
// !* by Hadi Zolfaghari, ARTORG Center, University of Bern, (zolfaghari.haadi@gmail.com)            *      
// !* October 2015 - March 2020                                                                      *            
// !**************************************************************************************************
extern "C" __global__ void product_dg_relax_kernel(int grid_id, int nx, int ny, int nz, int ny_cdg1, int ny_cdg2, int ny_cdg3, int nx_cdg, double* cdg1, double* cdg2, double* cdg3, double* rel, double* comp)



{
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

                  comp[ (i + 4) + (j + 4) * nx + (k + 4) * nx * ny] =   cdg1[(i) * nx_cdg + (grid_id - 1) * ny_cdg1 * nx_cdg] * rel[(i + 4) - 1 + (j + 4) * nx + (k + 4) * nx * ny]
                                                                     +  cdg1[(i) * nx_cdg + 2 + (grid_id - 1) * ny_cdg1 * nx_cdg] * rel[(i + 4) + 1 + (j + 4) * nx + (k + 4) * nx * ny]
                                                                     +  cdg2[(j) * nx_cdg + (grid_id - 1) * ny_cdg2 * nx_cdg] * rel[(i + 4)  + ((j + 4) - 1) * nx + (k + 4) * nx * ny]
                                                                     +  cdg2[(j) * nx_cdg + 2 + (grid_id - 1) * ny_cdg2 * nx_cdg] * rel[(i + 4) + ((j + 4) + 1) * nx + (k + 4) * nx * ny]
                                                                     +  cdg3[(k) * nx_cdg + (grid_id - 1) * ny_cdg3 * nx_cdg] * rel[(i + 4)  + (j + 4) * nx + ((k + 4) - 1) * nx * ny]
                                                                     +  cdg3[(k) * nx_cdg + 2 + (grid_id - 1) * ny_cdg3 * nx_cdg] * rel[(i + 4) + (j + 4) * nx + ((k + 4) + 1) * nx * ny]
                                                                     +  (cdg1[(i) * nx_cdg + (grid_id - 1) * ny_cdg1 * nx_cdg + 1]
                                                                     + cdg2[(j) * nx_cdg + (grid_id - 1) * ny_cdg2 * nx_cdg + 1]
                                                                     + cdg3[(k) * nx_cdg + (grid_id - 1) * ny_cdg3 * nx_cdg + 1]) * rel[ i + 4 + (j + 4) * nx + (k + 4) * nx * ny];



}



