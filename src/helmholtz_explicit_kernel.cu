// !**************************************************************************************************                                    // !* by Hadi Zolfaghari, ARTORG Center, University of Bern, (zolfaghari.haadi@gmail.com)            *      
// !* October 2015 - March 2020                                                                      *            
// !* Modified by Hadi Zolfaghari, University of Cambridge (hz382@damtp.cam.ac.uk)                   *      
// !* April 2020 -                                                                                   *            
// !**************************************************************************************************
//#include <stdio.h>
extern "C" __global__ void helmholtz_kernel(int component_id, int nx, int ny, int nz, int nxS1, int nxS2,int nxS3, double multL, double* vel, double* nl, double* cS1, double* cS2, double* cS3)

{
        int i = threadIdx.x + blockIdx.x * blockDim.x;
        int j = threadIdx.y + blockIdx.y * blockDim.y;
        int k = threadIdx.z + blockIdx.z * blockDim.z;



//!!//!!// The serial loop for testing the copy

//          for (int ii= 4; ii <= nx-5; ii++){  
//              for (int jj= 4; jj <= ny-5; jj++){
//                  for (int kk= 4; kk <= nz-5; kk++){

//                      output[ii + jj * nx + kk * nx * ny] = cS1[0 + (ii-4) * nxS] *phi[ii-2 + jj * nx + kk * nx * ny];                       
//                      for (int ll=1; ll<= 5; ll++){ 
                          
 //                          output[ii + jj * nx + kk * nx * ny] +=  phi[ii -2 + ll + jj * nx + kk * nx * ny]
//                                                                 * cS1[ll + (ii-4) * nxS];
//                      }
//                  }
//              }
//          }


//!//!//!the shifted indecies! TESTED! 


//    if ((i < 4) || (i > nx-5) || (j < 4) || (j > ny-5) || (k < 4) || (k > nz-5)) return;
//
//    if (component_id == 1) {
//
//        if ((i >= 4) && (i <= nx-5)) {
//            if ((j >= 4) && (j <= ny-5)) {     
//                if ((k >= 4) && (k <= nz-5)) {
//
//
//        int index1 = i + j * (nx) + k * (nx) * (ny);
//
//        output[index1] = cS[(i-4)*nxS] * phi[index1-2] + phi[ ( k * nx * ny + j * nx + i + 1 - 2) ] * cS[(i-4) * nxS + 1]
//                                                      + phi[ ( k * nx * ny + j * nx + i + 2 - 2) ] * cS[(i-4) * nxS + 2]
//                                                      + phi[ ( k * nx * ny + j * nx + i + 3 - 2) ] * cS[(i-4) * nxS + 3]
//                                                      + phi[ ( k * nx * ny + j * nx + i + 4 - 2) ] * cS[(i-4) * nxS + 4]
//                                                      + phi[ ( k * nx * ny + j * nx + i + 5 - 2) ] * cS[(i-4) * nxS + 5]; 
////        for (int l=1; l<= 5; l++){
//
//
////         output[index1] += phi[ ( k * nx * ny + j * nx + i + l - 2) ] * cS[(i-4) * nxS + l] ;
//
////                     }
//                  }
//               }
//            }
//
//
//
//}


//if ((i == 16) && (j == 16) && (k == 16)) {
//   int index1 = (i+4) + (j+4) * (nx) + (k+4) * (nx) * (ny);
//   printf("vel=%f", vel[index1-3]);
//};

//if (component_id == 2) {
//   printf("Value of dims=%d,%d,%d", nxS1,nxS2,nxS3);
//};



//!//!//! the warp-oriented implemtation

    if (component_id >= 1) {

////        if ((i >= 4) && (i <= nx-5)) {
////           if ((j >= 4) && (j <= ny-5)) {
////                if ((k >= 4) && (k <= nz-5)) {

        int index1 = (i+4) + (j+4) * (nx) + (k+4) * (nx) * (ny);

        double dd1 = cS1[(1 + i) * nxS1] * vel[index1 - 3] ;
////                                                      + phi[ ( (k+4) * nx * ny + (j+4) * nx + (i+4) + 1 - 2) ] * cS[(i) * nxS + 1]
////                                                      + phi[ ( (k+4) * nx * ny + (j+4) * nx + (i+4) + 2 - 2) ] * cS[(i) * nxS + 2]
////                                                      + phi[ ( (k+4) * nx * ny + (j+4) * nx + (i+4) + 3 - 2) ] * cS[(i) * nxS + 3]
////                                                      + phi[ ( (k+4) * nx * ny + (j+4) * nx + (i+4) + 4 - 2) ] * cS[(i) * nxS + 4]
////                                                      + phi[ ( (k+4) * nx * ny + (j+4) * nx + (i+4) + 5 - 2) ] * cS[(i) * nxS + 5];

////-- looped implementation 
        for (int l=1; l<= 6; l++){
             dd1  = dd1 + vel[ ( (k+4) * nx * ny + (j+4) * nx + (i+4) + l - 3) ] * cS1[(i + 1) * nxS1 + l];
        }

        for (int l=0; l<= 6; l++){
            dd1 = dd1 + vel[ ( (k+4) * nx * ny + ((j+4) + l -3) * nx + (i + 4)) ] * cS2[(j + 1) * nxS2 + l];
        }

         for (int l=0; l<= 6; l++){
             dd1 = dd1 + vel[ ( (k  + 4 + l -3) * nx * ny + (j+4) * nx + (i+4) ) ] * cS3[(k + 1) * nxS3 + l]; 
        }

        nl[index1] = nl[index1] - multL*dd1;

////         output[index1] += phi[ ( k * nx * ny + j * nx + i + l - 2) ] * cS[(i-4) * nxS + l] ;

////                     }
////                  }
////               }
////            }


}


}



