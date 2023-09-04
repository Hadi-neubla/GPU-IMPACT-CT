// !**************************************************************************************************                                    // !* by Hadi Zolfaghari, ARTORG Center, University of Bern, (zolfaghari.haadi@gmail.com)            *      
// !* October 2015 - March 2020                                                                      *            
// !* Modified by Hadi Zolfaghari, University of Cambridge (hz382@damtp.cam.ac.uk)                   *      
// !* April 2020 -                                                                                   *            
// !**************************************************************************************************
extern "C" __global__ void interpolate_kernel2(int dimmension_id, int nx, int ny, int nz, int nxI, int nyI, double* phi, double* cI, double* inter)
//for vel_pre


{
        int i = threadIdx.x + blockIdx.x * blockDim.x;
        int j = threadIdx.y + blockIdx.y * blockDim.y;
        int k = threadIdx.z + blockIdx.z * blockDim.z;



//!!//!!// The serial loop for testing the copy

//          for (int ii= 4; ii <= nx-5; ii++){  
//              for (int jj= 4; jj <= ny-5; jj++){
//                  for (int kk= 4; kk <= nz-5; kk++){

//                      inter[ii + jj * nx + kk * nx * ny] = cI1[0 + (ii-4) * nxI] *phi[ii-2 + jj * nx + kk * nx * ny];                       
//                      for (int ll=1; ll<= 5; ll++){ 
                          
 //                          inter[ii + jj * nx + kk * nx * ny] +=  phi[ii -2 + ll + jj * nx + kk * nx * ny]
//                                                                 * cI1[ll + (ii-4) * nxI];
//                      }
//                  }
//              }
//          }


//!//!//!the shifted indecies! TESTED! 


//    if ((i < 4) || (i > nx-5) || (j < 4) || (j > ny-5) || (k < 4) || (k > nz-5)) return;
//
//    if (dimmension_id == 1) {
//
//        if ((i >= 4) && (i <= nx-5)) {
//            if ((j >= 4) && (j <= ny-5)) {     
//                if ((k >= 4) && (k <= nz-5)) {
//
//
//        int index1 = i + j * (nx) + k * (nx) * (ny);
//
//        inter[index1] = cI[(i-4)*nxI] * phi[index1-2] + phi[ ( k * nx * ny + j * nx + i + 1 - 2) ] * cI[(i-4) * nxI + 1]
//                                                      + phi[ ( k * nx * ny + j * nx + i + 2 - 2) ] * cI[(i-4) * nxI + 2]
//                                                      + phi[ ( k * nx * ny + j * nx + i + 3 - 2) ] * cI[(i-4) * nxI + 3]
//                                                      + phi[ ( k * nx * ny + j * nx + i + 4 - 2) ] * cI[(i-4) * nxI + 4]
//                                                      + phi[ ( k * nx * ny + j * nx + i + 5 - 2) ] * cI[(i-4) * nxI + 5]; 
////        for (int l=1; l<= 5; l++){
//
//
////         inter[index1] += phi[ ( k * nx * ny + j * nx + i + l - 2) ] * cI[(i-4) * nxI + l] ;
//
////                     }
//                  }
//               }
//            }
//
//
//
//}


//!//!//! the warp-oriented implemtation


    if (dimmension_id == 1) {

////        if ((i >= 4) && (i <= nx-5)) {
////           if ((j >= 4) && (j <= ny-5)) {
////                if ((k >= 4) && (k <= nz-5)) {


        int index1 = (i+4) + (j+4) * (nx) + (k+4) * (nx) * (ny);

        inter[index1] = cI[(i+1)*nxI] * phi[index1-3] ;
//+ phi[ ( (k+4) * nx * ny + (j+4) * nx + (i+4) + 1 - 2) ] * cI[(i) * nxI + 1]
//                                                      + phi[ ( (k+4) * nx * ny + (j+4) * nx + (i+4) + 2 - 2) ] * cI[(i) * nxI + 2]
//                                                      + phi[ ( (k+4) * nx * ny + (j+4) * nx + (i+4) + 3 - 2) ] * cI[(i) * nxI + 3]
//                                                      + phi[ ( (k+4) * nx * ny + (j+4) * nx + (i+4) + 4 - 2) ] * cI[(i) * nxI + 4]
//                                                      + phi[ ( (k+4) * nx * ny + (j+4) * nx + (i+4) + 5 - 2) ] * cI[(i) * nxI + 5];
        for (int l=1; l<= 5; l++){
              inter[index1]  = inter[index1] + phi[ ( (k+4) * nx * ny + (j+4) * nx + (i+4) + l - 3) ] * cI[(i+1) * nxI + l];
        }

//         inter[index1] += phi[ ( k * nx * ny + j * nx + i + l - 2) ] * cI[(i-4) * nxI + l] ;

//                     }
////                  }
////               }
////            }



}



//!//!//! the shifted entries-- TESTED!




//   if (dimmension_id ==2) {
//
//        if ((i >= 4) && (i <= nx-5)) {
//            if ((j >= 4) && (j <= ny-5)) {     
//                if ((k >= 4) && (k <= nz-5)) {
//
//
//        int index1 = i + j * (nx) + k * (nx) * (ny);
//
//        inter[index1] = cI[(j-4)*nxI] * phi[index1- 2 * nx] + phi[ ( k * nx * ny + (j + 1 - 2) * nx + i) ] * cI[(j-4) * nxI + 1]
//                                                            + phi[ ( k * nx * ny + (j + 2 - 2) * nx + i) ] * cI[(j-4) * nxI + 2]
//                                                            + phi[ ( k * nx * ny + (j + 3 - 2) * nx + i) ] * cI[(j-4) * nxI + 3]
//                                                            + phi[ ( k * nx * ny + (j + 4 - 2) * nx + i) ] * cI[(j-4) * nxI + 4]
//                                                            + phi[ ( k * nx * ny + (j + 5 - 2) * nx + i) ] * cI[(j-4) * nxI + 5]; 
////        for (int l=1; l<= 5; l++){
//
////         inter[index1] += phi[ ( k * nx * ny + (j + l - 2) * nx + i) ] * cI[(j-4) * nxI + l] ;
//
////                     }
//                  }
//               }
//            }
//}



//!//!//! the warp oriented implementation


   if (dimmension_id ==2) {

////        if ((i >= 4) && (i <= nx-5)) {
////            if ((j >= 4) && (j <= ny-5)) {
////                if ((k >= 4) && (k <= nz-5)) {


        int index1 = (i+4) + (j+4) * (nx) + (k+4) * (nx) * (ny);

        inter[index1] = cI[(j+1)*nxI] * phi[index1- 3 * nx];
// + phi[ ( (k+4) * nx * ny + ((j+4) + 1 - 2) * nx + (i+4)) ] * cI[(j) * nxI + 1]
//                                                            + phi[ ( (k+4) * nx * ny + ((j+4) + 2 - 2) * nx + (i+4)) ] * cI[(j) * nxI + 2]
//                                                            + phi[ ( (k+4) * nx * ny + ((j+4) + 3 - 2) * nx + (i+4)) ] * cI[(j) * nxI + 3]
//                                                            + phi[ ( (k+4) * nx * ny + ((j+4) + 4 - 2) * nx + (i+4)) ] * cI[(j) * nxI + 4]
//                                                            + phi[ ( (k+4) * nx * ny + ((j+4) + 5 - 2) * nx + (i+4)) ] * cI[(j) * nxI + 5];
        for (int l=1; l<= 5; l++){

            inter[index1] = inter[index1] + phi[ ( (k+4) * nx * ny + ((j+4) + l - 3) * nx + (i+4)) ] * cI[(j+1) * nxI + l];
        }
//         inter[index1] += phi[ ( k * nx * ny + (j + l - 2) * nx + i) ] * cI[(j-4) * nxI + l] ;

//                     }
////                  }
////               }
////            }
}







//!//!//! the shifted entries !TESTED




//   if (dimmension_id == 3) { 
//
//        if ((i >= 4) && (i <= nx-5)) {
//            if ((j >= 4) && (j <= ny-5)) {     
//                if ((k >= 4) && (k <= nz-5)) {
//
//
//        int index1 = i + j * (nx) + k * (nx) * (ny);
//
//        inter[index1] = cI[(k-4)*nxI] * phi[index1- 2 * nx * ny] + phi[ ( (k + 1 -2) * nx * ny + j * nx + i ) ] * cI[(k-4) * nxI + 1] 
//                                                                 + phi[ ( (k + 2 -2) * nx * ny + j * nx + i ) ] * cI[(k-4) * nxI + 2] 
//                                                                 + phi[ ( (k + 3 -2) * nx * ny + j * nx + i ) ] * cI[(k-4) * nxI + 3] 
//                                                                 + phi[ ( (k + 4 -2) * nx * ny + j * nx + i ) ] * cI[(k-4) * nxI + 4] 
//                                                                 + phi[ ( (k + 5 -2) * nx * ny + j * nx + i ) ] * cI[(k-4) * nxI + 5] ; 
////        for (int l=1; l<= 5; l++){
//
//
////         inter[index1] += phi[ ( (k + l -2) * nx * ny + j * nx + i ) ] * cI[(k-4) * nxI + l] ;
//
//
////                     }
//                  }
//               }
//            }
//         
//}




//!//!//! the warp oriented implementation

   if (dimmension_id == 3) {

////        if ((i >= 4) && (i <= nx-5)) {
////            if ((j >= 4) && (j <= ny-5)) {
////                if ((k >= 4) && (k <= nz-5)) {

       i += 4;
       j += 4;
       k += 4;


        int index1 = i + j * (nx) + k * (nx) * (ny);



        inter[index1] = cI[(k+1-4)*nxI] * phi[index1- 3 * nx * ny];
// + phi[ ( (k + 1 -2) * nx * ny + j * nx + i ) ] * cI[(k-4) * nxI + 1]
//                                                                 + phi[ ( (k + 2 -2) * nx * ny + j * nx + i ) ] * cI[(k-4) * nxI + 2]
//                                                                 + phi[ ( (k + 3 -2) * nx * ny + j * nx + i ) ] * cI[(k-4) * nxI + 3]
//                                                                 + phi[ ( (k + 4 -2) * nx * ny + j * nx + i ) ] * cI[(k-4) * nxI + 4]
//                                                                 + phi[ ( (k + 5 -2) * nx * ny + j * nx + i ) ] * cI[(k-4) * nxI + 5] ;
        for (int l=1; l<= 5; l++){
         inter[index1] = inter[index1] + phi[ ( (k + l -3) * nx * ny + j * nx + i ) ] * cI[(k+1-4) * nxI + l]; 
        }

//         inter[index1] += phi[ ( (k + l -2) * nx * ny + j * nx + i ) ] * cI[(k-4) * nxI + l] ;


//                     }
////                  }
////               }
////            }

}








}



