// !**************************************************************************************************                                           
// !* by Hadi Zolfaghari, ARTORG Center, University of Bern, (zolfaghari.haadi@gmail.com)            *      
// !* October 2015 - March 2020                                                                      *            
// !**************************************************************************************************
extern "C" __global__ void stencil_kernel(int dimmension_id, int lower_ind, int upper_ind, int nx, int ny, int nz, int nxS, int nyS, double* phi, double* cS, double* output)



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
//    if (dimmension_id == 1) {
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


//!//!//! the warp-oriented implemtation


    if (dimmension_id == 1) {

////        if ((i >= 4) && (i <= nx-5)) {
////           if ((j >= 4) && (j <= ny-5)) {
////                if ((k >= 4) && (k <= nz-5)) {


        int index1 = (i+4) + (j+4) * (nx) + (k+4) * (nx) * (ny);

        output[index1] = cS[(-lower_ind -2 + i) * nxS] * phi[index1 + lower_ind] ;
//+ phi[ ( (k+4) * nx * ny + (j+4) * nx + (i+4) + 1 - 2) ] * cS[(i) * nxS + 1]
//                                                      + phi[ ( (k+4) * nx * ny + (j+4) * nx + (i+4) + 2 - 2) ] * cS[(i) * nxS + 2]
//                                                      + phi[ ( (k+4) * nx * ny + (j+4) * nx + (i+4) + 3 - 2) ] * cS[(i) * nxS + 3]
//                                                      + phi[ ( (k+4) * nx * ny + (j+4) * nx + (i+4) + 4 - 2) ] * cS[(i) * nxS + 4]
//                                                      + phi[ ( (k+4) * nx * ny + (j+4) * nx + (i+4) + 5 - 2) ] * cS[(i) * nxS + 5];
        for (int l=1; l<= (upper_ind - lower_ind); l++){
              output[index1]  = output[index1] + phi[ ( (k+4) * nx * ny + (j+4) * nx + (i+4) + l + lower_ind) ] * cS[(i - lower_ind + 2) * nxS + l];
        }

//         output[index1] += phi[ ( k * nx * ny + j * nx + i + l - 2) ] * cS[(i-4) * nxS + l] ;

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
//        output[index1] = cS[(j-4)*nxS] * phi[index1- 2 * nx] + phi[ ( k * nx * ny + (j + 1 - 2) * nx + i) ] * cS[(j-4) * nxS + 1]
//                                                            + phi[ ( k * nx * ny + (j + 2 - 2) * nx + i) ] * cS[(j-4) * nxS + 2]
//                                                            + phi[ ( k * nx * ny + (j + 3 - 2) * nx + i) ] * cS[(j-4) * nxS + 3]
//                                                            + phi[ ( k * nx * ny + (j + 4 - 2) * nx + i) ] * cS[(j-4) * nxS + 4]
//                                                            + phi[ ( k * nx * ny + (j + 5 - 2) * nx + i) ] * cS[(j-4) * nxS + 5]; 
////        for (int l=1; l<= 5; l++){
//
////         output[index1] += phi[ ( k * nx * ny + (j + l - 2) * nx + i) ] * cS[(j-4) * nxS + l] ;
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

        output[index1] = cS[(j - lower_ind + 2) * nxS] * phi[index1 + lower_ind * nx];
// + phi[ ( (k+4) * nx * ny + ((j+4) + 1 - 2) * nx + (i+4)) ] * cS[(j) * nxS + 1]
//                                                            + phi[ ( (k+4) * nx * ny + ((j+4) + 2 - 2) * nx + (i+4)) ] * cS[(j) * nxS + 2]
//                                                            + phi[ ( (k+4) * nx * ny + ((j+4) + 3 - 2) * nx + (i+4)) ] * cS[(j) * nxS + 3]
//                                                            + phi[ ( (k+4) * nx * ny + ((j+4) + 4 - 2) * nx + (i+4)) ] * cS[(j) * nxS + 4]
//                                                            + phi[ ( (k+4) * nx * ny + ((j+4) + 5 - 2) * nx + (i+4)) ] * cS[(j) * nxS + 5];
        for (int l=1; l<= (upper_ind - lower_ind); l++){

            output[index1] = output[index1] + phi[ ( (k+4) * nx * ny + ((j+4) + l + lower_ind) * nx + (i+4)) ] * cS[(j - lower_ind + 2) * nxS + l];
        }
//         output[index1] += phi[ ( k * nx * ny + (j + l - 2) * nx + i) ] * cS[(j-4) * nxS + l] ;

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
//        output[index1] = cS[(k-4)*nxS] * phi[index1- 2 * nx * ny] + phi[ ( (k + 1 -2) * nx * ny + j * nx + i ) ] * cS[(k-4) * nxS + 1] 
//                                                                 + phi[ ( (k + 2 -2) * nx * ny + j * nx + i ) ] * cS[(k-4) * nxS + 2] 
//                                                                 + phi[ ( (k + 3 -2) * nx * ny + j * nx + i ) ] * cS[(k-4) * nxS + 3] 
//                                                                 + phi[ ( (k + 4 -2) * nx * ny + j * nx + i ) ] * cS[(k-4) * nxS + 4] 
//                                                                 + phi[ ( (k + 5 -2) * nx * ny + j * nx + i ) ] * cS[(k-4) * nxS + 5] ; 
////        for (int l=1; l<= 5; l++){
//
//
////         output[index1] += phi[ ( (k + l -2) * nx * ny + j * nx + i ) ] * cS[(k-4) * nxS + l] ;
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



        output[index1] = cS[(k - (lower_ind + 2) - 4)*nxS] * phi[index1 + lower_ind * nx * ny];
// + phi[ ( (k + 1 -2) * nx * ny + j * nx + i ) ] * cS[(k-4) * nxS + 1]
//                                                                 + phi[ ( (k + 2 -2) * nx * ny + j * nx + i ) ] * cS[(k-4) * nxS + 2]
//                                                                 + phi[ ( (k + 3 -2) * nx * ny + j * nx + i ) ] * cS[(k-4) * nxS + 3]
//                                                                 + phi[ ( (k + 4 -2) * nx * ny + j * nx + i ) ] * cS[(k-4) * nxS + 4]
//                                                                 + phi[ ( (k + 5 -2) * nx * ny + j * nx + i ) ] * cS[(k-4) * nxS + 5] ;
        for (int l=1; l<= (upper_ind - lower_ind); l++){
         output[index1] = output[index1] + phi[ ( (k + l + lower_ind) * nx * ny + j * nx + i ) ] * cS[(k -(lower_ind + 2) -4) * nxS + l]; 
        }

//         output[index1] += phi[ ( (k + l -2) * nx * ny + j * nx + i ) ] * cS[(k-4) * nxS + l] ;


//                     }
////                  }
////               }
////            }

}








}



