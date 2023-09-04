// !**************************************************************************************************                                           
// !* by Hadi Zolfaghari, ARTORG Center, University of Bern, (zolfaghari.haadi@gmail.com)            *      
// !* October 2015 - March 2020                                                                      *            
// !**************************************************************************************************
extern "C" __global__ void gradient_kernel(int dimmension_id, int nx, int ny, int nz, int nxg, int nyg, double* pre, double* cGp, double* grad)



{
        int i = threadIdx.x + blockIdx.x * blockDim.x;
        int j = threadIdx.y + blockIdx.y * blockDim.y;
        int k = threadIdx.z + blockIdx.z * blockDim.z;



//!!//!!// The serial loop for testing the copy

//          for (int ii= 4; ii <= nx-5; ii++){  
//              for (int jj= 4; jj <= ny-5; jj++){
//                  for (int kk= 4; kk <= nz-5; kk++){

//                      grad[ii + jj * nx + kk * nx * ny] = cGp1[0 + (ii-4) * nxg] *pre[ii-2 + jj * nx + kk * nx * ny];                       
//                      for (int ll=1; ll<= 5; ll++){ 
                          
 //                          grad[ii + jj * nx + kk * nx * ny] +=  pre[ii -2 + ll + jj * nx + kk * nx * ny]
//                                                                 * cGp1[ll + (ii-4) * nxg];
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
//        grad[index1] = cGp[(i-4)*nxg] * pre[index1-2] + pre[ ( k * nx * ny + j * nx + i + 1 - 2) ] * cGp[(i-4) * nxg + 1]
//                                                      + pre[ ( k * nx * ny + j * nx + i + 2 - 2) ] * cGp[(i-4) * nxg + 2]
//                                                      + pre[ ( k * nx * ny + j * nx + i + 3 - 2) ] * cGp[(i-4) * nxg + 3]
//                                                      + pre[ ( k * nx * ny + j * nx + i + 4 - 2) ] * cGp[(i-4) * nxg + 4]
//                                                      + pre[ ( k * nx * ny + j * nx + i + 5 - 2) ] * cGp[(i-4) * nxg + 5]; 
////        for (int l=1; l<= 5; l++){
//
//
////         grad[index1] += pre[ ( k * nx * ny + j * nx + i + l - 2) ] * cGp[(i-4) * nxg + l] ;
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

        grad[index1] = cGp[(i)*nxg] * pre[index1-2] ;
//+ pre[ ( (k+4) * nx * ny + (j+4) * nx + (i+4) + 1 - 2) ] * cGp[(i) * nxg + 1]
//                                                      + pre[ ( (k+4) * nx * ny + (j+4) * nx + (i+4) + 2 - 2) ] * cGp[(i) * nxg + 2]
//                                                      + pre[ ( (k+4) * nx * ny + (j+4) * nx + (i+4) + 3 - 2) ] * cGp[(i) * nxg + 3]
//                                                      + pre[ ( (k+4) * nx * ny + (j+4) * nx + (i+4) + 4 - 2) ] * cGp[(i) * nxg + 4]
//                                                      + pre[ ( (k+4) * nx * ny + (j+4) * nx + (i+4) + 5 - 2) ] * cGp[(i) * nxg + 5];
        for (int l=1; l<= 5; l++){
              grad[index1]  = grad[index1] + pre[ ( (k+4) * nx * ny + (j+4) * nx + (i+4) + l - 2) ] * cGp[(i) * nxg + l];
        }

//         grad[index1] += pre[ ( k * nx * ny + j * nx + i + l - 2) ] * cGp[(i-4) * nxg + l] ;

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
//        grad[index1] = cGp[(j-4)*nxg] * pre[index1- 2 * nx] + pre[ ( k * nx * ny + (j + 1 - 2) * nx + i) ] * cGp[(j-4) * nxg + 1]
//                                                            + pre[ ( k * nx * ny + (j + 2 - 2) * nx + i) ] * cGp[(j-4) * nxg + 2]
//                                                            + pre[ ( k * nx * ny + (j + 3 - 2) * nx + i) ] * cGp[(j-4) * nxg + 3]
//                                                            + pre[ ( k * nx * ny + (j + 4 - 2) * nx + i) ] * cGp[(j-4) * nxg + 4]
//                                                            + pre[ ( k * nx * ny + (j + 5 - 2) * nx + i) ] * cGp[(j-4) * nxg + 5]; 
////        for (int l=1; l<= 5; l++){
//
////         grad[index1] += pre[ ( k * nx * ny + (j + l - 2) * nx + i) ] * cGp[(j-4) * nxg + l] ;
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

        grad[index1] = cGp[(j)*nxg] * pre[index1- 2 * nx];
// + pre[ ( (k+4) * nx * ny + ((j+4) + 1 - 2) * nx + (i+4)) ] * cGp[(j) * nxg + 1]
//                                                            + pre[ ( (k+4) * nx * ny + ((j+4) + 2 - 2) * nx + (i+4)) ] * cGp[(j) * nxg + 2]
//                                                            + pre[ ( (k+4) * nx * ny + ((j+4) + 3 - 2) * nx + (i+4)) ] * cGp[(j) * nxg + 3]
//                                                            + pre[ ( (k+4) * nx * ny + ((j+4) + 4 - 2) * nx + (i+4)) ] * cGp[(j) * nxg + 4]
//                                                            + pre[ ( (k+4) * nx * ny + ((j+4) + 5 - 2) * nx + (i+4)) ] * cGp[(j) * nxg + 5];
        for (int l=1; l<= 5; l++){

            grad[index1] = grad[index1] + pre[ ( (k+4) * nx * ny + ((j+4) + l - 2) * nx + (i+4)) ] * cGp[(j) * nxg + l];
        }
//         grad[index1] += pre[ ( k * nx * ny + (j + l - 2) * nx + i) ] * cGp[(j-4) * nxg + l] ;

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
//        grad[index1] = cGp[(k-4)*nxg] * pre[index1- 2 * nx * ny] + pre[ ( (k + 1 -2) * nx * ny + j * nx + i ) ] * cGp[(k-4) * nxg + 1] 
//                                                                 + pre[ ( (k + 2 -2) * nx * ny + j * nx + i ) ] * cGp[(k-4) * nxg + 2] 
//                                                                 + pre[ ( (k + 3 -2) * nx * ny + j * nx + i ) ] * cGp[(k-4) * nxg + 3] 
//                                                                 + pre[ ( (k + 4 -2) * nx * ny + j * nx + i ) ] * cGp[(k-4) * nxg + 4] 
//                                                                 + pre[ ( (k + 5 -2) * nx * ny + j * nx + i ) ] * cGp[(k-4) * nxg + 5] ; 
////        for (int l=1; l<= 5; l++){
//
//
////         grad[index1] += pre[ ( (k + l -2) * nx * ny + j * nx + i ) ] * cGp[(k-4) * nxg + l] ;
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



        grad[index1] = cGp[(k-4)*nxg] * pre[index1- 2 * nx * ny];
// + pre[ ( (k + 1 -2) * nx * ny + j * nx + i ) ] * cGp[(k-4) * nxg + 1]
//                                                                 + pre[ ( (k + 2 -2) * nx * ny + j * nx + i ) ] * cGp[(k-4) * nxg + 2]
//                                                                 + pre[ ( (k + 3 -2) * nx * ny + j * nx + i ) ] * cGp[(k-4) * nxg + 3]
//                                                                 + pre[ ( (k + 4 -2) * nx * ny + j * nx + i ) ] * cGp[(k-4) * nxg + 4]
//                                                                 + pre[ ( (k + 5 -2) * nx * ny + j * nx + i ) ] * cGp[(k-4) * nxg + 5] ;
        for (int l=1; l<= 5; l++){
         grad[index1] = grad[index1] + pre[ ( (k + l -2) * nx * ny + j * nx + i ) ] * cGp[(k-4) * nxg + l]; 
        }

//         grad[index1] += pre[ ( (k + l -2) * nx * ny + j * nx + i ) ] * cGp[(k-4) * nxg + l] ;


//                     }
////                  }
////               }
////            }

}








}



