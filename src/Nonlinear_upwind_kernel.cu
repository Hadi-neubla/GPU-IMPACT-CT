// !**************************************************************************************************                                    // !* by Hadi Zolfaghari, ARTORG Center, University of Bern, (zolfaghari.haadi@gmail.com)            *      
// !* October 2015 - March 2020                                                                      *            
// !* Modified by Hadi Zolfaghari, University of Cambridge (hz382@damtp.cam.ac.uk)                   *      
// !* April 2020 -                                                                                   *            
// !**************************************************************************************************
extern "C" __global__ void nonlinear_upwind_kernel(int component_id, int nx, int ny, int nz, int nxS, int nyS, double* pp, double* vel, double* nl, double* cS, double* cS2)



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


//!//!//! the warp-oriented implemtation

     ////--- u.du/dx
    if (component_id == 1) {

////        if ((i >= 4) && (i <= nx-5)) {
////           if ((j >= 4) && (j <= ny-5)) {
////                if ((k >= 4) && (k <= nz-5)) {

        int index1 = (i+4) + (j+4) * (nx) + (k+4) * (nx) * (ny);

        if (pp[index1] >= 0.) {

        double dd1 = cS[(1 + i) * nxS] * vel[index1 - 3] ;
////                                                      + phi[ ( (k+4) * nx * ny + (j+4) * nx + (i+4) + 1 - 2) ] * cS[(i) * nxS + 1]
////                                                      + phi[ ( (k+4) * nx * ny + (j+4) * nx + (i+4) + 2 - 2) ] * cS[(i) * nxS + 2]
////                                                      + phi[ ( (k+4) * nx * ny + (j+4) * nx + (i+4) + 3 - 2) ] * cS[(i) * nxS + 3]
////                                                      + phi[ ( (k+4) * nx * ny + (j+4) * nx + (i+4) + 4 - 2) ] * cS[(i) * nxS + 4]
////                                                      + phi[ ( (k+4) * nx * ny + (j+4) * nx + (i+4) + 5 - 2) ] * cS[(i) * nxS + 5];

////-- looped implementation 
        for (int l=1; l<= 6; l++){
             dd1  = dd1 + vel[ ( (k+4) * nx * ny + (j+4) * nx + (i+4) + l - 3) ] * cS[(i + 1) * nxS + l];
        }

             nl[index1] = dd1*pp[index1];
                      }
        else {

        double dd1 = cS2[(1 + i) * nxS] * vel[index1 - 3] ;
////                                                      + phi[ ( (k+4) * nx * ny + (j+4) * nx + (i+4) + 1 - 2) ] * cS[(i) * nxS + 1]
////                                                      + phi[ ( (k+4) * nx * ny + (j+4) * nx + (i+4) + 2 - 2) ] * cS[(i) * nxS + 2]
////                                                      + phi[ ( (k+4) * nx * ny + (j+4) * nx + (i+4) + 3 - 2) ] * cS[(i) * nxS + 3]
////                                                      + phi[ ( (k+4) * nx * ny + (j+4) * nx + (i+4) + 4 - 2) ] * cS[(i) * nxS + 4]
////                                                      + phi[ ( (k+4) * nx * ny + (j+4) * nx + (i+4) + 5 - 2) ] * cS[(i) * nxS + 5];

////-- looped implementation
        for (int l=1; l<= 6; l++){
             dd1  = dd1 + vel[ ( (k+4) * nx * ny + (j+4) * nx + (i+4) + l - 3) ] * cS2[(i + 1) * nxS + l];
        }

             nl[index1] = dd1*pp[index1];
                      }

       

////         output[index1] += phi[ ( k * nx * ny + j * nx + i + l - 2) ] * cS[(i-4) * nxS + l] ;

////                     }
////                  }
////               }
////            }



}



//!//!//! the shifted entries-- TESTED!




//   if (component_id ==2) {
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

   ////--v.du/dy 
   if (component_id == 2) {

////        if ((i >= 4) && (i <= nx-5)) {
////            if ((j >= 4) && (j <= ny-5)) {
////                if ((k >= 4) && (k <= nz-5)) {


        int index1 = (i+4) + (j+4) * (nx) + (k+4) * (nx) * (ny);

    if (pp[index1] >= 0.){
        double dd1 = cS[(j + 1) * nxS] * vel[index1 -3 * nx];
////                                                            + phi[ ( (k+4) * nx * ny + ((j+4) + 1 - 2) * nx + (i+4)) ] * cS[(j) * nxS + 1]
////                                                            + phi[ ( (k+4) * nx * ny + ((j+4) + 2 - 2) * nx + (i+4)) ] * cS[(j) * nxS + 2]
////                                                            + phi[ ( (k+4) * nx * ny + ((j+4) + 3 - 2) * nx + (i+4)) ] * cS[(j) * nxS + 3]
////                                                            + phi[ ( (k+4) * nx * ny + ((j+4) + 4 - 2) * nx + (i+4)) ] * cS[(j) * nxS + 4]
////                                                            + phi[ ( (k+4) * nx * ny + ((j+4) + 5 - 2) * nx + (i+4)) ] * cS[(j) * nxS + 5];
        for (int l=1; l<= 6; l++){

            dd1 = dd1 + vel[ ( (k+4) * nx * ny + ((j+4) + l -3) * nx + (i + 4)) ] * cS[(j + 1) * nxS + l];
        }

            nl[index1] = nl[index1] + dd1 * pp[index1];

                          }

    else {
        double dd1 = cS2[(j + 1) * nxS] * vel[index1 -3 * nx];
////                                                            + phi[ ( (k+4) * nx * ny + ((j+4) + 1 - 2) * nx + (i+4)) ] * cS[(j) * nxS + 1]
////                                                            + phi[ ( (k+4) * nx * ny + ((j+4) + 2 - 2) * nx + (i+4)) ] * cS[(j) * nxS + 2]
////                                                            + phi[ ( (k+4) * nx * ny + ((j+4) + 3 - 2) * nx + (i+4)) ] * cS[(j) * nxS + 3]
////                                                            + phi[ ( (k+4) * nx * ny + ((j+4) + 4 - 2) * nx + (i+4)) ] * cS[(j) * nxS + 4]
////                                                            + phi[ ( (k+4) * nx * ny + ((j+4) + 5 - 2) * nx + (i+4)) ] * cS[(j) * nxS + 5];
        for (int l=1; l<= 6; l++){

            dd1 = dd1 + vel[ ( (k+4) * nx * ny + ((j+4) + l -3) * nx + (i + 4)) ] * cS2[(j + 1) * nxS + l];
        }

            nl[index1] = nl[index1] + dd1 * pp[index1];

                          }


 
////         output[index1] += phi[ ( k * nx * ny + (j + l - 2) * nx + i) ] * cS[(j-4) * nxS + l] ;

////                    }
////                  }
////               }
////            }
}







//!//!//! the shifted entries !TESTED




//   if (component_id == 3) { 
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

   if (component_id == 3) {

////        if ((i >= 4) && (i <= nx-5)) {
////            if ((j >= 4) && (j <= ny-5)) {
////                if ((k >= 4) && (k <= nz-5)) {

////       i += 4;
////       j += 4;
////       k += 4;


        int index1 = i + 4 + (j + 4) * (nx) + (k + 4) * (nx) * (ny);


    if (pp[index1] >= 0.0) {
        double dd1 = cS[(k +1)*nxS] * vel[index1 -3 * nx * ny];
////                                                                 + phi[ ( (k + 1 -2) * nx * ny + j * nx + i ) ] * cS[(k-4) * nxS + 1]
////                                                                 + phi[ ( (k + 2 -2) * nx * ny + j * nx + i ) ] * cS[(k+1) * nxS + 2]
////                                                                 + phi[ ( (k + 3 -2) * nx * ny + j * nx + i ) ] * cS[(k+1) * nxS + 3]
////                                                                 + phi[ ( (k + 4 -2) * nx * ny + j * nx + i ) ] * cS[(k+1) * nxS + 4]
////                                                                 + phi[ ( (k + 5 -2) * nx * ny + j * nx + i ) ] * cS[(k+1) * nxS + 5] ;
        for (int l=1; l<= 6; l++){
             dd1 = dd1 + vel[ ( (k  + 4 + l -3) * nx * ny + (j+4) * nx + (i+4) ) ] * cS[(k +1) * nxS + l]; 
        }

             nl[index1] = nl[index1] + dd1 * pp[index1];

                           }
    else {
        double dd1 = cS2[(k +1)*nxS] * vel[index1 -3 * nx * ny];
////                                                                 + phi[ ( (k + 1 -2) * nx * ny + j * nx + i ) ] * cS[(k-4) * nxS + 1]
////                                                                 + phi[ ( (k + 2 -2) * nx * ny + j * nx + i ) ] * cS[(k+1) * nxS + 2]
////                                                                 + phi[ ( (k + 3 -2) * nx * ny + j * nx + i ) ] * cS[(k+1) * nxS + 3]
////                                                                 + phi[ ( (k + 4 -2) * nx * ny + j * nx + i ) ] * cS[(k+1) * nxS + 4]
////                                                                 + phi[ ( (k + 5 -2) * nx * ny + j * nx + i ) ] * cS[(k+1) * nxS + 5] ;
        for (int l=1; l<= 6; l++){
             dd1 = dd1 + vel[ ( (k  + 4 + l -3) * nx * ny + (j+4) * nx + (i+4) ) ] * cS2[(k +1) * nxS + l];
        }

             nl[index1] = nl[index1] + dd1 * pp[index1];

                           }
//         output[index1] += phi[ ( (k + l -2) * nx * ny + j * nx + i ) ] * cS[(k-4) * nxS + l] ;

//                     }
////                  }
////               }
////            }

}






   ////--- u.dv/dx
    if (component_id == 4) {

////        if ((i >= 4) && (i <= nx-5)) {
////           if ((j >= 4) && (j <= ny-5)) {
////                if ((k >= 4) && (k <= nz-5)) {

        int index1 = (i+4) + (j+4) * (nx) + (k+4) * (nx) * (ny);


     if (pp[index1] >= 0.0) {

        double dd1 = cS[(1 + i) * nxS] * vel[index1 - 3] ;
////                                                      + phi[ ( (k+4) * nx * ny + (j+4) * nx + (i+4) + 1 - 2) ] * cS[(i) * nxS + 1]
////                                                      + phi[ ( (k+4) * nx * ny + (j+4) * nx + (i+4) + 2 - 2) ] * cS[(i) * nxS + 2]
////                                                      + phi[ ( (k+4) * nx * ny + (j+4) * nx + (i+4) + 3 - 2) ] * cS[(i) * nxS + 3]
////                                                      + phi[ ( (k+4) * nx * ny + (j+4) * nx + (i+4) + 4 - 2) ] * cS[(i) * nxS + 4]
////                                                      + phi[ ( (k+4) * nx * ny + (j+4) * nx + (i+4) + 5 - 2) ] * cS[(i) * nxS + 5];

////-- looped implementation 
        for (int l=1; l<= 6; l++){
             dd1  = dd1 + vel[ ( (k+4) * nx * ny + (j+4) * nx + (i+4) + l - 3) ] * cS[(i + 1) * nxS + l];
        }
             nl[index1] = dd1*pp[index1];

                          }


     else{

        double dd1 = cS2[(1 + i) * nxS] * vel[index1 - 3] ;
////                                                      + phi[ ( (k+4) * nx * ny + (j+4) * nx + (i+4) + 1 - 2) ] * cS[(i) * nxS + 1]
////                                                      + phi[ ( (k+4) * nx * ny + (j+4) * nx + (i+4) + 2 - 2) ] * cS[(i) * nxS + 2]
////                                                      + phi[ ( (k+4) * nx * ny + (j+4) * nx + (i+4) + 3 - 2) ] * cS[(i) * nxS + 3]
////                                                      + phi[ ( (k+4) * nx * ny + (j+4) * nx + (i+4) + 4 - 2) ] * cS[(i) * nxS + 4]
////                                                      + phi[ ( (k+4) * nx * ny + (j+4) * nx + (i+4) + 5 - 2) ] * cS[(i) * nxS + 5];

////-- looped implementation
        for (int l=1; l<= 6; l++){
             dd1  = dd1 + vel[ ( (k+4) * nx * ny + (j+4) * nx + (i+4) + l - 3) ] * cS2[(i + 1) * nxS + l];
        }
             nl[index1] = dd1*pp[index1];

                          }
////         output[index1] += phi[ ( k * nx * ny + j * nx + i + l - 2) ] * cS[(i-4) * nxS + l] ;

////                     }
////                  }
////               }
////            }



}



//!//!//! the shifted entries-- TESTED!




//   if (component_id ==2) {
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

   ////--v.du/dy 
   if (component_id == 5) {

////        if ((i >= 4) && (i <= nx-5)) {
////            if ((j >= 4) && (j <= ny-5)) {
////                if ((k >= 4) && (k <= nz-5)) {


        int index1 = (i+4) + (j+4) * (nx) + (k+4) * (nx) * (ny);

   if (pp[index1] >= 0.0){

        double dd1 = cS[(j + 1) * nxS] * vel[index1 -3 * nx];
////                                                            + phi[ ( (k+4) * nx * ny + ((j+4) + 1 - 2) * nx + (i+4)) ] * cS[(j) * nxS + 1]
////                                                            + phi[ ( (k+4) * nx * ny + ((j+4) + 2 - 2) * nx + (i+4)) ] * cS[(j) * nxS + 2]
////                                                            + phi[ ( (k+4) * nx * ny + ((j+4) + 3 - 2) * nx + (i+4)) ] * cS[(j) * nxS + 3]
////                                                            + phi[ ( (k+4) * nx * ny + ((j+4) + 4 - 2) * nx + (i+4)) ] * cS[(j) * nxS + 4]
////                                                            + phi[ ( (k+4) * nx * ny + ((j+4) + 5 - 2) * nx + (i+4)) ] * cS[(j) * nxS + 5];
        for (int l=1; l<= 6; l++){

            dd1 = dd1 + vel[ ( (k+4) * nx * ny + ((j+4) + l -3) * nx + (i + 4)) ] * cS[(j + 1) * nxS + l];
        }

            nl[index1] = nl[index1] + dd1 * pp[index1];

                        }

   else {

        double dd1 = cS2[(j + 1) * nxS] * vel[index1 -3 * nx];
////                                                            + phi[ ( (k+4) * nx * ny + ((j+4) + 1 - 2) * nx + (i+4)) ] * cS[(j) * nxS + 1]
////                                                            + phi[ ( (k+4) * nx * ny + ((j+4) + 2 - 2) * nx + (i+4)) ] * cS[(j) * nxS + 2]
////                                                            + phi[ ( (k+4) * nx * ny + ((j+4) + 3 - 2) * nx + (i+4)) ] * cS[(j) * nxS + 3]
////                                                            + phi[ ( (k+4) * nx * ny + ((j+4) + 4 - 2) * nx + (i+4)) ] * cS[(j) * nxS + 4]
////                                                            + phi[ ( (k+4) * nx * ny + ((j+4) + 5 - 2) * nx + (i+4)) ] * cS[(j) * nxS + 5];
        for (int l=1; l<= 6; l++){

            dd1 = dd1 + vel[ ( (k+4) * nx * ny + ((j+4) + l -3) * nx + (i + 4)) ] * cS2[(j + 1) * nxS + l];
        }

            nl[index1] = nl[index1] + dd1 * pp[index1];

                        }
 
////         output[index1] += phi[ ( k * nx * ny + (j + l - 2) * nx + i) ] * cS[(j-4) * nxS + l] ;

////                    }
////                  }
////               }
////            }
}







//!//!//! the shifted entries !TESTED




//   if (component_id == 3) { 
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

   if (component_id == 6) {

////        if ((i >= 4) && (i <= nx-5)) {
////            if ((j >= 4) && (j <= ny-5)) {
////                if ((k >= 4) && (k <= nz-5)) {

////       i += 4;
////       j += 4;
////       k += 4;


        int index1 = i + 4 + (j + 4) * (nx) + (k + 4) * (nx) * (ny);

     if (pp[index1] >= 0.0) {

        double dd1 = cS[(k +1)*nxS] * vel[index1 -3 * nx * ny];
////                                                                 + phi[ ( (k + 1 -2) * nx * ny + j * nx + i ) ] * cS[(k-4) * nxS + 1]
////                                                                 + phi[ ( (k + 2 -2) * nx * ny + j * nx + i ) ] * cS[(k+1) * nxS + 2]
////                                                                 + phi[ ( (k + 3 -2) * nx * ny + j * nx + i ) ] * cS[(k+1) * nxS + 3]
////                                                                 + phi[ ( (k + 4 -2) * nx * ny + j * nx + i ) ] * cS[(k+1) * nxS + 4]
////                                                                 + phi[ ( (k + 5 -2) * nx * ny + j * nx + i ) ] * cS[(k+1) * nxS + 5] ;
        for (int l=1; l<= 6; l++){
             dd1 = dd1 + vel[ ( (k  + 4 + l -3) * nx * ny + (j+4) * nx + (i+4) ) ] * cS[(k +1) * nxS + l]; 
        }

             nl[index1] = nl[index1] + dd1 * pp[index1];

                             }

     else {

        double dd1 = cS2[(k +1)*nxS] * vel[index1 -3 * nx * ny];
////                                                                 + phi[ ( (k + 1 -2) * nx * ny + j * nx + i ) ] * cS[(k-4) * nxS + 1]
////                                                                 + phi[ ( (k + 2 -2) * nx * ny + j * nx + i ) ] * cS[(k+1) * nxS + 2]
////                                                                 + phi[ ( (k + 3 -2) * nx * ny + j * nx + i ) ] * cS[(k+1) * nxS + 3]
////                                                                 + phi[ ( (k + 4 -2) * nx * ny + j * nx + i ) ] * cS[(k+1) * nxS + 4]
////                                                                 + phi[ ( (k + 5 -2) * nx * ny + j * nx + i ) ] * cS[(k+1) * nxS + 5] ;
        for (int l=1; l<= 6; l++){
             dd1 = dd1 + vel[ ( (k  + 4 + l -3) * nx * ny + (j+4) * nx + (i+4) ) ] * cS2[(k +1) * nxS + l];
        }

             nl[index1] = nl[index1] + dd1 * pp[index1];

                             }

//         output[index1] += phi[ ( (k + l -2) * nx * ny + j * nx + i ) ] * cS[(k-4) * nxS + l] ;

//                     }
////                  }
////               }
////            }

}




    ////--- u.du/dx
    if (component_id == 7) {

////        if ((i >= 4) && (i <= nx-5)) {
////           if ((j >= 4) && (j <= ny-5)) {
////                if ((k >= 4) && (k <= nz-5)) {

        int index1 = (i+4) + (j+4) * (nx) + (k+4) * (nx) * (ny);

     if (pp[index1] >= 0.0) {

        double dd1 = cS[(1 + i) * nxS] * vel[index1 - 3] ;
////                                                      + phi[ ( (k+4) * nx * ny + (j+4) * nx + (i+4) + 1 - 2) ] * cS[(i) * nxS + 1]
////                                                      + phi[ ( (k+4) * nx * ny + (j+4) * nx + (i+4) + 2 - 2) ] * cS[(i) * nxS + 2]
////                                                      + phi[ ( (k+4) * nx * ny + (j+4) * nx + (i+4) + 3 - 2) ] * cS[(i) * nxS + 3]
////                                                      + phi[ ( (k+4) * nx * ny + (j+4) * nx + (i+4) + 4 - 2) ] * cS[(i) * nxS + 4]
////                                                      + phi[ ( (k+4) * nx * ny + (j+4) * nx + (i+4) + 5 - 2) ] * cS[(i) * nxS + 5];

////-- looped implementation 
        for (int l=1; l<= 6; l++){
             dd1  = dd1 + vel[ ( (k+4) * nx * ny + (j+4) * nx + (i+4) + l - 3) ] * cS[(i + 1) * nxS + l];
        }
             nl[index1] = dd1*pp[index1];

                            }

     else {

        double dd1 = cS2[(1 + i) * nxS] * vel[index1 - 3] ;
////                                                      + phi[ ( (k+4) * nx * ny + (j+4) * nx + (i+4) + 1 - 2) ] * cS[(i) * nxS + 1]
////                                                      + phi[ ( (k+4) * nx * ny + (j+4) * nx + (i+4) + 2 - 2) ] * cS[(i) * nxS + 2]
////                                                      + phi[ ( (k+4) * nx * ny + (j+4) * nx + (i+4) + 3 - 2) ] * cS[(i) * nxS + 3]
////                                                      + phi[ ( (k+4) * nx * ny + (j+4) * nx + (i+4) + 4 - 2) ] * cS[(i) * nxS + 4]
////                                                      + phi[ ( (k+4) * nx * ny + (j+4) * nx + (i+4) + 5 - 2) ] * cS[(i) * nxS + 5];

////-- looped implementation
        for (int l=1; l<= 6; l++){
             dd1  = dd1 + vel[ ( (k+4) * nx * ny + (j+4) * nx + (i+4) + l - 3) ] * cS2[(i + 1) * nxS + l];
        }
             nl[index1] = dd1*pp[index1];

                            }

////         output[index1] += phi[ ( k * nx * ny + j * nx + i + l - 2) ] * cS[(i-4) * nxS + l] ;

////                     }
////                  }
////               }
////            }



}



//!//!//! the shifted entries-- TESTED!




//   if (component_id ==2) {
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

   ////--v.du/dy 
   if (component_id == 8) {

////        if ((i >= 4) && (i <= nx-5)) {
////            if ((j >= 4) && (j <= ny-5)) {
////                if ((k >= 4) && (k <= nz-5)) {


        int index1 = (i+4) + (j+4) * (nx) + (k+4) * (nx) * (ny);

    if (pp[index1] >= 0.0) {

        double dd1 = cS[(j + 1) * nxS] * vel[index1 -3 * nx];
////                                                            + phi[ ( (k+4) * nx * ny + ((j+4) + 1 - 2) * nx + (i+4)) ] * cS[(j) * nxS + 1]
////                                                            + phi[ ( (k+4) * nx * ny + ((j+4) + 2 - 2) * nx + (i+4)) ] * cS[(j) * nxS + 2]
////                                                            + phi[ ( (k+4) * nx * ny + ((j+4) + 3 - 2) * nx + (i+4)) ] * cS[(j) * nxS + 3]
////                                                            + phi[ ( (k+4) * nx * ny + ((j+4) + 4 - 2) * nx + (i+4)) ] * cS[(j) * nxS + 4]
////                                                            + phi[ ( (k+4) * nx * ny + ((j+4) + 5 - 2) * nx + (i+4)) ] * cS[(j) * nxS + 5];
        for (int l=1; l<= 6; l++){

            dd1 = dd1 + vel[ ( (k+4) * nx * ny + ((j+4) + l -3) * nx + (i + 4)) ] * cS[(j + 1) * nxS + l];
        }

            nl[index1] = nl[index1] + dd1 * pp[index1];

                          }

    else {

        double dd1 = cS2[(j + 1) * nxS] * vel[index1 -3 * nx];
////                                                            + phi[ ( (k+4) * nx * ny + ((j+4) + 1 - 2) * nx + (i+4)) ] * cS[(j) * nxS + 1]
////                                                            + phi[ ( (k+4) * nx * ny + ((j+4) + 2 - 2) * nx + (i+4)) ] * cS[(j) * nxS + 2]
////                                                            + phi[ ( (k+4) * nx * ny + ((j+4) + 3 - 2) * nx + (i+4)) ] * cS[(j) * nxS + 3]
////                                                            + phi[ ( (k+4) * nx * ny + ((j+4) + 4 - 2) * nx + (i+4)) ] * cS[(j) * nxS + 4]
////                                                            + phi[ ( (k+4) * nx * ny + ((j+4) + 5 - 2) * nx + (i+4)) ] * cS[(j) * nxS + 5];
        for (int l=1; l<= 6; l++){

            dd1 = dd1 + vel[ ( (k+4) * nx * ny + ((j+4) + l -3) * nx + (i + 4)) ] * cS2[(j + 1) * nxS + l];
        }

            nl[index1] = nl[index1] + dd1 * pp[index1];

                          } 
////         output[index1] += phi[ ( k * nx * ny + (j + l - 2) * nx + i) ] * cS[(j-4) * nxS + l] ;

////                    }
////                  }
////               }
////            }
}







//!//!//! the shifted entries !TESTED




//   if (component_id == 3) { 
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

   if (component_id == 9) {

////        if ((i >= 4) && (i <= nx-5)) {
////            if ((j >= 4) && (j <= ny-5)) {
////                if ((k >= 4) && (k <= nz-5)) {

////       i += 4;
////       j += 4;
////       k += 4;


        int index1 = i + 4 + (j + 4) * (nx) + (k + 4) * (nx) * (ny);

   if ( pp[index1] >= 0.0) {

        double dd1 = cS[(k +1)*nxS] * vel[index1 -3 * nx * ny];
////                                                                 + phi[ ( (k + 1 -2) * nx * ny + j * nx + i ) ] * cS[(k-4) * nxS + 1]
////                                                                 + phi[ ( (k + 2 -2) * nx * ny + j * nx + i ) ] * cS[(k+1) * nxS + 2]
////                                                                 + phi[ ( (k + 3 -2) * nx * ny + j * nx + i ) ] * cS[(k+1) * nxS + 3]
////                                                                 + phi[ ( (k + 4 -2) * nx * ny + j * nx + i ) ] * cS[(k+1) * nxS + 4]
////                                                                 + phi[ ( (k + 5 -2) * nx * ny + j * nx + i ) ] * cS[(k+1) * nxS + 5] ;
        for (int l=1; l<= 6; l++){
             dd1 = dd1 + vel[ ( (k  + 4 + l -3) * nx * ny + (j+4) * nx + (i+4) ) ] * cS[(k + 1 ) * nxS + l]; 
        }

             nl[index1] = nl[index1] + dd1 * pp[index1];
                             }

   else{

        double dd1 = cS2[(k +1)*nxS] * vel[index1 -3 * nx * ny];
////                                                                 + phi[ ( (k + 1 -2) * nx * ny + j * nx + i ) ] * cS[(k-4) * nxS + 1]
////                                                                 + phi[ ( (k + 2 -2) * nx * ny + j * nx + i ) ] * cS[(k+1) * nxS + 2]
////                                                                 + phi[ ( (k + 3 -2) * nx * ny + j * nx + i ) ] * cS[(k+1) * nxS + 3]
////                                                                 + phi[ ( (k + 4 -2) * nx * ny + j * nx + i ) ] * cS[(k+1) * nxS + 4]
////                                                                 + phi[ ( (k + 5 -2) * nx * ny + j * nx + i ) ] * cS[(k+1) * nxS + 5] ;
        for (int l=1; l<= 6; l++){
             dd1 = dd1 + vel[ ( (k  + 4 + l -3) * nx * ny + (j+4) * nx + (i+4) ) ] * cS2[(k + 1 ) * nxS + l];
        }

             nl[index1] = nl[index1] + dd1 * pp[index1];
                             }

//         output[index1] += phi[ ( (k + l -2) * nx * ny + j * nx + i ) ] * cS[(k-4) * nxS + l] ;

//                     }
////                  }
////               }
////            }

}










}



