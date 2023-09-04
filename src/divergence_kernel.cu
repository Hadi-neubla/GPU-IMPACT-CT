// !**************************************************************************************************                                    // !* by Hadi Zolfaghari, ARTORG Center, University of Bern, (zolfaghari.haadi@gmail.com)            *      
// !* October 2015 - March 2020                                                                      *            
// !* Modified by Hadi Zolfaghari, University of Cambridge (hz382@damtp.cam.ac.uk)                   *      
// !* April 2020 -                                                                                   *            
// !**************************************************************************************************
extern "C" __global__ void divergence_kernel(int dimension_id, int nx, int ny, int nz, int nxg, int nyg, double* phi, double* div, double* cD)



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


//    if (dimension_id ==1) {
//
//
//
//          for (int ii= 3; ii <= nx-5; ii++){  
//              for (int jj= 4; jj <= ny-5; jj++){
//                  for (int kk= 4; kk <= nz-5; kk++){
//
//                      div[ii + jj * nx + kk * nx * ny] = cD[0 + (ii-4) * nxg] *phi[ii-3 + jj * nx + kk * nx * ny];                       
//                      for (int ll=1; ll<= 5; ll++){ 
//                            div[ii + jj * nx + kk * nx * ny] +=  phi[ii -3 + ll + jj * nx + kk * nx * ny]
//                                                                 * cD[ll + (ii-4) * nxg];
//                      }
//                  }
//              }
//          }
//
//
//
//
//    }


//        if ((i <= 3) || (i > nx-5)) {
//            if ((j <= 3) || (j > ny-5)) {
//                if ((k <= 3) || (k > nz-5)) {
//                           
//                return;
//        
//                }
//            }
//        }






    if (dimension_id == 1) {

//        if ((i >= 4) && (i <= nx-5)) {
//            if ((j >= 4) && (j <= ny-5)) {     
//                if ((k >= 4) && (k <= nz-5)) {
//

//        int index1 = i + j * (nx) + k * (nx) * (ny);
        div[ (i+4) + (j+4) * (nx) + (k+4) * (nx) * (ny)] =                     phi[ ( (k+4) * nx * ny + (j+4) * nx + (i+4) + 0 - 3) ] * cD[(i+1) * nxg]
                                                                             + phi[ ( (k+4) * nx * ny + (j+4) * nx + (i+4) + 1 - 3) ] * cD[(i+1) * nxg + 1]
                                                                             + phi[ ( (k+4) * nx * ny + (j+4) * nx + (i+4) + 2 - 3) ] * cD[(i+1) * nxg + 2]
                                                                             + phi[ ( (k+4) * nx * ny + (j+4) * nx + (i+4) + 3 - 3) ] * cD[(i+1) * nxg + 3]
                                                                             + phi[ ( (k+4) * nx * ny + (j+4) * nx + (i+4) + 4 - 3) ] * cD[(i+1) * nxg + 4]
                                                                             + phi[ ( (k+4) * nx * ny + (j+4) * nx + (i+4) + 5 - 3) ] * cD[(i+1) * nxg + 5];
        
//the check

                                                                
                                                                
                                                                
                                                                
                                                                


//                  }
//               }
//            }
//


}



    if (dimension_id ==2) {

//        div[i + j * (nx) + k * (nx) * (ny)] = 0;
//        if ((i >= 4) && (i <= nx-5)) {
//            if ((j >= 4) && (j <= ny-5)) {     
//                if ((k >= 4) && (k <= nz-5)) {


        int index1 = (i+4) + (j+4) * (nx) + (k+4) * (nx) * (ny);

        div[index1] =    div[index1]               + phi[ ( (k+4) * nx * ny + (j+4 + 0 - 3) * nx + i+4) ] * cD[(j+1) * nxg + 0]
                                                            + phi[ ( (k+4) * nx * ny + (j+4 + 1 - 3) * nx + i+4) ] * cD[(j+1) * nxg + 1]
                                                            + phi[ ( (k+4) * nx * ny + (j+4 + 2 - 3) * nx + i+4) ] * cD[(j+1) * nxg + 2]
                                                            + phi[ ( (k+4) * nx * ny + (j+4 + 3 - 3) * nx + i+4) ] * cD[(j+1) * nxg + 3]
                                                            + phi[ ( (k+4) * nx * ny + (j+4 + 4 - 3) * nx + i+4) ] * cD[(j+1) * nxg + 4]
                                                            + phi[ ( (k+4) * nx * ny + (j+4 + 5 - 3) * nx + i+4) ] * cD[(j+1) * nxg + 5];
//        for (int l=1; l<= 5; l++){

//         grad[index1] += pre[ ( k * nx * ny + (j + l - 2) * nx + i) ] * cGp[(j-4) * nxg + l] ;

//                     }
//                  }
//               }
//            }
}

   if (dimension_id == 3) { 

//        if ((i >= 4) && (i <= nx-5)) {
//            if ((j >= 4) && (j <= ny-5)) {     
//                if ((k >= 4) && (k <= nz-5)) {
//

        int index1 = i+4 + (j+4) * (nx) + (k+4) * (nx) * (ny);

        div[index1] =  div[index1]                 + phi[ ( (k+4 + 0 -3) * nx * ny + (j+4) * nx + i+4 ) ] * cD[(k+1) * nxg + 0]
                                                            + phi[ ( (k+4 + 1 -3) * nx * ny + (j+4) * nx + i+4 ) ] * cD[(k+1) * nxg + 1] 
                                                            + phi[ ( (k+4 + 2 -3) * nx * ny + (j+4) * nx + i+4 ) ] * cD[(k+1) * nxg + 2] 
                                                            + phi[ ( (k+4 + 3 -3) * nx * ny + (j+4) * nx + i+4 ) ] * cD[(k+1) * nxg + 3] 
                                                            + phi[ ( (k+4 + 4 -3) * nx * ny + (j+4) * nx + i+4 ) ] * cD[(k+1) * nxg + 4] 
                                                            + phi[ ( (k+4 + 5 -3) * nx * ny + (j+4) * nx + i+4 ) ] * cD[(k+1) * nxg + 5]; 
//        for (int l=1; l<= 5; l++){


//         grad[index1] += pre[ ( (k + l -2) * nx * ny + j * nx + i ) ] * cGp[(k-4) * nxg + l] ;


//                     }
//                  }
//               }
//            }
         
}

}



