// !**************************************************************************************************                                    // !* by Hadi Zolfaghari, ARTORG Center, University of Bern, (zolfaghari.haadi@gmail.com)            *      
// !* October 2015 - March 2020                                                                      *            
// !* Modified by Hadi Zolfaghari, University of Cambridge (hz382@damtp.cam.ac.uk)                   *      
// !* April 2020 -                                                                                   *            
// !**************************************************************************************************
extern "C" __global__ void axpby_kernel(int nx, int ny, int nz, double a, double b, double* x, double* y, double* z)


{
        int i = threadIdx.x + blockIdx.x * blockDim.x;
        int j = threadIdx.y + blockIdx.y * blockDim.y;
        int k = threadIdx.z + blockIdx.z * blockDim.z;

//!//!//! the warp-oriented implemtation

        int index = (i+4) + (j+4) * (nx) + (k+4) * (nx) * (ny);
        z[index] = a * x[index] + b * y[index];

 //       z[index] = a + 5*b ;

}

