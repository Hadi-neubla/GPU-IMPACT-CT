// !**************************************************************************************************                                           
// !* by Hadi Zolfaghari, ARTORG Center, University of Bern, (zolfaghari.haadi@gmail.com)            *      
// !* October 2015 - March 2020                                                                      *            
// !**************************************************************************************************
extern "C" __global__ void bicgstab_kernel(int nx, int ny, int nz, double omega, double alpha, double* phi, double* z2, double* z1, double* Ar, double* rr)


{
        int i = threadIdx.x + blockIdx.x * blockDim.x;
        int j = threadIdx.y + blockIdx.y * blockDim.y;
        int k = threadIdx.z + blockIdx.z * blockDim.z;

//!//!//! the warp-oriented implemtation

        int index = (i+4) + (j+4) * (nx) + (k+4) * (nx) * (ny);
        phi[index] = phi[index] + omega * z2[index] + alpha * z2[index];
        rr[index]  = rr[index] - omega * Ar[index];
}

