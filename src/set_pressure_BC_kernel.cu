extern "C" __global__ void set_pressure_BC_kernel(int type_id1, int type_id2, int type_id3, int k_start_index, int nx, int ny, int nz,  double pre_inlet, double pre_outlet, double pre_outlet_2, double* ct_geometry,  double* rel)

// type id1-3 take values 0 or 1 and identifies whether (1) or not (0) pressure Dirichlet bc is set at inlet (type_id1), descending aorta outlet (type_id2), and/or supra-aortic artery outlets (type_id3)    

{

        int i = threadIdx.x + blockIdx.x * blockDim.x;
        int j = threadIdx.y + blockIdx.y * blockDim.y;
        int k = threadIdx.z + blockIdx.z * blockDim.z;



//!!//!!// The serial loop for testing the copy





//!//!//! this is the shifted operation -- to be converted to power 2 operation for coalesed memory call on GPU



//        if ((i == 0) && (j == 0) && (k == 0)) {printf("Value of k index for 0,0,0=%d", k_start_index);}

//       if ((ct_geometry[(i+4) + (j+4) * nx + (k+4) * nx * ny] == 5.)  ) { printf("i index=%d", i);}


//!//!//! this is the power 2 implementation
       if ((type_id1 == 1) && (ct_geometry[(i+4) + (j+4) * nx + (k+4) * nx * ny] == 5.) && ((k+4) == (k_start_index - 1)) ) { 
                  rel[ (i + 4) + (j + 4) * nx + (k + 4) * nx * ny] =  pre_inlet;
       }
       if ((type_id2 == 1) && (ct_geometry[(i+4) + (j+4) * nx + (k+4) * nx * ny] == 9.) && ((k+4) == (k_start_index - 1))) { 
                  rel[ (i + 4) + (j + 4) * nx + (k + 4) * nx * ny] =  pre_outlet;
       }
       if ((type_id3 == 1) && (ct_geometry[(i+4) + (j+4) * nx + (k+4) * nx * ny] == 11.)) { 
                  rel[ (i + 4) + (j + 4) * nx + (k + 4) * nx * ny] =  pre_outlet_2;
       }
			  //cudaThreadSynchronize(); 
}



