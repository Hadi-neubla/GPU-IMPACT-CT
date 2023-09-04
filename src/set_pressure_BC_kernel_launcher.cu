// !**************************************************************************************************                                    // !* by Hadi Zolfaghari, ARTORG Center, University of Bern, (zolfaghari.haadi@gmail.com)            *      
// !* October 2015 - March 2020                                                                      *            
// !* Modified by Hadi Zolfaghari, University of Cambridge (hz382@damtp.cam.ac.uk)                   *      
// !* April 2020 -                                                                                   *            
// !**************************************************************************************************
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "set_pressure_BC_kernel.cu"
extern "C" void pressure_bc_launcher_(int *type_id1, int *type_id2, int *type_id3, int *k_start_index, int *nx_dev, int *ny_dev, int *nz_dev, double *pre_inlet,  double *pre_outlet, double *pre_outlet_2, double** ct_geometry_dev, double** rel_dev)
{

/*dim3    blocks((nx_dev-8)/8, (ny_dev-8)/8, (nz_dev-8)/8);
*/

int nx = (*nx_dev-8)/8;
int ny = (*ny_dev-8)/8;
int nz = (*nz_dev-8)/8;


dim3 blocks(nx,ny,nz);

int nxt = (*nx_dev-8)/blocks.x;
int nyt = (*ny_dev-8)/blocks.y;
int nzt = (*nz_dev-8)/blocks.z;

/*printf("Value of n=%d", nxt);
printf("Value of n=%d", nx);*/

 
dim3 threads(nxt,nyt,nzt);

/*int nxg = *nx_dev;
int nyg = *ny_dev;
int nzg = *nz_dev; 

cudaMemcpy(*phi_grad_dev, pre, sizeof(double) * nxg * nyg * nzg, cudaMemcpyHostToDevice ); */

/*int nxgc1 = *nx_grad_coef1;
int nygc1 = *ny_grad_coef1;


cudaMemcpy(*cGp1_dev, cGp1, sizeof(double) * nxgc1 * nygc1, cudaMemcpyHostToDevice ); */


//int row = sizeof(*phi_grad_dev) / sizeof(*phi_grad_dev[0]);
//int column = sizeof(*phi_grad_dev[0])/row;
//printf ("pressure host pointer:, %f\n", pre[0]);
//printf ("pressure host:, %f\n", phi_grad_dev[0]);

/*    if (mod(nx_dev, szblock) .ne. 0) blocks%x = blocks%x + 1 */
/*dim3    threads((nx_dev-8)/blocks.x, (ny_dev-8)/blocks.y, (nz_dev-8)/blocks.z);
*/

set_pressure_BC_kernel<<< blocks,threads >>>(*type_id1, *type_id2, *type_id3, *k_start_index, *nx_dev, *ny_dev, *nz_dev, *pre_inlet, *pre_outlet, *pre_outlet_2, *ct_geometry_dev, *rel_dev);
return;
}
