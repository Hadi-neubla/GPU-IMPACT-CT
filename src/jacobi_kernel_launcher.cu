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
#include "jacobi_kernel.cu"
extern "C" void jacobi_launcher_(int *grid_id, int *nx_dev_g, int *ny_dev_g, int *nz_dev_g, int *ny_cdg1_g, int *ny_cdg2_g, int* ny_cdg3_g, int *nx_cdg1_g, double* omega_g, double** cdg1_dev3, double** cdg2_dev3, double** cdg3_dev3, double** rel_dev3, double** bb_dev3, double** comp_dev3)
{

/*dim3    blocks((nx_dev-8)/8, (ny_dev-8)/8, (nz_dev-8)/8);
*/

//printf("kernel launch check = done") ;


int nx = (*nx_dev_g-8)/8;
int ny = (*ny_dev_g-8)/8;
int nz = (*nz_dev_g-8)/8;


dim3 blocks(nx,ny,nz);

int nxt = (*nx_dev_g-8)/blocks.x;
int nyt = (*ny_dev_g-8)/blocks.y;
int nzt = (*nz_dev_g-8)/blocks.z;

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

jacobi_kernel<<< blocks,threads >>>(*grid_id, *nx_dev_g, *ny_dev_g, *nz_dev_g, *ny_cdg1_g, *ny_cdg2_g, *ny_cdg3_g, *nx_cdg1_g, *omega_g, *cdg1_dev3, *cdg2_dev3, *cdg3_dev3, *rel_dev3, *bb_dev3, *comp_dev3);
return;
}
