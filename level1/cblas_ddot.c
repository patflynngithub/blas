// Adapted by Patrick Flynn : 5/16/19

/*  Example of using CBLAS BLAS Level 1 cblas ddot to get dot product of two double vectors
  
    ./ddot N          N = size of the two vectors
    
    Outputs timing to file
*/

#include "stdio.h"
#include "stdlib.h"
#include "sys/time.h"
#include "time.h"

#include <cblas.h>

// =============================================================================================================

/*  Print a vector
 
    Input:  vec - n x 1 double vector
            n   - # of elements in vector
    Output: nothing
*/
void print_vector(double *vec, int n)
{
    int i;
    
    for (i = 0; i < n; i++) {
        printf("%.1f ", *vec);
        vec++;
    }
    printf("\n");
}

// -------------------------------------------------------------------------------------------

int main(int argc, char* argv[])
{
  if(argc<2){
    printf("Input Error\n\n");
    printf("Proper Format: ./cblas_ddot vector_size\n\n");
    return 1;
  }

  // vector size
  int n = atoi(argv[1]);

  // va and vb are n x 1 vectors
  double* va = (double*)malloc(sizeof(double) * n);
  double* vb = (double*)malloc(sizeof(double) * n);

  // initialize vectors
  int i;
  for (i=0; i<n; i++)
    va[i] = i+1;
  for (i=0; i<n; i++)
    vb[i] = i+3;

//   srand((unsigned)time(NULL));
// 
//   for (i=0; i<n; i++)
//     va[i] = i%3+1;//(rand()%100)/10.0;

  // ddot_: increments between vector elements
  int inca=1, incb=1;

  // timing setup
  struct timeval start,finish;
  double duration;

  gettimeofday(&start, NULL);

  // ***************************************************************************

  double dp = cblas_ddot(n, va, inca, vb, incb);

  // ***************************************************************************

  gettimeofday(&finish, NULL);

  // calculate timing, flops
  duration = ((double)(finish.tv_sec-start.tv_sec)*1000000 + (double)(finish.tv_usec-start.tv_usec)) / 1000000;
  double flops  = 2.0 * n;
  double mflops = flops/duration*1.0e-6;

  // output timing, megaFLOPS to a file
  FILE *fp;
  fp = fopen("timingCBLAS_DDOT.txt", "a");
  fprintf(fp, "%d element vector\t%lf s\t%lf MFLOPS\n",n, duration, mflops);
  fclose(fp);

  // print vectors if they are small enough  
  if (n <= 10) {
    printf("va =\n"); print_vector(va, n);
    printf("vb =\n"); print_vector(vb, n);
  }

  printf("vector size = %d, dot product = %f\n", n, dp);

  free(va);
  free(vb);
  return 0;
}

