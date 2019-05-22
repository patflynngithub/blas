// Adapted by Patrick Flynn : 5/16/19

/*  Example of using BLAS Level 2 dtrsv_ to get solution x to Ax=b
    where A is a triangular square matrix (upper in this case)
  
    ./dtrsv N          size of Ax=b: A is NxN, b is Nx1
    
    Outputs timing to file
*/

#include "stdio.h"
#include "stdlib.h"
#include "sys/time.h"
#include "time.h"

extern double dtrsv_(char*, char*, char*, int*, double*, int*, double*, int*);

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
        printf("%.3f ", *vec);
        vec++;
    }
    printf("\n");
}

// -------------------------------------------------------------------------------------------

/*  Print a matrix
 
    Input:  mat - m x n double matrix
            m   - # of rows
            n   - # of columns
    Output: nothing
*/
void print_matrix(double *mat, int m, int n)
{
    int i,j;
    
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            printf("%.3f ", *mat);
            mat++;
        }
        printf("\n");
    }
}

// -------------------------------------------------------------------------------------------

int main(int argc, char* argv[])
{
  if(argc<2){
    printf("Input Error\n");
    return 1;
  }

  // array and vector size
  int n = atoi(argv[1]);

  // A is n x n matrix
  double* A = (double*)malloc(sizeof(double) * n * n);

  // b_x is n x 1 vector
  // input b_x is b before calling dtrsv_
  // output b_x is x after calling dtrsv_
  double* b_x = (double*)malloc(sizeof(double) * n);

  // initialize matrix and vector
  int i,j;
  double* B = A;
  for (i=0; i<n; i++)
    for (j=0; j<n; j++) {
        if (j<i) *B = 0.0;
        else     *B = i*n + j+1;
        B++;
    }

  for (i=0; i<n; i++)
    b_x[i] = i+10;

  if (n <= 10) {
    printf("A=\n"); print_matrix(A, n, n);
    printf("b =\n"); print_vector(b_x, n);
  }

  // setup for calling dtrsv_
  char uplo  = 'l';  // lower triangular in terms of Fortran
  char trans = 't';  // solve A**T*x = b
  char diag  = 'n';  // does A have a unit diagonal
  int  lda   =  n;   // first dimension of A (from Fortran point of view)
  int  incx  =  1;   // increment for elements of the vector

  // timing setup
  struct timeval start,finish;
  double duration;

  gettimeofday(&start, NULL);

  // ***************************************************************************

  dtrsv_(&uplo, &trans, &diag, &n, A, &n, b_x, &incx);

  // ***************************************************************************

  gettimeofday(&finish, NULL);

  // calculate timing, flops
  duration = ((double)(finish.tv_sec-start.tv_sec)*1000000 + (double)(finish.tv_usec-start.tv_usec)) / 1000000;

  // output timing, flops to a file
  FILE *fp;
  fp = fopen("timingDTRSV.txt", "a");
  fprintf(fp, "square matrix n = %d\t\t%lf s\n",n, duration);
  fclose(fp);

  // print solution
  if (n <= 10)
     printf("x =\n"); print_vector(b_x, n);

  free(b_x);
  free(A);
  return 0;
}

