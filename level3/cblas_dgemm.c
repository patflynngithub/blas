// Adapted by Patrick Flynn : 5/16/19

/*  Example of using CBLAS BLAS Level 3 cblas_dgemm to multiply two matrices
  
    ./dgemm  m k n          m = # of rows in first matrix
                            k = # of columns in first matrix and
                                # of rows in second matrix
                            n = # of columns in second matrix 
    
    Outputs timing to a file
*/

#include "stdio.h"
#include "stdlib.h"
#include "sys/time.h"
#include "time.h"

#include <cblas.h>

// =============================================================================================================

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
            printf("%f ", *mat);
            mat++;
        }
        printf("\n");
    }
}

// -------------------------------------------------------------------------------------------

int main(int argc, char* argv[])
{
    if(argc<4){
        printf("Input Error\n\n");
        printf("Proper Format: ./cblas_dgemm m k n   (note: A:mxk, B:kxn)\n\n");
        return 1;
    }

    // A matrix is m x k
    // B marix is k x n
    // C matrix will be m * n
    int m = atoi(argv[1]);
    int k = atoi(argv[2]);
    int n = atoi(argv[3]);

    int sizeofA = m * k;
    int sizeofB = k * n;
    int sizeofC = m * n;

    // timing setup
    struct timeval start,finish;
    double duration;

    // source and target matrices
    double* A = (double*)malloc(sizeof(double) * sizeofA);
    double* B = (double*)malloc(sizeof(double) * sizeofB);
    double* C = (double*)malloc(sizeof(double) * sizeofC);

    // initialize matrices
    int i;
    for (i=0; i<sizeofA; i++)
    A[i] = i+1;;

    for (i=0; i<sizeofB; i++)
    B[i] = i+1;

    for (i=0; i<sizeofC; i++)
    C[i] = 0.0;

    //   srand((unsigned)time(NULL));
    // 
    //   for (i=0; i<sizeofA; i++)
    //     A[i] = i%3+1;//(rand()%100)/10.0;
    // 
    //   for (i=0; i<sizeofB; i++)
    //     B[i] = i%3+1;//(rand()%100)/10.0;
    // 
    //   for (i=0; i<sizeofC; i++)
    //     C[i] = i%3+1;//(rand()%100)/10.0;


    /*
    // Goal: C = AB (in C-code)
    
    //  dgem:  C <- alpha * op(A) * op_A(B) + beta*C

    void cblas_dgemm(CBLAS_LAYOUT layout, CBLAS_TRANSPOSE TransA,
                 CBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const double alpha, const double *A,
                 const int lda, const double *B, const int ldb,
    */
    
    // BLAS cblas_dgemm_ setup
    CBLAS_LAYOUT    layout = CblasRowMajor;
    CBLAS_TRANSPOSE transA = CblasNoTrans;
    CBLAS_TRANSPOSE transB = CblasNoTrans;
    double alpha = 1.0;
    double beta = 0.0;

    printf("m=%d, k=%d, n=%d, alpha=%lf, beta=%lf, sizeofC=%d\n",m,k,n,alpha,beta,sizeofC);

    gettimeofday(&start, NULL);

    // ***************************************************************************

    cblas_dgemm(layout,transA,transB, m, n, k, alpha, A,  
                k ,B, n, beta ,C, n);                                                                    

    // ***************************************************************************

    gettimeofday(&finish, NULL);

    // calculate timing, flops
    duration = ((double)(finish.tv_sec-start.tv_sec)*1000000 + (double)(finish.tv_usec-start.tv_usec)) / 1000000;
    double flops  = 2.0 * m *n*k;
    double mflops = flops/duration*1.0e-6;  // megaFlOPs/sec
    
    // output timing, MFLOPS to a file
    FILE *fp;
    fp = fopen("timingCBLAS_DGEMM.txt", "a");
    fprintf(fp, "mxkxn: %dx%dx%d, size(C)=%d\t%lf s\t%lf MFLOPS/sec\n", m, k, n, sizeofC, duration, mflops);
    fclose(fp);

    // FOR DEBUGGING: print source/target matrices if they are small enough
    if (m < 10 && k < 10 && n < 10) {
    printf("A =\n"); print_matrix(A, m, k);
    printf("B =\n"); print_matrix(B, k, n);
    printf("C =\n"); print_matrix(C, m, n);
    }

    free(A);
    free(B);
    free(C);

    return 0;
}
