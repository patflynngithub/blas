// Adapted by Patrick Flynn : 5/16/19

/*  Example of using BLAS Level 3 dgemm_ to multiply two matrices

    Note: C calling Fortran BLAS routine directly

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

// BLAS
extern void dgemm_(char*, char*, int*, int*,int*, double*, double*, int*, double*, int*, double*, double*, int*);

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
        printf("Proper Format: ./dgemm m k n   (note: A:mxk, B:kxn)\n\n");
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


    // Goal: C = AB (in C-code)
    //
    // Note: Fortran column major matrix setup in BLAS forces C-code (row-major) adjustment to use dgemm_.
    //       From Fortran point of view, C-code matrices can be though of as tranposes:
    //
    //       i.e., C-Code A, B, C are Fortran A' B' C'
    //
    //       So want C' result
    //       C' = (AB)' = B'A'
    //
    // So dgemm_ used the following way in C-code
    //
    //    C <- alpha * op_B(B) * op_A(A) + beta*C

    // dgemm_ parameters:
    //
    // ta, tb  : indicate whether the input matrices are to be transposed or not 
    // n       : # of rows of op(first matrix=B) and output matrix
    // m       : # of columns of op(second matrix=A) and output matrix
    // k       : # of columns of op(first matrix=B) and rows of op(second matrix=A)
    // alpha   : scalar in front of op(first) * op(second)
    // B       : is first matrix (because in c-code using Fortran BLAS)
    // n (LDB) : # of increments to next column entry (considering col-major interpretation and possible padding)
    // A       : is second matrx (because in c-code using Fortran BLAS)
    // k (LDA) : # of increments to next column entry (considering col-major interpretation and possible padding)
    // beta    : scalar factor for C on rhs
    // C       : output matrix
    // n       : # of increments to next column entry (considering col-major interpretation and possible padding)

    // dgemm_ setup
    char ta = 'N';
    char tb = 'N';
    double alpha = 1.0;
    double beta = 0.0;

    printf("m=%d, k=%d, n=%d, alpha=%lf, beta=%lf, sizeofC=%d\n",m,k,n,alpha,beta,sizeofC);

    gettimeofday(&start, NULL);

    // ***************************************************************************

    dgemm_(&ta, &tb, &n, &m, &k, &alpha, B, &n, A, &k, &beta, C, &n); // this way there is no moving data around 
                                                                      // b/c of tranposing, etc.

    // ***************************************************************************

    gettimeofday(&finish, NULL);

    // calculate timing, flops
    duration = ((double)(finish.tv_sec-start.tv_sec)*1000000 + (double)(finish.tv_usec-start.tv_usec)) / 1000000;
    double flops  = 2.0 * m *n*k;
    double mflops = flops/duration*1.0e-6;  // megaFlOPs/sec
    
    // output timing, MFLOPS to a file
    FILE *fp;
    fp = fopen("timingDGEMM.txt", "a");
    fprintf(fp, "mxkxn: %dx%dx%d, size(C)=%d\t%lf s\t%lf MFLOPS/sec\n", m, k, n, sizeofC, duration, mflops);
    fclose(fp);

    // FOR DEBUGGING: print source/target matrices if they are small enough
    if (m <= 10 && k <= 10 && n <= 10) {
    printf("A =\n"); print_matrix(A, m, k);
    printf("B =\n"); print_matrix(B, k, n);
    printf("C =\n"); print_matrix(C, m, n);
    }

    free(A);
    free(B);
    free(C);

    return 0;
}

