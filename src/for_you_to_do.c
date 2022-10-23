#include "../include/for_you_to_do.h"

/**
 * 
 * this function computes LU factorization
 * for a square matrix
 * 
 * syntax 
 *  
 *  input : 
 *      A     n by n , square matrix
 *      ipiv  1 by n , vector
 *      n            , length of vector / size of matrix
 *  
 *  output :
 *      return -1 : if the matrix A is singular (max pivot == 0)
 *      return  0 : return normally 
 * 
 **/

int mydgetrf(double *A, int *ipiv, int n) 
{
    int i, t, j, k, maxIndex, temps;
    double max;
    double *temp = (double*)malloc(sizeof(double) * n);
    //outer loop to control the row
    for (i = 0; i < n - 1; i++)
    {
        // pivoting
        // select the row that has max value in a column
        maxIndex = i;
        max = fabs(A[i * n + i]);
        for (t = i + 1; t < n; t++){
            if (fabs(A[t * n + i]) > max){
                maxIndex = t;
                max = fabs(A[t * n + i]);
            }
        }
        if (max == 0){
            // cannot operate a singular matrix
            return -1;
        }
        else{
            if (maxIndex != i){
                // save pivoting information
                temps = ipiv[i];
                ipiv[i] = ipiv[maxIndex];
                ipiv[maxIndex] = temps;
                // swap rows
                memcpy(temp, A + i * n, n * sizeof(double));
                memcpy(A + i * n, A + maxIndex * n, n * sizeof(double));
                memcpy(A + maxIndex * n, temp, n * sizeof(double));
            }
        }

        // factorization
        for (j = i + 1; j < n; j++){
            A[j * n + i] = A[j * n + i] / A[i * n + i];
            for (k = i + 1; k < n; k++){
                A[j * n + k] -= A[j  *n + i] * A[i * n + k];
            }
        }
    }
    free(temp);
    return 0;
}

/**
 * 
 * this function computes triangular matrix - vector solver
 * for a square matrix . according to lecture slides, this
 * function computes forward AND backward subtitution in the
 * same function.
 * 
 * syntax 
 *  
 *  input :
 *      UPLO  'L' or 'U' , denotes whether input matrix is upper
 *                         lower triangular . ( forward / backward
 *                         substitution )
 * 
 *      A     n by n     , square matrix
 * 
 *      B     1 by n     , vector
 * 
 *      ipiv  1 by n     , vector , denotes interchanged index due
 *                                  to pivoting by mydgetrf()
 * 
 *      n                , length of vector / size of matrix
 *  
 *  output :
 *      none
 * 
 **/
void mydtrsv(char UPLO, double *A, double *B, int n, int *ipiv)
{
    int i, j;
    double sum;
    double *y = (double*)malloc(n * sizeof(double));
    if (UPLO == 'L'){
        y[0] = B[ipiv[0]];
        for (i = 1; i < n; i++){
            sum = 0;
            for (j = 0; j < i; j++){
                sum += y[j] * A[i*n + j];
            }
            y[i] = B[ipiv[i]] - sum;
        }
    }
    else if (UPLO == 'U'){
        y[n - 1] = B[n - 1] / A[(n - 1)*n + n - 1];
        for (i = n - 2; i >= 0; i--){
            sum = 0;
            for (j = i + 1; j < n; j++){
                sum += y[j] * A[i*n + j];
            }
            y[i] = (B[i] - sum) / A[i*n + i];
        }
    }
    memcpy(B, y, sizeof(double) * n);
    free(y);
}


int get_block_size(){
    //return the block size you use in your matrix multiplication code.
    /*add your code here, 128 is an example and can be modified. */
    // With the returned block size the test code will adaptively adjust the input sizes to avoid corner cases.
    return 128;

}

//The matrix multiplication function used in blocked GEPP.
// You need to let the mydgemm adapt to non-square inputs.
void mydgemm(double *A, double *B, double *C, int n, int i, int j, int k, int b)
{
    register int i1, j1, k1;
    for (i1 = i; (i1 < i + b) && (i1 < n); i1 += 3)
        for (j1 = j; (j1 < j + b) && (j1 < n); j1 += 3) {
            int c0 = i1 * n + j1;
            int c1 = c0 + n;
            int c2 = c1 + n;
            register double c00 = C[c0];
            register double c01 = C[c0 + 1];
            register double c02 = C[c0 + 2];
            register double c10 = C[c1];
            register double c11 = C[c1 + 1];
            register double c12 = C[c1 + 2];
            register double c20 = C[c2];
            register double c21 = C[c2 + 1];
            register double c22 = C[c2 + 2];
            for (k1 = k; (k1 < k + b) && (k1 < n); k1 += 3) {
                int a0 = i1 * n + k1;
                int a1 = a0 + n;
                int a2 = a1 + n;
                int b0 = k1 * n + j1;
                int b1 = b0 + n;
                int b2 = b1 + n;
                register double a00 = A[a0];
                register double a10 = A[a1];
                register double a20 = A[a2];
                register double b00 = B[b0];
                register double b01 = B[b0 + 1];
                register double b02 = B[b0 + 2];

                c00 -= a00 * b00;
                c01 -= a00 * b01;
                c02 -= a00 * b02;
                c10 -= a10 * b00;
                c11 -= a10 * b01;
                c12 -= a10 * b02;
                c20 -= a20 * b00;
                c21 -= a20 * b01;
                c22 -= a20 * b02;

                a00 = A[a0 + 1];
                a10 = A[a1 + 1];
                a20 = A[a2 + 1];
                b00 = B[b1];
                b01 = B[b1 + 1];
                b02 = B[b1 + 2];

                c00 -= a00 * b00;
                c01 -= a00 * b01;
                c02 -= a00 * b02;
                c10 -= a10 * b00;
                c11 -= a10 * b01;
                c12 -= a10 * b02;
                c20 -= a20 * b00;
                c21 -= a20 * b01;
                c22 -= a20 * b02;

                a00 = A[a0 + 2];
                a10 = A[a1 + 2];
                a20 = A[a2 + 2];
                b00 = B[b2];
                b01 = B[b2 + 1];
                b02 = B[b2 + 2];

                c00 -= a00 * b00;
                c01 -= a00 * b01;
                c02 -= a00 * b02;
                c10 -= a10 * b00;
                c11 -= a10 * b01;
                c12 -= a10 * b02;
                c20 -= a20 * b00;
                c21 -= a20 * b01;
                c22 -= a20 * b02;

            }
            C[c0] = c00;
            C[c0 + 1] = c01;
            C[c0 + 2] = c02;
            C[c1] = c10;
            C[c1 + 1] = c11;
            C[c1 + 2] = c12;
            C[c2] = c20;
            C[c2 + 1] = c21;
            C[c2 + 2] = c22;
        }
}

/**
 * 
 * this function computes LU decomposition
 * for a square matrix using block gepp introduced in course
 * lecture .
 * 
 * just implement the block algorithm you learned in class.
 * 
 * syntax 
 *  
 *  input :
 *      
 * 
 *      A     n by n     , square matrix
 * 
 *    
 * 
 *      ipiv  1 by n     , vector , denotes interchanged index due
 *                                  to pivoting by mydgetrf()
 * 
 *      n                , length of vector / size of matrix
 *     
 *      b                , block size   
 *  output :
 *      return -1 : if the matrix A is singular (max pivot == 0)
 *      return  0 : return normally 
 * 
 **/
int mydgetrf_block(double *A, int *ipiv, int n, int b) 
{
    int ib, i, j, k, maxIndex, temps;
    double max, sum;
    double *temp = (double*)malloc(sizeof(double) * n);

    for (ib = 0; ib < (n - 1); ib += b){
        for (i = ib; i < ib + b && i < n; i++){
            maxIndex = i;
            max = fabs(A[i*n + i]);

            for (j = i + 1; j < n; j++){
                if (fabs(A[j*n + i]) > max){
                    maxIndex = j;
                    max = fabs(A[j*n + i]);
                }
            }
            if (max == 0){
                return -1;
            }
            else{
                if (maxIndex != i){
                    temps = ipiv[i];
                    ipiv[i] = ipiv[maxIndex];
                    ipiv[maxIndex] = temps;
                    memcpy(temp, A + i * n, n * sizeof(double));
                    memcpy(A + i * n, A + maxIndex * n, n * sizeof(double));
                    memcpy(A + maxIndex * n, temp, n * sizeof(double));
                }
            }

            // factorization
            for (j = i + 1; j < n; j++){
                A[j*n + i] = A[j*n + i] / A[i*n + i];
                for (k = i + 1; k < ib + b && k < n; k++){// pay attention that bound of k
                    A[j*n + k] -= A[j*n + i] * A[i*n + k];
                }
            }
        }
        // update A Computing A=LL-1 A equals solving A_new from LL*A_new = A_old.
        for (i = ib; i < ib + b && i < n; i++){
            //don't know how to do it
        }
        // calculate the lower right conner matrix use BLAS3
        for (i = ib + b; i < n; i += b){
            for (j = ib + b; j < n; j += b){
                mydgemm(A, A, A, n, i, j, ib, b);
            }
        }
    }
    return 0;
}

