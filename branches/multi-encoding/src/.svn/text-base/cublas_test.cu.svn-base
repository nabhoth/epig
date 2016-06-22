/* Includes, cuda */
#include "cublas.h"



/* Matrix size */
#define N  (1275)

/* Main */
int main(int argc, char** argv)
{
    float* h_A;
    float* h_B;
    float* h_C;
    float* d_A = 0;
    float* d_B = 0;
    float* d_C = 0;
    float alpha = 1.0f;
    float beta = 0.0f;
    int n2 = N * N;
    int i;

    /* Initialize CUBLAS */
    cublasInit();

    /* Allocate host memory for the matrices */
    h_A = (float*)malloc(n2 * sizeof(h_A[0]));
    h_B = (float*)malloc(n2 * sizeof(h_B[0]));
    h_C = (float*)malloc(n2 * sizeof(h_C[0]));

    /* Fill the matrices with test data */
    for (i = 0; i < n2; i++)
    {
        h_A[i] = rand() / (float)RAND_MAX;
        h_B[i] = rand() / (float)RAND_MAX;
        h_C[i] = rand() / (float)RAND_MAX;
    }

    /* Allocate device memory for the matrices */
    cublasAlloc(n2, sizeof(d_A[0]), (void**)&d_A);
    cublasAlloc(n2, sizeof(d_B[0]), (void**)&d_B);
    cublasAlloc(n2, sizeof(d_C[0]), (void**)&d_C);

    /* Initialize the device matrices with the host matrices */
    cublasSetVector(n2, sizeof(h_A[0]), h_A, 1, d_A, 1);
    cublasSetVector(n2, sizeof(h_B[0]), h_B, 1, d_B, 1);
    cublasSetVector(n2, sizeof(h_C[0]), h_C, 1, d_C, 1);

    /* Performs operation using cublas */
    cublasSgemm('n', 'n', N, N, N, alpha,
                d_A, N, d_B, N, beta, d_C, N);

    /* Read the result back */
    cublasGetVector(n2, sizeof(h_C[0]), d_C, 1, h_C, 1);

    /* Memory clean up */
    free(h_A);
    free(h_B);
    free(h_C);
    cublasFree(d_A);
    cublasFree(d_B);
    cublasFree(d_C);

    /* Shutdown */
    cublasShutdown();

    return EXIT_SUCCESS;
}



