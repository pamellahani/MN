#include "mnblas.h"
#include "complexe2.h"
#include <stdio.h>

/*
const int M: This parameter represents the number of rows in matrix A.

const int N: This parameter represents the number of columns in matrix A.

const float alpha: This parameter is a scalar value used to multiply matrix A and vector X.

const float *A: This parameter is a pointer to matrix A, which is of size M x N (if row major) or N x M (if column major).

const int lda: This parameter is the leading dimension of matrix A. It represents the number of elements between the start of consecutive rows (if row major) or columns (if column major).

const float *X: This parameter is a pointer to the input vector X of size N (if row major) or M (if column major).

const int incX: This parameter is the increment for elements in vector X.

const float beta: This parameter is a scalar value used to multiply vector Y.

float *Y: This parameter is a pointer to the output vector Y of size M (if row major) or N (if column major).

const int incY: This parameter is the increment for elements in vector Y.


The function performs the operation Y = alpha * A * X + beta * Y in single precision.
*/

void mncblas_sgemv(const MNCBLAS_LAYOUT layout,
                   const MNCBLAS_TRANSPOSE TransA, const int M, const int N,
                   const float alpha, const float *A, const int lda,
                   const float *X, const int incX, const float beta,
                   float *Y, const int incY)
{

    int i, j;
    float temp;

    int num_threads = 4;
    omp_set_num_threads(num_threads);

    double start_time = omp_get_wtime();

#pragma omp parallel {}

    double end_time = omp_get_wtime();
    double elapsed_time = end_time - start_time;

    printf("Elapsed time: %f seconds\n", elapsed_time);

    if (TransA == MNCblasNoTrans)
    {
        for (i = 0; i < M; i++)
        {
            temp = 0.0;
            for (j = 0; j < N; j++)
            {
                temp += alpha * A[i * lda + j] * X[j * incX];
            }
            Y[i * incY] = beta * Y[i * incY] + temp;
        }
    }
    else if (TransA == MNCblasTrans)
    {
        for (i = 0; i < N; i++)
        {
            temp = 0.0;
            for (j = 0; j < M; j++)
            {
                temp += alpha * A[i * lda + j] * X[j * incX];
            }
            Y[i * incY] = beta * Y[i * incY] + temp;
        }
    }
}

void mncblas_dgemv(MNCBLAS_LAYOUT layout,
                   MNCBLAS_TRANSPOSE TransA, const int M, const int N,
                   const double alpha, const double *A, const int lda,
                   const double *X, const int incX, const double beta,
                   double *Y, const int incY)
{

    int i, j;
    double temp;

    if (TransA == MNCblasNoTrans)
    {
        for (i = 0; i < M; i++)
        {
            temp = 0.0;
            for (j = 0; j < N; j++)
            {
                temp += alpha * A[i * lda + j] * X[j * incX];
            }
            Y[i * incY] = beta * Y[i * incY] + temp;
        }
    }
    else if (TransA == MNCblasTrans)
    {
        for (i = 0; i < N; i++)
        {
            temp = 0.0;
            for (j = 0; j < M; j++)
            {
                temp += alpha * A[i * lda + j] * X[j * incX];
            }
            Y[i * incY] = beta * Y[i * incY] + temp;
        }
    }
}

void mncblas_cgemv(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA, const int M, const int N,
                   const void *alpha, const void *A, const int lda, const void *X, const int incX,
                   const void *beta, void *Y, const int incY)
{
    if (layout == MNCblasRowMajor)
    {
        if (TransA == MNCblasNoTrans)
        {
            int i, j, ix, iy;
            float alpha_real = ((complexe_float_t *)alpha)->real;
            float alpha_imag = ((complexe_float_t *)alpha)->imaginary;
            float beta_real = ((complexe_float_t *)beta)->real;
            float beta_imag = ((complexe_float_t *)beta)->imaginary;
            float *X_real = (float *)X;
            float *X_imag = X_real + 1;
            float *Y_real = (float *)Y;
            float *Y_imag = Y_real + 1;
            float *A_real = (float *)A;
            float *A_imag = A_real + 1;

            for (i = 0; i < M; i++)
            {
                ix = 0;
                iy = 0;

                Y_real[i * incY] *= beta_real;
                Y_imag[i * incY] *= beta_imag;

                for (j = 0; j < N; j++)
                {
                    float Aij_real = A_real[i * lda + j * 2];
                    float Aij_imag = A_imag[i * lda + j * 2];
                    float Xj_real = X_real[j * incX];
                    float Xj_imag = X_imag[j * incX];

                    Y_real[i * incY] += alpha_real * (Aij_real * Xj_real - Aij_imag * Xj_imag);
                    Y_imag[i * incY] += alpha_imag * (Aij_imag * Xj_real + Aij_real * Xj_imag);
                }
            }
        }
    }
}

void mncblas_zgemv(MNCBLAS_LAYOUT layout,
                   MNCBLAS_TRANSPOSE TransA, const int M, const int N,
                   const void *alpha, const void *A, const int lda,
                   const void *X, const int incX, const void *beta,
                   void *Y, const int incY)
{

    if (layout == MNCblasRowMajor)
    {
        if (TransA == MNCblasNoTrans)
        {
            int i, j, ix, iy;
            double alpha_real = ((complexe_double_t *)alpha)->real;
            double alpha_imag = ((complexe_double_t *)alpha)->imaginary;
            double beta_real = ((complexe_double_t *)beta)->real;
            double beta_imag = ((complexe_double_t *)beta)->imaginary;
            double *X_real = (float *)X;
            double *X_imag = X_real + 1;
            double *Y_real = (float *)Y;
            double *Y_imag = Y_real + 1;
            double *A_real = (float *)A;
            double *A_imag = A_real + 1;

            for (i = 0; i < M; i++)
            {
                ix = 0;
                iy = 0;

                Y_real[i * incY] *= beta_real;
                Y_imag[i * incY] *= beta_imag;

                for (j = 0; j < N; j++)
                {
                    double Aij_real = A_real[i * lda + j * 2];
                    double Aij_imag = A_imag[i * lda + j * 2];
                    double Xj_real = X_real[j * incX];
                    double Xj_imag = X_imag[j * incX];

                    Y_real[i * incY] += alpha_real * (Aij_real * Xj_real - Aij_imag * Xj_imag);
                    Y_imag[i * incY] += alpha_imag * (Aij_imag * Xj_real + Aij_real * Xj_imag);
                }
            }
        }
    }
}