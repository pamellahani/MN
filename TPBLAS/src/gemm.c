#include "mnblas.h"
#include "complexe2.h"
#include <stdio.h>

/*
The matrix operation performed is C = alpha * A * B + beta * C, where the multiplication is 
performed taking into account the transpose flags of A and B.
*/

void mncblas_sgemm(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA,
                   MNCBLAS_TRANSPOSE TransB, const int M, const int N,
                   const int K, const float alpha, const float *A,
                   const int lda, const float *B, const int ldb,
                   const float beta, float *C, const int ldc)
{
    int i, j, l;
    float temp;

    if (layout == MNCblasRowMajor)
    {
        if (TransA == MNCblasNoTrans && TransB == MNCblasNoTrans)
        {
            for (i = 0; i < M; i++)
            {
                for (j = 0; j < N; j++)
                {
                    temp = 0;
                    for (l = 0; l < K; l++)
                    {
                        temp += A[i * lda + l] * B[l * ldb + j];
                    }
                    C[i * ldc + j] = alpha * temp + beta * C[i * ldc + j];
                }
            }
        }
        else if (TransA == MNCblasNoTrans && TransB == MNCblasTrans)
        {
            for (i = 0; i < M; i++)
            {
                for (j = 0; j < N; j++)
                {
                    temp = 0;
                    for (l = 0; l < K; l++)
                    {
                        temp += A[i * lda + l] * B[j * ldb + l];
                    }
                    C[i * ldc + j] = alpha * temp + beta * C[i * ldc + j];
                }
            }
        }
        else if (TransA == MNCblasTrans && TransB == MNCblasNoTrans)
        {
            for (i = 0; i < M; i++)
            {
                for (j = 0; j < N; j++)
                {
                    temp = 0;
                    for (l = 0; l < K; l++)
                    {
                        temp += A[l * lda + i] * B[l * ldb + j];
                    }
                    C[i * ldc + j] = alpha * temp + beta * C[i * ldc + j];
                }
            }
        }
        else // TransA == MNCblasTrans && TransB == MNCblasTrans
        {
            for (i = 0; i < M; i++)
            {
                for (j = 0; j < N; j++)
                {
                    temp = 0;
                    for (l = 0; l < K; l++)
                    {
                        temp += A[l * lda + i] * B[j * ldb + l];
                    }
                    C[i * ldc + j] = alpha * temp + beta * C[i * ldc + j];
                }
            }
        }
    }
}

void mncblas_dgemm(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA,
                   MNCBLAS_TRANSPOSE TransB, const int M, const int N,
                   const int K, const double alpha, const double *A,
                   const int lda, const double *B, const int ldb,
                   const double beta, double *C, const int ldc){
                    {
    int i, j, l;
    double temp;

    if (layout == MNCblasRowMajor)
    {
        if (TransA == MNCblasNoTrans && TransB == MNCblasNoTrans)
        {
            for (i = 0; i < M; i++)
            {
                for (j = 0; j < N; j++)
                {
                    temp = 0;
                    for (l = 0; l < K; l++)
                    {
                        temp += A[i * lda + l] * B[l * ldb + j];
                    }
                    C[i * ldc + j] = alpha * temp + beta * C[i * ldc + j];
                }
            }
        }
        else if (TransA == MNCblasNoTrans && TransB == MNCblasTrans)
        {
            for (i = 0; i < M; i++)
            {
                for (j = 0; j < N; j++)
                {
                    temp = 0;
                    for (l = 0; l < K; l++)
                    {
                        temp += A[i * lda + l] * B[j * ldb + l];
                    }
                    C[i * ldc + j] = alpha * temp + beta * C[i * ldc + j];
                }
            }
        }
        else if (TransA == MNCblasTrans && TransB == MNCblasNoTrans)
        {
            for (i = 0; i < M; i++)
            {
                for (j = 0; j < N; j++)
                {
                    temp = 0;
                    for (l = 0; l < K; l++)
                    {
                        temp += A[l * lda + i] * B[l * ldb + j];
                    }
                    C[i * ldc + j] = alpha * temp + beta * C[i * ldc + j];
                }
            }
        }
        else // TransA == MNCblasTrans && TransB == MNCblasTrans
        {
            for (i = 0; i < M; i++)
            {
                for (j = 0; j < N; j++)
                {
                    temp = 0;
                    for (l = 0; l < K; l++)
                    {
                        temp += A[l * lda + i] * B[j * ldb + l];
                    }
                    C[i * ldc + j] = alpha * temp + beta * C[i * ldc + j];
                }
            }
        }
    }
}
                   }

void mncblas_cgemm(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA,
                   MNCBLAS_TRANSPOSE TransB, const int M, const int N,
                   const int K, const void *alpha, const void *A,
                   const int lda, const void *B, const int ldb,
                   const void *beta, void *C, const int ldc){
    int i, j, k;
    float temp;

    complexe_float_t *a = (complexe_float_t*)A; 
    complexe_float_t *b = (complexe_float_t*)B; 
    complexe_float_t *c = (complexe_float_t*)C; 
    complexe_float_t *alpha_c = (complexe_float_t*)alpha; 
    complexe_float_t *beta_c = (complexe_float_t*)beta; 
    complexe_float_t res; 

if (TransA == MNCblasNoTrans){
        if (TransB == MNCblasNoTrans){
            for(i = 0 ; i < M ; i++){
                for(j = 0 ; j < N ; j++){
                    res.real = 0.0f;
                    res.imaginary = 0.0f;
                    for(k = 0 ; k < K ; k++){
                        res = add_complexe_float(mult_complexe_float(a[k*M+i], b[j+k*N]), res);
                    }
                    c[j*M+i].real = beta_c->real * c[j*M+i].real + alpha_c->real * res.real - alpha_c->imaginary* res.imaginary;
                    c[j*M+i].imaginary = beta_c->real * c[j*M+i].imaginary + alpha_c->real * res.imaginary + alpha_c->imaginary* res.real;
                }
            }
        }else{
            for(i = 0 ; i < M ; i++){
                for(j = 0 ; j < N ; j++){
                    res.real = 0.0f;
                    res.imaginary = 0.0f;
                    for(k = 0 ; k < K ; k++){
                        res = add_complexe_float(mult_complexe_float(a[k*M+i], b[k+j*K]), res);
                    }
                    c[j*M+i].real = beta_c->real * c[j*M+i].real + alpha_c->real * res.real - alpha_c->imaginary* res.imaginary;
                    c[j*M+i].imaginary = beta_c->real * c[j*M+i].imaginary + alpha_c->real * res.imaginary + alpha_c->imaginary* res.real;
                }
            }
        }
    }else{
        if (TransB == MNCblasNoTrans){
            for(i = 0 ; i < M ; i++){
                for(j = 0 ; j < N ; j++){
                    res.real = 0.0f;
                    res.imaginary = 0.0f;
                    for(k = 0 ; k < K ; k++){
                        res = add_complexe_float(mult_complexe_float(a[i+k*M], b[j+k*N]), res);
                    }
                    c[j*M+i].real = beta_c->real * c[j*M+i].real + alpha_c->real * res.real - alpha_c->imaginary* res.imaginary;
                    c[j*M+i].imaginary = beta_c->real * c[j*M+i].imaginary + alpha_c->real * res.imaginary + alpha_c->imaginary* res.real;
                }
            }
        }else{
            for(i = 0 ; i < M ; i++){
                for(j = 0 ; j < N ; j++){
                    res.real = 0.0f;
                    res.imaginary = 0.0f;
                    for(k = 0 ; k < K ; k++){
                        res = add_complexe_float(mult_complexe_float(a[i+k*M], b[k+j*K]), res);
                    }
                    c[j*M+i].real = beta_c->real * c[j*M+i].real + alpha_c->real * res.real - alpha_c->imaginary* res.imaginary;
                    c[j*M+i].imaginary = beta_c->real * c[j*M+i].imaginary + alpha_c->real * res.imaginary + alpha_c->imaginary* res.real;
                }
            }
        }
    }
}


void mncblas_zgemm(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA,
                   MNCBLAS_TRANSPOSE TransB, const int M, const int N,
                   const int K, const void *alpha, const void *A,
                   const int lda, const void *B, const int ldb,
                   const void *beta, void *C, const int ldc){
    int i, j, k;
    double temp;

    complexe_double_t *a = (complexe_double_t*)A; 
    complexe_double_t *b = (complexe_double_t*)B; 
    complexe_double_t *c = (complexe_double_t*)C; 
    complexe_double_t *alpha_c = (complexe_double_t*)alpha; 
    complexe_double_t *beta_c = (complexe_double_t*)beta; 
    complexe_double_t res; 

if (TransA == MNCblasNoTrans){
        if (TransB == MNCblasNoTrans){
            for(i = 0 ; i < M ; i++){
                for(j = 0 ; j < N ; j++){
                    res.real = 0.0f;
                    res.imaginary = 0.0f;
                    for(k = 0 ; k < K ; k++){
                        res = add_complexe_double(mult_complexe_double(a[k*M+i], b[j+k*N]), res);
                    }
                    c[j*M+i].real = beta_c->real * c[j*M+i].real + alpha_c->real * res.real - alpha_c->imaginary* res.imaginary;
                    c[j*M+i].imaginary = beta_c->real * c[j*M+i].imaginary + alpha_c->real * res.imaginary + alpha_c->imaginary* res.real;
                }
            }
        }else{
            for(i = 0 ; i < M ; i++){
                for(j = 0 ; j < N ; j++){
                    res.real = 0.0f;
                    res.imaginary = 0.0f;
                    for(k = 0 ; k < K ; k++){
                        res = add_complexe_double(mult_complexe_double(a[k*M+i], b[k+j*K]), res);
                    }
                    c[j*M+i].real = beta_c->real * c[j*M+i].real + alpha_c->real * res.real - alpha_c->imaginary* res.imaginary;
                    c[j*M+i].imaginary = beta_c->real * c[j*M+i].imaginary + alpha_c->real * res.imaginary + alpha_c->imaginary* res.real;
                }
            }
        }
    }
    else{
        if (TransB == MNCblasNoTrans){
            for(i = 0 ; i < M ; i++){
                for(j = 0 ; j < N ; j++){
                    res.real = 0.0f;
                    res.imaginary = 0.0f;
                    for(k = 0 ; k < K ; k++){
                        res = add_complexe_double(mult_complexe_double(a[i+k*M], b[j+k*N]), res);
                    }
                    c[j*M+i].real = beta_c->real * c[j*M+i].real + alpha_c->real * res.real - alpha_c->imaginary* res.imaginary;
                    c[j*M+i].imaginary = beta_c->real * c[j*M+i].imaginary + alpha_c->real * res.imaginary + alpha_c->imaginary* res.real;
                }
            }
        }else{
            for(i = 0 ; i < M ; i++){
                for(j = 0 ; j < N ; j++){
                    res.real = 0.0f;
                    res.imaginary = 0.0f;
                    for(k = 0 ; k < K ; k++){
                        res = add_complexe_double(mult_complexe_double(a[i+k*M], b[k+j*K]), res);
                    }
                    c[j*M+i].real = beta_c->real * c[j*M+i].real + alpha_c->real * res.real - alpha_c->imaginary* res.imaginary;
                    c[j*M+i].imaginary = beta_c->real * c[j*M+i].imaginary + alpha_c->real * res.imaginary + alpha_c->imaginary* res.real;
                }
            }
        }
    }
}
