#include "mnblas.h"
#include "complexe.h"
#include "math.h" 

/*
 * Compute the Euclidean norm of a float vector.
 */
float MNCBLAS_SNRM2(const int N, const float *X, const int incX) {
    float norm = 0.0f, scale = 0.0f;
    for (int i = 0; i < N * incX; i += incX) {
        scale = fmaxf(scale, fabsf(X[i]));
    }
    if (scale == 0.0f) {
        return 0.0f;
    }
    for (int i = 0; i < N * incX; i += incX) {
        norm += powf(X[i] / scale, 2);
    }
    return scale * sqrtf(norm);
}

/*
 * Compute the Euclidean norm of a double vector.
 */
double MNCBLAS_DNRM2(const int N, const double *X, const int incX) {
    double norm = 0.0, scale = 0.0;
    for (int i = 0; i < N * incX; i += incX) {
        scale = fmax(scale, fabs(X[i]));
    }
    if (scale == 0.0) {
        return 0.0;
    }
    for (int i = 0; i < N * incX; i += incX) {
        norm += pow(X[i] / scale, 2);
    }
    return scale * sqrt(norm);
}

/*
 * Compute the Euclidean norm of a complex float vector.
 */
float mnblas_scnrm2(const int N, const void *X, const int incX) {
    float norm = 0.0f, scale = 0.0f;
    const complexe_float_t *cx = X;
    for (int i = 0; i < N * incX; i += incX) {
        scale = fmaxf(scale, fabsf(cx[i].real));
        scale = fmaxf(scale, fabsf(cx[i].imaginary));
    }
    if (scale == 0.0f) {
        return 0.0f;
    }
    for (int i = 0; i < N * incX; i += incX) {
        norm += powf(cx[i].real / scale, 2); 
        norm += powf(cx[i].imaginary / scale, 2);
    }
    return scale * sqrtf(norm);
}

/*
 * Compute the Euclidean norm of a complex double vector.
 */
double mnblas_dznrm2(const int N, const void *X, const int incX) {
    double norm = 0.0, scale = 0.0;
    const complexe_double_t *zx = X;
    for (int i = 0; i < N * incX; i += incX) {
        scale = fmax(scale, fabs(zx[i].real));
        scale = fmax(scale, fabs(zx[i].imaginary));
    }
    if (scale == 0.0) {
        return 0.0;
    }
    for (int i = 0; i < N * incX; i += incX) {
        norm += pow(zx[i].real / scale, 2); 
        norm += pow(zx[i].imaginary / scale, 2);
    }
    return scale * sqrt(norm);
}