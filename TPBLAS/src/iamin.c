#include "mnblas.h"
#include "complexe.h"
#include <stdio.h>

float fabsf(float x)
{
    if (x < 0.0f)
    {
        return -x;
    }
    else
    {
        return x;
    }
}

double fabs(double x)
{
    if (x < 0.0)
    {
        return -x;
    }
    else
    {
        return x;
    }
}

// Finds the index of the element with the smallest absolute value.

CBLAS_INDEX cblas_isamin(const int n, const float *x, const int incx)
{
    CBLAS_INDEX min_index = 0;
    float min_value = fabsf(x[0]);
    int i;
    for (i = 1; i < n; i += incx)
    {
        float abs_value = fabsf(x[i]);
        if (abs_value < min_value)
        {
            min_value = abs_value;
            min_index = i;
        }
    }
    return min_index;
}

CBLAS_INDEX cblas_idamin(const int n, const double *x, const int incx)
{
    CBLAS_INDEX min_index = 0;
    double min_value = fabs(x[0]);
    int i;
    for (i = 1; i < n; i+=incx)
    {
        double abs_value = fabs(x[i]);
        if (abs_value < min_value)
        {
            min_value = abs_value;
            min_index = i;
        }
    }
    return min_index;
}

CBLAS_INDEX cblas_icamin(const int n, const void *x, const int incx)
{
    CBLAS_INDEX min_index = 0;
    float min_value = fabsf(((complexe_float_t *)x)[0].real) + fabsf(((complexe_float_t *)x)[0].imaginary);
    int i;
    for (i = 1; i < n; i+=incx)
    {
        float abs_value = fabsf(((complexe_float_t *)x)[i].real) + fabsf(((complexe_float_t *)x)[i].imaginary);
        if (abs_value < min_value)
        {
            min_value = abs_value;
            min_index = i;
        }
    }
    return min_index;
}

CBLAS_INDEX cblas_izamin(const int n, const void *x, const int incx)
{
    CBLAS_INDEX min_index = 0;
    double min_value = fabs(((complexe_double_t *)x)[0].real) + fabs(((complexe_double_t *)x)[0].imaginary);
    int i;
    for (i = 1; i < n; i+=incx)
    {
        double abs_value = fabs(((complexe_double_t *)x)[i].real) + fabs(((complexe_double_t *)x)[i].imaginary);
        if (abs_value < min_value)
        {
            min_value = abs_value;
            min_index = i;
        }
    }
    return min_index;
}