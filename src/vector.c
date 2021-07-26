#include "vector.h"
#include <stdlib.h>
#include <stdio.h>

/** \brief Return fresh allocated real(double) array */
Rarray
alloc_rarr(unsigned int n)
{
    Rarray ptr = (Rarray) malloc(n * sizeof(double));
    if (ptr == NULL)
    {
        printf("\n\nProblem in Rarray allocation\n\n");
        exit(EXIT_FAILURE);
    }
    return ptr;
}


/** \brief Return fresh allocated complex(double) array */
Carray
alloc_carr(unsigned int n)
{
    Carray ptr = (Carray) malloc(n * sizeof(double complex));
    if (ptr == NULL)
    {
        printf("\n\nProblem in Carray allocation\n\n");
        exit(EXIT_FAILURE);
    }
    return ptr;
}


/** \brief Return fresh allocated RealVector struct address */
RealVector
alloc_rvec(unsigned int n)
{
    RealVector ptr = (RealVector) malloc(sizeof(_RealVector));
    if (ptr == NULL)
    {
        printf("\n\nProblem in RealVector allocation\n\n");
        exit(EXIT_FAILURE);
    }
    ptr->dim = n;
    ptr->vals = alloc_rarr(n);
    return ptr;
}


/** \brief Return fresh allocated ComplexVector struct address */
ComplexVector
alloc_cvec(unsigned int n)
{
    ComplexVector ptr = (ComplexVector) malloc(sizeof(_ComplexVector));
    if (ptr == NULL)
    {
        printf("\n\nProblem in ComplexVector allocation\n\n");
        exit(EXIT_FAILURE);
    }
    ptr->dim = n;
    ptr->vals = alloc_carr(n);
    return ptr;
}


/** \brief Free allocated RealVector struct and its internal array */
void
free_rvec(RealVector vec)
{
    free(vec->vals);
    free(vec);
}


/** \brief Free allocated ComplexVector struct and its internal array */
void
free_cvec(ComplexVector vec)
{
    free(vec->vals);
    free(vec);
}


/** \brief Copy values from the first array to the second */
void
carr_copy_values(unsigned int n, Carray from, Carray to)
{
    for (unsigned int i = 0; i < n; i++) to[i] = from[i];
}


/** \brief Copy values from the first array to the second */
void
rarr_copy_values(unsigned int n, Rarray from, Rarray to)
{
    for (unsigned int i = 0; i < n; i++) to[i] = from[i];
}


/** \brief Linear combination of two arrays with a constant scalar
 *
 * Set in out array the values : a0 + a1 * inp1 + a2 * inp2
 */
void
carr_linear_comb(
        unsigned int n,
        double complex a0,
        double complex a1,
        double complex a2,
        Carray inp1,
        Carray inp2,
        Carray out
)
{
    for (unsigned int i = 0; i < n; i++)
    {
        out[i] = a0 + a1 * inp1[i] + a2 * inp2[i];
    }
}


/** \brief Linear combination of two arrays with a constant scalar
 *
 * Set in out array the values : a0 + a1 * inp1 + a2 * inp2
 */
void
rarr_linear_comb(
        unsigned int n,
        double a0,
        double a1,
        double a2,
        Rarray inp1,
        Rarray inp2,
        Rarray out
)
{
    for (unsigned int i = 0; i < n; i++)
    {
        out[i] = a0 + a1 * inp1[i] + a2 * inp2[i];
    }
}
