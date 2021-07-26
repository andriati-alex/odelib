/**
 * \file vector.h
 * \author Alex Andriati
 * \date July 2021
 * \brief File containing suitable datatype for n-dimensional vectors
 *
 * A n-dimensional vector can be implemented as array of values
 * using pointers or struct, to also parse the vector dimension
 */

#ifndef VECTOR_H
#define VECTOR_H

#include <complex.h>
#include <stdio.h>
#include <stdlib.h>

/** \brief Array of real numbers as double pointers */
typedef double * Rarray;

/** \brief Array of complex numbers as double complex pointers */
typedef double complex * Carray;

/** \brief Struct for real vector definition with a given dimension */
typedef struct{
    unsigned int dim;   /// space dimension, also number of components
    Rarray vals;        /// array with values. Size of `dim`
} _RealVector;

/** \brief Real vector struct address **/
typedef _RealVector * RealVector;

/** \brief Struct for complex vector definition with a given dimension */
typedef struct{
    unsigned int dim;   /// space dimension, also number of components
    Carray vals;        /// array with values. Size of `dim`
} _ComplexVector;

/** \brief Complex vector struct address **/
typedef _ComplexVector * ComplexVector;


/** \brief Return fresh allocated real(double) array */
Rarray alloc_rarr(unsigned int array_size);


/** \brief Return fresh allocated complex(double) array */
Carray alloc_carr(unsigned int array_size);


/** \brief Return fresh allocated RealVector struct address */
RealVector alloc_rvec(unsigned int vec_size);


/** \brief Return fresh allocated ComplexVector struct address */
ComplexVector alloc_cvec(unsigned int vec_size);


/** \brief Free allocated RealVector struct and its internal array */
void free_rvec(RealVector vec);


/** \brief Free allocated ComplexVector struct and its internal array */
void free_cvec(ComplexVector vec);


/** \brief Copy values from the first array to the second */
void carr_copy_values(unsigned int array_size, Carray from, Carray to);


/** \brief Copy values from the first array to the second */
void rarr_copy_values(unsigned int array_size, Rarray from, Rarray to);


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
);


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
);

#endif
