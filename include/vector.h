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

/** \brief Array of real numbers as double pointers */
typedef double * Rarray;

/** \brief Array of complex numbers as double complex pointers */
typedef double complex * Carray;

/** \brief Struct for real vector definition with a given dimension */
typedef struct{
    unsigned int
        dim;    /// space dimension, also number of components
    Rarray
        vals;   /// array with values. Size of `dim`
} RealVector;

/** \brief Struct for complex vector definition with a given dimension */
typedef struct{
    unsigned int
        dim;    /// space dimension, also number of components
    Carray
        vals;   /// array with values. Size of `dim`
} ComplexVector;

#endif
