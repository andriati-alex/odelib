/**
 * \file vector.h
 * \author Alex Andriati
 * \date July 2021
 * \brief File containing suitable datatype for vectors as arrays
 *
 * A n-dimensional vector can be implemented as array of values
 * using pointers to real or complex numbers, here as double or
 * double complex pointer
 *
 * \note Developers whishing to contrinute see `arrays_assistant.h`
 *       which provide `static` implementation of basic functions
 */

#ifndef ARRAYS_H
#define ARRAYS_H

#include <complex.h>

/** \brief Array of real numbers as double pointers */
typedef double * Rarray;

/** \brief Array of complex numbers as double complex pointers */
typedef double complex * Carray;


#endif
