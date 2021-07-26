/**
 * \file derivative_signature.h
 * \author Alex Andriati
 * \date July 2021
 * \brief Function signatures to evaluate derivatives of ODE systems
 *
 * An Ordinary Differential Equation(ODE) system is defined by a set
 * of coupled first order differential equations. Such system always
 * can be made linear in function derivatives, where derivatives can
 * be written as y' = f(x, y) with y a vector of functions values at
 * grid point x. These signatures standardize user input function as
 * well provide general parameters struct for extra parameters
 */

#ifndef DERIVATIVE_SIGNATURE_H
#define DERIVATIVE_SIGNATURE_H

#include "vector.h"

/** \brief Struct with input parameters for derivatives computation */
typedef struct{
    unsigned int system_size;   /// number of equations in the system
    double x;                   /// grid point of the known solution
    Rarray y;
    void * extra_args;          /// user-defined external arguments
} _RealODEInputParameters;

/** \brief Input parameters struct address needed in function signature */
typedef _RealODEInputParameters * RealODEInputParameters;

/** \brief Struct with input parameters for derivatives computation */
typedef struct{
    unsigned int system_size;   /// number of equations in the system
    double x;                   /// grid point of the known solution
    Carray y;
    void * extra_args;          /// user-defined external arguments
} _ComplexODEInputParameters;

/** \brief Input parameters struct address needed in function signature */
typedef _ComplexODEInputParameters * ComplexODEInputParameters;

/**
 * \brief Function signature to compute derivatives of real ODE system
 *
 * \param 1 : optional arguments
 * \param 2 : vector with all functions evaluated at current grid point
 * \param 3 : (OUTPUT) function derivatives at current grid point
 */
typedef void (*real_sys_der)(RealODEInputParameters, Rarray);

/**
 * \brief Function signature to compute derivatives of complex ODE system
 *
 * \param 1 : optional arguments
 * \param 2 : vector with all functions evaluated at current grid point
 * \param 3 : (OUTPUT) function derivatives at current grid point
 */
typedef void (*cmplx_sys_der)(ComplexODEInputParameters, Carray);


#endif
