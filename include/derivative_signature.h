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
 * well provide general void pointer for extra optional arguments
 */

#ifndef DERIVATIVE_SIGNATURE_H
#define DERIVATIVE_SIGNATURE_H

#include "vector.h"

/**
 * \brief Compute derivative of real ODE system
 *
 * \param 1 : size of the system
 * \param 2 : grid point `x`
 * \param 3 : vector with all functions computed at grid point `x`
 * \param 4 : (OUTPUT) function derivatives at grid point `x`
 * \param 5 : optional arguments
 */
typedef void (*real_sys_der)(int, double, Rarray, Rarray, void *);

/**
 * \brief Compute derivative of complex ODE system
 *
 * \param 1 : size of the system
 * \param 2 : grid point `x`
 * \param 3 : vector with all functions computed at grid point `x`
 * \param 4 : (OUTPUT) function derivatives at grid point `x`
 * \param 5 : optional arguments
 */
typedef void (*cmplx_sys_der)(int, double, Carray, Carray, void *);

#endif
