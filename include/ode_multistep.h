/**
 * \file ode_multistep.h
 * \author Alex Andriati
 * \brief ODE integration routines with multi step methods
 *
 * These numerical routines require mutiple known steps simultaneously
 * to compute a new one. Generally these methods depend on single step
 * methods to initialize the first few known steps
 */

#ifndef ODE_MULTISTEP_H
#define ODE_MULTISTEP_H

#include "derivative_signature.h"

/** \brief Struct to provide workspace for multistep methods */
typedef struct{
    int
        ms_order,       /// number of previous steps required
        system_size;    /// number of equations in ODE system
    Carray
        prev_der;       /// Hold all required previous derivatives
} _ComplexWorkspaceMS;

/** \brief Struct address with working array for multistep methods */
typedef _ComplexWorkspaceMS * ComplexWorkspaceMS;

/** \brief Struct to provide workspace for multistep methods */
typedef struct{
    int
        ms_order,       /// number of previous steps required
        system_size;    /// number of equations in ODE system
    Rarray
        prev_der;       /// Hold all required previous derivatives
} _RealWorkspaceMS;

/** \brief Struct address with working array for multistep methods */
typedef _RealWorkspaceMS * RealWorkspaceMS;


/**
 * \brief General multistep basic operation
 *
 * Multistep methods have a general form separating functions to the
 * left-hand-side and their derivatives to the right-hand-side using
 * coefficients as
 * `a[0] * y_j+1 + a[1] * y_j + ... + a[m] * y_j+1-m =
 *          h * (b[0] * y'_j+1 + b[1] * y'_j + ... + b[m] * y'_j+1-m)`
 *
 * \param 1 : grid spacing `h`
 * \param 2 : grid point correspongind to function values `x`
 * \param 3 : function pointer to routine that compute derivatives
 * \param 4 : extra arguments required in `cmplx_sys_der` function
 * \param 5 : Workspace struct address to avoid memory allocation.
 *            The array within this struct must have at least size
 *            `m * n` if the method is explicit, and `(m + 1) * n`
 *            if the method implicit. In both cases, it shall have
 *            derivative of previous steps concatenated as
 *            `[y'_j y'_j-1 ...  y'_j+1-m]`
 * \param 6 : Concatenated function steps `y_j` of all previous steps
 *            required depending on the multistep order. Within system
 *            of size `n` and multistep order `m`, the array must have
 *            size `m * n` with `[y_j y_j-1 ...  y_j+1-m]` concatenated
 * \param 7 : Function weights `a` as array of `m + 1` elements. It is
 *            used in left-hand-side of the method, ignoring a[0] given
 *            by `y_j+1 + a[1] * y_j + ... + a[m] * y_j+1-m`
 * \param 8 : Derivative weights `b` as array of `m + 1` elements. Used
 *            in right-hand-side `h * (b[0] * y'_j+1 + b[1] * y'_j + ...
 *            b[m] * y'_j+1-m`. `b[0]` will be used depending on param 9
 * \param 9 : Number of iterations for implicit method if greater than 0
 *            If zero ignore first element of param 8 and apply explicit
 *            scheme, equivalent to using zero and set this param to 1
 * \param 10: (OUTPUT) solution at next grid step
 *            (INPUT)  use as predictor if parameter 9 is greater than 0
 */
void
cmplx_general_multistep
(
        double,
        double,
        cmplx_sys_der,
        void *,
        ComplexWorkspaceMS,
        Carray,
        Rarray,
        Rarray,
        unsigned int,
        Carray
);


/**
 * \brief General multistep basic operation
 *
 * Multistep methods have a general form separating functions to the
 * left-hand-side and their derivatives to the right-hand-side using
 * coefficients as
 * `a[0] * y_j+1 + a[1] * y_j + ... + a[m] * y_j+1-m =
 *          h * (b[0] * y'_j+1 + b[1] * y'_j + ... + b[m] * y'_j+1-m)`
 *
 * \param 1 : grid spacing `h`
 * \param 2 : grid point correspongind to function values `x`
 * \param 3 : function pointer to routine that compute derivatives
 * \param 4 : extra arguments required in `cmplx_sys_der` function
 * \param 5 : Workspace struct address to avoid memory allocation.
 *            The array within this struct must have at least size
 *            `m * n` if the method is explicit, and `(m + 1) * n`
 *            if the method implicit. In both cases, it shall have
 *            derivative of previous steps concatenated as
 *            `[y'_j y'_j-1 ...  y'_j+1-m]`
 * \param 6 : Concatenated function steps `y_j` of all previous steps
 *            required depending on the multistep order. Within system
 *            of size `n` and multistep order `m`, the array must have
 *            size `m * n` with `[y_j y_j-1 ...  y_j+1-m]` concatenated
 * \param 7 : Function weights `a` as array of `m + 1` elements. It is
 *            used in left-hand-side of the method, ignoring a[0] given
 *            by `y_j+1 + a[1] * y_j + ... + a[m] * y_j+1-m`
 * \param 8 : Derivative weights `b` as array of `m + 1` elements. Used
 *            in right-hand-side `h * (b[0] * y'_j+1 + b[1] * y'_j + ...
 *            b[m] * y'_j+1-m`. `b[0]` will be used depending on param 9
 * \param 9 : Number of iterations for implicit method if greater than 0
 *            If zero ignore first element of param 8 and apply explicit
 *            scheme, equivalent to using zero and set this param to 1
 * \param 10: (OUTPUT) solution at next grid step
 *            (INPUT)  use as predictor if parameter 9 is greater than 0
 */
void
real_general_multistep
(
        double,
        double,
        real_sys_der,
        void *,
        RealWorkspaceMS,
        Rarray,
        Rarray,
        Rarray,
        unsigned int,
        Rarray
);


#endif
