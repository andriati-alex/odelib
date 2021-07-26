/**
 * \file ode_multistep.h
 * \author Alex Andriati
 * \brief ODE integration routines with multistep methods
 *
 * These numerical routines require mutiple known steps simultaneously
 * to compute a new one. Generally these methods depend on single step
 * methods to initialize the first few known steps
 */

#ifndef ODE_MULTISTEP_H
#define ODE_MULTISTEP_H

#include "derivative_signature.h"

/** \brief Struct to provide complex workspace for multistep methods
 *
 * Provide basic data to simplify the general multistep methods API,
 * with max number of previous steps required, system size and array
 * with all previous derivatives pre-computed
 */
typedef struct{
    int
        ms_order,       /// number of previous steps required
        system_size;    /// number of equations in ODE system
    Carray
        prev_der;       /// Hold all required previous derivatives
} _ComplexWorkspaceMS;

/** \brief Workspace struct address for multistep methods */
typedef _ComplexWorkspaceMS * ComplexWorkspaceMS;

/** \brief Struct to provide real workspace for multistep methods
 *
 * Provide basic data to simplify the general multistep methods API,
 * with max number of previous steps required, system size and array
 * with all previous derivatives pre-computed
 */
typedef struct{
    int
        ms_order,       /// number of previous steps required
        system_size;    /// number of equations in ODE system
    Rarray
        prev_der;       /// Hold all required previous derivatives
} _RealWorkspaceMS;

/** \brief Struct address with working array for multistep methods */
typedef _RealWorkspaceMS * RealWorkspaceMS;


/** \brief Alloc struct internal array based on its integer fields */
void
alloc_cmplx_multistep_array(ComplexWorkspaceMS);

/** \brief Alloc struct internal array based on its integer fields */
void
alloc_real_multistep_array(RealWorkspaceMS);

/** \brief Free internal struct pointer */
void
free_cmplx_multistep_array(ComplexWorkspaceMS);

/** \brief Free internal struct pointer */
void
free_real_multistep_array(RealWorkspaceMS);

/** \brief Return fresh allocated struct address with internal fields set
 *
 * \param 1 : system size
 * \param 2 : multistep order (number of previous steps required)
 */
ComplexWorkspaceMS
get_cmplx_multistep_ws(unsigned int, unsigned int);

/** \brief Return fresh allocated struct address with internal fields set
 *
 * \param 1 : system size
 * \param 2 : multistep order (number of previous steps required)
 */
RealWorkspaceMS
get_real_multistep_ws(unsigned int, unsigned int);


/** \brief Free allocated complex workspace struct and its internal pointer */
void
free_cmplx_multistep_ws(ComplexWorkspaceMS);


/** \brief Free allcated real workspace struct and its internal pointer */
void
free_real_multistep_ws(RealWorkspaceMS);

/** \brief Prepare derivatives and solution to propagate next step */
void
cmplx_set_next_step(
        double,
        cmplx_sys_der,
        void *,
        ComplexWorkspaceMS,
        Carray,
        Carray
);

/** \brief Prepare derivatives and solution to propagate next step */
void
real_set_next_step(
        double,
        real_sys_der,
        void *,
        RealWorkspaceMS,
        Rarray,
        Rarray
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


/** \brief 4th order Adams-Bashforth(P)-Moulton(C) Predictor-Corrector step
 *
 * This routine carry out one step evolution of an ODE system using the 4th
 * order Adams-Bashforth for predictor and Adams-Moulton for corrector with
 * particular values for `real_general_multistep` routine
 *
 * \param 1 : grid spacing `h`
 * \param 2 : grid point correspongind to function values `x`
 * \param 3 : function pointer to routine that compute derivatives
 * \param 4 : extra arguments required in `cmplx_sys_der` function
 * \param 5 : Workspace struct address to avoid memory allocation.
 *            The array within this struct must have at least size
 *            `(m + 1) * n` due to implicit computation, must have
 *            derivative of previous steps concatenated as
 *            `[y'_j y'_j-1 ...  y'_j-3]`
 * \param 6 : Concatenated function steps required: `[y_j y_j-1 ...  y_j-3]`
 * \param 7 : Number of iterations for implicit part (Moulton), if
 *            zero, perform only explicit part (Bashforth)
 * \param 8: (OUTPUT) solution at next grid step
 */
void
real_adams4pc(
        double,
        double,
        real_sys_der,
        void *,
        RealWorkspaceMS,
        Rarray,
        unsigned int,
        Rarray
);


/** \brief 4th order Adams-Bashforth(P)-Moulton(C) Predictor-Corrector step
 *
 * This routine carry out one step evolution of an ODE system using the 4th
 * order Adams-Bashforth for predictor and Adams-Moulton for corrector with
 * particular values for `real_general_multistep` routine
 *
 * \param 1 : grid spacing `h`
 * \param 2 : grid point correspongind to function values `x`
 * \param 3 : function pointer to routine that compute derivatives
 * \param 4 : extra arguments required in `cmplx_sys_der` function
 * \param 5 : Workspace struct address to avoid memory allocation.
 *            The array within this struct must have at least size
 *            `(m + 1) * n` due to implicit computation, must have
 *            derivative of previous steps concatenated as
 *            `[y'_j y'_j-1 ...  y'_j-3]`
 * \param 6 : Concatenated function steps required: `[y_j y_j-1 ...  y_j-3]`
 * \param 7 : Number of iterations for implicit part (Moulton), if
 *            zero, perform only explicit part (Bashforth)
 * \param 8: (OUTPUT) solution at next grid step
 */
void
cmplx_adams4pc(
        double,
        double,
        cmplx_sys_der,
        void *,
        ComplexWorkspaceMS,
        Carray,
        unsigned int,
        Carray
);


#endif
