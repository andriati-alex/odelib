/**
 * \file ode_singlestep.h
 * \author Alex Andriati
 * \brief ODE integration routines with explicit single step methods
 *
 * Numerical integration methods of ODE systems that require only two
 * adjacent steps are named single step methods. Moreover, a subclass
 * of these methods are the explicit single step methods, a.k.a Runge-
 * Kutta methods, which are treated in this file
 */

#ifndef ODE_EXPLICIT_SINGLESTEP_H
#define ODE_EXPLICIT_SINGLESTEP_H

#include "derivative_signature.h"

/** \brief Struct to provide workspace for single step methods */
typedef struct{
    int system_size;
    Carray
        work1,
        work2,
        work3,
        work4,
        work5;
} _ComplexWorkspaceRK;

/** \brief Struct address working arrays for complex integration routines */
typedef _ComplexWorkspaceRK * ComplexWorkspaceRK;

/** \brief Struct to provide workspace for single step methods */
typedef struct{
    int system_size;
    Rarray
        work1,
        work2,
        work3,
        work4,
        work5;
} _RealWorkspaceRK;

/** \brief Struct address working arrays for real integration routines */
typedef _RealWorkspaceRK * RealWorkspaceRK;

/**
 * \brief 4th order Runge-Kutta method step integration
 *
 * \param 1 : grid spacing `h`
 * \param 2 : grid point correspongind to function values `x`
 * \param 3 : function pointer to routine that compute derivatives
 * \param 4 : extra arguments required in `cmplx_sys_der` function
 * \param 5 : Workspace struct address to avoid memory allocation
 * \param 6 : function values `y` computed at grid point `x`
 * \param 7 : (OUTPUT) function values at next grid point `x + h`
 */
void
cmplx_rungekutta4
(
        double,
        double,
        cmplx_sys_der,
        void *,
        ComplexWorkspaceRK,
        Carray,
        Carray
);


/**
 * \brief 4th order Runge-Kutta method step integration
 *
 * \param 1 : grid spacing `h`
 * \param 2 : grid point correspongind to function values `x`
 * \param 3 : function pointer to routine that compute derivatives
 * \param 4 : extra arguments required in `cmplx_sys_der` function
 * \param 5 : Workspace struct address to avoid memory allocation
 * \param 6 : function values `y` computed at grid point `x`
 * \param 7 : (OUTPUT) function values at next grid point `x + h`
 */
void
real_rungekutta4
(
        double dx,
        double x,
        real_sys_der yprime,
        void * args,
        RealWorkspaceRK ws,
        Rarray y,
        Rarray ynext
);


/**
 * \brief 2nd order Runge-Kutta method step integration
 *
 * \param 1 : grid spacing `h`
 * \param 2 : grid point correspongind to function values `x`
 * \param 3 : function pointer to routine that compute derivatives
 * \param 4 : extra arguments required in `cmplx_sys_der` function
 * \param 5 : Workspace struct address to avoid memory allocation
 * \param 6 : function values `y` computed at grid point `x`
 * \param 7 : (OUTPUT) function values at next grid point `x + h`
 */
void
cmplx_rungekutta2
(
        double,
        double,
        cmplx_sys_der,
        void *,
        ComplexWorkspaceRK,
        Carray,
        Carray
);


/**
 * \brief 2nd order Runge-Kutta method step integration
 *
 * \param 1 : grid spacing `h`
 * \param 2 : grid point correspongind to function values `x`
 * \param 3 : function pointer to routine that compute derivatives
 * \param 4 : extra arguments required in `cmplx_sys_der` function
 * \param 5 : Workspace struct address to avoid memory allocation
 * \param 6 : function values `y` computed at grid point `x`
 * \param 7 : (OUTPUT) function values at next grid point `x + h`
 */
void
real_rungekutta2
(
        double dx,
        double x,
        real_sys_der yprime,
        void * args,
        RealWorkspaceRK ws,
        Rarray y,
        Rarray ynext
);


#endif
