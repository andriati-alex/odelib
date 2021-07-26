/**
 * \file ode_singlestep.h
 * \author Alex Andriati
 * \brief ODE integration routines with explicit single step methods
 *
 * Numerical integration methods of ODE systems that require only one
 * previous adjacent step are named single step methods. Moreover, a
 * subclass of these methods are the explicit single step methods, a.k.a
 * Runge-Kutta methods, which are treated in this file
 */

#ifndef ODE_SINGLESTEP_H
#define ODE_SINGLESTEP_H

#include "derivative_signature.h"

/** \brief Struct to provide workspace for single step methods
 *
 * The arrays inside this structure hold intermediate steps in
 * Runge-Kutta methods which require derivative evaluations
 */
typedef struct{
    int
        system_size;
    Carray
        work1,
        work2,
        work3,
        work4,
        work5;
} _ComplexWorkspaceRK;

/** \brief Struct workspace address for single step methods */
typedef _ComplexWorkspaceRK * ComplexWorkspaceRK;

/** \brief Struct to provide workspace for single step methods
 *
 * The arrays inside this structure hold intermediate steps in
 * Runge-Kutta methods which require derivative evaluations
 */
typedef struct{
    int system_size;
    Rarray
        work1,
        work2,
        work3,
        work4,
        work5;
} _RealWorkspaceRK;

/** \brief Struct workspace address for single step methods */
typedef _RealWorkspaceRK * RealWorkspaceRK;


/** \brief Alloc internal arrays in struct address given */
void
alloc_cplx_rungekutta_wsarrays(ComplexWorkspaceRK);


/** \brief Alloc internal arrays in struct address given */
void
alloc_real_rungekutta_wsarrays(RealWorkspaceRK);


/** \brief Free internal pointers of struct given */
void
free_cplx_rungekutta_wsarrays(ComplexWorkspaceRK);


/** \brief Free internal pointers of struct given */
void
free_real_rungekutta_wsarrays(RealWorkspaceRK);


/** \brief Return fresh allocated struct address with internal fields set */
ComplexWorkspaceRK
get_cplx_rungekutta_ws(int sys_size);


/** \brief Return fresh allocated struct address with internal fields set */
RealWorkspaceRK
get_real_rungekutta_ws(int sys_size);


/** \brief Free allocated workspace struct and its internal arrays */
void
destroy_real_rungekutta_ws(RealWorkspaceRK);


/** \brief Free allocated workspace struct and its internal arrays */
void
destroy_cplx_rungekutta_ws(ComplexWorkspaceRK);


/**
 * \brief 4th order Runge-Kutta method step integration
 *
 * \param 1 : grid spacing `h`
 * \param 2 : current grid point `x`
 * \param 3 : function pointing to routine that compute derivatives
 * \param 4 : extra arguments (void pointer in _ComplexWorkspaceRK)
 * \param 5 : Workspace struct address for internal derivative computation
 * \param 6 : function values `y` computed at current grid point
 * \param 7 : (OUTPUT) function values at next grid point `x + h`
 */
void
cplx_rungekutta4
(
        double,
        double,
        cplx_odesys_der,
        void *,
        ComplexWorkspaceRK,
        Carray,
        Carray
);


/**
 * \brief 4th order Runge-Kutta method step integration
 *
 * \param 1 : grid spacing `h`
 * \param 2 : current grid point `x`
 * \param 3 : function pointing to routine that compute derivatives
 * \param 4 : extra arguments (void pointer in _RealWorkspaceRK)
 * \param 5 : Workspace struct address for internal derivative computation
 * \param 6 : function values `y` computed at current grid point `x`
 * \param 7 : (OUTPUT) function values at next grid point `x + h`
 */
void
real_rungekutta4
(
        double,
        double,
        real_odesys_der,
        void *,
        RealWorkspaceRK,
        Rarray,
        Rarray
);


/**
 * \brief 2nd order (simple)Runge-Kutta method step integration
 *
 * \param 1 : grid spacing `h`
 * \param 2 : current grid point `x`
 * \param 3 : function pointing to routine that compute derivatives
 * \param 4 : extra arguments (void pointer in _ComplexWorkspaceRK)
 * \param 5 : Workspace struct address to avoid memory allocation
 * \param 6 : function values `y` computed at current grid point `x`
 * \param 7 : (OUTPUT) function values at next grid point `x + h`
 */
void
cplx_rungekutta2
(
        double,
        double,
        cplx_odesys_der,
        void *,
        ComplexWorkspaceRK,
        Carray,
        Carray
);


/**
 * \brief 2nd order Runge-Kutta method step integration
 *
 * \param 1 : grid spacing `h`
 * \param 2 : current grid point `x`
 * \param 3 : function pointing to routine that compute derivatives
 * \param 4 : extra arguments (void pointer in _RealWorkspaceRK)
 * \param 5 : Workspace struct address to avoid memory allocation
 * \param 6 : function values `y` computed at current grid point `x`
 * \param 7 : (OUTPUT) function values at next grid point `x + h`
 */
void
real_rungekutta2
(
        double,
        double,
        real_odesys_der,
        void *,
        RealWorkspaceRK,
        Rarray,
        Rarray
);


#endif
