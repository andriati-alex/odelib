/**
 * \file ode_explicit_singlestep.c
 * \author Alex Andriati
 * \brief Source code for explicit single step ODE integration routines
 *
 * See function signature and description in header ode_singlestep.h
 * The algorithms are implemented providing the most possible similarity with
 * ref. [1], adapting to ODE system language whenever needed. Some references
 *
 * [1] Douglas Quinney, An introduction to the numerical solution of
 * differential equations, Revised Edition, 1987, cap. 2
 * [2] William H. Press et. al., Numerical Recipes in C, 2nd Edition
 * cap. 16
 * [3] Arieh Iserles, A first course in the numerical analysis of
 * differential equations, Cambridge, 2nd Edition, cap. 3
 */

#include "ode_singlestep.h"
#include <stdlib.h>


void
alloc_cmplx_rkwsarrays(ComplexWorkspaceRK ws)
{
    ws->work1 = alloc_carr(ws->system_size);
    ws->work2 = alloc_carr(ws->system_size);
    ws->work3 = alloc_carr(ws->system_size);
    ws->work4 = alloc_carr(ws->system_size);
    ws->work5 = alloc_carr(ws->system_size);
}


void
alloc_real_rkwsarrays(RealWorkspaceRK ws)
{
    ws->work1 = alloc_rarr(ws->system_size);
    ws->work2 = alloc_rarr(ws->system_size);
    ws->work3 = alloc_rarr(ws->system_size);
    ws->work4 = alloc_rarr(ws->system_size);
    ws->work5 = alloc_rarr(ws->system_size);
}


void
free_cmplx_rkwsarrays(ComplexWorkspaceRK ws)
{
    free(ws->work1);
    free(ws->work2);
    free(ws->work3);
    free(ws->work4);
    free(ws->work5);
}


void
free_real_rkwsarrays(RealWorkspaceRK ws)
{
    free(ws->work1);
    free(ws->work2);
    free(ws->work3);
    free(ws->work4);
    free(ws->work5);
}


ComplexWorkspaceRK
get_cmplx_rkws(int sys_size)
{
    ComplexWorkspaceRK
        ws = (ComplexWorkspaceRK) malloc(sizeof(_ComplexWorkspaceRK));
    if (ws == NULL)
    {
        printf("\n\nProblem in ComplexWorkspaceRK allocation\n\n");
        exit(EXIT_FAILURE);
    }
    ws->system_size = sys_size;
    alloc_cmplx_rkwsarrays(ws);
    return ws;
}


RealWorkspaceRK
get_real_rkws(int sys_size)
{
    RealWorkspaceRK
        ws = (RealWorkspaceRK) malloc(sizeof(_RealWorkspaceRK));
    if (ws == NULL)
    {
        printf("\n\nProblem in RealWorkspaceRK allocation\n\n");
        exit(EXIT_FAILURE);
    }
    ws->system_size = sys_size;
    alloc_real_rkwsarrays(ws);
    return ws;
}


void
free_real_rkws(RealWorkspaceRK ws)
{
    free_real_rkwsarrays(ws);
    free(ws);
}


void
free_cmplx_rkws(ComplexWorkspaceRK ws)
{
    free_cmplx_rkwsarrays(ws);
    free(ws);
}


void
cmplx_rungekutta4(
        double h,
        double x,
        cmplx_sys_der yprime,
        void * args,
        ComplexWorkspaceRK ws,
        Carray y,
        Carray ynext
)
{
    int
        i,
        sys_size;
    Carray
        k1,
        k2,
        k3,
        k4,
        karg;
    _ComplexODEInputParameters
        sys_params;

    sys_size = ws->system_size;
    k1 = ws->work1;
    k2 = ws->work2;
    k3 = ws->work3;
    k4 = ws->work4;
    karg = ws->work5;
    carr_copy_values(sys_size, y, karg);

    sys_params.y = karg;
    sys_params.extra_args = args;
    sys_params.system_size = sys_size;

    /* Start 4-th order Runge-Kutta algorithm as in [1] Equation (2.11.5) */
    sys_params.x = x;
    yprime(&sys_params, k1);
    for (i = 0; i < sys_size; i++)
    {
        karg[i] = 0.5 * h * k1[i] + y[i];
    }
    sys_params.x = x + 0.5 * h;
    yprime(&sys_params, k2);
    for (i = 0; i < sys_size; i++)
    {
        karg[i] = 0.5 * h * k2[i] + y[i];
    }
    sys_params.x = x + 0.5 * h;
    yprime(&sys_params, k3);
    for (i = 0; i < sys_size; i++)
    {
        karg[i] = h * k3[i] + y[i];
    }
    sys_params.x = x + h;
    yprime(&sys_params, k4);
    for (i = 0; i < sys_size; i++)
    {
        ynext[i] = y[i] + (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) * h / 6;
    }
}


void
real_rungekutta4(
        double h,
        double x,
        real_sys_der yprime,
        void * args,
        RealWorkspaceRK ws,
        Rarray y,
        Rarray ynext
)
{
    int
        i,
        sys_size;
    Rarray
        k1,
        k2,
        k3,
        k4,
        karg;
    _RealODEInputParameters
        sys_params;

    sys_size = ws->system_size;
    k1 = ws->work1;
    k2 = ws->work2;
    k3 = ws->work3;
    k4 = ws->work4;
    karg = ws->work5;
    rarr_copy_values(sys_size, y, karg);

    sys_params.y = karg;
    sys_params.extra_args = args;
    sys_params.system_size = sys_size;

    /* Start 4-th order Runge-Kutta algorithm as in [1] Equation (2.11.5) */
    sys_params.x = x;
    yprime(&sys_params, k1);
    for (i = 0; i < sys_size; i++)
    {
        karg[i] = 0.5 * h * k1[i] + y[i];
    }
    sys_params.x = x + 0.5 * h;
    yprime(&sys_params, k2);
    for (i = 0; i < sys_size; i++)
    {
        karg[i] = 0.5 * h * k2[i] + y[i];
    }
    sys_params.x = x + 0.5 * h;
    yprime(&sys_params, k3);
    for (i = 0; i < sys_size; i++)
    {
        karg[i] = h * k3[i] + y[i];
    }
    sys_params.x = x + h;
    yprime(&sys_params, k4);
    for (i = 0; i < sys_size; i++)
    {
        ynext[i] = y[i] + (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) * h / 6;
    }
}


void
cmplx_rungekutta2(
        double h,
        double x,
        cmplx_sys_der yprime,
        void * args,
        ComplexWorkspaceRK ws,
        Carray y,
        Carray ynext
)
{
    int
        i,
        sys_size;
    Carray
        k1,
        k2,
        karg;
    _ComplexODEInputParameters
        sys_params;

    sys_size = ws->system_size;
    k1 = ws->work1;
    k2 = ws->work2;
    karg = ws->work3;
    carr_copy_values(sys_size, y, karg);

    sys_params.y = karg;
    sys_params.extra_args = args;
    sys_params.system_size = sys_size;

    /* start 2nd order Runge-Kutta scheme as in [1] Equation (2.5.2) */
    sys_params.x = x;
    yprime(&sys_params, k1);
    for(i = 0; i < sys_size; i++)
    {
        karg[i] = h * k1[i] + y[i];
    }
    sys_params.x = x + h;
    yprime(&sys_params, k2);
    for (i = 0; i < sys_size; i++)
    {
        ynext[i] = y[i] + 0.5 * h * (k1[i] + k2[i]);
    }
}


void
real_rungekutta2(
        double h,
        double x,
        real_sys_der yprime,
        void * args,
        RealWorkspaceRK ws,
        Rarray y,
        Rarray ynext
)
{
    int
        i,
        sys_size;
    Rarray
        k1,
        k2,
        karg;
    _RealODEInputParameters
        sys_params;

    sys_size = ws->system_size;
    k1 = ws->work1;
    k2 = ws->work2;
    karg = ws->work3;
    rarr_copy_values(sys_size, y, karg);

    sys_params.y = karg;
    sys_params.extra_args = args;
    sys_params.system_size = sys_size;

    /* start 2nd order Runge-Kutta scheme as in [1] Equation (2.5.2) */
    sys_params.x = x;
    yprime(&sys_params, k1);
    for(i = 0; i < sys_size; i++)
    {
        karg[i] = h * k1[i] + y[i];
    }
    sys_params.x = x + h;
    yprime(&sys_params, k2);
    for (i = 0; i < sys_size; i++)
    {
        ynext[i] = y[i] + 0.5 * h * (k1[i] + k2[i]);
    }
}
