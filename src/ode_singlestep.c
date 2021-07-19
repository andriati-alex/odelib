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
    ws->work1 = (Carray) malloc(ws->system_size * sizeof(double complex));
    ws->work2 = (Carray) malloc(ws->system_size * sizeof(double complex));
    ws->work3 = (Carray) malloc(ws->system_size * sizeof(double complex));
    ws->work4 = (Carray) malloc(ws->system_size * sizeof(double complex));
    ws->work5 = (Carray) malloc(ws->system_size * sizeof(double complex));
}


void
alloc_real_rkwsarrays(RealWorkspaceRK ws)
{
    ws->work1 = (Rarray) malloc(ws->system_size * sizeof(double));
    ws->work2 = (Rarray) malloc(ws->system_size * sizeof(double));
    ws->work3 = (Rarray) malloc(ws->system_size * sizeof(double));
    ws->work4 = (Rarray) malloc(ws->system_size * sizeof(double));
    ws->work5 = (Rarray) malloc(ws->system_size * sizeof(double));
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
    ws->system_size = sys_size;
    alloc_cmplx_rkwsarrays(ws);
    return ws;
}


RealWorkspaceRK
get_real_rkws(int sys_size)
{
    RealWorkspaceRK
        ws = (RealWorkspaceRK) malloc(sizeof(_RealWorkspaceRK));
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


void free_cmplx_rkws(ComplexWorkspaceRK ws)
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

    sys_size = ws->system_size;
    k1 = ws->work1;
    k2 = ws->work2;
    k3 = ws->work3;
    k4 = ws->work4;
    karg = ws->work5;

    /* Start 4-th order Runge-Kutta algorithm as in [1] Equation (2.11.5) */
    yprime(sys_size, x, y, k1, args);
    for (i = 0; i < sys_size; i++)
    {
        karg[i] = 0.5 * h * k1[i] + y[i];
    }
    yprime(sys_size, x + 0.5 * h, karg, k2, args);
    for (i = 0; i < sys_size; i++)
    {
        karg[i] = 0.5 * h * k2[i] + y[i];
    }
    yprime(sys_size, x + 0.5 * h, karg, k3, args);
    for (i = 0; i < sys_size; i++)
    {
        karg[i] = h * k3[i] + y[i];
    }
    yprime(sys_size, x + h, karg, k4, args);
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

    sys_size = ws->system_size;
    k1 = ws->work1;
    k2 = ws->work2;
    k3 = ws->work3;
    k4 = ws->work4;
    karg = ws->work5;

    /* Start 4-th order Runge-Kutta algorithm as in [1] Equation (2.11.5) */
    yprime(sys_size, x, y, k1, args);
    for (i = 0; i < sys_size; i++)
    {
        karg[i] = 0.5 * h * k1[i] + y[i];
    }
    yprime(sys_size, x + 0.5 * h, karg, k2, args);
    for (i = 0; i < sys_size; i++)
    {
        karg[i] = 0.5 * h * k2[i] + y[i];
    }
    yprime(sys_size, x + 0.5 * h, karg, k3, args);
    for (i = 0; i < sys_size; i++)
    {
        karg[i] = h * k3[i] + y[i];
    }
    yprime(sys_size, x + h, karg, k4, args);
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

    sys_size = ws->system_size;
    k1 = ws->work1;
    k2 = ws->work2;
    karg = ws->work3;

    /* start 2nd order Runge-Kutta scheme as in [1] Equation (2.5.2) */
    yprime(sys_size, x, y, k1, args);
    for(i = 0; i < sys_size; i++)
    {
        karg[i] = h * k1[i] + y[i];
    }
    yprime(sys_size, x + h, karg, k2, args);
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

    sys_size = ws->system_size;
    k1 = ws->work1;
    k2 = ws->work2;
    karg = ws->work3;

    /* start 2nd order Runge-Kutta scheme as in [1] Equation (2.5.2) */
    yprime(sys_size, x, y, k1, args);
    for(i = 0; i < sys_size; i++)
    {
        karg[i] = h * k1[i] + y[i];
    }
    yprime(sys_size, x + h, karg, k2, args);
    for (i = 0; i < sys_size; i++)
    {
        ynext[i] = y[i] + 0.5 * h * (k1[i] + k2[i]);
    }
}
