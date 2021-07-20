#include "ode_multistep.h"
#include <stdio.h>


void
alloc_cmplx_multistep_array(ComplexWorkspaceMS ws)
{
    unsigned int
        full_size = (ws->ms_order + 1) * ws->system_size;
    ws->prev_der = (Carray) malloc(full_size * sizeof(double complex));
}


void
alloc_real_multistep_array(RealWorkspaceMS ws)
{
    unsigned int
        full_size = (ws->ms_order + 1) * ws->system_size;
    ws->prev_der = (Rarray) malloc(full_size * sizeof(double));
}


void
free_cmplx_multistep_array(ComplexWorkspaceMS ws)
{
    free(ws->prev_der);
}


void
free_real_multistep_array(RealWorkspaceMS ws)
{
    free(ws->prev_der);
}


ComplexWorkspaceMS
get_cmplx_multistep_ws(unsigned int ms_order, unsigned int sys_size)
{
    ComplexWorkspaceMS
        ws;
    ws = (ComplexWorkspaceMS) malloc(sizeof(_ComplexWorkspaceMS));
    ws->ms_order = ms_order;
    ws->system_size = sys_size;
    alloc_cmplx_multistep_array(ws);
    return ws;
}


RealWorkspaceMS
get_real_multistep_ws(unsigned int ms_order, unsigned int sys_size)
{
    RealWorkspaceMS
        ws;
    ws = (RealWorkspaceMS) malloc(sizeof(_RealWorkspaceMS));
    ws->ms_order = ms_order;
    ws->system_size = sys_size;
    alloc_real_multistep_array(ws);
    return ws;
}


void
free_cmplx_multistep_ws(ComplexWorkspaceMS ws)
{
    free(ws->prev_der);
    free(ws);
}


void
free_real_multistep_ws(RealWorkspaceMS ws)
{
    free(ws->prev_der);
    free(ws);
}


void
cmplx_set_next_step(
        double x,
        cmplx_sys_der yprime,
        void * args,
        ComplexWorkspaceMS ws,
        Carray y,
        Carray ynext
)
{
    int
        i,
        j,
        m,
        s;
    Carray
        der;

    m = ws->ms_order;
    s = ws->system_size;
    der = ws->prev_der;

    for (j = m - 1; j > 0; j--)
    {
        for (i = 0; i < s; i++)
        {
            y[i + j * s] = y[i + (j - 1) * s];
            der[i + j * s] = der[i + (j - 1) * s];
        }
    }
    for (i = 0; i < s; i++)
    {
        y[i] = ynext[i];
    }
    yprime(s, x, ynext, der, args);
}


void
real_set_next_step(
        double x,
        real_sys_der yprime,
        void * args,
        RealWorkspaceMS ws,
        Rarray y,
        Rarray ynext
)
{
    int
        i,
        j,
        m,
        s;
    Rarray
        der;

    m = ws->ms_order;
    s = ws->system_size;
    der = ws->prev_der;

    for (j = m - 1; j > 0; j--)
    {
        for (i = 0; i < s; i++)
        {
            y[i + j * s] = y[i + (j - 1) * s];
            der[i + j * s] = der[i + (j - 1) * s];
        }
    }
    for (i = 0; i < s; i++)
    {
        y[i] = ynext[i];
    }
    yprime(s, x, ynext, der, args);
}


void
cmplx_general_multistep(
        double h,
        double x,
        cmplx_sys_der yprime,
        void * args,
        ComplexWorkspaceMS ws,
        Carray y,
        Rarray a,
        Rarray b,
        unsigned int iter,
        Carray ynext
)
{
    int
        i,
        j,
        m,
        s,
        stride;
    double complex
        summ;
    Carray
        der;

    m = ws->ms_order;
    s = ws->system_size;
    der = ws->prev_der;

    if (!iter)
    {
        for (i = 0; i < s; i++)
        {
            summ = 0;
            for (j = 1; j <= m; j++)
            {
                stride = i + (j - 1) * s;
                summ = summ + h * b[j] * der[stride] - a[j] * y[stride];
            }
            ynext[i] = summ;
        }
    }
    /* Implicit scheme used as corrector *
     * `ynext` must provide a prediction */
    while (iter > 0)
    {
        yprime(s, x, ynext, &der[m * s], args);
        for (i = 0; i < s; i++)
        {
            summ = h * b[0] * der[i + m * s];
            for (j = 1; j <= m; j++)
            {
                stride = i + (j - 1) * s;
                summ = summ + h * b[j] * der[stride] - a[j] * y[stride];
            }
            ynext[i] = summ;
        }
        iter--;
    }
}


void
real_general_multistep(
        double h,
        double x,
        real_sys_der yprime,
        void * args,
        RealWorkspaceMS ws,
        Rarray y,
        Rarray a,
        Rarray b,
        unsigned int iter,
        Rarray ynext
)
{
    int
        i,
        j,
        m,
        s,
        stride;
    double complex
        summ;
    Rarray
        der;

    m = ws->ms_order;
    s = ws->system_size;
    der = ws->prev_der;

    if (!iter)
    {
        for (i = 0; i < s; i++)
        {
            summ = 0;
            for (j = 1; j <= m; j++)
            {
                stride = i + (j - 1) * s;
                summ = summ + h * b[j] * der[stride] - a[j] * y[stride];
            }
            ynext[i] = summ;
        }
    }
    /* Implicit scheme used as corrector *
     * `ynext` must provide a prediction */
    while (iter > 0)
    {
        yprime(s, x, ynext, &der[m * s], args);
        for (i = 0; i < s; i++)
        {
            summ = h * b[0] * der[i + m * s];
            for (j = 1; j <= m; j++)
            {
                stride = i + (j - 1) * s;
                summ = summ + h * b[j] * der[stride] - a[j] * y[stride];
            }
            ynext[i] = summ;
        }
        iter--;
    }
}
