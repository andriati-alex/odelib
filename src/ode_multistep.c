#include "ode_multistep.h"


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
            summ = b[0] * der[i + m * s];
            for (j = 1; j <= m; j++)
            {
                stride = i + (j - 1) * s;
                summ = summ + h * b[j] * der[stride] - a[j] * y[stride];
            }
            ynext[i] = summ;
        }
        iter--;
    }
    cmplx_set_next_step(x, yprime, args, ws, y, ynext);
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
            summ = b[0] * der[i + m * s];
            for (j = 1; j <= m; j++)
            {
                stride = i + (j - 1) * s;
                summ = summ + h * b[j] * der[stride] - a[j] * y[stride];
            }
            ynext[i] = summ;
        }
        iter--;
    }
    real_set_next_step(x, yprime, args, ws, y, ynext);
}
