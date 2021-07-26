#include "ode_multistep.h"
#include "arrays_assistant.h"


void
alloc_cplx_multistep_wsarray(ComplexWorkspaceMS ws)
{
    unsigned int
        full_size = (ws->ms_order + 1) * ws->system_size;
    ws->prev_der = alloc_carr(full_size);
}


void
alloc_real_multistep_wsarray(RealWorkspaceMS ws)
{
    unsigned int
        full_size = (ws->ms_order + 1) * ws->system_size;
    ws->prev_der = alloc_rarr(full_size);
}


void
free_cplx_multistep_wsarray(ComplexWorkspaceMS ws)
{
    free(ws->prev_der);
}


void
free_real_multistep_wsarray(RealWorkspaceMS ws)
{
    free(ws->prev_der);
}


ComplexWorkspaceMS
get_cplx_multistep_ws(unsigned int ms_order, unsigned int sys_size)
{
    ComplexWorkspaceMS
        ws;
    ws = (ComplexWorkspaceMS) malloc(sizeof(_ComplexWorkspaceMS));
    if (ws == NULL)
    {
        printf("\n\nProblem in ComplexWorkspaceMS allocation\n\n");
        exit(EXIT_FAILURE);
    }
    ws->ms_order = ms_order;
    ws->system_size = sys_size;
    alloc_cplx_multistep_wsarray(ws);
    return ws;
}


RealWorkspaceMS
get_real_multistep_ws(unsigned int ms_order, unsigned int sys_size)
{
    RealWorkspaceMS
        ws;
    ws = (RealWorkspaceMS) malloc(sizeof(_RealWorkspaceMS));
    if (ws == NULL)
    {
        printf("\n\nProbelm in RealWorkspaceMS allocation\n\n");
        exit(EXIT_FAILURE);
    }
    ws->ms_order = ms_order;
    ws->system_size = sys_size;
    alloc_real_multistep_wsarray(ws);
    return ws;
}


void
destroy_cplx_multistep_ws(ComplexWorkspaceMS ws)
{
    free(ws->prev_der);
    free(ws);
}


void
destroy_real_multistep_ws(RealWorkspaceMS ws)
{
    free(ws->prev_der);
    free(ws);
}


void
cplx_set_next_multistep(
        double xnext,
        cplx_odesys_der yprime,
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
    _ComplexODEInputParameters
        sys_params;

    m = ws->ms_order;
    s = ws->system_size;
    der = ws->prev_der;

    sys_params.x = xnext;
    sys_params.y = ynext;
    sys_params.system_size = s;
    sys_params.extra_args = args;

    for (j = m - 1; j > 0; j--)
    {
        for (i = 0; i < s; i++)
        {
            y[i + j * s] = y[i + (j - 1) * s];
            der[i + j * s] = der[i + (j - 1) * s];
        }
    }
    carr_copy_values(s, ynext, y);
    yprime(&sys_params, der);
}


void
real_set_next_multistep(
        double xnext,
        real_odesys_der yprime,
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
    _RealODEInputParameters
        sys_params;

    m = ws->ms_order;
    s = ws->system_size;
    der = ws->prev_der;

    sys_params.x = xnext;
    sys_params.y = ynext;
    sys_params.system_size = s;
    sys_params.extra_args = args;

    for (j = m - 1; j > 0; j--)
    {
        for (i = 0; i < s; i++)
        {
            y[i + j * s] = y[i + (j - 1) * s];
            der[i + j * s] = der[i + (j - 1) * s];
        }
    }
    rarr_copy_values(s, ynext, y);
    yprime(&sys_params, der);
}


void
cplx_general_multistep(
        double h,
        double x,
        cplx_odesys_der yprime,
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
        return;
    }

    /* Implicit scheme used as corrector *
     * `ynext` must provide a prediction */
    _ComplexODEInputParameters
        sys_params;
    sys_params.x = x + h;
    sys_params.y = ynext;
    sys_params.extra_args = args;
    sys_params.system_size = ws->system_size;
    while (iter > 0)
    {
        yprime(&sys_params, &der[m * s]);
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
        real_odesys_der yprime,
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
    double
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
        return;
    }

    /* Implicit scheme used as corrector *
     * `ynext` must provide a prediction */
    _RealODEInputParameters
        sys_params;
    sys_params.x = x + h;
    sys_params.y = ynext;
    sys_params.extra_args = args;
    sys_params.system_size = ws->system_size;
    while (iter > 0)
    {
        yprime(&sys_params, &der[m * s]);
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
real_adams4pc(
        double h,
        double x,
        real_odesys_der yprime,
        void * args,
        RealWorkspaceMS ws,
        Rarray y,
        unsigned int iter,
        Rarray ynext
)
{
    double
        ap[5] = {1.0, -1.0,   0.0,  0.0,  0.0},
        bp[5] = {0.0, 55.0, -59.0, 37.0, -9.0},
        ac[5] = {1.0, -1.0,   0.0,  0.0,  0.0},
        bc[5] = {9.0, 19.0,  -5.0,  1.0,  0.0};
    for (int i = 0; i < 5; i++)
    {
        bp[i] = bp[i] / 24;
        bc[i] = bc[i] / 24;
    }
    real_general_multistep(h, x, yprime, args, ws, y, ap, bp, 0, ynext);
    if (iter == 0) return;
    real_general_multistep(h, x, yprime, args, ws, y, ac, bc, iter, ynext);
}


void
cplx_adams4pc(
        double h,
        double x,
        cplx_odesys_der yprime,
        void * args,
        ComplexWorkspaceMS ws,
        Carray y,
        unsigned int iter,
        Carray ynext
)
{
    double
        ap[5] = {1.0, -1.0,   0.0,  0.0,  0.0},
        bp[5] = {0.0, 55.0, -59.0, 37.0, -9.0},
        ac[5] = {1.0, -1.0,   0.0,  0.0,  0.0},
        bc[5] = {9.0, 19.0,  -5.0,  1.0,  0.0};
    for (int i = 0; i < 5; i++)
    {
        bp[i] = bp[i] / 24;
        bc[i] = bc[i] / 24;
    }
    cplx_general_multistep(h, x, yprime, args, ws, y, ap, bp, 0, ynext);
    if (iter == 0) return;
    cplx_general_multistep(h, x, yprime, args, ws, y, ac, bc, iter, ynext);
}
