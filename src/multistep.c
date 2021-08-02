/**
 * \file multistep.c
 * \author Alex Andriati
 * \brief Source code for multistep methods
 *
 * See function signature and description in header multistep.h
 * Multistep methods have a common general form, using some coefficients
 * to combine previous steps and their derivatives. Therefore, there are
 * general routines for real and complex systems which consume as input
 * these coefficients, as well directed routines for specific approaches
 * which have specific coefficients. Among these routines with ready to
 * use coefficients are the Adams predictor-corrector of order 4 and 6.
 *
 * [1] Douglas Quinney, An introduction to the numerical solution of
 * differential equations, Revised Edition, 1987, cap. 2
 * [2] J.C. Butcher, Numerical methods for ordinary differential equations,
 * Wiley, 3rd Edition
 * [3] William H. Press et. al., Numerical Recipes in C, 2nd Edition
 * cap. 16
 * [4] Arieh Iserles, A first course in the numerical analysis of
 * differential equations, Cambridge, 2nd Edition, cap. 3
 */

#include "multistep.h"
#include "arrays_assistant.h"


static double
    ADAMS4_LEFT[5] = {1.0, -1.0, 0.0, 0.0, 0.0},
    ADAMS4_PRED[5] = {0.0, 55.0 / 24, -59.0 / 24, 37.0 / 24, -9.0 / 24},
    ADAMS4_CORR[5] = {9.0 / 24, 19.0 / 24, -5.0 / 24, 1.0 / 24, 0.0};

static double
    ADAMS6_LEFT[7] = {1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    ADAMS6_PRED[7] = {0.0, 4277.0 / 1440, -7923.0 / 1440, 9982.0 / 1440, -7298.0 / 1440, 2877.0 / 1440, -475.0 / 1440},
    ADAMS6_CORR[7] = {475.0 / 1440, 1427.0 / 1440, -798.0 / 1440, 482.0 / 1440, -173.0 / 1440, 27.0 / 1440, 0.0};


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


void
init_real_multistep(
        double h,
        real_odesys_der yprime,
        void * args,
        RealWorkspaceMS ws,
        Rarray y0,
        real_rk_routine rk,
        Rarray yms_init
)
{
    int
        i,
        j,
        sys_size;
    Rarray
        ywork;
    RealWorkspaceRK
        wsrk;
    _RealODEInputParameters
        inp;

    sys_size = ws->system_size;
    ywork = alloc_rarr(sys_size);
    wsrk = get_real_rungekutta_ws(sys_size);
    rarr_copy_values(sys_size, y0, ywork);

    inp.x = 0;
    inp.y = ywork;
    inp.extra_args = args;
    inp.system_size = sys_size;

    j = (ws->ms_order - 1) * sys_size;
    rarr_copy_values(sys_size, ywork, &yms_init[j]);
    yprime(&inp, &ws->prev_der[j]);

    for (i = 1; i < ws->ms_order; i++)
    {
        j = (ws->ms_order - 1 - i) * sys_size;
        (*rk)(h, inp.x, yprime, args, wsrk, ywork, &yms_init[j]);
        rarr_copy_values(sys_size, &yms_init[j], ywork);
        inp.x = i * h;
        yprime(&inp, &ws->prev_der[j]);
    }

    free(ywork);
    destroy_real_rungekutta_ws(wsrk);
}


void
init_cplx_multistep(
        double h,
        cplx_odesys_der yprime,
        void * args,
        ComplexWorkspaceMS ws,
        Carray y0,
        cplx_rk_routine rk,
        Carray yms_init
)
{
    int
        i,
        j,
        sys_size;
    Carray
        ywork;
    ComplexWorkspaceRK
        wsrk;
    _ComplexODEInputParameters
        inp;

    sys_size = ws->system_size;
    ywork = alloc_carr(sys_size);
    wsrk = get_cplx_rungekutta_ws(sys_size);
    carr_copy_values(sys_size, y0, ywork);

    inp.x = 0;
    inp.y = ywork;
    inp.extra_args = args;
    inp.system_size = sys_size;

    j = (ws->ms_order - 1) * sys_size;
    carr_copy_values(sys_size, ywork, &yms_init[j]);
    yprime(&inp, &ws->prev_der[j]);

    for (i = 1; i < ws->ms_order; i++)
    {
        j = (ws->ms_order - 1 - i) * sys_size;
        (*rk)(h, inp.x, yprime, args, wsrk, ywork, &yms_init[j]);
        carr_copy_values(sys_size, &yms_init[j], ywork);
        inp.x = i * h;
        yprime(&inp, &ws->prev_der[j]);
    }

    free(ywork);
    destroy_cplx_rungekutta_ws(wsrk);
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

    /* shift chunks representing concatenated previous steps */
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

    /* shift chunks representing concatenated previous steps */
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
    real_general_multistep(
            h, x, yprime, args, ws, y, ADAMS4_LEFT, ADAMS4_PRED, 0, ynext
    );
    if (iter == 0) return;
    real_general_multistep(
            h, x, yprime, args, ws, y, ADAMS4_LEFT, ADAMS4_CORR, iter, ynext
    );
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
    cplx_general_multistep(
            h, x, yprime, args, ws, y, ADAMS4_LEFT, ADAMS4_PRED, 0, ynext
    );
    if (iter == 0) return;
    cplx_general_multistep(
            h, x, yprime, args, ws, y, ADAMS4_LEFT, ADAMS4_CORR, iter, ynext
    );
}


void
cplx_adams6pc(
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
    cplx_general_multistep(
            h, x, yprime, args, ws, y, ADAMS6_LEFT, ADAMS6_PRED, 0, ynext
    );
    if (iter == 0) return;
    cplx_general_multistep(
            h, x, yprime, args, ws, y, ADAMS6_LEFT, ADAMS6_CORR, iter, ynext
    );
}


void
real_adams6pc(
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
    real_general_multistep(
            h, x, yprime, args, ws, y, ADAMS6_LEFT, ADAMS6_PRED, 0, ynext
    );
    if (iter == 0) return;
    real_general_multistep(
            h, x, yprime, args, ws, y, ADAMS6_LEFT, ADAMS6_CORR, iter, ynext
    );
}
