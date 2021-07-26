/**
 * \file adams4order_demo.c
 * \author Alex Andriati
 * \brief Run simple cases comparing Adams-Predictor-Corretor with RK4
 *
 * This example show how to use 4th order Adams-Predictor-Corretor
 * Some of the examples were taken from:
 * [1]  Douglas Quinney, An introduction to the numerical solution
 * of differential equations, Revised Edition, 1987, cap. 2
 * See examples 2.2.2, 2.2.3, 2.2.4 and 2.5.1. Another example was
 * taken from:
 * [2]  https://labmathdu.wordpress.com/solving-ivp-by-adams-fourth-order-predictor-corrector-method/
 * which correspond to the last equation in the system defined. The
 * values can be compared up to 7 decimal places with ref. [2] due
 * to REAL precision of Fortran that is smaller than C double type
 *
 * After build the application, run:
 * $ ./adams4order_demo <grid_step> <corrector_iterations>
 * where the final propagation point is 1. Use this example to see
 * the interplay of <grid_step> and <corrector_iterations> needed
 * for convergence and note that for sufficient small <grid_step>
 * (typically 0.005) only one iteration is needed
 */

#include <stdio.h>
#include <stdlib.h>
#include "ode_singlestep.h"
#include "ode_multistep.h"


/** \brief System derivatives with 4 equations */
void sys_der(RealODEInputParameters inp_params, Rarray yprime)
{
    double x = inp_params->x;
    Rarray y = inp_params->y;
    yprime[0] = y[0] + x;
    yprime[1] = y[1] / (1 + x * x);
    yprime[2] = y[2] * y[2] * x;
    yprime[3] = y[3] - x * x + 1; // taken from ref [2]
}


int main(int argc, char * argv[])
{
    int
        i,
        j,
        niter,
        nsteps;
    double
        h,
        y0[4],
        yrk4[4],
        yabm[16],
        yabm_next[4];
    _RealWorkspaceRK
        wsrk;
    _RealWorkspaceMS
        wsms;
    _RealODEInputParameters
        sys_params;

    if (argc > 3)
    {
        printf("\nMax 2 arguments accepted. %d given\n\n", argc - 1);
        exit(EXIT_FAILURE);
    }

    h = 0.1;
    niter = 1;
    switch (argc - 1)
    {
        case 2:
            sscanf(argv[1], "%lf", &h);
            sscanf(argv[2], "%d", &niter);
            break;
        case 1:
            sscanf(argv[1], "%lf", &h);
            break;
    }

    /* default final propagation grid point */
    nsteps = ((int) (1.0 + h / 2) / h);

    /* set initial condition */
    y0[0] = 1.0; y0[1] = 1.0; y0[2] = 1.0; y0[3] = 0.5;

    wsrk.system_size = 4;
    alloc_real_rkwsarrays(&wsrk);
    wsms.ms_order = 4;
    wsms.system_size = 4;
    alloc_real_multistep_array(&wsms);
    sys_params.system_size = 4;
    sys_params.extra_args = NULL;

    rarr_copy_values(wsrk.system_size, y0, yrk4);
    for (i = 0; i < wsms.ms_order; i++)
    {
        rarr_copy_values(wsrk.system_size, yrk4, y0);
        sys_params.x = i * h;
        sys_params.y = yrk4;
        j = (wsms.ms_order - 1 - i) * wsms.system_size;
        sys_der(&sys_params, &wsms.prev_der[j]);
        rarr_copy_values(wsms.system_size, yrk4, &yabm[j]);
        printf("\n%6.3lf", i * h);
        for (j = 0; j < wsrk.system_size; j++) printf(" %11.8lf", yrk4[j]);
        for (j = 0; j < wsms.system_size; j++) printf(" %17.14lf", yrk4[j]);
        real_rungekutta4(h, i * h, &sys_der, NULL, &wsrk, y0, yrk4);
    }

    /* recede the extra step unduly advanced in rk4 */
    rarr_copy_values(wsrk.system_size, y0, yrk4);

    for (i = 3; i < nsteps; i++)
    {
        real_adams4pc(h, i * h, &sys_der, NULL, &wsms, yabm, niter, yabm_next);
        real_set_next_step((i + 1) * h, &sys_der, NULL, &wsms, yabm, yabm_next);
        real_rungekutta4(h, i * h, &sys_der, NULL, &wsrk, y0, yrk4);
        rarr_copy_values(wsrk.system_size, yrk4, y0);
        printf("\n%6.3lf", (i + 1) * h);
        for (j = 0; j < wsrk.system_size; j++)
            printf(" %11.8lf", yrk4[j]);
        for (j = 0; j < wsms.system_size; j++)
            printf(" %17.14lf", yabm_next[j]);
    }

    free_real_rkwsarrays(&wsrk);
    free_real_multistep_array(&wsms);

    printf("\n\n");
    return 0;
}
