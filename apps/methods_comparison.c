/**
 * \file adams4order_demo.c
 * \author Alex Andriati
 * \brief Run simple cases comparing Adams-Predictor-Corretor with RK
 *
 * Use an analytic solvable equation to compare the convergence of
 * methods. The most accurate methods (RK5 and Adams6) converge in
 * all 12 decimals places shown using step <= 0.005
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "singlestep.h"
#include "multistep.h"


/** \brief Copy values from the first array to the second */
void
rarr_copy_values(unsigned int array_size, Rarray from, Rarray to)
{
    for (unsigned int i = 0; i < array_size; i++) to[i] = from[i];
}


/** \brief System derivatives */
void
sys_der(RealODEInputParameters inp_params, Rarray yprime)
{
    double x = inp_params->x;
    Rarray y = inp_params->y;
    yprime[0] = y[0] - x * x + 1;
}


double
analytic(double x, double y0)
{
    return (y0 - 1) * exp(x) + (1 + x) * (1 + x);
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
        y0[1],
        y00[1],
        yrk4[1],
        yrk5[1],
        yabm4[4],
        yabm6[6],
        yabm_next[1];
    _RealWorkspaceRK
        wsrk;
    _RealWorkspaceMS
        wsms4order,
        wsms6order;

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
    nsteps = ((int) (4.0 + h / 2) / h);

    /* set initial condition */
    y0[0] = 0.5;
    y00[0] = 0.5;

    wsrk.system_size = 1;
    alloc_real_rungekutta_wsarrays(&wsrk);
    wsms4order.ms_order = 4;
    wsms4order.system_size = 1;
    alloc_real_multistep_wsarray(&wsms4order);
    wsms6order.ms_order = 6;
    wsms6order.system_size = 1;
    alloc_real_multistep_wsarray(&wsms6order);

    init_real_multistep(h, &sys_der, NULL, &wsms4order, y0, &real_rungekutta4, yabm4);
    init_real_multistep(h, &sys_der, NULL, &wsms6order, y0, &real_rungekutta5, yabm6);

    printf("\ngrid x     Analytic     RungeKutta4     RungeKutta5     ");
    printf("Adams4step      Adams6step\n");
    printf("-------------------------------------------");
    printf("-------------------------------------------");

    /* print initial condition */
    printf("\n%6.3lf", 0.0);
    printf(" %15.12lf", analytic(0.0, 0.5));
    printf(" %15.12lf", y0[0]);
    printf(" %15.12lf", y00[0]);

    /* print first few rungekutta common steps */
    for (i = 0; i < 3; i++)
    {
        real_rungekutta4(h, i * h, &sys_der, NULL, &wsrk, y0, yrk4);
        rarr_copy_values(wsrk.system_size, yrk4, y0);
        real_rungekutta5(h, i * h, &sys_der, NULL, &wsrk, y00, yrk5);
        rarr_copy_values(wsrk.system_size, yrk5, y00);
        printf("\n%6.3lf", (i + 1) * h);
        printf(" %15.12lf", analytic((i + 1) * h, 0.5));
        for (j = 0; j < wsrk.system_size; j++) printf(" %15.12lf", yrk4[j]);
        for (j = 0; j < wsrk.system_size; j++) printf(" %15.12lf", yrk5[j]);
    }

    /* from this point multistep 4th order adams(PC) is also available */
    for (i = 3; i < 5; i++)
    {
        real_adams4pc(h, i * h, &sys_der, NULL, &wsms4order, yabm4, niter, yabm_next);
        real_set_next_multistep((i + 1) * h, &sys_der, NULL, &wsms4order, yabm4, yabm_next);
        real_rungekutta4(h, i * h, &sys_der, NULL, &wsrk, y0, yrk4);
        rarr_copy_values(wsrk.system_size, yrk4, y0);
        real_rungekutta5(h, i * h, &sys_der, NULL, &wsrk, y00, yrk5);
        rarr_copy_values(wsrk.system_size, yrk5, y00);
        printf("\n%6.3lf", (i + 1) * h);
        printf(" %15.12lf", analytic((i + 1) * h, 0.5));
        for (j = 0; j < wsrk.system_size; j++) printf(" %15.12lf", yrk4[j]);
        for (j = 0; j < wsrk.system_size; j++) printf(" %15.12lf", yrk5[j]);
        for (j = 0; j < wsms4order.system_size; j++)
            printf(" %15.12lf", yabm4[j]);
    }

    /* from this point multistep 4th order adams(PC) is also available */
    for (i = 5; i < nsteps; i++)
    {
        real_adams4pc(h, i * h, &sys_der, NULL, &wsms4order, yabm4, niter, yabm_next);
        real_set_next_multistep((i + 1) * h, &sys_der, NULL, &wsms4order, yabm4, yabm_next);
        real_adams6pc(h, i * h, &sys_der, NULL, &wsms6order, yabm6, niter, yabm_next);
        real_set_next_multistep((i + 1) * h, &sys_der, NULL, &wsms6order, yabm6, yabm_next);
        real_rungekutta4(h, i * h, &sys_der, NULL, &wsrk, y0, yrk4);
        rarr_copy_values(wsrk.system_size, yrk4, y0);
        real_rungekutta5(h, i * h, &sys_der, NULL, &wsrk, y00, yrk5);
        rarr_copy_values(wsrk.system_size, yrk5, y00);
        printf("\n%6.3lf", (i + 1) * h);
        printf(" %15.12lf", analytic((i + 1) * h, 0.5));
        for (j = 0; j < wsrk.system_size; j++) printf(" %15.12lf", yrk4[j]);
        for (j = 0; j < wsrk.system_size; j++) printf(" %15.12lf", yrk5[j]);
        for (j = 0; j < wsms4order.system_size; j++)
            printf(" %15.12lf", yabm4[j]);
        for (j = 0; j < wsms6order.system_size; j++)
            printf(" %15.12lf", yabm6[j]);
    }

    free_real_rungekutta_wsarrays(&wsrk);
    free_real_multistep_wsarray(&wsms4order);
    free_real_multistep_wsarray(&wsms6order);

    printf("\n\n");
    return 0;
}
