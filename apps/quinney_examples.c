/**
 * \file quinney_examples.c
 * \author Alex Andriati
 * \brief Run simple cases taken from Quinney's book
 *
 * Executable to provice a guide on the workflow to use the library.
 * All the elements from the library are explored, such as optional
 * arguments passed as struct, user-defined function for system derivatives,
 * how to cast the optional parameters from void pointer, and usage
 * of auxiliar function for workspace struct handling. The examples
 * were taken from:
 * [1] Douglas Quinney, An introduction to the numerical solution of
 * differential equations, Revised Edition, 1987, cap. 2
 * See examples 2.2.2, 2.2.3, 2.2.4 and 2.5.1. This last one, we can
 * check with the values given in the 2.5.1 for simple RungeKutta of
 * order 2 (aka simple RungeKutta)
 *
 * After build the application, run:
 * ./quinney_examples <grid_step>
 * where `grid_step` is a float value and the (max) final point is 1
 */

#include <stdio.h>
#include <stdlib.h>
#include "ode_singlestep.h"


/** \brief Extra parameters for derivatives computation */
struct sys_param_set{
    double
        coef1,
        coef2,
        coef3;
};


/** \brief Auxiliar function to copy values and update initial input **/
void copyvalues(int s, Rarray from, Rarray to)
{
    for (int i = 0; i < s; i++) to[i] = from[i];
}


/** \brief System derivatives computation with 4 equations */
void sys_der(int s, double x, Rarray y, Rarray yprime, void * args)
{
    struct sys_param_set * p = (struct sys_param_set *) args;
    yprime[0] = p->coef1 * y[0] + x;
    yprime[1] = p->coef2 * y[1] / (1 + x * x);
    yprime[2] = y[2] * y[2] * x;
    yprime[3] = p->coef3 * y[s - 1];
}


int main(int argc, char * argv[])
{
    int
        nsteps;
    double
        h;
    struct sys_param_set
        p = { .coef1 = 1.0, .coef2 = 1.0, .coef3 = -1.0 };
    Rarray
        yrk2,
        yrk4,
        ynext;
    _RealWorkspaceRK
        wsrk;

    if (argc > 2)
    {
        printf("\nMax 1 argument accepted. %d given\n\n", argc - 1);
        exit(EXIT_FAILURE);
    }
    if (argc == 2) sscanf(argv[1], "%lf", &h);
    else           h = 0.1;

    nsteps = ((int) 1.0005 / h);

    wsrk.system_size = 4;
    alloc_real_rkwsarrays(&wsrk);
    yrk2  = (double *) malloc(wsrk.system_size * sizeof(double));
    yrk4  = (double *) malloc(wsrk.system_size * sizeof(double));
    ynext = (double *) malloc(wsrk.system_size * sizeof(double));

    for (int i = 0; i < wsrk.system_size; i++)
    {
        yrk2[i] = 1.0;
        yrk4[i] = 1.0;
    }

    for (int i = 0; i <= nsteps; i++)
    {
        real_rungekutta2(h, i * h, &sys_der, &p, &wsrk, yrk2, ynext);
        printf("\n%8.3lf", i * h);
        for (int j = 0; j < wsrk.system_size; j++) printf(" %.5lf", yrk2[j]);
        copyvalues(wsrk.system_size, ynext, yrk2);
        real_rungekutta4(h, i * h, &sys_der, &p, &wsrk, yrk4, ynext);
        for (int j = 0; j < wsrk.system_size; j++) printf(" %.5lf", yrk4[j]);
        copyvalues(wsrk.system_size, ynext, yrk4);
    }

    free_real_rkwsarrays(&wsrk);
    free(yrk2);
    free(yrk4);
    free(ynext);

    printf("\n\n");
    return 0;
}
