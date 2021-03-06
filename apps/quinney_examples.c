/**
 * \file quinney_examples.c
 * \author Alex Andriati
 * \brief Run simple cases taken from Quinney's book
 *
 * Executable to provice a guide on the workflow to use the library.
 * All the elements from the library are explored, such as optional
 * arguments passed as struct, user-defined system of derivatives,
 * how to cast the optional parameters from void pointer, and usage
 * of auxiliar function for workspace struct handling. The examples
 * were taken from:
 * [1] Douglas Quinney, An introduction to the numerical solution of
 * differential equations, Revised Edition, 1987, cap. 2
 * See examples 2.2.2, 2.2.3, 2.2.4 and 2.5.1. This last one, we can
 * check with the values given in the book for RungeKutta of order 2
 * (aka simple RungeKutta). The multistep order is set to 1 to yield
 * the Euler method providing the lowest accuracy, therefore, we can
 * compare with values in table of example 2.2.3.
 *
 * After build the application, run:
 * $ ./quinney_examples <grid_step>
 * where `grid_step` is a float value and the final grid point is 1
 */

#include <stdio.h>
#include <stdlib.h>
#include "singlestep.h"
#include "multistep.h"


/** \brief Extra parameters for derivatives computation */
struct sys_extra_param{
    double
        coef1,
        coef2,
        coef3;
};


/** \brief Copy values from the first array to the second */
void
rarr_copy_values(unsigned int array_size, Rarray from, Rarray to)
{
    for (unsigned int i = 0; i < array_size; i++) to[i] = from[i];
}


/** \brief System derivatives with 4 equations */
void
sys_der(RealODEInputParameters inp_params, Rarray yprime)
{
    struct sys_extra_param * p = (
            (struct sys_extra_param *) inp_params->extra_args
    );
    double x = inp_params->x;
    Rarray y = inp_params->y;
    yprime[0] = p->coef1 * y[0] + x;
    yprime[1] = p->coef2 * y[1] / (1 + x * x);
    yprime[2] = y[2] * y[2] * x;
    yprime[3] = p->coef3 * y[3];
}


int main(int argc, char * argv[])
{
    int
        nsteps;
    double
        h,
        a[2],
        b[2];
    Rarray
        yrk2,   /* solution of 2nd order RungeKutta */
        yrk4,   /* solution of 4th order RungeKutta */
        yms2,   /* solution of 1st order multistep  */
        ynext;  /* Hold temporary integrator result */
    _RealWorkspaceRK
        wsrk;
    _RealWorkspaceMS
        wsms;
    struct sys_extra_param
        p = { .coef1 = 1.0, .coef2 = 1.0, .coef3 = -1.0 };
    _RealODEInputParameters
        sys_params;

    sys_params.system_size = 4;
    sys_params.extra_args = &p;

    /* These coefficients in multistep yield the Euler's method */
    a[1] = -1.0;
    b[1] =  1.0;

    /* Set input grid step and compute number of steps to finish at 1 */
    h = 0.1;
    if (argc > 2)
    {
        printf("\nMax 1 argument accepted. %d given\n\n", argc - 1);
        exit(EXIT_FAILURE);
    }
    if (argc == 2) sscanf(argv[1], "%lf", &h);
    if (h > 0.5)
    {
        printf("\nMax value for grid step is 0.5 but %.1lf given\n", h);
        exit(EXIT_FAILURE);
    }
    nsteps = ((int) (1.0 + h / 2) / h);

    /* workspace and arrays needed in integrators allocation */
    wsrk.system_size = 4;
    alloc_real_rungekutta_wsarrays(&wsrk);
    wsms.ms_order = 1;
    wsms.system_size = 4;
    alloc_real_multistep_wsarray(&wsms);
    yrk2  = (double *) malloc(wsrk.system_size * sizeof(double));
    yrk4  = (double *) malloc(wsrk.system_size * sizeof(double));
    yms2  = (double *) malloc(wsms.system_size * sizeof(double));
    ynext = (double *) malloc(wsrk.system_size * sizeof(double));

    /* set initial condition for each integrator used */
    for (int i = 0; i < wsrk.system_size; i++)
    {
        yrk2[i] = 1.0;
        yrk4[i] = 1.0;
        yms2[i] = 1.0; // In this case only one step is required for multistep
    }
    /* compute first derivative required for multistep method */
    sys_params.x = 0;
    sys_params.y = yms2;
    sys_der(&sys_params, wsms.prev_der);

    /* header for screen output */
    printf("\nstep x");
    printf("               Euler             ");
    printf("            Rungekutta2          ");
    printf("            Rungekutta4");
    printf("\n----------------------------------------------------");
    printf("------------------------------------------------------");

    for (int i = 0; i <= nsteps; i++)
    {
        printf("\n%6.3lf ", i * h);
        for (int j = 0; j < wsrk.system_size; j++) printf(" %.5lf", yms2[j]);
        real_general_multistep(
                h, i * h, &sys_der, &p, &wsms, yms2, a, b, 0, ynext
        );
        real_set_next_multistep((i + 1) * h, &sys_der, &p, &wsms, yms2, ynext);
        real_rungekutta2(h, i * h, &sys_der, &p, &wsrk, yrk2, ynext);
        printf(" ");
        for (int j = 0; j < wsrk.system_size; j++) printf(" %.5lf", yrk2[j]);
        rarr_copy_values(wsrk.system_size, ynext, yrk2);
        real_rungekutta4(h, i * h, &sys_der, &p, &wsrk, yrk4, ynext);
        printf(" ");
        for (int j = 0; j < wsrk.system_size; j++) printf(" %.5lf", yrk4[j]);
        rarr_copy_values(wsrk.system_size, ynext, yrk4);
    }

    free_real_rungekutta_wsarrays(&wsrk);
    free_real_multistep_wsarray(&wsms);
    free(yrk2);
    free(yrk4);
    free(yms2);
    free(ynext);

    printf("\n\n");
    return 0;
}
