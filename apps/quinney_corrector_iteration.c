/**
 * \file quinney_corretor_iteration.c
 * \author Alex Andriati
 * \brief Run specific predictor corretor iterations example
 *
 * Example 2.8.2 taken from:
 * [1] Douglas Quinney, An introduction to the numerical solution of
 * differential equations, Revised Edition, 1987, cap. 2
 * In this example, one of the simplest predictor-corretor scheme is
 * used in the equation y' = y^2 to show the effect of iterations in
 * the corrector part. Only one grid step is advanced.
 *
 * After build the application, run:
 * $ ./quinney_ex282 <grid_step>
 * Show in screen the result of 10 iterations of the correcotr part.
 * Within grid step smaller than 0.14 these 10 iterations are enough
 * to converge all decimal places shown
 *
 * WARNING: The book states that 4th order RungeKutta was used to get
 * the first step required in multistep but the result with all those
 * decimal places presented is only achieved from analytical solution
 * From RK4 one shall obtain 1.11111049 instead of 1.11111111 for 0.1
 * as grid step
 */

#include <stdio.h>
#include <stdlib.h>
#include "ode_multistep.h"


/** \brief System derivatives with 1 equation */
void sys_der(int s, double x, Rarray y, Rarray yprime, void * args)
{
    yprime[0] = y[0] * y[0];
}


int main(int argc, char * argv[])
{
    double
        h,
        a[3],
        b[3],
        aa[3],
        bb[3],
        y0[1],
        y1[1],
        yms[2];
    _RealWorkspaceMS
        wsms;

    /* See Quinney's book predictor example 2.8.1 */
    a[1] = -1.0;
    a[2] =  0.0;
    b[1] =  1.5;
    b[2] = -0.5;
    /* See Quinney's book corrector example 2.8.1 */
    aa[1] = -1.0;
    aa[2] =  0.0;
    bb[0] =  0.5;
    bb[1] =  0.5;
    bb[2] =  0.0;

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

    wsms.ms_order = 2;
    wsms.system_size = 1;
    alloc_real_multistep_array(&wsms);

    y0[0] = 1.0;
    y1[0] = 1.0 / (1 - h);  /* Exact solution instead of RK4 */
    yms[0] = y1[0];
    yms[1] = y0[0];
    sys_der(wsms.system_size, 0.0, y0, &wsms.prev_der[1], NULL);
    sys_der(wsms.system_size, h, y1, &wsms.prev_der[0], NULL);

    printf("\n%6.3lf  %11.8lf", 0.0, y0[0]);
    printf("\n%6.3lf  %11.8lf", h, y1[0]);

    real_general_multistep(h, h, &sys_der, NULL, &wsms, yms, a, b, 0, y0);

    printf("\n%6.3lf  %11.8lf  (predictor)", 2 * h, y0[0]);

    for (int i = 0; i < 10; i++)
    {
        real_general_multistep(
                h, 1 * h, &sys_der, NULL, &wsms, yms, aa, bb, 1, y0
        );
        printf("\n%6.3lf  %11.8lf  ", 2 * h, y0[0]);
        printf("(corrector %2d)", i + 1);
    }

    free_real_multistep_array(&wsms);

    printf("\n\n");
    return 0;
}
