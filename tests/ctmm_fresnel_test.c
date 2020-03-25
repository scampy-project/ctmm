/**
 * @file ctmm_test.c
 * @author Angus Bridges
 * @brief A simple test file for the ctmm library.
 * @version 0.1
 * @date 2020-01-23
 * 
 * This file should run all functions included in the ctmm library. If all tests
 * are succesful CTMM TEST SUCCESS will be printed to stdout. Failures should
 * specify the reason if possible. Pass -v as the first argument for more
 * detailed output.
 * 
 * Functions currently tested:
 *      ctmm_create_stack()
 *      ctmm_set_ind()
 *      ctmm_evaluate()
 *          evaluate_prop()
 *              kz()
 *          evaluate_tran()
 *              fresnel_coefs()
 *      ctmm_rts()
 *      ctmm_free_stack()
 *      ctmm_rtps()
 *      ctmm_rtc()
 * Functions to be added:
 * 
 * @copyright Copyright (c) 2019
 * 
 */

#include <stdio.h>
#include <math.h>
#include "ctmm.h"

int main()
{
    double rtps[8];

    double fresnel_n1 = 1;
    double fresnel_n2 = 1.515;
    double fresnel_t_in = 45*M_PI/180;
    double fresnel_t_out = asin(fresnel_n1*sin(fresnel_t_in)/fresnel_n2);
    //Component of wavevector parallel to the interface.
    double k_par = (2*M_PI/633e-9)*sin(fresnel_t_in);

    //Amplitude reflectivity and transmission coefficients (S/P polarisation)
    ctmm_complex fresnel_ctmm[4];
    double fresnel_calc[4];

    ctmm_stack stack;

    stack = ctmm_create_stack(2, 633e-9, fresnel_t_in);

    ctmm_set_ind(stack, 0, fresnel_n1, 0);
    ctmm_set_ind(stack, 1, fresnel_n2, 0);

    ctmm_set_d(stack, 0, 0);
    ctmm_set_d(stack, 1, 0);

    ctmm_evaluate(stack);

    ctmm_rtps(stack, rtps);

    fresnel_coefs(fresnel_ctmm, ctmm_complex_set(fresnel_n1, 0),
        ctmm_complex_set(fresnel_n2, 0), ctmm_complex_set(k_par, 0), 633e-9);

    /**
     * Frensel coefficient calculations see Born & Wolf (6th ed.)
     * p. 40, eqs. (20), (21).
     */
    fresnel_calc[0] = (fresnel_n2*cos(fresnel_t_in) - fresnel_n1*cos(fresnel_t_out))
        /(fresnel_n2*cos(fresnel_t_in) + fresnel_n1*cos(fresnel_t_out));
    fresnel_calc[1] = (2*fresnel_n1*cos(fresnel_t_in))
        /(fresnel_n2*cos(fresnel_t_in) + fresnel_n1*cos(fresnel_t_out));
    fresnel_calc[2] = (fresnel_n1*cos(fresnel_t_in) - fresnel_n2*cos(fresnel_t_out))
        /(fresnel_n1*cos(fresnel_t_in) + fresnel_n2*cos(fresnel_t_out));
    fresnel_calc[3] = (2*fresnel_n1*cos(fresnel_t_in))
        /(fresnel_n1*cos(fresnel_t_in) + fresnel_n2*cos(fresnel_t_out));

    for(int i=0; i<4; i++) {
        printf("calc: %f, ctmm: %f + %f",
            fresnel_calc[i],
            creal(fresnel_ctmm[i]),
            cimag(fresnel_ctmm[i]));
    }

    ctmm_free_stack(stack);

    return 0;
}