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
    unsigned int fail = 0;

    double rts[4];
    double rtps[8];
    double eps = 1e-13;

    double fresnel_n1 = 1;
    double fresnel_n2 = 1.515;
    double fresnel_t_in = 45*M_PI/180;
    double fresnel_t_out = asin(fresnel_n1*sin(fresnel_t_in)/fresnel_n2);

    //Amplitude reflectivity and transmission coefficients (S/P polarisation)
    double fresnel_rp;
    double fresnel_tp;
    double fresnel_rs;
    double fresnel_ts;

    //Power reflectivity and transmission coefficients
    double fresnel_power_rp;
    double fresnel_power_tp;
    double fresnel_power_rs;
    double fresnel_power_ts;

    double rtsAR[4] = {0.00882229765226, 0.99117770234774, 0.00882229765226, 0.99117770234774};
    double rtsAu[4] = {0.54794668295920, 0.02240833831642, 0.54794668295920, 0.02240833831642};

    ctmm_stack stack;

    stack = ctmm_create_stack(2, 633e-9, fresnel_t_in);

    ctmm_set_ind(stack, 0, fresnel_n1, 0);
    ctmm_set_ind(stack, 1, fresnel_n2, 0);

    ctmm_set_d(stack, 0, 0);
    ctmm_set_d(stack, 1, 0);

    ctmm_evaluate(stack);

    ctmm_rtps(stack, rtps);

    /**
     * Frensel coefficient calculations see Born & Wolf (6th ed.)
     * p. 40, eqs. (20), (21).
     */
    fresnel_rp = (fresnel_n2*cos(fresnel_t_in) - fresnel_n1*cos(fresnel_t_out))
        /(fresnel_n2*cos(fresnel_t_in) + fresnel_n1*cos(fresnel_t_out));
    fresnel_tp = (2*fresnel_n1*cos(fresnel_t_in))
        /(fresnel_n2*cos(fresnel_t_in) + fresnel_n1*cos(fresnel_t_out));
    fresnel_rs = (fresnel_n1*cos(fresnel_t_in) - fresnel_n2*cos(fresnel_t_out))
        /(fresnel_n1*cos(fresnel_t_in) + fresnel_n2*cos(fresnel_t_out));
    fresnel_ts = (2*fresnel_n1*cos(fresnel_t_in))
        /(fresnel_n1*cos(fresnel_t_in) + fresnel_n2*cos(fresnel_t_out));

    fresnel_power_rp = fresnel_rp*fresnel_rp;
    fresnel_power_tp = fresnel_n2*cos(fresnel_t_out)*fresnel_tp*fresnel_tp
        /(fresnel_n1*cos(fresnel_t_in));
    fresnel_power_rs = fresnel_rs*fresnel_rs;
    fresnel_power_ts = fresnel_n2*cos(fresnel_t_out)*fresnel_ts*fresnel_ts
        /(fresnel_n1*cos(fresnel_t_in));

    if (fabs(fresnel_power_rp - rtps[0]) > eps) {
        printf("\nFrensel coefficient failure: coefficient power rp.\n");
        fail = 1;
    } else if (fabs(fresnel_power_tp - rtps[1]) > eps) {
        printf("\nFrensel coefficient failure: coefficient power tp.\n");
        printf("%f, %f", fresnel_power_tp, rtps[1]);
        fail = 1;
    } else if (fabs(fresnel_power_rs - rtps[2]) > eps) {
        printf("\nFrensel coefficient failure: coefficient power rs.\n");
        fail = 1;
    } else if (fabs(fresnel_power_ts - rtps[3]) > eps) {
        printf("\nFrensel coefficient failure: coefficient phase ts.\n");
        fail = 1;
    }

    if (fail) {
        printf("\nFrensel coefficient calculation test failed.\n");
        return 1;
    }

    ctmm_free_stack(stack);
    
    stack = ctmm_create_stack(5, 633e-9, 0);
    
    ctmm_set_ind(stack, 0, 1, 0);
    ctmm_set_ind(stack, 1, 1.38, 0);
    ctmm_set_ind(stack, 2, 1.6, 0);
    ctmm_set_ind(stack, 3, 1.38, 0);
    ctmm_set_ind(stack, 4, 1.5, 0);

    ctmm_set_d(stack, 0, ctmm_get_vwl(stack)/2);
    ctmm_set_d(stack, 1, ctmm_get_vwl(stack)/2);
    ctmm_set_d(stack, 2, ctmm_get_vwl(stack)/2);
    ctmm_set_d(stack, 3, ctmm_get_vwl(stack)/2);
    ctmm_set_d(stack, 4, ctmm_get_vwl(stack)/2);

    ctmm_evaluate(stack);

    ctmm_rts(stack, rts);

    for (unsigned int i=0; i<4; i++) {
        if (fabs(rts[i] - rtsAR[i]) > eps) {
            fail = 1;
            printf("\nAR test failure, coefficient %d.", i);
        }
    }
    if (fail) {
        printf("\nAnti-reflection coating test failed.\n");
        return 2;
    }

    ctmm_free_stack(stack);

    stack = ctmm_create_stack(3, 633e-9, 0);
    
    ctmm_set_ind(stack, 0, 1, 0);
    ctmm_set_ind(stack, 1, 3, 3);
    ctmm_set_ind(stack, 2, 1, 0);

    ctmm_set_d(stack, 0, 0);
    ctmm_set_d(stack, 1, 50e-9);
    ctmm_set_d(stack, 2, 0);

    ctmm_evaluate(stack);

    ctmm_rts(stack, rts);

    for (unsigned int i=0; i<4; i++) {
        if (fabs(rts[i] - rtsAu[i]) > eps) {
            fail = 1;
            printf("\nAu test failure, coefficient %d.", i);
        }
    }
    if (fail) {
        printf("\nGold (metallic) coating test failed.\n");
        return 3;
    }

    ctmm_free_stack(stack);

    return 0;
}