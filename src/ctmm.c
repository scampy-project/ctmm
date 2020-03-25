/**
 * @file ctmm.c
 * @author Angus Bridges
 * @brief An optical transfer matrix modelling library.
 * @version 0.3.0
 * @date 2020-02-11
 *
 * For function documentation see ctmm.h.
 *
 * @copyright Copyright (c) 2020
 * 
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "ctmm.h"
#include "DLL_defines.h"

CTMM_EXPORT double test_func()
{
        printf("\ntest_func\n");
        return 12;
}

/**
 * Defining complex numbers in this way allows C99 complex operations to be
 * performed on ctmm variables in code that will only ever be compiled on a
 * standards compliant compiler, whilst allowing operations internal to ctmm to
 * remain cross platform for the MSVC compiler. Internal complex functions are
 * not optimized, but should compile more or less anywhere.
 */
#if __STDC_VERSION__ >= 199901L
        //check for C99 compiler
        #include <complex.h>

        ctmm_complex ctmm_complex_set(double re, double im)
        {
                return re + im*I;
        }
#else
        //custom definitions for non-standards compliant compilers (ie msvc)
        ctmm_complex ctmm_complex_set(double re, double im)
        {
                ctmm_complex z;
                z.re = re;
                z.im = im;

                return z;
        }
#endif

/**
 * @brief Struct holding stack parameters.
 */
struct ctmm_stack {
        unsigned int nlyrs;

        double vwl; //vacuum wavelength
        double t_in; //theta in, the incident angle
        ctmm_complex kx; //x component of wavevector in zeroth layer

        ctmm_complex *lyr_inds; //layer complex refractive indexes
        double *lyr_ds; //layer thicknesses

        ctmm_matrix tmp1;
        ctmm_matrix tmp2;
        ctmm_matrix tmp3;

        ctmm_matrix matrix; //stack matrix
};

ctmm_complex ctmm_complex_mul_real(ctmm_complex z, double x)
{
        return ctmm_complex_set(creal(z)*x, cimag(z)*x);
}

ctmm_complex ctmm_complex_add_real(ctmm_complex z, double x)
{
        return ctmm_complex_set(creal(z) + x, cimag(z));
}

ctmm_complex ctmm_complex_mul(ctmm_complex z1, ctmm_complex z2)
{
        return ctmm_complex_set(creal(z1)*creal(z2) - cimag(z1)*cimag(z2),
                creal(z1)*cimag(z2) + cimag(z1)*creal(z2));
}

ctmm_complex ctmm_complex_sub(ctmm_complex z1, ctmm_complex z2)
{
        return ctmm_complex_set(creal(z1) - creal(z2), cimag(z1) - cimag(z2));
}

ctmm_complex ctmm_complex_add(ctmm_complex z1, ctmm_complex z2)
{
        return ctmm_complex_set(creal(z1) + creal(z2), cimag(z1) + cimag(z2));
}

ctmm_complex ctmm_complex_div(ctmm_complex z1, ctmm_complex z2)
{
        double denom = creal(z2)*creal(z2) + cimag(z2)*cimag(z2);

        return ctmm_complex_set((creal(z1)*creal(z2) + cimag(z1)*cimag(z2))
        /denom, (cimag(z1)*creal(z2) - creal(z1)*cimag(z2))/denom);
}

double ctmm_complex_abs(ctmm_complex z)
{
        return sqrt(ctmm_complex_abs2(z));
}

ctmm_complex ctmm_complex_sqrt(ctmm_complex z)
{
        return ctmm_complex_mul_real(ctmm_complex_set(cos(0.5
                *ctmm_complex_arg(z)), sin(0.5*ctmm_complex_arg(z))),
                pow(ctmm_complex_abs2(z), 0.25));
}

ctmm_complex ctmm_complex_exp(ctmm_complex z)
{
        return ctmm_complex_set(exp(creal(z))*cos(cimag(z)),
                (exp(creal(z))*sin(cimag(z))));
}

ctmm_complex ctmm_complex_conj(ctmm_complex z)
{
        return ctmm_complex_set(creal(z), -cimag(z));
}

ctmm_complex ctmm_complex_inv(ctmm_complex z)
{
        double denom = creal(z)*creal(z) + cimag(z)*cimag(z);

        return ctmm_complex_set(creal(z)/denom, -cimag(z)/denom);
}

double ctmm_complex_abs2(ctmm_complex z)
{
        return creal(z)*creal(z) + cimag(z)*cimag(z);
}

double ctmm_complex_arg(ctmm_complex z)
{
        return atan2(cimag(z), creal(z));
}

void ctmm_matrix_set(ctmm_matrix *mat, unsigned int row, unsigned int col,
        ctmm_complex val)
{
        mat->data[row*4 + col] = val;
}

ctmm_complex ctmm_matrix_get(ctmm_matrix *mat, unsigned int row,
        unsigned int col)
{
        return mat->data[row*4 + col];
}

void ctmm_matrix_mul(ctmm_matrix *mat1, ctmm_matrix *mat2, ctmm_matrix *res)
{
        //Unrolled naive matrix multiplication. This could probably be
        //put in a loop and optimised out by the compiler anyway.
        res->data[0] = ctmm_complex_add(ctmm_complex_add(ctmm_complex_add(
                ctmm_complex_mul(mat1->data[0], mat2->data[0]),
                ctmm_complex_mul(mat1->data[1], mat2->data[4])),
                ctmm_complex_mul(mat1->data[2], mat2->data[8])),
                ctmm_complex_mul(mat1->data[3], mat2->data[12]));
        res->data[4] = ctmm_complex_add(ctmm_complex_add(ctmm_complex_add(
                ctmm_complex_mul(mat1->data[4], mat2->data[0]),
                ctmm_complex_mul(mat1->data[5], mat2->data[4])),
                ctmm_complex_mul(mat1->data[6], mat2->data[8])),
                ctmm_complex_mul(mat1->data[7], mat2->data[12]));
        res->data[8] = ctmm_complex_add(ctmm_complex_add(ctmm_complex_add(
                ctmm_complex_mul(mat1->data[8], mat2->data[0]),
                ctmm_complex_mul(mat1->data[9], mat2->data[4])),
                ctmm_complex_mul(mat1->data[10], mat2->data[8])),
                ctmm_complex_mul(mat1->data[11], mat2->data[12]));
        res->data[12] = ctmm_complex_add(ctmm_complex_add(ctmm_complex_add(
                ctmm_complex_mul(mat1->data[12], mat2->data[0]),
                ctmm_complex_mul(mat1->data[13], mat2->data[4])),
                ctmm_complex_mul(mat1->data[14], mat2->data[8])),
                ctmm_complex_mul(mat1->data[15], mat2->data[12]));
        res->data[1] = ctmm_complex_add(ctmm_complex_add(ctmm_complex_add(
                ctmm_complex_mul(mat1->data[0], mat2->data[1]),
                ctmm_complex_mul(mat1->data[1], mat2->data[5])),
                ctmm_complex_mul(mat1->data[2], mat2->data[9])),
                ctmm_complex_mul(mat1->data[3], mat2->data[13]));
        res->data[5] = ctmm_complex_add(ctmm_complex_add(ctmm_complex_add(
                ctmm_complex_mul(mat1->data[4], mat2->data[1]),
                ctmm_complex_mul(mat1->data[5], mat2->data[5])),
                ctmm_complex_mul(mat1->data[6], mat2->data[9])),
                ctmm_complex_mul(mat1->data[7], mat2->data[13]));
        res->data[9] = ctmm_complex_add(ctmm_complex_add(ctmm_complex_add(
                ctmm_complex_mul(mat1->data[8], mat2->data[1]),
                ctmm_complex_mul(mat1->data[9], mat2->data[5])),
                ctmm_complex_mul(mat1->data[10], mat2->data[9])),
                ctmm_complex_mul(mat1->data[11], mat2->data[13]));
        res->data[13] = ctmm_complex_add(ctmm_complex_add(ctmm_complex_add(
                ctmm_complex_mul(mat1->data[12], mat2->data[1]),
                ctmm_complex_mul(mat1->data[13], mat2->data[5])),
                ctmm_complex_mul(mat1->data[14], mat2->data[9])),
                ctmm_complex_mul(mat1->data[15], mat2->data[13]));
        res->data[2] = ctmm_complex_add(ctmm_complex_add(ctmm_complex_add(
                ctmm_complex_mul(mat1->data[0], mat2->data[2]),
                ctmm_complex_mul(mat1->data[1], mat2->data[6])),
                ctmm_complex_mul(mat1->data[2], mat2->data[10])),
                ctmm_complex_mul(mat1->data[3], mat2->data[14]));
        res->data[6] = ctmm_complex_add(ctmm_complex_add(ctmm_complex_add(
                ctmm_complex_mul(mat1->data[4], mat2->data[2]),
                ctmm_complex_mul(mat1->data[5], mat2->data[6])),
                ctmm_complex_mul(mat1->data[6], mat2->data[10])),
                ctmm_complex_mul(mat1->data[7], mat2->data[14]));
        res->data[10] = ctmm_complex_add(ctmm_complex_add(ctmm_complex_add(
                ctmm_complex_mul(mat1->data[8], mat2->data[2]),
                ctmm_complex_mul(mat1->data[9], mat2->data[6])),
                ctmm_complex_mul(mat1->data[10], mat2->data[10])),
                ctmm_complex_mul(mat1->data[11], mat2->data[14]));
        res->data[14] = ctmm_complex_add(ctmm_complex_add(ctmm_complex_add(
                ctmm_complex_mul(mat1->data[12], mat2->data[2]),
                ctmm_complex_mul(mat1->data[13], mat2->data[6])),
                ctmm_complex_mul(mat1->data[14], mat2->data[10])),
                ctmm_complex_mul(mat1->data[15], mat2->data[14]));
        res->data[3] = ctmm_complex_add(ctmm_complex_add(ctmm_complex_add(
                ctmm_complex_mul(mat1->data[0], mat2->data[3]),
                ctmm_complex_mul(mat1->data[1], mat2->data[7])),
                ctmm_complex_mul(mat1->data[2], mat2->data[11])),
                ctmm_complex_mul(mat1->data[3], mat2->data[15]));
        res->data[7] = ctmm_complex_add(ctmm_complex_add(ctmm_complex_add(
                ctmm_complex_mul(mat1->data[4], mat2->data[3]),
                ctmm_complex_mul(mat1->data[5], mat2->data[7])),
                ctmm_complex_mul(mat1->data[6], mat2->data[11])),
                ctmm_complex_mul(mat1->data[7], mat2->data[15]));
        res->data[11] = ctmm_complex_add(ctmm_complex_add(ctmm_complex_add(
                ctmm_complex_mul(mat1->data[8], mat2->data[3]),
                ctmm_complex_mul(mat1->data[9], mat2->data[7])),
                ctmm_complex_mul(mat1->data[10], mat2->data[11])),
                ctmm_complex_mul(mat1->data[11], mat2->data[15]));
        res->data[15] = ctmm_complex_add(ctmm_complex_add(ctmm_complex_add(
                ctmm_complex_mul(mat1->data[12], mat2->data[3]),
                ctmm_complex_mul(mat1->data[13], mat2->data[7])),
                ctmm_complex_mul(mat1->data[14], mat2->data[11])),
                ctmm_complex_mul(mat1->data[15], mat2->data[15]));
}

ctmm_stack ctmm_create_stack(unsigned int nlyrs, double vwl, double t_in)
{
        ctmm_stack stack;
        stack = (ctmm_stack)malloc(sizeof(struct ctmm_stack));

        stack->vwl = vwl;
        stack->t_in = t_in;

        stack->nlyrs = nlyrs;
        stack->lyr_inds =
                (ctmm_complex*) malloc(nlyrs * sizeof * stack->lyr_inds);
        stack->lyr_ds = (double*) malloc(nlyrs * sizeof * stack->lyr_ds);

        if (nlyrs < 1) {
                fprintf(stderr,
                        "\nError, stack initialised with no layers.\n");
        }

        return stack;
}

void ctmm_free_stack(ctmm_stack stack)
{
        free(stack->lyr_inds);
        free(stack->lyr_ds);
        free(stack);
}

void ctmm_set_t_in(ctmm_stack stack, double t_in)
{
        stack->t_in = t_in;
}

void ctmm_set_ind(ctmm_stack stack, unsigned int lyr_n, double x, double y)
{
        stack->lyr_inds[lyr_n] = ctmm_complex_set(x, y);
}

void ctmm_set_d(ctmm_stack stack, unsigned int lyr_n, double d)
{
        stack->lyr_ds[lyr_n] = d;
}

ctmm_complex kz(ctmm_complex n, ctmm_complex k, double vwl)
{
        n = ctmm_complex_mul_real(n, 2*M_PI/vwl);
        n = ctmm_complex_mul(n, n);
        k = ctmm_complex_mul(k, k);
        n = ctmm_complex_sub(n, k);
        n = ctmm_complex_sqrt(n);

        return n;
}

ctmm_complex ctmm_get_ind(ctmm_stack stack, unsigned int lyr_n)
{
        return stack->lyr_inds[lyr_n];
}

double ctmm_get_d(ctmm_stack stack, unsigned int lyr_n)
{
        return stack->lyr_ds[lyr_n];
}

double ctmm_get_vwl(ctmm_stack stack)
{
        return stack->vwl;
}

double ctmm_get_t_in(ctmm_stack stack)
{
        return stack->t_in;
}

unsigned int ctmm_get_nlyrs(ctmm_stack stack)
{
        return stack->nlyrs;
}

ctmm_matrix *ctmm_get_matrix(ctmm_stack stack)
{
        return &stack->matrix;
}

void evaluate_prop(ctmm_matrix *pmat, double d, ctmm_complex n,
                                        ctmm_complex k, double vwl)
{
        k = kz(n, k, vwl);

        ctmm_complex expn;

        expn = ctmm_complex_set(-d, 0);

        expn = ctmm_complex_mul(expn, k);
        expn = ctmm_complex_mul(expn, I);

        ctmm_matrix_set(pmat, 0, 1, ctmm_complex_set(0, 0));
        ctmm_matrix_set(pmat, 0, 2, ctmm_complex_set(0, 0));
        ctmm_matrix_set(pmat, 0, 3, ctmm_complex_set(0, 0));

        ctmm_matrix_set(pmat, 1, 0, ctmm_complex_set(0, 0));
        ctmm_matrix_set(pmat, 1, 2, ctmm_complex_set(0, 0));
        ctmm_matrix_set(pmat, 1, 3, ctmm_complex_set(0, 0));

        ctmm_matrix_set(pmat, 2, 1, ctmm_complex_set(0, 0));
        ctmm_matrix_set(pmat, 2, 0, ctmm_complex_set(0, 0));
        ctmm_matrix_set(pmat, 2, 3, ctmm_complex_set(0, 0));

        ctmm_matrix_set(pmat, 3, 1, ctmm_complex_set(0, 0));
        ctmm_matrix_set(pmat, 3, 2, ctmm_complex_set(0, 0));
        ctmm_matrix_set(pmat, 3, 0, ctmm_complex_set(0, 0));

        ctmm_matrix_set(pmat, 0, 0, ctmm_complex_exp(expn));
        ctmm_matrix_set(pmat, 1, 1,
                ctmm_complex_exp(ctmm_complex_mul_real(expn, -1)));
        ctmm_matrix_set(pmat, 2, 2, ctmm_complex_exp(expn));
        ctmm_matrix_set(pmat, 3, 3,
                ctmm_complex_exp(ctmm_complex_mul_real(expn, -1)));
}

void fresnel_coefs(ctmm_complex *coefs, ctmm_complex n0, ctmm_complex n1,
                                        ctmm_complex k, double vwl)
{
        ctmm_complex kz0, kz1;
        kz0 = kz(n0, k, vwl);
        kz1 = kz(n1, k, vwl);

        coefs[0] = ctmm_complex_mul(kz0, ctmm_complex_mul(n1, n1));
        coefs[1] = ctmm_complex_mul(kz1, ctmm_complex_mul(n0, n0));
        coefs[2] = ctmm_complex_sub(coefs[0], coefs[1]);
        coefs[3] = ctmm_complex_add(coefs[0], coefs[1]);
        coefs[0] = ctmm_complex_div(coefs[2], coefs[3]); //rP (TM)
        coefs[1] = ctmm_complex_add_real(coefs[0], 1);
        coefs[1] = ctmm_complex_mul(coefs[1], n0);
        coefs[1] = ctmm_complex_div(coefs[1], n1); //tP (TM)

        coefs[2] = ctmm_complex_sub(kz0, kz1);
        coefs[3] = ctmm_complex_add(kz0, kz1);
        coefs[2] = ctmm_complex_div(coefs[2], coefs[3]); //rS (TE)
        coefs[3] = ctmm_complex_add_real(coefs[2], 1); //tS (TE)
}

void evaluate_tran(ctmm_matrix *tmat, ctmm_complex k, ctmm_complex n0,
                                        ctmm_complex n1, double vwl)
{
        ctmm_complex coefs[4];
        fresnel_coefs(coefs, k, n0, n1, vwl);

        coefs[1] = ctmm_complex_inv(coefs[1]);
        coefs[3] = ctmm_complex_inv(coefs[3]);

        ctmm_matrix_set(tmat, 0, 2, ctmm_complex_set(0, 0));
        ctmm_matrix_set(tmat, 0, 3, ctmm_complex_set(0, 0));
        ctmm_matrix_set(tmat, 1, 2, ctmm_complex_set(0, 0));
        ctmm_matrix_set(tmat, 1, 3, ctmm_complex_set(0, 0));
        ctmm_matrix_set(tmat, 2, 0, ctmm_complex_set(0, 0));
        ctmm_matrix_set(tmat, 2, 1, ctmm_complex_set(0, 0));
        ctmm_matrix_set(tmat, 3, 0, ctmm_complex_set(0, 0));
        ctmm_matrix_set(tmat, 3, 1, ctmm_complex_set(0, 0));

        ctmm_matrix_set(tmat, 0, 0, coefs[1]);
        ctmm_matrix_set(tmat, 1, 1, coefs[1]);
        ctmm_matrix_set(tmat, 2, 2, coefs[3]);
        ctmm_matrix_set(tmat, 3, 3, coefs[3]);

        ctmm_matrix_set(tmat, 1, 0, ctmm_complex_mul(coefs[0], coefs[1]));
        ctmm_matrix_set(tmat, 0, 1, ctmm_complex_mul(coefs[0], coefs[1]));

        ctmm_matrix_set(tmat, 2, 3, ctmm_complex_mul(coefs[2], coefs[3]));
        ctmm_matrix_set(tmat, 3, 2, ctmm_complex_mul(coefs[2], coefs[3]));
}

void ctmm_evaluate(ctmm_stack stack)
{
        stack->kx = ctmm_complex_set(creal(stack->lyr_inds[0])*sin(stack->t_in)
                *2*M_PI/stack->vwl, 0);

        if (stack->nlyrs > 1) {
                evaluate_prop(&stack->matrix, stack->lyr_ds[0],
                        stack->lyr_inds[0], stack->kx, stack->vwl);

                for (unsigned int i=1; i<stack->nlyrs; i++) {
                        evaluate_prop(&stack->tmp1, stack->lyr_ds[i],
                                stack->lyr_inds[i], stack->kx, stack->vwl);
                        evaluate_tran(&stack->tmp2, stack->lyr_inds[i-1],
                                stack->lyr_inds[i], stack->kx, stack->vwl);
                        ctmm_matrix_mul(&stack->matrix, &stack->tmp2,
                                &stack->tmp3);
                        ctmm_matrix_mul(&stack->tmp3, &stack->tmp1,
                                &stack->matrix);
                }
        } else if (stack->nlyrs == 1) {
                evaluate_prop(&stack->matrix, stack->lyr_ds[0],
                        stack->lyr_inds[0], stack->kx, stack->vwl);
        }
}

void ctmm_rts(ctmm_stack stack, double *coefs)
{
        double t_fin, t_mult;

        t_fin = atan(creal(stack->kx)
                /creal(kz(stack->lyr_inds[stack->nlyrs - 1],
                stack->kx, stack->vwl)));
        t_mult = (cos(t_fin)*creal(stack->lyr_inds[stack->nlyrs - 1]))
                /(cos(stack->t_in)*creal(stack->lyr_inds[0]));

        coefs[0] = ctmm_complex_abs2(ctmm_complex_div(ctmm_matrix_get(
                &stack->matrix, 1, 0), ctmm_matrix_get(&stack->matrix, 0, 0)));
        coefs[1] = t_mult*ctmm_complex_abs2(ctmm_complex_inv(
                ctmm_matrix_get(&stack->matrix, 0, 0)));
        coefs[2] = ctmm_complex_abs2(ctmm_complex_div(ctmm_matrix_get(
                &stack->matrix, 3, 2), ctmm_matrix_get(&stack->matrix, 2, 2)));
        coefs[3] = t_mult*ctmm_complex_abs2(ctmm_complex_inv(
                ctmm_matrix_get(&stack->matrix, 2, 2)));
}

void ctmm_rtps(ctmm_stack stack, double *coefs)
{
        double t_fin, t_mult;

        t_fin = atan(creal(stack->kx)
                /creal(kz(stack->lyr_inds[stack->nlyrs - 1],
                stack->kx, stack->vwl)));
        t_mult = (cos(t_fin)*creal(stack->lyr_inds[stack->nlyrs - 1]))
                /(cos(stack->t_in)*creal(stack->lyr_inds[0]));

        coefs[0] = ctmm_complex_abs2(ctmm_complex_div(ctmm_matrix_get(
                &stack->matrix, 1, 0), ctmm_matrix_get(&stack->matrix, 0, 0)));
        coefs[1] = t_mult*ctmm_complex_abs2(ctmm_complex_inv(
                ctmm_matrix_get(&stack->matrix, 0, 0)));
        coefs[2] = ctmm_complex_abs2(ctmm_complex_div(ctmm_matrix_get(
                &stack->matrix, 3, 2), ctmm_matrix_get(&stack->matrix, 2, 2)));
        coefs[3] = t_mult*ctmm_complex_abs2(ctmm_complex_inv(
                ctmm_matrix_get(&stack->matrix, 2, 2)));
        coefs[4] = ctmm_complex_arg(ctmm_complex_div(ctmm_matrix_get(
                &stack->matrix, 1, 0), ctmm_matrix_get(&stack->matrix, 0, 0)));
        coefs[5] = ctmm_complex_arg(ctmm_complex_inv(
                ctmm_matrix_get(&stack->matrix, 0, 0)));
        coefs[6] = ctmm_complex_arg(ctmm_complex_div(ctmm_matrix_get(
                &stack->matrix, 3, 2), ctmm_matrix_get(&stack->matrix, 2, 2)));
        coefs[7] = ctmm_complex_arg(ctmm_complex_inv(
                ctmm_matrix_get(&stack->matrix, 2, 2)));
}

void ctmm_rtc(ctmm_stack stack, ctmm_complex *coefs)
{
        coefs[0] = ctmm_complex_div(ctmm_matrix_get(&stack->matrix, 1, 0),
                ctmm_matrix_get(&stack->matrix, 0, 0));
        coefs[1] = ctmm_complex_inv(ctmm_matrix_get(&stack->matrix, 0, 0));
        coefs[2] = ctmm_complex_div(ctmm_matrix_get(&stack->matrix, 3, 2),
                ctmm_matrix_get(&stack->matrix, 2, 2));
        coefs[3] = ctmm_complex_inv(ctmm_matrix_get(&stack->matrix, 2, 2));
}