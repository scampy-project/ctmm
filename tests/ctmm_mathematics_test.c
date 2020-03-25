/**
 * @file ctmm_test.c
 * @author Angus Bridges
 * @brief A test for the mathematica functions implemented in the ctmm library.
 * @version 0.1
 * @date 2020-01-20
 *
 * This file should run all mathematical functions included in the ctmm library.
 * If all tests are succesful CTMM MATHEMATICS TEST SUCCESS will be printed to
 * stdout. Failures should specify the reason if possible. Pass -v as the first
 * argument for more detailed output.
 *
 * Functions currently tested:
 *      ctmm_matrix_set()
 *      ctmm_complex_set()
 *      ctmm_matrix_get()
 *      ctmm_matrix_mul()
 *          ctmm_complex_add()
 *          ctmm_complex_mul()
 *      ctmm_complex_exp()
 *      ctmm_complex_mul_real()
 *      ctmm_complex_div()
 * Functions to be added:
 *      ctmm_complex_add_real()
 *      ctmm_complex_sub()
 *      ctmm_complex_sqrt()
 *      ctmm_complex_conj()
 *      ctmm_complex_inv()
 *      ctmm_complex_abs2()
 *      ctmm_complex_arg()
 *
 * To compile this test run...
 *
 * @copyright Copyright (c) 2019
 *
 */

#include <stdio.h>
#include "ctmm.h"

int main()
{
    unsigned int fail = 0;

    double eps = 1e-13;

    double element_test_values_real[16] =
        {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};
    double element_test_values_imag[16] =
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    double multiplication_test_values_real[16] =
        {0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, -3, -3, -3, -3};
    double multiplication_test_values_imag[16] =
        {0, 1, 2, 3, -0, -1, -2, -3, 0, 1, 2, 3, 0, 1, 2, 3};

    ctmm_matrix mat1;
    ctmm_matrix mat2;
    ctmm_matrix mat3;

    ctmm_complex cval1;

    ctmm_matrix_set(&mat1, 0, 0, ctmm_complex_set(1, 0));
    ctmm_matrix_set(&mat1, 0, 1, ctmm_complex_set(0, 0));
    ctmm_matrix_set(&mat1, 0, 2, ctmm_complex_set(0, 0));
    ctmm_matrix_set(&mat1, 0, 3, ctmm_complex_set(0, 0));

    ctmm_matrix_set(&mat1, 1, 0, ctmm_complex_set(0, 0));
    ctmm_matrix_set(&mat1, 1, 1, ctmm_complex_set(1, 0));
    ctmm_matrix_set(&mat1, 1, 2, ctmm_complex_set(0, 0));
    ctmm_matrix_set(&mat1, 1, 3, ctmm_complex_set(0, 0));

    ctmm_matrix_set(&mat1, 2, 0, ctmm_complex_set(0, 0));
    ctmm_matrix_set(&mat1, 2, 1, ctmm_complex_set(0, 0));
    ctmm_matrix_set(&mat1, 2, 2, ctmm_complex_set(1, 0));
    ctmm_matrix_set(&mat1, 2, 3, ctmm_complex_set(0, 0));

    ctmm_matrix_set(&mat1, 3, 0, ctmm_complex_set(0, 0));
    ctmm_matrix_set(&mat1, 3, 1, ctmm_complex_set(0, 0));
    ctmm_matrix_set(&mat1, 3, 2, ctmm_complex_set(0, 0));
    ctmm_matrix_set(&mat1, 3, 3, ctmm_complex_set(1, 0));

    for (int i=0; i<4; i++) {
            for (int j=0; j<4; j++) {
                if (fabs(creal(ctmm_matrix_get(&mat1, i, j))
                    - element_test_values_real[i*4 + j]) > eps) {
                        fail = 1;
                        printf("\n\n\t\t\tFail element (%d, %d), real part.\n",
                            i, j);
                } else if (fabs(cimag(ctmm_matrix_get(&mat1, i, j))
                    - element_test_values_imag[i*4 + j]) > eps) {
                        fail = 1;
                        printf("\n\n\t\t\tFail element (%d, %d), imaginary part.\n",
                            i, j);
                }
            }
    }

    if (fail) {
        printf("Failure on matrix set/get test.");
        return 1;
    }

    //Matrix multiplication test also tests complex addition and multiplication.

    ctmm_matrix_set(&mat2, 0, 0, ctmm_complex_set(0, 0));
    ctmm_matrix_set(&mat2, 0, 1, ctmm_complex_set(0, 1));
    ctmm_matrix_set(&mat2, 0, 2, ctmm_complex_set(0, 2));
    ctmm_matrix_set(&mat2, 0, 3, ctmm_complex_set(0, 3));

    ctmm_matrix_set(&mat2, 1, 0, ctmm_complex_set(1, -0));
    ctmm_matrix_set(&mat2, 1, 1, ctmm_complex_set(1, -1));
    ctmm_matrix_set(&mat2, 1, 2, ctmm_complex_set(1, -2));
    ctmm_matrix_set(&mat2, 1, 3, ctmm_complex_set(1, -3));

    ctmm_matrix_set(&mat2, 2, 0, ctmm_complex_set(2, 0));
    ctmm_matrix_set(&mat2, 2, 1, ctmm_complex_set(2, 1));
    ctmm_matrix_set(&mat2, 2, 2, ctmm_complex_set(2, 2));
    ctmm_matrix_set(&mat2, 2, 3, ctmm_complex_set(2, 3));

    ctmm_matrix_set(&mat2, 3, 0, ctmm_complex_set(-3, 0));
    ctmm_matrix_set(&mat2, 3, 1, ctmm_complex_set(-3, 1));
    ctmm_matrix_set(&mat2, 3, 2, ctmm_complex_set(-3, 2));
    ctmm_matrix_set(&mat2, 3, 3, ctmm_complex_set(-3, 3));

    ctmm_matrix_mul(&mat1, &mat2, &mat3);

    for (int i=0; i<4; i++) {
            for (int j=0; j<4; j++) {
                if (fabs(creal(ctmm_matrix_get(&mat3, i, j))
                    - multiplication_test_values_real[i*4 + j]) > eps) {
                        fail = 1;
                        printf("\n\n\t\t\tFail element (%d, %d), real part.\n",
                            i, j);
                } else if (fabs(cimag(ctmm_matrix_get(&mat3, i, j))
                    - multiplication_test_values_imag[i*4 + j]) > eps) {
                        fail = 1;
                        printf("\n\n\t\t\tFail element (%d, %d), imaginary part.\n",
                            i, j);
                }
            }
    }

    if (fail) {
        printf("Failure on matrix multiplication test.");
        return 1;
    }

    cval1 = ctmm_complex_exp(ctmm_complex_mul_real(I, M_PI));

    if (fabs(creal(cval1) + 1) > eps) {
        fail = 1;
    } else if (fabs(cimag(cval1)) > eps) {
        fail = 1;
    }

    if (fail) {
        printf("Failure on exponention test.");
        return 1;
    }

    cval1 = ctmm_complex_div(ctmm_complex_set(1, 3), ctmm_complex_set(1, 2));

    if (fabs(creal(cval1) - 1.4) > eps) {
        fail = 1;
    } else if (fabs(cimag(cval1) - 0.2) > eps) {
        fail = 1;
    }

    if (fail) {
        printf("Failure on division test.");
        return 1;
    }

    cval1 = ctmm_complex_set(0, 1);
    cval1 = ctmm_complex_sqrt(cval1);

    if (fabs(creal(cval1) - sqrt(0.5)) > eps) {
        fail = 1;
    } else if (fabs(cimag(cval1) - sqrt(0.5)) > eps) {
        fail = 1;
    }
    if (fail) {
        printf("Failure on sqrt test.");
        return 1;
    }

    return 0;
}