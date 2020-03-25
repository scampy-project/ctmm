/**
 * @file ctmm.h
 * @author Angus Bridges
 * @brief An optical transfer matrix modelling library.
 * @version 0.3.0
 * @date 2020-02-11
 *
 * ctmm is an optical thin film transfer matrix modelling library written in C.
 * It is primarily designed to provide a lightweight and efficient backend for
 * tspy, a python package for optical network modelling. A 4x4 transfer matrix
 * methodology is implemented, treating both polarisations simultaneously.
 * Support for birefringent materials may be added in the future.
 *
 * @copyright Copyright (c) 2020
 */

#ifndef _CTMM_H_
#define _CTMM_H_

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "DLL_defines.h"

#if __STDC_VERSION__ >= 199901L
    //check for C99 compiler
    #include <complex.h>
    typedef double _Complex ctmm_complex;
#else
    //custom definitions for non-standards compliant compilers (ie msvc)

    typedef struct {
        double re, im;
    } ctmm_complex;

    //redefine macros from <complex.h> for compatability with C99 complex.
    #define creal(z) ((z).re)
    #define cimag(z) ((z).im)

    #ifndef im_i
    const ctmm_complex im_i = {0, 1};
    #define I im_i
    #endif
#endif

/**
 * @brief returns ctmm_complex number re + im*I.
 * 
 * @param re real part of complex number re + im*I
 * @param im imaginary part of complex number re + im*I
 * 
 * The exact definition of this function varies depending on if a C99 compliant
 * compiler is available.
 */
CTMM_EXPORT ctmm_complex ctmm_complex_set(double re, double im);

/**
 * Fixed size 4x4 complex matrix.
 */
CTMM_EXPORT typedef struct ctmm_matrix {
    ctmm_complex data[16];
} ctmm_matrix;

/**
 * @breif multiplies ctmm_complex by real number
 * 
 * @param z complex number
 * @param x real number
 * 
 * @return ctmm_complex z*x
 */
CTMM_EXPORT ctmm_complex ctmm_complex_mul_real(ctmm_complex z, double x);

/**
 * @breif add ctmm_complex to real number
 * 
 * @param z complex number
 * @param x real number
 * 
 * @return ctmm_complex z + x
 */
CTMM_EXPORT ctmm_complex ctmm_complex_add_real(ctmm_complex z, double x);

/**
 * @breif multiplies two ctmm_complex numbers
 * 
 * @param z1 complex number
 * @param z2 complex number
 * 
 * @return ctmm_complex z1*z2
 */
CTMM_EXPORT ctmm_complex ctmm_complex_mul(ctmm_complex z1, ctmm_complex z2);

/**
 * @breif subtracts one ctmm_complex from another
 * 
 * @param z1 complex number
 * @param z2 complex number
 * 
 * @return ctmm_complex z1 - z2
 */
CTMM_EXPORT ctmm_complex ctmm_complex_sub(ctmm_complex z1, ctmm_complex z2);

/**
 * @breif adds one ctmm_complex to another
 * 
 * @param z1 complex number
 * @param z2 complex number
 * 
 * @return ctmm_complex z1 + z2
 */
CTMM_EXPORT ctmm_complex ctmm_complex_add(ctmm_complex z1, ctmm_complex z2);

/**
 * @breif divides one ctmm_complex by another
 * 
 * @param z1 complex number
 * @param z2 complex number
 * 
 * @return ctmm_complex z1/z2
 */
CTMM_EXPORT ctmm_complex ctmm_complex_div(ctmm_complex z1, ctmm_complex z2);

/**
 * @breif calculates square root of ctmm_complex
 * 
 * @param z complex number
 * 
 * @return ctmm_complex sqrt(z)
 */
CTMM_EXPORT ctmm_complex ctmm_complex_sqrt(ctmm_complex z);

/**
 * @breif calculates e to the power ctmm_complex
 * 
 * @param z complex number
 * 
 * @return ctmm_complex e^z
 */
CTMM_EXPORT ctmm_complex ctmm_complex_exp(ctmm_complex z);

/**
 * @breif calculates complex conjugate of ctmm_complex
 * 
 * @param z complex number
 * 
 * @return ctmm_complex conj(z)
 */
CTMM_EXPORT ctmm_complex ctmm_complex_conj(ctmm_complex z);

/**
 * @breif calculates inverse of ctmm_complex
 * 
 * @param z complex number
 * 
 * @return ctmm_complex 1/z
 */
CTMM_EXPORT ctmm_complex ctmm_complex_inv(ctmm_complex z);

/**
 * @breif calculates absolute value squared of ctmm_complex
 * 
 * @param z complex number
 * 
 * @return ctmm_complex |z|^2
 */
CTMM_EXPORT double ctmm_complex_abs2(ctmm_complex z);

/**
 * @breif calculates argument of ctmm_complex
 * 
 * @param z complex number
 * 
 * @return ctmm_complex atan2(Im(z), Re(z))
 */
CTMM_EXPORT double ctmm_complex_arg(ctmm_complex z);

/**
 * @breif sets matrix value mat[row][col] to val
 * 
 * @param mat pointer to matrix
 * @param row row
 * @param col column
 * @param val complex value of mat[row][col]
 */
CTMM_EXPORT void ctmm_matrix_set(ctmm_matrix *mat, unsigned int row,
    unsigned int col, ctmm_complex val);

/**
 * @breif gets matrix value mat[row][col]
 * 
 * @param mat pointer to matrix
 * @param row row
 * @param col column
 * 
 * @return ctmm_complex value of mat[row][col]
 */
CTMM_EXPORT ctmm_complex ctmm_matrix_get(ctmm_matrix *mat, unsigned int row,
    unsigned int col);

/**
 * @breif multiplies two ctmm matrices and sets result res = mat1*mat2
 * 
 * @param mat1 pointer to matrix1
 * @param mat1 pointer to matrix2
 * @param res pointer to result matrix
 */
CTMM_EXPORT void ctmm_matrix_mul(ctmm_matrix *mat1, ctmm_matrix *mat2,
    ctmm_matrix *res);

/**
 * @brief Pointer to stack parameters.
 */
CTMM_EXPORT typedef struct ctmm_stack *ctmm_stack;

/**
 * @brief Initialises a new stack.
 * 
 * The stack should be deallocated with ctmm_free_stack(stack).
 * 
 * @param nlyrs Number of layers in stack.
 * @param vwl Vacuum wavelength of illuminating light (metres).
 * @param t_in Angle of incidence (radians).
 * 
 * @return ctmm_stack with initialised arrays and matricies.
 */
CTMM_EXPORT ctmm_stack ctmm_create_stack(unsigned int nlyrs, double vwl,
    double t_in);

/**
 * @brief Frees memory allocated upon stack creation.
 * 
 * @param stack 
 */
CTMM_EXPORT void ctmm_free_stack(ctmm_stack stack);

/**
 * @brief Sets the angle of incidence.
 * 
 * @param stack ctmm_stack to set index in.
 * @param t_in Angle of incidence.
 */
CTMM_EXPORT void ctmm_set_t_in(ctmm_stack stack, double t_in);

/**
 * @brief Sets the refractive index of layer lyr_n to (x + y*i).
 * 
 * @param stack ctmm_stack to set index in.
 * @param lyr_n Layer number to set the index of.
 * @param x Real component of refractive index.
 * @param y Imaginary component of refractive index.
 */
CTMM_EXPORT void ctmm_set_ind(ctmm_stack stack, unsigned int lyr_n, double x,
    double y);

/**
 * @brief Sets the thickness of layer lyr_n to d.
 * 
 * @param stack ctmm_stack to set thickness in.
 * @param lyr_n Layer number to set the index of.
 * @param d Thickness of layer lyr_n.
 */
CTMM_EXPORT void ctmm_set_d(ctmm_stack stack, unsigned int lyr_n, double d);

/**
 * @brief Returns the refractive index of layer lyr_n to (x + y*i).
 * 
 * @param stack ctmm_stack to get index from.
 * @param lyr_n Layer number to get the index of.
 */
CTMM_EXPORT ctmm_complex ctmm_get_ind(ctmm_stack stack, unsigned int lyr_n);

/**
 * @brief Returns the thickness of layer lyr_n.
 * 
 * @param stack ctmm_stack to get thickness from.
 * @param lyr_n Layer number to get the index of.
 */
CTMM_EXPORT double ctmm_get_d(ctmm_stack stack, unsigned int lyr_n);

/**
 * @brief Returns the vacuum wavelength used by stack.
 * 
 * @param stack ctmm_stack to get vwl from.
 */
CTMM_EXPORT double ctmm_get_vwl(ctmm_stack stack);

/**
 * @brief Returns the angle of incidence used by stack.
 * 
 * @param stack ctmm_stack to get t_in from.
 */
CTMM_EXPORT double ctmm_get_t_in(ctmm_stack stack);

/**
 * @brief Returns the number of layers in stack.
 * 
 * @param stack ctmm_stack to get nlyrs from.
 */
CTMM_EXPORT unsigned int ctmm_get_nlyrs(ctmm_stack stack);

/**
 * @brief Returns a pointer to the first element of the stack matrix.
 * 
 * @param stack ctmm_stack to get matrix from.
 */
CTMM_EXPORT ctmm_matrix *ctmm_get_matrix(ctmm_stack stack);

/**
 * @brief Calculate kz.
 * 
 * Calculate the z (along stack) component of the wavevector from the
 * x (parallel to interfaces) component of the wavevector in the zeroth
 * layer.
 * 
 * @param n Refractive index of layer.
 * @param k Wavevector x component in zeroth layer.
 * @param vwl Vacuum wavelength of illuminating light.
 * 
 * @return Wavevector z component.
 */
CTMM_EXPORT ctmm_complex kz(ctmm_complex n, ctmm_complex k, double vwl);

/**
 * @brief Evaluate propigation matrix, for internal use.
 * 
 * @param pmat Matrix to be filled.
 * @param d Layer thickness.
 * @param n Layer refractive index.
 * @param k Wavevector component parallel to interfaces.
 * @param vwl Vacuum wavelength of illuminating light.
 */
CTMM_EXPORT void evaluate_prop(ctmm_matrix *pmat, double d, ctmm_complex n,
    ctmm_complex k, double vwl);

/**
 * @brief Calculate Fresnel coefficients.
 * 
 * @param coefs Array to be filled with coefficients.
 * @param n0 Refractive index in incident layer.
 * @param n1 Refractive index in exident layer.
 * @param k Wavevector component parallel to interfaces.
 * @param vwl Vacuum wavelength of illuminating light.
 */
CTMM_EXPORT void fresnel_coefs(ctmm_complex *coefs, ctmm_complex n0,
    ctmm_complex n1, ctmm_complex k, double vwl);

/**
 * @brief Calculate transfer matrix, for internal use.
 * 
 * @param tmat Matrix to be filled.
 * @param k Wavevector component parallel to interfaces.
 * @param n0 Refractive index in incident layer.
 * @param n1 Refractive index in exident layer.
 * @param vwl Vacuum wavelength of illuminating light.
 */
CTMM_EXPORT void evaluate_tran(ctmm_matrix *tmat, ctmm_complex k,
    ctmm_complex n0, ctmm_complex n1, double vwl);

/**
 * @brief Calculate the stack transfer matrix.
 * 
 * @param stack 
 */
CTMM_EXPORT void ctmm_evaluate(ctmm_stack stack);

/**
 * @brief Calculates power reflectivities and transmissions.
 * 
 * Reflectivity and transmission coefficients are calculated for both S and P
 * polarised light. The calculated coefficients are stored in the 'coefs' array.
 * The stack matrix must have already been evaluated.
 * 
 * This function may return incorrect results if the first or final layer in
 * the stack is absorbing (that is, have a non-zero imaginary component of the
 * refractive index).
 * 
 * @param coefs Array of coefficients ordered [rP, tP, rS, tS].
 * @param stack Stack struct.
 */
CTMM_EXPORT void ctmm_rts(ctmm_stack stack, double *coefs);

/**
 * @brief Calculates power and phase reflectivities and transmissions.
 * 
 * Reflectivity and transmission coefficients and phases are calculated for both
 * S and P polarised light. The calculated coefficients are stored in the
 * 'coefs' array. The stack matrix must have already been evaluated.
 * 
 * This function may return incorrect results if the first or final layer in
 * the stack is absorbing (that is, have a non-zero imaginary component of the
 * refractive index).
 * 
 * @param coefs Pointer to length 8 array of coefficients ordered [rP, tP, rS,
 *      tS, prP, ptP, prS, ptS].
 * @param stack Stack struct.
 */
CTMM_EXPORT void ctmm_rtps(ctmm_stack stack, double *coefs);

/**
 * @brief calculates fresnel coefficients of stack
 * 
 * @param coefs pointer to length 4 ctnn_complex array to hold coefficients
 *      ordered [r_p, t_p, r_s, t_s].
 * @param stack pointer to stack struct
 */
CTMM_EXPORT void ctmm_rtc(ctmm_stack stack, ctmm_complex *coefs);

#endif
