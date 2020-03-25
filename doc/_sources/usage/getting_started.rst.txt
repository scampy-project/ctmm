Getting started
===============

`ctmm` is intended to form a backend to other software, rather than be used
directly, in particular via the python bindings. The following brief example
demonstrates how `ctmm` can be used directly from C to model a single layer
of glass surrounded by air at normal incidence.

In order to use `ctmm` the header file `ctmm.h` must be included, and the
shared library produced during the compilation step linked against. ::

    #include "ctmm.h"

A `ctmm_stack` can then be created and initialised with ::

    ctmm_stack stack;
    stack = ctmm_create_stack(3, 633e-9, 0);

`ctmm_create_stack` takes three arguments, the number of layers (3), the
illuminating wavelength (`633e-9 nm` HeNe laser here) and the angle of incidence
(in radians).

Layer properties are set with the :ref:`ctmm_set_ind <ctmm_set_ind>` and
:ref:`ctmm_set_d <ctmm_set_d>` functions. All `ctmm` functions take a
`ctmm_stack` as the first argument. For both the set index and set thickness (d)
functions, the second argument is the layer number to set. For
:ref:`ctmm_set_ind <ctmm_set_ind>` this is follow by the real and imaginary
components of the refractive index, likewise for :ref:`ctmm_set_d <ctmm_set_d>`
this is followed by the layer thickness. For our air-glass-air example, the
properties will be set with ::

    ctmm_set_ind(stack, 0, 1, 0);
    ctmm_set_ind(stack, 1, 1.5, 0);
    ctmm_set_ind(stack, 2, 1, 0);

    ctmm_set_d(stack, 0, 0);
    ctmm_set_d(stack, 1, 0.25*(633e-9));
    ctmm_set_d(stack, 2, 0);

The stack transfer matrix can then be evaluated with ::

    ctmm_evaluate(stack);

If only the matrix is needed, a pointer to the `ctmm_matrix` transfer matrix
can be obtained with ::

    ctmm_matrix *tmat;
    tmat = ctmm_get_matrix(stack);

The `ctmm_matrix` type is a typedefed struct with a single member, a fixed
size array of `ctmm_complex` represeting a 4x4 matrix in row major order.

The `ctmm_complex` type is another typedefed struct that mimics the structure
of the C99 complex type, and exists purely to simplify the process of writing
`ctmm` for both C99 and non-C99 compliant compilers. On systems where the C99
complex type is available, the `ctmm_complex` type should be compatable.

Individual elements of a `ctmm_matrix` can be accessed with ::

    ctmm_complex elm;
    elm = ctmm_matrix_get(tmat, 0, 0);

In this case accessing the element `(0, 0)`.

`ctmm` can also calculate the amplitude and power reflectivity and transmission
coefficients for each polarisations, as well as the phase change on reflection
and transmission with the following :ref:`ctmm_rtc <ctmm_rtc>`,
:ref:`ctmm_rts <ctmm_rts>` and :ref:`ctmm_rtps <ctmm_rtps>` functions ::

    ctmm_complex rtc[4];
    double rts[4];
    double rtps[8];

    ctmm_rtc(stack, rtc);
    ctmm_rts(stack, rts);
    ctmm_rtps(stack, rtps);

where `rtc` is an array of the amplitude reflectivity and transmission
coefficients, `rts` is an array of the power reflectivity and transmission
coefficients, and `rtps` is an array of the power coefficients, and the phase
changes on reflection and transmission for each polarisation.

Once the stack is no longer needed, the allocated memory must be freed with ::

    ctmm_free_stack(stack);