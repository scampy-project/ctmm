API documentation
=================

`ctmm` consists of a single shared library, and an associated header file. The
functions in this header file can be broken down into three main groups, only
one of which should need to be used by the end user.

The core functions group documents all functions needed to model thin film
stacks and calculate their properties.

Internal functions are documented in the internal functions group, these are
used to calculate optical coefficients and intermediate matrices that are
required for the evaluation of the stack matrix. These are exposed in the header
mainly to allow testing, although calculating individual propagation and
transfer matrices may be useful.

Complex mathematical functions redefined in `ctmm` purely to allow use with
non-C99 compliant compilers asre documented in the mathematical functions group.
Again, these are exposed largely for testing purposes.

.. toctree::
   :maxdepth: 1
   :caption: Function groups:

   ctmm
   ctmm_internal
   ctmm_mathematical