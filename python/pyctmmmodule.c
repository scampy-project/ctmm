#define PY_SSIZE_T_CLEAN
#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "numpy/arrayobject.h"
#include "ctmm.h"

void pyctmm_free_stack(PyObject *capsule)
{
    ctmm_stack stack_ptr = (ctmm_stack) PyCapsule_GetPointer(capsule, NULL);
    ctmm_free_stack(stack_ptr);
}

static PyObject *method_create_stack(PyObject *self,  PyObject *args)
{
    unsigned int nlyrs;
    double vwl;
    double t_in;
    ctmm_stack stack;
    PyObject *stack_ptr;

    if (!PyArg_ParseTuple(args, "Idd", &nlyrs, &vwl, &t_in)) {
        return NULL;
    }

    stack = ctmm_create_stack(nlyrs, vwl, t_in);
    stack_ptr = PyCapsule_New(stack, NULL, pyctmm_free_stack);

    return stack_ptr;
}

static PyObject *method_set_t_in(PyObject *self,  PyObject *args)
{
    PyObject *stack_ptr;
    ctmm_stack stack;
    double t_in;

    if (!PyArg_ParseTuple(args, "Od", &stack_ptr, &t_in)) {
        return NULL;
    }

    stack = (ctmm_stack) PyCapsule_GetPointer(stack_ptr, NULL);
    ctmm_set_t_in(stack, t_in);

    Py_RETURN_NONE;
}

static PyObject *method_set_ind(PyObject *self,  PyObject *args)
{
    PyObject *stack_ptr;
    ctmm_stack stack;
    unsigned int lyr_n;
    double re;
    double im;

    if (!PyArg_ParseTuple(args, "OIdd", &stack_ptr, &lyr_n, &re, &im)) {
        return NULL;
    }

    stack = (ctmm_stack) PyCapsule_GetPointer(stack_ptr, NULL);
    ctmm_set_ind(stack, lyr_n, re, im);

    Py_RETURN_NONE;
}

static PyObject *method_set_d(PyObject *self,  PyObject *args)
{
    PyObject *stack_ptr;
    ctmm_stack stack;
    unsigned int lyr_n;
    double d;

    if (!PyArg_ParseTuple(args, "OId", &stack_ptr, &lyr_n, &d)) {
        return NULL;
    }

    stack = (ctmm_stack) PyCapsule_GetPointer(stack_ptr, NULL);
    ctmm_set_d(stack, lyr_n, d);

    Py_RETURN_NONE;
}

static PyObject *method_get_ind(PyObject *self,  PyObject *args)
{
    PyObject *stack_ptr;
    ctmm_stack stack;
    unsigned int lyr_n;
    ctmm_complex ind;

    if (!PyArg_ParseTuple(args, "OI", &stack_ptr, &lyr_n)) {
        return NULL;
    }

    stack = (ctmm_stack) PyCapsule_GetPointer(stack_ptr, NULL);
    ind = ctmm_get_ind(stack, lyr_n);

    return Py_BuildValue("dd", creal(ind), cimag(ind));
}

static PyObject *method_get_d(PyObject *self,  PyObject *args)
{
    PyObject *stack_ptr;
    ctmm_stack stack;
    unsigned int lyr_n;
    double d;

    if (!PyArg_ParseTuple(args, "OI", &stack_ptr, &lyr_n)) {
        return NULL;
    }

    stack = (ctmm_stack) PyCapsule_GetPointer(stack_ptr, NULL);
    d = ctmm_get_d(stack, lyr_n);

    return Py_BuildValue("d", d);
}

static PyObject *method_evaluate(PyObject *self,  PyObject *args)
{
    PyObject *stack_ptr;
    ctmm_stack stack;

    if (!PyArg_ParseTuple(args, "O", &stack_ptr)) {
        return NULL;
    }

    stack = (ctmm_stack) PyCapsule_GetPointer(stack_ptr, NULL);
    ctmm_evaluate(stack);

    Py_RETURN_NONE;
}

static PyObject *method_get_matrix(PyObject *self,  PyObject *args)
{
    PyObject *stack_ptr;
    ctmm_stack stack;
    PyObject *stack_ndarray;
    npy_intp dims[2];

    if (!PyArg_ParseTuple(args, "O", &stack_ptr)) {
        return NULL;
    }

    stack = (ctmm_stack) PyCapsule_GetPointer(stack_ptr, NULL);

    dims[0] = 4;
    dims[1] = 4;

    stack_ndarray = PyArray_SimpleNewFromData(2, dims, NPY_COMPLEX128,
        &(ctmm_get_matrix(stack)->data));
    Py_INCREF(stack_ptr);

    return Py_BuildValue("O", stack_ndarray);
}

static PyObject *method_get_amplitude(PyObject *self,  PyObject *args)
{
    PyObject *stack_ptr;
    ctmm_stack stack;
    PyObject *rtc_ndarray;
    npy_intp dims[1];

    if (!PyArg_ParseTuple(args, "O", &stack_ptr)) {
        return NULL;
    }

    stack = (ctmm_stack) PyCapsule_GetPointer(stack_ptr, NULL);

    dims[0] = 4;
    rtc_ndarray = PyArray_SimpleNew(1, dims, NPY_COMPLEX128);
    ctmm_rtc(stack,
        (ctmm_complex *) PyArray_BYTES((PyArrayObject *) rtc_ndarray));

    return Py_BuildValue("O", rtc_ndarray);
}

static PyObject *method_get_power(PyObject *self,  PyObject *args)
{
    PyObject *stack_ptr;
    ctmm_stack stack;
    PyObject *rts_ndarray;
    npy_intp dims[1];

    if (!PyArg_ParseTuple(args, "O", &stack_ptr)) {
        return NULL;
    }

    stack = (ctmm_stack) PyCapsule_GetPointer(stack_ptr, NULL);

    dims[0] = 4;
    rts_ndarray = PyArray_SimpleNew(1, dims, NPY_FLOAT64);
    ctmm_rts(stack, (double *) PyArray_BYTES((PyArrayObject *) rts_ndarray));

    return Py_BuildValue("O", rts_ndarray);
}

static PyObject *method_get_power_phase(PyObject *self,  PyObject *args)
{
    PyObject *stack_ptr;
    ctmm_stack stack;
    PyObject *rtps_ndarray;
    npy_intp dims[1];

    if (!PyArg_ParseTuple(args, "O", &stack_ptr)) {
        return NULL;
    }

    stack = (ctmm_stack) PyCapsule_GetPointer(stack_ptr, NULL);

    dims[0] = 8;
    rtps_ndarray = PyArray_SimpleNew(1, dims, NPY_FLOAT64);
    ctmm_rtps(stack, (double *) PyArray_BYTES((PyArrayObject *) rtps_ndarray));

    return Py_BuildValue("O", rtps_ndarray);
}

static PyObject *method_get_fresnel(PyObject *self,  PyObject *args)
{
    PyObject *fresnel_ndarray;
    npy_intp dims[1];

    double n0_re;
    double n0_im;
    double n1_re;
    double n1_im;
    double k_re;
    double k_im;
    double vwl;

    if (!PyArg_ParseTuple(args, "ddddddd",
        &n0_re,
        &n0_im,
        &n1_re,
        &n1_im,
        &k_re,
        &k_im,
        &vwl)) {
        return NULL;
    }

    dims[0] = 8;
    fresnel_ndarray = PyArray_SimpleNew(1, dims, NPY_COMPLEX128);
    fresnel_coefs((ctmm_complex *) PyArray_BYTES((PyArrayObject *)
        fresnel_ndarray),
        ctmm_complex_set(n0_re, n0_im),
        ctmm_complex_set(n1_re, n1_im),
        ctmm_complex_set(k_re, k_im),
        vwl);

    return Py_BuildValue("O", fresnel_ndarray);
}

static PyMethodDef PyctmmMethods[] = {
    {"create_stack", method_create_stack, METH_VARARGS,
        "Initialised pyctmm stack."},
    {"set_t_in", method_set_t_in, METH_VARARGS,
        "Set angle of incidence."},
    {"set_ind", method_set_ind, METH_VARARGS,
        "Set layer index."},
    {"set_d", method_set_d, METH_VARARGS,
        "Set layer thickness."},
    {"get_ind", method_get_ind, METH_VARARGS,
        "Get layer index."},
    {"get_d", method_get_d, METH_VARARGS,
        "Get layer thickness."},
    {"evaluate", method_evaluate, METH_VARARGS,
        "Evaluate transfer matrix."},
    {"get_matrix", method_get_matrix, METH_VARARGS,
        "Return numpy array of type NP_COMPLEX128 containing the stack matrix."},
    {"get_amplitude", method_get_amplitude, METH_VARARGS,
        "Return numpy array of type NP_COMPLEX128 containing the amplitude reflectivity and transmission coefficients."},
    {"get_power", method_get_power, METH_VARARGS,
        "Return numpy array of type NP_FLOAT64 containing the power reflectivity and transmission coefficients."},
    {"get_power_phase", method_get_power_phase, METH_VARARGS,
        "Return numpy array of type NP_FLOAT64 containing the power and phase reflectivity and transmission coefficients."},
    {"get_fresnel", method_get_fresnel, METH_VARARGS,
        "Return numpy array of type NP_COMPLEX128 containing calculated Fresnel coefficients."},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef pyctmmmodule = {
    PyModuleDef_HEAD_INIT,
    "pyctmm",
    "Python interface for ctmm transfer matrix modelling library",
    -1,
    PyctmmMethods
};

PyMODINIT_FUNC PyInit_pyctmm(void)
{
    PyObject *m;

    m = PyModule_Create(&pyctmmmodule);
    if (m == NULL) {
        return NULL;
    }

    import_array();

    return m;
}