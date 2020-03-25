Python bindings
===============

`ctmm` can be accessed from python through the `pyctmm` module, a CPython
extension module exposing the core functions of `ctmm`.

Installation
------------

The python bindings can be installed, after having compiled `ctmm`, with ::

    cd python
    python setup.py install --user

A C compiler must be available on the system to compile the CPython extension,
on Microsoft Windows this means the installation commands must be run from an
MSVC Developer Command Prompt.

Usage
-----

The `pyctmm` module can be imported into python with ::

    import pyctmm

Core functions are accessed from this module, for example a three layer stack
illuminated with light at `633 nm` at normal incidence is created with ::

    pyctmm.create_stack(3, 633e-9, 0);