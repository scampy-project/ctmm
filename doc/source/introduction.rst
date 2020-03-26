.. |br| raw:: html

    <br />

Introduction
============

`ctmm` (C transfer matrix model) is a software library for modelling optical
multilayer thin films in C. It is designed to be as simple as possible, whilst
modelling both S and P polarisations simultaneously and correctly handling
absorbing (metallic) optical materials with complex refractive indices.

`ctmm` is primarily intended to be used as a backend for other software, rather
than as a standalone tool, in particular `ctmm` acts as a transfer matrix
modelling tool for the `strapy` optical modelling package. As such python
bindings are provided in the form of a CPython extension module.

For a detailed introduction to the transfer matrix modelling technique for
thin film modelling see:

    * P. Yeh, “Optical Waves in Layered Media,” (John Wiley & Sons, Inc.,Hoboken, New Jersey, 1998), pp. 62–63.

As per the terms of the MIT licence, this software is provided "as is", without warranty of any kind. Please see the `licence file <https://github.com/strapy-project/strapy/blob/master/LICENSE>`_ for full details.

Licencing
---------

`ctmm` is licenced under the MIT licence.

Reporting bugs
--------------

If you find a bug in `ctmm` please either create an issue on the GitHub
repository, or contact one of the authors (see below).

Authors
-------

    * `Angus Bridges <https://github.com/AngusBridges>`_:sup:`1, 2` (initial development) 
    * Andrew Yacoot :sup:`1`
    * Thomas Kissinger :sup:`2`
    * Ralph P. Tatam :sup:`2`

:sup:`1` National Physical Laboratory, Teddington, Middlesex, TW11 0LW, United Kingdom |br|
:sup:`2` Centre for Engineering Photonics, Cranfield University, MK43 0AL, United Kingdom

Bug reports and queries regarding the code should be directed to Angus Bridges
in the first instance.
