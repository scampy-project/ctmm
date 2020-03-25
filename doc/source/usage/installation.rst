Installation
============

ctmm uses CMake for build file generation, generating makefiles for either GNU
Make or Microsoft NMAKE. Building should be possible on any system that provides
the C standard libraries, however the process has only been
tested on Windows (using NMAKE and MSVC build tools) and Linux (using GNU Make
and gcc) systems. Once the source has been downloaded and CMake and a C compiler
installed, the following commands build ctmm.

GNU Make
--------
From a terminal open in the ctmm root directory: ::

    mkdir build
    cd build
    cmake -B. .. -G"Unix Makefiles"
    make all

Microsoft NMAKE
---------------

Before building `ctmm` `Build Tools for Visual Studio
<https://visualstudio.microsoft.com/downloads>`_ (scroll down to Tools for
Visual Studio, Build Tools for Visual Studio) must be installed to provide a
python compatable C compiler, and `CMake <https://cmake.org/download/>`_
(select 'Add CMake to the system path for current user' during the installation
process) must be installed to generate makefiles. Once the build tools have been
installed, from an MSVC Developer Command Prompt open in the ctmm root
directory: ::

    mkdir build
    cd build
    cmake -B. .. -G"NMake Makefiles"
    nmake all

If the error
"`'cmake' is not recognized as an internal or external command, operable program or batch file.`"
is returned, ensure CMake is installed, and has been added to the system PATH
environment variable.

If CMake cannot find the compiler
("`cl is not a full path and was not found in the PATH`" or similar error) check
that an Developer Command Prompt has been used. A suitable command prompt can
be started from the Windows start menu under
`Developer Command Prompt for VS 20XX`.

Testing
-------

Two basic test programs are provided, `ctmm_test.c` and
`ctmm_mathematics_test.c`. `ctmm_test.c` is intended to test all functions
related to transfer matrix modelling in the `ctmm` library.
`ctmm_mathematics_test.c` tests the mathematical functions that have been
redefined in ctmm to provide compatability with non-C99 standard compliant
compilers. After building the project the compiled test binaries can be found
under `tests/bin`. On systems using GNU Make these tests can be run with ::

    make test

from the build folder. If ctmm has been build with NMAKE the tests can similarly
be run with the command ::

    nmake test

Python bindings
---------------

To install the python bindings (`pyctmm`) following the build step: ::

    cd ../python
    pip install .

To use `pyctmm` a python 3 installation (see installation instructions for
`Windows <https://docs.python.org/3/using/windows.html>`_, or
`Unix <https://docs.python.org/3/using/unix.html>`_) is required.
