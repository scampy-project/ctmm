# ctmm

ctmm is an optical thin film transfer matrix modelling library written in C. It is primarily designed to provide a lightweight and efficient backend for tspy, a python package for optical network modelling. A 4x4 transfer matrix methodology is implemented, treating both polarisations simultaneously.

For full documentation, see [ctmm.readthedocs.io](https://ctmm.readthedocs.io/)

## Getting Started

For useage of the library the test programs (see below) provide an example of applying the library to thin film optical modelling. A basic outline of the purpose and useage of each function may be found in ctmm.h.

A detailed treatment of the transfer matrix modelling technique employed by ctmm can be found in P. Yeh, “Optical Waves in Layered Media,” (John Wiley & Sons, Inc., Hoboken, New Jersey, 1998), pp. 62–63.

### Building from source

ctmm uses CMake for build file generation. Building should be possible on any system that provides the C standard libraries, however the process has only been tested on Windows (using NMAKE and MSVC build tools) and Linux (using GNU Make and gcc) systems. Once the source has been downloaded and CMake and a C compiler installed, the following commands build ctmm.

#### GNU Make

1. `mkdir build`
2. `cd build`
3. `cmake -B. .. -G"Unix Makefiles"`
4. `make all`

#### Microsoft NMAKE

From a MSVC Developer Command Prompt:
1. `mkdir build`
2. `cd build`
3. `cmake -B. .. -G"NMake Makefiles"`
4. `nmake all`

## Testing

Two basic test programs are provided, ctmm_test.c and ctmm_mathematics_test.c. ctmm_test.c is intended to test all functions related to transfer matrix modelling in the ctmm library. ctmm_mathematics_test.c tests the mathematical functions that have been redefined in ctmm to provide compatability with non-C99 standard compliant compilers. After building the project the compiled test binaries can be found under tests/bin. On systems using GNU Make these tests can be run with `make test` from the build folder. If ctmm has been build with NMAKE the tests can similarly be run with the command `nmake test`.

## Python bindings

A basic python interface to ctmm is included in the python subdirectory. pyctmm is reliant on a CPython 3 interpreter and numpy being installed. To install this python extension (after having built ctmm):

1. `cd python`
2. `pip install .`

## Versioning

[SemVer](http://semver.org/) is used for versioning. For the versions available, see the [tags on this repository](https://github.com/scampy-project/ctmm).

## Authors

* **Angus Bridges** - *Initial development* - [AngusBridges](https://github.com/AngusBridges)

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details
