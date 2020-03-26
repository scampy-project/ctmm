import setuptools
import os
import numpy as np

def main():
    library_dirs = ['./pyctmm', '.']
    include_dirs = ['./include', np.get_include()]

    if os.name == 'nt':
        libraries = ['ctmm']
        package_data = {'pyctmm': ['ctmm.dll']}

        modules = setuptools.Extension('pyctmm.pyctmm',
            sources = ['pyctmmmodule.c'],
            include_dirs = include_dirs,
            libraries = libraries,
            library_dirs = library_dirs)
    else:
        libraries = ['ctmm', 'm']
        package_data = {'pyctmm': ['libctmm.so']}

        modules = setuptools.Extension('pyctmm.pyctmm',
            sources = ['pyctmmmodule.c'],
            include_dirs = include_dirs,
            libraries = libraries,
            library_dirs = library_dirs,
            runtime_library_dirs=['$ORIGIN'])

    setuptools.setup(name = 'pyctmm',
        version = '1.0.1',
        description = 'Python interface for the ctmm optical transfer matrix modelling library',
        author = 'Angus Bridges',
        author_email = 'angus.bridges@npl.co.uk',
        ext_modules = [modules],
        packages = ['pyctmm'],
        package_data = package_data,
        include_package_data = True,
        install_requires=['numpy'],
        zip_safe = False
    )

if __name__ == '__main__':
    main()