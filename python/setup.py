#===================================================================================================
#! /usr/bin/env python
#===================================================================================================

# System imports
from distutils.core import *
from distutils      import sysconfig

# Third-party modules - we depend on numpy for everything
import numpy

# Obtain the numpy include directory.  This logic works across numpy versions.
try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()

# sz_wrap extension module
_SZpack = Extension("_SZpack",
                   ["SZpack.i","SZpack.python.cpp"],
                   include_dirs = [numpy_include, "./.", "../.", "../include",
                                   "/usr/local/include","/opt/local/include"],
                   libraries = ['gsl', 'gslcblas', 'SZpack'],
                   library_dirs=['./.','../.',"/usr/local/lib","/opt/local/lib"],
                   extra_compile_args = ['-fPIC'],
                   swig_opts=['-modern', '-I../include', '-c++']
                   )

# NumyTypemapTests setup
setup(  name        = "SZpack functions",
        description = "This package provides the main SZpack functions for Python",

        author      = "E. Switzer & J. Chluba",
        version     = "1.0",
        ext_modules = [_SZpack]
        )
        
#===================================================================================================
#===================================================================================================
