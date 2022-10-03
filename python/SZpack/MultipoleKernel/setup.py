from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext as build_pyx
from Cython.Build import cythonize
import numpy as np

dirs = [np.get_include(), ".","..", "../../..", "../../../include", "../../../src", "../../../src/database", "../../../src/database.opt", 
        "../../../StartUp", "../../../Tools/Definitions", "../../../Tools/Cosmology", "../../../Tools/Integration", 
        "../../../Tools/Simple_routines","/opt/local/include", "/usr/local/include"]

ext = [
    Extension("check", ["check.pyx","pyCallFunction.pyx"],
        include_dirs=dirs,
        libraries=['gsl', 'gslcblas'], language="c++",
        library_dirs=['.','..',"/opt/local/lib","/usr/local/lib"],
        extra_compile_args=['-O3','-fPIC','-std=c++11','-stdlib=libc++'],
        extra_link_args=['-O3']),
    Extension("params", ["params.pyx"],
        include_dirs=dirs,
        libraries=['gsl', 'gslcblas'], language="c++",
        library_dirs=['.','..',"/opt/local/lib","/usr/local/lib"],
        extra_compile_args=['-O3','-fPIC','-std=c++11','-stdlib=libc++'],
        extra_link_args=['-O3']),
    Extension("PySZ_electronDistribution", ["PySZ_electronDistribution.pyx"],
        include_dirs=dirs,
        libraries=['gsl', 'gslcblas'], language="c++",
        library_dirs=['.','..',"/opt/local/lib","/usr/local/lib"],
        extra_compile_args=['-O3','-fPIC','-std=c++11','-stdlib=libc++'],
        extra_link_args=['-O3']),
    Extension("PySZ_multipoleKernel", ["PySZ_multipoleKernel.pyx","pyCallFunction.pyx"],
        include_dirs=dirs,
        libraries=['gsl', 'gslcblas'], language="c++",
        library_dirs=['.','..',"/opt/local/lib","/usr/local/lib"],
        extra_compile_args=['-O3','-fPIC','-std=c++11','-stdlib=libc++'],
        extra_link_args=['-O3']),
]

setup(  name        = "PyMultipoleKernel",
        ext_modules = cythonize(ext,
                    compiler_directives={'language_level': 3})
        )