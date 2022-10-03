from glob import glob
from os import path
from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension, build_ext

home = path.realpath("..")

places = sorted(glob(path.join(home,"python/*.cpp"))+
                glob(path.join(home,"python/MultipoleKernel/*.cpp"))+
                glob(path.join(home,"python/NonThermalAsymptotic/*.cpp"))+
                glob(path.join(home,"StartUp/global_variables.cpp"))+
                glob(path.join(home,"Tools/Cosmology/*.cpp"))+
                glob(path.join(home,"Tools/Simple_routines/*.cpp"))+
                glob(path.join(home,"Tools/Integration/*.cpp"))+
                glob(path.join(home,"src/*.cpp"))+
                glob(path.join(home,"SZpack.cpp")))

dirs = [".", "..", "../include", "../src", "../src/database", "../src/database.opt", 
        "../StartUp", "../Tools/Definitions", "../Tools/Cosmology", "../Tools/Integration", 
        "../Tools/Simple_routines","/opt/local/include", "/usr/local/include"]

ext_modules = [
    Pybind11Extension(
        "SZpack", places,
        include_dirs=dirs,
        libraries=['gsl', 'gslcblas'], language="c++",
        library_dirs=["/opt/local/lib","/usr/local/lib"],
        define_macros=[('UsePybind11', None)]
    ),
]

setup(  name        = "SZpack",
        description = "This package provides the main SZpack functions for Python",
        author      = "E. Switzer, J. Chluba & E. Lee",
        version     = "2.0",
        cmdclass={"build_ext": build_ext},
        ext_modules=ext_modules,
        zip_safe=False
)