import numpy as np
cimport numpy as np
cimport cython
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool

#Import parameters class firstly
cdef extern from "physical_consts.h" :
    pass
cdef extern from "Definitions.h" :
    pass
cdef extern from "routines.cpp":
    pass
cdef extern from "routines.h" :
    pass

cdef extern from "parser.cpp":
    pass
cdef extern from "parser.h" :
    pass
cdef extern from "global_functions.h":
    pass

cdef extern from "Parameters.cpp":
    pass
cdef extern from "Parameters.h" :
    pass

cdef extern from "Relativistic_MB.cpp" :
    pass
cdef extern from "Relativistic_MB.h" :
    pass
cdef extern from "nPl_derivatives.cpp" :
    pass
cdef extern from "nPl_derivatives.h" :
    pass

cdef extern from "Chebyshev_Int.cpp" :
    pass
cdef extern from "Chebyshev_Int.h" :
    pass
cdef extern from "Patterson.cpp" :
    pass
cdef extern from "Patterson.h" :
    pass
cdef extern from "Integration_routines.cpp" :
    pass
cdef extern from "Integration_routines.h" :
    pass
    

cdef extern from "SZ_multipoleKernel.cpp" :
    pass

from params cimport Parameters

#Definitions
cdef extern from "pyFunctionWrapper.h":
    cdef cppclass PyElectronDist:
        PyElectronDist()
        PyElectronDist(object)
        double run()

cdef extern from "SZ_multipoleKernel.h":
    cdef cppclass MultipoleKernel :
        MultipoleKernel() except +
        MultipoleKernel(int l_i, double s_i, double eta_i) except +
        MultipoleKernel(int l_i, double s_i, double eta_i, double Int_eps_i)  except +

        void Update_s(double s_i)
        double Calculate_integrated()
        double Calculate_formula()
        double Calculate_stable()
        double tMax, tMin

    cdef cppclass BeamKernel :
        BeamKernel() except +
        BeamKernel(int l_i, double s_i, double eta_i, double mup_i) except +

        void Update_s(double s_i)
        void Update_l(int l_i)
        double Calculate_formula()
        double tMax, tMin, mupMax, mupMin

    cdef cppclass IntegralKernel :
        IntegralKernel() except +
        IntegralKernel(double x_i, double betac_i, double muc_i, double eps_Int_i) except +
        IntegralKernel(int k, Parameters fp) except +

        void Update_x()
        double compute_distortion(string mode, PyElectronDist)
        double compute_kernel(int l_i, double s_i, PyElectronDist)
        double compute_distortion_fixed_eta(string mode, double eta_i)

        double compute_beam_distortion(string mode, PyElectronDist)
        double compute_beam_kernel(int l_i, double s_i, PyElectronDist)
        double compute_beam_distortion_fixed_eta(string mode, double eta_i)
