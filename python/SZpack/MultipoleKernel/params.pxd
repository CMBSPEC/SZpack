import numpy as np
cimport numpy as np
cimport cython
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool

#Import parameters class firstly
cdef extern from "parser.cpp":
    pass
cdef extern from "routines.cpp":
    pass
cdef extern from "Parameters.cpp":
    pass
cdef extern from "parser.h" :
    pass
cdef extern from "global_functions.h":
    pass
cdef extern from "routines.h" :
    pass
cdef extern from "Definitions.h" :
    pass
cdef extern from "physical_consts.h" :
    pass

cdef extern from "Parameters.h":
    cdef cppclass  MeansParameters : #Contains the Parameters needed for running means
        double Omega, Sigma, kappa    #For temperature dispersion
        vector[double] omegas, sigmas #For higher order temperature dispersion

        void assignOmegas(double omega0, double omega1, double omega2)
        void assignSigmas(double sigma0, double sigma1, double sigma2)

    cdef cppclass KernelParameters : #Contains the Parameters needed for running the kernel calculations
        double smin, smax
        vector[double] srange
        int l

    cdef cppclass DcomputeValues : #The parameters to calculate the derivatives for means
        int dThe, dbeta_para, dbeta2_perp
        vector[double] dDn_dThe

        void setValues(int DThe, int Dbeta_para, int Dbeta2_perp)

    cdef cppclass RareParameters :
        string RunMode #Takes values such as "all", "kin", "monopole", "dipole", "quadrupole", "monopole_corr" or "full".
        #If left blank or unrecognised, defaults to containing "all" (or "full", for 5D)

        void setTCMB(double T0)
        double TCMB()
        double Dn_DI_conversion() #conversion factor x^3 Dn --> DI in MJy / sr

    cdef cppclass CalculatedParameters :
        double The, gammao, gammac, mucc, xfac, xfacCMB, CNSN_Dtau
        double betac_para, betac2_perp

    cdef cppclass Parameters:
        double xmin, xmax, Dtau, Te, betac, muc, betao, muo
        int gridpoints
        vector[double] xcmb
        int T_order, beta_order  #only for Asymptotic & CNSN
        double relative_accuracy  #only for 5D & 3D
        int kmax, accuracy_level #only for CNSNopt

        double Te_max #only for the moment method running

        CalculatedParameters calc
        KernelParameters kernel
        MeansParameters means
        RareParameters rare

        DcomputeValues D

        #default values set in constructor
        Parameters() except +
        void copyParameters(Parameters copyParameters)
        void SetParametersFromFile(string filename)

        void Set_x(double x_min, double x_max, int gridpoints_i)
        void Set_x_from_nu(double nu_min, double nu_max, int gridpoints)
        void updateT(double T)
        void PopulateRunningArrays() #Call to update vectors if you change min/max or gridpoints
        void setCalcValues() #needs to be called if muc, muo, betac or betao are changed!

        bool CheckValues() #Check values are correct

cdef vec_to_py(vector[double] v)
cdef py_to_vec(np.ndarray v_py)
cdef vec_to_py_int(vector[int] v)

cdef class Params:
    cdef Parameters c_params