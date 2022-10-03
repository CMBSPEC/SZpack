# distutils: language = c++

import numpy as np
from libcpp.vector cimport vector 
cimport cython

from params cimport Parameters
from params cimport Params
from params cimport vec_to_py
from params cimport py_to_vec
from PySZ_multipoleKernel cimport MultipoleKernel
from PySZ_multipoleKernel cimport BeamKernel
from PySZ_multipoleKernel cimport IntegralKernel

cdef extern from "SZ_multipoleKernel.h" :
    cdef void compute_SZ_distortion_kernel(vector[double] &Dn, vector[double] x, double betac, double muc, 
                                           double eps_Int, bool DI, PyElectronDist, string mode)
    cdef void compute_SZ_distortion_kernel(vector[double] &Dn, Parameters &fp, bool DI, PyElectronDist)

    cdef void compute_averaged_kernel(vector[double] &Dn, int l, vector[double] s, double betac, double muc, 
                                      double eps_Int, PyElectronDist)
    cdef void compute_averaged_kernel(vector[double] &Dn, Parameters &fp, PyElectronDist)

    cdef void compute_SZ_distortion_fixed_eta(vector[double] &Dn, vector[double] x, double eta, double betac, double muc,
                                              double eps_Int, bool DI, string mode)
    cdef void compute_SZ_distortion_fixed_eta(vector[double] &Dn, Parameters &fp, bool DI, double eta)

    cdef void compute_SZ_distortion_beam_kernel(vector[double] &Dn, vector[double] x, double mup, double betac, double muc, 
                                                double eps_Int, bool DI, PyElectronDist, string mode)
    cdef void compute_SZ_distortion_beam_kernel(vector[double] &Dn, Parameters &fp, bool DI, double mup, PyElectronDist)

    cdef void compute_averaged_beam_kernel(vector[double] &Dn, int l, vector[double] s, double mup, double betac, double muc, 
                                           double eps_Int, PyElectronDist)
    cdef void compute_averaged_beam_kernel(vector[double] &Dn, Parameters &fp, double mup, PyElectronDist)

    cdef void compute_SZ_distortion_beam_kernel_fixed_eta(vector[double] &Dn, vector[double] x, double mup, double eta, 
                                                          double betac, double muc, double eps_Int, bool DI, string mode)
    cdef void compute_SZ_distortion_beam_kernel_fixed_eta(vector[double] &Dn, Parameters &fp, bool DI, double mup, double eta)

# Multipole Kernel Implementation
def kernel_formula(int l, np.ndarray ss, double p0):
    results = np.zeros(len(ss))
    cdef MultipoleKernel c_MultipoleKernel = MultipoleKernel(l, ss[0], p0)
    results[0] = c_MultipoleKernel.Calculate_formula()
    for runvar in range(1, len(ss)):
        c_MultipoleKernel.Update_s(ss[runvar])
        results[runvar] = c_MultipoleKernel.Calculate_formula()
    return results

def kernel_integrated(int l, np.ndarray ss, double p0, double Int_eps = 1.0e-4):
    results = np.zeros(len(ss))
    cdef MultipoleKernel c_MultipoleKernel = MultipoleKernel(l, ss[0], p0, Int_eps)
    results[0] = c_MultipoleKernel.Calculate_integrated()
    for runvar in range(1, len(ss)):
        c_MultipoleKernel.Update_s(ss[runvar])
        results[runvar] = c_MultipoleKernel.Calculate_integrated()
    return results

def kernel_stable(int l, np.ndarray ss, double p0, double Int_eps = 1.0e-4):
    results = np.zeros(len(ss))
    cdef MultipoleKernel c_MultipoleKernel = MultipoleKernel(l, ss[0], p0, Int_eps)
    results[0] = c_MultipoleKernel.Calculate_stable()
    for runvar in range(1, len(ss)):
        c_MultipoleKernel.Update_s(ss[runvar])
        results[runvar] = c_MultipoleKernel.Calculate_stable()
    return results

def s_limits(double p0):
    gamma0 = np.sqrt(1+p0*p0)
    return np.log(np.array([(gamma0-p0)/(gamma0+p0),(gamma0+p0)/(gamma0-p0)]))


#Beam Kernel Implementation
def beam_kernel_formula(int l, np.ndarray ss, double p0, double mup):
    results = np.zeros(len(ss))
    cdef BeamKernel c_BeamKernel = BeamKernel(l, ss[0], p0, mup)
    results[0] = c_BeamKernel.Calculate_formula()
    for runvar in range(1, len(ss)):
        c_BeamKernel.Update_s(ss[runvar])
        results[runvar] = c_BeamKernel.Calculate_formula()
    return results

def beam_s_limits(double p0, double mup):
    cdef BeamKernel c_BeamKernel = BeamKernel(0, 0.0, p0, mup)
    return np.log(np.array([c_BeamKernel.tMin, c_BeamKernel.tMax]))

def beam_mup_limits(double p0, double s):
    cdef BeamKernel c_BeamKernel = BeamKernel(0, s, p0, 0.0)
    return np.array([c_BeamKernel.mupMin, c_BeamKernel.mupMax])

#The Integrals of the Kernel Implementation
def distortion(Params functionParameters, electronDist, bool DI=True):
    cdef Parameters fp = functionParameters.c_params
    cdef PyElectronDist func = PyElectronDist(electronDist)
    cdef vector[double] Dn
    compute_SZ_distortion_kernel(Dn, fp, DI, func)
    return vec_to_py(Dn)

def distortion_from_variables(np.ndarray xcmb, double betac, double muc, electronDist, str mode, double eps_Int, bool DI=True):
    b_mode = mode.encode('utf-8')
    cdef PyElectronDist func = PyElectronDist(electronDist)
    cdef vector[double] Dn
    compute_SZ_distortion_kernel(Dn, py_to_vec(xcmb), betac, muc, eps_Int, DI, func, b_mode)
    return vec_to_py(Dn)

def averaged_kernel(Params functionParameters, electronDist, int l = -1):
    cdef Parameters fp = functionParameters.c_params
    if (l>=0): fp.kernel.l = l
    cdef PyElectronDist func = PyElectronDist(electronDist)
    cdef vector[double] Dn
    compute_averaged_kernel(Dn, fp, func)
    return vec_to_py(Dn)

def averaged_kernel_from_variables(np.ndarray s_array, int l, double betac, double muc, electronDist, double eps_Int):
    cdef PyElectronDist func = PyElectronDist(electronDist)
    cdef vector[double] Dn
    compute_averaged_kernel(Dn, l, py_to_vec(s_array), betac, muc, eps_Int, func)
    return vec_to_py(Dn)

def distortion_fixed_eta(Params functionParameters, double eta, bool DI=True):
    cdef Parameters fp = functionParameters.c_params
    cdef vector[double] Dn
    compute_SZ_distortion_fixed_eta(Dn, fp, DI, eta)
    return vec_to_py(Dn)

def distortion_fixed_eta_from_variables(np.ndarray xcmb, double eta, double betac, double muc, str mode, double eps_Int, bool DI=True):
    b_mode = mode.encode('utf-8')
    cdef vector[double] Dn
    compute_SZ_distortion_fixed_eta(Dn, py_to_vec(xcmb), eta, betac, muc, eps_Int, DI, b_mode)
    return vec_to_py(Dn)

def distortion_beam(Params functionParameters, bool DI, electronDist, double mup):
    cdef Parameters fp = functionParameters.c_params
    cdef PyElectronDist func = PyElectronDist(electronDist)
    cdef vector[double] Dn
    compute_SZ_distortion_beam_kernel(Dn, fp, DI, mup, func)
    return vec_to_py(Dn)

def distortion_beam_from_variables(np.ndarray xcmb, double mup, double betac, double muc, electronDist, str mode, double eps_Int, bool DI=True):
    b_mode = mode.encode('utf-8')
    cdef PyElectronDist func = PyElectronDist(electronDist)
    cdef vector[double] Dn
    compute_SZ_distortion_beam_kernel(Dn, py_to_vec(xcmb), mup, betac, muc, eps_Int, DI, func, b_mode)
    return vec_to_py(Dn)

def averaged_kernel_beam(Params functionParameters, electronDist, double mup, int l = -1):
    cdef Parameters fp = functionParameters.c_params
    cdef PyElectronDist func = PyElectronDist(electronDist)
    cdef vector[double] Dn
    compute_averaged_beam_kernel(Dn, fp, mup, func)
    return vec_to_py(Dn)

def averaged_kernel_beam_from_variables(np.ndarray s_array, int l, double mup, double betac, double muc, 
                                        electronDist, double eps_Int):
    cdef PyElectronDist func = PyElectronDist(electronDist)
    cdef vector[double] Dn
    compute_averaged_beam_kernel(Dn, l, py_to_vec(s_array), mup, betac, muc, eps_Int, func)
    return vec_to_py(Dn)

def distortion_beam_fixed_eta(Params functionParameters, double mup, double eta, bool DI=True):
    cdef Parameters fp = functionParameters.c_params
    cdef vector[double] Dn
    compute_SZ_distortion_beam_kernel_fixed_eta(Dn, fp, DI, mup, eta)
    return vec_to_py(Dn)

def distortion_beam_fixed_eta_from_variables(np.ndarray xcmb, double mup, double eta, double betac, double muc, 
                                             str mode, double eps_Int, bool DI=True):
    b_mode = mode.encode('utf-8')
    cdef vector[double] Dn
    compute_SZ_distortion_beam_kernel_fixed_eta(Dn, py_to_vec(xcmb), mup, eta, betac, muc, eps_Int, DI, b_mode)
    return vec_to_py(Dn)