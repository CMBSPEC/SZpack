# distutils: language = c++

import numpy as np
from libcpp.vector cimport vector 

from params cimport Parameters

cdef vec_to_py(vector[double] v):
    cdef int length = v.size()
    cdef np.ndarray v_py = np.zeros(length)
    for i in range(length):
        v_py[i] = v[i]
    return v_py

cdef py_to_vec(np.ndarray v_py):
    cdef int length = len(v_py)
    cdef vector[double] v
    v.resize(length)
    for i in range(length):
        v[i] = v_py[i]
    return v

cdef vec_to_py_int(vector[int] v):
    cdef int length = v.size()
    cdef np.ndarray v_py = np.zeros(length)
    for i in range(length):
        v_py[i] = v[i]
    return v_py

cdef class Params:
    def __cinit__(self):
        self.c_params = Parameters()

    def copy_parameters(self, Params copyParameters) :
        self.c_params.copyParameters(copyParameters.c_params)

    def set_parameters_from_file(self, str filename) :
        b_filename= filename.encode('utf-8')
        self.c_params.SetParametersFromFile(b_filename)

    def set_x_array(self, float x_min, float x_max, int gridpoints) :
        self.c_params.Set_x(x_min, x_max, gridpoints)

    def set_s_array(self, float s_min, float s_max, int gridpoints) :
        self.c_params.kernel.smin = s_min
        self.c_params.kernel.smax = s_max
        self.c_params.gridpoints = gridpoints
        self.c_params.PopulateRunningArrays()

    def set_x_from_nu(self, float nu_min, float nu_max, int gridpoints) :
        self.c_params.Set_x_from_nu(nu_min, nu_max, gridpoints)

    def check_values(self) :
        return self.c_params.CheckValues()

    def means_assign_omegas(self, float omega0, float omega1, float omega2) :
        self.c_params.means.assignOmegas(omega0, omega1, omega2)

    def means_assign_sigmas(self, float sigma0, float sigma1, float sigma2) :
        self.c_params.means.assignSigmas(sigma0, sigma1, sigma2)

    # Attribute access
    @property
    def xmin(self):
        return self.c_params.xmin

    @property
    def xmax(self):
        return self.c_params.xmax

    @property
    def gridpoints(self):
        return self.c_params.gridpoints

    @property
    def xcmb(self):
        return vec_to_py(self.c_params.xcmb)

    @property
    def nucmb(self):
        const_h_kb_e9 = 1.0e9*6.62606896e-27/(1.3806504e-23*1.0e+7)
        return vec_to_py(self.c_params.xcmb)*self.c_params.rare.TCMB()/const_h_kb_e9

    @property
    def Te(self):
        return self.c_params.Te
    @Te.setter
    def Te(self, Te):
        self.c_params.updateT(Te)

    @property
    def Dtau(self):
        return self.c_params.Dtau
    @Dtau.setter
    def Dtau(self, Dtau):
        self.c_params.Dtau = Dtau
        self.c_params.setCalcValues()

    @property
    def betac(self):
        return self.c_params.betac
    @betac.setter
    def betac(self, betac):
        self.c_params.betac = betac
        self.c_params.setCalcValues()

    @property
    def muc(self):
        return self.c_params.muc
    @muc.setter
    def muc(self, muc):
        self.c_params.muc = muc
        self.c_params.setCalcValues()

    @property
    def betao(self):
        return self.c_params.betao
    @betao.setter
    def betao(self, betao):
        self.c_params.betao = betao
        self.c_params.setCalcValues()

    @property
    def muo(self):
        return self.c_params.muo
    @muo.setter
    def muo(self, muo):
        self.c_params.muo = muo
        self.c_params.setCalcValues()

    #Less common attributes
    @property
    def T_order(self):
        return self.c_params.T_order
    @T_order.setter
    def T_order(self, T_order):
        self.c_params.T_order = T_order
    
    @property
    def beta_order(self):
        return self.c_params.beta_order
    @beta_order.setter
    def beta_order(self, beta_order):
        self.c_params.beta_order = beta_order

    @property
    def relative_accuracy(self):
        return self.c_params.relative_accuracy
    @relative_accuracy.setter
    def relative_accuracy(self, relative_accuracy):
        self.c_params.relative_accuracy = relative_accuracy

    @property
    def kmax(self):
        return self.c_params.kmax
    @kmax.setter
    def kmax(self, kmax):
        self.c_params.kmax = kmax
    
    @property
    def accuracy_level(self):
        return self.c_params.accuracy_level
    @accuracy_level.setter
    def accuracy_level(self, accuracy_level):
        self.c_params.accuracy_level = accuracy_level
    
    @property
    def Te_max(self):
        return self.c_params.Te_max
    @Te_max.setter
    def Te_max(self, Te_max):
        self.c_params.Te_max = Te_max

    @property
    def means_omega(self):
        return self.c_params.means.Omega
    @means_omega.setter
    def means_omega(self, Omega):
        self.c_params.means.Omega = Omega

    @property
    def means_sigma(self):
        return self.c_params.means.Sigma
    @means_sigma.setter
    def means_sigma(self, Sigma):
        self.c_params.means.Sigma = Sigma

    @property
    def means_kappa(self):
        return self.c_params.means.kappa
    @means_kappa.setter
    def means_kappa(self, kappa):
        self.c_params.means.kappa = kappa

    @property
    def means_omegas(self):
        return vec_to_py(self.c_params.means.omegas)

    @property
    def means_sigmas(self):
        return vec_to_py(self.c_params.means.sigmas)

    @property
    def kernel_l(self):
        return self.c_params.kernel.l
    @kernel_l.setter
    def kernel_l(self, l):
        self.c_params.kernel.l = l

    @property
    def kernel_smin(self):
        return self.c_params.kernel.smin

    @property
    def kernel_smax(self):
        return self.c_params.kernel.smax

    @property
    def kernel_srange(self):
        return vec_to_py(self.c_params.kernel.srange)
    
    @property
    def rare_RunMode(self):
        return self.c_params.rare.RunMode
    @rare_RunMode.setter
    def rare_RunMode(self, RunMode):
        b_RunMode = RunMode.encode('utf-8')
        self.c_params.rare.RunMode = b_RunMode

    @property
    def rare_TCMB(self):
        return self.c_params.rare.TCMB()
    @rare_TCMB.setter
    def rare_TCMB(self, rare_TCMB):
        self.c_params.rare.setTCMB(rare_TCMB)

    @property
    def rare_Dn_DI_conversion(self):
        return self.c_params.rare.Dn_DI_conversion()

    @property
    def calc_The(self):
        return self.c_params.calc.The

    @property
    def calc_gammao(self):
        return self.c_params.calc.gammao

    @property
    def calc_gammac(self):
        return self.c_params.calc.gammac

    @property
    def calc_mucc(self):
        return self.c_params.calc.mucc

    @property
    def calc_xfac(self):
        return self.c_params.calc.xfac

    @property
    def calc_xfacCMB(self):
        return self.c_params.calc.xfacCMB

    @property
    def calc_CNSN_Dtau(self):
        return self.c_params.calc.CNSN_Dtau

    @property
    def calc_betac_para(self):
        return self.c_params.calc.betac_para

    @property
    def calc_betac2_perp(self):
        return self.c_params.calc.betac2_perp