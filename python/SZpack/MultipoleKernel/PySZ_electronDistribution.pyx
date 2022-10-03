cimport numpy as np
cimport cython

cdef extern from "SZ_electron_distributions.h" :
    cdef double Boltzmann_Dist(double eta, double Te)

    cdef double CosmicRay_Dist(double eta, double alpha, double p1, double p2)
    cdef double ThermalCosmicRay_Dist(double eta, double Te, double alpha, double p1, double p2)
    cdef double DoublePower_Dist(double eta, double alpha1, double alpha2, double p1, double pcr, double p2)
    cdef double kappa_Dist(double eta, double Te, double kappa)
    cdef double MultiMaxwellian_Dist(double eta, double Te)

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def boltzmann(double eta, double Te):
    return Boltzmann_Dist(eta, Te)

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def cosmic_ray(double eta, double alpha = 2.5, double p1 = 0.1, double p2 = 10.0):
    return CosmicRay_Dist(eta, alpha, p1, p2)

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def thermal_cosmic_ray(double eta, double Te, double alpha = 2.5, double p1 = 0.2, double p2 = 10.0):
    return ThermalCosmicRay_Dist(eta, Te, alpha, p1, p2)

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def double_power_law(double eta, double alpha1 = 0.5, double alpha2 = 2.5, 
                   double p1 = 0.01, double pcr = 0.2, double p2 = 10.0):
    return DoublePower_Dist(eta, alpha1, alpha2, p1, pcr, p2)

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def kappa(double eta, double Te, double kappa = 2.0):
    return kappa_Dist(eta, Te, kappa)

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def multi_maxwellian(double eta, double Te):
    return MultiMaxwellian_Dist(eta, Te)


##----------------------------------------------------------------------------------
## End User Functions
##----------------------------------------------------------------------------------
def boltz_p5(double eta):
    return boltzmann(eta, 0.5)

#TODO: Figure out a way to fix the integration routine so this still gives meaningful results even when the distribution is big or small
def boltz_2(double eta):
    return boltzmann(eta, 2.0)

def boltz_5(double eta):
    return boltzmann(eta, 5.0)

def boltz_10(double eta):
    return boltzmann(eta, 10.0)

def CR_default(double eta):
    return cosmic_ray(eta, 2.5, 0.1, 10.0)

def ThCR_T5_a2p5(double eta):
    return thermal_cosmic_ray(eta, 5.0, 2.5, 0.2, 10.0)

def ThCR_T5_a2(double eta):
    return thermal_cosmic_ray(eta, 5.0, 2.0, 0.2, 10.0)

def ThCR_T5_a3(double eta):
    return thermal_cosmic_ray(eta, 5.0, 3.0, 0.2, 10.0)

def ThCR_T5_a2p5_p0p3(double eta):
    return thermal_cosmic_ray(eta, 5.0, 2.5, 0.3, 10.0)

def kappa_T5_k2(double eta):
    return kappa(eta, 5.0, 2.0)

def kappa_T5_k5(double eta):
    return kappa(eta, 5.0, 5.0)

def kappa_T5_k10(double eta):
    return kappa(eta, 5.0, 10.0)

def multi_T5(double eta):
    return multi_maxwellian(eta, 5.0)