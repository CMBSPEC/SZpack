#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/functional.h>
#include "SZ_electron_distributions.h"

using namespace pybind11::literals;
namespace py = pybind11;

vector<double> c_default{0.755, 0.0609, 2.54, 13.3, 17.58};
vector<double> a_default{0.483, 152, 6843, 57.4, 12.4};

void init_ex_e_dists(py::module_ &m){
    m.def("boltzmann", &Boltzmann_Dist, "eta"_a, "Te"_a=5.0, 
            "The maxwell-juttner, 'relativistic boltzmann' distribution. This is the standard distribution used in SZ calculations. "
            "See the documentation for more details.");
    m.def("cosmic_ray", &CosmicRay_Dist, "eta"_a, "alpha"_a=2.5, "p1"_a=0.1, "p2"_a=10.0,
            "This is a power-law distribution, with -alpha as the exponent and p1 and p2 as the momentum cut offs, i.e., for p1 < eta < p2. "
            "Everywhere else the distribution returns 0. See the documentation or Ensslin & Kaiser (2000) for more details.");
    m.def("thermal_cosmic_ray", &ThermalCosmicRay_Dist, "eta"_a, "Te"_a=5.0, "alpha"_a=2.5, 
            "p1"_a = 0.2, "p2"_a=10.0, "This is a maxwell-juttner distribution for 0 < eta < p1, and then an power-law decay with exponent "
            "-alpha for p1 < eta < p2. This isn't normalised, please use the thermal_cosmic_ray_norm function to correct this.See the "
            "documentation or Ensslin & Kaiser (2000) for more details.");
    m.def("thermal_cosmic_ray_norm", &ThermalCosmicRay_Norm, "Te"_a=5.0, "alpha"_a=2.5, 
            "p1"_a = 0.2, "p2"_a=10.0, 
            "This is the normalising factor for the thermal cosmic ray distribution. This is a multiplicative factor which can only be "
            "calculated numerically. See the documentation for more details.");
    m.def("double_power_law", &DoublePower_Dist, "eta"_a, "alpha1"_a=0.5, "alpha2"_a=2.5,
            "p1"_a=0.01, "pcr"_a=0.2, "p2"_a=10.0, 
            "This is a double power-law i.e., two power-laws with a cross over at pcr. For p1 < eta < pcr the exponent is -alpha1 and "
            "for pcr < eta < p2 the exponent is alpha 2. See the documentation or Colfrancesco et al (2003) for more details.");
    m.def("kappa", &kappa_Dist, "eta"_a, "Te"_a=5.0, "kappa"_a = 2.0, 
            "The kappa distribution. This tends to a boltzmann distribution as kappa -> infinity. It has a power law tail with an exponent of "
            "-(kappa+1) for high momenta. See the documentation or Kaastra et al. (2009) for more details.");
    m.def("multimaxwellian", [](double eta, double Te, vector<double> c, vector<double> a) { return MultiMaxwellian_Dist(eta, Te, c, a);}, 
            "eta"_a, "Te"_a=5.0, "c"_a=c_default, "a"_a=a_default,
            "A distribution to mimic multiple overlaying thermal distributions. The arrays c and a represent the expansion coefficients of this "
            "method. See the documentation or Kaastra et al. (2009) for more details.");
}