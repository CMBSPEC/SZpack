#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "Parameters.h"

using namespace pybind11::literals;
namespace py = pybind11;

void init_ex_params(py::module_ &m){
    py::class_<Parameters>(m, "parameters")
        .def(py::init<>())
        .def("copy_parameters", &Parameters::copyParameters, "A function to deep copy a parameters instance",
                "parameters_to_copy"_a)
        .def("set_parameters_from_file", &Parameters::SetParametersFromFile, "filename"_a,
                "A function to set the parameters from a given file, formatted as in runfiles/parameters.dat")
        .def("set_x_array", &Parameters::Set_x, "A function to set the xcmb array, x = (h*nu)/(kB*TCMB).", "x_min"_a, "x_max"_a, "gridpoints"_a)
        .def("set_x_from_nu", &Parameters::Set_x_from_nu, "A function to set the xcmb array from frequency values in GHz", 
                "nu_min"_a, "nu_max"_a, "gridpoints"_a)
        .def("check_values", &Parameters::CheckValues, "A function that checks all the values set for Parameters. "
                "Returns false if there are values that will not function elsewhere. This function will fix some rarer values to default values. "
                "Pay attention to the error messages.")
        .def_property_readonly("xmin", [](const Parameters &a) {return a.xmin;},
                "The minimum value of (observer frame) dimensionless frequency, x = (h*nu)/(kB*TCMB). xmin>=0.1.")
        .def_property_readonly("xmax", [](const Parameters &a) {return a.xmax;},
                "The maximum value of (observer frame) dimensionless frequency, x = (h*nu)/(kB*TCMB). xmin<=50.0.")
        .def_property_readonly("gridpoints", [](const Parameters &a) {return a.gridpoints;}, "The number of values in the x array.")
        .def_property_readonly("xcmb", [](const Parameters &a) {return py::array(a.xcmb.size(), a.xcmb.data());},
                "The array of values of (observer frame) dimensionless frequency, x = (h*nu)/(kB*TCMB).")
        .def_property_readonly("nucmb", [](Parameters &a) {vector<double> nucmb = a.get_nucmb(); return py::array(nucmb.size(), nucmb.data());},
                "The array of values of (observer frame) frequency [in GHz].")
        .def_property("Te", [](const Parameters &a) {return a.Te;}, &Parameters::updateT,
                "The cluster/line of sight temperature [in keV] in the cluster frame. Te > 0. Note we define y = Te*Dtau, so this also helps "
                "determine the amplitude of the signal.")
        .def_property("Dtau", [](const Parameters &a) {return a.Dtau;}, [](Parameters &a, double Dtau) {a.Dtau = Dtau; a.setCalcValues();},
                "The optical depth in the cluster frame. Dtau > 0.")
        .def_property("betac", [](const Parameters &a) {return a.betac;}, [](Parameters &a, double betac) {a.betac = betac; a.setCalcValues();},
                "The peculiar velocity of cluster with respect to CMB frame (betac = v/c).")
        .def_property("muc", [](const Parameters &a) {return a.muc;}, [](Parameters &a, double muc) {a.muc = muc; a.setCalcValues();},
                "The direction cosine of the cluster velocity with respect to the line-of-sight in the CMB rest frame.")
        .def_property("betao", [](const Parameters &a) {return a.betao;}, [](Parameters &a, double betao) {a.betao = betao; a.setCalcValues();},
                "The peculiar velocity of observer with respect to CMB frame (betao = v/c).")
        .def_property("muo", [](const Parameters &a) {return a.muo;}, [](Parameters &a, double muo) {a.muo = muo; a.setCalcValues();},
                "The direction cosine for the line-of-sight with respect to the observers velocity. The angle is measured in the observer frame.")
        .def_readwrite("T_order", &Parameters::T_order, "The order of expansion with respect to Te. Used only for CNSN and aymptotic methods.")
        .def_readwrite("beta_order", &Parameters::beta_order, 
                "The order of expansion with respect to betac. Used only for CNSN, CNSNopt and aymptotic methods.")
        .def_readwrite("relative_accuracy", &Parameters::relative_accuracy, 
                "The relative accuracy for numerical integration. Used only for the 3D and 5D methods.")
        .def_readwrite("kmax", &Parameters::kmax, "Used only for CNSNopt mode; depends on accuracy_level. See Table 1 of Chluba et al. 2013.")
        .def_readwrite("accuracy_level", &Parameters::accuracy_level, 
                "Used only for CNSNopt mode. See Table 1 of Chluba et al. 2013. accuracy_level = {0,1,2,3}.")
        .def_readwrite("Te_max", &Parameters::Te_max, "Used only in CNSNopt mode (for moment method). See Chluba et al. 2013.")
        .def_property("means_omega", [](const Parameters &a) {return a.means.Omega;}, [](Parameters &a, double Omega) {a.means.Omega = Omega;},
                "The second dimensionless temperature moment, <(Te-Te0)^2>/<Te>^2. i.e., the temperature dispersion term. Used only in means mode.")
        .def_property("means_sigma", [](const Parameters &a) {return a.means.Sigma;}, [](Parameters &a, double Sigma) {a.means.Sigma = Sigma;},
                "The first temperature-velocity moment, <(Te-Te0)(beta_para-beta_para0)>/<Te>. i.e., velocity-temperature cross term. "
                "Used only in means mode.")
        .def_property("means_kappa", [](const Parameters &a) {return a.means.kappa;}, [](Parameters &a, double kappa) {a.means.kappa = kappa;},
                "The second velocity moment, <(beta_para-beta_para0)^2>. i.e., the line-of-sight velocity dispersion term. "
                "Used only in means and means_ex modes.")
        .def_property_readonly("means_omegas", [](const Parameters &a) {return a.means.omegas;}, "Higher order dimensionless temperature moments. "
                "<(Te-Te0)^k>/<Te>^k for k = 2, 3 and 4 respectively. Used only in means_ex modes.")
        .def_property_readonly("means_sigmas", [](const Parameters &a) {return a.means.sigmas;}, "Higher order temperature-velocity moments. "
                "<(beta_para-beta_para0)*(Te-Te0)^k>/<Te>^k for k = 1, 2 and 3 respectively. Used only in means_ex modes.")
        .def("means_assign_omegas", [](Parameters &a, double omega0, double omega1, double omega2) {a.means.assignOmegas(omega0, omega1, omega2);}, 
                "omega0"_a, "omega1"_a, "omega2"_a, "A function to set the values for means_omegas.")
        .def("means_assign_sigmas", [](Parameters &a, double sigma0, double sigma1, double sigma2) {a.means.assignSigmas(sigma0, sigma1, sigma2);}, 
                "sigma0"_a, "sigma1"_a, "sigma2"_a, "A function to set the values for means_sigmas.")
        .def_property("runmode", [](const Parameters &a) {return a.rare.RunMode;}, [](Parameters &a, string runMode) {a.rare.RunMode = runMode;},
                "The runmode the distortion is calculated under. There are 7 modes recognised: \n"
                "\"monopole\": Only scattering of the monopole without the second order kinetic correction. i.e., O(betac^0 muc^0).\n"
                "\"dipole\": Only scattering of the dipole. That is, the first order kinematic correction. i.e., O(betac muc).\n"
                "\"quadrupole\": Only scattering of the quadrupole. That is, the second order kinematic correction. i.e., O(betac^2 muc^2).\n"
                "\"monopole_corr\": Only scattering of the second order kinematic correction to monopole. i.e., O(betac^2 muc^0).\n"
                "\"all\": All of the above terms added together.\n"
                "\"kin\": Only the kinematic terms. i.e., everything except \"monopole\" added together.\n"
                "\"full\": Only used for 5D method. Calculates the integral fully numerically with no approximations or analytic precomputation. "
                "This can be very slow. In other methods \"all\" is used when this is called.")
        .def_property("rare_TCMB", [](Parameters &a) {return a.rare.TCMB();}, [](Parameters &a, double TCMB) {a.rare.setTCMB(TCMB);},
                "The CMB temperature assumed in the calculations [in K]. By default this is set to T0_CMB = 2.726 K.")
        .def_property_readonly("rare_Dn_DI_conversion", [](Parameters &a) {return a.rare.Dn_DI_conversion();},
                "This is the value, conv, such that DI = conv*(x^3)*Dn, where DI is the distortion to the CMB in MJy/sr, and Dn is the change in "
                "photon density. In general, the DI parameter in functions means this won't be necessary. If returned true, the functions return DI "
                "and if false, Dn. This value depends on rare_TCMB. It is equal to 2*(kB*TCMB)^3/(h*c)^2.")
        .def_property_readonly("rare_x_nu_conversion", [](Parameters &a) {return a.rare.x_nu_conversion();},
                "This is the value, conv, such that nu = conv*x, where nu is the frequency in GHz and x the dimensionless frequency. This value "
                "depends on rare_TCMB. It is equal to TCMB/(h*kB)/1e9.")
        .def_property_readonly("calc_The", [](const Parameters &a) {return a.calc.The;}, 
                "The cluster/line of sight dimensionless temperature in the cluster frame. The = kB*Te/(me*c^2).")
        .def_property_readonly("calc_gammac", [](const Parameters &a) {return a.calc.gammac;},
                "The dimensionless energy associated with betac. i.e., gammac = 1/sqrt(1-betac^2).")
        .def_property_readonly("calc_gammao", [](const Parameters &a) {return a.calc.gammao;},
                "The dimensionless energy associated with betao. i.e., gammao = 1/sqrt(1-betao^2).")
        .def_property_readonly("calc_mucc", [](const Parameters &a) {return a.calc.mucc;},
                "muc in the cluster frame. Depends on betac and muc.")
        .def_property_readonly("calc_xfac", [](const Parameters &a) {return a.calc.xfac;},
                "The transformation factor to convert x into the cluster frame. Depends on betac, muc, betao and muo.")
        .def_property_readonly("calc_xfacCMB", [](const Parameters &a) {return a.calc.xfacCMB;},
                "The transformation factor to convert x into the CMB frame. Depends on betao and muo.")
        .def_property_readonly("calc_CNSN_Dtau", [](const Parameters &a) {return a.calc.CNSN_Dtau;},
                "Dtau but using the Nozawa 2006 convention instead of Chluba 2012. Depends on betac, muc and Dtau.")
        .def_property_readonly("calc_betac_para", [](const Parameters &a) {return a.calc.betac_para;},
                "betac in the direction of the line-of-sight. i.e., betac_para = betac*muc.")
        .def_property_readonly("calc_betac2_perp", [](const Parameters &a) {return a.calc.betac2_perp;},
                "The square of betac perpendicular to the line-of-sight. i.e., betac2_perp = (beta^2)*(1-muc^2).")
        .def("__repr__", [](const Parameters &a) { return "<SZpack.parameters>";})
        .def("__str__", [](const Parameters &a) {
                return "A summary of the most commonly used parameters:\n"
                   "xmin = "+to_string(a.xmin)+"; xmax = "+to_string(a.xmax)+"; gridpoints = "+to_string(a.gridpoints)+"\n"
                   "runmode = "+a.rare.RunMode+"; Dtau = "+to_string(a.Dtau)+"\n"
                   "Te = "+to_string(a.Te)+" keV (The ~ "+to_string(a.calc.The)+")\n"
                   "betac = "+to_string(a.betac)+"; muc = "+to_string(a.muc)+"\n"
                   "betao = "+to_string(a.betao)+"; muo = "+to_string(a.muo)+"\n"
                   "Runmode specific parameters:\n"
                   "relative_accuracy = "+to_string(a.relative_accuracy)+"\n"
                   "beta_order = "+to_string(a.beta_order)+"; T_order = "+to_string(a.T_order)+"\n"
                   "accuracy_level = "+to_string(a.accuracy_level)+"; kmax = "+to_string(a.kmax)+"; Te_max = "+to_string(a.Te_max);})
        ;
}