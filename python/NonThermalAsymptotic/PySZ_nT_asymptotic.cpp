#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/functional.h>
#include "SZ_nT_asymptotic.h"
#include "routines.h"
#include "Parameters.h"

using namespace pybind11::literals;
namespace py = pybind11;

void init_ex_nT_asym(py::module_ &m){
    m.def("radio_distortion", [](Parameters fp, bool DI, bool CMBframe){
                    vector<double> Dn; compute_radio_distortion(Dn, fp, DI, CMBframe);
                    return py::array(Dn.size(), Dn.data());}, "Params"_a, "DI"_a=true, "CMBframe"_a=false,
                    "A function to calculate the SZ signal from an asymptotic expansion around Te -- for a radio signal as described in "
                    "Lee & Chluba 2021. See the documentation for more details.");

    m.def("distortion", [](const std::function<double(int, double)> &dDist, Parameters fp, bool DI, bool CMBframe){
                    vector<double> Dn; compute_nonthermal_distortion_asymptotic(Dn, fp, dDist, DI, CMBframe);
                    return py::array(Dn.size(), Dn.data());}, "dist_derivs"_a, "Params"_a, "DI"_a=true, "CMBframe"_a=false,
                    "A function to calculate the SZ signal from an asymptotic expansion around Te, for a given set of background photon "
                    "distribution derivatives i.e. x^k x^k/dx^k(background light). See the documentation or Lee & Chluba 2021 for more details.");

    m.def("distortion_from_variables", [](const std::function<double(int, double)> &dDist, vector<double> xcmb, double The, double betac, 
                    double muc, int Te_order, int betac_order, string mode, bool DI, bool CMBframe){
                    vector<double> Dn; compute_nonthermal_distortion_asymptotic(Dn, xcmb, The, betac, muc, Te_order, betac_order, dDist, 
                    mode, DI, CMBframe); return py::array(Dn.size(), Dn.data());}, "dist_derivs"_a, "xcmb"_a, "The"_a, "betac"_a, "muc"_a, 
                    "Te_order"_a=10, "betac_order"_a=2, "mode"_a="all", "DI"_a=true, "CMBframe"_a=true,
                    "A function to calculate the SZ signal from an asymptotic expansion around Te, for a given set of background "
                    "photon distribution derivatives i.e. x^k x^k/dx^k(background light) from variables. See the documentation or "
                    "Lee & Chluba 2021 for more details.");

    m.def("null_of_combined_signal", [](const std::function<double(int, double)> &Ddist, Parameters fp)
                    {return compute_null_of_combined_signal(fp, Ddist);}, "dist_derivs"_a, "Params"_a, "A function to calculate the null of the "
                    "combined typical CMB backlit asymptotic signals and the imputted photon distribution. The distribution is imputted by derivatives "
                    "i.e. x^k x^k/dx^k(background light)).");

    m.def("null_of_combined_signal_radio", [](Parameters fp) {return compute_null_of_combined_signal(fp);}, "Params"_a, 
                    "A function to calculate the null of the combined radio distortion and the typical CMB backlit asymptotic signals.");

    m.def("PowerLawDerivs", &PowerLawDerivs, "k"_a, "x"_a, "gamma"_a, "The relevant derivatives for a power law, i.e., x^k d^k/dx^k (x^-gamma)");
}