#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "SZpack.h"

using namespace pybind11::literals;
namespace py = pybind11;

void init_ex_params(py::module_ &);
void init_ex_asymptotic(py::module_ &);
void init_ex_CNSNopt(py::module_ &);
void init_ex_CNSN(py::module_ &);
void init_ex_3d(py::module_ &);
void init_ex_5d(py::module_ &);
void init_ex_nonrel(py::module_ &);

PYBIND11_MODULE(SZpack, m) {
    m.doc() = "The python wrapper for SZpack.";
    init_ex_params(m);

    m.def("compute_non_relativistic", [](Parameters fp, bool DI) {vector<double> Dn; compute_signal_nonrelativistic(Dn, fp, DI); 
                        return py::array(Dn.size(), Dn.data());}, "Params"_a, "DI"_a=true,
                        "A function to calculate the signal with no relativistic corrections. "
                        "Note setting Te = 0 in general will return a result of 0, since Te is being wrapped into y.");
    m.def("compute_5d", [](Parameters fp, bool DI) {vector<double> Dn; compute_signal_5D(Dn, fp, DI);
                        return py::array(Dn.size(), Dn.data());}, "Params"_a, "DI"_a=true,
                        "A function to calculate the signal by computing all 5 integrals numerically.");
    m.def("compute_3d", [](Parameters fp, bool DI) {vector<double> Dn; compute_signal_3D(Dn, fp, DI);
                        return py::array(Dn.size(), Dn.data());}, "Params"_a, "DI"_a=true,
                        "A function to calculate the signal signal by precomputing two of the integrals analytically, "
                        "under an assumption of low values for betac.");
    m.def("compute_asymptotic", [](Parameters fp, bool DI, bool CMBframe) {vector<double> Dn; compute_signal_asymptotic(Dn, fp, DI, CMBframe);
                        return py::array(Dn.size(), Dn.data());}, "Params"_a, "DI"_a=true, "CMBframe"_a=false,
                        "A function to calculate the signal from an asymptotic expansion around Te. This setting is best for Te < 2 keV. "
                        "It will throw an error for Te > 20 keV. CMBframe determines if the signal is calculated in the CMB frame or the observer frame.");
    m.def("compute_CNSN", [](Parameters fp, bool DI, bool CMBframe) {vector<double> Dn; compute_signal_CNSN(Dn, fp, DI, CMBframe);
                        return py::array(Dn.size(), Dn.data());}, "Params"_a, "DI"_a=true, "CMBframe"_a=false,
                        "A function to calculate the signal from a pivot based interpolation scheme described in Chluba et al. 2012. "
                        "This method will only calculate for 2 keV < Te < 75 keV. CMBframe determines if the signal is calculated in the "
                        "CMB frame or the observer frame.");
    m.def("compute_CNSN_opt", [](Parameters fp, bool DI, bool CMBframe) {vector<double> Dn; compute_signal_CNSN_opt(Dn, fp, DI, CMBframe);
                        return py::array(Dn.size(), Dn.data());}, "Params"_a, "DI"_a=true, "CMBframe"_a=false,
                        "A function to calculate the signal from a slightly more involved form of the pivot based interpolation, summarised in "
                        "Chluba et al. 2013. The precision of this method depends on setting for kmax and accuracy level. "
                        "CMBframe determines if the signal is calculated in the CMB frame or the observer frame.");
    m.def("compute_combo", [](Parameters fp, bool DI, bool CMBframe) {vector<double> Dn; compute_signal_combo(Dn, fp, DI, CMBframe);
                        return py::array(Dn.size(), Dn.data());}, "Params"_a, "DI"_a=true, "CMBframe"_a=false,
                        "A function to calculate the signal using the asymptotic mode for Te < 2 keV and CNSN above this. "
                        "CMBframe determines if the signal is calculated in the CMB frame or the observer frame.");
    m.def("compute_precise", [](Parameters fp, bool DI) {vector<double> Dn; compute_signal_precise(Dn, fp, DI);
                        return py::array(Dn.size(), Dn.data());}, "Params"_a, "DI"_a=true,
                        "A function to calculate the signal using the asymptotic mode for Te < 2 keV, the 3d method for Te < 75 keV and "
                        "the CNSN method otherwise. CMBframe determines if the signal is calculated in the CMB frame or the observer frame.");
    m.def("compute_means", [](Parameters fp, bool yw, bool DI) {vector<double> Dn; compute_signal_means(Dn, fp, DI, yw);
                        return py::array(Dn.size(), Dn.data());}, "Params"_a, "yw"_a, "DI"_a=true,
                        "A function to calculate the signal using moments to account for variations of temperature and velocity within the line of sight. "
                        "This function depends on params.means_omega, params.means_sigma and params.means_kappa. The method is described in detail in "
                        "Chluba et al. 2013. yw determines whether a the moments are calculated using a y- or tau-weighted method.");
    m.def("compute_means_ex", [](Parameters fp, bool yw, bool DI) {vector<double> Dn; compute_signal_means_ex(Dn, fp, DI);
                        return py::array(Dn.size(), Dn.data());}, "Params"_a, "yw"_a, "DI"_a=true,
                        "A function to calculate the signal using moments to account for variations of temperature and velocity within the line of sight. "
                        "This goes to higher order than the means method -- up to The^4 terms for the temperature dispersion and The^3 for the temperature-"
                        "velcoity cross terms. This function depends on params.means_omegas, params.means_sigmas and params.means_kappa. The method is "
                        "described in detail in Chluba et al. 2013. yw determines whether a the moments are calculated using a y- or tau-weighted method.");
    m.def("compute_two_temperatures", [](double ftau, double DT_T, Parameters fp, bool DI, bool CMBframe){vector<double> Dn; 
                        compute_signal_TwoTemperatures(Dn, ftau, DT_T, fp, DI, CMBframe); return py::array(Dn.size(), Dn.data());}, 
                        "ftau"_a, "DT_T"_a, "Params"_a, "DI"_a=true, "CMBframe"_a=false,
                        "A function to calculate the combined signal of two overlaying thermal distributions of electrons. ftau determines the fraction of" 
                        "optical depth contributed by the second temperature, i.e., DTau1/(Dtau0+Dtau1) and DT_T the proportional difference between the "
                        "temperatures. i.e., (T1-T0)/T0");
    m.def("compute_relativistic_corrections", [](Parameters fp, bool DI) {vector<double> Dn; compute_signal_RelCorrs(Dn, fp, DI);
                        return py::array(Dn.size(), Dn.data());}, "Params"_a, "DI"_a=true,
                        "A function to calculate the relativistic corrections. This calculates the differences between the non-relativistic "
                        "method and the combo method, in the CMB frame.");
    m.def("compute_temperature_dispersion", [](Parameters fp, bool DI) {vector<double> Dn; compute_signal_TDispersion(Dn, fp, DI);
                        return py::array(Dn.size(), Dn.data());}, "Params"_a, "DI"_a=true,
                        "A function to calculate the temperature dispersion. This calculates the difference between the means signal given the omega, "
                        "sigma and kappa values that have been set, and the signal without these variables, all calculated in the CMB frame");
    m.def("compute_null", [](Parameters fp, bool inGHz) {double null = (inGHz?fp.rare.x_nu_conversion():1.0)*compute_null_of_SZ_signal(fp); return null;},
                        "Params"_a, "inGHz"_a=false, "A function to find the null crossing point for any given parameter set");


    m.def("compute_combo_from_variables", [](vector<double> xcmb, double The, double betac, double muc, int Te_order, int betac_order, 
                        string mode, bool DI, bool CMBframe){
                            Parameters fp = Parameters();
                            fp.xcmb = xcmb; fp.xmin = xcmb[0]; fp.xmax = xcmb[-1]; fp.gridpoints = xcmb.size();
                            fp.updateT(The*const_me);
                            fp.betac = betac; fp.muc = muc;
                            fp.Dtau = 1.0;
                            fp.T_order = Te_order; fp.beta_order = betac_order;
                            fp.rare.RunMode = mode;

                            vector<double> Dn; compute_signal_combo(Dn, fp, DI, CMBframe); return py::array(Dn.size(), Dn.data());
                        }, "xcmb"_a, "The"_a, "betac"_a, "muc"_a, "Te_order"_a, "betac_order"_a, "mode"_a, "DI"_a=true, "CMBframe"_a=true,
                        "A function to calculate the signal from the given variables using the asymptotic mode for Te < 2 keV and CNSN above this.");

    m.def("compute_precise_from_variables", [](vector<double> xcmb, double The, double betac, double muc, string mode, bool DI){
                            Parameters fp = Parameters();
                            fp.xcmb = xcmb; fp.xmin = xcmb[0]; fp.xmax = xcmb[-1]; fp.gridpoints = xcmb.size();
                            fp.updateT(The*const_me);
                            fp.betac = betac; fp.muc = muc;
                            fp.Dtau = 1.0;
                            fp.T_order = 10; fp.beta_order = 2; fp.relative_accuracy = 1e-6;
                            fp.rare.RunMode = mode;

                            vector<double> Dn; compute_signal_precise(Dn, fp, DI); return py::array(Dn.size(), Dn.data());
                        }, "xcmb"_a, "The"_a, "betac"_a, "muc"_a, "mode"_a, "DI"_a=true, "A function to calculate the signal from the given "
                        "variables using the asymptotic mode for Te < 2 keV, the 3d method for Te < 75 keV and the CNSN method otherwise.");

    m.def("compute_means_from_variables", [](vector<double> xcmb, double The, double betac, double muc, string mode, 
                        double omega, double sigma, double kappa, bool yw, bool DI){
                            Parameters fp = Parameters();
                            fp.xcmb = xcmb; fp.xmin = xcmb[0]; fp.xmax = xcmb[-1]; fp.gridpoints = xcmb.size();
                            fp.updateT(The*const_me);
                            fp.betac = betac; fp.muc = muc;
                            fp.Dtau = 1.0;
                            fp.means.Omega = omega; fp.means.Sigma = sigma; fp.means.kappa = kappa;
                            fp.rare.RunMode = mode;

                            vector<double> Dn; compute_signal_means(Dn, fp, DI, yw); return py::array(Dn.size(), Dn.data());
                        }, "xcmb"_a, "The"_a, "betac"_a, "muc"_a, "mode"_a, "omega"_a, "sigma"_a, "kappa"_a, "yw"_a, "DI"_a=true, 
                        "A function to calculate the signal from the given variables using moments to account for variations of temperature and "
                        "velocity within the line of sight. The method is described in detail in Chluba et al. 2013. yw determines whether a the "
                        "moments are calculated using a y- or tau-weighted method.");

    m.def("compute_means_ex_from_variables", [](vector<double> xcmb, double The, double betac, double muc, string mode, 
                        vector<double> omegas, vector<double> sigmas, double kappa, bool yw, bool DI){
                            if ((omegas.size() != 3)||(sigmas.size() != 3)){
                                print_error("Error: Must input omegas and sigmas as np arrays of length 3.");
                                vector<double> Dn = {}; return py::array(Dn.size(), Dn.data());
                            }
                                
                            Parameters fp = Parameters();
                            fp.xcmb = xcmb; fp.xmin = xcmb[0]; fp.xmax = xcmb[-1]; fp.gridpoints = xcmb.size();
                            fp.updateT(The*const_me);
                            fp.betac = betac; fp.muc = muc;
                            fp.Dtau = 1.0;
                            fp.means.assignOmegas(omegas[0], omegas[1], omegas[2]);
                            fp.means.assignSigmas(sigmas[0], sigmas[1], sigmas[2]); 
                            fp.means.kappa = kappa;
                            fp.rare.RunMode = mode;

                            vector<double> Dn; compute_signal_means_ex(Dn, fp, DI, yw); return py::array(Dn.size(), Dn.data());
                        }, "xcmb"_a, "The"_a, "betac"_a, "muc"_a, "mode"_a, "omega"_a, "sigma"_a, "kappa"_a, "yw"_a, "DI"_a=true,
                        "A function to calculate the signal from the given variables using moments to account for variations of temperature and velocity "
                        "within the line of sight. This goes to higher order than the means method -- up to The^4 terms for the temperature dispersion "
                        "and The^3 for the temperature-velcoity cross terms. The method is described in detail in Chluba et al. 2013. yw determines "
                        "whether a the moments are calculated using a y- or tau-weighted method.");


    m.def("Dcompute_combo_CMB", [](Parameters fp, int dThe, int dbeta_para, int dbeta2_perp, bool yw, bool DI) {
                        vector<double> Dn; Dcompute_signal_combo_CMB(Dn, fp, dThe, dbeta_para, dbeta2_perp, yw, DI); 
                        return py::array(Dn.size(), Dn.data());}, "Params"_a, "dThe"_a, "dbeta_para"_a, "dbeta2_perp"_a,
                        "yw"_a=false, "DI"_a=true,
                        "A function to calculate, for every point along for the associated signal, the derivatives with respect to The, betac_para and "
                        "betac2_perp. betac_para is beta in the direction of the cluster movement i.e., betac*muc and betac2_perp is beta "
                        "perpendicular to this. i.e., betac^2*(1-muc^2). The input integers give the order of the d^i/dThe^i, d^j/dbetac_para^j and "
                        "d^k/dbetac_perp^k derivatives. dThe<=4, dbeta_para<=2 and dbeta2_perp<=1. If dThe=dbeta_para=dbeta2_perp=0 the orignal signal "
                        "is returned.");
    m.def("Dcompute_combo_for_x", [](Parameters fp, int dx) { vector<double> Dn; Dcompute_signal_combo_for_x(Dn, fp, dx); 
                        return py::array(Dn.size(), Dn.data());}, "Params"_a, "dx"_a, 
                        "A function to calculate, for every point along the associated signal, the derivatives with respect to x. i.e., d^i/dx^i of the "
                        "signal. 0<=dx<=4. If dx=0 the original signal is returned.");


    m.def("Dcompute_from_variables", [](vector<double> xcmb, int dThe, int dbeta_para, int dbeta2_perp, 
                                        double The, double betac, double muc, bool yw, bool DI){
                            Parameters fp = Parameters();
                            fp.xcmb = xcmb; fp.xmin = xcmb[0]; fp.xmax = xcmb[-1]; fp.gridpoints = xcmb.size();
                            fp.updateT(The*const_me);
                            fp.betac = betac; fp.muc = muc;
                            fp.Dtau = 1.0;
                            vector<double> Dn; Dcompute_signal_combo_CMB(Dn, fp, dThe, dbeta_para, dbeta2_perp, yw, DI);
                            return py::array(Dn.size(), Dn.data());
                        }, "xcmb"_a, "dThe"_a, "dbeta_para"_a, "dbeta2_perp"_a, "The"_a, "betac"_a, "muc"_a, "yw"_a=false, "DI"_a=true,
                        "A function to calculate for the given variables, for every point along for the associated signal, the derivatives with respect "
                        "to The, betac_para and betac2_perp. betac_para is beta in the direction of the cluster movement i.e., betac*muc and betac2_perp "
                        "is beta perpendicular to this. i.e., betac^2*(1-muc^2). The input integers give the order of the d^i/dThe^i, d^j/dbetac_para^j "
                        "and d^k/dbetac_perp^k derivatives. dThe<=4, dbeta_para<=2 and dbeta2_perp<=1. If dThe=dbeta_para=dbeta2_perp=0 the orignal "
                        "signal is returned.");

    m.def("Dcompute_for_x_from_variables", [](vector<double> xcmb, int dx, double The, double betac, double muc){
                            Parameters fp = Parameters();
                            fp.xcmb = xcmb; fp.xmin = xcmb[0]; fp.xmax = xcmb[-1]; fp.gridpoints = xcmb.size();
                            fp.updateT(The*const_me);
                            fp.betac = betac; fp.muc = muc;
                            fp.Dtau = 1.0;
                            vector<double> Dn; Dcompute_signal_combo_for_x(Dn, fp, dx); return py::array(Dn.size(), Dn.data());
                        }, "xcmb"_a, "dx"_a, "The"_a, "betac"_a, "muc"_a, 
                        "A function to calculate for the given variables, for every point along the associated signal, the derivatives with respect to x. "
                        "i.e., d^i/dx^i of the signal. 0<=dx<=4. If dx=0 the original signal is returned.");

    m.def("convert_signal_DT", [](vector<double> xcmb, vector<double> Dn){vector<double> DT; convert_signal_DT(DT, xcmb, Dn);
                        return py::array(DT.size(), DT.data());}, "xcmb"_a, "Dn"_a,
                        "A function to convert the Dn signal to DT/TCMB i.e., (the effective change CMB temperature)/(CMB temperature) at a given "
                        "frequency. This function inverts the boltzmann distribution to calculate the conversion accurately. " 
                        "Note this takes in Dn and not DI.");
    
    m.def("convert_signal_DT_approx", [](vector<double> xcmb, vector<double> Dn){vector<double> DT; convert_signal_DT_approx(DT, xcmb, Dn);
                        return py::array(DT.size(), DT.data());}, "xcmb"_a, "Dn"_a,
                        "A function to convert the Dn signal to DT/TCMB i.e., (the effective change CMB temperature)/(CMB temperature) at a given "
                        "frequency. This function uses a conventional approximation from the first derivative of the boltzmann distribution. "
                        "Note this takes in Dn and not DI.");

    py::module asym = m.def_submodule("Asymptotic", "The submodule for some specific asymptotic method functions.");
    init_ex_asymptotic(asym);
    py::module cnsnopt = m.def_submodule("CNSN_basis_opt", "The submodule for some specific CNSNopt method functions.");
    init_ex_CNSNopt(cnsnopt);
    py::module cnsn = m.def_submodule("CNSN_basis", "The submodule for some specific CNSN method functions.");
    init_ex_CNSN(cnsn);
    py::module int3d = m.def_submodule("Integral3D", "The submodule for some specific 3D method functions.");
    init_ex_3d(int3d);
    py::module int5d = m.def_submodule("Integral5D", "The submodule for some specific 5D method functions.");
    init_ex_5d(int5d);
    py::module nonrel = m.def_submodule("Nonrelativistic", "The submodule for some specific nonrelativistic method functions.");
    init_ex_nonrel(nonrel);
}