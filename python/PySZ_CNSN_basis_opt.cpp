#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "SZ_CNSN_basis.opt.h"

using namespace pybind11::literals;
namespace py = pybind11;

void init_ex_CNSNopt(py::module_ &m){
    py::class_<IntegralCNSNopt>(m, "class_CNSNopt", py::module_local())
        .def(py::init<>())
        .def("set_from_params", [](IntegralCNSNopt &a, double x, Parameters fp, bool CMBframe){
                        a = IntegralCNSNopt(x, fp, CMBframe);}, "x"_a, "Params"_a, "CMBframe"_a, 
                        "A function to initialise the class from a parameters class.")
        .def("set_from_variables", [](IntegralCNSNopt &a, double x_i, double The_i, double betac_i, double muc_i, 
                        int kmax_i, int accuracy_level_i, int betac_order_i, bool CMBframe){
                        a = IntegralCNSNopt(x_i, The_i, betac_i, muc_i, kmax_i, accuracy_level_i, betac_order_i, CMBframe);},
                        "x"_a, "The"_a, "betac"_a, "muc"_a, "kmax"_a, "accuracy_level"_a, "betac_order"_a, "CMBframe"_a, 
                        "A function to initialise the class from variables.")
        .def("set_for_basis_functions_only", [](IntegralCNSNopt &a, double x_i, int region, int Te_order_i) {
                        a = IntegralCNSNopt(x_i, region, Te_order_i);}, "x"_a, "region"_a, "Te_order"_a,
                        "A function to initialise the class when it is only wanted for calculating the basis functions, Y, D, Q and Mcorr.")
                        
        .def("update_x", &IntegralCNSNopt::Update_x, "x"_a, "A function to set the value of x used in the function.")
        .def("compute_distortion", &IntegralCNSNopt::compute_distortion, "runmode"_a, 
                        "A function to calculate the CNSNopt signal for the initialised values.")

        .def_property_readonly("Y", [](const IntegralCNSNopt &a) {return a.Y;}, "These are the basis functions for the monopole, Y_n as "
                        "defined in Chluba et al. 2012. These are exposed for completeness, but will rarely be necessary to use.")
        .def_property_readonly("D", [](const IntegralCNSNopt &a) {return a.D;}, "These are the basis functions for the dipole, D^{kin}_n as "
                        "defined in Chluba et al. 2012. These are exposed for completeness, but will rarely be necessary to use.")
        .def_property_readonly("Q", [](const IntegralCNSNopt &a) {return a.Q;}, "These are the basis functions for the quadrupole, Q^{kin}_n as "
                        "defined in Chluba et al. 2012. These are exposed for completeness, but will rarely be necessary to use.")
        .def_property_readonly("Mcorr", [](const IntegralCNSNopt &a) {return a.Mcorr;}, "These are the basis functions for the monopole correction"
                        ", Y^{kin}_n as defined in Chluba et al. 2012. These are exposed for completeness, but will rarely be necessary to use.")

        .def("__repr__", [](const Parameters &a) { return "<SZpack.CNSNopt.class_CNSNopt>";})
    ;

    m.def("compute_from_variables", [](vector<double> xcmb, double The, double betac, double muc, int kmax, int accuracy_level, 
                        int betac_order, string mode, bool DI, bool CMBframe) {vector<double> Dn; 
                        compute_SZ_distortion_CNSN_basis_opt(Dn, xcmb, The, betac, muc, kmax, betac_order, mode, accuracy_level, DI, CMBframe); 
                        return py::array(Dn.size(), Dn.data());}, "xcmb"_a, "The"_a, "betac"_a, "muc"_a, "kmax"_a, "accuracy_leve"_a,
                        "betac_order"_a, "mode"_a, "DI"_a=true, "CMBframe"_a=true, 
                        "A function to calculate the CNSNopt signal for the given values.");

    m.def("find_optimal_kmax", [](int accuracy_level, double Te_max) { int kmax = 0, iregmax = 0; 
                        determine_optimal_kmax(accuracy_level, Te_max, kmax, iregmax);
                        return kmax;}, "accuracy_level"_a, "Te_max"_a,
                        "A function to find the optimal value for kmax given an accuracy level goal and a maximum temperature "
                        "(Te_max) that will be used.");
}