#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "SZ_Integral.5D.h"
#include "Parameters.h"

using namespace pybind11::literals;
namespace py = pybind11;

void init_ex_5d(py::module_ &m){
    py::class_<Integral5D>(m, "class_integral_5D", py::module_local())
        .def(py::init<>())
        .def("set_from_params", [](Integral5D &a, double x, Parameters fp){
                        a = Integral5D(x, fp);}, "x"_a, "Params"_a, 
                        "A function to initialise the class from a parameters class.")
        .def("set_from_variables", [](Integral5D &a, double x_i, double The_i, double betac_i, double muc_i, double eps_Int){
                        a = Integral5D(x_i, The_i, betac_i, muc_i, eps_Int);},
                        "x"_a, "The"_a, "betac"_a, "muc"_a, "relative_accuracy"_a, 
                        "A function to initialise the class from variables.")

        .def("update_x", &Integral5D::Update_x, "x"_a, "A function to set the value of x used in the function.")
        .def("compute_distortion", &Integral5D::compute_distortion, "runmode"_a, 
                        "A function to calculate the 5D integrated signal for the initialised values.")

        .def("__repr__", [](const Parameters &a) { return "<SZpack.integral_5D.class_integral_5D>";})
    ;

    m.def("compute_from_variables", [](vector<double> xcmb, double The, double betac, double muc, double eps_Int, string mode,
                        bool DI) {vector<double> Dn; 
                        compute_SZ_distortion_Patterson_5D(Dn, xcmb, The, betac, muc, eps_Int, mode, DI); 
                        return py::array(Dn.size(), Dn.data());}, "xcmb"_a, "The"_a, "betac"_a, "muc"_a, "relative_accuracy"_a,
                        "mode"_a, "DI"_a=true, "A function to calculate the 5D integrated signal for the given values.");
}