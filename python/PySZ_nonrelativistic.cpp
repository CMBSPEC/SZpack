#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "SZ_nonrelativistic.h"
#include "Parameters.h"

using namespace pybind11::literals;
namespace py = pybind11;

void init_ex_nonrel(py::module_ &m){
    py::class_<IntegralNonRelativistic>(m, "class_nonrelativistic", py::module_local())
        .def(py::init<>())
        .def("set_from_params", [](IntegralNonRelativistic &a, double x, Parameters fp){
                        a = IntegralNonRelativistic(x, fp);}, "x"_a, "Params"_a, 
                        "A function to initialise the class from a parameters class.")
        .def("set_from_variables", [](IntegralNonRelativistic &a, double x_i, double The_i, double betac_i, double muc_i){
                        a = IntegralNonRelativistic(x_i, The_i, betac_i, muc_i);},
                        "x"_a, "The"_a, "betac"_a, "muc"_a, "A function to initialise the class from variables.")

        .def("update_x", &IntegralNonRelativistic::Update_x, "x"_a, "A function to set the value of x used in the function.")
        .def("compute_distortion", &IntegralNonRelativistic::compute_distortion, "runmode"_a, 
                        "A function to calculate the nonrelativistic signal for the initialised values.")

        .def("__repr__", [](const Parameters &a) { return "<SZpack.nonrelativistic.class_nonrelativistic>";})
    ;

    m.def("compute_from_variables", [](vector<double> xcmb, double The, double betac, double muc, string mode, bool DI) 
                        {vector<double> Dn; compute_SZ_distortion_nonrelativistic(Dn, xcmb, The, betac, muc, mode, DI); 
                        return py::array(Dn.size(), Dn.data());}, "xcmb"_a, "The"_a, "betac"_a, "muc"_a, "mode"_a, "DI"_a=true, 
                        "A function to calculate the nonrelativistic signal for the given values.");
}