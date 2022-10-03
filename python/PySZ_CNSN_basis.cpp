#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "SZ_CNSN_basis.h"

using namespace pybind11::literals;
namespace py = pybind11;

void init_ex_CNSN(py::module_ &m){
    py::class_<IntegralCNSN>(m, "class_CNSN", py::module_local())
        .def(py::init<>())
        .def("set_from_params", [](IntegralCNSN &a, int k, Parameters fp, bool CMBframe){
                        a = IntegralCNSN(k, fp, CMBframe, true);}, "k"_a, "Params"_a, "CMBframe"_a, 
                        "A function to initialise the class from a parameters class. k is the index of xcmb that will be calculated.")
        .def("set_from_variables", [](IntegralCNSN &a, double x_i, double The_i, double betac_i, double muc_i, 
                        int Te_order_i, int betac_order_i, bool CMBframe){
                        a = IntegralCNSN(x_i, The_i, betac_i, muc_i, Te_order_i, betac_order_i, CMBframe);},
                        "x"_a, "The"_a, "betac"_a, "muc"_a, "Te_order"_a, "betac_order"_a, "CMBframe"_a, 
                        "A function to initialise the class from variables.")
        .def("set_from_variables_highest_order", [](IntegralCNSN &a, double x_i, double The_i, double betac_i, 
                        double muc_i, bool CMBframe){ a = IntegralCNSN(x_i, The_i, betac_i, muc_i, CMBframe);},
                        "x"_a, "The"_a, "betac"_a, "muc"_a,  "CMBframe"_a, 
                        "A function to initialise the class from variables. Te_order and betac_order are set to their highest values.")

                        
        .def("update_x", &IntegralCNSN::Update_x, "x"_a, "A function to set the value of x used in the function.")
        .def("compute_distortion", &IntegralCNSN::compute_distortion, "runmode"_a, 
                        "A function to calculate the CNSN signal for the initialised values.")
        .def("Dcompute_distortion", [](IntegralCNSN &a, int dThe, int dbeta_para, int dbeta2_perp){
                        vector<double> dDn; a.Dcompute_distortion(dThe, dbeta_para, dbeta2_perp, dDn); 
                        return py::array(dDn.size(), dDn.data());}, "dThe"_a, "dbeta_para"_a, "dbeta2_perp"_a,
                        "A function to calculate the derivatives at this point with respect to The, betac_para and betac2_perp. "
                        "betac_para is beta in the direction of the cluster movement i.e., betac*muc and betac2_perp is beta "
                        "perpendicular to this. i.e., betac^2*(1-muc^2). The input integers give the order of d/dThe that the will "
                        "be calculated up to, and then the order of the d^j/dbetac_para^j and d^k/dbetac_perp^k derivatives.")

        .def("__repr__", [](const Parameters &a) { return "<SZpack.CNSN.class_CNSN>";})
    ;

    m.def("compute_from_variables", [](vector<double> xcmb, double The, double betac, double muc, int Te_order, int betac_order,
                        string mode, bool DI, bool CMBframe) {vector<double> Dn; 
                        compute_SZ_distortion_CNSN_basis(Dn, xcmb, The, betac, muc, Te_order, betac_order, mode, DI, CMBframe); 
                        return py::array(Dn.size(), Dn.data());}, "xcmb"_a, "The"_a, "betac"_a, "muc"_a, "Te_order"_a, "betac_order"_a,
                        "mode"_a, "DI"_a=true, "CMBframe"_a=true, 
                        "A function to calculate the CNSN signal for the given values.");

    m.def("Dcompute_from_variables", [](double x, int dThe, int dbeta_para, int dbeta2_perp, 
                        double The, double betac, double muc, bool CMBframe) { vector<double> dDn_dThe; 
                        Dcompute_SZ_distortion_CNSN_basis(x, dThe, dbeta_para,dbeta2_perp, The, betac, muc, dDn_dThe, CMBframe); 
                        return py::array(dDn_dThe.size(), dDn_dThe.data());}, "x"_a, "dThe"_a, "dbeta_para"_a, "dbeta2_perp"_a,
                        "The"_a, "betac"_a, "muc"_a, "CMBframe"_a=true,
                        "A function to calculate the derivatives for the given values with respect to The, betac_para and betac2_perp. "
                        "betac_para is beta in the direction of the cluster movement i.e., betac*muc and betac2_perp is beta "
                        "perpendicular to this. i.e., betac^2*(1-muc^2). The input integers give the order of d/dThe that the will "
                        "be calculated up to, and then the order of the d^j/dbetac_para^j and d^k/dbetac_perp^k derivatives. "
                        "dThe<=4, dbeta_para<=2 and dbeta2_perp<=1.");

    m.def("Dcompute_from_params", [](double x, Parameters fp, int dThe, int dbeta_para, int dbeta2_perp, bool CMBframe) { 
                            fp.D.dThe = dThe; fp.D.dbeta_para = dbeta_para; fp.D.dbeta2_perp = dbeta2_perp;  
                            Dcompute_SZ_distortion_CNSN_basis(x, fp, CMBframe); 
                            return py::array(fp.D.dDn_dThe.size(), fp.D.dDn_dThe.data());
                        }, "x"_a, "Params"_a, "dThe"_a, "dbeta_para"_a, "dbeta2_perp"_a, "CMBframe"_a=false,
                        "A function to calculate the derivatives for the given values with respect to The, betac_para and betac2_perp. "
                        "betac_para is beta in the direction of the cluster movement i.e., betac*muc and betac2_perp is beta "
                        "perpendicular to this. i.e., betac^2*(1-muc^2). The input integers give the order of d/dThe that the will "
                        "be calculated up to, and then the order of the d^j/dbetac_para^j and d^k/dbetac_perp^k derivatives. "
                        "dThe<=4, dbeta_para<=2 and dbeta2_perp<=1.");
}