#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/functional.h>
#include "SZ_multipoleKernel.h"

using namespace pybind11::literals;
namespace py = pybind11;

void init_ex_kernel(py::module_ &m){
//Calculating distortions through the multipole kernel
    m.def("distortion", [](const std::function<double(double)> &eDist, Parameters fp, bool DI){
                    vector<double> Dn; compute_SZ_distortion_kernel(Dn, fp, DI, eDist);
                    return py::array(Dn.size(), Dn.data());}, "electronDist"_a, "Params"_a, "DI"_a=true,
                    "A function to calculate the SZ signal from a given electron distribution.");
    
    m.def("averaged_kernel", [](const std::function<double(double)> &eDist, Parameters fp, int l){
                    if (l>= 0) {fp.kernel.l = l;}
                    vector<double> K; compute_averaged_kernel(K, fp, eDist);
                    return py::array(K.size(), K.data());}, "electronDist"_a, "Params"_a, "l"_a=-1,
                    "A function to calculate electron distribution averaged kernel. If l is not set, l is used from Params.");
    
    m.def("distortion_fixed_momentum", [](double eta, Parameters fp, bool DI){
                    vector<double> Dn; compute_SZ_distortion_fixed_eta(Dn, fp, DI, eta);
                    return py::array(Dn.size(), Dn.data());}, "eta"_a, "Params"_a, "DI"_a=true,
                    "A function to calculate the SZ signal for a fixed momentum electron population.");

                    
    m.def("distortion_from_variables", [](const std::function<double(double)> &eDist, vector<double> xcmb, double betac, double muc, 
                    string mode, double eps_Int, bool DI){
                    vector<double> Dn; compute_SZ_distortion_kernel(Dn, xcmb, betac, muc, eps_Int, DI, eDist, mode);
                    return py::array(Dn.size(), Dn.data());}, "electronDist"_a, "xcmb"_a, "betac"_a, "muc"_a, 
                    "mode"_a, "relative_accuracy"_a=1e-4, "DI"_a=true,
                    "A function to calculate the SZ signal from a given electron distribution from variables.");
    
    m.def("averaged_kernel_from_variables", [](const std::function<double(double)> &eDist, vector<double> s_array, int l, double eps_Int){
                    vector<double> K; compute_averaged_kernel(K, l, s_array, eps_Int , eDist);
                    return py::array(K.size(), K.data());}, "electronDist"_a, "s_array"_a, "l"_a, "relative_accuracy"_a=1e-4,
                    "A function to calculate electron distribution averaged kernel from variables.");
    
    m.def("distortion_fixed_momentum_from_variables", [](double eta, vector<double> xcmb, double betac, double muc, string mode, double eps_Int, 
                    bool DI) {vector<double> Dn; compute_SZ_distortion_fixed_eta(Dn, xcmb, eta, betac, muc, eps_Int, DI, mode);
                    return py::array(Dn.size(), Dn.data());}, "eta"_a, "xcmb"_a, "betac"_a, "muc"_a, "mode"_a, "relative_accuracy"_a=1e-4, 
                    "DI"_a=true, "A function to calculate the SZ signal for a fixed momentum electron population from variables.");

//Calculating distortions through the beam kernel
    m.def("beam_distortion", [](const std::function<double(double)> &eDist, double mup, Parameters fp, bool DI){
                    vector<double> Dn; compute_SZ_distortion_beam_kernel(Dn, fp, DI, mup, eDist);
                    return py::array(Dn.size(), Dn.data());}, "electronDist"_a, "mup"_a, "Params"_a, "DI"_a=true,
                    "A function to calculate the SZ signal for a beam-like given electron distribution.");

    m.def("beam_averaged_kernel", [](const std::function<double(double)> &eDist, double mup, Parameters fp, int l){
                    if (l>= 0) {fp.kernel.l = l;}
                    vector<double> K; compute_averaged_beam_kernel(K, fp, mup, eDist);
                    return py::array(K.size(), K.data());}, "electronDist"_a, "mup"_a, "Params"_a, "l"_a=-1,
                    "A function to calculate electron distribution averaged beam kernel. If l is not set, l is used from Params.");

    m.def("distortion_fixed_momentum", [](double eta, double mup, Parameters fp, bool DI){
                    vector<double> Dn; compute_SZ_distortion_beam_kernel_fixed_eta(Dn, fp, DI, mup, eta);
                    return py::array(Dn.size(), Dn.data());}, "eta"_a, "mup"_a, "Params"_a, "DI"_a=true,
                    "A function to calculate the SZ signal for a beam-like fixed momentum electron population.");


    m.def("beam_distortion_from_variables", [](const std::function<double(double)> &eDist, double mup, vector<double> xcmb, 
                    double betac, double muc, string mode, double eps_Int, bool DI){
                    vector<double> Dn; compute_SZ_distortion_beam_kernel(Dn, xcmb, mup, betac, muc, eps_Int, DI, eDist, mode);
                    return py::array(Dn.size(), Dn.data());}, "electronDist"_a, "mup"_a, "xcmb"_a, "betac"_a, "muc"_a, "mode"_a,
                    "relative_accuracy"_a=1e-4, "DI"_a=true,
                    "A function to calculate the SZ signal for a beam-like given electron distribution from variables.");

    m.def("beam_averaged_kernel_from_variables", [](const std::function<double(double)> &eDist, double mup, vector<double> s_array, 
                    int l, double eps_Int){ 
                    vector<double> K; compute_averaged_beam_kernel(K, l, s_array, mup, eps_Int, eDist);
                    return py::array(K.size(), K.data());}, "electronDist"_a, "mup"_a, "s_array"_a, "l"_a, "relative_accuracy"_a=1e-4,
                    "A function to calculate electron distribution averaged beam kernel from variables.");

    m.def("distortion_fixed_momentum", [](double eta, double mup, vector<double> xcmb, double betac, double muc, string mode, 
                    double eps_Int, bool DI){
                    vector<double> Dn; compute_SZ_distortion_beam_kernel_fixed_eta(Dn, xcmb, mup, eta, betac, muc, eps_Int, DI, mode);
                    return py::array(Dn.size(), Dn.data());}, "eta"_a, "mup"_a, "xcmb"_a, "betac"_a, "muc"_a, "mode"_a, 
                    "relative_accuracy"_a=1e-4, "DI"_a=true,
                    "A function to calculate the SZ signal for a beam-like fixed momentum electron population from variables.");

// Multipole Kernel Implementation
    m.def("kernel_formula", [](int l, vector<double> s_array, double eta) {
                    int np = s_array.size(); vector<double> K; K.resize(np); 
                    MultipoleKernel MK = MultipoleKernel(l, s_array[0], eta);
                    K[0] = MK.Calculate_formula();
                    for (int var = 1; var < np; var++){
                        MK.Update_s(s_array[var]); K[var] = MK.Calculate_formula();
                    } return py::array(K.size(), K.data());},
                    "l"_a, "s_array"_a, "eta"_a, "A function to calculate the multipole kernel from the analytic "
                    "formula. (See Lee et al. 2021 for more details)");
    
    m.def("kernel_integrated", [](int l, vector<double> s_array, double eta, double Int_eps) {
                    int np = s_array.size(); vector<double> K; K.resize(np); 
                    MultipoleKernel MK = MultipoleKernel(l, s_array[0], eta, Int_eps);
                    K[0] = MK.Calculate_integrated();
                    for (int var = 1; var < np; var++){
                        MK.Update_s(s_array[var]); K[var] = MK.Calculate_integrated();
                    } return py::array(K.size(), K.data());},
                    "l"_a, "s_array"_a, "eta"_a, "relative_accuracy"_a=1.0e-4, "A function to calculate the multipole "
                    "kernel through integration. (See Lee et al. 2021 for more details)");

    m.def("kernel_stable", [](int l, vector<double> s_array, double eta, double Int_eps) {
                    int np = s_array.size(); vector<double> K; K.resize(np); 
                    MultipoleKernel MK = MultipoleKernel(l, s_array[0], eta, Int_eps);
                    K[0] = MK.Calculate_stable();
                    for (int var = 1; var < np; var++){
                        MK.Update_s(s_array[var]); K[var] = MK.Calculate_stable();
                    } return py::array(K.size(), K.data());},
                    "l"_a, "s_array"_a, "eta"_a, "relative_accuracy"_a=1.0e-4, "A function to calculate the multipole "
                    "kernel. At low eta this function calculates through the integrated method and above that the "
                    "formula method to ensure numerical stability up to l = 10. (See Lee et al. 2021 for more details)");
                
    m.def("s_limits", [](double eta){ MultipoleKernel MK = MultipoleKernel(0, 0.0, eta); return MK.s_limits();},
                    "eta"_a, "A function to calculate the minimum and maximum s values that can be scattered to, given an "
                    "input electron energy.");

//Beam Kernel Implementation
    m.def("beam_kernel_formula", [](int l, vector<double> s_array, double eta, double mup){
                    int np = s_array.size(); vector<double> K; K.resize(np);
                    BeamKernel BK = BeamKernel(l, s_array[0], eta, mup);
                    K[0] = BK.Calculate_formula();
                    for (int var = 1; var < np; var++){
                        BK.Update_s(s_array[var]); K[var] = BK.Calculate_formula();
                    } return py::array(K.size(), K.data());},
                    "l"_a, "s_array"_a, "eta"_a, "mup"_a, "A function to calculate the beam kernel from the analytic "
                    "formula. (See Lee et al. 2021 for more details)");

    m.def("beam_s_limits", [](double eta, double mup){ BeamKernel BK = BeamKernel(0, 0.0, eta, mup); 
                    return BK.s_limits();}, "eta"_a, "mup"_a,
                    "A function to calculate the minimum and maximum s values that can be scattered to, given an input "
                    "electron energy and angle.");

    m.def("beam_mup_limits", [](double eta, double s){ BeamKernel BK = BeamKernel(0, s, eta, 0.0);
                    return BK.mup_limits();}, "eta"_a, "s"_a,
                    "A function to calculate the minimum and maximum mup values that can be scattered to, given an "
                    "input energy and energy ratio.");
}