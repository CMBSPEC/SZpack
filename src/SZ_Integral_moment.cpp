//==================================================================================================
//
// //TODO: Fill this box
//
//==================================================================================================
//
// Author: Elizabeth Lee
// Based off work by Jens Chluba
//
//==================================================================================================

#include "SZ_Integral_moment.h"

using namespace std;

normFunc NormFunc;
integrationValues IValues;

void setNormFunc(normFunc NF){
    NormFunc = NF;
}

double func_1(double lz){ return 1.0; }
double func_1(double lz, double yk){ return 1.0; }

//==================================================================================================
double compute_Integral(double (*func)(double l), double Int_eps, double lmin, double lmax)
{
    double eps_abs = 1.0e-300;
    double r = Integrate_using_Patterson_adaptive(max(-1.0, lmin), min(-0.1, lmax), 
                                                Int_eps, eps_abs, func);
    r += Integrate_using_Patterson_adaptive(max(-0.1, lmin), min(-0.01, lmax), 
                                          Int_eps, eps_abs, func);
    r += Integrate_using_Patterson_adaptive(max(-0.01, lmin), min(-0.001, lmax), 
                                          Int_eps, eps_abs, func);
    r += Integrate_using_Patterson_adaptive(max(-0.001, lmin), min(-0.0001, lmax), 
                                          Int_eps, eps_abs, func);
    r += Integrate_using_Patterson_adaptive(max(-0.0001, lmin), min(-0.00001, lmax), 
                                          Int_eps, eps_abs, func);
    r += Integrate_using_Patterson_adaptive(max(-0.00001, lmin), min(0.00001, lmax), 
                                          Int_eps, eps_abs, func);
    r += Integrate_using_Patterson_adaptive(max(0.00001, lmin), min(0.0001, lmax), 
                                          Int_eps, eps_abs, func);
    r += Integrate_using_Patterson_adaptive(max(0.0001, lmin), min(0.001, lmax), 
                                          Int_eps, eps_abs, func);
    r += Integrate_using_Patterson_adaptive(max(0.001, lmin), min(0.01, lmax), 
                                          Int_eps, eps_abs, func);
    r += Integrate_using_Patterson_adaptive(max(0.01, lmin), min(0.1, lmax), 
                                          Int_eps, eps_abs, func);
    r += Integrate_using_Patterson_adaptive(max(0.1, lmin), min(1.0, lmax), 
                                          Int_eps, eps_abs, func); 
    return r;
}

//==================================================================================================
//
// functions to compute y, z, b, and c integrals according to CSNN 2012
//
//==================================================================================================
integrationValues::integrationValues(){
    lx0 = ly0 = lsc = lang = 0;
    oneD = true;
    usek = false;
    lx = ly = 0;
    k = 0;
    integrand = func_1;
    integrand_yk = func_1;
}

integrationValues::integrationValues(bool OneD){
    integrationValues();
    oneD = OneD;
}

integrationValues::integrationValues(double Lx0, double Ly0, double Lsc, double Lang, Integrand_1D Integrand, bool OneD){
    lx0 = Lx0;
    ly0 = Ly0;
    lsc = Lsc;
    lang = Lang;
    integrand = Integrand;
    oneD = OneD;
    usek = false;
}

integrationValues::integrationValues(double Lx0, double Ly0, double Lsc, double Lang, Integrand_1D_yk Integrand, int K, bool OneD){
    lx0 = Lx0;
    ly0 = Ly0;
    lsc = Lsc;
    lang = Lang;
    k = K;
    integrand_yk = Integrand;
    oneD = OneD;
    usek = true;
}

//==================================================================================================
double func_dtau(double l)
{
    if (IValues.oneD){
        return clusterF_1D.Ne(l, PSP);
    }
    return clusterF_3D.Ne(IValues.lx, IValues.ly, l, PSP);
}

double compute_tau(double zsc)
{
    return zsc*const_sigT*compute_Integral(func_dtau, g_eps_Int, -1.0, 1.0);
}

//==================================================================================================
double func_dySZ(double l)
{
    if (IValues.oneD){
        return clusterF_1D.Ne(l, PSP)*clusterF_1D.Te(l, PSP)/const_me;
    }
    return clusterF_3D.Ne(IValues.lx, IValues.ly, l, PSP)*clusterF_3D.Te(IValues.lx, IValues.ly, l, PSP)/const_me;
}

double compute_ySZ(double zsc)
{ 
    return zsc*const_sigT*compute_Integral(func_dySZ, g_eps_Int, -1.0, 1.0);
}

double func_dyk(double l)
{ 
    double Te;
    if (IValues.oneD){
        Te = clusterF_1D.Te(l, PSP);
        return clusterF_1D.Ne(l, PSP)*pow(Te/const_me, IValues.k + 1);
    }
    Te = clusterF_3D.Te(IValues.lx, IValues.ly, l, PSP);
    return clusterF_3D.Ne(IValues.lx, IValues.ly, l, PSP)*pow(Te/const_me, IValues.k + 1); 
}

double compute_yk(double zsc, int k)
{ 
    IValues.k = k;
    return zsc*const_sigT*compute_Integral(func_dyk, g_eps_Int, -1.0, 1.0);
}

//==================================================================================================
// low temperature moments
//==================================================================================================
double func_dy_k_cond(double l)
{ 
    double Te, D;
    if (IValues.oneD){
        Te = clusterF_1D.Te(l, PSP);
        D = Te/const_me-The_ref_global;
        return NormFunc(Te/const_me)*clusterF_1D.Ne(l, PSP)*pow(D, IValues.k + 1); 
    }
    Te = clusterF_3D.Te(IValues.lx, IValues.ly, l, PSP);
    D = Te/const_me-The_ref_global;
    return NormFunc(Te/const_me)*clusterF_3D.Ne(IValues.lx, IValues.ly, l, PSP)*pow(D, IValues.k + 1); 
}

double compute_y_k_cond(int k, int ir, double zsc, SZ_moment_method SZM)
{ 
    IValues.k = k;
    double r = 0.0;
    for(int z = 0; z<(int) SZM.Te_zeros[ir].size(); z++)
    {
        r += compute_Integral(func_dy_k_cond, g_eps_Int,
                              SZM.Te_zeros[ir][z].lmin, SZM.Te_zeros[ir][z].lmax);
    }
    return zsc*const_sigT*r;
}

//==================================================================================================
double func_db_k_0_cond(double l)
{
    double Te, D, b;
    if (IValues.oneD){
        Te = clusterF_1D.Te(l, PSP);
        b = clusterF_1D.betac(l, PSP);
        D = Te/const_me-The_ref_global;
        return NormFunc(Te/const_me)*clusterF_1D.Ne(l, PSP)*b*b*pow(D, IValues.k); 
    }
    Te = clusterF_3D.Te(IValues.lx, IValues.ly, l, PSP);
    b = clusterF_3D.betac(IValues.lx, IValues.ly, l, PSP);
    D = Te/const_me-The_ref_global;
    return NormFunc(Te/const_me)*clusterF_3D.Ne(IValues.lx, IValues.ly, l, PSP)*b*b*pow(D, IValues.k); 
}

double compute_b_k_0_cond(int k, int ir, double zsc, SZ_moment_method SZM)
{ 
    IValues.k = k;
    double r=0.0;
    
    for(int z=0; z<(int)SZM.Te_zeros[ir].size(); z++)
        
        r+=compute_Integral(func_db_k_0_cond, g_eps_Int,
                            SZM.Te_zeros[ir][z].lmin, SZM.Te_zeros[ir][z].lmax);
    
    return zsc*const_sigT*r;
}

double compute_b_0(double zsc)
{ 
    IValues.k = 0;
    return zsc*const_sigT*compute_Integral(func_db_k_0_cond, g_eps_Int, -1.0, 1.0);
}

//==================================================================================================
double func_db_k_1_cond(double l)
{ 
    double Te, b, muc, D;
    if (IValues.oneD){
        Te = clusterF_1D.Te(l, PSP);
        b = clusterF_1D.betac(l, PSP);
        muc = clusterF_1D.muc(l, PSP);
        D = Te/const_me-The_ref_global;
        return NormFunc(Te/const_me)*clusterF_1D.Ne(l, PSP)*b*muc*pow(D, IValues.k);
    }
    Te = clusterF_3D.Te(IValues.lx, IValues.ly, l, PSP);
    b = clusterF_3D.betac(IValues.lx, IValues.ly, l, PSP);
    muc = clusterF_3D.muc(IValues.lx, IValues.ly, l, PSP);
    D = Te/const_me-The_ref_global;
    return NormFunc(Te/const_me)*clusterF_3D.Ne(IValues.lx, IValues.ly, l, PSP)*b*muc*pow(D, IValues.k);
}

double compute_b_k_1_cond(int k, int ir, double zsc, SZ_moment_method SZM)
{ 
    IValues.k = k;
    double r=0.0;
    
    for(int z=0; z<(int)SZM.Te_zeros[ir].size(); z++)
        
        r+=compute_Integral(func_db_k_1_cond, g_eps_Int,
                            SZM.Te_zeros[ir][z].lmin, SZM.Te_zeros[ir][z].lmax);
    
    return zsc*const_sigT*r;
}

double compute_b_1(double zsc)
{ 
    IValues.k = 0;
    return zsc*const_sigT*compute_Integral(func_db_k_1_cond, g_eps_Int, -1.0, 1.0);
}

//==================================================================================================
double func_db_k_2_cond(double l)
{ 
    double Te, b, muc, D;
    if (IValues.oneD){
        Te = clusterF_1D.Te(l, PSP);
        b = clusterF_1D.betac(l, PSP);
        muc = clusterF_1D.muc(l, PSP);
        D = Te/const_me-The_ref_global;
        return NormFunc(Te/const_me)*clusterF_1D.Ne(l, PSP)*b*b*(1.5*muc*muc-0.5)*pow(D, IValues.k);
    }
    Te = clusterF_3D.Te(IValues.lx, IValues.ly, l, PSP);
    b = clusterF_3D.betac(IValues.lx, IValues.ly, l, PSP);
    muc = clusterF_3D.muc(IValues.lx, IValues.ly, l, PSP);
    D = Te/const_me-The_ref_global;
    return NormFunc(Te/const_me)*clusterF_3D.Ne(IValues.lx, IValues.ly, l, PSP)*b*b*(1.5*muc*muc-0.5)*pow(D, IValues.k);
}

double compute_b_k_2_cond(int k, int ir, double zsc, SZ_moment_method SZM)
{ 
    IValues.k = k;
    double r=0.0;
    
    for(int z=0; z<(int)SZM.Te_zeros[ir].size(); z++)
        
        r+=compute_Integral(func_db_k_2_cond, g_eps_Int,
                            SZM.Te_zeros[ir][z].lmin, SZM.Te_zeros[ir][z].lmax);
    
    return zsc*const_sigT*r;
}

double compute_b_2(double zsc)
{ 
    IValues.k = 0;
    return zsc*const_sigT*compute_Integral(func_db_k_2_cond, g_eps_Int, -1.0, 1.0);
}

//==================================================================================================
//
// additional average over x and y
// 
//==================================================================================================
double beam_function(double Dlx, double Dly, double lang)
{
    double l2=Dlx*Dlx+Dly*Dly;
    return exp(-l2/2.0/lang/lang)/TWOPI/lang/lang;
}

//==================================================================================================
double compute_Int_xy(double ly)
{
    IValues.ly = ly;
    double fac = beam_function(IValues.lx-IValues.lx0, IValues.ly-IValues.ly0, IValues.lang);
    
    if (IValues.usek) {
        return fac*IValues.integrand_yk(IValues.lsc, IValues.k);
    }
    return fac*IValues.integrand(IValues.lsc); 
}

double compute_Int_x(double lx)
{ 
    double a = IValues.ly0 - 6*IValues.lang;
    double b = IValues.ly0 + 6*IValues.lang;
    IValues.lx = lx;
    
    return Integrate_using_Patterson_adaptive(a, b, 1.0e-4, 1.0e-300, compute_Int_xy);
}

//==================================================================================================
double compute_Int_3D(double lx0, double ly0, double lsc, double lang, Integrand_1D integrand)
{ 
    IValues = integrationValues(lx0, ly0, lsc, lang, integrand);
    double a = IValues.lx0 - 6*IValues.lang;
    double b = IValues.lx0 + 6*IValues.lang;
    
    return Integrate_using_Patterson_adaptive(a, b, 2.0e-4, 1.0e-300, compute_Int_x);
}

double compute_Int_3D_yk(double lx0, double ly0, double lsc, double lang, int k,
                         Integrand_1D_yk integrand)
{ 
    IValues = integrationValues(lx0, ly0, lsc, lang, integrand, k);
    double a = IValues.lx0 - 6*IValues.lang;
    double b = IValues.lx0 + 6*IValues.lang;
    
    return Integrate_using_Patterson_adaptive(a, b, 2.0e-4, 1.0e-300, compute_Int_x);
}
//==================================================================================================
//
// computing temperature moments using VIK
//
//==================================================================================================
double dy_k_VIK(double l)
{
    double dtau = Ne_func(l, PSP)*const_sigT;
    double Te = Te_func(l, PSP), The = Te/const_me;
    return dtau*pow(The, PSP.k+1);
}

double integrate_yk_VIK(int k)
{ 
    PSP.k = k;
    double r = compute_Integral(dy_k_VIK, 1.0e-5);
    return PSP.zsc*PSP.Cluster.Get_rc_cm()*r;
}

//==================================================================================================
//
// computing SZ signal by explicitly integrating along the line of sight for each
// frequency bin. This routine is just for comparisons & rather slow.
//
//==================================================================================================
double dSZ_VIK(double l)
{
    Parameters temp = Parameters();
    temp.copyParameters(parameters);
    temp.xcmb[0] = PSP.xcmb;
    temp.Dtau = Ne_func(l, PSP)*const_sigT;
    temp.Te = Te_func(l, PSP);
    temp.betac = PSP.betac;
    temp.muc = PSP.muc;
    temp.betao = temp.muo = 0;
    temp.setCalcValues();

    return compute_signal_combo(0, temp);
}

double integrate_SZ_VIK(double xcmb)
{ 
    PSP.xcmb=xcmb;
    double r = compute_Integral(dSZ_VIK, 1.0e-5);
   
    return PSP.zsc*PSP.Cluster.Get_rc_cm()*r;
}
