//==================================================================================================
//
// This program allows computing the thermal and kinematic SZ effect using the basis functions of 
// Chluba, Nagai, Sazonov & Nelson, 2012. Up to 20th order temperature corrections are allowed. 
// Kinematic corrections can be included up to second order in the clusters peculiar velocity. 
// For 0.01 < x < 30, 2keV < Te < 75keV and beta < 0.01 the obtained SZ signal should be accurate 
// at the level of 0.001% when 20 temperature order are included. This precision is acheived using 
// three temperature pivots, The=0.01, 0.03, and 0.1. For Te < 2keV the asymptotic expansions 
// should be used.
//
//==================================================================================================
//
// Author: Jens Chluba  & Elizabeth Lee
//
// first implementation: April 2012
// last modification   : February 2020
//
//==================================================================================================
// 18th Apr  2018: improved precision of second derivative with Te
// 12th Sept 2013: fixed bug because of betac^2 P_2(muc) == beta_para^2 - beta_perp^2/2
//  4st  Aug: added derivatives of basis functions in the CMB rest frame
// 22th July: added low and high temperature expansions. Now 2keV < Te < 75keV is covered
// 21th July: added headers with required data; this avoids time taken for reading the files
// 10th July: added basis functions in the CMB rest frame
//  8th July: changed definition of S^kin; temperature-independent terms are fully canceled
// 02/2020: Added Class structure

#include "SZ_CNSN_basis.h"

//==================================================================================================
// these are the reference temperatures that were used to compute basis;
//==================================================================================================
const double Tref_basis_low=0.01; 
const double Tref_basis_mid=0.03; 
const double Tref_basis_high=0.1; 

#include "./src/database/SZ_basis.The_0.01.h"
#include "./src/database/SZ_basis.The_0.03.h"
#include "./src/database/SZ_basis.The_0.1.h" 
// CMB frame data
#include "./src/database/SZ_basis.The_0.01.CMB.h"
#include "./src/database/SZ_basis.The_0.03.CMB.h"
#include "./src/database/SZ_basis.The_0.1.CMB.h" 


//==================================================================================================
CNSNsplineMembers::CNSNsplineMembers(){
    loaded = false;
}

void CNSNsplineMembers::setup_expansion_splines(vector<int> &spline_mem_indices, const double D[][22])
{
    int cols=22;
    int np=1000;

    spline_mem_indices.resize(cols, -1);
    vector<double> xarr(np), yarr(np);

    for(int i=0; i<np; i++) { xarr[i] = D[i][0]; }
    
    for(int c=1; c<cols; c++)
    {
        for(int i=0; i<np; i++) { yarr[i] = D[i][c]; }
        
        spline_mem_indices[c-1]=calc_spline_coeffies_JC(np, &xarr[0], &yarr[0], 
                                                        " setup_expansion_splines ");
    }
}

void CNSNsplineMembers::loadSplines(double refThe, const double dY[][22],const double dD[][22],const double dQ[][22],const double dMcorr[][22]){
    if(!loaded){
        setup_expansion_splines(Y, dY);
        setup_expansion_splines(D, dD);
        setup_expansion_splines(Q, dQ);
        setup_expansion_splines(Mcorr, dMcorr);
        loaded = true;
        The_ref = refThe;
    }
}

//==================================================================================================

namespace CNSNsplines{
    CNSNsplineMembers lowT, midT, highT, lowT_CMB, midT_CMB, highT_CMB;

    //loading expansion splines
    void load_lowT() {lowT.loadSplines(Tref_basis_low, Data_thSZ_low, Data_kSZ_low, Data_qkSZ_low, Data_mkSZ_low);}
    void load_midT() {midT.loadSplines(Tref_basis_mid, Data_thSZ, Data_kSZ, Data_qkSZ, Data_mkSZ);}
    void load_highT() {highT.loadSplines(Tref_basis_high, Data_thSZ_high, Data_kSZ_high, Data_qkSZ_high, Data_mkSZ_high);}
    void load_lowT_CMB() {lowT_CMB.loadSplines(Tref_basis_low, Data_thSZ_I, Data_kSZ_I, Data_qkSZ_I, Data_mkSZ_I);}
    void load_midT_CMB() {midT_CMB.loadSplines(Tref_basis_mid, Data_thSZ_II, Data_kSZ_II, Data_qkSZ_II, Data_mkSZ_II);}
    void load_highT_CMB() {highT_CMB.loadSplines(Tref_basis_high, Data_thSZ_III, Data_kSZ_III, Data_qkSZ_III, Data_mkSZ_III);}
}

//==================================================================================================
//
// compute SZ effect using interpolation of derivative terms
//
//==================================================================================================

IntegralCNSN::IntegralCNSN(double xfac_i, double x_i, double The_i, double betac_i, double muc_i, int Te_order_i, int betac_order_i, bool CMB){
    xfac = xfac_i;
    x = xfac*x_i;
    The = The_i;
    betac = betac_i;
    muc = muc_i;
    CMBframe = CMB;
    Te_order = Te_order_i+1;
    betac_order = betac_order_i;

    // These are only used in the dcompute methods
    beta_perp=betac*sqrt(1.0-muc*muc);
    beta_para=betac*muc;

    x3=pow(x, 3);
    exp_mx=exp(-x);
    dex=one_minus_exp_mx(x, exp_mx);

    setSpline();
}

IntegralCNSN::IntegralCNSN(double x_i, double The_i, double betac_i, double muc_i, int Te_order_i, int betac_order_i, bool CMB)
    : IntegralCNSN(1.0, x_i, The_i, betac_i, muc_i, Te_order_i, betac_order_i, CMB) {}

IntegralCNSN::IntegralCNSN()
    : IntegralCNSN(0.0, 0.0, 0.0, 0.0, 20, 2, false) {}

IntegralCNSN::IntegralCNSN(double x_i, double The_i, double betac_i, double muc_i, bool CMB)
    : IntegralCNSN(x_i, The_i, betac_i, muc_i, 20, 2, CMB) {}


IntegralCNSN::IntegralCNSN(int k, Parameters fp, bool CMB, bool inputOrders)
    : IntegralCNSN((inputOrders ? (CMB ? fp.calc.xfacCMB : fp.calc.xfac) : 1.0), fp.xcmb[k], fp.calc.The, fp.betac, 
                   CMB ? fp.muc : fp.calc.mucc, inputOrders ? fp.T_order : 20, inputOrders ? fp.beta_order : 2, CMB) {}

void IntegralCNSN::Update_x(double x_i){
    x = xfac*x_i;
    x3=pow(x, 3);
    exp_mx=exp(-x);
    dex=one_minus_exp_mx(x, exp_mx);
}

void IntegralCNSN::setSpline(){
    if(The<=0.015){
        if (CMBframe){
            CNSNsplines::load_lowT_CMB();
            spline = CNSNsplines::lowT_CMB;
        }
        else {
            CNSNsplines::load_lowT();
            spline = CNSNsplines::lowT;
        }
    }
    else if(The>=0.045){
        if (CMBframe){
            CNSNsplines::load_highT_CMB();
            spline = CNSNsplines::highT_CMB;
        }
        else {
            CNSNsplines::load_highT();
            spline = CNSNsplines::highT;
        }
    }
    else {
        if (CMBframe){
            CNSNsplines::load_midT_CMB();
            spline = CNSNsplines::midT_CMB;
        }
        else {
            CNSNsplines::load_midT();
            spline = CNSNsplines::midT;
        }
    }
}

// generic function for interpolation (x^3 Dn was used)
double IntegralCNSN::Dn_SZ_splines(double The_ref, vector<int> &spline_mem){
    double lx=log(x);
    double r=calc_spline_JC(lx, spline_mem[0]);

    for(int i=1; i<Te_order; i++) 
        r+=pow(The/The_ref-1.0, i)*calc_spline_JC(lx, spline_mem[i])/factorial(i);
    
    r*=norm_df_RM_dTheta(The)/norm_df_RM_dTheta(The_ref)*The_ref/The;

    return r/x3;    
}

//Needs spline defining first!
double IntegralCNSN::Calculate_monopole(){
    return Dn_SZ_splines(spline.The_ref, spline.Y);
}

double IntegralCNSN::Calculate_dipole(){
    if(betac_order < 1) return 0.0;
    double rnorm=x*exp_mx/dex/dex;

    if(Te_order==0) return betac*muc*rnorm;
    
    double r = Dn_SZ_splines(spline.The_ref, spline.D);
    r+=rnorm;
    return r*betac*muc;
}

double IntegralCNSN::Calculate_quadrupole(){
    if(betac_order < 2) return 0.0;
    double rnorm = CMBframe ? 11.0/30.0 : -3.0/10.0;
    rnorm *= x*exp_mx/dex/dex*x*(1.0+exp_mx)/dex;

    if(Te_order==0) return betac*betac*(1.5*muc*muc-0.5)*rnorm;
    
    double r = Dn_SZ_splines(spline.The_ref, spline.Q);
    r+=rnorm;
    return r*betac*betac*(1.5*muc*muc-0.5);
}

double IntegralCNSN::Calculate_monopole_correction(){
    if(betac_order < 2) return 0.0;
    double rnorm = CMBframe ? x*exp_mx/dex/dex*(x*(1.0+exp_mx)/dex-3.0)/3.0 : 0.0;

    if(Te_order==0) return betac*betac*rnorm;
    
    double r = Dn_SZ_splines(spline.The_ref, spline.Mcorr);
    r+=rnorm;
    return r*betac*betac;
}

double IntegralCNSN::Calculate_kinetic_correction(){
    return Calculate_dipole()+Calculate_quadrupole()+Calculate_monopole_correction();
}

double IntegralCNSN::Calculate_All(){
    return Calculate_kinetic_correction()+Calculate_monopole();
}

// analytic derivatives in betac_parallel and betac_perp
double IntegralCNSN::Dn_for_The(){
    double Dn_th  = Calculate_monopole();
    double Dn_k1  = beta_para*Calculate_dipole();
    double Dn_k2  = (beta_para*beta_para - 0.5*beta_perp*beta_perp)*Calculate_quadrupole(); 
    double Dn_k2m = (beta_para*beta_para + beta_perp*beta_perp)*Calculate_monopole_correction();
    
    return Dn_th+Dn_k1+Dn_k2+Dn_k2m;
}

double IntegralCNSN::Dn_dbeta_para(){
    double Dn_k1  = Calculate_dipole();
    double Dn_k2  = 2.0*beta_para*Calculate_quadrupole();
    double Dn_k2m = 2.0*beta_para*Calculate_monopole_correction();
    
    return Dn_k1+Dn_k2+Dn_k2m;
}

double IntegralCNSN::Dn_d2beta_para(){
    double Dn_k2  = Calculate_quadrupole(); 
    double Dn_k2m = Calculate_monopole_correction(); 
    
    return Dn_k2+Dn_k2m;
}

double IntegralCNSN::Dn_dbeta2_perp(){
    double Dn_k2  = -0.5*Calculate_quadrupole();
    double Dn_k2m = Calculate_monopole_correction(); 
    
    return Dn_k2+Dn_k2m;
}

void IntegralCNSN::compute_all_Te_derivatives_upto_dThe(int dThe, std::function<double()> f){
    double d0, dp, dm, dp2, dm2;
    double eps=0.01;
    
    double The_fid = The;
    Te_order = 21;
    betac_order = 2;
    betac = muc = 1.0;
    
    dDn_dThe.resize(dThe+1);
    d0=f();
    dDn_dThe[0]=d0;
    
    if(dThe>0){
        The = The_fid*(1.0+eps);
        dp = f();
        The = The_fid*(1.0-eps);
        dm = f();
        The = The_fid*(1.0+2*eps);
        dp2= f();
        The = The_fid*(1.0-2*eps);
        dm2= f();
        The = The_fid;
        
        dDn_dThe[1]=(-dp2+8.0*dp-8.0*dm+dm2)/12.0/eps;

        if(dThe>1){
            dDn_dThe[2]=(-dp2+16.0*dp-30.0*d0+16.0*dm-dm2)/12.0/pow(eps, 2) / (2.0);

            if(dThe>2){
                dDn_dThe[3]=(dp2-2.0*dp+2.0*dm-dm2)/2.0/pow(eps, 3) / (6.0);
            
                if(dThe>3){
                    dDn_dThe[4]=(dp2-4.0*dp+6.0*d0-4.0*dm+dm2)/pow(eps, 4) / (24.0);
                }
            }
        }
    }
}

double IntegralCNSN::compute_distortion(string mode){
    run_mode="all";
    if(mode=="monopole" || mode=="dipole" || mode=="quadrupole" || mode=="monopole_corr" ||mode=="kin"){
        run_mode = mode;
    }

    if(run_mode=="monopole"){
        return Calculate_monopole();
    } 
    else if(run_mode=="dipole"){ 
        return Calculate_dipole();
    }
    else if(run_mode=="quadrupole"){
        return Calculate_quadrupole();
    }
    else if(run_mode=="monopole_corr"){ 
        return Calculate_monopole_correction();
    }
    else if(run_mode=="kin"){
        return Calculate_kinetic_correction();
    }
    else{
        return Calculate_All();
    }
}

void IntegralCNSN::Dcompute_distortion(int dThe, int dbeta_para, int dbeta2_perp, vector<double> &dDn){
    dDn_dThe.resize(dThe+1);
    for(int k=0; k<dThe+1; k++) dDn_dThe[k]=0.0;
    
    if(dThe>4) return;
    if(dbeta_para>2) return;
    if(dbeta2_perp>1) return;
    if(dbeta2_perp==1 && dbeta_para!=0) return;
    
    if(dbeta2_perp==1){
        compute_all_Te_derivatives_upto_dThe(dThe, [this]() { return this->Dn_dbeta2_perp();});
    }
    else if(dbeta_para==1){
        compute_all_Te_derivatives_upto_dThe(dThe, [this]() { return this->Dn_dbeta_para();});
    }
    else if(dbeta_para==2){
        compute_all_Te_derivatives_upto_dThe(dThe, [this]() { return this->Dn_d2beta_para();});
    }
    //==============================================================================================
    // The derivatives only
    //==============================================================================================
    else{    
        compute_all_Te_derivatives_upto_dThe(dThe, [this]() { return this->Dn_for_The();});
    }
    dDn = dDn_dThe;
}



//==================================================================================================
//
// compute Dn using improved expansion in cluster frame (Chluba, Nagai, Sazonov, Nelson, 2012)
//
// mode == "monopole"      --> only scattering of monopole without second order kinematic corr
// mode == "dipole"        --> only scattering of dipole     (first order kinematic correction)
// mode == "quadrupole"    --> only scattering of quadrupole (second order kinematic correction)
// mode == "monopole_corr" --> only scattering of second order kinematic correction to monopole
// mode == "all"           --> all terms added
// mode == "kin"           --> only kinematic terms
//
//==================================================================================================
double compute_SZ_distortion_CNSN_basis(double x, 
                                        double The, double betac, double muc, 
                                        int Te_order, int betac_order, 
                                        string mode, bool CMBframe){
    IntegralCNSN szDistortion = IntegralCNSN(x,The,betac,muc,Te_order,betac_order, CMBframe);
    return szDistortion.compute_distortion(mode);
}

void compute_SZ_distortion_CNSN_basis(vector<double> &Dn, vector<double> x,
                                      double The, double betac, double muc, int Te_order, int betac_order,
                                      string mode, bool DI, bool CMBframe){
    int gridpoints = x.size();
    Dn.resize(gridpoints);
    Parameters fp = Parameters(); //This is just to get a value for the Dn_DI conversion 
    IntegralCNSN szDistortion = IntegralCNSN(x[0],The,betac,muc,Te_order,betac_order,CMBframe);
    Dn[0] = (DI ? pow(x[0],3.0)*fp.rare.Dn_DI_conversion() : 1.0)*szDistortion.compute_distortion(mode);
    for(int k = 1; k < gridpoints; k++){
        szDistortion.Update_x(x[k]);
        Dn[k] = szDistortion.compute_distortion(mode);
        if (DI) { Dn[k] *= pow(x[k],3.0)*fp.rare.Dn_DI_conversion(); }
    }
}

double compute_SZ_distortion_CNSN_basis(int k, Parameters fp, bool CMBframe){
    IntegralCNSN szDistortion = IntegralCNSN(fp.k_inRange(k), fp, CMBframe);
    return szDistortion.compute_distortion(fp.rare.RunMode);
}

void compute_SZ_distortion_CNSN_basis(vector<double> &Dn, Parameters fp, bool DI, bool CMBframe){
    Dn.resize(fp.gridpoints);
    IntegralCNSN szDistortion = IntegralCNSN(0, fp, CMBframe);
    Dn[0] = (DI ? pow(fp.xcmb[0],3.0)*fp.rare.Dn_DI_conversion() : 1.0)*fp.Dtau*szDistortion.compute_distortion(fp.rare.RunMode);
    for(int k = 1; k < fp.gridpoints; k++){
        szDistortion.Update_x(fp.xcmb[k]);
        Dn[k] = fp.Dtau*szDistortion.compute_distortion(fp.rare.RunMode);
        if (DI) { Dn[k] *= pow(fp.xcmb[k],3.0)*fp.rare.Dn_DI_conversion(); }
    }
}

//==================================================================================================
// Derivatives (The^k d^k_dThe /k!) (betapara^m d^m_dbetapara /m!) (beta2perp^l d^l_beta2perp /l!) S
// in the CMB frame for a resting observer. Maximal orders in The and betac are used to compute 
// the derivatives.
//
// constraints: dThe<=4; dbeta_para<=2; dbeta2_perp<=1;
//==================================================================================================
void Dcompute_SZ_distortion_CNSN_basis(double x, 
                                           int dThe, int dbeta_para, int dbeta2_perp,
                                           double The, double betac, double muc,
                                           vector<double> &dDn_dThe, bool CMBframe){
    IntegralCNSN szDistortion = IntegralCNSN(x,The,betac,muc,CMBframe);
    szDistortion.Dcompute_distortion(dThe, dbeta_para,dbeta2_perp, dDn_dThe);
}

void Dcompute_SZ_distortion_CNSN_basis(double x, Parameters &fp, bool CMBframe){
    IntegralCNSN szDistortion = IntegralCNSN(0, fp, CMBframe, false);
    szDistortion.Update_x(x);
    szDistortion.Dcompute_distortion(fp.D.dThe, fp.D.dbeta_para, fp.D.dbeta2_perp, fp.D.dDn_dThe);
}
//==================================================================================================
//==================================================================================================