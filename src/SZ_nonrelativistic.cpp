//==================================================================================================
//
// This program allows computing the traditional non-relativistic thermal SZ and kinematic effect. 
// i.e, for example, as described in the review Mroczkowski et al., 2018.
//
//==================================================================================================
//
// Author: Elizabeth Lee
//
// first implementation: August 2020
//
//==================================================================================================

#include "SZ_nonrelativistic.h"

IntegralNonRelativistic::IntegralNonRelativistic(){
    betac=muc=The=x=exp_mx=dex=0.0;
    run_mode="";
}

IntegralNonRelativistic::IntegralNonRelativistic(double x_i, double The_i, double betac_i, double muc_i){
    betac = betac_i;
    muc = muc_i;
    The = The_i;
    x = x_i;
    exp_mx = exp(-x);
    dex = one_minus_exp_mx(x, exp_mx); 
    run_mode = "";
}

IntegralNonRelativistic::IntegralNonRelativistic(int k, Parameters fp)
    : IntegralNonRelativistic(fp.xcmb[k], fp.calc.The, fp.betac, fp.muc) {}

void IntegralNonRelativistic::Update_x(double x_i){
    x = x_i;
    exp_mx = exp(-x);
    dex = one_minus_exp_mx(x, exp_mx); 
}

double IntegralNonRelativistic::Calculate_monopole(){
    double Gx = x*exp_mx/pow(dex, 2);
    double F = (x*(exp_mx+1.0)/dex-4.0);
    return Gx*F*The;
}

double IntegralNonRelativistic::Calculate_dipole(){
    double Gx = x*exp_mx/pow(dex, 2);
    return betac*muc*Gx;
}

double IntegralNonRelativistic::Calculate_All(){
    return Calculate_dipole()+Calculate_monopole();
}

double IntegralNonRelativistic::compute_distortion(string mode){
    double r = 0.0; 
    run_mode = mode;

    if(run_mode=="monopole"){
        r = Calculate_monopole();
    } 
    else if(run_mode=="dipole"){ 
        r = Calculate_dipole();
    }
    else if(run_mode=="quadrupole"){
        r = 0.0;
    }
    else if(run_mode=="monopole_corr"){ 
        r = 0.0;
    }
    else if(run_mode=="kin"){
        r = Calculate_dipole();
    }
    else {
        r = Calculate_All();
    }
    return r;
}

//==================================================================================================
//
// 3D integration carried out using Patterson scheme
//
// eps_Int : relative accuracy for numerical integration (lower than 10^-6 is hard to achieve)
//
// mode == "monopole"      --> only scattering of monopole without second order kinematic corr
// mode == "dipole"        --> only scattering of dipole     (first order kinematic correction)
// mode == "quadrupole"    --> only scattering of quadrupole (second order kinematic correction)
// mode == "monopole_corr" --> only scattering of second order kinematic correction to monopole
// mode == "all"           --> all terms added
// mode == "kin"           --> only kinematic terms
//
//==================================================================================================

double compute_SZ_distortion_nonrelativistic(double x, double The, double betac, double muc, string mode){
    IntegralNonRelativistic szDistortion = IntegralNonRelativistic(x,The,betac,muc);
    return szDistortion.compute_distortion(mode);
}

void compute_SZ_distortion_nonrelativistic(vector<double> &Dn, vector<double> &x,
                                           double The, double betac, double muc, string mode, bool DI){
    int gridpoints = x.size();
    Dn.resize(gridpoints);
    Parameters fp = Parameters(); //This is just to get a value for the Dn_DI conversion 
    IntegralNonRelativistic szDistortion = IntegralNonRelativistic(x[0],The,betac,muc);
    Dn[0] = (DI ? pow(x[0],3.0)*fp.rare.Dn_DI_conversion() : 1.0)*szDistortion.compute_distortion(mode);
    for(int k = 1; k < gridpoints; k++){
        szDistortion.Update_x(x[k]);
        Dn[k] = szDistortion.compute_distortion(mode);
        if (DI) { Dn[k] *= pow(x[k],3.0)*fp.rare.Dn_DI_conversion(); }
    }
}

double compute_SZ_distortion_nonrelativistic(int k, Parameters fp){
    IntegralNonRelativistic szDistortion = IntegralNonRelativistic(fp.k_inRange(k), fp);
    return szDistortion.compute_distortion(fp.rare.RunMode);
}

void compute_SZ_distortion_nonrelativistic(vector<double> &Dn, Parameters &fp, bool DI){
    Dn.resize(fp.gridpoints);
    IntegralNonRelativistic szDistortion = IntegralNonRelativistic(0, fp);
    Dn[0] = (DI ? pow(fp.xcmb[0],3.0)*fp.rare.Dn_DI_conversion() : 1.0)*fp.Dtau*szDistortion.compute_distortion(fp.rare.RunMode);
    for(int k = 1; k < fp.gridpoints; k++){
        szDistortion.Update_x(fp.xcmb[k]);
        Dn[k] = fp.Dtau*szDistortion.compute_distortion(fp.rare.RunMode);
        if (DI) { Dn[k] *= pow(fp.xcmb[k],3.0)*fp.rare.Dn_DI_conversion(); }
    }
}

//==================================================================================================
//==================================================================================================