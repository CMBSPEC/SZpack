//==================================================================================================
//
// Computation of the thermal and k-SZ effect using explicit integration of the collision term in
// different orders of beta. The cluster is assumed to be isothermal and moving at a given speed 
// betac with direction muc relative to the line of sight. The integrals are carried out in the 
// cluster frame. 
//
//==================================================================================================
//
// Author: Jens Chluba & Elizabeth Lee
//
// first implementation: April 2012
// last modification   : February 2020
//
//==================================================================================================
// 8th July: changed definition of S^kin, so that the temperature-independent term is fully canceled
// 02/2020: Recreated into class form

#include "SZ_Integral.3D.h"

Integral3D::Integral3D(){
    Int_eps=betac=muc=The=x=mu=mup=xi=dx=xp=dsig0=dsig1=dsig2=nfac=exp_mx=dex=exp_mxp=dexp=0.0;
    run_mode="";
    xfac=1.0;
}

Integral3D::Integral3D(double x_i, double The_i, double betac_i, double muc_i, double eps_Int_i){
    Int_eps = eps_Int_i;
    betac = betac_i;
    muc = muc_i;
    The = The_i;
    x = x_i;
    mu=mup=xi=dx=xp=dsig0=dsig1=dsig2=nfac=exp_mx=dex=exp_mxp=dexp=0.0;
    run_mode = "";
    xfac = 1.0;
}

Integral3D::Integral3D(double x_i, Parameters fp)
    : Integral3D(fp.calc.xfac*x_i, fp.calc.The, fp.betac, fp.calc.mucc, fp.relative_accuracy) {
        xfac = fp.calc.xfac;
    }

void Integral3D::Update_x(double x_i){
    x = xfac*x_i;
}

void Integral3D::Calculate_shared_variables(){
    //Standard definitions
    double eta=sqrt(2.0*The*xi);
    //==============================================================================================
    // mu == cosine of angle of gamma  and beta 
    // mup== cosine of angle of gamma' and beta
    //==============================================================================================    
    double beta=eta/sqrt(1.0+eta*eta);
    double gamma2=1.0+eta*eta;
    //
    double alpha_sc=1.0-mu*mup;
    double gmumup=0.5*(1.0-mu*mu)*(1.0-mup*mup);
    double beta_sc=alpha_sc*alpha_sc+gmumup;
    //
    double kappa =1.0-beta*mu ;
    double kappap=1.0-beta*mup;

    double zeta=1.0/gamma2/kappa/kappap;

    //The variables based on x
    dx=x*beta*(mup-mu)/kappap;
    xp=x+dx;

    exp_mx=exp(-x);
    dex=one_minus_exp_mx(x, exp_mx); 
    exp_mxp=exp(-xp);
    dexp=one_minus_exp_mx(xp, exp_mxp); 

    //d2sigma/dmu/dmup for the 3 multipoles as in Appendix B CNSN
    nfac = 3.0/8.0/PI*kappa/pow(kappap, 2)/gamma2;

    dsig0 = 1.0-zeta*(alpha_sc-0.5*beta_sc*zeta);

    dsig1 = dsig0*mup*mu;
    dsig1 += zeta*(1.0-zeta*alpha_sc)*gmumup;

    double P2c=0.5*(3.0*mu*mu*mup*mup-1.0);
    double F0=dsig0*P2c;
    double F1=zeta*zeta*(P2c-1.0)*gmumup;
    double F2=( 1.0+zeta*(2.0+0.75*zeta*gmumup-3.0*alpha_sc*(1.0-zeta)) )*gmumup;
    dsig2 = F0 + 2.5*F1 + 1.5*F2;
}

double Integral3D::Calculate_monopole(){
    double F = exp_mx/dex*(exp(-dx)*dex/dexp-1.0);
    return dsig0*F;
}

double Integral3D::Calculate_dipole(){
    double Gx = x*exp_mx/pow(dex, 2);
    double Gp = xp*exp_mxp/pow(dexp, 2);
    
    return betac*muc*(-Gp*dsig1 + Gx*dsig0);
}

double Integral3D::Calculate_quadrupole(){
    double Qp = xp*xp*exp_mxp*(exp_mxp+1.0)/pow(dexp, 3);
    double Qx = x*x*exp_mx*(exp_mx+1.0)/pow(dex, 3);
    double P2m = 0.5*(3.0*muc*muc-1.0);

    return betac*betac*P2m*(Qp*dsig2-Qx*dsig0)/3.0;
}

double Integral3D::Calculate_monopole_correction(){
    double r=xp/x;
    double F=x*(exp_mx+1.0)/dex*( r*r*exp(-dx)*pow(dex/dexp, 3)*(exp_mxp+1.0)/(exp_mx+1.0)-1.0 );    
    F+=3.0*(1.0-r*exp(-dx)*pow(dex/dexp, 2));    
    F*=x*exp_mx/pow(dex, 2); //This is equivalent to (Qp-Qx)-3(Gp-Gx)
    
    return betac*betac*dsig0*F/6.0;
}

double Integral3D::Calculate_kinetic_correction(){
    return Calculate_dipole()+Calculate_quadrupole()+Calculate_monopole_correction();
}

double Integral3D::Calculate_All(){
    return Calculate_kinetic_correction()+Calculate_monopole();
}

double Integral3D::sig_Boltzmann_Compton(double int_mup){
    mup = int_mup;
    Calculate_shared_variables();
    
    double r = 0.0; 

    if(run_mode=="monopole"){
        r = Calculate_monopole();
    } 
    else if(run_mode=="dipole"){ 
        r = Calculate_dipole();
    }
    else if(run_mode=="quadrupole"){
        r = Calculate_quadrupole();
    }
    else if(run_mode=="monopole_corr"){ 
        r = Calculate_monopole_correction();
    }
    else if(run_mode=="kin"){
        r = Calculate_kinetic_correction();
    }
    else if(run_mode=="all"){
        r = Calculate_All();
    }
    return nfac*r;
}

double Integral3D::mup_Int(double int_mu){
    mu = int_mu;

    double a=-1.0, b=1.0;
    double epsrel=Int_eps*0.8, epsabs=1.0e-300;
    
    return Integrate_using_Patterson_adaptive(a, b, epsrel, epsabs, [this](double int_var) { return this->sig_Boltzmann_Compton(int_var);});
}

double Integral3D::mue_Int(double int_xi){
    xi = int_xi;

    double a=-1.0, b=1.0;
    double epsrel=Int_eps*0.9, epsabs=1.0e-300;

    double integral = Integrate_using_Patterson_adaptive(a, b, epsrel, epsabs, [this](double int_var) { return this->mup_Int(int_var);});
    return integral*df_RM_dTheta(xi, The, 0);
}

double Integral3D::xi_Int(){
    double a=0.0, lim=30.0, b=lim*(1.0+0.5*lim*The);
    double epsrel=Int_eps, epsabs=1.0e-300;
    
    return Integrate_using_Patterson_adaptive(a, b, epsrel, epsabs, [this](double int_var) { return this->mue_Int(int_var);});
}

double Integral3D::compute_distortion(string mode){
    run_mode="all";
    if(mode=="monopole" || mode=="dipole" || mode=="quadrupole" || mode=="monopole_corr" ||mode=="kin"){
        run_mode = mode;
    }
    double norm_f=norm_df_RM_dTheta(The);
    return TWOPI*TWOPI*norm_f*xi_Int();
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

double compute_SZ_distortion_Patterson_3D(double x, 
                                          double The, double betac, double muc,
                                          double eps_Int, string mode){
    Integral3D szDistortion = Integral3D(x,The,betac,muc,eps_Int);
    return szDistortion.compute_distortion(mode);
}

void compute_SZ_distortion_Patterson_3D(vector<double> &Dn, vector<double> x,
                                        double The, double betac, double muc, double eps_Int,
                                        string mode, bool DI){
    int gridpoints = x.size();
    Dn.resize(gridpoints);
    Parameters fp = Parameters(); //This is just to get a value for the Dn_DI conversion 
    Integral3D szDistortion = Integral3D(x[0],The,betac,muc,eps_Int);
    Dn[0] = (DI ? pow(x[0],3.0)*fp.rare.Dn_DI_conversion() : 1.0)*szDistortion.compute_distortion(mode);
    for(int k = 1; k < gridpoints; k++){
        szDistortion.Update_x(x[k]);
        Dn[k] = szDistortion.compute_distortion(mode);
        if (DI) { Dn[k] *= pow(x[k],3.0)*fp.rare.Dn_DI_conversion(); }
    }
}

double compute_SZ_distortion_Patterson_3D(double x, Parameters fp){
    Integral3D szDistortion = Integral3D(x, fp);
    return szDistortion.compute_distortion(fp.rare.RunMode);
}

void compute_SZ_distortion_Patterson_3D(vector<double> &Dn, Parameters &fp, bool DI){
    Dn.resize(fp.gridpoints);
    Integral3D szDistortion = Integral3D(fp.xcmb[0], fp);
    Dn[0] = (DI ? pow(fp.xcmb[0],3.0)*fp.rare.Dn_DI_conversion() : 1.0)*fp.Dtau*szDistortion.compute_distortion(fp.rare.RunMode);
    for(int k = 1; k < fp.gridpoints; k++){
        szDistortion.Update_x(fp.xcmb[k]);
        Dn[k] = fp.Dtau*szDistortion.compute_distortion(fp.rare.RunMode);
        if (DI) { Dn[k] *= pow(fp.xcmb[k],3.0)*fp.rare.Dn_DI_conversion(); }
    }
}

//==================================================================================================
//==================================================================================================
