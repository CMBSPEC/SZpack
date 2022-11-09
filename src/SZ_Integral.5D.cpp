//==================================================================================================
//
// Computation of the thermal and k-SZ effect using explicit integration of the 5D collision term. 
// The cluster is assumed to be isothermal and moving at a given speed betac with direction muc 
// relative to the line of sight. The integrals are carried out in the cluster frame. Optionally 
// all higher order terms in betac can be included. Use of this routine is only recommended for 
// checks of the results, as it is rather slow.
//
//==================================================================================================
//
// Author: Jens Chluba & Elizabeth Lee
//
// first implementation: May 2012
// last modification: February 2020
//
//==================================================================================================
// 02/2020: Recreated into class form
//
//==================================================================================================

#include "SZ_Integral.5D.h"

#include <cmath>

#include "Parameters.h"
#include "Patterson.h"
#include "Relativistic_MB.h"
#include "nPl_derivatives.h"
#include "Definitions.h"

Integral5D::Integral5D(){
    Int_eps=betac=muc=The=x=mu=mup=xi=phi=phip=dx=xp=dsig=mucp=0.0;
    exp_mx=dex=exp_mxp=dexp=nx=nxp=G=Gp=0.0;
    run_mode="";
    xfac = 1.0;
}

Integral5D::Integral5D(double x_i, double The_i, double betac_i, double muc_i, double eps_Int_i){
    Int_eps = eps_Int_i;
    betac = betac_i;
    muc = muc_i;
    The = The_i;
    x = x_i;
    mu=mup=xi=phi=phip=dx=xp=dsig=mucp=exp_mx=dex=exp_mxp=dexp=nx=nxp=G=Gp=0.0;
    run_mode = "";
    xfac = 1.0;
}

Integral5D::Integral5D(double x_i, Parameters fp)
    : Integral5D(fp.calc.xfac*x_i, fp.calc.The, fp.betac, fp.calc.mucc, fp.relative_accuracy) {
        xfac = fp.calc.xfac;
    }

void Integral5D::Update_x(double x_i){
    x = xfac*x_i;
}

void Integral5D::Calculate_shared_variables(){
    //Standard definitions
    double eta=sqrt(2.0*The*xi);
    //==============================================================================================
    // mu == cosine of angle of gamma  and beta 
    // mup== cosine of angle of gamma' and beta
    //==============================================================================================    
    double beta=eta/sqrt(1.0+eta*eta);
    double gamma2=1.0+eta*eta;
    //
    mucp=muc*mup+cos(0.0 -phip)*sqrt( (1.0-muc*muc)*(1.0-mup*mup) );
    double muep=mu*mup+cos(phi-phip)*sqrt( (1.0-mu*mu)*(1.0-mup*mup) );
    //
    double kappa =1.0-beta*mu;
    double kappap=1.0-beta*muep;
    double zeta=1.0/gamma2/kappa/kappap;
    //
    double alpha_sc=1.0-mup;

    //The variables based on x
    dx=x*beta*(muep-mu)/kappap;
    xp=x+dx;

    //d2sigma/dmu/dmup as defined equation 2 CNSN
    dsig = 3.0/8.0/PI*kappa/pow(kappap, 2)/gamma2*(1.0-zeta*alpha_sc*(1.0-0.5*zeta*alpha_sc));
}

double Integral5D::Calculate_full(){
    double gammac=1.0/sqrt(1.0-betac*betac);
    double xv =gammac*x *(1.0+betac*muc );
    double xvp=gammac*xp*(1.0+betac*mucp);
    //
    exp_mx=exp(-xv);
    dex=one_minus_exp_mx(xv, exp_mx);
    dexp=one_minus_exp_mx(xvp); 
    double nx=exp_mx/dex;
    //
    return nx*(exp(xv-xvp)*dex/dexp-1.0);
}
//The cut higher order version from before
void Integral5D::Calculate_HigherOrder_Variables(){
    exp_mx=exp(-x);
    dex=one_minus_exp_mx(x, exp_mx); 
    exp_mxp=exp(-xp);
    dexp=one_minus_exp_mx(xp, exp_mxp);
    nx=exp_mx/dex;
    nxp=exp_mxp/dexp;
    G=x*nx/dex;
    Gp=xp*nxp/dexp;
}

double Integral5D::Calculate_monopole(){
    return nxp-nx;
}
double Integral5D::Calculate_dipole(){
    return betac*(muc*G-mucp*Gp);
}
double Integral5D::Calculate_quadrupole(){
    double Q=G*x*(1.0+exp_mx)/dex;
    double Qp=Gp*xp*(1.0+exp_mxp)/dexp;
    double P2 =0.5*(3.0*muc *muc -1.0);
    double P2p=0.5*(3.0*mucp*mucp-1.0);
    return betac*betac*(Qp*P2p-Q*P2)/3.0;
}
double Integral5D::Calculate_monopole_correction(){
    double M=G*(x*(1.0+exp_mx)/dex-3.0);
    double Mp=Gp*(xp*(1.0+exp_mxp)/dexp-3.0);
    return betac*betac*(Mp-M)/6.0;
}
double Integral5D::Calculate_kinetic_correction(){
    return Calculate_dipole()+Calculate_quadrupole()+Calculate_monopole_correction();
}

double Integral5D::Calculate_All(){
    return Calculate_kinetic_correction()+Calculate_monopole();
}

double Integral5D::sig_Boltzmann_Compton(double int_phip){
    phip = int_phip;
    Calculate_shared_variables();
    
    double r = 0.0; 

    if(run_mode=="full"){
        r = Calculate_full();
    }
    else {
        Calculate_HigherOrder_Variables();
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
    }
    return dsig*r;

}

double Integral5D::phip_Int(double int_phi){
    phi = int_phi;
    double a=0.0, b=TWOPI;
    double epsrel=Int_eps*0.6, epsabs=1.0e-300;
    
    return Integrate_using_Patterson_adaptive(a, b, epsrel, epsabs, [this](double int_var) { return this->sig_Boltzmann_Compton(int_var);});;
}

double Integral5D::phi_Int(double int_mup){
    mup = int_mup;
    double a=0.0, b=TWOPI;
    double epsrel=Int_eps*0.7, epsabs=1.0e-300;
    
    return Integrate_using_Patterson_adaptive(a, b, epsrel, epsabs, [this](double int_var) { return this->phip_Int(int_var);});
}

double Integral5D::mup_Int(double mu_int){
    mu = mu_int;
    double a=-1.0, b=1.0;
    double epsrel=Int_eps*0.8, epsabs=1.0e-300;
    
    return Integrate_using_Patterson_adaptive(a, b, epsrel, epsabs, [this](double int_var) { return this->phi_Int(int_var);});
}

double Integral5D::mue_Int(double xi_int){
    xi = xi_int;
    double a=-1.0, b=1.0;
    double epsrel=Int_eps*0.9, epsabs=1.0e-300;

    double integral = Integrate_using_Patterson_adaptive(a, b, epsrel, epsabs, [this](double int_var) { return this->mup_Int(int_var);});
    return integral*f_RM(xi, The);
}

double Integral5D::xi_Int(){
    double a=0.0, lim=30.0, b=lim*(1.0+0.5*lim*The);
    double epsrel=Int_eps, epsabs=1.0e-300;
    
    return Integrate_using_Patterson_adaptive(a, b, epsrel, epsabs, [this](double int_var) { return this->mue_Int(int_var);});
}

double Integral5D::compute_distortion(string mode){
    run_mode="full";
    if(mode=="all" || mode=="monopole" || mode=="dipole" || mode=="quadrupole" || mode=="monopole_corr" ||mode=="kin"){
        run_mode = mode;
    }
    return norm_f_RM(The)*xi_Int();
}

//==================================================================================================
//
// 5D integration carried out using Patterson scheme. The SZ signal is obtained in the cluster frame
//
//==================================================================================================
double compute_SZ_distortion_Patterson_5D(double x, double The, double betac, double muc,
                                          double eps_Int, string mode){
    Integral5D szDistortion = Integral5D(x,The,betac,muc,eps_Int);
    return szDistortion.compute_distortion(mode);
}

void compute_SZ_distortion_Patterson_5D(vector<double> &Dn, vector<double> x,
                                        double The, double betac, double muc, double eps_Int,
                                        string mode, bool DI){
    int gridpoints = x.size();
    Dn.resize(gridpoints);
    Parameters fp = Parameters(); //This is just to get a value for the Dn_DI conversion 
    Integral5D szDistortion = Integral5D(x[0],The,betac,muc,eps_Int);
    Dn[0] = (DI ? pow(x[0],3.0)*fp.rare.Dn_DI_conversion() : 1.0)*szDistortion.compute_distortion(mode);
    for(int k = 1; k < gridpoints; k++){
        szDistortion.Update_x(x[k]);
        Dn[k] = szDistortion.compute_distortion(mode);
        if (DI) { Dn[k] *= pow(x[k],3.0)*fp.rare.Dn_DI_conversion(); }
    }
}

double compute_SZ_distortion_Patterson_5D(double x, Parameters fp){
    Integral5D szDistortion = Integral5D(x, fp);
    return szDistortion.compute_distortion(fp.rare.RunMode);
}

void compute_SZ_distortion_Patterson_5D(vector<double> &Dn, Parameters fp, bool DI){
    Dn.resize(fp.gridpoints);
    Integral5D szDistortion = Integral5D(fp.xcmb[0], fp);
    Dn[0] = (DI ? pow(fp.xcmb[0],3.0)*fp.rare.Dn_DI_conversion() : 1.0)*fp.Dtau*szDistortion.compute_distortion(fp.rare.RunMode);
    for(int k = 1; k < fp.gridpoints; k++){
        szDistortion.Update_x(fp.xcmb[k]);
        Dn[k] = fp.Dtau*szDistortion.compute_distortion(fp.rare.RunMode);
        if (DI) { Dn[k] *= pow(fp.xcmb[k],3.0)*fp.rare.Dn_DI_conversion(); }
    }
}

//==================================================================================================
//==================================================================================================
