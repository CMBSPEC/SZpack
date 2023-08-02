//==================================================================================================
//
// Electron distributions for high energy non thermal electron contributions to the SZ signal. These
// are used by SZ_Integral.Kernel. Distributions taken from TODO:
//
//==================================================================================================
//==================================================================================================
//
// Author: Elizabeth Lee
//
// first implementation: September 2020
// last modification   : September 2020
//
//==================================================================================================

#include "SZ_electron_distributions.h"
#include "global_functions.h"

//TODO: Fill in the sources for these
//==================================================================================================
// The electron momentum distribution for thermal electrons
//==================================================================================================
double Boltzmann_Dist(double eta, double Te) {
    double gamma = sqrt(1.0+eta*eta);
    return Boltzmann_Dist_gamma(eta, gamma, Te);
}

double Boltzmann_Dist_gamma(double eta, double gamma, double Te) {
    double The = Te/const_me;
    double exp_th=f_RM_gamma(gamma, The);
    double norm = norm_f_RM_gamma(The);
    return eta*eta*exp_th*norm;
}

//TODO: add non-relativistic Maxwell boltzmann also

//==================================================================================================
// kinematic boost model
//==================================================================================================
// The full form
double FullKinematicBoost_Dist(double eta, double mup, double phip, double Te, double betac, double muc){
    double The = Te/const_me;
    double gammac = 1.0/sqrt(1.0-betac*betac);
    double gamma = sqrt(1.0+eta*eta);
    double fth = Boltzmann_Dist_gamma(eta, gamma*gammac,Te);
    double mug = muc*mup+cos(0.0-phip)*sqrt((1.0-muc*muc)*(1.0-mup*mup));
    double factor = exp(-mug*gammac*betac*eta/The);
    return factor*fth;
}
// The full form using the modified Juttner distribution
double FullKinematicModifiedJuttner(double eta, double mup, double phip, double Te, double betac, double muc){
    double The = Te/const_me;
    double gammac = 1.0/sqrt(1.0-betac*betac);
    double gamma = sqrt(1.0+eta*eta);
    double mug = muc*mup+cos(0.0-phip)*sqrt((1.0-muc*muc)*(1.0-mup*mup));
    double cmbgamma = gammac*(gamma+mug*betac*eta);
    double fth = exp(-cmbgamma/The);
    return eta*eta*fth/cmbgamma;
}

double intfunc_FullKinematicModifiedJuttner(double eta, double Te, double betac, double muc){
    double a=-1.0, b=1.0;
    double epsrel=1.0e-6, epsabs=1.0e-300;
    double norm = Integrate_using_Patterson_adaptive(a, b, epsrel, epsabs, 
        [eta, Te, betac, muc](double int_var) {return FullKinematicModifiedJuttner(eta, int_var, 0.0, Te, betac, muc);});
    return norm;
}

double Norm_FullKinematicModifiedJuttner(double Te, double betac, double muc){
    double a=0.0, lim=30.0, b=lim*(1.0+0.5*lim*0.05);
    double epsrel=1.0e-6, epsabs=1.0e-300;
    double norm = Integrate_using_Patterson_adaptive(a, b, epsrel, epsabs, 
        [Te, betac, muc](double int_var) {return intfunc_FullKinematicModifiedJuttner(int_var, Te, betac, muc);});
    return 2.0/norm;
}

// The multipole expansion
double KinematicBoost_Dist(double eta, double Te, double betac, double muc, int l){
    double The = Te/const_me;
    double gammac = 1/sqrt(1-betac*betac);
    double gamma = sqrt(1.0+eta*eta);
    double fth = Boltzmann_Dist_gamma(eta, gamma*gammac,Te);
    double Pl = LegendreP(l, muc);
    double A = gammac*betac*eta/The;
    double norm = sqrt(PI/2/A);
    return gsl_sf_bessel_Inu(l+0.5,A)*Pl*fth*(2.0*l+1.0)*norm*pow(-1,l);
}

// The betac expansion of the multipole expansion
double KinematicBoost_Dist_exp(double eta, double Te, double betac, double muc, int l, int betac_order){
    if (betac_order > 5 || betac_order<0) {
        print_message("betac_order must be between 0 and 5. Function returning the full kinematic boost dist");
        return KinematicBoost_Dist(eta, Te, betac, muc, l);
    }
    double fth = Boltzmann_Dist(eta,Te);
    double The = Te/const_me;
    double gam0 = sqrt(eta*eta+1);
    if (l==0) {
        double factor = 1.0;
        if (betac_order <2) { return factor*fth; }
        factor += betac*betac*(eta*eta-3.0*gam0*The)/The/The/6.0;
        if (betac_order <4) { return factor*fth; }
        factor += pow(betac,4)*(pow(eta,4)-10.0*The*gam0*eta*eta+5.0*The*The*(7.0*eta*eta+3.0-9.0*gam0*The))/The/The/The/The/120.0;
        return factor*fth;
    }
    if (l==1) {
        if (betac_order == 0) { return 0; }
        double multfactor = -betac*muc*eta/The;
        double factor = 1.0;
        if (betac_order < 3) { return factor*multfactor*fth; }
        factor += betac*betac*(eta*eta+5.0*The*(The-gam0))/The/The/10.0;
        if (betac_order < 5) { return factor*multfactor*fth; }
        factor += pow(betac,4)*(pow(eta,4)-14.0*The*gam0*eta*eta+7.0*The*The*(11.0*eta*eta+5.0*(1.0-5.0*The*gam0+3.0*The*The)))/pow(The,4)/280.0;
        return factor*multfactor*fth;
    }
    if (l==2) {
        if (betac_order < 2) { return 0; }
        double multfactor = betac*betac*(3.0*muc*muc-1.0)*eta*eta/The/The/6.0;
        double factor = 1.0;
        if (betac_order < 4) { return factor*multfactor*fth; }
        factor += betac*betac*(eta*eta-7.0*gam0*The+14.0*The*The)/The/The/14.0;
        return factor*multfactor*fth; 
    }
    if (l==3) {
        if (betac_order < 3) { return 0; }
        double multfactor = -pow(betac,3)*eta*eta*eta*muc*(-3.0+5.0*muc*muc)/The/The/The/30.0;
        double factor = 1.0;
        if (betac_order < 5) { return factor*multfactor*fth; }
        factor += betac*betac*(eta*eta-9.0*gam0*The+27.0*The*The)/The/The/18.0;
        return factor*multfactor*fth;
    }
    print_message("Expansion only implemented for l<=3. Function returning the full kinematic boost dist");
    return KinematicBoost_Dist(eta, Te, betac, muc, l);
}

//==================================================================================================
// Different models for temperature distributions
//==================================================================================================
double CosmicRay_Dist(double eta, double alpha, double p1, double p2) {
    if (eta < p1 || eta > p2){
        return 0;
    }
    return (alpha-1)*pow(eta,-alpha)/(pow(p1,1-alpha)-pow(p2,1-alpha));
}

double ThermalCosmicRay_Dist(double eta, double Te, double alpha, double p1, double p2) {
    if (eta < p1){
        return Boltzmann_Dist(eta, Te);
    }
    if (eta < p2){
        return Boltzmann_Dist(p1, Te)*pow((eta/p1),-alpha);
    }
    return 0;
}

double ThermalCosmicRay_Norm(double Te, double alpha, double p1, double p2) {
    double a=0.0, b=p2;
    double epsrel=1.0e-6, epsabs=1.0e-300;
    double norm = Integrate_using_Patterson_adaptive(a, b, epsrel, epsabs, 
        [Te, alpha, p1, p2](double int_var) {return ThermalCosmicRay_Dist(int_var, Te, alpha, p1, p2);});
    return 1.0/norm;
}

double DoublePower_Dist(double eta, double alpha1, double alpha2, double p1, double pcr, double p2) {
    if (eta < p1 || eta > p2){
        return 0;
    }

    double prefactor = (pow(p1,1-alpha1)-pow(pcr,1-alpha1))/(alpha1-1);
    prefactor += pow(pcr,-alpha1+alpha2)*(pow(pcr,1-alpha2)-pow(p2,1-alpha2))/(alpha2-1);

    if (eta < pcr){
        return pow(eta,-alpha1)/prefactor;
    }
    return pow(pcr,-alpha1+alpha2)*pow(eta,-alpha2)/prefactor;
}

double kappa_Dist(double eta, double Te, double kappa) {
    //kappa an integer >=2
    double The = Te/const_me;

    //double dbdeta = pow(eta*eta+1.0,-1.5);
    double factor = 2.0/sqrt(TWOPI*The)/The;
    double kappa_part = Gamma_JC(kappa+1)/Gamma_JC(kappa-0.5)/pow(kappa-1.5,1.5);
    double E_part = pow(1+eta*eta/2.0/The/(kappa-1.5),-(kappa+1));
    return factor*kappa_part*E_part*eta*eta;
}

double MultiMaxwellian_Dist(double eta, double Te, vector<double> c, vector<double> a) {
    if (a.size() != c.size()){
        exit_error("MultiMaxwellian_Dist: Vectors c and a must be the same length.");
    }

    double The = Te/const_me;

    int np = a.size();
    double x = eta/sqrt(2*The);
    double f = c[0]*exp(-a[0]*x*x);
    for (int loopvar = 1; loopvar < np; loopvar++){
        f += c[loopvar]/pow(a[loopvar]+x*x,2+loopvar/2.0);
    }
    return f*x/2/The*eta;
}

double MultiMaxwellian_Dist(double eta, double Te) {
    vector<double> c{0.755, 0.0609, 2.54, 13.3, 17.58};
    vector<double> a{0.483, 152, 6843, 57.4, 12.4};
    return MultiMaxwellian_Dist(eta, Te, c, a);
}