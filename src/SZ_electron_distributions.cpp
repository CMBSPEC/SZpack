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
    double The = Te/const_me;
    double gamma = sqrt(1.0+eta*eta);
    double exp_th=f_RM_gamma(gamma, The);
    double norm = norm_f_RM_gamma(The);
    return eta*eta*exp_th*norm;
}

//TODO: add non-relativistic Maxwell boltzmann also

//==================================================================================================
// kinematic boost model
//==================================================================================================
double KinematicBoost_Dist(double eta, double Te, double betac, double muc, int l){
    double The = Te/const_me;
    double gammac = 1/sqrt(1-betac*betac);
    double fth = Boltzmann_Dist(eta*gammac,Te);
    double Pl = LegendreP(l, muc);
    double A = gammac*betac*eta/The;
    double norm = sqrt(PI/2/A);
    return gsl_sf_bessel_Inu(l+0.5,A)*Pl*fth*(2.0*l+1.0)*norm*pow(-1,l);
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

    double dbdeta = pow(eta*eta+1.0,-1.5);
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