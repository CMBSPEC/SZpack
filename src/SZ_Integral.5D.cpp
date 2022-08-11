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
// Author: Jens Chluba (CITA, University of Toronto)
//
// first implementation: May 2012
// last modification: June 2012
//
//==================================================================================================

//==================================================================================================
// Standards
//==================================================================================================
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include <iomanip>
#include <cmath>
#include <limits.h>
#include <vector>

//==================================================================================================
// required libs
//==================================================================================================
#include "physical_consts.h"
#include "routines.h"
#include "Relativistic_MB.h"
#include "nPl_derivatives.h"
#include "Integration_routines.h"
#include "SZ_Integral.5D.h"

//==================================================================================================
//
// namespaces
//
//==================================================================================================
using namespace std;

double Int_eps; // relative accuracy for numerical integration (lower than 10^-6 is hard to achieve)

#define CUT_HIGHER_ORDERS  // Defining this means that higher order corrections in beta are omitted

//==================================================================================================
// 
// cross section times statistical factor
//
//==================================================================================================
double sig_Boltzmann_Compton(double x, double The, double betac, double muc, 
                             double xi, double mue, double mup, double phie, double phip)
{
    double eta=sqrt(2.0*The*xi);

    //==============================================================================================
    // mue== cosine of angle of gamma  and e 
    // mup== cosine of angle of gamma  and gamma'
    //==============================================================================================    
    double beta=eta/sqrt(1.0+eta*eta);
    double gamma2=1.0+eta*eta;
    //
    double muep=mue*mup+cos(phie-phip)*sqrt( (1.0-mue*mue)*(1.0-mup*mup) );
    double mucp=muc*mup+cos(0.0 -phip)*sqrt( (1.0-muc*muc)*(1.0-mup*mup) );
    //
    double kappa =1.0-beta*mue;
    double kappap=1.0-beta*muep;
    //
    double dx=x*beta*(muep-mue)/kappap;
    double xp=x+dx;
    double zeta=1.0/gamma2/kappa/kappap;
    //
    double alpha_sc=1.0-mup;
    //
    double gammac=1.0/sqrt(1.0-betac*betac);
    double xv =gammac*x *(1.0+betac*muc );
    double xvp=gammac*xp*(1.0+betac*mucp);
    //
    double exp_mx=exp(-xv);
    double dex=one_minus_exp_mx(xv, exp_mx); 
    double nx=exp_mx/dex;
    //
    double dexp=one_minus_exp_mx(xvp); 
    double F=nx*(exp(xv-xvp)*dex/dexp-1.0);

#ifdef CUT_HIGHER_ORDERS
    //------------------------------------------------------------
    // cut series at after quadrupole term
    //------------------------------------------------------------
    exp_mx=exp(-x);
    dex=one_minus_exp_mx(x, exp_mx); 
    nx=exp_mx/dex;
    //
    double exp_mxp=exp(-xp);
    dexp=one_minus_exp_mx(xp);
    double nxp=exp_mxp/dexp;
    
    double G=x*nx/dex;
    double Gp=xp*nxp/dexp;
    
    double Q=G*x*(1.0+exp_mx)/dex;
    double Qp=Gp*xp*(1.0+exp_mxp)/dexp;
    double P2 =0.5*(3.0*muc *muc -1.0);
    double P2p=0.5*(3.0*mucp*mucp-1.0);
    
    double M=G*(x*(1.0+exp_mx)/dex-3.0);
    double Mp=Gp*(xp*(1.0+exp_mxp)/dexp-3.0);
    
    F =betac*betac*((Qp*P2p-Q*P2)/3.0+(Mp-M)/6.0);
    F+=betac*(muc*G-mucp*Gp);
    F+=nxp-nx;
    //------------------------------------------------------------
#endif    
    
    return 3.0/8.0/PI*kappa/pow(kappap, 2)/gamma2*( 1.0-zeta*alpha_sc*(1.0-0.5*zeta*alpha_sc) )*F;
}

//==================================================================================================
//
// integrations using Patterson scheme
//
//==================================================================================================
double Boltzmann_I(double phip, void *userdata)
{ 
    double *d=(double *)userdata;
    
    return sig_Boltzmann_Compton(d[0], d[1], d[2], d[3], 
                                 d[4], d[5], d[6], d[7], phip);
}

//==================================================================================================
double phip_Int(double x, double The, double betac, double muc, 
                double xi, double mue, double mup, double phie)
{
    double a=0.0, b=TWOPI;
    double epsrel=Int_eps*0.6, epsabs=1.0e-300;
    
    double d[8]={x, The, betac, muc, xi, mue, mup, phie};
    void *userdata=(void *)d;
    
    return Integrate_using_Patterson_adaptive(a, b, epsrel, epsabs, Boltzmann_I, userdata);
}

//==================================================================================================
double Boltzmann_II(double phie, void *userdata)
{ 
    double *d=(double *)userdata;
    
    return phip_Int(d[0], d[1], d[2], d[3], d[4], d[5], d[6], phie);
}

//==================================================================================================
double phie_Int(double x, double The, double betac, double muc, 
                double xi, double mue, double mup)
{
    double a=0.0, b=TWOPI;
    double epsrel=Int_eps*0.7, epsabs=1.0e-300;
    
    double d[7]={x, The, betac, muc, xi, mue, mup};
    void *userdata=(void *)d;
    
    return Integrate_using_Patterson_adaptive(a, b, epsrel, epsabs, Boltzmann_II, userdata);
}

//==================================================================================================
double Boltzmann_III(double mup, void *userdata)
{ 
    double *d=(double *)userdata;
    
    return phie_Int(d[0], d[1], d[2], d[3], d[4], d[5], mup);
}

//==================================================================================================
double mup_Int(double x, double The, double betac, double muc, double xi, double mue)
{
    double a=-1.0, b=1.0;
    double epsrel=Int_eps*0.8, epsabs=1.0e-300;
    
    double d[6]={x, The, betac, muc, xi, mue};
    void *userdata=(void *)d;
    
    return Integrate_using_Patterson_adaptive(a, b, epsrel, epsabs, Boltzmann_III, userdata);
}

//==================================================================================================
double Boltzmann_IV(double mue, void *userdata)
{ 
    double *d=(double *)userdata;
    
    return mup_Int(d[0], d[1], d[2], d[3], d[4], mue);
}

//==================================================================================================
double mue_Int(double x, double The, double betac, double muc, double xi)
{
    double a=-1.0, b=1.0;
    double epsrel=Int_eps*0.9, epsabs=1.0e-300;
    
    double d[5]={x, The, betac, muc, xi};
    void *userdata=(void *)d;
    
    return Integrate_using_Patterson_adaptive(a, b, epsrel, epsabs, Boltzmann_IV, userdata);
}

//==================================================================================================
double Boltzmann_V(double xi, void *userdata)
{ 
    double *d=(double *)userdata;
    
    return mue_Int(d[0], d[1], d[2], d[3], xi)*f_RM(xi, d[1]);
}

//==================================================================================================
double xi_Int(double x, double The, double betac, double muc)
{
    double a=0.0, lim=30.0, b=lim*(1.0+0.5*lim*The);
    double epsrel=Int_eps, epsabs=1.0e-300;
    
    double d[4]={x, The, betac, muc};
    void *userdata=(void *)d;
    
    return Integrate_using_Patterson_adaptive(a, b, epsrel, epsabs, Boltzmann_V, userdata);
}

//==================================================================================================
//
// 5D integration carried out using Patterson scheme. The SZ signal is obtained in the cluster frame
//
//==================================================================================================
double compute_SZ_distortion_Patterson_5D(double x, 
                                          double The, double betac, double muc, 
                                          double eps_Int)
{ 
    Int_eps=eps_Int;
    
    return norm_f_RM(The)*xi_Int(x, The, betac, muc); 
}

//==================================================================================================
//==================================================================================================
