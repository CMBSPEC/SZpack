//==================================================================================================
//
// Computation of the thermal and k-SZ effect using explicit integration of the collision term in
// different orders of beta. The cluster is assumed to be isothermal and moving at a given speed 
// betac with direction muc relative to the line of sight. The integrals are carried out in the 
// cluster frame. 
//
//==================================================================================================
//
// Author: Jens Chluba  (CITA, University of Toronto)
//
// first implementation: April 2012
// last modification   : July  2012
//
//==================================================================================================
// 8th July: changed definition of S^kin, so that the temperature-independent term is fully canceled

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
#include "SZ_Integral.3D.h"

//==================================================================================================
//
// namespaces
//
//==================================================================================================
using namespace std;


//==================================================================================================
namespace SZIntegral3D
{
    
double Int_eps_3D; // relative accuracy for numerical integration (< 10^-6 is hard to achieve)

string run_mode="all";

double glob_betac, glob_muc;

//==================================================================================================
// 
// cross section for monopole, dipole and quadrupole scattering
//
//==================================================================================================
double sig_Boltzmann_Compton(double x, double The, double xi, double mu, double mug)
{
    double eta=sqrt(2.0*The*xi);

    //==============================================================================================
    // mu == cosine of angle of gamm   and beta 
    // mug== cosine of angle of gamma' and beta
    //==============================================================================================    
    double beta=eta/sqrt(1.0+eta*eta);
    double gamma2=1.0+eta*eta;
    //
    double alpha_sc=1.0-mu*mug;
    double gmumup=0.5*(1.0-mu*mu)*(1.0-mug*mug);
    double beta_sc=alpha_sc*alpha_sc+gmumup;
    //
    double kappa =1.0-beta*mu ;
    double kappap=1.0-beta*mug;
    double dx=x*beta*(mug-mu)/kappap;
    double xp=x+dx;
    double zeta=1.0/gamma2/kappa/kappap;
    //
    double exp_mxp=exp(-xp);
    double dexp=one_minus_exp_mx(xp, exp_mxp); 
    //
    double nfac=3.0/8.0/PI*kappa/pow(kappap, 2)/gamma2;
    //
    double dsig_0_dmudmup=1.0-zeta*(alpha_sc-0.5*beta_sc*zeta);
    
    double exp_mx, dex;
    double F, F0, F1, F2, r;
    
    //==============================================================================================
    // monopole terms only
    //==============================================================================================
    if(run_mode=="monopole")
    {
        exp_mx=exp(-x);
        dex=one_minus_exp_mx(x, exp_mx); 
        F=exp_mx/dex*(exp(-dx)*dex/dexp-1.0);
        
        return nfac*dsig_0_dmudmup*F;
    }
    
    //==============================================================================================
    // dipole terms
    //==============================================================================================
    else if(run_mode=="dipole")
    {
        double G=xp*exp_mxp/pow(dexp, 2);
        F0=dsig_0_dmudmup*mug*mu;
        F1=zeta*(1.0-zeta*alpha_sc)*gmumup;
        
        return -nfac*G*( F0 + F1 );
    }
    
    //==============================================================================================
    // quadrupole terms
    //==============================================================================================
    else if(run_mode=="quadrupole")
    {
        double P2c=0.5*(3.0*mu*mu*mug*mug-1.0);
        double Q=xp*xp*exp_mxp*(exp_mxp+1.0)/pow(dexp, 3);
        F0=dsig_0_dmudmup*P2c;
        F1=zeta*zeta*(P2c-1.0)*gmumup;
        F2=( 1.0+zeta*(2.0+0.75*zeta*gmumup-3.0*alpha_sc*(1.0-zeta)) )*gmumup;
        F=F0 + 2.5*F1 + 1.5*F2;
        
        //JC08.07: take out term for beta = 0
        exp_mx=exp(-x);
        double Qx=x*x*exp_mx*(exp_mx+1.0)/pow(one_minus_exp_mx(x, exp_mx), 3);
        F0=(1.0-(alpha_sc-0.5*beta_sc))*P2c;
        F1=(P2c-1.0)*gmumup;
        F2=(3.0+0.75*gmumup)*gmumup;
        //Qx=0; //JC08.07: old definition of S^kin
        
        return nfac/3.0*( Q*F - Qx*(F0+2.5*F1+1.5*F2) );
    }
    
    //==============================================================================================
    // monopole correction terms
    //==============================================================================================
    else if(run_mode=="monopole_corr")
    {
        exp_mx=exp(-x);
        dex=one_minus_exp_mx(x, exp_mx); 
        r=xp/x;
        F=x*(exp_mx+1.0)/dex*( r*r*exp(-dx)*pow(dex/dexp, 3)*(exp_mxp+1.0)/(exp_mx+1.0)-1.0 );    
        F+=3.0*(1.0-r*exp(-dx)*pow(dex/dexp, 2));    
        F*=x*exp_mx/pow(dex, 2);
        
        return nfac/6.0*dsig_0_dmudmup*F;
    }
    
    //==============================================================================================
    // all terms
    //==============================================================================================
    exp_mx=exp(-x);
    dex=one_minus_exp_mx(x, exp_mx); 
    F=exp_mx/dex*(exp(-dx)*dex/dexp-1.0);
    double r0=dsig_0_dmudmup*F;
    
    double G=x*exp_mx/pow(dex, 2);
    double Gp=xp*exp_mxp/pow(dexp, 2);
    F0=dsig_0_dmudmup*mug*mu;
    F1=zeta*(1.0-zeta*alpha_sc)*gmumup;
    double r1=-Gp*( F0 + F1 ) + G*dsig_0_dmudmup;
    
    double P2c=0.5*(3.0*mu*mu*mug*mug-1.0);
    double Q=x*x*exp_mx*(exp_mx+1.0)/pow(dex, 3);
    double Qp=xp*xp*exp_mxp*(exp_mxp+1.0)/pow(dexp, 3);
    F0=dsig_0_dmudmup*P2c;
    F1=zeta*zeta*(P2c-1.0)*gmumup;
    F2=( 1.0+zeta*(2.0+0.75*zeta*gmumup-3.0*alpha_sc*(1.0-zeta)) )*gmumup;
    double r2=(Qp*( F0 + 2.5*F1 + 1.5*F2)-Q*dsig_0_dmudmup)/3.0;
    
    r=xp/x;
    double DM=x*(exp_mx+1.0)/dex*( r*r*exp(-dx)*pow(dex/dexp, 3)*(exp_mxp+1.0)/(exp_mx+1.0)-1.0 );    
    DM+=3.0*(1.0-r*exp(-dx)*pow(dex/dexp, 2));    
    DM*=x*exp_mx/pow(dex, 2);
    double r0c=dsig_0_dmudmup*DM/6.0;
    
    r =glob_betac*glob_betac*(0.5*(3.0*glob_muc*glob_muc-1.0)*r2+r0c);
    r+=glob_muc*glob_betac*r1;
    if(run_mode=="all") r+=r0; // not added if only kinematic terms wanted
    
    return nfac*r;
}

//==================================================================================================
//
// integrations using Patterson
//
//==================================================================================================
double Boltzmann_I(double mup, void *userdata)
{ 
    double *d=(double *)userdata;
    
    return sig_Boltzmann_Compton(d[0], d[1], d[2], d[3], mup);
}

//==================================================================================================
double mup_Int(double x, double The, double xi, double mu)
{
    double a=-1.0, b=1.0;
    double epsrel=Int_eps_3D*0.8, epsabs=1.0e-300;
    
    double d[4]={x, The, xi, mu};
    void *userdata=(void *)d;
    
    return Integrate_using_Patterson_adaptive(a, b, epsrel, epsabs, Boltzmann_I, userdata);
}

//==================================================================================================
double Boltzmann_II(double mu, void *userdata)
{ 
    double *d=(double *)userdata;
    
    return mup_Int(d[0], d[1], d[2], mu);
}

//==================================================================================================
double mu_Int(double x, double The, double xi)
{
    double a=-1.0, b=1.0;
    double epsrel=Int_eps_3D*0.9, epsabs=1.0e-300;
    
    double d[3]={x, The, xi};
    void *userdata=(void *)d;
    
    return Integrate_using_Patterson_adaptive(a, b, epsrel, epsabs, Boltzmann_II, userdata);
}

//==================================================================================================
double Boltzmann_III(double xi, void *userdata)
{ 
    double *d=(double *)userdata;
    
    return mu_Int(d[0], d[1], xi)*df_RM_dTheta(xi, d[1], (int)d[2]);
}

//==================================================================================================
double xi_Int(double x, double The, int k)
{
    double a=0.0, lim=30.0, b=lim*(1.0+0.5*lim*The);
    double epsrel=Int_eps_3D, epsabs=1.0e-300;
    
    double d[3]={x, The, (double)k};
    void *userdata=(void *)d;
    
    return Integrate_using_Patterson_adaptive(a, b, epsrel, epsabs, Boltzmann_III, userdata);
}

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
using namespace SZIntegral3D;

double compute_SZ_distortion_Patterson_3D(double x, 
                                          double The, double betac, double muc, 
                                          string mode,
                                          double eps_Int)
{
    glob_betac=betac; 
    glob_muc=muc;
    Int_eps_3D=eps_Int;
    run_mode=mode;
    
    double norm_f=norm_df_RM_dTheta(The);
    double r=TWOPI*TWOPI*norm_f*xi_Int(x, The, 0);
    
    if(run_mode=="dipole") 
    {
        r+=x*exp(-x)/pow(one_minus_exp_mx(x), 2);
        r*=betac*muc;
    }
    else if(run_mode=="quadrupole") 
    {
        //r-=x*x*exp(-x)*(exp(-x)+1.0)/pow(one_minus_exp_mx(x), 3)/3.0; //JC08.07: old def of S^kin
        r-=0.3*x*x*exp(-x)*(exp(-x)+1.0)/pow(one_minus_exp_mx(x), 3);
        r*=betac*betac*(3.0*muc*muc-1.0)*0.5;
    }  
    else if(run_mode=="monopole_corr") r*=betac*betac;
    
    return r;
}

//==================================================================================================
//==================================================================================================
