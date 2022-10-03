//==================================================================================================
//
// run_SZ_moment_method routine
// 
//==================================================================================================
//
// purpose: give some explicit examples of how to use the SZpack temperature-velocity moment method
// to calculate the SZ signals for simple cluster profiles. The code can be run with default 
// parameters that are set in the function
//
//      int main(int narg, char *args[])
//
// below. Alternatively one can call it with a start file.
//
//==================================================================================================
//
// Author: Jens Chluba (CITA, University of Toronto)
//
// first implementation: July 2012
// last modification   : July 2012
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
#include "Integration_routines.h"
#include "routines.h"
#include "nPl_derivatives.h"
#include "SZpack.h"

#include "gsl/gsl_multimin.h"
#include <gsl/gsl_linalg.h>

//==================================================================================================
//
// namespaces
//
//==================================================================================================
using namespace std;

bool show_mess=1;

//==================================================================================================
//
// print message to screen
//
//==================================================================================================
void print_message(string mess)
{
    if(show_mess)
    {
        cout << "\n " << setfill('=') << setw(90) << "=" << endl;  
        cout << " || " << mess << endl;
        cout << " " << setfill('=') << setw(90) << "=" << endl << endl;   
    }
    
    return;
}

//==================================================================================================
double LegendreP(int l, double x)
{
    vector<double> P(l+2);
    
    P[0]=1;
    P[1]=x;
    
    for(int k=2; k<=l; k++) P[k]=((2.0*k-1)*x*P[k-1]-(k-1.0)*P[k-2])/k;
    
    return P[l];
}

double LegendreP_m1(int l, double x)
{
    vector<double> P(l+2);
    
    P[0]=0;
    P[1]=-sqrt(1.0-x*x);
    
    for(int k=2; k<=l; k++) P[k]=((2.0*k-1)*x*P[k-1]-k*P[k-2])/(k-1);
    
    return P[l];
}

//==================================================================================================
//
// const density sphere
//
//==================================================================================================
double tau0sphere(double z0, double mur)
{
    return sqrt(1.0-z0*z0*(1.0-mur*mur))-z0*mur;
}

double tau0sphere(double b)
{
    return 2.0*sqrt(1.0-b*b);
}

//==================================================================================================
double dtau_dmur_sphere(double mur, void *p)
{
    double *d=(double *)p;
    double z0=d[0];
    int l=d[1];
    
    return LegendreP(l, mur)*( tau0sphere(z0, mur)+pow(-1.0, l)*tau0sphere(z0, -mur) );
}

double taul_hat_Int_sphere(double z0, int l)
{
    double a=0.0, b=1.0;
    double epsrel=1.0e-5, epsabs=1.0e-16;
    
    double d[2]={z0, l};
    void *userdata=(void *)d;
    
    return Integrate_using_Patterson_adaptive(a, b, epsrel, epsabs, dtau_dmur_sphere, userdata)/2.0;
}

//==================================================================================================
double dtau_ds_sphere(double s, void *p)
{
    double *d=(double *)p;
    double b=d[0];
    int l=d[1];
    
    double z0eff=sqrt(b*b+s*s)+1.0e-300;
    double mureff=s/z0eff;
    
    return LegendreP(l, mureff)*taul_hat_Int_sphere(z0eff, l);
}

double taul_av_Int_sphere(double b, int l)
{
    if(l>0 && l%2) return 0;
    
    double av=0, bv=sqrt(1.0-b*b);
    double epsrel=1.0e-5, epsabs=1.0e-16;
    
    double d[2]={b, l};
    void *userdata=(void *)d;
    
    return 2.0*(2.0*l+1.0)*Integrate_using_Patterson_adaptive(av, bv, epsrel, epsabs,
                                                              dtau_ds_sphere, userdata);
}

//==================================================================================================
//
// plot <tau_l>
//
//==================================================================================================
void output_plot_tau_l_sphere(string filename, int lmax, int np)
{
    vector<double> xc(np);
    init_xarr(0.001, 1.0-1.0e-4, &xc[0], np, 0, 0);
    
    ofstream ofile(filename.c_str());
    ofile.precision(10);
    
    double r=0.0;
    
    vector<double> taulc(lmax+2);
    vector<double> taulv(lmax+2);
    
    for(int k=0; k<=lmax; k+=1)
    {
        taulc[k]=taul_av_Int_sphere(0.0, k)/(0.5*pow(tau0sphere(0.0), 2));
        r+=taulc[k];
        
        cout << k << " " << taulc[k] << " " << r << endl;
    }
    
    for(int m=0; m<np; m++)
    {
        cout << scientific << xc[m] << endl;
        ofile << scientific << xc[m] << " ";
        
        double tauscale=0.5*pow(tau0sphere(0), 2);
        
        r=0;
        for(int k=0; k<=lmax; k+=2)
        {
            taulv[k]=taul_av_Int_sphere(xc[m], k);
            
            ofile << scientific
                  << taulv[k]/tauscale/taulc[k]
                  << " ";
            
            r+=taulv[k]/(0.5*pow(tau0sphere(xc[m]), 2));
        }
        
        ofile << r << " "; // check sum == 1
        
        double tauscale2=tau0sphere(0)*tau0sphere(xc[m]);
        double kappa_E_iso=0.5*pow(tau0sphere(xc[m]), 2);
        double kappa_T=taulv[0]+taulv[2]*0.1-kappa_E_iso;
        double kappa_S=-0.6*taulv[2];
        double kappa_E=taulv[0]+taulv[2]*0.1;
        
        ofile << kappa_E_iso/tauscale << " " << kappa_T/tauscale << " "
              <<  kappa_S/tauscale << " " << kappa_E/tauscale << " ";

        ofile << kappa_E_iso/tauscale2 << " " << kappa_T/tauscale2 << " "
              <<  kappa_S/tauscale2 << " " << kappa_E/tauscale2 << " ";

        ofile << endl;
    }
    
    ofile.close();
    
    return;
}

//==================================================================================================
//==================================================================================================


//==================================================================================================
//
// isothermal beta profile
//
//==================================================================================================
double Ne_beta(double r, double beta)
{
    return pow(1.0+r*r, -1.5*beta);
}

double tau0(double beta, double b)
{
    return sqrt(PI)*Gamma_JC(1.5*beta-0.5)/Gamma_JC(1.5*beta)*pow(1.0+b*b, 0.5-1.5*beta);
}

//==================================================================================================
double dtau_ds(double ls, void *p)
{
    double s=exp(ls);
    double *d=(double *)p;
    double z0=d[0], mur=d[1], beta=d[2];
    
    double Ne=Ne_beta(sqrt(s*s+z0*z0+2.0*s*z0*mur), beta);
    
    return s*Ne;
}

double s_Int(double z0, double mur, double beta)
{
    double a=log(1.0e-6), b, r=0.0;
    double epsrel=1.0e-4;
    
    double d[3]={z0, mur, beta};
    void *userdata=(void *)d;
    
    for(int k=0; k<=20; k++)
    {
        b=log(10.0)+a;
        double r1=Integrate_using_Patterson_adaptive(a, b, epsrel, epsrel*(1.0e-100+r)*0.01, dtau_ds, userdata);
        r+=r1;
        
        a=b;
        
        if(fabs(r1/r)<epsrel*0.01 && a>10.0) break;
    }
    
    return r;
}

//==================================================================================================
double dtau_dsdmu(double mur, void *p)
{
    double *d=(double *)p;
    double Pl=LegendreP((int)d[2], mur);
    return Pl*( s_Int(d[0], mur, d[1])+pow(-1.0, (int)d[2])*s_Int(d[0], -mur, d[1]) );
}

double mur_Int_explicit(double z0, double beta, int l)
{
    double a=0, b=1.0;
    double epsrel=1.0e-4, epsabs=1.0e-10;
    
    double d[3]={z0, beta, l};
    void *userdata=(void *)d;

    return Integrate_using_Patterson_adaptive(a, b, epsrel, epsabs, dtau_dsdmu, userdata)/2.0;
}

//==================================================================================================
vector<bool> tau_spline_set(100, 0);
vector<int> memindex_tau_spline(100, 0);

double mur_Int(double z0, double beta, int l)
{
    if(!tau_spline_set[l])
    {
        int np=1000;
        double a=log(1.0e-6), b=log(1.0e+4);

        cout << " setting up sline for l = " << l << " (# points = " << np << ")" << endl;
        
        vector<double> xc(np), yc(np);
        init_xarr(a, b, &xc[0], np, 0, 0);
        
        for(int k=0; k<np; k++) yc[k]=mur_Int_explicit(exp(xc[k]), beta, l);
        
        memindex_tau_spline[l]=calc_spline_coeffies_JC(np, &xc[0], &yc[0]);
        tau_spline_set[l]=1;
    }
    
    z0=min(max(2.0e-6, z0), 9.9e+3);

    return calc_spline_JC(log(z0), memindex_tau_spline[l]);
}

//==================================================================================================
double dtau_dz0(double lz0, void *p)
{
    double *d=(double *)p;
    double z0=exp(lz0), b=d[0], beta=d[1];
    int l=(int)d[2];

    double z0eff=sqrt(b*b+z0*z0)+1.0e-300;
    double mureff=z0/z0eff;
 
    double Ne=Ne_beta(z0eff, beta);
    double Pl=LegendreP(l,  mureff);
    
    return z0*Ne*Pl*mur_Int(z0eff, beta, l);
}

double taul_av_Int_beta(double b0, double beta, int l)
{
    if(l>0 && l%2) return 0;

    double epsrel=2.0e-4, epsabs=1.0e-10;
    
    double d[3]={b0, beta, l};
    void *userdata=(void *)d;
    
    double a=log(1.0e-6), b=log(1.0e+4);
    return 2.0*(2.0*l+1.0)*Integrate_using_Patterson_adaptive(a, b, epsrel, epsabs, dtau_dz0, userdata);
}

double z0_Int(double b0, double beta, int l)
{
    return taul_av_Int_beta(b0, beta, l)/(2.0*l+1.0);
}

//==================================================================================================
double dtau_dz0_m1(double lz0, void *p)
{
    double *d=(double *)p;
    double z0=exp(lz0), b=d[0], beta=d[1];
    int l=(int)d[2];

    double z0eff=sqrt(b*b+z0*z0)+1.0e-300;
    double mureff=z0/z0eff;
    
    double Ne=Ne_beta(z0eff, beta);
    double Plp=LegendreP_m1(l,  mureff);
    
    return z0*Ne*Plp*mur_Int(z0eff, beta, l);
}

double z0_Int_m1(double b0, double beta, int l)
{
    if(l>0 && (l+1)%2) return 0;

    double epsrel=2.0e-4, epsabs=1.0e-10;
    
    double d[3]={b0, beta, l};
    void *userdata=(void *)d;
    
    double a=log(1.0e-6), b=log(1.0e+4);
    return 2.0*Integrate_using_Patterson_adaptive(a, b, epsrel, epsabs, dtau_dz0_m1, userdata);
}

//==================================================================================================
double lambda_function(double b0, double beta, int l)
{
    double dum1=0.0;

    if(l>=1) dum1=l*z0_Int(b0, beta, l-1);
    
    dum1+=(l+1.0)*z0_Int(b0, beta, l+1);
    
    return dum1;
}

double kappa_function(double b0, double beta, int l)
{
    double dum1=0.0;
    
    if(l>=1) dum1=z0_Int_m1(b0, beta, l-1);
    
    dum1-=z0_Int_m1(b0, beta, l+1);
    
    return dum1;
}

//==================================================================================================
//
// plot <tau_l>
//
//==================================================================================================
void output_plot_tau_l(string filename, int lmax, int np)
{
    vector<double> xc(np);
    init_xarr(1.0e-1, 1.0e+2, &xc[0], np, 1, 0);
    
    ofstream ofile(filename.c_str());
    ofile.precision(10);
    
    double beta=2.0/3.0, r=0.0;
    
    vector<double> taulc(lmax+2);
    vector<double> taulv(lmax+2);
    
    for(int k=0; k<=lmax; k+=2)
    {
        taulc[k]=taul_av_Int_beta(0.0, beta, k)/(0.5*pow(tau0(beta, 0), 2));
        r+=taulc[k];
        
        cout << k << " " << taulc[k] << " " << r << endl;
    }
    
    for(int m=0; m<np; m++)
    {
        cout << scientific << xc[m] << endl;
        ofile << scientific << xc[m] << " ";
        
        double tauscale=0.5*pow(tau0(beta, 0), 2);
        
        r=0;
        for(int k=0; k<=lmax; k+=2)
        {
            taulv[k]=taul_av_Int_beta(xc[m], beta, k);
            
            ofile << scientific
                  << taulv[k]/tauscale/taulc[k]
                  << " ";
            
            r+=taulv[k]/(0.5*pow(tau0(beta, xc[m]), 2));
        }
        
        ofile << r << " ";
        
        double tauscale2=tau0(beta, 0)*tau0(beta, xc[m]);
        double kappa_E_iso=0.5*pow(tau0(beta, xc[m]), 2);
        double kappa_T=taulv[0]+taulv[2]*0.1-kappa_E_iso;
        double kappa_S=-0.6*taulv[2];
        double kappa_E=taulv[0]+taulv[2]*0.1;
        
        ofile << tau0(beta, xc[m])/tau0(beta, 0) << " ";
        ofile << kappa_E_iso/tauscale << " " << kappa_T/tauscale << " "
              << kappa_S/tauscale << " " << kappa_E/tauscale << " ";
        
        ofile << kappa_E_iso/tauscale2 << " " << kappa_T/tauscale2 << " "
              << kappa_S/tauscale2 << " " << kappa_E/tauscale2 << " ";

        ofile << endl;
    }
    
    ofile.close();
    
    return;
}

//==================================================================================================
//
// CMB anisotropies signal
//
//==================================================================================================
void output_CMB_isotropy_signals(string filename, int np, double Te)
{
    vector<double> xc(np);
    init_xarr(1.0e-1, 40.0, &xc[0], np, 1, 0);
    
    ofstream ofile(filename.c_str());
    ofile.precision(8);
    
    for(int m=0; m<np; m++)
    {
        cout << xc[m] << endl;
        
        ofile << xc[m] << " ";
        
       for(int l=0; l<=5; l++)
            ofile << pow(xc[m], 3)*compute_SZ_distortion_Patterson_multiple(xc[m], l, Te/const_me, "CMB")
                         /(Te/const_me) << " ";
        
        ofile << endl;
    }
    
    ofile.close();
   
    return;
}

//==================================================================================================
//
// plot first few correction functions
//
//==================================================================================================
void output_CMB_isotropy_Y_functions(string filename, int np)
{
    vector<double> xc(np);
    init_xarr(1.0e-1, 40.0, &xc[0], np, 1, 0);
    
    ofstream ofile(filename.c_str());
    ofile.precision(8);
    vector<double> Y(3);
    
    for(int m=0; m<np; m++)
    {
        cout << xc[m] << endl;
        
        ofile << xc[m] << " ";
        
        compute_Y(xc[m], Y);
        for(int l=0; l<(int)Y.size(); l++) ofile << pow(xc[m], 3)*Y[l] << " ";
        
        compute_Yl0_k(xc[m], 0, Y);
        for(int l=0; l<(int)Y.size(); l++) ofile << pow(xc[m], 3)*Y[l] << " ";

        compute_Yl0_k(xc[m], 1, Y);
        for(int l=0; l<(int)Y.size(); l++) ofile << pow(xc[m], 3)*Y[l] << " ";

        compute_Yl0_k(xc[m], 2, Y);
        for(int l=0; l<(int)Y.size(); l++) ofile << pow(xc[m], 3)*Y[l] << " ";

        compute_Yl0_k(xc[m], 3, Y);
        for(int l=0; l<(int)Y.size(); l++) ofile << pow(xc[m], 3)*Y[l] << " ";

        ofile << endl;
    }
    
    ofile.close();
    
    return;
}

void compute_Y(double x, vector<double> &Y);

//==================================================================================================
//
// plot lowest order signal
//
//==================================================================================================
void output_lowest_order_signal(string filename, int np, double b)
{
    vector<double> xc(np);
    init_xarr(1.0e-1, 40.0, &xc[0], np, 1, 0);
    
    ofstream ofile(filename.c_str());
    ofile.precision(8);
    vector<double> Y(3);
    
    double The=0.01, tau=0.01;
    double beta=2.0/3.0;

    double ISA=tau*0.5;
    double scale=500;
    
    for(int m=0; m<np; m++)
    {
        cout << xc[m] << endl;
        
        ofile << xc[m] << " ";
        
        compute_Y(xc[m], Y);
        double x3=pow(xc[m], 3), Y0=Y[0];
        
        ofile << x3*Y0 << " ";
        
        compute_Yl0_k(xc[m], 0, Y);
        ofile << scale*ISA*x3*The*Y[0] << " ";

        // constant density sphere
        double tauscale=0.5*tau0sphere(b)*tau0sphere(b);
        double tau_0=taul_av_Int_sphere(b, 0)/tauscale;
        double tau_2=taul_av_Int_sphere(b, 2)/tauscale;
        double kappa_T=tau_0+0.1*tau_2-1.0;
        double kappa_S=-0.6*tau_2;
        double kappa_E=tau_0+0.1*tau_2;
        
        cout << tau_0 << " " << tau_2 << endl;
        
        ofile << scale*ISA*x3*kappa_T*Y0 << " ";
        ofile << scale*ISA*x3*kappa_S*The*Y0 << " ";
        ofile << scale*ISA*x3*kappa_E*The*Y[0] << " ";
        
        ofile << scale*ISA*x3*( (kappa_T+kappa_S*The)*Y0 + kappa_E*The*Y[0] ) << " ";

        // isothermal beta model
        tauscale=0.5*tau0(beta, b)*tau0(beta, b);
        tau_0=taul_av_Int_beta(b, beta, 0)/tauscale;
        tau_2=taul_av_Int_beta(b, beta, 2)/tauscale;
        kappa_T=tau_0+0.1*tau_2-1.0;
        kappa_S=-0.6*tau_2;
        kappa_E=tau_0+0.1*tau_2;
  
        cout << tau_0 << " " << tau_2 << endl;
        
        ofile << scale*ISA*x3*kappa_T*Y0 << " ";
        ofile << scale*ISA*x3*kappa_S*The*Y0 << " ";
        ofile << scale*ISA*x3*kappa_E*The*Y[0] << " ";

        ofile << scale*ISA*x3*( (kappa_T+kappa_S*The)*Y0 + kappa_E*The*Y[0] ) << " ";

        ofile << endl;
    }

    ofile.close();
    
    return;
}

//==================================================================================================
//
// plot emission/absorption terms
//
//==================================================================================================
void output_emission_absorption(string filename, int npp)
{
    vector<double> rc(npp);
    init_xarr(1.0e-1, 1.0e+2, &rc[0], npp, 1, 0);
    
    ofstream ofile1(filename.c_str());
    ofile1.precision(10);
    
    double beta=2.0/3.0;
    double tau0v=tau0(beta, 0);
    
    for(int m=0; m<npp; m++)
    {
        double tauv=tau0(beta, rc[m]);
        double dum=z0_Int(rc[m], beta, 0)+0.5*z0_Int(rc[m], beta, 2);
        
        cout << rc[m] << endl;
        
        ofile1 << scientific << rc[m] << " " << tauv/tau0v << " " << pow(tauv/tau0v, 2) << " ";
        
        ofile1 << scientific << (dum-0.5*tauv*tauv)/(0.5*tau0v*tau0v) << " " << dum/(0.5*tau0v*tau0v) << endl;
    }
    
    ofile1.close();
    
    return;
}

//==================================================================================================
void output_emission_absorption_kinetic(string filename, int npp)
{
    vector<double> rc(npp);
    init_xarr(1.0e-1, 1.0e+2, &rc[0], npp, 1, 0);
    
    ofstream ofile1(filename.c_str());
    ofile1.precision(10);
    
    double beta=2.0/3.0;
    double tau0v=tau0(beta, 0);
    
    for(int m=0; m<npp; m++)
    {
        double tauv=tau0(beta, rc[m]);
        
        cout << rc[m] << endl;
        
        ofile1 << scientific << rc[m] << " " << tauv/tau0v << " " << pow(tauv/tau0v, 2) << " ";
        
        for(int l=0; l<=3; l++)
        {
            double dum =lambda_function(rc[m], beta, l);
            double dum1= kappa_function(rc[m], beta, l);
            
            cout << scientific << l << " " << dum/(0.5*tau0v*tau0v) << " " << dum1/(0.5*tau0v*tau0v) << endl;
            ofile1 << scientific << dum/(0.5*tau0v*tau0v) << " " << dum1/(0.5*tau0v*tau0v) << " ";
        }
        
        cout << endl;
        ofile1 << endl;
    }
    
    ofile1.close();
    
    return;
}

//==================================================================================================
//
// plot of Y_l0 functions for different temperatures
//
//==================================================================================================
void output_DI2_E_temperature(string filename, int np, double Te)
{
    vector<double> xc(np);
    init_xarr(1.0e-1, 40.0, &xc[0], np, 1, 0);
    
    ofstream ofile(filename.c_str());
    ofile.precision(8);
    
    double a0=pow(Te/const_me, 2);
    double a1=-0.4*pow(Te/const_me, 2);
    double a2= 0.1*pow(Te/const_me, 2);
    double a3=-3.0/70.0*pow(Te/const_me, 2);
    
    for(int m=0; m<np; m++)
    {
        cout << xc[m] << endl;
        double x3=pow(xc[m], 3);
        
        ofile << xc[m] << " ";
        
        ofile << x3*compute_SZ_distortion_asym_multiple(xc[m], 0, Te/const_me, 0 )/a0 << " ";
        //
        ofile << x3*compute_SZ_distortion_asym_multiple(xc[m], 0, Te/const_me, 10)/a0 << " ";
        ofile << x3*compute_SZ_distortion_asym_multiple(xc[m], 1, Te/const_me, 10)/a1 << " ";
        ofile << x3*compute_SZ_distortion_asym_multiple(xc[m], 2, Te/const_me, 10)/a2 << " ";
        ofile << x3*compute_SZ_distortion_asym_multiple(xc[m], 3, Te/const_me, 10)/a3 << " ";
        //
        ofile << x3*compute_SZ_distortion_Patterson_multiple_Kernel_2D(xc[m], 0, Te/const_me)/a0 << " ";
        ofile << x3*compute_SZ_distortion_Patterson_multiple_Kernel_2D(xc[m], 1, Te/const_me)/a1 << " ";
        ofile << x3*compute_SZ_distortion_Patterson_multiple_Kernel_2D(xc[m], 2, Te/const_me)/a2 << " ";
        ofile << x3*compute_SZ_distortion_Patterson_multiple_Kernel_2D(xc[m], 3, Te/const_me)/a3 << " ";

        ofile << endl;
    }
    
    ofile.close();
    
    return;
}

//==================================================================================================
//
// plot of Y_l0 functions for different temperatures
//
//==================================================================================================
void output_SZ_signal(string filename, int np, double tau, double Te)
{
    vector<double> xc(np);
    init_xarr(1.0e-1, 40.0, &xc[0], np, 1, 0);
    
    ofstream ofile(filename.c_str());
    ofile.precision(8);
    
    for(int m=0; m<np; m++)
    {
        cout << xc[m] << endl;

        ofile << xc[m] << " ";
        
        double Dn0=compute_SZ_signal_combo(xc[m], tau, Te, 0, 0, 0, 0);
        double x3=pow(xc[m], 3), fac=500.0*tau*tau/2.0;
        
        ofile << Dn_DI_conversion*x3*Dn0 << " ";
        
        double S0=compute_SZ_distortion_Patterson_multiple_Kernel_2D(xc[m], 0, Te/const_me);
        double S2=compute_SZ_distortion_Patterson_multiple_Kernel_2D(xc[m], 2, Te/const_me);
        
        ofile << Dn_DI_conversion*x3*fac*S0 << " ";

        double kappa_T=0.69+0.1*0.14-1.0;
        ofile << Dn_DI_conversion*x3*fac*( kappa_T*Dn0/tau + 0.69*S0+0.1*0.14*S2 ) << " ";
        
        kappa_T=0.87+0.1*0.15-1.0;
        ofile << Dn_DI_conversion*x3*fac*( kappa_T*Dn0/tau + 0.87*S0+0.1*0.15*S2 ) << " ";

        ofile << endl;
    }
    
    ofile.close();
    
    return;
}

//==================================================================================================
//
// main call of routine. Here parameters are handled.
//
//==================================================================================================
int main(int narg, char *args[])
{
    output_SZ_signal("./DI_signal.5keV.dat", 400, 0.01, 5.0);
    output_SZ_signal("./DI_signal.10keV.dat", 400, 0.01, 10.0);
    output_SZ_signal("./DI_signal.15keV.dat", 400, 0.01, 15.0);
    output_SZ_signal("./DI_signal.25keV.dat", 400, 0.01, 25.0);
    exit(0);
    
    output_DI2_E_temperature("./DI2_E.5keV.dat", 400, 5.0);
//    output_DI2_E_temperature("./DI2_E.5keV.dat", 400, 5.0);
//    output_DI2_E_temperature("./DI2_E.10keV.dat", 400, 10.0);
//    output_DI2_E_temperature("./DI2_E.15keV.dat", 400, 15.0);
//    output_DI2_E_temperature("./DI2_E.25keV.dat", 400, 25.0);
    exit(0);
    
    output_lowest_order_signal("./lowest_order.dat", 400, 0);
//    exit(0);

//    output_plot_tau_l_sphere("./tau_l_sphere.dat", 12, 400);
//    exit(0);
    
//    output_plot_tau_l("./tau_l_data.dat", 30, 50);
//    exit(0);
    
//    output_emission_absorption("./spatial_data.dat", 50);
    output_emission_absorption_kinetic("./spatial_data.kinetic.dat", 50);
    exit(0);
/*    
    int lmax=100;
    double beta=2.0/3.0, r=0.0;
    
    vector<double> taulc(lmax+2);
    vector<double> taulv(lmax+2);
    double b=1.0;
    
    for(int k=0; k<=lmax; k+=2)
    {
        taulc[k]=z0_Int(b, beta, k)/(0.5*pow(tau0(beta, b), 2));
        r+=(2.0*k+1.0)*taulc[k];
        cout << k << " " << (2.0*k+1.0)*taulc[k] << " " << r << endl;
    }
    exit(0);
*/    
    output_CMB_isotropy_Y_functions("./CMB_scattering.Y.dat", 400);
    exit(0);
    
    output_CMB_isotropy_signals("./CMB_scattering.Te_20keV.dat", 400, 20.0);
    exit(0);
    
    
    //==============================================================================================
    // second scattering to thermal SZ and CMB scattering
    //==============================================================================================
    int np=400;
    vector<double> xc(np);
    init_xarr(1.0e-1, 40.0, &xc[0], np, 1, 0);
    
    for(int l=0; l<=4; l++) cout << P_l_Kernel_Norm(l, 0.01, 1e-5) << endl;
    
    cout << P_l_Kernel(4, 0, 0.01) << endl;
        
    ofstream ofileK("./Kernel.dat");
    ofileK.precision(8);
    
    double Thev=50.0/const_me;
    double sig[4];
    
    sig[0]=P_l_Kernel(0, 0, Thev);
    sig[1]=P_l_Kernel(1, 0, Thev);
    sig[2]=P_l_Kernel(2, 0, Thev);
    sig[3]=P_l_Kernel(3, 0, Thev);
    
    double sFWHM=2.0*atanh(sqrt(2.0*log(2.0)*Thev));
    np++;
    init_xarr(-4.0*sFWHM, 4.0*sFWHM, &xc[0], np, 0, 0);

    for(int m=0; m<0*np; m++)
    {
        cout << xc[m] << endl;
        
        ofileK << xc[m] << " ";
        
        for(int l=0; l<=4; l++)
            ofileK << P_l_Kernel(l, xc[m], Thev) << " "
                   << Compute_Kernel_interpol(l, xc[m], Thev) << " ";

 /*       ofileK << xc[m]/sFWHM << " ";
        
        for(int l=0; l<=4; l++)
            ofileK << P_l_Kernel(l, xc[m], Thev)/sig[l] << " "
                   << Compute_Kernel_interpol(l, xc[m], Thev)/sig[l] << " ";
 */       
        ofileK << endl;
    }
    ofileK.close();
    
//    exit(0);
    
    init_xarr(1.0e-1, 40.0, &xc[0], np, 1, 0);

    double Te=25.0, tau=0.01;
    ofstream ofile("./CMB.data.25keV.compare.dat");
    ofile.precision(8);
    
    for(int m=0; m<np; m++)
    {
        cout << xc[m] << endl;
        
        ofile << xc[m] << " ";
        
//        for(int l=0; l<=3; l++)
//            ofile << pow(xc[m], 3)*compute_SZ_distortion_Patterson_multiple(xc[m], l, Te/const_me, "SEC") << " ";

/*        for(int l=0; l<=3; l++)
            ofile << pow(xc[m], 3)*compute_SZ_distortion_Patterson_multiple_Kernel(xc[m], l, Te/const_me) << " ";
        
        for(int l=0; l<=3; l++)
            ofile << pow(xc[m], 3)*compute_SZ_distortion_asym_multiple(xc[m], l, Te/const_me, 10) << " ";
*/        
        double Dn0=compute_SZ_signal_combo(xc[m], tau, Te, 0, 0, 0, 0);
        double x3=pow(xc[m], 3), fac=Dn_DI_conversion*tau*tau/2.0;
        
        ofile << Dn_DI_conversion*x3*Dn0 << " ";
        
        ofile << fac*x3*compute_SZ_distortion_asym_multiple(xc[m], 0, Te/const_me, 10) << " ";

        ofile << fac*x3*(
                     0.61*compute_SZ_distortion_asym_multiple(xc[m], 0, Te/const_me, 10)
                    +0.087/10.0*compute_SZ_distortion_asym_multiple(xc[m], 2, Te/const_me, 10)
                  ) << " ";

        ofile << Dn_DI_conversion*x3*tau/2.0*(0.61+0.087/10.0-1.0)*Dn0 << " ";
        
        double Dsig2=compute_SZ_distortion_Patterson_multiple(1.0, 2, Te/const_me, "SIG")-0.1;
        
        ofile << fac*x3*Dsig2*Dn0 << " ";
 
        double S2_0=compute_SZ_distortion_Patterson_multiple_Kernel_2D(xc[m], 0, Te/const_me);
        double S2_2=compute_SZ_distortion_Patterson_multiple_Kernel_2D(xc[m], 2, Te/const_me);
        
        ofile << fac*x3*S2_0 << " ";

        ofile << fac*x3*(0.61*S2_0+0.087/10.0*S2_2) << " ";
        
        ofile << fac*x3*(0.61*S2_0+0.087/10.0*S2_2+Dsig2*Dn0)
                 +Dn_DI_conversion*x3*tau/2.0*(0.61+0.087/10.0-1.0)*Dn0 << " ";
        
        //------------------------------------------------------------
        // x-derivatives
        //------------------------------------------------------------
        double Ix =compute_SZ_signal_combo(xc[m], 1.0, Te, 0, 0, 0, 0);
        double Ixp=compute_SZ_signal_combo(xc[m]*(1.0+0.001), 1.0, Te, 0, 0, 0, 0);
        double Ixm=compute_SZ_signal_combo(xc[m]*(1.0-0.001), 1.0, Te, 0, 0, 0, 0);
        double dI =(Ixp-Ixm)/(2.0*0.001);
        double d2I=(Ixp-2.0*Ix+Ixm)/pow(0.001, 2);

        ofile << fac*x3*Te/const_me*xc[m]*(4.0*dI+xc[m]*d2I)*0.62 << endl;
    }
    
    ofile.close();
    
    cout << Dn_DI_conversion << endl;
    
    return 0;

/*
    //==============================================================================================
    // compute scattering cross section correction
    //==============================================================================================
    int np=50;
    vector<double> Thev(np);
    init_xarr(1.0e-3, 1.0, &Thev[0], np, 1, 0);
   
    ofstream ofile("./sig.data.dat");
    ofile.precision(8);
    
    for(int m=0; m<np; m++)
    {
        ofile << Thev[m] << " ";
        
        for(int l=0; l<=4; l++)
            ofile << compute_SZ_distortion_Patterson_multiple(1.0, l, Thev[m], "SIG") << " ";
        
        ofile << endl;
    }
    
    ofile.close();
*/
    return 0;
}

//==================================================================================================
//==================================================================================================
