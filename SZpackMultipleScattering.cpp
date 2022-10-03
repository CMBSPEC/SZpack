//======================================================================================
// Author: Elizabeth Lee
// Based on work by Jens Chluba 
//
//======================================================================================

#include "SZpackMultipleScattering.h"

//==================================================================================================
//
// const density sphere
//
//==================================================================================================

namespace constDensitySphere{
    double z0, mur, b0;
    int l;
    
    double tau0(double B0) {
        return 2.0*sqrt(1.0-B0*B0);
    }

    double tau0(double Z0, double Mur) {
        return sqrt(1.0-Z0*Z0*(1.0-Mur*Mur))-Z0*Mur;
    }

    double dtau_dmu(double Mur) {
        return LegendreP(l, Mur)*( tau0(z0, Mur)+pow(-1.0, l)*tau0(z0, -Mur) );
    }

    double dtau_ds(double S) {
        double z0eff  = sqrt(b0*b0+S*S)+1.0e-300;
        double mureff = S/z0eff;
        
        return LegendreP(l, mureff)*tau_hat(z0eff, l);
    }

    double tau_hat(double Z0, int L) //Tau averaged over mu_r
    {
        z0 = Z0;
        l  = L;
        
        double a = 0.0, b = 1.0;
        double epsrel = 1.0e-5, epsabs = 1.0e-16;
        
        return Integrate_using_Patterson_adaptive(a, b, epsrel, epsabs, dtau_dmu)/2.0;
    }

    double tau_av(double B0, int L) // Tau averaged over mu_r and s
    {
        if(L > 0 && L % 2){
            return 0;
        }

        b0 = B0;
        l = L;

        double a = 0, b = sqrt(1.0-b0*b0);
        double epsrel = 1.0e-5, epsabs = 1.0e-16;
        
        return 2.0*(2.0*l+1.0)*Integrate_using_Patterson_adaptive(a, b, epsrel, epsabs, dtau_ds);
    }
}

//==================================================================================================
//
// isothermal beta profile
//
//==================================================================================================

namespace isothermalBeta{
    double beta = 2.0/3.0;
    double z0, mur;
    double b0;
    int l;

    void setBeta(double Beta) {
        beta = Beta;
    }

    double Ne(double r2) {
        return pow(1.0+r2, -1.5*beta);
    }

    double tau0(double b) {
        return sqrt(PI)*Gamma_JC(1.5*beta-0.5)/Gamma_JC(1.5*beta)*pow(1.0+b*b, 0.5-1.5*beta);
    }

    //==================================================================================================
    double dtau_dlns(double lns) {
        double s = exp(lns);
        double N = Ne(s*s+z0*z0+2.0*s*z0*mur);
        return s*N;
    }

    double dtau_dsdmu(double mur) {
        double Pl = LegendreP(l, mur);
        return Pl*( s_Int(z0, mur)+pow(-1.0, l)*s_Int(z0, -mur) );
    }

    double dtau_dlnz0(double lnz0) {
        double Z0 = exp(lnz0);
        double z0eff  = sqrt(b0*b0+Z0*Z0)+1.0e-300;
        double mureff = Z0/z0eff;
    
        double N  = Ne(z0eff);
        double Pl = LegendreP(l,  mureff);
        
        return Z0*N*Pl*mur_Int(z0eff, l);
    }

    double dtau_dlnz0_m1(double lnz0) {
        double Z0 = exp(lnz0);
        double z0eff = sqrt(b0*b0+Z0*Z0)+1.0e-300;
        double mureff = Z0/z0eff;
        
        double N   = Ne(z0eff);
        double Plm = LegendreP_m1(l,  mureff);
        
        return Z0*N*Plm*mur_Int(z0eff, l);
    }

    //==================================================================================================

    double s_Int(double Z0, double Mur) //Tau averaged over s
    {
        double a = log(1.0e-6), b, r = 0.0;
        double epsrel = 1.0e-4;
        
        z0 = Z0;
        mur = Mur;
        
        for(int k=0; k<=20; k++)
        {
            b=log(10.0)+a;
            double r1=Integrate_using_Patterson_adaptive(a, b, epsrel, epsrel*(1.0e-100+r)*0.01, dtau_dlns);
            r+=r1;
            
            a=b;
            
            if(fabs(r1/r)<epsrel*0.01 && a>10.0) break;
        }
        
        return r;
    }

    double mur_Int_explicit(double Z0, int L) //Tau explicitly averaged over mu and s
    {
        double a=0, b=1.0;
        double epsrel=1.0e-4, epsabs=1.0e-10;
        
        z0 = Z0;
        l = L;

        return Integrate_using_Patterson_adaptive(a, b, epsrel, epsabs, dtau_dsdmu)/2.0;
    }

    //==================================================================================================
    vector<bool> tau_spline_set(100, 0);
    vector<int> memindex_tau_spline(100, 0);

    double mur_Int(double Z0, int L)  //Tau averaged over mu and s splined between calculated values
    {
        z0 = Z0;
        l = L;

        if(!tau_spline_set[l])
        {
            int np = 1000;
            double a = log(1.0e-6), b = log(1.0e+4);

            cout << " setting up spline for l = " << l << " (# points = " << np << ")" << endl;
            
            vector<double> xc(np), yc(np);
            init_xarr(a, b, &xc[0], np, 0, 0);
            
            for(int k=0; k<np; k++) yc[k]=mur_Int_explicit(exp(xc[k]), l);
            
            memindex_tau_spline[l]=calc_spline_coeffies_JC(np, &xc[0], &yc[0]);
            tau_spline_set[l]=1;
        }
        
        z0 = min(max(2.0e-6, z0), 9.9e+3);

        return calc_spline_JC(log(z0), memindex_tau_spline[l]);
    }

    //==================================================================================================
    double z0_Int(double B0, int L)//Tau averaged over mu, s and z.
    {
        b0 = B0;
        l = L;
        
        if(l>0 && l%2) return 0;

        double a=log(1.0e-6), b=log(1.0e+4);
        double epsrel=2.0e-4, epsabs=1.0e-10;

        return 2.0*Integrate_using_Patterson_adaptive(a, b, epsrel, epsabs, dtau_dlnz0);
    }

    double z0_Int_m1(double B0, int L) //Tau averaged over mu, s and z using the second multipole.
    {
        b0 = B0;
        l = L;

        if(l>0 && (l+1)%2) return 0;

        double epsrel=2.0e-4, epsabs=1.0e-10;
        double a=log(1.0e-6), b=log(1.0e+4);

        return 2.0*Integrate_using_Patterson_adaptive(a, b, epsrel, epsabs, dtau_dlnz0_m1);
    }

    //==================================================================================================
    double tau_av(double B0, int L)  //Tau averaged over mu, s and z, with prefactor
    {
        b0 = B0;
        l = L;
        return (2.0*l+1.0)*z0_Int(b0, l);
    }

    //==================================================================================================
    double lambda_function(double B0, int L)
    {
        b0 = B0;
        l = L;

        double dum = 0.0;
        if(l>=1) dum = l*z0_Int(b0, l-1);
        
        dum += (l+1.0)*z0_Int(b0, l+1);
        
        return dum;
    }

    double kappa_function(double B0, int L)
    {
        b0 = B0;
        l = L;

        double dum = 0.0;
        if(l>=1) dum = z0_Int_m1(b0, l-1);
        
        dum-=z0_Int_m1(b0, l+1);
        
        return dum;
    }
}
