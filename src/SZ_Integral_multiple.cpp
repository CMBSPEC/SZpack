//==================================================================================================
//
// Computation multiple scattering correction to thermal SZ effect according to Chluba et al. 2013
//
//==================================================================================================
//
// Author: Jens Chluba  (Johns Hopkins University)
//
// first implementation: Jan 2013
// last modification   : Jan 2013
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
#include "Patterson.h"
#include "SZ_Integral_multiple.h"
#include "SZ_Integral_multiple_data.h"
#include "SZpack.h"
#include "Parameters.h"

//==================================================================================================
//
// namespaces
//
//==================================================================================================
using namespace std;


//==================================================================================================
namespace SZIntegralmultiple
{

double Int_eps_3D; // relative accuracy for numerical integration (< 10^-6 is hard to achieve)

string run_mode="SIG";

//==================================================================================================
//
// azimuthally averaged cross section for multipole l<=3
//
//==================================================================================================
double sig_Boltzmann_Compton(int l, double The, double xi, double mu, double mug)
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
    double zeta=1.0/gamma2/kappa/kappap;
    //
    double nfac=3.0/8.0/PI*kappa/pow(kappap, 2)/gamma2;
    //
    double dsig_0_dmudmup=1.0-zeta*(alpha_sc-0.5*beta_sc*zeta);
    
    if(l==0) return nfac*dsig_0_dmudmup;
    
    else if(l==1) return nfac*(mu*mug*dsig_0_dmudmup + (1.0-alpha_sc*zeta)*zeta*gmumup);
                               
    else if(l==2)
    {
        double P2mu=0.5*(3.0*mu*mu-1.0);
        double P2mug=0.5*(3.0*mug*mug-1.0);
        //
        double rt=16.0*mu*mug+(1.0-16.0*mu*mug-mug*mug-mu*mu*(1.0-17.0*mug*mug))*zeta;
        //
        return nfac*(P2mu*P2mug*dsig_0_dmudmup + 3.0/16.0*rt*zeta*gmumup);
    }
    
    else if(l==3)
    {
        double mu2=mu*mu, mug2=mug*mug;
        double P3mu =0.5*(5.0*mu *mu -3.0)*mu;
        double P3mug=0.5*(5.0*mug*mug-3.0)*mug;
        //
        double rt1=2.0*(1.0-5.0*mu2)*(1.0-5.0*mug2)*(1.0-zeta);
        double rt2=mu*mug*(7.0-15.0*mug2-5.0*mu2*(3.0-11.0*mug2))*zeta;
        
        return nfac*(P3mu*P3mug*dsig_0_dmudmup + 3.0/16.0*(rt1+rt2)*zeta*gmumup);
    }
    
    else if(l==4)
    {
        double mu2=mu*mu, mug2=mug*mug;
        double P4mu =(3.0-5.0*mu2 *(6.0-7.0*mu2 ))/8.0;
        double P4mug=(3.0-5.0*mug2*(6.0-7.0*mug2))/8.0;
        //
        double rt1=mu*(3.0-7.0*mu2)*mug*(3.0-7.0*mug2)*(1.0-zeta);
        double rt2=(1.0-mu2*(8.0-7.0*mu2)-8.0*mug2+8.0*(17.0-28.0*mu2)*mu2*mug2
                            +7.0*(1.0-mu2*(32.0-63.0*mu2))*mug2*mug2)*zeta/8.0;
        
        return nfac*(P4mu*P4mug*dsig_0_dmudmup + 5.0/8.0*(rt1+rt2)*zeta*gmumup);
    }
    
    return 0.0;
}

//==================================================================================================
// 
// collision terms for different runmodes
//
//==================================================================================================
double sig_Boltzmann_Compton(double x, double The, double xi, double mu, double mug, int l)
{
    double eta=sqrt(2.0*The*xi);

    //==============================================================================================
    // mu == cosine of angle of gamm   and beta 
    // mug== cosine of angle of gamma' and beta
    //==============================================================================================    
    double beta=eta/sqrt(1.0+eta*eta);
    //
    double kappap=1.0-beta*mug;
    double dx=x*beta*(mug-mu)/kappap;
    double xp=x+dx;
    //
    double exp_mx=exp(-x);
    double dex=one_minus_exp_mx(x, exp_mx);
    double exp_mxp=exp(-xp);
    double dexp=one_minus_exp_mx(xp, exp_mxp); 
    
    //==============================================================================================
    // compute cross section
    //==============================================================================================
    if(run_mode=="SIG") return sig_Boltzmann_Compton(l, The, xi, mu, mug);
    
    //==============================================================================================
    // scattering of CMB anisotropies; relevant spectral redistribution
    //==============================================================================================
    if(run_mode=="CMB")
    {
        double G=x*exp_mx/pow(dex, 2);
        double Gp=xp*exp_mxp/pow(dexp, 2);
        return sig_Boltzmann_Compton(l, The, xi, mu, mug)*(Gp-G);
    }
    
    //==============================================================================================
    // second scattering correction full
    //==============================================================================================
    else if(run_mode=="SEC")
    {
        Parameters temp = Parameters();
        temp.xcmb[0] = x;
        temp.Dtau = 1.0;
        temp.Te = The*const_me;
        temp.betac = temp.muc = temp.betao = temp.muo = 0;
        temp.relative_accuracy = Int_eps_3D*0.75;
        temp.setCalcValues();
        double S=compute_signal_3D (0, temp);
        temp.xcmb[0] = xp;
        double Sp=compute_signal_3D(0, temp);
        return sig_Boltzmann_Compton(l, The, xi, mu, mug)*(Sp-S);
        
        if(xp>=30.0 || x>=30.0 || x<=0.01 || xp<=0.01) return 0.0;
        
//        double S =compute_SZ_signal_combo(x , 1.0, The*const_me, 0, 0, 0, 0);
//        double Sp=compute_SZ_signal_combo(xp, 1.0, The*const_me, 0, 0, 0, 0);
//        return sig_Boltzmann_Compton(l, The, xi, mu, mug)*(Sp-S);
    }
        
    return 0.0;
}

//==================================================================================================
//
// integrations using Patterson
//
//==================================================================================================
double Boltzmann_I(double mup, void *userdata)
{ 
    double *d=(double *)userdata;
    
    return sig_Boltzmann_Compton(d[0], d[1], d[2], d[3], mup, (int)d[4]);
}

//==================================================================================================
double mup_Int(double x, double The, double xi, double mu, int l)
{
    double a=-1.0, b=1.0;
    double epsrel=Int_eps_3D*0.8, epsabs=1.0e-300;
    
    double d[5]={x, The, xi, mu, (double)l};
    void *userdata=(void *)d;
    
    return Integrate_using_Patterson_adaptive(a, b, epsrel, epsabs, Boltzmann_I, userdata);
}

//==================================================================================================
double Boltzmann_II(double mu, void *userdata)
{ 
    double *d=(double *)userdata;
    
    return mup_Int(d[0], d[1], d[2], mu, (int)d[3]);
}

//==================================================================================================
double mu_Int(double x, double The, double xi, int l)
{
    double a=-1.0, b=1.0;
    double epsrel=Int_eps_3D*0.9, epsabs=1.0e-300;
    
    double d[4]={x, The, xi, (double)l};
    void *userdata=(void *)d;
    
    return Integrate_using_Patterson_adaptive(a, b, epsrel, epsabs, Boltzmann_II, userdata);
}

//==================================================================================================
double Boltzmann_III(double xi, void *userdata)
{ 
    double *d=(double *)userdata;
    
    return mu_Int(d[0], d[1], xi, (int)d[2])*df_RM_dTheta(xi, d[1], 0);
}

//==================================================================================================
double xi_Int(double x, double The, int l)
{
    double a=0.0, lim=30.0, b=lim*(1.0+0.5*lim*The);
    double epsrel=Int_eps_3D, epsabs=1.0e-300;
    
    double d[3]={x, The, (double)l};
    void *userdata=(void *)d;
    
    return Integrate_using_Patterson_adaptive(a, b, epsrel, epsabs, Boltzmann_III, userdata);
}


//==================================================================================================
//
// functions for Compton scattering kernel
//
//==================================================================================================
double Boltzmann_I_Kernel(double mup, void *userdata)
{
    double *d=(double *)userdata;
    
    double beta=sqrt(2.0*d[1]*d[2])/sqrt(1.0+2.0*d[1]*d[2]);
    
    double mu=(1.0-exp(d[0])*(1.0-beta*mup))/beta;
    
    return (1.0-beta*mu)/beta*sig_Boltzmann_Compton((int)d[3], d[1], d[2], mu, mup);
}

double mu_Int_Kernel(double s, double The, double xi, int l)
{
    double beta=sqrt(2.0*The*xi)/sqrt(1.0+2.0*The*xi);
    
    double a=max(-1.0, (1.0-exp(-s)*(1.0+beta))/beta);
    double b=min( 1.0, (1.0-exp(-s)*(1.0-beta))/beta);
    
    double epsrel=Int_eps_3D*0.9, epsabs=1.0e-300;
    
    double d[4]={s, The, xi, (double)l};
    void *userdata=(void *)d;
    
    return Integrate_using_Patterson_adaptive(a, b, epsrel, epsabs, Boltzmann_I_Kernel, userdata);
}

double Boltzmann_II_Kernel(double xi, void *userdata)
{
    double *d=(double *)userdata;
    
    return mu_Int_Kernel(d[0], d[1], xi, (int)d[2])*df_RM_dTheta(xi, d[1], 0);
}

double xi_Int_Kernel(int l, double s, double The)
{
    double a=pow(sinh(fabs(s)*0.5), 2)/2.0/The, lim=30.0, b=lim*(1.0+0.5*lim*The);
    
    if(b<a) return 0.0; 
    
    double epsrel=Int_eps_3D, epsabs=1.0e-300;
    
    double d[3]={s, The, (double)l};
    void *userdata=(void *)d;
    
    return Integrate_using_Patterson_adaptive(a, b, epsrel, epsabs, Boltzmann_II_Kernel, userdata);
}

double Boltzmann_III_Kernel(double s, void *userdata)
{
    double *d=(double *)userdata;

    return xi_Int_Kernel((int)d[1], s, d[0]);
}

double s_Int_Kernel(int l, double The)
{
    double lim=30.0, slim=2.0*asinh(sqrt(The*lim*(2.0+lim*The)));
    double a=-slim, b=slim;
    
    double epsrel=Int_eps_3D*1.05, epsabs=1.0e-300;
    
    double d[2]={The, (double)l};
    void *userdata=(void *)d;
    
    return Integrate_using_Patterson_adaptive(a, b, epsrel, epsabs, Boltzmann_III_Kernel, userdata);
}

}

//==================================================================================================
//
// global access to functions
//
//==================================================================================================
using namespace SZIntegralmultiple;


//==================================================================================================
//
// Compton scattering kernel
//
//==================================================================================================
double P_l_Kernel(int l, double s, double The, double eps_Int)
{
    Int_eps_3D=eps_Int;
    double norm_f=norm_df_RM_dTheta(The)/The;
    //TODO: This /The scaling is because I changed this function many many commits before
    double r=TWOPI*TWOPI*norm_f*xi_Int_Kernel(l, s, The);
    
    return r;
}

double P_l_Kernel_Norm(int l, double The, double eps_Int)
{
    Int_eps_3D=eps_Int;
    
    double norm_f=norm_df_RM_dTheta(The)/The;
    //TODO: This /The scaling is because I changed this function many many commits before
    double r=TWOPI*TWOPI*norm_f*s_Int_Kernel(l, The);
    
    return r;
}

//==================================================================================================
double sFWHM(double The){ return 2.0*atanh(sqrt(2.0*log(2.0)*The)); }

void create_2D_Kernel_Lib(string filename, int lmax, int nps, int npT,
                          double The_min, double The_max, double eps_Int)
{
    vector<double> xc(nps);
    init_xarr(1.0e-4, 4.0, &xc[0], nps, 1, 0);

    vector<double> Tc(npT);
    init_xarr(The_min, The_max, &Tc[0], npT, 1, 0);

    ofstream ofile(filename.c_str());
    ofile.precision(10);
    
    ofile << lmax << " " << 2*nps+1 << " " << npT << endl;

    for(int k=nps-1; k>=0; k--) ofile << scientific << -xc[k] << " ";
    ofile << 0 << " ";
    for(int k=0; k<nps; k++) ofile << scientific << xc[k] << " ";
    ofile << endl;
    
    for(int k=0; k<npT; k++) ofile << scientific << Tc[k] << " ";
    ofile << endl;
    
    for(int l=0; l<=lmax; l++)
    {
        for(int kT=0; kT<npT; kT++)
        {
            double The=Tc[kT];
            double sval=sFWHM(The);
            
            for(int k=nps-1; k>=0; k--)
                ofile << scientific << P_l_Kernel(l, -xc[k]*sval, The, eps_Int) << " ";
            
            ofile << scientific << P_l_Kernel(l, 0.0, The, eps_Int) << " ";
            
            for(int k=0; k<nps; k++)
                ofile << scientific << P_l_Kernel(l,  xc[k]*sval, The, eps_Int) << " ";
            
            ofile << endl;
        }
    }
    
    ofile.close();
    
    return;
}

//==================================================================================================
bool Lib_is_loaded=0;
vector<double> The_lib;
vector<double> s_scale_lib;
vector<vector<vector<double> > > Kernel_lib;

void read_2D_Kernel_Lib(string filename)
{
    The_lib.clear();
    s_scale_lib.clear();
    Kernel_lib.clear();
    ifstream ifile(filename.c_str());
    
    if(!ifile)
    {
        create_2D_Kernel_Lib(filename, 5, 300, 300, 0.002, 0.2, 1.0e-5);
        read_2D_Kernel_Lib(filename);
        
        return;
    }
    
    int lmax, nps, npT;
    double dum;
    vector<double> dumv;
    vector<vector<double> > dumvv;
    
    ifile >> lmax;
    ifile >> nps;
    ifile >> npT;
    
    for(int k=0; k<nps; k++)
    {
        ifile >> dum;
        s_scale_lib.push_back(dum);
    }

    for(int k=0; k<npT; k++)
    {
        ifile >> dum;
        The_lib.push_back(log(dum));
    }

    for(int l=0; l<=lmax; l++)
    {
        dumvv.clear();
        
        for(int kT=0; kT<npT; kT++)
        {
            dumv.clear();
            
            for(int k=0; k<nps; k++)
            {
                ifile >> dum;
                dumv.push_back(dum);
            }
            
            dumvv.push_back(dumv);
        }
        
        Kernel_lib.push_back(dumvv);
    }
    
    ifile.close();
    
    return;
}

//==================================================================================================
double Compute_Kernel_interpol(int l, double s, double The)
{
    if(l>3) return 0;
    
    if(Lib_is_loaded==0)
    {
        read_2D_Kernel_Lib("./Kernel_Lib.dat");
        Lib_is_loaded=1;
    }

    double s_scale=s/sFWHM(The);
    
    if(s_scale < s_scale_lib[0] || s_scale_lib[s_scale_lib.size()-1]< s_scale) return 0.0;
    if(log(The) < The_lib[0] || The_lib[The_lib.size()-1]< log(The)) return 0.0;

    unsigned long j=0;
        
    locate_JC(&s_scale_lib[0], s_scale_lib.size(), s_scale, &j);
    
    int nT=The_lib.size();
    int index_T=floor( (log(The)-The_lib[0])/((The_lib[nT-1]-The_lib[0])/nT) );
    
    j++;
    index_T++;
    
    j=min(s_scale_lib.size()-4, j);
    index_T=min(nT-4, index_T);
    
    int jmid=(s_scale_lib.size()-1)/2;
    if((int)j<=jmid) j=min(jmid-4, j);
        
    double x = s_scale;
    double y = log(The);
    double x1 = s_scale_lib[j+0];
    double x2 = s_scale_lib[j+1];
    double x3 = s_scale_lib[j+2];
    double x4 = s_scale_lib[j+3];
    //
    double y1 = The_lib[index_T+0];
    double y2 = The_lib[index_T+1];
    double y3 = The_lib[index_T+2];
    double y4 = The_lib[index_T+3];
    
    /* ------------------------------------------- */
    /* Coefficients for 1-D interpolations         */
    /* ------------------------------------------- */
    double a[4], b[4];
    //
    a[0]= (x - x2)*(x - x3)*(x - x4)/((x1 - x2)*(x1 - x3)*(x1 - x4));
    a[1]= (x - x1)*(x - x3)*(x - x4)/((x1 - x2)*(x2 - x3)*(x4 - x2));
    a[2]= (x - x1)*(x - x2)*(x - x4)/((x1 - x3)*(x2 - x3)*(x3 - x4));
    a[3]= (x - x1)*(x - x2)*(x - x3)/((x3 - x4)*(x4 - x1)*(x2 - x4));
    //
    b[0]= (y - y2)*(y - y3)*(y - y4)/((y1 - y2)*(y1 - y3)*(y1 - y4));
    b[1]= (y - y1)*(y - y3)*(y - y4)/((y1 - y2)*(y2 - y3)*(y4 - y2));
    b[2]= (y - y1)*(y - y2)*(y - y4)/((y1 - y3)*(y2 - y3)*(y3 - y4));
    b[3]= (y - y1)*(y - y2)*(y - y3)/((y3 - y4)*(y4 - y1)*(y2 - y4));
    
    double fxy=0.0, fx;
    
    for(int i=0; i<4; i++)
    {
        fx=0.0;
        for(int k=0; k<4; k++) fx+=a[k]*Kernel_lib[l][i+index_T][k+j];
        fxy+=b[i]*fx;                                                          // interpolation in y
    }
    
    return fxy;
}


//==================================================================================================
//
// multiple scattering correction using the scattering kernel library
//
// eps_Int : relative accuracy for numerical integration (lower than 10^-6 is hard to achieve)
//
//==================================================================================================
double Boltzmann_DI_Kernel(double s, void *userdata)
{
    double *d=(double *)userdata;
    
    double x=d[0], xp=d[0]*exp(s);
    
    if(xp>=30.0 || x>=30.0 || x<=0.01 || xp<=0.01) return 0.0;

    Parameters temp = Parameters();
    temp.Dtau = 1.0;
    temp.Te = d[1]*const_me;
    temp.betac = temp.betao = temp.muc = temp.muo = 0;
    temp.setCalcValues();

    temp.xcmb[0] = x;
    double S =compute_signal_combo(0, temp);
    temp.xcmb[0] = xp;
    double Sp=compute_signal_combo(0, temp);
    
    return Compute_Kernel_interpol((int)d[2], s, d[1])*(Sp-S);
}

double DI_Int_Kernel(int l, double x, double The)
{
    double a=-4.0*sFWHM(The), b=4.0*sFWHM(The);
    
    double epsrel=Int_eps_3D*1.05, epsabs=1.0e-300;
    
    double d[3]={x, The, (double)l};
    void *userdata=(void *)d;
    
    return Integrate_using_Patterson_adaptive(a, b, epsrel, epsabs, Boltzmann_DI_Kernel, userdata);
}

//==================================================================================================
double Boltzmann_DI_Kernel_single(double s, void *userdata)
{
    double *d=(double *)userdata;
    
    double x=d[0], xp=d[0]*exp(s);
    
    double S =exp(-x) /(1.0-exp(-x) );
    double Sp=exp(-xp)/(1.0-exp(-xp));
    
    return Compute_Kernel_interpol(0, s, d[1])*(Sp-S);
}

double DI_Int_Kernel_single(double x, double The)
{
    double a=-4.0*sFWHM(The), b=4.0*sFWHM(The);
    
    double epsrel=Int_eps_3D*0.9, epsabs=1.0e-300;
    
    double d[2]={x, The};
    void *userdata=(void *)d;
    
    return Integrate_using_Patterson_adaptive(a, b, epsrel, epsabs, Boltzmann_DI_Kernel_single, userdata);
}

double Boltzmann_DI_Kernel_II(double s, void *userdata)
{
    double *d=(double *)userdata;
    
    double x=d[0], xp=d[0]*exp(s);
    
    double S =DI_Int_Kernel_single(x , d[1]);
    double Sp=DI_Int_Kernel_single(xp, d[1]);
    
    return Compute_Kernel_interpol((int)d[2], s, d[1])*(Sp-S);
}

double DI_Int_Kernel_2D(int l, double x, double The)
{
    double a=-4.0*sFWHM(The), b=4.0*sFWHM(The);
    
    double epsrel=Int_eps_3D, epsabs=1.0e-300;
    
    double d[3]={x, The, (double)l};
    void *userdata=(void *)d;
    
    return Integrate_using_Patterson_adaptive(a, b, epsrel, epsabs, Boltzmann_DI_Kernel_II, userdata);
}

//==================================================================================================
double compute_SZ_distortion_Patterson_multiple_Kernel(double x, int l,
                                                       double The, double eps_Int)
{
    Int_eps_3D=eps_Int;
    
    return DI_Int_Kernel(l, x, The);
}

double compute_SZ_distortion_Patterson_multiple_Kernel_2D(double x, int l,
                                                          double The, double eps_Int)
{
    Int_eps_3D=eps_Int;
    
    return DI_Int_Kernel_2D(l, x, The);
}

//==================================================================================================
//
// 3D integration for multiple scattering terms carried out using Patterson scheme
//
// eps_Int : relative accuracy for numerical integration (lower than 10^-6 is hard to achieve)
//
// mode == "SEC"      --> second scattering correction of thermal SZ effect (x is irrelevant)
// mode == "SIG"      --> compute scattering cross section correction
// mode == "CMB"      --> compute scattering of CMB multipoles
//
//==================================================================================================
double compute_SZ_distortion_Patterson_multiple(double x, int l,
                                                double The,
                                                string mode,
                                                double eps_Int)
{
    Int_eps_3D=eps_Int;
    run_mode=mode;
    
    double norm_f=norm_df_RM_dTheta(The)/The;
    //TODO: This /The scaling is because I changed this function many many commits before
    double r=TWOPI*TWOPI*norm_f*xi_Int(x, The, l);
    
    return r;
}

//==================================================================================================
//
// asymptotic expansion of multiple scattering terms. 
//
//==================================================================================================

//==================================================================================================
// auxillary function
//==================================================================================================
void compute_Yl0_k(double x, const double Mptr[][25], vector<double> &Yk)
{
    int ymax=Yk.size()-1;
    
    vector<double> Pk(2*(ymax+2)+1), beta(2*(ymax+2)+2);
    
    for(int k=0; k<=2*(ymax+2); k++) Pk[k]=Pfunc(k, x);
    
    double exp_mx=exp(-x);
    double fac=(-x)/(1.0-exp_mx);
    double nPl=exp_mx/(1.0-exp_mx);
    
    for(int k=0; k<=ymax; k++)
    {
        int kmax=2*(k+2);
        beta[kmax+1]=0.0;
        
        for(int l=kmax; l>=1; l--) beta[l]=fac*(Mptr[k][l]*Pk[l]+beta[l+1]);
        
        Yk[k]=nPl*beta[1];
    }
    
    return;
}

void compute_Yl0_k(double x, int l, vector<double> &Yk)
{
    if(l>3 || Yk.size()>10) return;

    if(l==0) compute_Yl0_k(x, M0_ij, Yk);
    else if(l==1) compute_Yl0_k(x, M1_ij, Yk);
    else if(l==2) compute_Yl0_k(x, M2_ij, Yk);
    else if(l==3) compute_Yl0_k(x, M3_ij, Yk);
    
    return;
}

//==================================================================================================
//
// asymptotic expansion of multiple scattering terms according to Chluba et al. 2013. A maximum
// of 9 temperature correction terms, i.e., O(The^11) can be included (0<=Te_corr_order<=9)
//
//==================================================================================================
double compute_SZ_distortion_asym_multiple(double x, int l, double The, int Te_corr_order)
{
    if(l>3 || Te_corr_order<0) return 0;
    if(Te_corr_order>9) Te_corr_order=9;
        
    vector<double> Yk(Te_corr_order+1);
    
    if(l==0) compute_Yl0_k(x, M0_ij, Yk);
    else if(l==1) compute_Yl0_k(x, M1_ij, Yk);
    else if(l==2) compute_Yl0_k(x, M2_ij, Yk);
    else if(l==3) compute_Yl0_k(x, M3_ij, Yk);
    
    double r=The*Yk[Te_corr_order];
    for(int k=Te_corr_order-1; k>=0; k--) r=The*(r+Yk[k]);
    
    return The*r;
}

//==================================================================================================
//==================================================================================================
