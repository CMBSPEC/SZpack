//==================================================================================================
//
// SZpack functions
//
//==================================================================================================
//
// purpose: computation of the SZ signal according to Chluba, Nagai, Sazonov, Nelson, 2012
//          and Chluba, Switzer, Nagai, Nelson, 2012.
//
// comments: - computations are performed in the single-scattering approximation
//           - polarization effects are neglected
//           - the electron distribution function is assumed to be thermal
//==================================================================================================
//
// Author: Jens Chluba (CITA, University of Toronto and Johns Hopkins University)
//
// first implementation: May  2012
// last modification   : Aug  2017
//
//==================================================================================================
// 28th Aug 2017: added y-weighted moment method for temperature corrections
// 13th May 2015: improved performance of combo-means-methods && fixed bug for asymptotic derivative
// 20th  Dec: added optimized SZ signal function to minimize number of temperature terms
// 29th  Aug: added expansion of SZ signal around mean values of tau, TeSZ, and betac*muc with
//            higher moments included. 
//  8th  Aug: added function to compute SZ null.
//  6th  Aug: added expansion of SZ signal around mean values of tau, TeSZ, and betac*muc. 
//  4th  Aug: added derivatives of basis functions in the CMB rest frame
// 22th July: added combination of asymptotic expansion + CNSN basis functions. For Te < 2keV 
//            the asymptotic expansion is used while for 2keV < Te < 75keV the basis of CNSN is
//            applied. At 0.01 < x < 30 the precision should be similar to 0.001% for beta < 0.01. 
//            Also added simple sanity checks for approximation functions.
// 10th July: added functions to compute the SZ signal using temperature-velocity moments. These   
//            routines allow taking into account the detailed temperature and velocity structure   
//            of the cluster medium along the line-of-sight.
//--------------------------------------------------------------------------------------------------

//==================================================================================================
//
//  Unless stated otherwise the main parameters are:
//  
//  xo          : observer frame photon frequency xo == h nu / k T0 (T0 == todays CMB temperature)
//  Dtau        : line-of-sight scattering optical depth in the cluster frame
//  Te          : temperature of the (cluster frame) electron gas in keV
//  betac       : velocity beta=v/c of the cluster in the CMB frame
//  muc         : direction cosine of the cluster velocity with respect to the line-of-sight in  
//                the CMB rest frame
//  betao       : velocity beta=v/c of the observer in the CMB frame 
//                (from CMB dipole betao=1.241e-3)
//  muo         : direction cosine for the line-of-sight with respect to the observers velocity. 
//                The angle is measured in the observer frame
//
//==================================================================================================

#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <vector>

#include "physical_consts.h"
#include "routines.h"
#include "SZpack.h"

using namespace std;


//==================================================================================================
//
// As discussed by Chluba et el. 2012, Nozawa 1998 and 2006 used another convention for the optical
// depth variable.
// 
//      Dtau --> gammac (1-betac muc) Dtau 
//
// for an observer at rest in the CMB frame. To use this convention call 
// 'use_Nozawa2006_convention();' before excecuting the other routines. To reset to the convention 
// of Chluba et al. call 'use_CNSN_convention();', which is the default. These functions only have
// to be called once. Subsequent calls of routines will use the convention that was set last.
//
//==================================================================================================
bool CNSN2012_convention=1;

void use_Nozawa2006_convention(){ CNSN2012_convention=0; }
void use_CNSN_convention(){ CNSN2012_convention=1; }


//==================================================================================================
//
//  Full numerical integration of the 5-dimensional Compton collision term. All terms in betac are
//  included. In principle here one can easily add the scattering of primordial CMB anisotropies.
//
//  eps_Int : relative accuracy for numerical integration (lower than 10^-6 is hard to achieve)
//
//==================================================================================================
double compute_SZ_signal_5D(double xo, 
                            double Dtau, double Te, 
                            double betac, double muc, 
                            double betao, double muo, 
                            double eps_Int)
{
    double The=Te/const_me; // conversion keV --> The = k Te / me c^2
    double gammao=1.0/sqrt(1.0-betao*betao);
    double gammac=1.0/sqrt(1.0-betac*betac);
    
    // transformation of muc to mucc
    double mucc=(muc-betac)/(1.0-betac*muc);
        
    // transformation of xo to xc
    double xc=gammac*(1.0-betac*mucc)*gammao*xo*(1.0+betao*muo);
    
    // change to Nozawa 2006 convention for scattering optical depth
    if(!CNSN2012_convention) Dtau*=gammac*(1.0-betac*mucc);
    
    return Dtau*compute_SZ_distortion_Patterson_5D(xc, The, betac, mucc, eps_Int);
}

//--------------------------------------------------------------------------------------------------
// vector versions (output is returned in xo-vector)
//--------------------------------------------------------------------------------------------------
void compute_SZ_signal_5D(double *xo, int np, 
                          double Dtau, double Te, 
                          double betac, double muc, 
                          double betao, double muo, 
                          double eps_Int)
{
    for(int m=0; m<np; m++) 
        xo[m]=compute_SZ_signal_5D(xo[m], Dtau, Te, betac, muc, betao, muo, eps_Int);
    
    return;
}



void compute_SZ_signal_5D(vector<double> &xo, 
                          double Dtau, double Te, 
                          double betac, double muc, 
                          double betao, double muo, 
                          double eps_Int)
{
    compute_SZ_signal_5D(&xo[0], xo.size(), Dtau, Te, betac, muc, betao, muo, eps_Int);
    
    return;
}


//==================================================================================================
//
//  Full numerical integration after reduction of Compton collision term to 3 dimensions. Only terms
//  up to betac^2 are retained. Execution of this routine is significantly faster than for the full 
//  5D integral.
//
//  eps_Int : relative accuracy for numerical integration (lower than 10^-6 is hard to achieve)
//
//==================================================================================================
double compute_SZ_signal_3D(double xo, 
                            double Dtau, double Te, 
                            double betac, double muc, 
                            double betao, double muo, 
                            double eps_Int)
{
    double The=Te/const_me; // conversion keV --> The = k Te / me c^2
    double gammao=1.0/sqrt(1.0-betao*betao);
    double gammac=1.0/sqrt(1.0-betac*betac);
    
    // transformation of muc to mucc
    double mucc=(muc-betac)/(1.0-betac*muc);
    
    // transformation of xo to xc
    double xc=gammac*(1.0-betac*mucc)*gammao*xo*(1.0+betao*muo);
    
    // change to Nozawa 2006 convention for scattering optical depth
    if(!CNSN2012_convention) Dtau*=gammac*(1.0-betac*mucc);

    return Dtau*compute_SZ_distortion_Patterson_3D(xc, The, betac, mucc, "all", eps_Int);
}

//--------------------------------------------------------------------------------------------------
// vector versions (output is returned in xo-vector)
//--------------------------------------------------------------------------------------------------
void compute_SZ_signal_3D(double *xo, int np, 
                          double Dtau, double Te, 
                          double betac, double muc, 
                          double betao, double muo, 
                          double eps_Int)
{
    for(int m=0; m<np; m++) 
        xo[m]=compute_SZ_signal_3D(xo[m], Dtau, Te, betac, muc, betao, muo, eps_Int);
    
    return;
}



void compute_SZ_signal_3D(vector<double> &xo, 
                          double Dtau, double Te, 
                          double betac, double muc, 
                          double betao, double muo, 
                          double eps_Int)
{
    compute_SZ_signal_3D(&xo[0], xo.size(), Dtau, Te, betac, muc, betao, muo, eps_Int);
    
    return;
}


//==================================================================================================
//
//  Asymptotic expansion of the Compton collision integral similar to Itoh et al. but up to 
//  10th order in temperature. Kinematic terms are computed according to Chluba et al. 2012.
//
//  The additional parameters are:
//  
//  Te_order    (<=10): maximal order of temperature corrections. Te_order==0 means that the  
//                      normal thSZ formula is used.
//  betac_order (<=2) : maximal order of the kinematic corrections in the cluster frame. 
//                      betac_order == 0 means only the thSZ effect is considered.
//
//==================================================================================================
double compute_SZ_signal_asymptotic(double xo, 
                                    double Dtau, double Te, 
                                    double betac, double muc, 
                                    double betao, double muo, 
                                    int Te_order, int betac_order)
{
    //==============================================================================================
    // simple sanity checks
    //==============================================================================================
    if(Te>20.0)
    { 
        cerr << " compute_SZ_signal_asymptotic :: Temperature really (!) high "
             << "(Te = " << Te << " keV)" << endl;
    }
    
    if(Te_order>10){ Te_order=10; }
    if(betac_order>2){ betac_order=2; }
    
    //==============================================================================================
    // computations
    //==============================================================================================
    double The=Te/const_me; // conversion keV --> The = k Te / me c^2
    double gammao=1.0/sqrt(1.0-betao*betao);
    double gammac=1.0/sqrt(1.0-betac*betac);
    
    // transformation of muc to mucc
    double mucc=(muc-betac)/(1.0-betac*muc);
    
    // transformation of xo to xc
    double xc=gammac*(1.0-betac*mucc)*gammao*xo*(1.0+betao*muo);
    
    // change to Nozawa 2006 convention for scattering optical depth
    if(!CNSN2012_convention) Dtau*=gammac*(1.0-betac*mucc);

    return Dtau*compute_SZ_distortion_asymptotic(xc, The, betac, mucc, Te_order, betac_order,"all");
}

//--------------------------------------------------------------------------------------------------
// vector versions (output is returned in xo-vector)
//--------------------------------------------------------------------------------------------------
void compute_SZ_signal_asymptotic(double *xo, int np, 
                                  double Dtau, double Te, 
                                  double betac, double muc, 
                                  double betao, double muo, 
                                  int Te_order, int betac_order)
{
    for(int m=0; m<np; m++) 
        xo[m]=compute_SZ_signal_asymptotic(xo[m], Dtau, Te, 
                                           betac, muc, betao, muo, 
                                           Te_order, betac_order);
    
    return;
}



void compute_SZ_signal_asymptotic(vector<double> &xo, 
                                  double Dtau, double Te, 
                                  double betac, double muc, 
                                  double betao, double muo, 
                                  int Te_order, int betac_order)
{
    compute_SZ_signal_asymptotic(&xo[0], xo.size(), Dtau, Te,                                           
                                 betac, muc, betao, muo, 
                                 Te_order, betac_order);
    
    return;
}


//==================================================================================================
//
//  Improved expansion of the Compton collision integral following the approach of 
//  Chluba, Nagai, Sazonov, Nelson, 2012 
//
//  The additional parameters are:
//  
//  Te_order    (<=20): maximal order of temperature corrections. Te_order<=4 is not recommended.
//  betac_order (<=2) : maximal order of the kinematic corrections in the cluster frame. 
//                      betac_order == 0 means only the thSZ effect is considered.
//
//==================================================================================================
// 22th July, 2012: The allowed temperature range has been increased to 2keV < Te < 75keV 
//                  (see SZ_CNSN_basis.cpp for details)
//--------------------------------------------------------------------------------------------------
double compute_SZ_signal_CNSN_basis(double xo, 
                                    double Dtau, double Te, 
                                    double betac, double muc, 
                                    double betao, double muo, 
                                    int Te_order, int betac_order)
{
    //==============================================================================================
    // simple sanity checks
    //==============================================================================================
    if(xo<0.01 || xo>50.0)
    { 
        cerr << " compute_SZ_signal_CNSN_basis :: You are outside of grid \n" << endl;
        exit(0);
    }
    
    if(Te<2.0 || Te>75.0)
    { 
        cerr << " compute_SZ_signal_CNSN_basis :: Temperature too high/low \n" << endl;
        exit(0);
    }

    if(Te_order>20){ Te_order=20; }
    if(betac_order>2){ betac_order=2; }
    
    //==============================================================================================
    // computations
    //==============================================================================================
    double The=Te/const_me; // conversion keV --> The = k Te / me c^2
    double gammao=1.0/sqrt(1.0-betao*betao);
    double gammac=1.0/sqrt(1.0-betac*betac);
    
    // transformation of muc to mucc
    double mucc=(muc-betac)/(1.0-betac*muc);
    
    // transformation of xo to xc
    double xc=gammac*(1.0-betac*mucc)*gammao*xo*(1.0+betao*muo);
    
    // change to Nozawa 2006 convention for scattering optical depth
    if(!CNSN2012_convention) Dtau*=gammac*(1.0-betac*mucc);

    return Dtau*compute_SZ_distortion_CNSN_basis(xc, The, betac, mucc, Te_order, betac_order,"all");
}

//--------------------------------------------------------------------------------------------------
// vector versions (output is returned in xo-vector)
//--------------------------------------------------------------------------------------------------
void compute_SZ_signal_CNSN_basis(double *xo, int np, 
                                  double Dtau, double Te, 
                                  double betac, double muc, 
                                  double betao, double muo, 
                                  int Te_order, int betac_order)
{
    for(int m=0; m<np; m++) 
        xo[m]=compute_SZ_signal_CNSN_basis(xo[m], Dtau, Te, 
                                           betac, muc, betao, muo, 
                                           Te_order, betac_order);
    
    return;
}



void compute_SZ_signal_CNSN_basis(vector<double> &xo, 
                                  double Dtau, double Te, 
                                  double betac, double muc, 
                                  double betao, double muo, 
                                  int Te_order, int betac_order)
{
    compute_SZ_signal_CNSN_basis(&xo[0], xo.size(), Dtau, Te,                                           
                                 betac, muc, betao, muo, 
                                 Te_order, betac_order);
    
    return;
}


//==================================================================================================
//
// Similar to compute_SZ_signal_CNSN_basis-function defined above but here kmax and accuracy_level
// can be set according to the definitions of CSNN2012. Note that here the SZ signal is computed 
// for an observer in the CMB rest frame, while the SZ signal in the observer frame is obtained by
// Lorentz-transformation.
//
// kmax: depends on accuracy_level (see Table 1 in CSNN2012)
// accuracy_level: 0, 1, 2, 3
//
// Dependening on kmax and accuracy_level the maximal allowed temperature varies slightly but is
// usually larger than 60 keV (Table 1 in CSNN2012).
//
// (added 20.12.2012)
//==================================================================================================
double compute_SZ_signal_CNSN_basis_opt(double xo,
                                        double Dtau, double Te,
                                        double betac, double muc,
                                        double betao, double muo,
                                        int kmax, int betac_order,
                                        int accuracy_level)
{
    //==============================================================================================
    // simple sanity checks
    //==============================================================================================
    if(xo<0.01 || xo>50.0)
    {
        cerr << " compute_SZ_signal_CNSN_basis_opt :: You are outside of grid \n" << endl;
        exit(0);
    }
    
    //==============================================================================================
    // computations
    //==============================================================================================
    double The=Te/const_me; // conversion keV --> The = k Te / me c^2
    double gammao=1.0/sqrt(1.0-betao*betao);
    
    // transformation of xo to xc
    double xc=gammao*xo*(1.0+betao*muo);
    
    return Dtau*compute_SZ_distortion_CNSN_basis_opt(xc, The, betac, muc, kmax,
                                                     betac_order, "all", accuracy_level);
}

//--------------------------------------------------------------------------------------------------
// vector versions (output is returned in xo-vector)
//--------------------------------------------------------------------------------------------------
void compute_SZ_signal_CNSN_basis_opt(double *xo, int np,
                                      double Dtau, double Te,
                                      double betac, double muc,
                                      double betao, double muo,
                                      int kmax, int betac_order,
                                      int accuracy_level)
{
    for(int m=0; m<np; m++)
        xo[m]=compute_SZ_signal_CNSN_basis_opt(xo[m], Dtau, Te,
                                               betac, muc, betao, muo,
                                               kmax, betac_order, accuracy_level);
    
    return;
}



void compute_SZ_signal_CNSN_basis_opt(vector<double> &xo,
                                      double Dtau, double Te,
                                      double betac, double muc,
                                      double betao, double muo,
                                      int kmax, int betac_order,
                                      int accuracy_level)
{
    compute_SZ_signal_CNSN_basis_opt(&xo[0], xo.size(), Dtau, Te,
                                     betac, muc, betao, muo,
                                     kmax, betac_order, accuracy_level);
    
    return;
}


//==================================================================================================
// 
// Combination of asymptotic expansion and CNSN basis functions. This routine should reproduce the 
// full numerical result for 0.01 < x < 30, Te < 75keV, and beta < 0.01 with precision similar
// to 0.001%. (added 22th July, 2012, by JC)
//
//==================================================================================================
double compute_SZ_signal_combo(double xo, 
                               double Dtau, double Te, 
                               double betac, double muc, 
                               double betao, double muo)
{
    if(Te<2.0) return compute_SZ_signal_asymptotic(xo, Dtau, Te, betac, muc, betao, muo, 10, 2); 
    return compute_SZ_signal_CNSN_basis(xo, Dtau, Te, betac, muc, betao, muo, 20, 2); 
}

//--------------------------------------------------------------------------------------------------
// vector versions (output is returned in xo-vector)
//--------------------------------------------------------------------------------------------------
void compute_SZ_signal_combo(double *xo, int np, 
                             double Dtau, double Te, 
                             double betac, double muc, 
                             double betao, double muo)
{
    for(int m=0; m<np; m++) 
        xo[m]=compute_SZ_signal_combo(xo[m], Dtau, Te, betac, muc, betao, muo);

    return;
}
    


void compute_SZ_signal_combo(vector<double> &xo, 
                             double Dtau, double Te, 
                             double betac, double muc, 
                             double betao, double muo)
{
    compute_SZ_signal_combo(&xo[0], xo.size(), Dtau, Te, betac, muc, betao, muo);

    return;
}


//==================================================================================================
// Derivatives (The^k d^k_dThe /k!) (d^m_dbetapara /m!) (d^l_beta2perp /l!) S(...)
// in the CMB frame for a resting observer. Maximal orders in The and betac are used to compute 
// the derivatives.
//
// constraints: dThe<=4; dbeta_para<=2; dbeta2_perp<=1;
//==================================================================================================
void Dcompute_SZ_signal_combo_CMB(double x, int k, int m, int l,
                                  double Dtau, double Te, double betac, double muc,
                                  vector<double> &dDn_dThe)
{
    if(Te<2.0) Dcompute_SZ_distortion_asymptotic_CMB(x, k, m, l, Te/const_me, betac, muc, dDn_dThe); 
    else
    {
        //==========================================================================================
        // simple sanity checks
        //==========================================================================================
        if(x<0.01 || x>50.0)
        { 
            cerr << " Dcompute_SZ_signal_combo_CMB :: You are outside of grid \n" << endl;
            exit(0);
        }
        
        if(Te>75.0)
        { 
            cerr << " Dcompute_SZ_signal_combo_CMB :: Temperature too high/low \n" << endl;
            exit(0);
        }

        Dcompute_SZ_distortion_CNSN_basis_CMB(x, k, m, l, Te/const_me, betac, muc, dDn_dThe); 
    }
    
    for(int g=0; g<=k; g++) dDn_dThe[g]*=Dtau;
    
    return;
}


//==================================================================================================
// 
// Expansion of SZ signal around mean values of tau, TeSZ, and betac*muc. This routine should 
// reproduce the full numerical result for 0.01 < x < 30, Te < 75keV, and beta < 0.01 with 
// precision similar to 0.001%, assuming smooth cluster profile. (added 6th Aug, 2012, JC)
//
//==================================================================================================
double compute_SZ_signal_combo_means(double xo, 
                                     // mean parameters
                                     double tau, double TeSZ, double betac_para,
                                     // variances
                                     double omega, double sigma, 
                                     double kappa, double betac2_perp)
{
    vector<double> dDn_dThe;
    
    // mean signal + temperature dispersion term
    Dcompute_SZ_signal_combo_CMB(xo, 2, 0, 0, tau, TeSZ, betac_para, 1.0, dDn_dThe);
    double r=dDn_dThe[0]+dDn_dThe[2]*omega;
    
    // velocity - temperature cross term
    if(sigma!=0.0)
    {
        Dcompute_SZ_signal_combo_CMB(xo, 1, 1, 0, tau, TeSZ, betac_para, 1.0, dDn_dThe);
        r+=dDn_dThe[1]*sigma;
    }
    
    // betac_parallel dispersion term
    if(kappa!=0.0)
    {
        Dcompute_SZ_signal_combo_CMB(xo, 0, 2, 0, tau, TeSZ, betac_para, 1.0, dDn_dThe);
        r+=dDn_dThe[0]*kappa;
    }

    // betac_perp dispersion term
    if(betac2_perp!=0.0)
    {
        Dcompute_SZ_signal_combo_CMB(xo, 0, 0, 1, tau, TeSZ, betac_para, 1.0, dDn_dThe);
        r+=dDn_dThe[0]*betac2_perp;
    }
    
    return r;
}

//--------------------------------------------------------------------------------------------------
// vector versions (output is returned in xo-vector)
//--------------------------------------------------------------------------------------------------
void compute_SZ_signal_combo_means(double *xo, int np, 
                                   // mean parameters
                                   double tau, double TeSZ, double betac_para,
                                   // variances
                                   double omega, double sigma, 
                                   double kappa, double betac2_perp)
{
    for(int m=0; m<np; m++) 
        xo[m]=compute_SZ_signal_combo_means(xo[m], tau, TeSZ, betac_para, 
                                            omega, sigma, kappa, betac2_perp);
    
    return;
}


void compute_SZ_signal_combo_means(vector<double> &xo, 
                                   // mean parameters
                                   double tau, double TeSZ, double betac_para,
                                   // variances
                                   double omega, double sigma, 
                                   double kappa, double betac2_perp)
{
    compute_SZ_signal_combo_means(&xo[0], xo.size(), tau, TeSZ, betac_para, 
                                  omega, sigma, kappa, betac2_perp);
    
    return;
}


//==================================================================================================
// 
// Compute null of SZ signal using expansion around mean values of tau, TeSZ, and betac*muc. This  
// routine should reproduce the full numerical result for Te < 75keV, and beta < 0.01 with high 
// precision, assuming smooth cluster profile. (added 8th Aug, 2012, JC)
//
//==================================================================================================
double func_SZ_null(double *x, void *p)
{
    double *d=(double *)p;
    return compute_SZ_signal_combo_means(*x, d[0], d[1], d[2], d[3], d[4], d[5], d[6]);
}

double compute_null_of_SZ_signal(double tau, double TeSZ, double betac_para,
                                 double omega, double sigma, 
                                 double kappa, double betac2_perp)
{
    double d[]={tau, TeSZ, betac_para, omega, sigma, kappa, betac2_perp};
    void *p=(void *)d;
    return find_root_brent(func_SZ_null, p, 0.1, 10.0, 1.0e-5);
}


//==================================================================================================
//
// Extended versions of expansions around average values. Up to The^4 terms can be included for
// the thSZ signal and up to The^3 for first order cross terms (temperature corrections to higher 
// order velocity terms can in principle be activated but have been omitted at this point). 
// The parameters are:
// 
// omega[0..2] == omega^(1..3)         [i.e., O(The), O(The^2), O(The^3) & O(The^4)]
// sigma[0..2] == sigma^(1..3, 1)      [i.e., O(betac The), O(betac The^2), O(betac The^3)]
// kappa       == kappa^(1)            [i.e., O(betac^2 muc^2)]
// betac2_perp == <beta_c^2 (1-muc^2)> 
//
// according to definitions of CSNN 2012. Values that are not need should be set to zero.
// (added 29th Aug, 2012, JC)
//
//==================================================================================================
double compute_SZ_signal_combo_means_ex(double xo, 
                                        // mean parameters
                                        double tau, double TeSZ, double betac_para,
                                        // variances
                                        double omega[3], double sigma[3], 
                                        double kappa, double betac2_perp)
{
    vector<double> dDn_dThe;
    
    // mean signal + temperature dispersion term
    Dcompute_SZ_signal_combo_CMB(xo, 4, 0, 0, tau, TeSZ, betac_para, 1.0, dDn_dThe);
    double r=dDn_dThe[0]+dDn_dThe[2]*omega[0]+dDn_dThe[3]*omega[1]+dDn_dThe[4]*omega[2];
    
    // velocity - temperature cross term
    if(sigma[0]!=0.0 || sigma[1]!=0.0 || sigma[2]!=0.0)
    {
        Dcompute_SZ_signal_combo_CMB(xo, 3, 1, 0, tau, TeSZ, betac_para, 1.0, dDn_dThe);
        r+=dDn_dThe[1]*sigma[0]+dDn_dThe[2]*sigma[1]+dDn_dThe[3]*sigma[2];
    }
    
    // betac_parallel dispersion term
    if(kappa!=0.0)
    {
        Dcompute_SZ_signal_combo_CMB(xo, 0, 2, 0, tau, TeSZ, betac_para, 1.0, dDn_dThe);
        r+=dDn_dThe[0]*kappa;
    }
    
    // betac_perp dispersion term
    if(betac2_perp!=0.0)
    {
        Dcompute_SZ_signal_combo_CMB(xo, 0, 0, 1, tau, TeSZ, betac_para, 1.0, dDn_dThe);
        r+=dDn_dThe[0]*betac2_perp;
    }
    
    return r;    
}

//--------------------------------------------------------------------------------------------------
// vector versions (output is returned in xo-vector)
//--------------------------------------------------------------------------------------------------
void compute_SZ_signal_combo_means_ex(double *xo, int np, 
                                      // mean parameters
                                      double tau, double TeSZ, double betac_para,
                                      // variances
                                      double omega[3], double sigma[3], 
                                      double kappa, double betac2_perp)
{
    for(int m=0; m<np; m++) 
        xo[m]=compute_SZ_signal_combo_means_ex(xo[m], tau, TeSZ, betac_para, 
                                               omega, sigma, kappa, betac2_perp);
    
    return;
}

void compute_SZ_signal_combo_means_ex(vector<double> &xo, 
                                      // mean parameters
                                      double tau, double TeSZ, double betac_para,
                                      // variances
                                      double omega[3], double sigma[3], 
                                      double kappa, double betac2_perp)
{
    compute_SZ_signal_combo_means_ex(&xo[0], xo.size(), tau, TeSZ, betac_para, 
                                     omega, sigma, kappa, betac2_perp);
    
    return;
}    

//==================================================================================================
//
// y-weighted moment expansion for up to O(The^4) for thSZ only. Using the definition
// dy = (kTe / mec^2) dtau in the expressions below, we have the parameters
//
// y      : y     = int dy               [dimensionless]
// TeSZ   : TeSZ  = y^-1 int Te dy       [keV  ]
// TeSZ2  : TeSZ2 = y^-1 int Te^2 dy     [keV^2]
// TeSZ3  : TeSZ3 = y^-1 int Te^3 dy     [keV^3]
// TeSZ4  : TeSZ4 = y^-1 int Te^4 dy     [keV^4]
//
// The integrals have to be taken along the considered line-of-sight.
//
// [added 28.08.2017 JC]
//==================================================================================================
double compute_SZ_signal_combo_means_yw(double xo,
                                        // mean parameters
                                        double y, double TeSZ,
                                        double TeSZ2, double TeSZ3, double TeSZ4)
{
    // convert moments to omega_y
    double The=TeSZ/const_me;
    double tau=y/The;
    double oy2=(TeSZ2==0           ? 0 : TeSZ2/pow(TeSZ, 2) - 1.0);
    double oy3=(TeSZ3==0 || oy2==0 ? 0 : TeSZ3/pow(TeSZ, 3) - 3.0*oy2 - 1.0);
    double oy4=(TeSZ4==0 || oy3==0 ? 0 : TeSZ4/pow(TeSZ, 4) - 4.0*oy3 - 6.0*oy2 - 1.0);

    vector<double> dDn_dThe;
    
    // mean signal + temperature dispersion term
    Dcompute_SZ_signal_combo_CMB(xo, 4, 0, 0, tau, TeSZ, 0.0, 1.0, dDn_dThe);
    
    double H0= dDn_dThe[0];
    double H2=(dDn_dThe[2]-dDn_dThe[1]+dDn_dThe[0]);
    double H3=(dDn_dThe[3]-dDn_dThe[2]+dDn_dThe[1]-dDn_dThe[0]);
    double H4=(dDn_dThe[4]-dDn_dThe[3]+dDn_dThe[2]-dDn_dThe[1]+dDn_dThe[0]);
    
    double r=H0+H2*oy2+H3*oy3+H4*oy4;
    
    return r;
}

//--------------------------------------------------------------------------------------------------
// vector versions (output is returned in xo-vector)
//--------------------------------------------------------------------------------------------------
void compute_SZ_signal_combo_means_yw(double *xo, int np,
                                      // mean parameters
                                      double y, double TeSZ,
                                      double TeSZ2, double TeSZ3, double TeSZ4)
{
    for(int m=0; m<np; m++)
        xo[m]=compute_SZ_signal_combo_means_yw(xo[m], y, TeSZ, TeSZ2, TeSZ3, TeSZ4);
    
    return;
}


void compute_SZ_signal_combo_means_yw(vector<double> &xo,
                                      // mean parameters
                                      double y, double TeSZ,
                                      double TeSZ2, double TeSZ3, double TeSZ4)
{
    compute_SZ_signal_combo_means_yw(&xo[0], xo.size(), y, TeSZ, TeSZ2, TeSZ3, TeSZ4);
    
    return;
}

//==================================================================================================
//==================================================================================================
