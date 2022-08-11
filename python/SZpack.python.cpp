//==================================================================================================
//
// basic SZpack python functions
// 
//==================================================================================================
//
// Author: Jens Chluba & Eric Switzer (CITA, University of Toronto)
//
// first implementation:  Aug 2012
// last modification   :  Feb 2013
//
//==================================================================================================

#include <iostream>
#include <cstdlib>

#include "SZpack.h"
#include "SZpack.python.h"

using namespace std;

//--------------------------------------------------------------------------------------------------
//
// Unless stated otherwise, the functions below all return Dn(x) [distortion in terms of photon
// occupation number] in the observer frame [defined by betao and muo]. The spectral intensity
// is DI = 2 h nu^3 / c^2 Dn = 2 h/ c^2 (kTgamma/h)^3 x^3 Dn(x), where the factor is
// 2 h/ c^2 (kTgamma/h)^3 ~ 270 MJy / sr. For more details about the parameters see 'SZpack.h'.
//
//--------------------------------------------------------------------------------------------------

//--------------------------------------------------------------------------------------------------
//
// 5D integral
//
//--------------------------------------------------------------------------------------------------
double compute_5d(double xo,
                  double Dtau, double Te, 
                  double betac, double muc, 
                  double betao, double muo, 
                  double eps_Int)
{ 
    return compute_SZ_signal_5D(xo, Dtau, Te, betac, muc, betao, muo, eps_Int); 
}


void compute_5d(double *xo, int np, 
                  double Dtau, double Te, 
                  double betac, double muc, 
                  double betao, double muo, 
                  double eps_Int)
{ 
    compute_SZ_signal_5D(xo, np, Dtau, Te, betac, muc, betao, muo, eps_Int); 
    return;
}

//--------------------------------------------------------------------------------------------------
//
// 3D integral
//
//--------------------------------------------------------------------------------------------------
double compute_3d(double xo,
                  double Dtau, double Te, 
                  double betac, double muc, 
                  double betao, double muo, 
                  double eps_Int)
{ 
    return compute_SZ_signal_3D(xo, Dtau, Te, betac, muc, betao, muo, eps_Int); 
}


void compute_3d(double *xo, int np, 
                  double Dtau, double Te, 
                  double betac, double muc, 
                  double betao, double muo, 
                  double eps_Int)
{ 
    compute_SZ_signal_3D(xo, np, Dtau, Te, betac, muc, betao, muo, eps_Int); 
    return;
}

//--------------------------------------------------------------------------------------------------
//
// asymptotic expansion
//
//--------------------------------------------------------------------------------------------------
double compute_asym(double xo,
                    double Dtau, double Te, 
                    double betac, double muc, 
                    double betao, double muo, 
                    int Te_order, int betac_order)
{ 
    return compute_SZ_signal_asymptotic(xo, Dtau, Te, betac, muc, betao, muo, 
                                        Te_order, betac_order); 
}


void compute_asym(double *xo, int np, 
                  double Dtau, double Te, 
                  double betac, double muc, 
                  double betao, double muo, 
                  int Te_order, int betac_order)
{ 
    compute_SZ_signal_asymptotic(xo, np, Dtau, Te, betac, muc, betao, muo, 
                                 Te_order, betac_order); 
    return;
}


//--------------------------------------------------------------------------------------------------
//
// basis functions according to CNSN2012 with extension of CSNN2012
//
//--------------------------------------------------------------------------------------------------
double compute_CNSN(double xo,
                    double Dtau, double Te, 
                    double betac, double muc, 
                    double betao, double muo, 
                    int Te_order, int betac_order)
{
    return compute_SZ_signal_CNSN_basis(xo, Dtau, Te, betac, muc, betao, muo, 
                                        Te_order, betac_order);     
}

void compute_CNSN(double *xo, int np, 
                  double Dtau, double Te, 
                  double betac, double muc, 
                  double betao, double muo, 
                  int Te_order, int betac_order)
{ 
    compute_SZ_signal_CNSN_basis(xo, np, Dtau, Te, betac, muc, betao, muo, 
                                 Te_order, betac_order); 
    return;
}


//--------------------------------------------------------------------------------------------------
//
// basis functions according to CNSN2012 with extension of CSNN2012 but with optimal settings
// in order to minimize the number of terms needed.
//
//--------------------------------------------------------------------------------------------------
double compute_CNSN_opt(double xo,
                        double Dtau, double Te,
                        double betac, double muc,
                        double betao, double muo,
                        int kmax, int betac_order, int accuracy_level)
{
    return compute_SZ_signal_CNSN_basis_opt(xo, Dtau, Te, betac, muc, betao, muo,
                                            kmax, betac_order, accuracy_level);
}

void compute_CNSN_opt(double *xo, int np,
                      double Dtau, double Te,
                      double betac, double muc,
                      double betao, double muo,
                      int kmax, int betac_order, int accuracy_level)
{
    compute_SZ_signal_CNSN_basis_opt(xo, np, Dtau, Te, betac, muc, betao, muo,
                                     kmax, betac_order, accuracy_level);
    return;
}

//--------------------------------------------------------------------------------------------------
//
// combination of asymptotic expansion and CNSN-basis to cover a wide range of temperatures
//
//--------------------------------------------------------------------------------------------------
double compute_combo(double xo, 
                     double Dtau, double Te, 
                     double betac, double muc, 
                     double betao, double muo)
{
    return compute_SZ_signal_combo(xo, Dtau, Te, betac, muc, betao, muo);     
}


void compute_combo(double *xo, int np, 
                   double Dtau, double Te, 
                   double betac, double muc, 
                   double betao, double muo)
{ 
    compute_SZ_signal_combo(xo, np, Dtau, Te, betac, muc, betao, muo); 
    return;
}


//--------------------------------------------------------------------------------------------------
// 
// These functions are expansion of the mean SZ signal around tau, TeSZ, betac_para. Again Dn is 
// returned. The variances omega1, sigma, kappa, betac2_perp are all defined in CSNN 2012.
//
//--------------------------------------------------------------------------------------------------
double compute_combo_means(double xo, 
                           // mean parameters
                           double tau, double TeSZ, double betac_para,
                           // variances
                           double omega1, double sigma,
                           double kappa, double betac2_perp)
{
    return compute_SZ_signal_combo_means(xo, tau, TeSZ, betac_para, 
                                         omega1, sigma, kappa, betac2_perp);
}

void compute_combo_means(double *xo, int np, 
                         // mean parameters
                         double tau, double TeSZ, double betac_para,
                         // variances
                         double omega1, double sigma,
                         double kappa, double betac2_perp)
{ 
    compute_SZ_signal_combo_means(xo, np, tau, TeSZ, betac_para, 
                                  omega1, sigma, kappa, betac2_perp);
    return;
}


//--------------------------------------------------------------------------------------------------
// 
// Extended expansions of the mean SZ signal around tau, TeSZ, betac_para. Dn is
// returned. The variances omega[3], sigma[3], kappa, betac2_perp are all defined in CSNN 2012. See
// SZpack.h for more information. Here nomega and nsigma must be 3.
//
//--------------------------------------------------------------------------------------------------
double compute_combo_means_ex(double xo, 
                              // mean parameters
                              double tau, double TeSZ, double betac_para,
                              // variances
                              double *omega, int nomega, double *sigma, int nsigma,
                              double kappa, double betac2_perp)
{
    if(nomega!=3 || nsigma!=3) 
    {
        cerr << " compute_combo_means_ex: check provided number of omega/sigma moments " << endl; 
        exit(0);
    }

    double og[3]={omega[0], omega[1], omega[2]};
    double sg[3]={sigma[0], sigma[1], sigma[2]};

    return compute_SZ_signal_combo_means_ex(xo, tau, TeSZ, betac_para, 
                                            og, sg, kappa, betac2_perp);     
}

void compute_combo_means_ex(double *xo, int np, 
                            // mean parameters
                            double tau, double TeSZ, double betac_para,
                            // variances
                            double *omega, int nomega, double *sigma, int nsigma,
                            double kappa, double betac2_perp)
{ 
    if(nomega!=3 || nsigma!=3) 
    {
        cerr << " compute_combo_means_ex: check provided number of omega/sigma moments " << endl; 
        exit(0);
    }
    
    double og[3]={omega[0], omega[1], omega[2]};
    double sg[3]={sigma[0], sigma[1], sigma[2]};
    
    compute_SZ_signal_combo_means_ex(xo, np, tau, TeSZ, betac_para, 
                                     og, sg, kappa, betac2_perp); 
    return;
}

//--------------------------------------------------------------------------------------------------
//
// y-weighted moment expansion
//
//--------------------------------------------------------------------------------------------------
double compute_combo_means_yw(double xo,
                              // mean parameters
                              double tau, double TeSZ,
                              double TeSZ2, double TeSZ3, double TeSZ4)
{
    return compute_SZ_signal_combo_means_yw(xo, tau, TeSZ, TeSZ2, TeSZ3, TeSZ4);
}

void compute_combo_means_yw(double *xo, int np,
                            // mean parameters
                            double tau, double TeSZ,
                            double TeSZ2, double TeSZ3, double TeSZ4)
{
    compute_SZ_signal_combo_means_yw(xo, np, tau, TeSZ, TeSZ2, TeSZ3, TeSZ4);
    return;
}

//--------------------------------------------------------------------------------------------------
// 
// Cross-over frequency of the SZ signal
//
//--------------------------------------------------------------------------------------------------
double compute_null(double tau, double TeSZ, double betac_para,
                    double omega1, double sigma,
                    double kappa, double betac2_perp)
{
    return compute_null_of_SZ_signal(tau, TeSZ, betac_para, 
                                     omega1, sigma, kappa, betac2_perp);
}


//==================================================================================================
//==================================================================================================
