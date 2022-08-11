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

#ifndef SZPACK_H
#define SZPACK_H

#include <cmath>
#include <vector>

#include "SZ_Integral.5D.h"
#include "SZ_Integral.3D.h"
#include "SZ_asymptotic.h"
#include "SZ_CNSN_basis.h"
#include "SZ_CNSN_basis.opt.h"
#include "SZ_Integral_multiple.h"

using namespace std;

const string SZpack_version="SZpack v1.1";

//==================================================================================================
// conversion factor x^3 Dn --> DI in MJy / sr
//==================================================================================================
const double T0_CMB=2.726;
const double Dn_DI_conversion=13.33914078*pow(T0_CMB, 3);

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
void use_Nozawa2006_convention();
void use_CNSN_convention();


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
                            double eps_Int=1.0e-4);

//--------------------------------------------------------------------------------------------------
// vector versions (output is returned in xo-vector)
//--------------------------------------------------------------------------------------------------
void compute_SZ_signal_5D(double *xo, int np, 
                          double Dtau, double Te, 
                          double betac, double muc, 
                          double betao, double muo, 
                          double eps_Int=1.0e-4);

void compute_SZ_signal_5D(vector<double> &xo, 
                          double Dtau, double Te, 
                          double betac, double muc, 
                          double betao, double muo, 
                          double eps_Int=1.0e-4);


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
                            double eps_Int=1.0e-4);

//--------------------------------------------------------------------------------------------------
// vector versions (output is returned in xo-vector)
//--------------------------------------------------------------------------------------------------
void compute_SZ_signal_3D(double *xo, int np, 
                          double Dtau, double Te, 
                          double betac, double muc, 
                          double betao, double muo, 
                          double eps_Int=1.0e-4);

void compute_SZ_signal_3D(vector<double> &xo, 
                          double Dtau, double Te, 
                          double betac, double muc, 
                          double betao, double muo, 
                          double eps_Int=1.0e-4);


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
                                    int Te_order, int betac_order);

//--------------------------------------------------------------------------------------------------
// vector versions (output is returned in xo-vector)
//--------------------------------------------------------------------------------------------------
void compute_SZ_signal_asymptotic(double *xo, int np, 
                                  double Dtau, double Te, 
                                  double betac, double muc, 
                                  double betao, double muo, 
                                  int Te_order, int betac_order);

void compute_SZ_signal_asymptotic(vector<double> &xo, 
                                  double Dtau, double Te, 
                                  double betac, double muc, 
                                  double betao, double muo, 
                                  int Te_order, int betac_order);


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
                                    int Te_order, int betac_order);

//--------------------------------------------------------------------------------------------------
// vector versions (output is returned in xo-vector)
//--------------------------------------------------------------------------------------------------
void compute_SZ_signal_CNSN_basis(double *xo, int np, 
                                  double Dtau, double Te, 
                                  double betac, double muc, 
                                  double betao, double muo, 
                                  int Te_order, int betac_order);

void compute_SZ_signal_CNSN_basis(vector<double> &xo, 
                                  double Dtau, double Te, 
                                  double betac, double muc, 
                                  double betao, double muo, 
                                  int Te_order, int betac_order);


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
                                        int accuracy_level);

//--------------------------------------------------------------------------------------------------
// vector versions (output is returned in xo-vector)
//--------------------------------------------------------------------------------------------------
void compute_SZ_signal_CNSN_basis_opt(double *xo, int np,
                                      double Dtau, double Te,
                                      double betac, double muc,
                                      double betao, double muo,
                                      int kmax, int betac_order,
                                      int accuracy_level);

void compute_SZ_signal_CNSN_basis_opt(vector<double> &xo,
                                      double Dtau, double Te,
                                      double betac, double muc,
                                      double betao, double muo,
                                      int kmax, int betac_order,
                                      int accuracy_level);

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
                               double betao, double muo);

//--------------------------------------------------------------------------------------------------
// vector versions (output is returned in xo-vector)
//--------------------------------------------------------------------------------------------------
void compute_SZ_signal_combo(double *xo, int np, 
                             double Dtau, double Te, 
                             double betac, double muc, 
                             double betao, double muo);

void compute_SZ_signal_combo(vector<double> &xo, 
                             double Dtau, double Te, 
                             double betac, double muc, 
                             double betao, double muo);


//==================================================================================================
// Derivatives (The^k d^k_dThe /k!) (d^m_dbetapara /m!) (d^l_beta2perp /l!) S(...)
// in the CMB frame for a resting observer. Maximal orders in The and betac are used to compute 
// the derivatives.
//
// constraints: dThe<=4; dbeta_para<=2; dbeta2_perp<=1;
//==================================================================================================
void Dcompute_SZ_signal_combo_CMB(double x, int k, int m, int l,
                                  double Dtau, double Te, double betac, double muc,
                                  vector<double> &dDn_dThe);


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
                                     double kappa, double betac2_perp);

//--------------------------------------------------------------------------------------------------
// vector versions (output is returned in xo-vector)
//--------------------------------------------------------------------------------------------------
void compute_SZ_signal_combo_means(double *xo, int np, 
                                   // mean parameters
                                   double tau, double TeSZ, double betac_para,
                                   // variances
                                   double omega, double sigma, 
                                   double kappa, double betac2_perp);

void compute_SZ_signal_combo_means(vector<double> &xo, 
                                   // mean parameters
                                   double tau, double TeSZ, double betac_para,
                                   // variances
                                   double omega, double sigma, 
                                   double kappa, double betac2_perp);


//==================================================================================================
// 
// Compute null of SZ signal using expansion around mean values of tau, TeSZ, and betac*muc. This  
// routine should reproduce the full numerical result for Te < 75keV, and beta < 0.01 with high 
// precision, assuming smooth cluster profile. (added 8th Aug, 2012, JC)
//
//==================================================================================================
double compute_null_of_SZ_signal(double tau, double TeSZ, double betac_para,
                                 double omega, double sigma, 
                                 double kappa, double betac2_perp);


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
                                        double kappa, double betac2_perp);

//--------------------------------------------------------------------------------------------------
// vector versions (output is returned in xo-vector)
//--------------------------------------------------------------------------------------------------
void compute_SZ_signal_combo_means_ex(double *xo, int np, 
                                      // mean parameters
                                      double tau, double TeSZ, double betac_para,
                                      // variances
                                      double omega[3], double sigma[3], 
                                      double kappa, double betac2_perp);

void compute_SZ_signal_combo_means_ex(vector<double> &xo, 
                                      // mean parameters
                                      double tau, double TeSZ, double betac_para,
                                      // variances
                                      double omega[3], double sigma[3], 
                                      double kappa, double betac2_perp);

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
                                        double TeSZ2, double TeSZ3, double TeSZ4);

//--------------------------------------------------------------------------------------------------
// vector versions (output is returned in xo-vector)
//--------------------------------------------------------------------------------------------------
void compute_SZ_signal_combo_means_yw(double *xo, int np,
                                      // mean parameters
                                      double y, double TeSZ,
                                      double TeSZ2, double TeSZ3, double TeSZ4);

void compute_SZ_signal_combo_means_yw(vector<double> &xo,
                                      // mean parameters
                                      double y, double TeSZ,
                                      double TeSZ2, double TeSZ3, double TeSZ4);

#endif

//==================================================================================================
//==================================================================================================
