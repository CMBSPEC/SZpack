//==================================================================================================
//
// basic SZpack python functions
// 
//==================================================================================================
//
// Author: Jens Chluba & Eric Switzer (CITA, University of Toronto)
//
// first implementation:  Aug 2012
// last modification   :  Aug 2012
//
//==================================================================================================

#ifndef SZPACK_PY_H
#define SZPACK_PY_H

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
                  double eps_Int=1.0e-4);

void compute_5d(double *xo, int np, 
                double Dtau, double Te, 
                double betac, double muc, 
                double betao, double muo, 
                double eps_Int=1.0e-4);

//--------------------------------------------------------------------------------------------------
//
// 3D integral
//
//--------------------------------------------------------------------------------------------------
double compute_3d(double xo,
                  double Dtau, double Te, 
                  double betac, double muc, 
                  double betao, double muo, 
                  double eps_Int=1.0e-4);

void compute_3d(double *xo, int np, 
                double Dtau, double Te, 
                double betac, double muc, 
                double betao, double muo, 
                double eps_Int=1.0e-4);

//--------------------------------------------------------------------------------------------------
//
// asymptotic expansion
//
//--------------------------------------------------------------------------------------------------
double compute_asym(double xo,
                    double Dtau, double Te, 
                    double betac, double muc, 
                    double betao, double muo, 
                    int Te_order, int betac_order);

void compute_asym(double *xo, int np, 
                  double Dtau, double Te, 
                  double betac, double muc, 
                  double betao, double muo, 
                  int Te_order, int betac_order);

//--------------------------------------------------------------------------------------------------
//
// basis functions according to CNSN2012 with extension of CSNN2012
//
//--------------------------------------------------------------------------------------------------
double compute_CNSN(double xo,
                    double Dtau, double Te, 
                    double betac, double muc, 
                    double betao, double muo, 
                    int Te_order, int betac_order);

void compute_CNSN(double *xo, int np, 
                  double Dtau, double Te, 
                  double betac, double muc, 
                  double betao, double muo, 
                  int Te_order, int betac_order);

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
                        int kmax, int betac_order, int accuracy_level);

void compute_CNSN_opt(double *xo, int np,
                      double Dtau, double Te,
                      double betac, double muc,
                      double betao, double muo,
                      int kmax, int betac_order, int accuracy_level);

//--------------------------------------------------------------------------------------------------
//
// combination of asymptotic expansion and CNSN-basis to cover a wide range of temperatures
//
//--------------------------------------------------------------------------------------------------
double compute_combo(double xo,
                     double Dtau, double Te, 
                     double betac, double muc, 
                     double betao, double muo);

void compute_combo(double *xo, int np,
                   double Dtau, double Te, 
                   double betac, double muc, 
                   double betao, double muo);

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
                           double kappa, double betac2_perp);

void compute_combo_means(double *xo, int np, 
                         // mean parameters
                         double tau, double TeSZ, double betac_para,
                         // variances
                         double omega1, double sigma,
                         double kappa, double betac2_perp);

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
                              double kappa, double betac2_perp);

void compute_combo_means_ex(double *xo, int np, 
                            // mean parameters
                            double tau, double TeSZ, double betac_para,
                            // variances
                            double *omega, int nomega, double *sigma, int nsigma,
                            double kappa, double betac2_perp);

//--------------------------------------------------------------------------------------------------
// 
// Cross-over frequency of the SZ signal
//
//--------------------------------------------------------------------------------------------------
double compute_null(double tau, double TeSZ, double betac_para,
                    double omega1, double sigma,
                    double kappa, double betac2_perp);

#endif

//==================================================================================================
//==================================================================================================
