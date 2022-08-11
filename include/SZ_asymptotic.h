//==================================================================================================
//
// This program allows computing the thermal SZ and kinematic effect. Up to 10th order temperature
// corrections are allowed. The basis functions are computed using recursion relations, based on 
// Eulerian numbers to determine the derivatives x^k d^k nPl(x)/ dx^k.
//
// The expressions for the thermal SZ effect are equivalent to those of Itoh et al., 1998 and 
// Shimon & Rephaeli, 2004, but to higher order kinematic corrections are computed according to 
// Chluba, Nagai, Sazonov & Nelson, 2012. 
//
//==================================================================================================
//
// Author: Jens Chluba  (CITA, University of Toronto)
//
// first implementation: April 2012
// last modification   : July  2012
//
//==================================================================================================
//  4st  Aug: added derivatives of basis functions in the CMB rest frame
// 10th July: added basis functions in the CMB rest frame

#ifndef SZ_ASYMPTOTIC_H
#define SZ_ASYMPTOTIC_H

using namespace std;

//==================================================================================================
//
// compute Dn using asymptotic expansion in cluster frame
//
// mode == "monopole"      --> only scattering of monopole without second order kinematic corr
// mode == "dipole"        --> only scattering of dipole     (first order kinematic correction)
// mode == "quadrupole"    --> only scattering of quadrupole (second order kinematic correction)
// mode == "monopole_corr" --> only scattering of second order kinematic correction to monopole
// mode == "all"           --> all terms added
// mode == "kin"           --> only kinematic terms
//
//==================================================================================================
double compute_SZ_distortion_asymptotic(double x, 
                                        double The, double betac, double muc, 
                                        int Torder, int betac_order, 
                                        string mode);

//==================================================================================================
//
// compute Dn using asymptotic expansion in CMB rest frame (added 10.07 by JC)
//
// mode == "monopole"      --> only monopole part without second order kinematic corr
// mode == "dipole"        --> only dipolar part     (first order kinematic correction)
// mode == "quadrupole"    --> only quadrupolar part (second order kinematic correction)
// mode == "monopole_corr" --> only second order kinematic correction to monopole part
// mode == "all"           --> all terms added
// mode == "kin"           --> only kinematic terms
//
//==================================================================================================
double compute_SZ_distortion_asymptotic_CMB(double x, 
                                            double The, double betac, double muc, 
                                            int Te_order, int betac_order,
                                            string mode);

//==================================================================================================
// Derivatives (The^k d^k_dThe /k!) (betapara^m d^m_dbetapara /m!) (beta2perp^l d^l_beta2perp /l!) S
// in the CMB frame for a resting observer. Maximal orders in The and betac are used to compute 
// the derivatives.
//
// constraints: dThe<=4; dbeta_para<=2; dbeta2_perp<=1;
//==================================================================================================
void Dcompute_SZ_distortion_asymptotic_CMB(double x, 
                                           int dThe, int dbeta_para, int dbeta2_perp,
                                           double The, double betac, double muc,
                                           vector<double> &dDn_dThe);

//==================================================================================================
//
// access basis functions (added 10.07 by JC)
//
//==================================================================================================
void compute_Y(double x, vector<double> &Y);
void compute_Ykin(double x, vector<double> &Ykin);
void compute_Dkin(double x, vector<double> &Dkin);
void compute_Qkin(double x, vector<double> &Qkin);
//
void compute_M_CMB(double x, vector<double> &MCMB);
void compute_D_CMB(double x, vector<double> &DCMB);
void compute_Q_CMB(double x, vector<double> &QCMB);

#endif

//==================================================================================================
//==================================================================================================
