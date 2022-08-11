//==================================================================================================
//
// This program allows computing the thermal and kinematic SZ effect using the basis functions of 
// Chluba, Nagai, Sazonov & Nelson, 2012. Up to 20th order temperature corrections are allowed. 
// Kinematic corrections can be included up to second order in the clusters peculiar velocity. 
// For 0.01 < x < 30, 2keV < Te < 75keV and beta < 0.01 the obtained SZ signal should be accurate 
// at the level of 0.001% when 20 temperature order are included. This precision is acheived using 
// three temperature pivots, The=0.01, 0.03, and 0.1. For Te < 2keV the asymptotic expansions 
// should be used.
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
// 22th July: added low and high temperature expansions. Now 2keV < Te < 75keV is covered
// 21th July: added headers with required data; this avoids time taken for reading the files
// 10th July: added basis functions in the CMB rest frame
//  8th July: changed definition of S^kin; temperature-independent terms are fully canceled

#ifndef SZ_CNSN_BASIS_H
#define SZ_CNSN_BASIS_H

using namespace std;

//==================================================================================================
//
// compute Dn using improved expansion in cluster frame (Chluba, Nagai, Sazonov, Nelson, 2012)
//
// mode == "monopole"      --> only scattering of monopole without second order kinematic corr
// mode == "dipole"        --> only scattering of dipole     (first order kinematic correction)
// mode == "quadrupole"    --> only scattering of quadrupole (second order kinematic correction)
// mode == "monopole_corr" --> only scattering of second order kinematic correction to monopole
// mode == "all"           --> all terms added
// mode == "kin"           --> only kinematic terms
//
//==================================================================================================
double compute_SZ_distortion_CNSN_basis(double x, 
                                        double The, double betac, double muc, 
                                        int Te_order, int betac_order,
                                        string mode);

//==================================================================================================
//
// compute Dn using expansion of CNSN2012 but in the CMB frame (added 10.07 by JC)
//
// mode == "monopole"      --> only scattering of monopole without second order kinematic corr
// mode == "dipole"        --> only scattering of dipole     (first order kinematic correction)
// mode == "quadrupole"    --> only scattering of quadrupole (second order kinematic correction)
// mode == "monopole_corr" --> only scattering of second order kinematic correction to monopole
// mode == "all"           --> all terms added
// mode == "kin"           --> only kinematic terms
//
// Note that since only terms up to O(betac^2) are included, the results obtained with this routine  
// are slightly different from those obtained by Lorentz transformation of the cluster-frame signal 
// into the CMB frame.
//
//==================================================================================================
double compute_SZ_distortion_CNSN_basis_CMB(double x, 
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
void Dcompute_SZ_distortion_CNSN_basis_CMB(double x, 
                                           int dThe, int dbeta_para, int dbeta2_perp,
                                           double The, double betac, double muc,
                                           vector<double> &dDn_dThe);

#endif

//==================================================================================================
//==================================================================================================
