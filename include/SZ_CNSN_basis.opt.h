//==================================================================================================
//
// This program allows computing the thermal SZ and kinematic effect using the improved basis of 
// Chluba, Nagai, Sazonov & Nelson, 2012. Everything is similar to the routines provided in 
// SZ_CNSN_basis.h, but here the computations are carried out in the CMB frame only. Since only terms
// up to O(betac^2) are included, this leads to some tiny differences of the results. However, for 
// the temperature-velocity moment method explained by Chluba, Switzer, Nelson & Nagai, 2012, this
// choice is beneficial. 
//
//==================================================================================================
//
// Author: Jens Chluba  (CITA, University of Toronto)
//
// first implementation: July 2012
// last modification   : July 2012
//
//==================================================================================================

#ifndef SZ_CNSN_BASIS_OPT_H
#define SZ_CNSN_BASIS_OPT_H

using namespace std;


//==================================================================================================
//
// access basis functions and Normalization, etc. These functions are needed by the SZ moment 
// method but are otherwise not too useful for general purpose applications.
//
//==================================================================================================
double N_func_CNSN(double The);
void compute_Y_CNSN    (double x, int region, vector<double> &Y);
void compute_M_CNSN_CMB(double x, int region, vector<double> &MCMB);
void compute_D_CNSN_CMB(double x, int region, vector<double> &DCMB);
void compute_Q_CNSN_CMB(double x, int region, vector<double> &QCMB);

//==================================================================================================
// computes the optimal value for kmax given the accuracy goal and required maximal temperature
//==================================================================================================
void determine_optimal_kmax(int accuracy_level, double Te_max, int &kmax, int &iregmax);

vector<double> Get_temperature_regions(int kmax, int accuracy_level);
vector<double> Get_temperature_pivots (int kmax, int accuracy_level);
vector<int> Get_region_indices (int kmax, int accuracy_level);


//==================================================================================================
//
// compute Dn using improved expansion in CMB rest frame
//
// mode == "monopole"      --> only monopole part without second order kinematic corr
// mode == "dipole"        --> only dipolar part     (first order kinematic correction)
// mode == "quadrupole"    --> only quadrupolar part (second order kinematic correction)
// mode == "monopole_corr" --> only second order kinematic correction to monopole part
// mode == "all"           --> all terms added
// mode == "kin"           --> only kinematic terms
//
//==================================================================================================
double compute_SZ_distortion_CNSN_basis_opt(double x, 
                                            double The, double betac, double muc, 
                                            int kmax, int betac_order,
                                            string mode, int accuracy_level);

#endif

//==================================================================================================
//==================================================================================================
