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
// Author: Jens Chluba & Elizabeth Lee
//
// first implementation: April 2012
// last modification   : February 2020
//
//==================================================================================================
// 18th Apr  2018: improved precision of second derivative with Te
// 12th Sept 2013: fixed bug because of betac^2 P_2(muc) == beta_para^2 - beta_perp^2/2
//  4st  Aug: added derivatives of basis functions in the CMB rest frame
// 22th July: added low and high temperature expansions. Now 2keV < Te < 75keV is covered
// 21th July: added headers with required data; this avoids time taken for reading the files
// 10th July: added basis functions in the CMB rest frame
//  8th July: changed definition of S^kin; temperature-independent terms are fully canceled
// 02/2020: Added Class structure

#ifndef SZ_CNSN_BASIS_H
#define SZ_CNSN_BASIS_H

using namespace std;

#include "physical_consts.h"
#include "routines.h"
#include "Relativistic_MB.h"
#include "nPl_derivatives.h"
#include "Parameters.h" 

class CNSNsplineMembers{
    public:
    vector<int> Y, D, Q, Mcorr;
    double The_ref;
    bool loaded;

    CNSNsplineMembers();

    private:
    void setup_expansion_splines(vector<int> &spline_mem_indices, const double D[][22]);

    public:
    void loadSplines(double refThe, const double dY[][22],const double dD[][22],const double dQ[][22],const double dMcorr[][22]);
};

//==================================================================================================
//
// compute SZ effect using interpolation of derivative terms
//
//==================================================================================================
class IntegralCNSN{
    private:
    double betac, muc, The, xfac, x;
    int Te_order, betac_order;
    string run_mode;
    bool CMBframe;
    double x3, exp_mx, dex;
    double beta_para, beta_perp;
    CNSNsplineMembers spline;
    vector<double> dDn_dThe;

    public:
    IntegralCNSN();
    IntegralCNSN(double xfac_i, double x_i, double The_i, double betac_i, double muc_i, int Te_order_i, int betac_order_i, bool CMB);
    IntegralCNSN(double x_i, double The_i, double betac_i, double muc_i, int Te_order_i, int betac_order_i, bool CMB);
    IntegralCNSN(double x_i, double The_i, double betac_i, double muc_i, bool CMB);
    IntegralCNSN(int k, Parameters fp, bool CMB, bool inputOrders=true);
    void Update_x(double x_i);

    private:
    void setSpline();

    // generic function for interpolation (x^3 Dn was used). Outputs Dn
    double Dn_SZ_splines(double The_ref, vector<int> &spline_mem);

    double Calculate_monopole();
    double Calculate_dipole();
    double Calculate_quadrupole();
    double Calculate_monopole_correction();

    double Calculate_kinetic_correction();
    double Calculate_All();

    // analytic derivatives in betac_parallel and betac_perp
    // Note since betac muc = beta_para, betac^2 P_2(muc) == beta_para^2 - 0.5*beta_perp^2
    // These are all calculated assuming it has already been set betac = muc = 1.
    double Dn_for_The();
    double Dn_dbeta_para();
    double Dn_d2beta_para();
    double Dn_dbeta2_perp();

    void compute_all_Te_derivatives_upto_dThe(int dThe, std::function<double()> f);

    public:
    double compute_distortion(string mode);
    void Dcompute_distortion(int dThe, int dbeta_para, int dbeta2_perp, vector<double> &dDn);
};

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
double compute_SZ_distortion_CNSN_basis(double x, double The, double betac, double muc, 
                                        int Te_order, int betac_order, 
                                        string mode, bool CMBframe = true);
void compute_SZ_distortion_CNSN_basis(vector<double> &Dn, vector<double> x,
                                      double The, double betac, double muc, int Te_order, int betac_order,
                                      string mode, bool DI, bool CMBframe = true); 

double compute_SZ_distortion_CNSN_basis(int k, Parameters fp, bool CMBframe = false);
void compute_SZ_distortion_CNSN_basis(vector<double> &Dn, Parameters fp, bool DI, bool CMBframe = false);

//==================================================================================================
// Derivatives (The^k d^k_dThe /k!) (betapara^m d^m_dbetapara /m!) (beta2perp^l d^l_beta2perp /l!) S
// in the CMB frame for a resting observer. Maximal orders in The and betac are used to compute 
// the derivatives.
//
// constraints: dThe<=4; dbeta_para<=2; dbeta2_perp<=1;
//==================================================================================================
void Dcompute_SZ_distortion_CNSN_basis(double x, int dThe, int dbeta_para, int dbeta2_perp,
                                       double The, double betac, double muc,
                                       vector<double> &dDn_dThe, bool CMBframe = true);

void Dcompute_SZ_distortion_CNSN_basis(double x, Parameters &fp, bool CMBframe);
#endif

//==================================================================================================
//==================================================================================================
