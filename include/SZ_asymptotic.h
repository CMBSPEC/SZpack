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
// Author: Jens Chluba & Elizabeth Lee
//
// first implementation: April 2012
// last modification   : February 2020
//
//==================================================================================================
//  4st  Aug: added derivatives of basis functions in the CMB rest frame
// 10th July: added basis functions in the CMB rest frame
// 02/2020: Added Class structure

#ifndef SZ_ASYMPTOTIC_H
#define SZ_ASYMPTOTIC_H

//==================================================================================================
// required libs
//==================================================================================================
#include "physical_consts.h"
#include "routines.h"
#include "nPl_derivatives.h"
#include "Parameters.h" 

using namespace std;

class IntegralAsymptotic{
    public:
    vector<double> Y, D, Q, Mcorr;

    private:
    double xfac, betac, muc, The, x;
    int Te_order, betac_order;
    double beta_para, beta_perp;
    string run_mode;
    bool CMBframe;
    double exp_mx, dex;
    double fac, nPl;
    vector<double> Pk, beta;
    vector<double> dDn_dThe;

    public:
    IntegralAsymptotic();
    IntegralAsymptotic(double xfac_i, double x_i, double The_i, double betac_i, double muc_i, int Te_order_i, int betac_order_i, bool CMB);
    IntegralAsymptotic(double x_i, double The_i, double betac_i, double muc_i, int Te_order_i, int betac_order_i, bool CMB);
    IntegralAsymptotic(double x_i, double The_i, double betac_i, double muc_i, bool CMB);
    IntegralAsymptotic(double x_i, bool CMB);
    IntegralAsymptotic(int k, Parameters fp, bool CMB, bool inputOrders=true);
    void Update_x(double x_i);
    
    private:
    void Calculate_shared_variables();

    void Compute_Y();
    void Compute_D();
    void Compute_Q();
    void Compute_Mcorr();

    // approximation for relativistic SZ effect
    // const temperature case y^(k)= Theta^k tau_e
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

    //==================================================================================================
    void compute_all_Te_derivatives_upto_dThe(int dThe, std::function<double()> f);

    public:
    double compute_distortion(string mode);
    void Dcompute_distortion(int dThe, int dbeta_para, int dbeta2_perp, vector<double> &dDn);
};


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
double compute_SZ_distortion_asymptotic(double x, double The, double betac, double muc, 
                                        int Te_order, int betac_order, 
                                        string mode, bool CMBframe = true);
void compute_SZ_distortion_asymptotic(vector<double> &Dn, vector<double> x,
                                      double The, double betac, double muc, int Te_order, int betac_order,
                                      string mode, bool DI, bool CMBframe = true);                                 

double compute_SZ_distortion_asymptotic(int k, Parameters fp, bool CMBframe = false);
void compute_SZ_distortion_asymptotic(vector<double> &Dn, Parameters fp, bool DI, bool CMBframe = false);

//==================================================================================================
// Derivatives (The^k d^k_dThe /k!) (betapara^m d^m_dbetapara /m!) (beta2perp^l d^l_beta2perp /l!) S
// in the CMB frame for a resting observer. Maximal orders in The and betac are used to compute 
// the derivatives.
//
// constraints: dThe<=4; dbeta_para<=2; dbeta2_perp<=1;
//==================================================================================================
void Dcompute_SZ_distortion_asymptotic(double x, int dThe, int dbeta_para, int dbeta2_perp,
                                       double The, double betac, double muc,
                                       vector<double> &dDn_dThe, bool CMBframe = true);

void Dcompute_SZ_distortion_asymptotic(double x, Parameters &fp, bool CMBframe);

#endif

//==================================================================================================
//==================================================================================================