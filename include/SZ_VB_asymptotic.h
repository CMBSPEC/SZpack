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

#ifndef SZ_VB_ASYMPTOTIC_H
#define SZ_VB_ASYMPTOTIC_H

#include <functional>
#include <string>
#include <vector>

class Parameters;

using namespace std;

typedef std::function<double (int, double)> distributionDerivs;

class NonThermalAsymptotic{
    public:
    vector<double> Y, D, Q, Mcorr;

    private:
    double xfac, betac, muc, The, x;
    int Te_order, betac_order;
    bool CMBframe;
    string run_mode;
    distributionDerivs Ddist;
    vector<double> Pk, beta;

    public:
    NonThermalAsymptotic();
    NonThermalAsymptotic(double xfac_i, double x_i, double The_i, double betac_i, double muc_i, 
                         int Te_order_i, int betac_order_i, distributionDerivs Ddist_i, bool CMB);
    NonThermalAsymptotic(double x_i, double The_i, double betac_i, double muc_i, 
                         int Te_order_i, int betac_order_i, distributionDerivs Ddist_i, bool CMB);
    NonThermalAsymptotic(double x_i, Parameters fp, distributionDerivs Ddist_i, bool CMB, bool inputOrders=true);
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

    public:
    double compute_distortion(string mode);
};

//==================================================================================================
//
// compute the nonthermal Dn using asymptotic expansion in cluster frame
//
// mode == "monopole"      --> only scattering of monopole without second order kinematic corr
// mode == "dipole"        --> only scattering of dipole     (first order kinematic correction)
// mode == "quadrupole"    --> only scattering of quadrupole (second order kinematic correction)
// mode == "monopole_corr" --> only scattering of second order kinematic correction to monopole
// mode == "all"           --> all terms added
// mode == "kin"           --> only kinematic terms
//
//==================================================================================================
double compute_nonthermal_distortion_asymptotic(double x, double The, double betac, double muc, 
                                                int Te_order, int betac_order, 
                                                std::function<double(int, double)> Ddist, 
                                                string mode, bool CMBframe = true);
void compute_nonthermal_distortion_asymptotic(vector<double> &Dn, vector<double> x,
                                              double The, double betac, double muc, int Te_order, 
                                              int betac_order, std::function<double(int, double)> Ddist,
                                              string mode, bool DI, bool CMBframe = true);

double compute_nonthermal_distortion_asymptotic(double x, Parameters fp, std::function<double(int, double)> Ddist, 
                                                bool CMBframe = false);
void compute_nonthermal_distortion_asymptotic(vector<double> &Dn, Parameters fp, 
                                              std::function<double(int, double)> Ddist, bool DI, 
                                              bool CMBframe = false);

//==================================================================================================
//
// Radio SZ functions -- in particular power law functions
//
//==================================================================================================

// The derivatives for the function x^k d^k/dx^k (x^-gamma)
double PowerLawDerivs(int k, double x, double gamma);

double RadioDerivs(int k, double x);
void compute_radio_distortion(vector<double> &Dn, Parameters fp, bool DI, bool CMBframe = false);
double compute_radio_distortion(double x, Parameters fp, bool CMBframe = false);

double compute_null_of_combined_signal(Parameters fp);
double compute_null_of_combined_signal(Parameters fp, std::function<double(int, double)> Ddist);


#endif

//==================================================================================================
//==================================================================================================