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
// Author: Jens Chluba & Elizabeth Lee
//
// first implementation: July 2012
// last modification   : February 2020
//
//==================================================================================================
// 02/2020: Added Class structure

#ifndef SZ_CNSN_BASIS_OPT_H
#define SZ_CNSN_BASIS_OPT_H

#include <vector>

using namespace std;

#include "physical_consts.h"
#include "routines.h"
#include "Relativistic_MB.h"
#include "nPl_derivatives.h"
#include "SZ_asymptotic.h"


class CNSNoptSplineMembers{
    public:
    vector<int> Y, D, Q, Mcorr;
    double The_ref;
    bool loaded;

    CNSNoptSplineMembers();

    private:
    void setup_expansion_splines(vector<int> &spline_mem_indices, const double D[][8]);

    public:
    void loadSplines(double refThe, const double dY[][8],const double dD[][8],const double dQ[][8],const double dMcorr[][8]);
};

//==================================================================================================

class precision_settings {
    private:
    int length;

    public:
    int region, kmax;
    vector<double> Te_regions;
    vector<double> Te_pivots;
    vector<int> pivot_indices;

    private:
    int get_pivot_index(int reg);
    
    public:
    void initialise_settings(int maxk, int len, const vector<double> Tregions, const vector<double> Tpivots);
    int determine_Te_region(double Te_max);
};

//==================================================================================================
//
// SZ effect using interpolation of derivative terms
//
//==================================================================================================

class IntegralCNSNopt{
    public:
    vector<double> Y, D, Q, Mcorr;
    precision_settings accuracy;
    
    private:
    double betac, muc, The, x, xfac;
    int kmax, Te_order, betac_order;
    string run_mode;
    bool CMBframe;
    double x3, exp_mx, dex;

    int region, piv;
    CNSNoptSplineMembers spline;

    public:
    IntegralCNSNopt();
    IntegralCNSNopt(double xfac_i, double x_i, double The_i, double betac_i, double muc_i, int kmax_i, int accuracy_level, int betac_order_i, bool CMB);
    IntegralCNSNopt(double x_i, double The_i, double betac_i, double muc_i, int kmax_i, int accuracy_level, int betac_order_i, bool CMB);
    IntegralCNSNopt(int kmax, int accuracy_level); //For Accuracy settings only
    IntegralCNSNopt(double x_i, int region, int Te_order_i); //for calculating the basis functions
    IntegralCNSNopt(double x_i, Parameters fp, bool CMB);
    void Update_x(double x_i);

    private:
    // Get reference to required accuracy setting
    void getAccuracySettings(int kmax, int accuracy_level);

    // generic function for interpolation (x^3 Dn was used). Outputs Dn
    double Dn_SZ_splines(double The_ref, vector<int> &spline_mem);

    // access basis functions and Normalization, etc. These functions are needed by the SZ moment 
    // function to derive the Y, D, Q and Mcorr vectors as in asymptotic
    void Compute_XX(double The_ref, vector<double> &XX, vector<int> &spline_mem);

    void Compute_Y();
    void Compute_D();
    void Compute_Q();
    void Compute_Mcorr();

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
double compute_SZ_distortion_CNSN_basis_opt(double x, double The, double betac, double muc, 
                                            int kmax, int betac_order,
                                            string mode, int accuracy_level, bool CMBframe = true);
void compute_SZ_distortion_CNSN_basis_opt(vector<double> &Dn, vector<double> x,
                                          double The, double betac, double muc, int kmax, int betac_order,
                                          string mode, int accuracy_level, bool DI, bool CMBframe = true); 

double compute_SZ_distortion_CNSN_basis_opt(double x, Parameters fp, bool CMBframe = true);
void compute_SZ_distortion_CNSN_basis_opt(vector<double> &Dn, Parameters fp, bool DI, bool CMBframe = true);

//==================================================================================================
//TODO: These are only used in SZ_moment_method. May be able to write this code better to remove the need for these functions.
//==================================================================================================
vector<double> Get_temperature_regions(int kmax, int accuracy_level);
vector<int> Get_region_indices (int kmax, int accuracy_level);
vector<double> Get_temperature_pivots (int kmax, int accuracy_level);

//==================================================================================================
// computes the optimal value for kmax given the accuracy goal and required maximal temperature
//==================================================================================================
void determine_optimal_kmax(int accuracy_level, double Te_max, int &kmax, int &iregmax);

//==================================================================================================
// access basis functions (always in CMB frame)
//==================================================================================================
void compute_Y_CNSNopt(double x, int region, vector<double> &Y);
void compute_D_CNSNopt(double x, int region, vector<double> &D);
void compute_Q_CNSNopt(double x, int region, vector<double> &Q);
void compute_Mcorr_CNSNopt(double x, int region, vector<double> &Mcorr);

#endif

//==================================================================================================
//==================================================================================================
