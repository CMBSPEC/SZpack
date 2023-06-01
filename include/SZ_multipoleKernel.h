//==================================================================================================
//
// Computation of the thermal and k-SZ effect using explicit integration collision term.
// Using a precomputed kernel over the angular directions. 
// TODO: Update this!!!!
// The cluster is assumed to be isothermal and moving at a given speed betac with direction muc 
// relative to the line of sight. The integrals are carried out in the cluster frame.
//
//==================================================================================================
//
// Author: Elizabeth Lee
// Based on work by Jens Chluba (CITA, University of Toronto)
//
// first implementation: October 2018
// last modification: March 2021
//
//==================================================================================================

#ifndef SZ_MULTIPOLE_KERNEL_H
#define SZ_MULTIPOLE_KERNEL_H

#include "nPl_derivatives.h"
#include "Parameters.h"
#include "Integration_routines.h"
#include "Patterson.h"
#include "routines.h"
#include "Definitions.h"
//#include "SZ_electron_distributions.h"

using namespace std;

//TODO: Check that this electron distribution stuff is correct
typedef std::function<double (double)> electronDistribution;
//this defines our electron distributions as functions of l and eta

static double plainDistribution(double eta){
    return 1.0;
}

class MultipoleKernel{
    private:
    int l;
    double s, t, eta, gamma0;
    double Int_eps, mus, mup, r; 
    double Lpart;
    vector<double> P, h0, h2, h4, K1, K2; 
    double s_low, s_high;

    public:
    MultipoleKernel();
    MultipoleKernel(int l_i, double s_i, double eta_i);
    MultipoleKernel(int l_i, double s_i, double eta_i, double Int_eps_i);
    void Update_s(double s_i);
    void Update_l(int l_i);
    bool electron_anisotropy;

    private:
    void Calculate_integral_variables();
    double sigl_Boltzmann_Compton(double mup_int);
    double mup_Int(double mus_int);
    double mus_Int();
    void Calculate_formula_variables();

    public:
    vector<double> s_limits();
    double Calculate_integrated();
    double Calculate_formula();
    double Calculate_stable();
    double Calculate_electron_multipoles();
};

class BeamKernel{
    private:
    int l;
    double s, t, mup, eta, gamma0;
    double f0, f1, f2; 

    public:
    BeamKernel();
    BeamKernel(int l_i, double s_i, double eta_i, double mup_i);
    void Update_s(double s_i);
    void Update_l(int l_i);

    private:
    void Calculate_monopole();
    void Calculate_dipole();
    void Calculate_quadrupole();

    public:
    vector<double> s_limits();
    vector<double> mup_limits();
    double Calculate_formula();
};

class IntegralKernel{
    private:
    double Int_eps, betac, muc, x;
    string run_mode;
    int l;
    double s, eta, xfac, xp, mup;

    bool beam_kernel, fixed_eta, electron_anisotropy;

    electronDistribution etaDistribution;
    MultipoleKernel MK;
    BeamKernel BK;

    public:
    IntegralKernel();
    IntegralKernel(double x_i, double betac_i, double muc_i, double eps_Int_i);
    IntegralKernel(double x_i, Parameters fp);
    void Update_x(double x_i);

    private:
    void Calculate_shared_variables();
    double Calculate_kernel(int l_i);

    double Calculate_monopole(int l_i);
    double Calculate_dipole(int l_i);
    double Calculate_quadrupole(int l_i);
    double Calculate_monopole_correction(int l_i);

    double sig_Boltzmann_Compton(double int_eta);

    double eta_Int(double int_s);
    double s_Int();

    public:
    double compute_distortion(string mode, electronDistribution eDistribution, int l_i=0, bool e_anis=false);
    double compute_kernel(int l_i, double s_i, electronDistribution eDistribution, bool e_anis=false);
    double compute_distortion_fixed_eta(string mode, double eta_, int l_i=0, bool e_anis=false);

    double compute_beam_distortion(double mup_i, string mode, electronDistribution eDistribution, int l_i=0);
    double compute_beam_kernel(double mup_i, int l_i, double s_i, electronDistribution eDistribution);
    double compute_beam_distortion_fixed_eta(double mup_i, string mode, double eta_i, int l_i=0);

    //TODO: function to compute distortion for fixed p0
    //TODO: All three functions but for the beam kernel instead of the multipole kernel
    //TODO: Look into integration routines to figure out if there is a better one for the variable momentum limits
};

//==================================================================================================
//
// 3D integration carried out using Patterson scheme
//
// eps_Int : relative accuracy for numerical integration (lower than 10^-6 is hard to achieve)
//
// These modes now refer to the background photon population used -- i.e., whether it is n or some derivative of n
// mode == "monopole"      --> only scattering of monopole without second order kinematic corr
// mode == "dipole"        --> only scattering of dipole     (first order kinematic correction)
// mode == "quadrupole"    --> only scattering of quadrupole (second order kinematic correction)
// mode == "monopole_corr" --> only scattering of second order kinematic correction to monopole
// mode == "all"           --> Not Implemented here!
// mode == "kin"           --> Not Implemented here!
// Note these cannot be simply combined to access the kinematic corrections!
//
//==================================================================================================

void compute_SZ_distortion_kernel(vector<double> &Dn, Parameters &fp, bool DI, std::function<double(double)> eDistribution, int l=0, bool e_anis=false);
void compute_averaged_kernel(vector<double> &Dn, Parameters &fp, std::function<double(double)> eDistribution, int l=0, bool e_anis=false);
void compute_SZ_distortion_fixed_eta(vector<double> &Dn, Parameters &fp, bool DI, double eta, int l=0,  bool e_anis=false);

void compute_SZ_distortion_beam_kernel(vector<double> &Dn, Parameters &fp, bool DI, double mup, 
                                        std::function<double(double)> eDistribution, int l=0);
void compute_averaged_beam_kernel(vector<double> &Dn, Parameters &fp, double mup, std::function<double(double)> eDistribution, int l=0);
void compute_SZ_distortion_beam_kernel_fixed_eta(vector<double> &Dn, Parameters &fp, bool DI, double mup, double eta, int l=0);

#endif