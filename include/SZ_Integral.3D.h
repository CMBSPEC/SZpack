//==================================================================================================
//
// Computation of the thermal and k-SZ effect using explicit integration of the collision term in
// different orders of beta. The cluster is assumed to be isothermal and moving at a given speed 
// betac with direction muc relative to the line of sight. The integrals are carried out in the 
// cluster frame. 
//
//==================================================================================================
//
// Author: Jens Chluba  & Elizabeth Lee
//
// first implementation: April 2012
// last modification   : February 2020
//
//==================================================================================================
// 8th July: changed definition of S^kin, so that the temperature-independent term is fully canceled
// 02/2020: Added Class structure

#ifndef SZ_INTEGRAL_3D_H
#define SZ_INTEGRAL_3D_H

#include "Relativistic_MB.h"
#include "nPl_derivatives.h"
#include "Integration_routines.h"
#include "Parameters.h"

using namespace std;

//==================================================================================================
//
// Class Structure
//
//==================================================================================================
class Integral3D{
    private:
    double Int_eps, betac, muc, The, x;
    string run_mode;
    double mu, mup, xi;
    double xfac, dx, xp;
    double dsig0, dsig1, dsig2, nfac;

    double exp_mx, dex, exp_mxp, dexp;

    public:
    Integral3D();
    Integral3D(double x_i, double The_i, double betac_i, double muc_i, double eps_Int_i);
    Integral3D(int k, Parameters fp);
    void Update_x(double x_i);

    private:
    void Calculate_shared_variables();

    double Calculate_monopole();
    double Calculate_dipole();
    double Calculate_quadrupole();
    double Calculate_monopole_correction();
    double Calculate_kinetic_correction();
    double Calculate_All();

    double sig_Boltzmann_Compton(double mug);

    double mup_Int(double mue);
    double mue_Int(double xi_int);
    double xi_Int();

    public:
    double compute_distortion(string mode);
};

//==================================================================================================
//
// 3D integration carried out using Patterson scheme
//
// eps_Int : relative accuracy for numerical integration (lower than 10^-6 is hard to achieve)
//
// mode == "monopole"      --> only scattering of monopole without second order kinematic corr
// mode == "dipole"        --> only scattering of dipole     (first order kinematic correction)
// mode == "quadrupole"    --> only scattering of quadrupole (second order kinematic correction)
// mode == "monopole_corr" --> only scattering of second order kinematic correction to monopole
// mode == "all"           --> all terms added
// mode == "kin"           --> only kinematic terms
//
//==================================================================================================

double compute_SZ_distortion_Patterson_3D(double x, 
                                          double The, double betac, double muc,
                                          double eps_Int, string mode = "all");
void compute_SZ_distortion_Patterson_3D(vector<double> &Dn, vector<double> x,
                                        double The, double betac, double muc, double eps_Int,
                                        string mode, bool DI);      

double compute_SZ_distortion_Patterson_3D(int k, Parameters fp);
void compute_SZ_distortion_Patterson_3D(vector<double> &Dn, Parameters &fp, bool DI);


#endif

//==================================================================================================
//==================================================================================================