//==================================================================================================
//
// Computation of the thermal and k-SZ effect using explicit integration of the 5D collision term. 
// The cluster is assumed to be isothermal and moving at a given speed betac with direction muc 
// relative to the line of sight. The integrals are carried out in the cluster frame. Optionally 
// all higher order terms in betac can be included. Use of this routine is only recommended for 
// checks of the results, as it is rather slow.
//
//==================================================================================================
//
// Author: Jens Chluba & Elizabeth Lee
//
// first implementation: May 2012
// last modification: February 2020
//
//==================================================================================================
// 02/2020: Recreated into class form
//
//==================================================================================================

//==================================================================================================

#ifndef SZ_INTEGRAL_5D_H
#define SZ_INTEGRAL_5D_H

#include "Relativistic_MB.h"
#include "nPl_derivatives.h"
#include "Integration_routines.h"
#include "Parameters.h" 

using namespace std;

class Integral5D{
    private:
    double Int_eps, betac, muc, The, x;
    string run_mode;
    double mu, mup, xi, phi, phip;
    double xfac, dx, xp;
    double dsig, mucp;

    double exp_mx, dex, exp_mxp, dexp;
    double nx, nxp, G, Gp;

    public:
    Integral5D();
    Integral5D(double x_i, double The_i, double betac_i, double muc_i, double eps_Int_i);
    Integral5D(double x_i, Parameters fp);
    void Update_x(double x_i);

    private:
    void Calculate_shared_variables();

    double Calculate_full();

    //The cut higher order version from before
    void Calculate_HigherOrder_Variables();

    double Calculate_monopole();
    double Calculate_dipole();
    double Calculate_quadrupole();
    double Calculate_monopole_correction();
    double Calculate_kinetic_correction();
    double Calculate_All();

    double sig_Boltzmann_Compton(double int_phip);

    double phip_Int(double int_phi);
    double phi_Int(double int_mup);
    double mup_Int(double mu_int);
    double mue_Int(double xi_int);
    double xi_Int();

    public:
    double compute_distortion(string mode);
};

//==================================================================================================
//
// 5D integration carried out using Patterson scheme. The SZ signal is obtained in the cluster frame
//
//==================================================================================================
double compute_SZ_distortion_Patterson_5D(double x, double The, double betac, double muc,
                                          double eps_Int, string mode = "full");
void compute_SZ_distortion_Patterson_5D(vector<double> &Dn, vector<double> x,
                                        double The, double betac, double muc, double eps_Int,
                                        string mode, bool DI);

double compute_SZ_distortion_Patterson_5D(double x, Parameters fp);
void compute_SZ_distortion_Patterson_5D(vector<double> &Dn, Parameters fp, bool DI);

#endif

//==================================================================================================
//==================================================================================================