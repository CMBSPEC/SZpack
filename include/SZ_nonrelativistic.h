//==================================================================================================
//
// This program allows computing the traditional non-relativistic thermal SZ and kinematic effect. 
// i.e, as described in Zeldovich & Sunyaev 1969 and Sunyaev & Zeldovich 1980
//
//==================================================================================================
//
// Author: Elizabeth Lee
//
// first implementation: August 2020
//
//==================================================================================================

#ifndef SZ_NONREL_H
#define SZ_NONREL_H

#include <string>
#include <vector>

class Parameters;

using namespace std;

//==================================================================================================
//
// Class Structure
//
//==================================================================================================
class IntegralNonRelativistic{
    private:
    double betac, muc, The, x;
    string run_mode;
    double exp_mx, dex;

    public:
    IntegralNonRelativistic();
    IntegralNonRelativistic(double x_i, double The_i, double betac_i, double muc_i);
    IntegralNonRelativistic(double x_i, Parameters fp);
    void Update_x(double x_i);

    private:
    double Calculate_monopole();
    double Calculate_dipole();
    double Calculate_All();

    public:
    double compute_distortion(string mode);
};

//==================================================================================================
//
// Numerical solution directly calculated
//
// mode == "monopole"      --> only scattering of monopole without second order kinematic corr
// mode == "dipole"        --> only scattering of dipole     (first order kinematic correction)
// mode == "quadrupole"    --> Automatically 0 within this regime
// mode == "monopole_corr" --> Automatically 0 within this regime
// mode == "all"           --> all terms added
// mode == "kin"           --> equivalent to dipole.
//
//==================================================================================================

double compute_SZ_distortion_nonrelativistic(double x, double The, double betac, double muc, string mode = "all");
void compute_SZ_distortion_nonrelativistic(vector<double> &Dn, vector<double> &x,
                                           double The, double betac, double muc, string mode, bool DI);

double compute_SZ_distortion_nonrelativistic(double x, Parameters fp);
void compute_SZ_distortion_nonrelativistic(vector<double> &Dn, Parameters &fp, bool DI);

#endif

//==================================================================================================
//==================================================================================================