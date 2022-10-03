//==================================================================================================
//
// Handling the output from SZpack, both to file and to screen
//
//==================================================================================================

#ifndef output_MultipleScattering_H
#define output_MultipleScattering_H

#include "output.h"
#include "SZ_moment_method.h"
#include "SZpackMultipleScattering.h"

using namespace std;

extern Parameters parameters;

//TODO: All of these functions cout the xcmb as they pass through, instead of anything to determine if it's done the correct thing
// Also could do with putting this output behind a show mess variable and such.

typedef double (*t_0)(double);
typedef double (*t_av)(double, int);

template <typename somestream>
void SetUpOutput(string description, string fileAddition, Parameters &fp, somestream &ofile, int precision, bool printHeader);

//==================================================================================================
//
// plot <tau_l>
//
//==================================================================================================
void output_tau_l(int lmax, int np, vector<double> &xc, t_0 tau0, t_av tau_av, ofstream &ofile);

void output_tau_l_sphere(int lmax, int np);

void output_tau_l_beta(int lmax, int np, double beta);

//==================================================================================================
//
// plot of Y_l0 functions for different temperatures
//
//==================================================================================================
void output_DI2_E_temperature(Parameters &fp);

//==================================================================================================
//
// plot of Y_l0 functions for different temperatures
//
//==================================================================================================
void output_SZ_signal_l(Parameters &fp);

//==================================================================================================
//
// plot emission/absorption terms
//
//==================================================================================================
void output_emission_absorption(Parameters &fp);

//==================================================================================================
void output_emission_absorption_kinetic(Parameters &fp);

//==================================================================================================
//
// CMB anisotropies signal
//
//==================================================================================================
void output_CMB_isotropy_signals(Parameters &fp);

//==================================================================================================
//
// plot first few correction functions
//
//==================================================================================================
void output_CMB_isotropy_Y_functions(Parameters &fp);

void output_lowest_order_signal_model(double b, double ISA, double scale, double The, double x3,
                                      double Y0, double Yl0, t_0 tau0, t_av tau_av, ofstream &ofile);

//==================================================================================================
//
// plot lowest order signal
//
//==================================================================================================
void output_lowest_order_signal(Parameters &fp, double b);

#endif
