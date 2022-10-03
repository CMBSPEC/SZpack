//==================================================================================================
//
// Handling the output from SZpack, both to file and to screen
//
//==================================================================================================

#ifndef output_Moment_H
#define output_Moment_H

#include "output.h"
#include "SZ_moment_method.h"
#include "SZpackMoment.h"

using namespace std;

extern Parameters parameters;

//==================================================================================================
template<typename somestream>
void output_moment_vector(somestream &output, const vector<double> &mv, Parameters fp, SZ_moment_method SZM);

void saveAndPrintMomentVector(string fileAddition, vector<double> &mv, Parameters fp, SZ_moment_method SZM);

//==================================================================================================
//
// header for matrix info
//
//==================================================================================================
template<typename somestream>
void header_for_moment_matrix(somestream &output, string description, Parameters fp, SZ_moment_method SZM);

void export_moment_matrix(string description, simpleMatrix &B, Parameters fp, SZ_moment_method SZM);

void export_vector(vector<double> A, string description = "unknown", Parameters fp = parameters);

template<typename somestream>
void output_profile_slice(somestream &output, double (*f)(double l), double lsc, int np);

void show_profile_slice(double (*f)(double l), double lsc, int np);
void export_profile_slice(string fname, double (*f)(double l), double lsc, int np);

//==================================================================================================
//
// This function gives an example of how the moment method and the class `SZ_moment_method' are used
// to compute the SZ signal. Two variables type for the high temperature moments are available, 
// which should all give practically the same SZ signal. The different steps in the calculation are
// commented below and one can play with things a bit. //TODO: Sort comments!
//
//==================================================================================================
void output_distortion_moments_isothermal(SZ_moment_method &SZM, string mess = "", Parameters fp = parameters);

//==================================================================================================
//
// Calculations for computing SZ signal for different cluster profile according to 
// fits of Vikhlinin et al 2006
//
//==================================================================================================
void output_distortion_moments_cluster_VIK(string fileAddition, SZ_cluster_profiles &CL, 
                                           double x_rc, double y_rc, SZ_moment_method &SZM,
                                           Parameters fp, string mess = "");
//--------------------------------------------------------------------------------------------------
void output_distortion_cluster_VIK_explicit(string fileAddition, SZ_cluster_profiles &CL,
                                               double x_rc, double y_rc, Parameters fp);

//==================================================================================================
//
// Calculates the moments for the VIK method above - and outputs the y_k, rho_k and omega_k associated.
//
//==================================================================================================
template<typename somestream>
void write_yk_ykiso(somestream &f0, somestream &f1, somestream &f2, int kmax);

void output_distortion_cluster_yk(double y_rc, SZ_cluster_profiles &CL, Parameters fp);

//==================================================================================================
//
// Compute degeneracy factors
//
//==================================================================================================
void compute_degeneracy_functions(Parameters fp = parameters);

//==================================================================================================
void output_TSZ_bestfit(SZ_cluster_profiles &CL, SZ_moment_method &SZM, Parameters fp, double x_rc_max, int np);
//TODO: This is not currently used!

#endif
