//==================================================================================================
//
// //TODO: Fill this box
//
//==================================================================================================
//
// Author: Elizabeth Lee
// Based off work by Jens Chluba
//
//==================================================================================================

#ifndef SZ_INTEGRAL_MOMENT_H
#define SZ_INTEGRAL_MOMENT_H

#include "SZ_moment_method.h"
#include "routines.h"
#include "Integration_routines.h"
#include "Patterson.h"
#include "nPl_derivatives.h"
#include "physical_consts.h"

using namespace std;

extern profile_slice_params PSP;
extern Parameters parameters;

//==================================================================================================
double compute_Integral(double (*func)(double l), double Int_eps, double lmin = -2.0, double lmax = 2.0);

//==================================================================================================
//
// functions to compute y, z, b, and c integrals according to CSNN 2012
//
//==================================================================================================
typedef double (*Integrand_1D)(double lz);
typedef double (*Integrand_1D_yk)(double lz, double k);
typedef double (*normFunc)(double The);

static double The_ref_global;
static double g_eps_Int = 1.0e-6;

struct integrationValues{
    double lx0, ly0, lsc, lang;
    Integrand_1D integrand;
    bool oneD;
    bool usek;

    double lx, ly;

    Integrand_1D_yk integrand_yk;
    int k;

    integrationValues();
    integrationValues(bool OneD);
    integrationValues(double Lx0, double Ly0, double Lsc, double Lang, Integrand_1D Integrand, bool OneD = false);
    integrationValues(double Lx0, double Ly0, double Lsc, double Lang, Integrand_1D_yk Integrand, int K, bool OneD = false);
};

extern integrationValues IValues;
void setNormFunc(normFunc NF);

//==================================================================================================
double compute_tau(double zsc);

double compute_ySZ(double zsc);
                              
double compute_yk(double zsc, int k);

//==================================================================================================
// low temperature moments
//==================================================================================================
double compute_y_k_cond(int k, int ir, double zsc, SZ_moment_method SZM);

double compute_b_k_0_cond(int k, int ir, double zsc, SZ_moment_method SZM);
double compute_b_0(double zsc);

double compute_b_k_1_cond(int k, int ir, double zsc, SZ_moment_method SZM);
double compute_b_1(double zsc);

double compute_b_k_2_cond(int k, int ir, double zsc, SZ_moment_method SZM);
double compute_b_2(double zsc);

//==================================================================================================
//
// additional average over x and y
// 
//==================================================================================================
double compute_Int_3D(double lx0, double ly0, double lsc, double lang, Integrand_1D integrand);

double compute_Int_3D_yk(double lx0, double ly0, double lsc, double lang, int k,
                         Integrand_1D_yk integrand);

//==================================================================================================
//
// computing temperature moments using VIK
//
//==================================================================================================
double integrate_yk_VIK(int k);

//==================================================================================================
//
// computing SZ signal by explicitly integrating along the line of sight for each
// frequency bin. This routine is just for comparisons & rather slow.
//
//==================================================================================================
double integrate_SZ_VIK(double xcmb);

//==================================================================================================
//==================================================================================================
#endif
