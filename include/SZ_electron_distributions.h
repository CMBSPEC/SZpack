//==================================================================================================
//
// Electron distributions for high energy non thermal electron contributions to the SZ signal. These
// are used by SZ_Integral.Kernel. Distributions taken from TODO:
//
//==================================================================================================
//==================================================================================================
//
// Author: Elizabeth Lee
//
// first implementation: September 2020
// last modification   : September 2020
//
//==================================================================================================

#ifndef SZ_ELECTRON_DISTRIBUTIONS_H
#define SZ_ELECTRON_DISTRIBUTIONS_H

#include <gsl/gsl_sf_bessel.h>
#include "physical_consts.h"
#include "routines.h"
#include "Relativistic_MB.h"
#include "Integration_routines.h"
#include "Patterson.h"

//TODO: Fill in the sources for these
//==================================================================================================
// The electron momentum distribution for thermal electrons
//==================================================================================================
double Boltzmann_Dist(double eta, double Te);
double Boltzmann_Dist_gamma(double eta, double gamma, double Te);

//==================================================================================================
// kinematic boost model
//==================================================================================================
double KinematicBoost_Dist(double eta, double Te, double betac, double muc, int l);
double KinematicBoost_Dist_exp(double eta, double Te, double betac, double muc, int l, int betac_order);


//==================================================================================================
// Different models for temperature distributions
//==================================================================================================
double CosmicRay_Dist(double eta, double alpha = 2.5, double p1 = 0.1, double p2 = 10.0);

double ThermalCosmicRay_Dist(double eta, double Te, double alpha = 2.5, double p1 = 0.2, double p2 = 10.0);
double ThermalCosmicRay_Norm(double Te, double alpha = 2.5, double p1 = 0.2, double p2 = 10.0);

double DoublePower_Dist(double eta, double alpha1 = 0.5, double alpha2 = 2.5, 
                        double p1 = 0.01, double pcr = 0.2, double p2 = 10.0);

double kappa_Dist(double eta, double Te, double kappa = 2.0);

double MultiMaxwellian_Dist(double eta, double Te, vector<double> c, vector<double> a);
double MultiMaxwellian_Dist(double eta, double Te);

#endif