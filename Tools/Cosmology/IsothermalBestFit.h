//==================================================================================================
//
// routine to find best-fit solutions using isothermal model
//
//==================================================================================================
//
// Author: Elizabeth Lee
// Based off work by Jens Chluba
//
//==================================================================================================

#ifndef SZ_ISOTHERMAL_H
#define SZ_ISOTHERMAL_H

#include "Parameters.h"
#include "SZpack.h"

#include <gsl/gsl_multimin.h>
#include <gsl/gsl_linalg.h>

extern Parameters parameters;

//--------------------------------------------------------------------------------------------------
struct minimizer_Data
{
    vector<double> xa, Ia;
    vector<double> Iapprox;
};

//--------------------------------------------------------------------------------------------------
void find_best_fit_isothermal_model(vector<double> &xa, vector<double> &Ia,
                                   vector<double> &solution, vector<double> &Isol, 
                                   int numpar = 2);

#endif
