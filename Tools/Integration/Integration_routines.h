//==================================================================================================
// Author Jens Chluba Jan 2010
//==================================================================================================
#ifndef INTEGRATION_ROUTINES_H
#define INTEGRATION_ROUTINES_H

#include <iostream>
#include <string>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <complex>
#include <vector>

#include "routines.h"
#include "Chebyshev_Int.h"
#include "Patterson.h"

using namespace std;

//==================================================================================================
// routines for integration
//==================================================================================================
double Integrate_gk15_GSL(double a, double b, double epsrel, double epsabs,
                          double (*fptr)(double, void *), void *para=NULL);

double Integrate_gk31_GSL(double a, double b, double epsrel, double epsabs, 
                          double (*fptr)(double, void *), void *para=NULL);

double Integrate_gk61_GSL(double a, double b, double epsrel, double epsabs, 
                          double (*fptr)(double, void *), void *para=NULL);

double Integrate_gk_GSL(int key, double a, double b, double epsrel, double epsabs, 
                        double (*fptr)(double, void *), void *para=NULL);

#endif

//==================================================================================================
//==================================================================================================
