//==================================================================================================
// Author Jens Chluba Jan 2010
// Last modification: April 2011
//==================================================================================================
#ifndef PATTERSON_H
#define PATTERSON_H

# include <cstdlib>
# include <cmath>
# include <iostream>
# include <iomanip>
# include <cstring>
# include <vector>

# include "routines.h"
# include "global_functions.h"
//==================================================================================================
// patterson formulae & integration
// Integral_-1^1 f(x) dx = f(0)+Sum _1^order w[i]*(f(0-x[i])+f(0+x[i]))
//==================================================================================================
int Integrate_using_Patterson(double a, double b, double epsrel, double epsabs, 
                              double (*fptr)(double), int *neval, double *r);

int Integrate_using_Patterson(double a, double b, double epsrel, double epsabs, 
                              double (*fptr)(double, void *p), int *neval, double *r, void *p);
    
double Integrate_using_Patterson_adaptive(double a, double b, double epsrel, double epsabs, 
                                          double (*fptr)(double));

double Integrate_using_Patterson_adaptive(double a, double b, double epsrel, double epsabs, 
                                          double (*fptr)(double, void *p), void *p);


double Integrate_using_Patterson_adaptive(double a, double b, double epsrel, double epsabs, 
                                          std::function<double(double)> f);

#endif  

//==================================================================================================
//==================================================================================================
