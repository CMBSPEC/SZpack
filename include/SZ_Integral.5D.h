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
// Author: Jens Chluba (CITA, University of Toronto)
//
// first implementation: May 2012
// last modification: June 2012
//
//==================================================================================================

#ifndef SZ_INTEGRAL_5D_H
#define SZ_INTEGRAL_5D_H

using namespace std;

//==================================================================================================
//
// 5D integration carried out using Patterson scheme. The SZ signal is obtained in the cluster frame
//
// eps_Int : relative accuracy for numerical integration (lower than 10^-6 is hard to achieve)
//
//==================================================================================================
double compute_SZ_distortion_Patterson_5D(double x, 
                                          double The, double betac, double muc, 
                                          double eps_Int=1.0e-4);

#endif

//==================================================================================================
//==================================================================================================
