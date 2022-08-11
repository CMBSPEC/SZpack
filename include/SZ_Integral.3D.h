//==================================================================================================
//
// Computation of the thermal and k-SZ effect using explicit integration of the collision term in
// different orders of beta. The cluster is assumed to be isothermal and moving at a given speed 
// betac with direction muc relative to the line of sight. The integrals are carried out in the 
// cluster frame. 
//
//==================================================================================================
//
// Author: Jens Chluba  (CITA, University of Toronto)
//
// first implementation: April 2012
// last modification   : July  2012
//
//==================================================================================================
// 8th July: changed definition of S^kin, so that the temperature-independent term is fully canceled

#ifndef SZ_INTEGRAL_3D_H
#define SZ_INTEGRAL_3D_H

using namespace std;

//==================================================================================================
//
// 3D integration carried out using Patterson scheme. The SZ signal is obtained in the cluster frame
//
// eps_Int : relative accuracy for numerical integration (lower than 10^-6 is hard to achieve)
//
// mode == "monopole"      --> only scattering of monopole without second order kinematic corr
// mode == "dipole"        --> only scattering of dipole     (first order kinematic correction)
// mode == "quadrupole"    --> only scattering of quadrupole (second order kinematic correction)
// mode == "monopole_corr" --> only scattering of second order kinematic correction to monopole
// mode == "all"           --> all terms added
// mode == "kin"           --> only kinematic terms
//
//==================================================================================================
double compute_SZ_distortion_Patterson_3D(double x, 
                                          double The, double betac, double muc, 
                                          string mode,
                                          double eps_Int=1.0e-4);

#endif

//==================================================================================================
//==================================================================================================
