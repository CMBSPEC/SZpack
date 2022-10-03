//======================================================================================
// Author: Elizabeth Lee
// Based on work by Jens Chluba 
//
// //TODO: Check that these are doing these things in sensible ways!
//======================================================================================

#ifndef SZPack_MultipleScattering_H
#define SZPack_MultipleScattering_H

#include <vector>
#include <cmath>
#include "physical_consts.h"
#include "routines.h"
#include "Patterson.h"

#include "SZ_Integral_multiple.h"

using namespace std;

//==================================================================================================
//
// const density sphere
//
//==================================================================================================

namespace constDensitySphere{
    double tau0(double B0);

    double tau_hat(double Z0, int L); //Tau averaged over mu_r

    double tau_av(double B0, int L); // Tau averaged over mu_r and s
}

//==================================================================================================
//
// isothermal beta profile
//
//==================================================================================================

namespace isothermalBeta{
    void setBeta(double Beta);
    
    double Ne(double r2);

    double tau0(double b);

    //==================================================================================================
    double s_Int(double Z0, double Mur); //Tau averaged over s

    double mur_Int_explicit(double Z0, int L); //Tau explicitly averaged over mu and s

    //==================================================================================================
    double mur_Int(double Z0, int L); //Tau averaged over mu and s splined between calculated values

    //==================================================================================================
    double z0_Int(double B0, int L); //Tau averaged over mu, s and z.

    double z0_Int_m1(double B0, int L); //Tau averaged over mu, s and z using the second multipole.

    //==================================================================================================
    double tau_av(double B0, int L); //Tau averaged over mu, s and z, with prefactor

    //==================================================================================================
    double lambda_function(double B0, int L);

    double kappa_function(double B0, int L);
}
#endif
