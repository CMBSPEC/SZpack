//==================================================================================================
//
// run_SZ_moment_method routine
// 
//==================================================================================================
//
// purpose: give some explicit examples of how to use the SZpack temperature-velocity moment method
// to calculate the SZ signals for simple cluster profiles. The code can be run with default 
// parameters that are set in the function
//
//      int main(int narg, char *args[])
//
// below. Alternatively one can call it with a start file.
//
//==================================================================================================
//
// Author: Jens Chluba (CITA, University of Toronto)
//
// first implementation: July 2012
// last modification   : July 2012
//
//==================================================================================================

//==================================================================================================
// Standards
//==================================================================================================
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include <iomanip>
#include <cmath>
#include <limits.h>
#include <vector>

//==================================================================================================
// required libs
//==================================================================================================
#include "physical_consts.h"
#include "Integration_routines.h"
#include "routines.h"
#include "nPl_derivatives.h"
#include "parser.h"
#include "global_functions.h"
#include "Parameters.h"
#include "SZpack.h"
#include "SZpackMultipleScattering.h"

#include "gsl/gsl_multimin.h"
#include <gsl/gsl_linalg.h>

//==================================================================================================
//
// namespaces
//
//==================================================================================================
using namespace std;

//==================================================================================================
// these are several Modules (i.e. plugin) that are directly loaded when compiling. If there are 
// changes to these files one should usually type "make clean" before "make"
//==================================================================================================
#include "./StartUp/global_variables.cpp"
#include "./StartUp/initial_setup.cpp"
#include "output.h"
#include "output_multiple_scattering.h"

extern Parameters parameters;

//==================================================================================================
//
// main call of routine. Here parameters are handled.
//
//==================================================================================================
int main(int narg, char *args[])
{
    parameters.gridpoints = 400;
    parameters.xmin = 1.0e-1;
    parameters.xmax = 40.0;
    parameters.PopulateRunningArrays();
    parameters.Dtau = 0.01;
    parameters.betac = parameters.muc = parameters.betao = parameters.muo = 0;

    double Trange[] = {5.0, 10.0, 15.0, 25.0};
    vector<double> temperature_range(begin(Trange), end(Trange));
    parameters.T = TemperatureIterators(temperature_range);

    for (int i = 0; i < parameters.T.np; i++){
        parameters.Te = parameters.T.Trange[i];
        output_SZ_signal_l(parameters);
    }
    //exit(0);
    
    for (int i = 0; i < parameters.T.np; i++){
        parameters.Te = parameters.T.Trange[i];
        output_DI2_E_temperature(parameters); //TODO: Rename this function
    }
    //exit(0);
    
    parameters.Te = 0.01*const_me;
    isothermalBeta::setBeta(2.0/3.0);
    output_lowest_order_signal(parameters, 0);
    //exit(0);

    output_tau_l_sphere(12, 400);
    //exit(0);
    
    output_tau_l_beta(30, 50, 2.0/3.0); //TODO: This is very slow
    //exit(0);
    
    parameters.gridpoints = 50; //TODO: set up a better initialisation of this probably
    parameters.xmin = 1.0e-1; //TODO: I also think maybe this should be a different set of running rather than xcmb tbh
    parameters.xmax = 1.0e+2;
    parameters.PopulateRunningArrays();

    isothermalBeta::setBeta(2.0/3.0);

    output_emission_absorption(parameters);

    output_emission_absorption_kinetic(parameters);
    //exit(0);
/*    
    int lmax=100;
    double beta=2.0/3.0, r=0.0;
    
    vector<double> taulc(lmax+2);
    vector<double> taulv(lmax+2);
    double b=1.0;
    
    for(int k=0; k<=lmax; k+=2)
    {
        taulc[k]=z0_Int(b, beta, k)/(0.5*pow(tau0(beta, b), 2));
        r+=(2.0*k+1.0)*taulc[k];
        cout << k << " " << (2.0*k+1.0)*taulc[k] << " " << r << endl;
    }
    exit(0);
*/    
    parameters.gridpoints = 400;
    parameters.xmin = 1.0e-1;
    parameters.xmax = 40.0;
    parameters.PopulateRunningArrays();

    output_CMB_isotropy_Y_functions(parameters);
    exit(0);

    parameters.Te = 20.0;
    output_CMB_isotropy_signals(parameters);
    exit(0);
    
    //TODO: Clean from here down into actual functions and such and such!
    
    //==============================================================================================
    // second scattering to thermal SZ and CMB scattering
    //==============================================================================================
    parameters.gridpoints = 400;
    parameters.xmin = 1.0e-1;
    parameters.xmax = 40.0;
    parameters.PopulateRunningArrays();
    parameters.Te = 50.0;
    parameters.setCalcValues();
    
    for(int l=0; l<=4; l++) cout << P_l_Kernel_Norm(l, 0.01, 1e-5) << endl;
    
    cout << P_l_Kernel(4, 0, 0.01) << endl;
        
    ofstream ofileK("./Kernel.dat");
    ofileK.precision(8);
    
    double sig[4];
    
    sig[0]=P_l_Kernel(0, 0, parameters.calc.The);
    sig[1]=P_l_Kernel(1, 0, parameters.calc.The);
    sig[2]=P_l_Kernel(2, 0, parameters.calc.The);
    sig[3]=P_l_Kernel(3, 0, parameters.calc.The);
    
    double sFWHM = 2.0*atanh(sqrt(2.0*log(2.0)*parameters.calc.The));

    parameters.xmin = -4.0*sFWHM;
    parameters.xmax =  4.0*sFWHM;
    parameters.gridpoints ++;
    //TODO: Make a function handle this!
    init_xarr(parameters.xmin, parameters.xmax, &parameters.xcmb[0], parameters.gridpoints, 0, 0);

    for(int m=0; m<0*parameters.gridpoints; m++)
    {
        cout << parameters.xcmb[m] << endl;
        
        ofileK << parameters.xcmb[m] << " ";
        
        for(int l=0; l<=4; l++)
            ofileK << P_l_Kernel(l, parameters.xcmb[m], parameters.calc.The) << " "
                   << Compute_Kernel_interpol(l, parameters.xcmb[m], parameters.calc.The) << " ";

 /*       ofileK << xc[m]/sFWHM << " ";
        
        for(int l=0; l<=4; l++)
            ofileK << P_l_Kernel(l, xc[m], Thev)/sig[l] << " "
                   << Compute_Kernel_interpol(l, xc[m], Thev)/sig[l] << " ";
 */       
        ofileK << endl;
    }
    ofileK.close();
    
//    exit(0);

    parameters.xmin = 1.0e-1;
    parameters.xmax =  40.0;
    parameters.PopulateRunningArrays();

    parameters.Te = 25.0;
    parameters.Dtau = 0.01;
    parameters.betac = parameters.muc = parameters.betao = parameters.muo = 0;
    parameters.setCalcValues();

    ofstream ofile("./CMB.data.25keV.compare.dat");
    ofile.precision(8);
    
    Parameters temp = Parameters();
    temp.copyParameters(parameters);
    temp.Dtau = 1.0;
    temp.setCalcValues();

    for(int m=0; m<parameters.gridpoints; m++)
    {
        cout << parameters.xcmb[m] << endl;
        
        ofile << parameters.xcmb[m] << " ";
        
//        for(int l=0; l<=3; l++)
//            ofile << pow(xc[m], 3)*compute_SZ_distortion_Patterson_multiple(xc[m], l, Te/const_me, "SEC") << " ";

/*        for(int l=0; l<=3; l++)
            ofile << pow(xc[m], 3)*compute_SZ_distortion_Patterson_multiple_Kernel(xc[m], l, Te/const_me) << " ";
        
        for(int l=0; l<=3; l++)
            ofile << pow(xc[m], 3)*compute_SZ_distortion_asym_multiple(xc[m], l, Te/const_me, 10) << " ";
*/        
        double Dn0=compute_signal_combo(m, parameters);
        double x3=pow(parameters.xcmb[m], 3), fac=parameters.rare.Dn_DI_conversion()*parameters.Dtau*parameters.Dtau/2.0;
        
        ofile << parameters.rare.Dn_DI_conversion()*x3*Dn0 << " ";
        
        ofile << fac*x3*compute_SZ_distortion_asym_multiple(parameters.xcmb[m], 0, parameters.calc.The, 10) << " ";

        ofile << fac*x3*(
                     0.61*compute_SZ_distortion_asym_multiple(parameters.xcmb[m], 0, parameters.calc.The, 10)
                    +0.087/10.0*compute_SZ_distortion_asym_multiple(parameters.xcmb[m], 2, parameters.calc.The, 10)
                  ) << " ";

        ofile << parameters.rare.Dn_DI_conversion()*x3*parameters.Dtau/2.0*(0.61+0.087/10.0-1.0)*Dn0 << " ";
        
        double Dsig2=compute_SZ_distortion_Patterson_multiple(1.0, 2, parameters.calc.The, "SIG")-0.1;
        
        ofile << fac*x3*Dsig2*Dn0 << " ";
 
        double S2_0=compute_SZ_distortion_Patterson_multiple_Kernel_2D(parameters.xcmb[m], 0, parameters.calc.The);
        double S2_2=compute_SZ_distortion_Patterson_multiple_Kernel_2D(parameters.xcmb[m], 2, parameters.calc.The);
        
        ofile << fac*x3*S2_0 << " ";

        ofile << fac*x3*(0.61*S2_0+0.087/10.0*S2_2) << " ";
        
        ofile << fac*x3*(0.61*S2_0+0.087/10.0*S2_2+Dsig2*Dn0)
                 +parameters.rare.Dn_DI_conversion()*x3*parameters.Dtau/2.0*(0.61+0.087/10.0-1.0)*Dn0 << " ";
        
        //------------------------------------------------------------
        // x-derivatives
        //------------------------------------------------------------
        
        double Ix =compute_signal_combo(m, temp);
        double xp = temp.xcmb[m]*(1.0+0.001), xm = temp.xcmb[m]*(1.0-0.001);
        temp.xcmb[m] = xp;
        double Ixp=compute_signal_combo(m, temp);
        temp.xcmb[m] = xm;
        double Ixm=compute_signal_combo(m, temp);
        double dI =(Ixp-Ixm)/(2.0*0.001);
        double d2I=(Ixp-2.0*Ix+Ixm)/pow(0.001, 2);

        ofile << fac*x3*temp.calc.The*parameters.xcmb[m]*(4.0*dI+parameters.xcmb[m]*d2I)*0.62 << endl;
    }
    
    ofile.close();
    
    cout << parameters.rare.Dn_DI_conversion() << endl;
    
    return 0;

/*
    //==============================================================================================
    // compute scattering cross section correction
    //==============================================================================================
    int np=50;
    vector<double> Thev(np);
    init_xarr(1.0e-3, 1.0, &Thev[0], np, 1, 0);
   
    ofstream ofile("./sig.data.dat");
    ofile.precision(8);
    
    for(int m=0; m<np; m++)
    {
        ofile << Thev[m] << " ";
        
        for(int l=0; l<=4; l++)
            ofile << compute_SZ_distortion_Patterson_multiple(1.0, l, Thev[m], "SIG") << " ";
        
        ofile << endl;
    }
    
    ofile.close();
*/
    return 0;
}

//==================================================================================================
//==================================================================================================
