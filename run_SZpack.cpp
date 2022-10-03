//==================================================================================================
//
// run_SZpack routine
// 
//==================================================================================================
//
// purpose: give some explicit examples of how to use SZpack. The code can be run with default 
// parameters that are set in the function
//
//      int main(int narg, char *args[])
//
// below. Alternatively one can call it with a start file. Not all functions of SZpack are called
// explicitly. For more options see 'SZpack.h'
//
//==================================================================================================
//
// Author: Jens Chluba & Elizabeth Lee
//
// first implementation: May 2012
// last modification   : March 2020
//
//==================================================================================================
// 28th Aug,  2017: added y-weighted moment method for temperature corrections
// 20th Dec,  2012: added runmodes CNSNopt and cleaned the routines up.
//  5th Nov,  2012: added examples for derivative outputs and to illustrate the
//                  temperature-velocity moment method for smooth cluster profiles.
// 10th Aug,  2012: output of distortion now also in MJy / sr
// 23th July, 2012: added function to compute the precision different expansions
// 22th July, 2012: added combo function option

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
#include "SZpack.h"
#include "physical_consts.h"
#include "routines.h"
#include "parser.h"
#include "global_functions.h"
#include "Parameters.h"

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
 
extern Parameters parameters;

//==================================================================================================
//
// Main call of routine.
//
//==================================================================================================
int main(int narg, char *args[])
{
    //==================================================================
    // startup with different number of arguments
    //==================================================================
    if(narg==1){}
    
    else if(narg==2) parameters.mode=args[narg-1];
    
    else if(narg==3)
    {
        parameters.mode=args[narg-2];
        string filename=args[narg-1];
        
        print_message("Reading parameters from parameter file: " + filename);
        read_startup_data_parser(filename, parameters);
    }

    else{ exit_error(" Too many/few parameters "); }

    distortionModes modes = distortionModes();

    //==================================================================
    // change this to use the convention of Nozawa et al. 2006 for 
    // the optical depth. Interpretation of the SZ signal is different
    // in this case. The convention of Chluba et al., 2012 is used  
    // otherwise, which allow a clean separation of kinematic and 
    // scattering effect.
    // UseNozawaConvention is defined in .StartUp/global_variables.cpp 
    // but could just be replaced with a boolean as wanted.
    //==================================================================
    setConvention(UseNozawaConvention);

    //==================================================================
    // call in the different runmodes
    //==================================================================
    if(parameters.mode=="5D")
    {
        output_SZ_distortion(modes.Int5D);
    }
    
    else if(parameters.mode=="3D")
    {
        output_SZ_distortion(modes.Int3D);
    }

    else if(parameters.mode=="Kernel")
    {
        //output_SZ_distortion(modes.Kernel);
        //TODO: Sort this
    }
    
    else if(parameters.mode=="ASYM")
    {
        output_SZ_distortion(modes.Asymptotic);
    }
    
    //==================================================================
    // Added the basis at additional temperature pivots. This function 
    // should give very precise results for 2keV < Te < 75keV
    //==================================================================
    else if(parameters.mode=="CNSN")
    {
        output_SZ_distortion(modes.CNSN);
    }
    
    //==================================================================
    // Added the basis at additional temperature pivots. This function
    // should give very precise results for 2keV < Te < 75keV
    //==================================================================
    else if(parameters.mode=="CNSNopt")
    {
        output_SZ_distortion(modes.CNSNopt);
    }
    //==================================================================
    // combination of asymptotic expansion + basis functions of CNSN2012
    // This function should give very precise results for Te < 75keV and
    // can be used for comparison.
    //==================================================================
    else if(parameters.mode=="COMBO")
    {
        output_SZ_distortion(modes.Combo);
    }
        
    //==================================================================
    // check the precision of the basis function. 
    //==================================================================
    else if(parameters.mode=="ACC")
    {
        double DI_s = 4.881e-5; // fiducial accuracy according to CSNN 2012

        method Method = modes.Asymptotic; // This can take any of the modes and check its accuracy!
        Parameters temp = Parameters();
        temp.copyParameters(parameters);
        temp.fileEnding = "_Accuracy"+temp.fileEnding;

        ofstream ofile;
        SetUpOutput("Computing the precision of ", Method, temp, ofile);
        ofile << "#\n# Output format: Te [K] | Max error/fiducial accuracy " << endl;

        temp.gridpoints = 300; //This overwrites the user set value to ensure better accuracy
        temp.PopulateRunningArrays();

        temp.T = TemperatureIterators(1.0, 15.0, 50); //originally 1, 4, 40 in moment method
        
        for(double Te : temp.T.Trange){
            compute_precision_of_basis(Te, Method, temp, ofile, DI_s);
        }
        
        ofile.close();
    }

    //==================================================================
    // expansion of SZ signal around mean values
    //==================================================================
    else if(parameters.mode=="MEANS")
    {
        output_SZ_distortion(modes.Means);
    }

    //==================================================================
    // expansion of SZ signal around mean values using y-weighted values
    //==================================================================
    else if(parameters.mode=="MEANSYW")
    {
        output_SZ_distortion(modes.Means_Yweighted);
    }

    //==================================================================
    // compute null of SZ signal. For details see function
    //==================================================================
    else if(parameters.mode=="NULL")
    {
        string fileAddition = "SZ";
        ofstream ofile;

        parameters.fileEnding = "_Null"+parameters.fileEnding;
        SetUpOutput("Computing null of SZ signal", fileAddition, parameters, ofile);
        ofile << "#\n# Output format: Te [K] | x0 = (h nu/k T0) :  x_null - x0 (for each element in array) ";
        ofile << "| -(d2/dx2(I) / d/dx(I)) at x0" << endl;

        Parameters temp = Parameters();
        temp.copyParameters(parameters);
        temp.T = TemperatureIterators(1.0, 70.0, 200, 1); //Change 1 to 0 to have a linear grid in T

        temp.Dtau = 1.0;
        temp.betac = temp.betao = temp.muc = temp.muo = 0.0;
        temp.setCalcValues();
        temp.calc.betac_para = temp.means.Omega = temp.means.Sigma = temp.means.kappa = temp.calc.betac2_perp = 0.0;

        vector<double> betac_para_array = {0.001, -0.001, 0.005, -0.005, 0.01, -0.01};
        vector<double> omega_array = {0.01, 0.05, 0.1, 0.2, 0.4};
        vector<double> sigma_array = {0.001, -0.001, 0.002, -0.002, 0.005, -0.005, 0.01, -0.01};
        vector<double> kappa_array = {0.0001};
        vector<double> betac2_perp_array = {0.0001}; //These are example values to iterate over!

        for(int k = 0; k < temp.T.np; k++)
        {
            temp.updateT(temp.T.Trange[k]);
            double x0=compute_null_of_SZ_signal(temp);
            
            ofile << temp.Te << " " << x0 << ":" << endl;
                //
            IterateForNull(ofile, betac_para_array, temp.calc.betac_para, temp, x0);
                //
            IterateForNull(ofile, omega_array, temp.means.Omega, temp, x0);
                //
            IterateForNull(ofile, sigma_array, temp.means.Sigma, temp, x0);
                //
            IterateForNull(ofile, kappa_array, temp.means.kappa, temp, x0);
                //
            IterateForNull(ofile, betac2_perp_array, temp.calc.betac2_perp, temp, x0);
            
            double dI = Dcompute_signal_combo_for_x(x0, temp, 1);
            double d2I = Dcompute_signal_combo_for_x(x0, temp, 2);

            ofile << -d2I/dI << endl;
            
            cout << temp.Te << endl;
        }
        
        ofile.close();
    }

    //==================================================================
    // compute derivatives of SZ signal 
    //==================================================================
    else if(parameters.mode=="DERIVS")
    {//outputs derivatives for a number of temperatures.
        int maxderiv=4;
        
        output_derivatives(5, maxderiv);
        output_derivatives(8, maxderiv);
        output_derivatives(10, maxderiv);
        output_derivatives(20, maxderiv);
        output_derivatives(30, maxderiv);
        output_derivatives(40, maxderiv);
        output_derivatives(50, maxderiv);
    }
    
    //==================================================================
    // expansion of SZ signal around mean values
    //==================================================================
    else if(parameters.mode=="RELCORR")
    {
        // for additional second order parameters see function itself
        output_SZ_distortion(modes.Means);
        output_SZ_distortion(modes.RelativisticCorrections);
        output_SZ_distortion(modes.TemperatureDispersion);
    }
    
    //==================================================================
    // SZ signal for two-temperature case
    //==================================================================
    else if(parameters.mode=="TWOT")
    {
        double ftau = 0.2;
        double DT_T = 0.5;

        Parameters fp = Parameters();
        fp.copyParameters(parameters);
        
        ofstream ofile;
        string fileAddition = "SZ_moments.TwoTemp.ftau_"+DoubletoString(ftau)+"-Te_"+DoubletoString(fp.Te)+"-Delta_"+DoubletoString(DT_T);
        SetUpOutput("Two-temperature case", fileAddition, fp, ofile);

        // Calculating the precise version
        vector<double> Dn;
        compute_signal_TwoTemperatures(Dn, ftau, DT_T, fp, true);

        // Setting up the details to calculate the same signal through the means method
        fp.updateT(fp.Te*(1.0+ftau*DT_T));

        double omega0 = ftau*(1.0-ftau)*pow(DT_T/(1.0+ftau*DT_T), 2);
        double omega1 = ftau*(1.0-ftau)*(1.0-2.0*ftau) * pow(DT_T/(1.0+ftau*DT_T), 3);
        fp.means.Omega = omega0;
        fp.means.assignOmegas(omega0, omega1, 0.0);

        cout << fp.Te << " " << fp.means.omegas[0] << " " << fp.means.omegas[1] << " " << (1.0-2.0*ftau)/ftau << endl;

        // Calculating the means and means expansion versions
        vector<double> Dn_means, Dn_ex;
        compute_signal_means(Dn_means, fp, true);
        compute_signal_means_ex(Dn_ex, fp, true);

        // Outputting all the details
        ofile << "#\n# Specific parameters:" << endl;
        ofile << "# ftau= " << ftau << " DT_T= " << DT_T << endl;
        ofile << "#\n# Calculated parameters:" << endl;
        ofile << "# averageT= " << fp.Te << " omega0= " << omega0 << " omega1= " << omega1 << endl;
        ofile << "#\n# Output format: x = (h nu/k T0) | nu [GHz] | DI(x) in MJy/sr - precise";
        ofile << " | DI(x) in MJy/sr - means | DI(x) in MJy/sr - means_ex " << endl;

        for(int k = 0; k < fp.gridpoints; k++){   
            ofile << fp.xcmb[k] << " " << fp.xcmb[k]*fp.rare.TCMB()/const_h_kb/1.0e+9<< " ";
            ofile << Dn[k] << " " << Dn_means[k] << " " << Dn_ex[k];
            ofile << endl;
        }
        ofile.close();
    }
    
    //==================================================================
    // average y relativistic corrections
    //==================================================================
    else if(parameters.mode=="AVYRELSZ")
    {
        //outputs a number nT of normalised SZ signals for temperatures Te+n*DTe * a multiplicative factor
        //fac = Io the observer intensity * chosen y parameter.
        string fname="./outputs.normalised/SZ_normalized.dat";
        ofstream ofile(fname.c_str());
        ofile.precision(16);
        
        double Te=0.001, DTe=3.0;
        int nT=10, nfreq=1000;
        double xmin = 0.01, xmax = 25.0;
    
        double xo=0.01760533;
        double y_param = 6.493939*1.0e-2;
        
        vector<double> norm(nT);
        vector<double> dum(nfreq);
        vector<vector<double> > SZsig(nT);
        
        //reset parameters to match what is wanted.
        Parameters temp = Parameters();
        temp.copyParameters(parameters);
        temp.Te = Te;
        temp.Dtau = 1.0;
        temp.betac = temp.betao = 0.0;
        temp.muc = temp.muo = 1.0;
        temp.setCalcValues();

        temp.xmin = xmin;
        temp.xmax = xmax;
        temp.gridpoints = nfreq;
        temp.xcmb.resize(nfreq);
        init_xarr(xmin, xmax, &temp.xcmb[0], nfreq, 0, 0); // linear frequency grid, to allow for easier normalisation.
        
        ofile << "DI*y/norm ";
        // compute signals for different temperatures
        for(int l=0; l<nT; l++){
            compute_signal_combo(SZsig[l], temp, true, true);
            
            temp.updateT(temp.Te+DTe);
            ofile << temp.Te << " ";
        }
        ofile << endl;
        
        // compute normalizations
        for(int l=0; l<10; l++) { //linear integration approximation.
            norm[l]=0.0;
            for(int k=0; k<nfreq; k++){
                norm[l]+=SZsig[l][k];
            }
            norm[l]*=(temp.xcmb[2]-temp.xcmb[1]);
        }

        // output to file
        for(int k=0; k<nfreq; k++) {
            ofile << temp.xcmb[k]/xo << " ";
            for(int l=0; l<10; l++){
                ofile << SZsig[l][k]*y_param/norm[l] << " ";
            }
            ofile << endl;
        }
        
        ofile.close();
    }

    //==================================================================
    else print_error("Unknown runmode: "+parameters.mode);
    
    return 0;
}

//==================================================================================================
//==================================================================================================
