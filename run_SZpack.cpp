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
// last modification   : February 2023
//
//==================================================================================================

#include <iostream>
#include <string>
#include <cmath>
#include <vector>

//==================================================================================================
// required libs
//==================================================================================================
#include "SZpack.h"
#include "physical_consts.h"
#include "routines.h"
#include "global_functions.h"
#include "Parameters.h"
#include "ModeFunctions.h"

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


// An optimisation to allow for the parameters object to be implicitly called instead of explicitly 
// passing it to each function call.
extern Parameters parameters;

//==================================================================================================
//
// Main call of routine.
//
//==================================================================================================
int main(int narg, char *args[])
{
    //==============================================================================================
    // startup with different number of arguments
    //==============================================================================================
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

    //==============================================================================================
    // Change this to use the convention of Nozawa et al. 2006 for the optical depth. Interpretation 
    // of the SZ signal is different in this case. The convention of Chluba et al., 2012 is used  
    // otherwise, which allow a clean separation of kinematic and scattering effect.
    // UseNozawaConvention is defined in .StartUp/global_variables.cpp but can be replaced with a 
    // boolean as wanted.
    //==============================================================================================
    setConvention(UseNozawaConvention);

    //==============================================================================================
    // Call in the different runmodes
    // This is a very large switch allowing for run_SZ_pack to calculate 
    // a number of different objects without too much thought. i.e., 
    // calling most of these will calculate the SZ effect from the 
    // parameters defined. 
    //==============================================================================================
    // Calculating the signal by a full 5-dimensional integral.
    // This is the slowest setting, and even slower with runmode of "full". It is also the most 
    // 'accurate', however the accuracy of the other methods are fully documented. Normally this is 
    // not the most relevant mode to use. See //TODO: for more details.
    //==============================================================================================
    if(parameters.mode=="5D")
    {
        output_SZ_distortion(modes.Int5D);
    }
    
    //==============================================================================================
    // Calculating the signal by precomputing two of the integrals analytically, under an assumption 
    // of low values for betac. This computes a 3-dimensional integral.
    //==============================================================================================
    else if(parameters.mode=="3D")
    {
        output_SZ_distortion(modes.Int3D);
    }
    
    //==============================================================================================
    // An asymptotic expansion around Te. This setting is best for Te < 2 keV. 
    // It will throw an error for Te > 20 keV.
    //==============================================================================================
    else if(parameters.mode=="ASYM")
    {
        output_SZ_distortion(modes.Asymptotic);
    }
    
    //==============================================================================================
    // A pivot based interpolation scheme described in Chluba et al. 2012.
    // This method will only calculate for 2 keV < Te < 75 keV. 
    //==============================================================================================
    else if(parameters.mode=="CNSN")
    {
        output_SZ_distortion(modes.CNSN);
    }
    
    //==============================================================================================
    // A slightly more involved form of the previous scheme, summarised in Chluba et al. 2013.
    // This should give very precise results for 2keV < Te < 75keV.
    //==============================================================================================
    else if(parameters.mode=="CNSNopt")
    {
        output_SZ_distortion(modes.CNSNopt);
    }
    //==============================================================================================
    // This setting will automatically calculate using the asymptotic mode for Te < 2 keV and CNSN 
    // above this. This function should give very precise results for Te < 75keV.
    //==============================================================================================
    else if(parameters.mode=="COMBO")
    {
        output_SZ_distortion(modes.Combo);
    }


    //==============================================================================================
    // This is the combo method, but above 75 keV it will default to the 3D method.
    //==============================================================================================
    else if(parameters.mode=="PRECISE")
    {
        output_SZ_distortion(modes.Precise);
    }
        
    //==============================================================================================
    // This is an example of a function to calculate the accuracy of a given method for calculating
    // the SZ effect. This compares to the precise method to determine the maximal absolute 
    // deviation between it and the given method. 
    //==============================================================================================
    else if(parameters.mode=="ACCURACY")
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
            // This function is declared in output.h and calculates and outputs the required value
            compute_precision_of_basis(Te, Method, temp, ofile, DI_s);
        }
        
        ofile.close();
    }

    //==============================================================================================
    // Calculating the signal using moments to account for variations of temperature and velocity 
    // within the line of sight. Described in detail in Chluba et al. 2013. i.e., the means method.
    //==============================================================================================
    else if(parameters.mode=="MEANS")
    {
        output_SZ_distortion(modes.Means);
    }

    //==============================================================================================
    // This is the means method as before, but now the weighting used in determining these values 
    // is y-weighted, rather than Dtau-weighted. See Lee et al. 2020 for more details.
    //==============================================================================================
    else if(parameters.mode=="MEANSYW")
    {
        output_SZ_distortion(modes.Means_Yweighted);
    }

    //==============================================================================================
    //TODO: fix this 
    //==============================================================================================
    else if(parameters.mode=="Kernel")
    {
        //output_SZ_distortion(modes.Kernel);
        //TODO: Sort this
    }

    //==============================================================================================
    // An example of a function using the null method to find the change in null crossing point for 
    // a variety of values of temperature, betac, omega, sigma, kappa and betac2_perp.
    //==============================================================================================
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
        //Setting a variety of temperatures to find the nulls in
        temp.T = TemperatureIterators(1.0, 70.0, 200, 1); //Change 1 to 0 to have a linear grid in T

        // Overwriting a number of the parameters set to be the specific values wanted in the function.
        temp.Dtau = 1.0;
        temp.betac = temp.betao = temp.muc = temp.muo = 0.0;
        temp.setCalcValues();
        temp.calc.betac_para = temp.means.Omega = temp.means.Sigma = temp.means.kappa = temp.calc.betac2_perp = 0.0;

        // Setting example arrays to iterate over, we will find the nulls for each of these values at
        // each of the temperatures defined above
        vector<double> betac_para_array = {0.001, -0.001, 0.005, -0.005, 0.01, -0.01};
        vector<double> omega_array = {0.01, 0.05, 0.1, 0.2, 0.4};
        vector<double> sigma_array = {0.001, -0.001, 0.002, -0.002, 0.005, -0.005, 0.01, -0.01};
        vector<double> kappa_array = {0.0001};
        vector<double> betac2_perp_array = {0.0001}; 

        for(int k = 0; k < temp.T.np; k++)
        {
            temp.updateT(temp.T.Trange[k]);
            double x0=compute_null_of_SZ_signal(temp);
            
            ofile << temp.Te << " " << x0 << ":" << endl;
            
            // This function is defined in output.h and calculates the change in null for each value 
            // in the given array, for the defined parameter being changed.
            IterateForNull(ofile, betac_para_array,
                [temp](const double &betac_para) mutable -> Parameters& {
                    temp.calc.betac_para = betac_para; return temp;
                }, x0);
            IterateForNull(ofile, omega_array,
                [temp](const double &Omega) mutable -> Parameters& {
                    temp.means.Omega = Omega; return temp;
                }, x0);
            IterateForNull(ofile, sigma_array,
                [temp](const double &Sigma) mutable -> Parameters& {
                    temp.means.Sigma = Sigma; return temp;
                }, x0);
            IterateForNull(ofile, kappa_array,
                [temp](const double &kappa) mutable -> Parameters& {
                    temp.means.kappa = kappa; return temp;
                }, x0);
            IterateForNull(ofile, betac2_perp_array,
                [temp](const double &betac2_perp) mutable -> Parameters& {
                    temp.calc.betac2_perp = betac2_perp; return temp;
                }, x0);
            
            // Calculating the derivatives defined in header
            double dI = Dcompute_signal_combo_for_x(x0, temp, 1);
            double d2I = Dcompute_signal_combo_for_x(x0, temp, 2);

            ofile << -d2I/dI << endl;
            
            // A terminal output to keep track of the progress of the function
            cout << temp.Te << endl;
        }
        
        ofile.close();
    }

    //==============================================================================================
    // This is an example function to computes various derivatives of the SZ signal.
    // The methods calculate the derivatives in the CMB frame, with respect to The, betac_para and 
    // betac2_perp. Derivatives can also be calculated with respect to x. 
    // Most of the function details are hidden within output_derivatives while this function 
    // calculates and outputs the derivatives:
    // 1/(i!j!k!)*(d^i/dThe^i)(d^j/dbetac_para^j)(d^k/dbetac2_perp^k)(signal) for i<=4, j<=2, k<=1
    // (where i=j=k=0 returns the combo signal itself)
    // It then additionally computes d^n/dx^n (signal) for n<=2 (n=0 returns the combo signal)
    //==============================================================================================
    else if(parameters.mode=="DERIVS")
    {
        int maxderiv=4; //i.e., the maximum value of i
        
        // This outputs the derivatives for a selection of temperatures.
        // output_derivatives can be found in output.h/output.cpp
        output_derivatives(5, maxderiv);
        output_derivatives(8, maxderiv);
        output_derivatives(10, maxderiv);
        output_derivatives(20, maxderiv);
        output_derivatives(30, maxderiv);
        output_derivatives(40, maxderiv);
        output_derivatives(50, maxderiv);
    }
    
    //==============================================================================================
    // This is an example function to calculate the diffeerences between the means method and the 
    // combo method variously.
    //==============================================================================================
    else if(parameters.mode=="CORR")
    {
        // The means method in question
        output_SZ_distortion(modes.Means);
        // This function calculates the relativistic corrections. i.e., the differences between the 
        // non-relativistic method and the combo method, in the CMB frame.
        output_SZ_distortion(modes.RelativisticCorrections);
        // This calculates the difference between the means signal given the omega, sigma and kappa 
        // values that have been set (i.e., the means method), and the signal without these 
        // variables (the combo method), all calculated in the CMB frame.
        output_SZ_distortion(modes.TemperatureDispersion);
    }
    
    //==============================================================================================
    // An example function using the two temperatures method to compare to the means method
    //==============================================================================================
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
    
    //==============================================================================================
    // An example function to calculate a number nT of normalised SZ signals for temperatures 
    // Te+n*DTe
    //==============================================================================================
    else if(parameters.mode=="AVYRELSZ")
    {
        string fname="./outputs.normalised/SZ_normalized.dat";
        ofstream ofile(fname.c_str());
        ofile.precision(16);
        
        // Setting up the specific values that are used in the function
        double Te=0.001, DTe=3.0;
        int nT=10, nfreq=1000;
        double xmin = 0.01, xmax = 25.0;
    
        // Specific y parameter for outputting in a normalised manner
        double y_param = 6.493939*1.0e-2;
        
        // vectors to be filled when iterated over
        vector<double> norm(nT);
        vector<vector<double> > SZsig(nT);
        
        // reset parameters to match what is wanted.
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
        vector<double> nus = temp.get_nucmb();
        
        ofile << "DI*y/norm ";
        // compute signals for different temperatures
        for(int l=0; l<nT; l++){
            compute_signal_combo(SZsig[l], temp, true, true);
            
            temp.updateT(temp.Te+DTe);
            ofile << temp.Te << " ";
        }
        ofile << endl;
        
        // compute normalizations
        for(int l=0; l<10; l++) { //linear integration approximation
            norm[l]=0.0;
            for(int k=0; k<nfreq; k++){
                norm[l]+=SZsig[l][k];
            }
            norm[l]*=(temp.xcmb[2]-temp.xcmb[1]);
        }

        // output to file
        for(int k=0; k<nfreq; k++) {
            ofile << nus[k] << " ";
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
