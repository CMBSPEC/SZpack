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
#include "ModeFunctions.h"
#include "SZpack.h"

#include "SZ_moment_method.h"
#include "SZ_cluster_profiles.h"
#include "ClusterFunctions.h"

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
#include "output_moment_method.h"

extern Parameters parameters;

//==================================================================================================
//
// examples from Vikhlinin et al. 2006
//
//==================================================================================================
Cluster_param_Ne A478N = {10.170, 155.5, 2928.9, 1.254, 0.704, 5.000, 0.762, 23.84, 1.000}; 
Cluster_param_Te A478T = {11.06, 0.27, 0.02, 5.00, 0.4, 0.38, 129, 1.60};

Cluster_param_Ne A2029N = {15.721, 84.2, 908.9, 1.164, 0.545, 1.669, 3.510, 5.00, 1.000}; 
Cluster_param_Te A2029T = {16.19, 3.04, -0.03, 1.57, 5.9, 0.1, 93, 0.48}; 

Cluster_param_Ne A2390N = {3.605, 308.2, 1200.0, 1.891, 0.658, 0.563, 0, 1, 1};
Cluster_param_Te A2390T = {19.34, 2.46, -0.1, 5.00, 10.0, 0.12, 214, 0.08}; 

//==================================================================================================
//
// main call of routine. Here parameters are handled.
//
//==================================================================================================
int main(int narg, char *args[])
{
    parameters.mode = "ISO";
    
    //===============================================================================
    // startup with different number of arguments
    //===============================================================================
    if(narg==1){} //TODO: Fix this so it is less yuck.
    
    else if(narg==2) parameters.mode=args[narg-1];
    
    else if(narg==3)
    {
        parameters.mode=args[narg-2];
        string filename=args[narg-1];
        
        print_message("Reading parameters from parameter file: " + filename);
        
        read_startup_data_parser(filename, parameters);
    }

    else{ exit_error("Too many/few parameters"); }
        
    distortionModes modes = distortionModes();

    //===============================================================================
    // combination of asymptotic expansion + basis functions of CNSN2012
    // This function should give very precise results for Te < 75keV and
    // can be used for comparison.
    //===============================================================================
    if(parameters.mode=="COMBO") //TODO: Switch switch switch!
    {
        output_SZ_distortion(modes.Combo);
        return 0;
    }
    
    //===============================================================================
    // check the precision of the basis function. Settings should be 
    // changed in compute_precision_of_basis() above
    //===============================================================================
    else if(parameters.mode=="ACC")
    {
        //compute_precision_of_basis(modes.CNSNopt); //TODO: fix this
        print_error("Currently broken.");
        return 0;
    }
    
    //===============================================================================
    // combination of asymptotic expansion + basis functions of CNSN2012
    // This function should give very precise results for Te < 75keV and
    // can be used for comparison.
    //===============================================================================
    else if(parameters.mode=="MEANS")
    {
        parameters.fileEnding = ".III" + parameters.fileEnding;
        output_SZ_distortion(modes.Means);
        return 0;
    }

    //===============================================================================
    //
    // matrix setup for moment method
    //
    //===============================================================================
    if(parameters.mode=="ISO" || parameters.mode=="ISOOPT" || parameters.mode=="VIK" || parameters.mode=="VIKOPT") 
        parameters.fileEnding=".acc_"+to_string(parameters.accuracy_level)+"."+parameters.mode+parameters.fileEnding;
    
    //===============================================================================
    // This version uses the accuracy settings according to CSNN 2012
    // for given kmax==Te_order
    //===============================================================================
    SZ_moment_method SZMoments = SZ_moment_method(parameters, true);
    
    //===============================================================================
    // In this case the optimal kmax is determined internally according 
    // to the required accuracy goal and maximal temperature. This 
    // version of the moment matrix minimizes the # of parameters needed.
    // Calling the functions of 'SZMomentsopt' might give exactly the  
    // same as those of 'SZMoments'! This depends on the settings.
    //===============================================================================
    SZ_moment_method SZMomentsopt = SZ_moment_method(parameters, false);

    //===============================================================================
    //
    // store the matrices corresponding to the different moment versions
    // of variables for later computations. These are actually all one 
    // needs to compute S = F * m !
    //
    //===============================================================================
    if(parameters.mode=="MATRIX")
    {
        print_message("Storing SZ signal matrices");
        
        export_vector(parameters.xcmb, "x", parameters);

        export_moment_matrix("M", SZMoments.M, parameters, SZMoments);
        export_moment_matrix("MT", SZMoments.MT, parameters, SZMoments);

        export_moment_matrix("M.opt", SZMomentsopt.M, parameters, SZMomentsopt);
        export_moment_matrix("MT.opt", SZMomentsopt.MT, parameters, SZMomentsopt);

        return 0;
    }
    
    //===============================================================================
    // call in the different runmodes
    //===============================================================================
    if(parameters.mode=="ISO")
    {
        output_distortion_moments_isothermal(SZMoments);
        return 0;
    }
    
    else if(parameters.mode=="ISOOPT")
    {
        parameters.fileEnding=".opt"+parameters.fileEnding;
        output_distortion_moments_isothermal(SZMomentsopt, " (opt)");
        return 0;
    }

    //==================================================================
    // compute null of SZ signal
    //==================================================================
    else if(parameters.mode=="NULL")
    {
        //output_SZ_null(); //TODO: fix this
        print_error("Currently broken.");
        return 0;
    }
    
    //===============================================================================
    //
    // Setting cluster parameters. Examples defined above.
    //
    //===============================================================================
    //SZ_cluster_profiles CL(A478N, A478T, "A478");
    SZ_cluster_profiles CL(A2029N, A2029T, "A2029");
    //SZ_cluster_profiles CL(A2390N, A2390T, "A2390");
    
    CL.show_cluster_parameters();
    parameters.fileEnding="."+CL.Get_label()+parameters.fileEnding;
    
    //===============================================================================
    // compute yk temperature moments for some slices. 
    //===============================================================================
    if(parameters.mode=="VIKMOM")
    {
        print_message("Computing yk-integrals");
        parameters.kmax = 10;
        parameters.gridpoints = 200;
        output_distortion_cluster_yk(0.0001/CL.Get_rc(), CL, parameters);
        parameters.gridpoints = 100;
        output_distortion_cluster_yk(10.0/CL.Get_rc(), CL, parameters);
        output_distortion_cluster_yk(50.0/CL.Get_rc(), CL, parameters);
        output_distortion_cluster_yk(100.0/CL.Get_rc(), CL, parameters);
        return 0;
    }
    
    //===============================================================================
    // compute derivatives of SZ signal 
    //===============================================================================
    if(parameters.mode=="DERIVS")
    {
        parameters.Dtau=1.0;
        parameters.betac=0.0;
        int maxderiv = 4;
        
        output_derivatives(10, maxderiv);
        output_derivatives(20, maxderiv);
        output_derivatives(30, maxderiv);
        output_derivatives(40, maxderiv);
        output_derivatives(50, maxderiv);

        return 0;
    }
    
    //===============================================================================
    // SZ signal for two-temperature case
    //===============================================================================
    if(parameters.mode=="TWOT")
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
    
    //===============================================================================
    // compute derivatives of SZ signal 
    //===============================================================================
    if(parameters.mode=="DEGEN")
    {
        parameters.xmin = 0.1; //This overwrites user input values so possibly shouldn't be here
        parameters.xmax = 30.0; 
        parameters.gridpoints = 200;
        parameters.PopulateRunningArrays();
        compute_degeneracy_functions();
        return 0;
    }

    //===============================================================================
    // chose particular line-of-sight through the cluster. xs and ys are in units of
    // the clusters core radius.
    //===============================================================================
    double xs=0.0;
//    double ys=500.0/CL.Get_rc();
//    add=".y_500kpc"+add; 
    double ys=1.0;
    parameters.fileEnding=".y_1rc.S1"+parameters.fileEnding; 
    
/*    output_TSZ_bestfit(outpath+"SZ_moments_", add, xcmb, CL, SZMoments, 
                       betac, muc, 10.0, 40);
    
    exit(0);
*/  
    //===============================================================================
    // compute SZ signal for different cluster profile according to 
    // fits of Vikhlinin et al 2006
    //===============================================================================
    if(parameters.mode=="VIK")
    {//TODO: CNSN expansion always returning 0. Need to check this is doing the right thing!!!!
        output_distortion_moments_cluster_VIK("SZ_moments_clusters", CL, xs, ys, SZMoments, parameters);
    }

    else if(parameters.mode=="VIKOPT")
    {
        parameters.fileEnding=".opt"+parameters.fileEnding;
        output_distortion_moments_cluster_VIK("SZ_moments_clusters", CL, xs, ys, SZMomentsopt, parameters, " (opt)");
    }

    //===============================================================================
    // computing SZ signal by explicitly integrating along the line of sight for each
    // frequency bin. This routine is just for comparisons & rather slow.
    //===============================================================================
    else if(parameters.mode=="VIKEXP")
    {
        parameters.fileEnding=".explicit"+parameters.fileEnding;
        output_distortion_cluster_VIK_explicit("SZ_moments_clusters", CL, xs, ys, parameters);
    }
    
    return 0;
}

//==================================================================================================
//==================================================================================================
