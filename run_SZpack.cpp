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
// Author: Jens Chluba (CITA, University of Toronto and Johns Hopkins University)
//
// first implementation: May 2012
// last modification   : Dec 2012
//
//==================================================================================================
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

//==================================================================================================
//
// namespaces
//
//==================================================================================================
using namespace std;

bool show_mess=1; // set to 0 to make things go silent

//==================================================================================================
//
// print message to screen
//
//==================================================================================================
void print_message(string mess)
{
    if(show_mess)
    {
        cout << "\n " << setfill('=') << setw(90) << "=" << endl;  
        cout << " || " << mess << endl;
        cout << " " << setfill('=') << setw(90) << "=" << endl << endl;   
    }
    
    return;
}

//==================================================================================================
//
// header for file
//
//==================================================================================================
template <typename somestream>
void output_header(somestream &ofile, 
                   double xmin, double xmax, int np, 
                   double Dtau, double Te, 
                   double betac, double muc, 
                   double betao, double muo, 
                   double eps_Int, int Te_order, int beta_order,
                   string mode)
{
    ofile.precision(6);

    ofile << setfill('#') << setw(90) << "#" << endl << "#" << endl;  
    ofile << "# SZ signal computed with " << SZpack_version << endl;
    ofile << "# Runmode: " << mode << endl;
    ofile << "#\n# General parameters:" << endl;
    ofile << "# xmin= " << xmin << " xmax= " << xmax << " np= " << np << endl;
    ofile << "# Dtau= " << Dtau << endl;
    ofile << "# Te= " << Te << "keV (The ~ " << Te/const_me << ")" << endl;
    ofile << "# betac= " << betac << " muc= " << muc << endl;
    ofile << "# betao= " << betao << " muo= " << muo << endl;
    ofile << "#\n# Runmode specific parameters:\n# eps= " << eps_Int << endl;
    ofile << "# Te_order= " << Te_order << " beta_order= " << beta_order << endl;
    ofile << "#\n" << setfill('#') << setw(90) << "#" << endl; 
    ofile << "#\n# Output format: x = (h nu/k T0) | x^3 Dn(x) | DI(x) in MJy/sr " << endl;
    ofile << "# Here a CMB temperature of T0 = " << T0_CMB << " K is assumed." << endl;
    ofile << "#\n" << setfill('#') << setw(90) << "#" << endl; 
    
    return;
}

//--------------------------------------------------------------------------------------------------
void output_to_file(ofstream &f, double x, double Dn)
{
    double x3=pow(x, 3);
    f << x << " " << x3*Dn << " " << Dn_DI_conversion*x3*Dn << endl;
    if(show_mess) cout  << " x= " << x << " " << x3*Dn << " " << Dn_DI_conversion*x3*Dn << endl;
    
    return;
}

//==================================================================================================
//
// compute distortion
//
//==================================================================================================
void output_SZ_distortion_5D(string fname, 
                             vector<double> &xa, 
                             double Dtau, double Te, 
                             double betac, double muc, 
                             double betao, double muo, 
                             double eps_Int=1.0e-4)
{
    print_message("Carrying out 5D integral");
    
    int np=xa.size();
    double xmin=xa[0], xmax=xa[np-1];

    ofstream ofile(fname.c_str());
    
    output_header(ofile, xmin, xmax, np, Dtau, Te, betac, muc, betao, muo, 
                  eps_Int, 0, 0, "5D integral");
    
    ofile.precision(16);

    for(int k=0; k<np; k++)
    {
        double dum=compute_SZ_signal_5D(xa[k], Dtau, Te, 
                                        betac, muc, betao, muo, 
                                        eps_Int);

        output_to_file(ofile, xa[k], dum);
    }
        
    ofile.close();
        
    return;
}

//==================================================================================================
void output_SZ_distortion_3D(string fname, 
                             vector<double> &xa, 
                             double Dtau, double Te, 
                             double betac, double muc, 
                             double betao, double muo, 
                             double eps_Int=1.0e-4)
{
    print_message("Carrying out 3D integral");

    int np=xa.size();
    double xmin=xa[0], xmax=xa[np-1];
    
    ofstream ofile(fname.c_str());

    output_header(ofile, xmin, xmax, np, Dtau, Te, betac, muc, betao, muo, 
                  eps_Int, 0, 0, "3D integral");

    ofile.precision(16);

    for(int k=0; k<np; k++)
    {
        double dum=compute_SZ_signal_3D(xa[k], Dtau, Te, 
                                        betac, muc, betao, muo, 
                                        eps_Int);

        output_to_file(ofile, xa[k], dum);
    }
    
    ofile.close();
    
    return;
}

//==================================================================================================
void output_SZ_distortion_asym(string fname, 
                               vector<double> &xa, 
                               double Dtau, double Te, 
                               double betac, double muc, 
                               double betao, double muo,
                               int Te_order, int beta_order)
{
    print_message("Using asymptotic expansion of collision integral (Te<13 keV at most)");
    
    int np=xa.size();
    double xmin=xa[0], xmax=xa[np-1];
    
    ofstream ofile(fname.c_str());
    
    output_header(ofile, xmin, xmax, np, Dtau, Te, betac, muc, betao, muo, 
                  0, Te_order, beta_order, "asymptotic");
    
    ofile.precision(16);
    
    for(int k=0; k<np; k++)
    {
        double dum=compute_SZ_signal_asymptotic(xa[k], Dtau, Te, 
                                                betac, muc, betao, muo, 
                                                Te_order, beta_order);
        
        output_to_file(ofile, xa[k], dum);
    }
    
    ofile.close();
    
    return;
}

//==================================================================================================
void output_SZ_distortion_CNSN(string fname, 
                               vector<double> &xa, 
                               double Dtau, double Te, 
                               double betac, double muc, 
                               double betao, double muo,
                               int Te_order, int beta_order)
{
    print_message("Using improved basis of CNSN 2012 (2keV < Te < 75 keV)");
    
    int np=xa.size();
    double xmin=xa[0], xmax=xa[np-1];
    
    ofstream ofile(fname.c_str());
    
    output_header(ofile, xmin, xmax, np, Dtau, Te, betac, muc, betao, muo, 
                  0, Te_order, beta_order, "CNSN-basis");
    
    ofile.precision(16);
    
    for(int k=0; k<np; k++)
    {
        double dum=compute_SZ_signal_CNSN_basis(xa[k], Dtau, Te, 
                                                betac, muc, betao, muo, 
                                                Te_order, beta_order);

        output_to_file(ofile, xa[k], dum);
    }
    
    ofile.close();
    
    return;
}

//==================================================================================================
void output_SZ_distortion_CNSN_opt(string fname,
                                   vector<double> &xa,
                                   double Dtau, double Te,
                                   double betac, double muc,
                                   double betao, double muo,
                                   int kmax, int beta_order, int accuracy_level)
{
    print_message("Using improved basis of CNSN 2012 with optimization of temperature terms");
    
    int np=xa.size();
    double xmin=xa[0], xmax=xa[np-1];
    
    ofstream ofile(fname.c_str());
    
    output_header(ofile, xmin, xmax, np, Dtau, Te, betac, muc, betao, muo,
                  0, kmax, beta_order, "CNSNopt-basis");
    
    ofile.precision(16);
    
    for(int k=0; k<np; k++)
    {
        double dum=compute_SZ_signal_CNSN_basis_opt(xa[k], Dtau, Te,
                                                    betac, muc, betao, muo,
                                                    kmax, beta_order, accuracy_level);
        
        output_to_file(ofile, xa[k], dum);
    }
    
    ofile.close();
    
    return;
}

//==================================================================================================
void output_SZ_distortion_combo(string fname, 
                                vector<double> &xa, 
                                double Dtau, double Te, 
                                double betac, double muc, 
                                double betao, double muo)
{
    print_message("Using combo of asymptotic expansion and CNSN basis (Te < 75 keV at most)");
    
    int np=xa.size();
    double xmin=xa[0], xmax=xa[np-1];
    
    ofstream ofile(fname.c_str());
    
    output_header(ofile, xmin, xmax, np, Dtau, Te, betac, muc, betao, muo, 
                  0, 0, 0, "combo");
    
    ofile.precision(16);
    
    for(int k=0; k<np; k++)
    {
        double dum=compute_SZ_signal_combo(xa[k], Dtau, Te, betac, muc, betao, muo);
        
        output_to_file(ofile, xa[k], dum);
    }
    
    ofile.close();
    
    return;
}

//==================================================================================================
//
// This function can be used to compute the precision of the different expansions. The results can 
// be compared with the full numerical result or the approximation obtained with the combo function.
// For Te < 75 eV the combo function should be sufficient as reference, but above currently the 
// integrals have to be carried out explicitly.
//
//==================================================================================================
void compute_precision_of_basis(string fname, double xmin, double xmax)
{
    print_message("Computing the precision of the basis");
    
    int np=300, nT=50;
    vector<double> xa(np), ya_3D(np), ya(np);
    vector<double> Tea(nT);

    init_xarr(xmin, xmax, &xa[0], np, 1, 0);  // change 1 --> 0 to have linear grid in x
    init_xarr(1.0, 15.0, &Tea[0], nT, 0, 0);  // change 1 --> 0 to have linear grid in x
    
    ofstream ofile(fname.c_str());
    ofile.precision(16);
    
    double Dtau=0.01, DI_s=4.881e-5; // fiducial accuracy according to CSNN 2012
    double betac=0.01, muc=1.0, betao=0.0, muo=0.0;
    
    for(int t=0; t<nT; t++)
    {
        double Te=Tea[t];

        ofile << Te << " ";        
        if(show_mess) cout << " Te= "<< Te << " keV done " << endl;
        
        //==================================================================
        // compute precise solution
        //==================================================================
        for(int k=0; k<np; k++) 
            
//            ya_3D[k]=pow(xa[k], 3)*compute_SZ_signal_3D(xa[k], Dtau, Te, betac, muc,  
//                                                        betao, muo, 1.0e-5);
        
            ya_3D[k]=pow(xa[k], 3)*compute_SZ_signal_combo(xa[k], Dtau, Te, 
                                                           betac, muc, betao, muo);
        
        //==================================================================
        // compute approximation
        //==================================================================
        for(int m=0; m<=10; m++)
        {
            for(int k=0; k<np; k++) 

                ya[k]=pow(xa[k], 3)*compute_SZ_signal_asymptotic(xa[k], Dtau, Te, betac, muc, 
                                                                 betao, muo, m, 2);

//                ya[k]=pow(xa[k], 3)*compute_SZ_signal_CNSN_basis(xa[k], Dtau, Te, betac, muc, 
//                                                                 betao, muo, m, 2);
            
            //==============================================================
            // find maximal absolute deviation and save to file
            //==============================================================
            double max=0.0;
            
            for(int k=0; k<np; k++)
                if(fabs(ya[k]-ya_3D[k])>=max) max=fabs(ya[k]-ya_3D[k]);
            
            ofile << max/DI_s << " ";            
        }
        
        ofile << endl;
    }
    
    ofile.close();
    
    return;
}

//==================================================================================================
//
// compute distortion using expansion around mean values. 
//
//==================================================================================================
void output_SZ_distortion_means(string fname, 
                                vector<double> &xa, 
                                double Dtau, double Te, 
                                double betac, double muc, double omega[3], double sig[3])
{
    print_message("Using expansion around mean values");
    
    int np=xa.size();
    ofstream ofile(fname.c_str());
    
    output_header(ofile, xa[0], xa[np-1], np, Dtau, Te, 
                  betac, muc, 0, 0, 0, 0, 0, "means");
    
    ofile.precision(16);

    // additional, second order parameters
    double betapar=betac*muc;
    double beta2perp=betac*betac*(1.0-muc*muc);
    double kappa=0.0;
    
    for(int k=0; k<np ; k++)
    {
        double dum=compute_SZ_signal_combo_means_ex(xa[k], Dtau, Te, betapar,
                                                    omega, sig, kappa, beta2perp);
        
        output_to_file(ofile, xa[k], dum);
    }
    
    ofile.close();
    
    return;
}

//==================================================================================================
//
// output for Fig. 8 & 9; this is meant to be an example for the use of this function
//
//==================================================================================================
void output_SZ_null(string fname)
{
    print_message("Computing null of SZ signal");
    
    ofstream ofile(fname.c_str());
    
    ofile.precision(16);
    
    int nT=200;
    vector<double> Tea(nT);
    
    init_xarr(1.0, 70.0, &Tea[0], nT, 1, 0);  // change 1 --> 0 to have linear grid in x
    
    for(int k=0; k<nT; k++)
    {
        double Te=Tea[k];
        double x0=compute_null_of_SZ_signal(1, Te, 0.0, 0, 0, 0, 0);
        
        ofile << Te << " " 
              << x0 << " " 
              << compute_null_of_SZ_signal(1, Te, 0.001, 0, 0, 0, 0) << " " 
              << compute_null_of_SZ_signal(1, Te, -0.001, 0, 0, 0, 0) << " " 
              << compute_null_of_SZ_signal(1, Te, 0.005, 0, 0, 0, 0) << " " 
              << compute_null_of_SZ_signal(1, Te, -0.005, 0, 0, 0, 0) << " " 
              << compute_null_of_SZ_signal(1, Te, 0.01, 0, 0, 0, 0) << " " 
              << compute_null_of_SZ_signal(1, Te, -0.01, 0, 0, 0, 0) << " " 
              //
              << compute_null_of_SZ_signal(1, Te, 0.0, 0.01, 0, 0, 0) - x0 << " " 
              << compute_null_of_SZ_signal(1, Te, 0.0, 0.05, 0, 0, 0) - x0 << " " 
              << compute_null_of_SZ_signal(1, Te, 0.0, 0.1, 0, 0, 0) - x0 << " " 
              << compute_null_of_SZ_signal(1, Te, 0.0, 0.2, 0, 0, 0) - x0 << " " 
              << compute_null_of_SZ_signal(1, Te, 0.0, 0.4, 0, 0, 0) - x0 << " " 
              //
              << compute_null_of_SZ_signal(1, Te, 0.0, 0.0, 0.001, 0, 0) - x0 << " " 
              << compute_null_of_SZ_signal(1, Te, 0.0, 0.0, -0.001, 0, 0) - x0 << " " 
              << compute_null_of_SZ_signal(1, Te, 0.0, 0.0, 0.002, 0, 0) - x0 << " " 
              << compute_null_of_SZ_signal(1, Te, 0.0, 0.0, -0.002, 0, 0) - x0 << " " 
              << compute_null_of_SZ_signal(1, Te, 0.0, 0.0, 0.005, 0, 0) - x0 << " " 
              << compute_null_of_SZ_signal(1, Te, 0.0, 0.0, -0.005, 0, 0) - x0 << " " 
              << compute_null_of_SZ_signal(1, Te, 0.0, 0.0, 0.01, 0, 0) - x0 << " " 
              << compute_null_of_SZ_signal(1, Te, 0.0, 0.0, -0.01, 0, 0) - x0 << " " 
              //
              << compute_null_of_SZ_signal(1, Te, 0.0, 0.0, 0.0, 0.0001, 0) - x0 << " " 
              << compute_null_of_SZ_signal(1, Te, 0.0, 0.0, 0.0, 0.0, 0.0001) - x0 << " ";
        
        double Ix=pow(x0, 3)*compute_SZ_signal_combo(x0, 1.0, Te, 0, 0, 0, 0);
        double xp=x0*(1.0+0.001), xm=x0*(1.0-0.001);
        double Ixp=pow(xp, 3)*compute_SZ_signal_combo(xp, 1.0, Te, 0, 0, 0, 0);
        double Ixm=pow(xm, 3)*compute_SZ_signal_combo(xm, 1.0, Te, 0, 0, 0, 0);
        double dI=(Ixp-Ixm)/(2.0*0.001);
        double d2I=(Ixp-2.0*Ix+Ixm)/pow(0.001, 2) / 2.0;
        
        ofile << -d2I/dI << endl;
        
        cout << Te << endl;
    }
    
    ofile.close();
    
    return;
}

//==================================================================================================
//
// output for Fig. 3; this is meant to be an example for the use of derivative functions
//
//==================================================================================================
void output_SZ_distortion_derivatives(string fname, 
                                      vector<double> &xa, 
                                      double Dtau, double Te, 
                                      double betac, double muc, int kmax)
{
    print_message("Computing derivatives of S in CMB rest frame");
    
    int np=xa.size();
    
    ofstream ofile(fname.c_str());
    
    ofile.precision(16);
    
    vector<double> dDn_dThe;
    
    for(int k=0; k<np; k++)
    {
        double x3=pow(xa[k], 3);
        ofile << xa[k] << " ";
        
        Dcompute_SZ_signal_combo_CMB(xa[k], kmax, 0, 0, Dtau, Te, betac, muc, dDn_dThe);
        for(int m=0; m<=kmax; m++) ofile << dDn_dThe[m]*x3 << " ";
        
        Dcompute_SZ_signal_combo_CMB(xa[k], kmax, 1, 0, Dtau, Te, betac, muc, dDn_dThe);
        for(int m=0; m<=kmax; m++) ofile << dDn_dThe[m]*x3 << " ";
        
        Dcompute_SZ_signal_combo_CMB(xa[k], kmax, 2, 0, Dtau, Te, betac, muc, dDn_dThe);
        for(int m=0; m<=kmax; m++) ofile << dDn_dThe[m]*x3 << " ";
        
        Dcompute_SZ_signal_combo_CMB(xa[k], kmax, 0, 1, Dtau, Te, betac, muc, dDn_dThe);
        for(int m=0; m<=kmax; m++) ofile << dDn_dThe[m]*x3 << " ";
        
        //------------------------------------------------------------
        // x-derivatives
        //------------------------------------------------------------
        double Ix=pow(xa[k], 3)*compute_SZ_signal_combo(xa[k], 1.0, Te, 0, 0, 0, 0);
        double xp=xa[k]*(1.0+0.001), xm=xa[k]*(1.0-0.001);
        double Ixp=pow(xp, 3)*compute_SZ_signal_combo(xp, 1.0, Te, 0, 0, 0, 0);
        double Ixm=pow(xm, 3)*compute_SZ_signal_combo(xm, 1.0, Te, 0, 0, 0, 0);
        double dI=(Ixp-Ixm)/(2.0*0.001);
        double d2I=(Ixp-2.0*Ix+Ixm)/pow(0.001, 2) / 2.0;
        
        ofile << dI << " " << d2I << endl;
        
        cout  << " x= " << xa[k] << " " << dDn_dThe[0]*x3 << endl;
    }
    
    ofile.close();
    
    return;
}


//==================================================================================================
//
// example for a two temperature case
//
//==================================================================================================
void output_SZ_distortion_two_temp(string fname, 
                                   vector<double> &xa, 
                                   double Dtau, double ftau, 
                                   double Te, double DT_T,
                                   double betac, double muc)
{
    print_message("Two-temperature case");
    
    int np=xa.size();
    
    ofstream ofile(fname.c_str());

    ofile.precision(16);
    
    double DI_s=4.881e-5;   // fiducial precision
    int betao=0.0, muo=0.0; // resting observer
    
    for(int k=0; k<np; k++)
    {
        double x3=pow(xa[k], 3);
        double Tb=Te*(1.0+DT_T);
        double S1=x3*compute_SZ_signal_combo(xa[k], Dtau*(1.0-ftau), Te, betac, muc, betao, muo);
        double S2=x3*compute_SZ_signal_combo(xa[k], Dtau*ftau, Tb, betac, muc, betao, muo);
        
        ofile << xa[k] << " ";
        
        double Tm=Te*(1.0+ftau*DT_T);
        double omega1=ftau*(1.0-ftau) * pow(DT_T/(1.0+ftau*DT_T), 2);
        double omega2=ftau*(1.0-ftau)*(1.0-2.0*ftau) * pow(DT_T/(1.0+ftau*DT_T), 3);
        
        vector<double> dDn_dThe(4);
        Dcompute_SZ_signal_combo_CMB(xa[k], 3, 0, 0, Dtau, Tm, betac, muc, dDn_dThe);
        
        double Sm=x3*dDn_dThe[0];
        double d2Sm=x3*dDn_dThe[2];
        double d3Sm=x3*dDn_dThe[3];
        
        if(k==0) cout << omega1 << " " << omega2 << " " << (1.0-2.0*ftau)/ftau << endl;
        
        ofile << S1+S2 << " " << Sm << " " 
              << (S1+S2-Sm)/DI_s << " " << (S1+S2-Sm)/d2Sm << " " 
              << (S1+S2-(Sm+d2Sm*omega1))/DI_s << " " 
              << d2Sm*omega1/DI_s << " " << (d2Sm*omega1 + d3Sm*omega2)/DI_s << " ";
        
        ofile << endl;
    }
    
    ofile.close();
    
    return;
}

//==================================================================================================
//
// compute importance of relativistic corrections
//
//==================================================================================================
void output_SZ_distortion_rel_corrs(string fname,
                                    vector<double> &xa,
                                    double Dtau, double Te,
                                    double betac, double muc)
{
    print_message("Using expansion around mean values and computing relativistic correction");
    
    int np=xa.size();
    ofstream ofile(fname.c_str());
    
    output_header(ofile, xa[0], xa[np-1], np, Dtau, Te,
                  betac, muc, 0, 0, 0, 0, 0, "rel_corrs");
    
    ofile.precision(16);
    
    // additional, second order parameters
    double betapar=betac*muc;
    double beta2perp=betac*betac*(1.0-muc*muc);
    
    for(int k=0; k<np ; k++)
    {
        double dum=compute_SZ_signal_combo_means(xa[k], Dtau, Te, betapar, 0, 0, 0, beta2perp);
        
        dum-=compute_SZ_signal_asymptotic(xa[k], Dtau, Te, 0, 0, 0, 0, 0, 0);
        
        // in kJy / sr
        output_to_file(ofile, xa[k], 1.0e+3*dum);
    }
    
    ofile.close();
    
    return;
}

//==================================================================================================
//
// compute importance of temperature dispersion
//
//==================================================================================================
void output_SZ_distortion_temp_dis(string fname,
                                    vector<double> &xa,
                                    double Dtau, double Te,
                                    double betac, double muc,
                                    double omega, double sig)
{
    print_message("Using expansion around mean values and computing temperature dispersion");
    
    int np=xa.size();
    ofstream ofile(fname.c_str());
    
    output_header(ofile, xa[0], xa[np-1], np, Dtau, Te,
                  betac, muc, 0, 0, 0, 0, 0, "temp_dis");
    
    ofile.precision(16);
    
    // additional, second order parameters
    double betapar=betac*muc;
    double beta2perp=betac*betac*(1.0-muc*muc);
    
    for(int k=0; k<np ; k++)
    {
        double dum=compute_SZ_signal_combo_means(xa[k], Dtau, Te, betapar,
                                                 omega, sig, 0, beta2perp);
        
        dum-=compute_SZ_signal_combo_means(xa[k], Dtau, Te, betapar, 0, 0, 0, beta2perp);
        
        // in kJy / sr
        output_to_file(ofile, xa[k], 1.0e+3*dum);
    }
    
    ofile.close();
    
    return;
}

//==================================================================================================
//
// read parameters
//
//==================================================================================================
void read_parameters(string fname,
                     double &xmin, double &xmax, int &np, 
                     double &Dtau, double &Te, 
                     double &betac, double &muc, 
                     double &betao, double &muo,
                     int &Te_order, int &beta_order, double &eps_Int,
                     int &kmax, int &accuracy_level,
                     string &outpath, string &add)
{
    ifstream ifile(fname.c_str());
    
    ifile >> xmin; ifile >> xmax; ifile >> np;
    ifile >> Dtau; ifile >> Te;
    ifile >> betac; ifile >> muc;
    ifile >> betao; ifile >> muo;
    ifile >> Te_order; ifile >> beta_order; ifile >> eps_Int;
    ifile >> kmax; ifile >> accuracy_level;
    ifile >> outpath; ifile >> add;
    
    //==================================================================
    // simple sanity checks
    //==================================================================
    if(xmin<0.01 || xmax>50.0){ cerr << " read_parameters :: change x-range " << endl; exit(0); }
    if(np<=0){ cerr << " read_parameters :: change number of frequency points " << endl; exit(0); }

    if(Dtau<0){ cerr << " read_parameters :: change Dtau " << endl; exit(0); }

    if(betac<0 || betac>0.1){ cerr << " read_parameters :: check betac " << endl; exit(0); }
    if(muc<-1.0 || muc>1){ cerr << " read_parameters :: check muc " << endl; exit(0); }
    if(betao<0 || betao>0.1){ cerr << " read_parameters :: check betao " << endl; exit(0); }
    if(muo<-1.0 || muo>1){ cerr << " read_parameters :: check muo " << endl; exit(0); }

    if(eps_Int<1e-8 || eps_Int>0.001){cerr << " read_parameters :: check eps_Int "<< endl; exit(0);}

    if(accuracy_level<0 || accuracy_level>3)
    {cerr << " read_parameters :: check accuracy_level "<< endl; exit(0);}
    
    return;
}
 
//==================================================================================================
//
// main call of routine. Here parameters are handled.
//
//==================================================================================================
int main(int narg, char *args[])
{
    //==================================================================
    // example parameters
    //==================================================================
    double xmin=0.1;
    double xmax=30.0;
    int np=100;
    //
    double Dtau=0.01;
    double Te=0.03*const_me;
    //
    double betac=0.01;
    double muc=1.0;
    //
    double betao=0.001241;
    double muo=0.0;
    //
    double eps_Int=1.0e-4;
    int Te_order=10;
    int beta_order=2;
    //
    int kmax=Te_order=5;
    int accuracy_level=1;
    //
    string outpath=SZPACKDIR+"./outputs/";
    string add=".dat";
    string mode="3D";
    
    //==================================================================
    // startup with different number of arguments
    //==================================================================
    if(narg==1){}
    
    else if(narg==2) mode=args[narg-1];
    
    else if(narg==3)
    {
        mode=args[narg-2];
        string fname=args[narg-1];
        
        print_message("Reading parameters from parameter file: " + fname);
        
        read_parameters(fname, xmin, xmax, np, Dtau, Te, 
                        betac, muc, betao, muo,
                        Te_order, beta_order, eps_Int,
                        kmax, accuracy_level,
                        outpath, add);
    }

    else{ cerr << " Too many/few parameters " << endl; exit(1); }
        
    //==================================================================
    // uncomment this to use the convention of Nozawa et al. 2006 for 
    // the optical depth. Interpretation of the SZ signal is different
    // in this case. The convention of Chluba et al., 2012 is used  
    // otherwise, which allow a clean separation of kinematic and 
    // scattering effect.
    //==================================================================
    //use_Nozawa2006_convention();
    //add=".Nozawa_conv"+add;

    //==================================================================
    // setting up frequency points. Here in principle any frequency 
    // distribution can be used.
    //==================================================================
    vector<double> xcmb(np); 
    init_xarr(xmin, xmax, &xcmb[0], np, 1, 0); // change 1 --> 0 to have 
                                               // linear grid in x
    
    //==================================================================
    // call in the different runmodes
    //==================================================================
    if(mode=="5D")
    {
        output_SZ_distortion_5D(outpath+"SZ_Integral.5D"+add, 
                                xcmb, Dtau, Te, 
                                betac, muc, betao, muo, 
                                eps_Int);
    }
    
    else if(mode=="3D")
    {
        output_SZ_distortion_3D(outpath+"SZ_Integral.3D"+add, 
                                xcmb, Dtau, Te, 
                                betac, muc, betao, muo, 
                                eps_Int);
    }
    
    else if(mode=="ASYM")
    {
        output_SZ_distortion_asym(outpath+"SZ_asymptotic"+add, 
                                  xcmb, Dtau, Te, 
                                  betac, muc, betao, muo, 
                                  Te_order, beta_order);
    }
    
    //==================================================================
    // Added the basis at additional temperature pivots. This function 
    // should give very precise results for 2keV < Te < 75keV
    //==================================================================
    else if(mode=="CNSN")
    {
        output_SZ_distortion_CNSN(outpath+"SZ_CNSN_basis"+add,
                                  xcmb, Dtau, Te, 
                                  betac, muc, betao, muo, 
                                  Te_order, beta_order);
    }
    
    //==================================================================
    // Added the basis at additional temperature pivots. This function
    // should give very precise results for 2keV < Te < 75keV
    //==================================================================
    else if(mode=="CNSNopt")
    {
        output_SZ_distortion_CNSN_opt(outpath+"SZ_CNSN_opt_basis"+add,
                                      xcmb, Dtau, Te,
                                      betac, muc, betao, muo,
                                      kmax, beta_order, accuracy_level);
    }
    //==================================================================
    // combination of asymptotic expansion + basis functions of CNSN2012
    // This function should give very precise results for Te < 75keV and
    // can be used for comparison.
    //==================================================================
    else if(mode=="COMBO")
    {
        output_SZ_distortion_combo(outpath+"SZ_combo"+add, 
                                   xcmb, Dtau, Te, 
                                   betac, muc, betao, muo);
    }
        
    //==================================================================
    // check the precision of the basis function. Settings should be 
    // changed directly in compute_precision_of_basis() above
    //==================================================================
    else if(mode=="ACC")
    {
        compute_precision_of_basis(outpath+"rec.moments.asym.beta"+add, 
                                   xmin, xmax);
    }

    //==================================================================
    // expansion of SZ signal around mean values
    //==================================================================
    else if(mode=="MEANS")
    {
        //=======================================================================
        // examples from CSNN 2012
        //=======================================================================
        //double Dtau=7.7563e-3, Te=4.1051; double omegas[3]={0.066421, 0.021254, 0};
        //double Dtau=3.1682e-3, Te=6.3163; double omegas[3]={0.14538, -0.021378, 0};
        //double Dtau=0.62251e-3, Te=1.1251; double omegas[3]={0.69288, 1.1797, 0};
        //double Dtau=1.3738e-3, Te=3.8462; double omegas[3]={0.29556, 0.17844, 0};
        //double sigmas[3]={0, 0, 0};
        
        //=======================================================================
        double omegas[3]={0.2, 0.1, 0.05};     // omega^(1), omega^(2), omega^(3)
        double sigmas[3]={0, 0, 0};            // sigma^(1), sigma^(2), sigma^(3)

        // for additional second order parameters see function itself
        output_SZ_distortion_means(outpath+"SZ_means"+add,
                                   xcmb, Dtau, Te, betac, muc, omegas, sigmas);
    }

    //==================================================================
    // compute null of SZ signal. For details see function
    //==================================================================
    else if(mode=="NULL")
    {
        output_SZ_null(outpath+"SZ_null"+add);
    }

    //==================================================================
    // compute derivatives of SZ signal 
    //==================================================================
    else if(mode=="DERIVS")
    {
        Dtau=1.0;
        betac=0.0;
        
        output_SZ_distortion_derivatives(outpath+"SZ_moments.derivs.10keV"+add, xcmb, 
                                         Dtau, 10, betac, muc, 4);
        output_SZ_distortion_derivatives(outpath+"SZ_moments.derivs.20keV"+add, xcmb, 
                                         Dtau, 20, betac, muc, 4);
        output_SZ_distortion_derivatives(outpath+"SZ_moments.derivs.30keV"+add, xcmb, 
                                         Dtau, 30, betac, muc, 4);
        output_SZ_distortion_derivatives(outpath+"SZ_moments.derivs.40keV"+add, xcmb, 
                                         Dtau, 40, betac, muc, 4);
        output_SZ_distortion_derivatives(outpath+"SZ_moments.derivs.50keV"+add, xcmb, 
                                         Dtau, 50, betac, muc, 4);
    }
    
    //==================================================================
    // SZ signal for two-temperature case
    //==================================================================
    else if(mode=="TWOT")
    {
        double ftau=0.2;
        double Te1=5;
        double Delta=0.5;
        
        output_SZ_distortion_two_temp(outpath+"SZ_moments.two-temp.ftau_0.2.Te_5.Delta_0.5"+add, 
                                      xcmb, Dtau, ftau, Te1, Delta, betac, muc);
    }
    
    //==================================================================
    // expansion of SZ signal around mean values
    //==================================================================
    else if(mode=="RELCORR")
    {
        //=======================================================================
        double omegas[3]={0.0, 0.0, 0.0};      // omega^(1), omega^(2), omega^(3)
        double sigmas[3]={0, 0, 0};            // sigma^(1), sigma^(2), sigma^(3)

        double PInt=0.57;
        
        Te=3.0;
        string adds=".3keV";
        betac=0.0;
        Dtau=0.00636608*pow(Te/5, 0.49)*PInt;
        
        // for additional second order parameters see function itself
        output_SZ_distortion_means(outpath+"SZ_means"+adds+add,
                                   xcmb, Dtau, Te, betac, muc, omegas, sigmas);
        
        output_SZ_distortion_rel_corrs(outpath+"SZ_rel_corrs"+adds+add,
                                       xcmb, Dtau, Te, betac, muc);

        output_SZ_distortion_temp_dis(outpath+"SZ_temp_dis"+adds+add,
                                      xcmb, Dtau, Te, betac, muc, 0.2, 0);
    }

    //==================================================================
    // expansion of SZ signal around mean values
    //==================================================================
    else if(mode=="RELCORRSUPERCORE")
    {
        //=======================================================================
        double omegas[3]={0.0, 0.0, 0.0};      // omega^(1), omega^(2), omega^(3)
        double sigmas[3]={0, 0, 0};            // sigma^(1), sigma^(2), sigma^(3)
        
        outpath="./outputs.SuperCOrE/";
        Te=3.0;
        string adds=".3keV";
        betac=0.0;
        Dtau=const_me/Te*1.0e-4;
        
        // for additional second order parameters see function itself
        output_SZ_distortion_means(outpath+"SZ_means"+adds+add,
                                   xcmb, Dtau, Te, betac, muc, omegas, sigmas);
        
        output_SZ_distortion_rel_corrs(outpath+"SZ_rel_corrs"+adds+add,
                                       xcmb, Dtau, Te, betac, muc);
        
        output_SZ_distortion_temp_dis(outpath+"SZ_temp_dis"+adds+add,
                                      xcmb, Dtau, Te, betac, muc, 0.2, 0);
    }
    
    //==================================================================
    else cerr << " Unknown runmode: " << mode << endl;
    
    return 0;
}

//==================================================================================================
//==================================================================================================
