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
#include "SZpack.h"

#include "SZ_moment_method.h"
#include "SZ_cluster_profiles.h"

#include "gsl/gsl_multimin.h"
#include <gsl/gsl_linalg.h>

//==================================================================================================
//
// namespaces
//
//==================================================================================================
using namespace std;

bool show_mess=1;

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
                   int Te_order, int beta_order,
                   double accuracy_goal, double Te_max,
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
    ofile << "#\n# Runmode specific parameters:" << endl;
    ofile << "# Te_order= " << Te_order << " beta_order= " << beta_order << endl;
    ofile << "# accuracy_goal= " << accuracy_goal << " Te_max= " << Te_max << endl;
    ofile << "#\n" << setfill('#') << setw(90) << "#" << endl; 
    ofile << "#\n# Output format: x = (h nu/k T0) | x^3 Dn(x) | DI(x) in MJy/sr " << endl;
    ofile << "# Here a CMB temperature of T0 = " << T0_CMB << " K is assumed." << endl;
    ofile << "#\n" << setfill('#') << setw(90) << "#" << endl; 
    
    return;
}

//==================================================================================================
//
// routine to find best-fit solutions using isothermal model
//
//==================================================================================================
struct minimizer_Data
{
    vector<double> xa, Ia;
    vector<double> Iapprox;
};

//--------------------------------------------------------------------------------------------------
double f_chisq_betac(const gsl_vector *v, void *params)
{
    double tau=gsl_vector_get(v, 0);
    double Teiso=gsl_vector_get(v, 1);
    double betac=gsl_vector_get(v, 2);
    double muc=1.0, betao=0.0, muo=0.0;
    
    minimizer_Data *p = (minimizer_Data *)params;
    
    double residual=0.0;
    for(int k=0; k<(int)p->xa.size(); k++)
    {
        p->Iapprox[k]=pow(p->xa[k], 3)*compute_SZ_signal_combo(p->xa[k], tau, Teiso, 
                                                               betac, muc, betao, muo);
        
        residual+=pow(p->Iapprox[k]-p->Ia[k], 2);
    }
    //    cout << residual << " " << tau << " " << Teiso << endl;
    
    return residual;
}

double f_chisq(const gsl_vector *v, void *params)
{
    double tau=gsl_vector_get(v, 0);
    double Teiso=gsl_vector_get(v, 1);
    double betac=0.0, muc=0.0, betao=0.0, muo=0.0;
    
    minimizer_Data *p = (minimizer_Data *)params;
    
    double residual=0.0;
    for(int k=0; k<(int)p->xa.size(); k++)
    {
        p->Iapprox[k]=pow(p->xa[k], 3)*compute_SZ_signal_combo(p->xa[k], tau, Teiso, 
                                                               betac, muc, betao, muo);
        
        residual+=pow(p->Iapprox[k]-p->Ia[k], 2);
    }
    //    cout << residual << " " << tau << " " << Teiso << endl;
    
    return residual;
}

//--------------------------------------------------------------------------------------------------
int find_best_fit_isothermal_model(vector<double> &xa, vector<double> &Ia,
                                   vector<double> &solution, vector<double> &Isol, 
                                   int numpar=2)
{
    minimizer_Data MinData;
    
    MinData.xa=xa;
    MinData.Ia=Ia;
    MinData.Iapprox.resize(Ia.size());
    
    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
    gsl_multimin_fminimizer *s = NULL;
    gsl_vector *ss, *x;
    gsl_multimin_function minex_func;
    
    size_t iter = 0;
    int status;
    double size;
    
    /* Starting point */
    x = gsl_vector_alloc (numpar);
    gsl_vector_set (x, 0, 0.01);
    gsl_vector_set (x, 1, 5.0);
    if(numpar>2) gsl_vector_set (x, 2, 0.01);
    
    /* Set initial step sizes to 1 */
    ss = gsl_vector_alloc (numpar);
    gsl_vector_set (ss, 0, 0.001);
    gsl_vector_set (ss, 1, 0.5);
    if(numpar>2) gsl_vector_set (ss, 2, 0.001);
    
    /* Initialize method and iterate */
    minex_func.n = numpar;
    
    if(numpar==2) minex_func.f = f_chisq;
    else if(numpar==3) minex_func.f = f_chisq_betac;
    else{ cerr << " find_best_fit_isothermal_model :: N/A " << endl; exit(0); }
    
    minex_func.params = (void *) &MinData;
    
    s = gsl_multimin_fminimizer_alloc (T, minex_func.n);
    gsl_multimin_fminimizer_set (s, &minex_func, x, ss);
    
    do
    {
        iter++;
        status = gsl_multimin_fminimizer_iterate(s);
        
        if(status) break;
        
        size = gsl_multimin_fminimizer_size (s);
        status = gsl_multimin_test_size (size, 1.0e-3);
        
        if (status == GSL_SUCCESS)
        {
            printf (" converged to minimum\n");
        }
    }
    while (status == GSL_CONTINUE && iter < 500);
    
    solution.clear();
    solution.push_back(gsl_vector_get (s->x, 0));
    solution.push_back(gsl_vector_get (s->x, 1));
    if(numpar>2) solution.push_back(gsl_vector_get (s->x, 2));
    Isol=MinData.Iapprox;
    
    for(int k=0; k<(int)minex_func.n; k++) cout << " " << k << " " << solution[k] << endl;
    
    gsl_vector_free(x);
    gsl_vector_free(ss);
    gsl_multimin_fminimizer_free (s);
    
    return status;
}

//==================================================================================================
//
// compute distortion using combo method. This is just to check the precision of the result obtain
// with the temperature-velocity moment method.
//
//==================================================================================================
void output_SZ_distortion_combo(string fname, 
                                vector<double> &xa, 
                                double Dtau, double Te, 
                                double betac, double muc, 
                                double betao, double muo)
{
    print_message("Using combo of asymptotic expansion and CNSN basis (Te < 75 keV at most)");
    
    int np=xa.size();
    ofstream ofile(fname.c_str());
    
    output_header(ofile, xa[0], xa[np-1], np, Dtau, Te, 
                  betac, muc, betao, muo, 
                  0, 0, 0, 0, "combo");
    
    ofile.precision(16);
    
    for(int k=0; k<np ; k++)
    {
        double dum=compute_SZ_signal_combo(xa[k], Dtau, Te, betac, muc, betao, muo);
        
        ofile << xa[k] << " " << dum*pow(xa[k], 3) << endl;
        cout  << " x= " << xa[k] << " " << dum*pow(xa[k], 3) << endl;
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
                                double betac, double muc, 
                                double betao, double muo)
{
    print_message("Using expansion around mean values");
    
    int np=xa.size();
    ofstream ofile(fname.c_str());
    
    output_header(ofile, xa[0], xa[np-1], np, Dtau, Te, 
                  betac, muc, betao, muo, 
                  0, 0, 0, 0, "means");
    
    ofile.precision(16);
    
//    double pars[2]={7.7563e-3, 4.1051}; double omega[3]={0.066421, 0.021254, 0};
//    double pars[2]={3.1682e-3, 6.3163}; double omega[3]={0.14538, -0.021378, 0};
//    double pars[2]={0.62251e-3, 1.1251}; double omega[3]={0.69288, 1.1797, 0};
    double pars[2]={1.3738e-3, 3.8462}; double omega[3]={0.29556, 0.17844, 0};
    
    double sig[3]={0, 0, 0};
    
    for(int k=0; k<np ; k++)
    {
        double dum=compute_SZ_signal_combo_means_ex(xa[k], pars[0], pars[1], 0, omega, sig, 0, 0);
        
        ofile << xa[k] << " " << dum*pow(xa[k], 3) << endl;
        cout  << " x= " << xa[k] << " " << dum*pow(xa[k], 3) << endl;
    }
    
    ofile.close();
    
    return;
}

//==================================================================================================
//
// compute distortions using temperature-velocity moments
//
//==================================================================================================
void compute_isothermal_moments(double Dtau, double Te, 
                                double betac, double muc, 
                                vector<double> &mv, SZ_moment_method &SZM,
                                int variable)
{
    //===============================================================================
    // setting up the moments in different forms
    //===============================================================================
    mv.resize(SZM.Get_total_number_of_moments());
    for(int k=0; k<SZM.Get_total_number_of_moments(); k++) mv[k]=0.0;
    //
    int kmax  =SZM.Get_kmax();
    int beta_order=SZM.Get_order_b();
    
    double The=Te/const_me;
    int dn=2+kmax, dm=1+kmax;
    int ir=SZM.Get_index_of_Te_region(Te);
    
    //===============================================================================
    // low temperature moments (isothermal)
    //===============================================================================
    if(ir==0)
    {
        for(int m=0; m<=kmax; m++) mv[m]=Dtau*pow(The, m+1);
        
        if(beta_order>0)
        {
            for(int m=0; m<=kmax; m++) mv[dn+m]=Dtau*betac*muc*pow(The, m+1);
            
            if(beta_order>1)
            {          
                double b2=betac*betac, P2=1.5*muc*muc-0.5;
                for(int m=0; m<=kmax; m++) mv[2*dn+m]=Dtau*b2   *pow(The, m+1);
                for(int m=0; m<=kmax; m++) mv[3*dn+m]=Dtau*b2*P2*pow(The, m+1);
            }
        }
    }
    
    //===============================================================================
    // pure velocity moments
    //===============================================================================
    if(beta_order>0)
    {
        mv[dn-1]=Dtau*betac*muc;
        
        if(beta_order>1)
        {          
            mv[2*dn-1]=Dtau*betac*betac;
            mv[3*dn-1]=Dtau*betac*betac*(1.5*muc*muc-0.5);
        }
    }
    
    //===============================================================================
    // high temperature moments
    //===============================================================================
    if(ir>0)
    {
        // lookup where the moments have to be written in the moment vector
        int is=SZM.Get_start_index_of_moment_set(ir);
        double The_r=SZM.Get_The_ref(ir), Nf=N_func_CNSN(The);
        
        // (var 0)
        if(variable==0)
        {
            for(int m=0; m<=kmax; m++) mv[is+m]=Dtau*Nf*pow(The-The_r, m);
            
            if(beta_order>0)
            {
                for(int m=0; m<=kmax; m++) mv[is+dm+m]=Dtau*betac*muc*Nf*pow(The-The_r, m);
                
                if(beta_order>1)
                {          
                    double b2=betac*betac, P2=1.5*muc*muc-0.5;
                    for(int m=0; m<=kmax; m++) mv[is+2*dm+m]=Dtau*b2   *Nf*pow(The-The_r, m);
                    for(int m=0; m<=kmax; m++) mv[is+3*dm+m]=Dtau*b2*P2*Nf*pow(The-The_r, m);
                }
            }
        }
    
        // (var 1) 
        if(variable==1)
        {
            for(int m=0; m<=kmax; m++) mv[is+m]=Dtau*Nf*pow(The, m);
            
            if(beta_order>0)
            {
                for(int m=0; m<=kmax; m++) mv[is+dm+m]=Dtau*betac*muc*Nf*pow(The, m);
                
                if(beta_order>1)
                {          
                    double b2=betac*betac, P2=1.5*muc*muc-0.5;
                    for(int m=0; m<=kmax; m++) mv[is+2*dm+m]=Dtau*b2   *Nf*pow(The, m);
                    for(int m=0; m<=kmax; m++) mv[is+3*dm+m]=Dtau*b2*P2*Nf*pow(The, m);
                }
            }
        }
    }
    
    return;
}    

//==================================================================================================
//
// This function gives an example of how the moment method and the class `SZ_moment_method' are used
// to compute the SZ signal. Two variables type for the high temperature moments are available, 
// which should all give practically the same SZ signal. The different steps in the calculation are
// commented below and one can play with things a bit.
//
//==================================================================================================
void output_SZ_distortion_moments_isothermal(string fname, string outpath, string add,
                                             vector<double> &xcmb, 
                                             double Dtau, double Te, 
                                             double betac, double muc, 
                                             double betao, double muo,
                                             SZ_moment_method &SZM,
                                             int Te_order, int beta_order,
                                             int accuracy_goal, double Te_max,
                                             string mess="")
{
    print_message("Using moment method for isothermal cluster to obtain SZ signal"+mess);
    
    //===============================================================================
    //
    // setting up the moments in different form and computing SZ signal
    //
    //===============================================================================
    vector<double> mv, xv;
    vector<vector<double> > DIv(2);

    //===============================================================================
    // Calling the functions of 'SZM' for the two variable. The results should 
    // be identical to machine precision. For details about the variables see the 
    // moment setup routine above and CSNN 2012.
    //===============================================================================
    for(int variable=0; variable<2; variable++)
    {
        compute_isothermal_moments(Dtau, Te, betac, muc, mv, SZM, variable);
        SZM.compute_SZ_signal_moments(mv, xv, DIv[variable], variable);

        //===========================================================================
        // output moment vector
        //===========================================================================
        SZM.show_moment_vector(mv);
        
        string name=outpath+"moments.var_"+int_to_string(variable)+add;
        SZM.export_moment_vector(name, mv);
    }
        
    //===============================================================================
    //
    // output SZ signal
    //
    //===============================================================================
    ofstream ofile(fname.c_str());
    
    output_header(ofile, xcmb[0], xcmb[xcmb.size()-1], xcmb.size(), 
                  Dtau, Te, betac, muc, betao, muo, 
                  Te_order, beta_order, accuracy_goal, Te_max, 
                  "moments isothermal"+mess);
    
    ofile.precision(16);
    
    for(int n=0; n<(int)xv.size(); n++)
        ofile << xv[n] << " " << DIv[0][n] << " " << DIv[1][n] << endl;
    
    ofile.close();

    mv.clear();
    xv.clear();
    DIv.clear();
    
    return;
}


//==================================================================================================
//
// Computation of SZ signal using the simple cluster profile functions of Vikhlinin et al. 2006.
// Here the moments are computed internally by the SZ_moment_method class functions. The profiles
// for Ne and Te are supplied at different slices through the cluster medium.
//
//==================================================================================================
struct profile_slice_params
{
    double xs, ys;
    double zsc;
    double betac, muc;
    double xcmb;
    int k;
    SZ_cluster_profiles *Cluster_ptr;
};
profile_slice_params PSP;

//--------------------------------------------------------------------------------------------------
double zero_func(double l){ return 0.0; }
double betac_func_3D(double lx, double ly, double lz){ return PSP.betac; }
double muc_func_3D(double lx, double ly, double lz){ return PSP.muc; }

double Ne_func_3D(double lx, double ly, double lz)
{ return PSP.Cluster_ptr->Ne(lx*PSP.zsc, ly*PSP.zsc, lz*PSP.zsc); }

double Te_func_3D(double lx, double ly, double lz)
{ return PSP.Cluster_ptr->Te(lx*PSP.zsc, ly*PSP.zsc, lz*PSP.zsc); }

//--------------------------------------------------------------------------------------------------
double betac_func(double l){ return PSP.betac; }
double muc_func(double l){ return PSP.muc; }

double Ne_func(double l)
{ return Ne_func_3D(PSP.xs/PSP.zsc, PSP.ys/PSP.zsc, l); }

double Te_func(double l)
{ return Te_func_3D(PSP.xs/PSP.zsc, PSP.ys/PSP.zsc, l); }

//--------------------------------------------------------------------------------------------------
void output_SZ_distortion_moments_cluster_VIK(string fname, string outpath, string add,
                                              vector<double> &xcmb, 
                                              SZ_cluster_profiles &CL, 
                                              double x_rc, double y_rc,
                                              double betac, double muc, 
                                              double betao, double muo,
                                              SZ_moment_method &SZM,
                                              int Te_order, int beta_order,
                                              int accuracy_goal, double Te_max,
                                              string mess="")
{
    print_message("Using moment method for cluster profiles (VIK) to obtain SZ signal"+mess);

    //===============================================================================
    PSP.betac=betac;
    PSP.muc=muc;

    PSP.Cluster_ptr=&CL;

    PSP.xs=x_rc;
    PSP.ys=y_rc;
    PSP.zsc=1000.0;
    
    //SZM.show_profile_slice(Ne_func, PSP.zsc, 200);
    //SZM.export_profile_slice(outpath+"Ne_slice"+add, Ne_func, PSP.zsc, 200);
    //SZM.show_profile_slice(Te_func, PSP.zsc, 200);
    //SZM.export_profile_slice(outpath+"Te_slice"+add, Te_func, PSP.zsc, 200);

    //===============================================================================
    //
    // setting up the moments in different form and computing SZ signal
    //
    //===============================================================================
    vector<double> mv, xv, dumv, DIsmooth;
    vector<vector<double> > DIv(2), DIviso(2), DInonrel(2);
    vector<vector<double> > TSZetc(2);

    //===============================================================================
    // Calling the functions of 'SZM' for the two variable. The results should 
    // be identical to machine precision. For details about the variables see the 
    // moment setup routine above and CSNN 2012.
    //===============================================================================
    for(int variable=0; variable<2; variable++)
    {
        SZM.compute_SZ_signal_moments(Ne_func, Te_func, betac_func, muc_func,  
                                      PSP.zsc*PSP.Cluster_ptr->Get_rc_cm(), mv, xv, 
                                      DIv[variable], TSZetc[variable], 
                                      variable, 2.0e-6);
        
        cout << " x_rc= " << x_rc << " y_rc= " << y_rc 
        << " tau= " << TSZetc[variable][0] 
        << " ySZ= " << TSZetc[variable][1] 
        << " TSZ= " << TSZetc[variable][2] << " keV" << endl;

        SZM.compute_SZ_signal_moments(Ne_func_3D, Te_func_3D, 
                                      betac_func_3D, muc_func_3D, 
                                      x_rc/PSP.zsc, y_rc/PSP.zsc, 
                                      PSP.zsc*PSP.Cluster_ptr->Get_rc_cm(), 
                                      1.0e-1/PSP.zsc, 
                                      mv, xv, DIv[variable], TSZetc[variable], 
                                      variable, 2.0e-6);
        
        wait_f_r(TSZetc[variable][0]);
        wait_f_r(TSZetc[variable][1]);
        wait_f_r(TSZetc[variable][2]);
               
        //===========================================================================
        // output moment vector
        //===========================================================================
        SZM.show_moment_vector(mv);
        
        cout << " x_rc= " << x_rc << " y_rc= " << y_rc 
             << " tau= " << TSZetc[variable][0] 
             << " ySZ= " << TSZetc[variable][1] 
             << " TSZ= " << TSZetc[variable][2] << " keV" << endl;
        
        string name=outpath+"moments.var_"+int_to_string(variable)+add;
        SZM.export_moment_vector(name, mv);
        
        //===========================================================================
        // compute SZ signal using isothermal moments
        //===========================================================================
        print_message("Now using isothermal approximation");

        compute_isothermal_moments(TSZetc[variable][0], TSZetc[variable][2], 
                                   betac, muc, mv, SZM, variable);
        SZM.compute_SZ_signal_moments(mv, xv, DIviso[variable], variable);
        
        //===========================================================================
        // output moment vector
        //===========================================================================
        SZM.show_moment_vector(mv);
        
        name=outpath+"moments.iso.var_"+int_to_string(variable)+add;
        SZM.export_moment_vector(name, mv);    

        if(variable==0)
        {
            //=======================================================================
            // compute non-relativistic case
            //=======================================================================
            mv[0]=TSZetc[variable][1];
            for(int m=1; m<(int) mv.size(); m++) mv[m]=0.0;

            SZM.compute_SZ_signal_moments(mv, xv, DInonrel[0], 0);
            
            //=======================================================================
            // compute using smooth profile approximation
            //=======================================================================  
            vector<double> omegak;
            
            SZM.compute_SZ_signal_moments_smooth(Ne_func, Te_func, betac_func, muc_func,  
                                                 PSP.zsc*PSP.Cluster_ptr->Get_rc_cm(), 
                                                 omegak, xv, DIsmooth, dumv, 1);    
            
            name=outpath+"moments.omegak"+add;
            ofstream f0(name.c_str());
            f0.precision(8);
            for(int k=0; k<(int) omegak.size(); k++) f0 << k << " " << omegak[k] << endl;
            f0.close();
        }
        
        //===========================================================================
        // compute high precision isothermal approximation
        //===========================================================================
        else if(variable==1)
        {
            DInonrel[1]=xv;
            for(int k=0; k<(int) xv.size(); k++)
                DInonrel[1][k]=pow(xv[k], 3)*compute_SZ_signal_combo(xv[k], 
                                                                     TSZetc[variable][0], 
                                                                     TSZetc[variable][2], 
                                                                     betac, muc, betao, muo);
        }
    }
    
    //===============================================================================
    // find best-fit isothermal model
    //===============================================================================
    print_message("Computing best-fit isothermal model");
    vector<double> solution, Isol;
    find_best_fit_isothermal_model(xv, DIv[1], solution, Isol, 3);
    
    cout << "\n tau= " << TSZetc[0][0] << " tau_bestfit= " << solution[0] << endl;
    cout << " TSZ= " << TSZetc[0][2] << " TSZ_bestfit= " << solution[1] << endl;
    
    //===============================================================================
    //
    // output SZ signal
    //
    //===============================================================================
    ofstream ofile(fname.c_str());
    
    output_header(ofile, xcmb[0], xcmb[xcmb.size()-1], xcmb.size(), 
                  TSZetc[0][0], TSZetc[0][2], betac, muc, betao, muo, 
                  Te_order, beta_order, accuracy_goal, Te_max, 
                  "moments VIK cluster profiles"+mess);
    
    ofile.precision(16);
        
    for(int n=0; n<(int)xv.size(); n++)
        ofile << xv[n] << " " 
              << DIv[0][n] << " " << DIv[1][n] << " " 
              << DIviso[0][n] << " " << DIviso[1][n] << " " 
              << DInonrel[0][n] << " " << DInonrel[1][n] << " " 
              << DIsmooth[n] << " " << Isol[n] << endl;
    
    ofile.close();
    
    //===============================================================================
    // clean up
    //===============================================================================
    mv.clear();
    xv.clear();
    DIv.clear();
    DIviso.clear();
    TSZetc.clear();
    
    return;
}

//==================================================================================================
//
// computing SZ signal by explicitly integrating along the line of sight for each
// frequency bin. This routine is just for comparisons & rather slow.
//
//==================================================================================================
double dSZ_VIK(double l)
{
    double dtau=Ne_func(l)*const_sigT;
    double Te=Te_func(l);
    return compute_SZ_signal_combo(PSP.xcmb, dtau, Te, PSP.betac, PSP.muc, 0, 0);
}

double integrate_SZ_VIK(double xcmb)
{ 
    PSP.xcmb=xcmb;
    
    double r=Integrate_using_Patterson_adaptive(-1.0, -0.1, 1.0e-5, 1.0e-300, dSZ_VIK);
    r+=Integrate_using_Patterson_adaptive(-0.1, -0.01, 1.0e-5, 1.0e-300, dSZ_VIK);
    r+=Integrate_using_Patterson_adaptive(-0.01, -0.001, 1.0e-5, 1.0e-300, dSZ_VIK);
    r+=Integrate_using_Patterson_adaptive(-0.001, -0.0001, 1.0e-5, 1.0e-300, dSZ_VIK);
    r+=Integrate_using_Patterson_adaptive(-0.0001, -0.00001, 1.0e-5, 1.0e-300, dSZ_VIK);
    r+=Integrate_using_Patterson_adaptive(-0.00001, 0.00001, 1.0e-5, 1.0e-300, dSZ_VIK);
    r+=Integrate_using_Patterson_adaptive(0.00001, 0.0001, 1.0e-5, 1.0e-300, dSZ_VIK);
    r+=Integrate_using_Patterson_adaptive(0.0001, 0.001, 1.0e-5, 1.0e-300, dSZ_VIK);
    r+=Integrate_using_Patterson_adaptive(0.001, 0.01, 1.0e-5, 1.0e-300, dSZ_VIK);
    r+=Integrate_using_Patterson_adaptive(0.01, 0.1, 1.0e-5, 1.0e-300, dSZ_VIK);
    r+=Integrate_using_Patterson_adaptive(0.1, 1.0, 1.0e-5, 1.0e-300, dSZ_VIK);
   
    return PSP.zsc*PSP.Cluster_ptr->Get_rc_cm()*r;
}

//==================================================================================================
//
// computing temperature moments
//
//==================================================================================================
double dy_k_VIK(double l)
{
    double dtau=Ne_func(l)*const_sigT;
    double Te=Te_func(l), The=Te/const_me;
    return dtau*pow(The, PSP.k+1);
}

double integrate_yk_VIK(double xs, int k)
{ 
    PSP.xs=xs;
    PSP.k=k;
    
    double r=Integrate_using_Patterson_adaptive(-1.0, -0.1, 1.0e-5, 1.0e-300, dy_k_VIK);
    r+=Integrate_using_Patterson_adaptive(-0.1, -0.01, 1.0e-5, 1.0e-300, dy_k_VIK);
    r+=Integrate_using_Patterson_adaptive(-0.01, -0.001, 1.0e-5, 1.0e-300, dy_k_VIK);
    r+=Integrate_using_Patterson_adaptive(-0.001, -0.0001, 1.0e-5, 1.0e-300, dy_k_VIK);
    r+=Integrate_using_Patterson_adaptive(-0.0001, -0.00001, 1.0e-5, 1.0e-300, dy_k_VIK);
    r+=Integrate_using_Patterson_adaptive(-0.00001, 0.00001, 1.0e-5, 1.0e-300, dy_k_VIK);
    r+=Integrate_using_Patterson_adaptive(0.00001, 0.0001, 1.0e-5, 1.0e-300, dy_k_VIK);
    r+=Integrate_using_Patterson_adaptive(0.0001, 0.001, 1.0e-5, 1.0e-300, dy_k_VIK);
    r+=Integrate_using_Patterson_adaptive(0.001, 0.01, 1.0e-5, 1.0e-300, dy_k_VIK);
    r+=Integrate_using_Patterson_adaptive(0.01, 0.1, 1.0e-5, 1.0e-300, dy_k_VIK);
    r+=Integrate_using_Patterson_adaptive(0.1, 1.0, 1.0e-5, 1.0e-300, dy_k_VIK);
    
    return PSP.zsc*PSP.Cluster_ptr->Get_rc_cm()*r;
}

//==================================================================================================
void output_SZ_distortion_cluster_VIK_explicit(string fname,
                                               vector<double> &xa, 
                                               SZ_cluster_profiles &CL, 
                                               double x_rc, double y_rc,
                                               double betac, double muc, 
                                               double betao, double muo)
{
    print_message("Computing SZ signal explicitly as line-of-sight integral (slow...) ");
    
    //===============================================================================
    PSP.betac=betac;
    PSP.muc=muc;
    PSP.Cluster_ptr=&CL;

    PSP.xs=x_rc;
    PSP.ys=y_rc;
    PSP.zsc=1000.0;
    
    //===============================================================================
    int np=xa.size();
    vector<double> Ia(np);

    for(int k=0; k<np ; k++) Ia[k]=pow(xa[k], 3)*integrate_SZ_VIK(xa[k]);
    
    print_message("Computing best-fit isothermal model");
    vector<double> solution, Isol;
    find_best_fit_isothermal_model(xa, Ia, solution, Isol);
    
    double tau=integrate_yk_VIK(PSP.xs, -1);
    double ySZ=integrate_yk_VIK(PSP.xs, 0);
    double TSZ=ySZ*const_me/tau;
    
    //===============================================================================
    ofstream ofile(fname.c_str());
    
    output_header(ofile, xa[0], xa[np-1], np, 0, 0, 
                  betac, muc, betao, muo, 
                  0, 0, 0, 0, "explicit");
    
    ofile.precision(16);
    
    for(int k=0; k<np ; k++)
    {
        ofile << xa[k] << " " << Ia[k] << " " << Isol[k] << endl;
        cout  << " x= " << xa[k] << " " << Ia[k] << " " << Isol[k] << endl;
    }
    
    ofile.close();

    cout << "\n tau= " << tau << " tau_bestfit= " << solution[0] << endl;
    cout << " TSZ= " << TSZ << " TSZ_bestfit= " << solution[1] << endl;
    
    return;
}

//==================================================================================================
void write_yk_ykiso(double xs, ofstream &f0, ofstream &f1, ofstream &f2, int kmax)
{
    double tau=integrate_yk_VIK(xs, -1);
    double ySZ=integrate_yk_VIK(xs, 0);
    double TSZ=ySZ*const_me/tau;

    f0 << xs << " " << tau << " " << TSZ << " " << ySZ << " ";
    f1 << xs << " " << tau << " " << TSZ << " " << ySZ << " ";
    f2 << xs << " " << tau << " " << TSZ << " " << ySZ << " ";

    vector<double> rhok(kmax+1, 1.0);
    
    for(int m=1; m<=kmax; m++)
    {
        double yk   =integrate_yk_VIK(xs, m);
        double ykiso=ySZ*pow(TSZ/const_me, m);
        rhok[m] =yk/ykiso;
        
        f0 << yk   << " ";
        f1 << rhok[m] << " ";
    }
    
    for(int m=0; m<=kmax; m++)
    {
        double omegak=pow(-1.0, m+1);
        for(int l=1; l<=m+1; l++)
            omegak+=Binomial_coeff(m+1, l)*pow(-1.0, m+1-l)*rhok[l-1];
        
        f2 << omegak << " ";
    }

    f0 << endl;
    f1 << endl;
    f2 << endl;

    return;
}

void output_SZ_distortion_cluster_yk(string root, string add,
                                     SZ_cluster_profiles &CL, 
                                     double y_rc, int kmax, int np)
{
    //===============================================================================
    PSP.Cluster_ptr=&CL;
    
    PSP.ys=y_rc;
    PSP.zsc=1000.0;
    
    //===============================================================================
    vector<double> lz(np); 
    init_xarr(1.0e-5, 1, &lz[0], np, 1, 0); 

    string fname=root+"yk"+add;
    ofstream ofile(fname.c_str());
    ofile.precision(16);

    fname=root+"rhok"+add;
    ofstream ofile1(fname.c_str());
    ofile1.precision(16);
    
    fname=root+"wk"+add;
    ofstream ofile2(fname.c_str());
    ofile2.precision(16);
    
    double xs;
    
/*    for(int k=np-1; k>=0; k--) 
    {
        xs=-PSP.zsc*lz[k];
        write_yk_ykiso(xs, ofile, ofile1, ofile2, kmax);
    }
*/    
    xs=0.0;
    write_yk_ykiso(xs, ofile, ofile1, ofile2, kmax);
    
    for(int k=0; k<np; k++)
    {
        xs=PSP.zsc*lz[k];
        write_yk_ykiso(xs, ofile, ofile1, ofile2, kmax);
    }
    
    ofile.close();
    ofile1.close();
    ofile2.close();
    
    return;
}


//==================================================================================================
void output_TSZ_bestfit(string root, string add,
                        vector<double> &xa, 
                        SZ_cluster_profiles &CL, 
                        SZ_moment_method &SZM,
                        double betac, double muc,
                        double x_rc_max, int np)
{
    //===============================================================================
    PSP.Cluster_ptr=&CL;
    PSP.ys=0.0;
    PSP.zsc=1000.0;
    PSP.betac=betac;
    PSP.muc=muc;
    
    //===============================================================================
    vector<double> lz(np); 
    init_xarr(1.0e-5, x_rc_max/PSP.zsc, &lz[0], np, 1, 0); 
    
    string fname=root+"TSZ.bestfit"+add;
    ofstream ofile(fname.c_str());
    ofile.precision(16);
    
    double xs, tau, ySZ, TSZ;

    int nx=xa.size();
    vector<double> Ia(nx), mv, TSZetc;
    vector<double> solution, Isol;
    
    for(int k=0; k<np; k++)
    {
        xs=PSP.zsc*lz[k];
        tau=integrate_yk_VIK(xs, -1);
        ySZ=integrate_yk_VIK(xs, 0);
        TSZ=ySZ*const_me/tau;
        
        ofile << xs << " " << tau << " " << TSZ << " ";

        //===========================================================================
//        for(int k=0; k<nx; k++)  Ia[k]=pow(xa[k], 3)*integrate_SZ_VIK(xa[k]);

        SZM.compute_SZ_signal_moments(Ne_func, Te_func, betac_func, muc_func,  
                                      PSP.zsc*PSP.Cluster_ptr->Get_rc_cm(), mv, xa, 
                                      Ia, TSZetc, 1, 2.0e-6);
        
        find_best_fit_isothermal_model(xa, Ia, solution, Isol);
        
        //===========================================================================
        // estimate best-fit case
        //===========================================================================
        double w1=integrate_yk_VIK(xs, 1)/ySZ/(TSZ/const_me)-1.0;
        double est_tau=tau*( 1.0-w1/( 1.0+0.027*pow(TSZ, 0.86)) );
        double est_Te =TSZ*( 1.0+w1*exp(-0.026*pow(TSZ, 0.86)) );
        
        ofile << " " << solution[0] << " " << solution[1] << " " << Te_func(lz[k])
              << " " << est_tau << " " << est_Te << endl;
    }
            
    ofile.close();
    
    return;
}

//==================================================================================================
void output_SZ_distortion_derivatives(string fname, 
                                      vector<double> &xa, 
                                      double Dtau, double Te, 
                                      double betac, double muc, int kmax)
{
    print_message("Computing derivatives of S in CMB rest frame");
    
    int np=xa.size();
    
    ofstream ofile(fname.c_str());
    
//    output_header(ofile, xmin, xmax, np, Dtau, Te, betac, muc, betao, muo, 0, 0, 0, "combo");
    
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
        double Ixp=pow(xa[k]*(1.0+0.001), 3)*compute_SZ_signal_combo(xa[k]*(1.0+0.001), 1.0, Te, 0, 0, 0, 0);
        double Ixm=pow(xa[k]*(1.0-0.001), 3)*compute_SZ_signal_combo(xa[k]*(1.0-0.001), 1.0, Te, 0, 0, 0, 0);
        double dI=(Ixp-Ixm)/(2.0*0.001);
        double d2I=(Ixp-2.0*Ix+Ixm)/pow(0.001, 2) / 2.0;
        
        ofile << dI << " " << d2I << endl;
        
        cout  << " x= " << xa[k] << " " << dDn_dThe[0]*x3 << endl;
    }
    
    ofile.close();
    
    return;
}


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
//    output_header(ofile, xmin, xmax, np, Dtau, Te, betac, muc, betao, muo, 0, 0, 0, "combo");
    ofile.precision(16);
    
    double DI_s=4.881e-5; // fiducial precision
    int betao=0, muo=0;
    
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
// Compute degeneracy factors
//
//==================================================================================================
void compute_degeneracy_functions(string fname, double xmin, double xmax, int np, 
                                  double betac, double muc)
{
    print_message("Computing the degeneracy coefficients");
    
    double Dtau=1.0;
    
    vector<double> xa(np);
    init_xarr(xmin, xmax, &xa[0], np, 1, 0); // change 1 --> 0 to have linear grid in x
    
    int nT=100;
    vector<double> Ta(nT);
    init_xarr(0.01, 50.0, &Ta[0], nT, 1, 0); // change 1 --> 0 to have linear grid in x

    ofstream ofile(fname.c_str());
    ofile.precision(16);
    
    vector<vector<double> > a_v(3, vector<double>(np, 0));
    vector<vector<double> > s_v(4, vector<double>(np, 0));
    vector<double> dDn_dThe, M(3*3, 0.0), A(4*3, 0.0);
        
    gsl_permutation * p = gsl_permutation_alloc (3);
    gsl_vector *b = gsl_vector_alloc (3);
    gsl_vector *x_GSL = gsl_vector_alloc (3);
    
    for(int m=0; m<nT; m++)
    {
        double Te=Ta[m];
        
        for(int k=0; k<np; k++)
        {
            double x3=pow(xa[k], 3);
            
            Dcompute_SZ_signal_combo_CMB(xa[k], 2, 0, 0, Dtau, Te, betac, muc, dDn_dThe);
            a_v[0][k]=dDn_dThe[0]*x3;
            a_v[1][k]=dDn_dThe[1]*x3;
            //
            s_v[0][k]=dDn_dThe[2]*x3;
            
            Dcompute_SZ_signal_combo_CMB(xa[k], 1, 1, 0, Dtau, Te, betac, muc, dDn_dThe);
            a_v[2][k]=dDn_dThe[0]*x3*0.01;
            //
            s_v[1][k]=dDn_dThe[1]*x3*0.01;

            Dcompute_SZ_signal_combo_CMB(xa[k], 0, 2, 0, Dtau, Te, betac, muc, dDn_dThe);
            s_v[2][k]=dDn_dThe[0]*x3*0.01*0.01;

            Dcompute_SZ_signal_combo_CMB(xa[k], 0, 0, 1, Dtau, Te, betac, muc, dDn_dThe);
            s_v[3][k]=dDn_dThe[0]*x3*0.01*0.01;
         }
        
        //===============================================================================
        // compute signal matrices
        //===============================================================================
        for(int c=0; c<3; c++) 
            for(int r=0; r<3; r++) 
                M[c+3*r]=DC_sumprod(a_v[c], a_v[r]);

        for(int c=0; c<4; c++) 
            for(int r=0; r<3; r++) 
                A[c+4*r]=DC_sumprod(s_v[c], a_v[r]);
        
        //===============================================================================
        // no velocity vector present
        //===============================================================================
        double det2 =    M[0+3*1]*M[0+3*1] - M[0+3*0]*M[1+3*1];
        double alpha= -( M[0+3*1]*A[0+4*1] - M[1+3*1]*A[0+4*0] ) / det2; 
        double beta =  ( M[0+3*1]*A[0+4*0] - M[0+3*0]*A[0+4*1] ) / det2; 
        
        ofile << Te << " " << alpha << " " << beta << " ";

        //===============================================================================
        // velocity vector present and only tau & Te varied
        //===============================================================================
        det2 =    M[0+3*1]*M[0+3*1] - M[0+3*0]*M[1+3*1];
        alpha= -( M[0+3*1]*M[2+3*1] - M[1+3*1]*M[2+3*0] ) / det2; 
        beta =  ( M[0+3*1]*M[2+3*0] - M[0+3*0]*M[2+3*1] ) / det2; 
                          
        ofile << alpha << " " << beta << " ";
        
        //===============================================================================
        // with velocity present
        //===============================================================================
        gsl_matrix_view mm=gsl_matrix_view_array (&M[0], 3, 3);
        int s;
        gsl_linalg_LU_decomp (&mm.matrix, p, &s);
        
        for(int c=0; c<4; c++) 
        {
            for(int r=0; r<3; r++) b->data[r]=A[c+4*r];

            gsl_linalg_LU_solve (&mm.matrix, p, b, x_GSL);
            
            for(int k=0; k<3; k++) ofile << x_GSL->data[k] << " ";
        }
        
        ofile << endl;
    }
    
    gsl_permutation_free (p);
    gsl_vector_free (x_GSL);    
    gsl_vector_free (b);    

    ofile.close();
    
    return;
}

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
        double Ixp=pow(x0*(1.0+0.001), 3)*compute_SZ_signal_combo(x0*(1.0+0.001), 1.0, Te, 0, 0, 0, 0);
        double Ixm=pow(x0*(1.0-0.001), 3)*compute_SZ_signal_combo(x0*(1.0-0.001), 1.0, Te, 0, 0, 0, 0);
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
// This function can be used to compute the precision of the different expansions. The results can 
// be compared with the full numerical result or the approximation obtained with the combo function.
// For Te < 75 eV the combo function should be sufficient as reference, but above currently the 
// integrals have to be carried out explicitly.
//
//==================================================================================================
void compute_precision_of_basis(string fname, double xmin, double xmax, int accuracy_goal)
{
    print_message("Computing the precision of the basis");
    
    int np=200, nT=40;
    vector<double> xa(np), ya_3D(np), ya(np);
    vector<double> Tea(nT);

    init_xarr(xmin, xmax, &xa[0], np, 1, 0);  // change 1 --> 0 to have linear grid in x
    init_xarr(1.0, 4.0, &Tea[0], nT, 0, 0);  // change 1 --> 0 to have linear grid in x
    
    ofstream ofile(fname.c_str());
    ofile.precision(16);
    
    double Dtau=0.01, DI_s=4.881e-5; // fiducial precision
    double betac=0.01, muc=1.0, betao=0.0, muo=1.0;
    
    for(int t=0; t<nT; t++)
    {
        double Te=Tea[t];

        ofile << Te << " ";
        
        //===========================================================================
        // compute precise solution
        //===========================================================================
        for(int k=0; k<np; k++)             
//            ya_3D[k]=pow(xa[k], 3)*compute_SZ_signal_3D(xa[k], Dtau, Te, betac, muc,  
//                                                        betao, muo, 1.0e-6);        
            ya_3D[k]=pow(xa[k], 3)*compute_SZ_signal_combo(xa[k], Dtau, Te, 
                                                           betac, muc, betao, muo);

        //===========================================================================
        // compute approximation
        //===========================================================================
        for(int m=6; m<=6; m++)
        {
            for(int k=0; k<np; k++) 

                ya[k]=pow(xa[k], 3)*Dtau*compute_SZ_distortion_CNSN_basis_opt(xa[k], Te/const_me, 
                                                                              betac, muc, 
                                                                              m, 2, "all", 
                                                                              accuracy_goal);

            //=======================================================================
            // find maximal absolute deviation and save to file
            //=======================================================================
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
// read parameters
//
//==================================================================================================
void read_parameters(string fname, 
                     double &xmin, double &xmax, int &np, 
                     double &Dtau, double &Te, 
                     double &betac, double &muc, 
                     double &betao, double &muo,
                     int &Te_order, int &beta_order, 
                     double &Te_max, int &accuracy_goal,
                     string &outpath, string &add)
{
    ifstream ifile(fname.c_str());
    
    ifile >> xmin; ifile >> xmax; ifile >> np;
    ifile >> Dtau; ifile >> Te;
    ifile >> betac; ifile >> muc;
    ifile >> betao; ifile >> muo;
    ifile >> Te_order; ifile >> beta_order; 
    ifile >> Te_max; ifile >> accuracy_goal;
    ifile >> outpath; ifile >> add;
    
    //===============================================================================
    // simple sanity checks
    //===============================================================================
    if(xmin<0.01 || xmax>50.0){ cerr << " read_parameters :: change x-range " << endl; exit(0); }
    if(np<=0){ cerr << " read_parameters :: change number of frequency points " << endl; exit(0); }

    if(Dtau<0){ cerr << " read_parameters :: change Dtau " << endl; exit(0); }

    if(betac<0 || betac>0.1){ cerr << " read_parameters :: check betac " << endl; exit(0); }
    if(muc<-1.0 || muc>1){ cerr << " read_parameters :: check muc " << endl; exit(0); }
    if(betao<0 || betao>0.1){ cerr << " read_parameters :: check betao " << endl; exit(0); }
    if(muo<-1.0 || muo>1){ cerr << " read_parameters :: check muo " << endl; exit(0); }

    if(Te_max<0){cerr << " read_parameters :: check Te_max "<< endl;exit(0);}

    return; 
}
 
//==================================================================================================
//
// main call of routine. Here parameters are handled.
//
//==================================================================================================
int main(int narg, char *args[])
{
    //===============================================================================
    // example parameters
    //===============================================================================
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
    int Te_order=3;
    int beta_order=2;
    //
    int accuracy_goal=0;
    double Te_max=30.0;
    //
    string outpath=SZPACKDIR+"./outputs/";
    string add=".dat";
    string mode="ISO";
    
    //===============================================================================
    // startup with different number of arguments
    //===============================================================================
    if(narg==1){}
    
    else if(narg==2) mode=args[narg-1];
    
    else if(narg==3)
    {
        mode=args[narg-2];
        string fname=args[narg-1];
        
        print_message("Reading parameters from parameter file: " + fname);
        
        read_parameters(fname, xmin, xmax, np, Dtau, Te, 
                        betac, muc, betao, muo,
                        Te_order, beta_order, 
                        Te_max, accuracy_goal, 
                        outpath, add);
    }

    else{ cerr << " Too many/few parameters " << endl; exit(1); }
        
    //===============================================================================
    // setting up frequency points. Here in principle any frequency 
    // distribution can be used.
    //===============================================================================
    vector<double> xcmb(np); 
    init_xarr(xmin, xmax, &xcmb[0], np, 1, 0); // change 1 --> 0 to have linear grid in x
    
    //===============================================================================
    // combination of asymptotic expansion + basis functions of CNSN2012
    // This function should give very precise results for Te < 75keV and
    // can be used for comparison.
    //===============================================================================
    if(mode=="COMBO")
    {
        output_SZ_distortion_combo(outpath+"SZ_combo"+add,
                                   xcmb, Dtau, Te, 
                                   betac, muc, betao, muo);
        return 0;
    }
    
    //===============================================================================
    // check the precision of the basis function. Settings should be 
    // changed in compute_precision_of_basis() above
    //===============================================================================
    else if(mode=="ACC")
    {
        compute_precision_of_basis(outpath+"rec.moments.acc_IV.beta"+add, 
                                   xmin, xmax, accuracy_goal=3);
        return 0;
    }
    
    //===============================================================================
    // combination of asymptotic expansion + basis functions of CNSN2012
    // This function should give very precise results for Te < 75keV and
    // can be used for comparison.
    //===============================================================================
    else if(mode=="MEANS")
    {
        output_SZ_distortion_means(outpath+"SZ_means.III"+add,
                                   xcmb, Dtau, Te, 
                                   betac, muc, betao, muo);
        return 0;
    }

    //===============================================================================
    //
    // matrix setup for moment method
    //
    //===============================================================================
    if(mode=="ISO" || mode=="ISOOPT") 
        add=".acc_"+int_to_string(accuracy_goal)+".ISO"+add;
    
    else if(mode=="VIK" || mode=="VIKOPT") 
        add=".acc_"+int_to_string(accuracy_goal)+".VIK"+add;
    
    //===============================================================================
    // This version uses the accuracy settings according to CSNN 2012
    // for given kmax==Te_order
    //===============================================================================
    SZ_moment_method SZMoments(xcmb, Te_order, accuracy_goal, beta_order);
    
    //===============================================================================
    // In this case the optimal kmax is determined internally according 
    // to the required accuracy goal and maximal temperature. This 
    // version of the moment matrix minimizes the # of parameters needed.
    // Calling the functions of 'SZMomentsopt' might give exactly the  
    // same as those of 'SZMoments'! This depends on the settings.
    //===============================================================================
    SZ_moment_method SZMomentsopt(xcmb, Te_max, accuracy_goal, beta_order);
    
    //===============================================================================
    //
    // store the matrices corresponding to the different moment versions
    // of variables for later computations. These are actually all one 
    // needs to compute S = F * m !
    //
    //===============================================================================
    if(mode=="MATRIX")
    {
        print_message("Storing SZ signal matrices");
        
        SZMoments.export_x_vector(outpath+"vector.x"+add);
    
        SZMoments.export_moment_matrix_M(outpath+"matrix.M"+add);
        SZMoments.export_moment_matrix_MT(outpath+"matrix.MT"+add);
    
        SZMomentsopt.export_moment_matrix_M(outpath+"matrix.M.opt"+add);
        SZMomentsopt.export_moment_matrix_MT(outpath+"matrix.MT.opt"+add);

        return 0;
    }
    
    //===============================================================================
    // call in the different runmodes
    //===============================================================================
    if(mode=="ISO")
    {
        output_SZ_distortion_moments_isothermal(outpath+"SZ_moments_isothermal"+add, 
                                                outpath, add, xcmb, Dtau, Te, 
                                                betac, muc, betao, muo, SZMoments,
                                                Te_order, beta_order, 
                                                accuracy_goal, Te_max);
        return 0;
    }
    
    else if(mode=="ISOOPT")
    {
        add=".opt"+add;
        output_SZ_distortion_moments_isothermal(outpath+"SZ_moments_isothermal"+add, 
                                                outpath, add, xcmb, Dtau, Te, 
                                                betac, muc, betao, muo, SZMomentsopt,
                                                Te_order, beta_order, 
                                                accuracy_goal, Te_max, " (opt)");
        return 0;
    }

    //==================================================================
    // compute null of SZ signal
    //==================================================================
    else if(mode=="NULL")
    {
        output_SZ_null(outpath+"SZ_null"+add);
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
    add="."+CL.Get_label()+add;
    
    //===============================================================================
    // compute yk temperature moments for some slices. 
    //===============================================================================
    if(mode=="VIKMOM")
    {
        print_message("Computing yk-integrals");
        
        output_SZ_distortion_cluster_yk(outpath+"SZ_moments_", ".y_0kpc"+add,
                                        CL, 0.0001/CL.Get_rc(), 10, 200);

        output_SZ_distortion_cluster_yk(outpath+"SZ_moments_", ".y_10kpc"+add,
                                        CL, 10.0/CL.Get_rc(), 10, 100);
        
        output_SZ_distortion_cluster_yk(outpath+"SZ_moments_", ".y_50kpc"+add,
                                        CL, 50.0/CL.Get_rc(), 10, 100);
        
        output_SZ_distortion_cluster_yk(outpath+"SZ_moments_", ".y_100kpc"+add,
                                        CL, 100.0/CL.Get_rc(), 10, 100);
                
        return 0;
    }
    
    //===============================================================================
    // compute derivatives of SZ signal 
    //===============================================================================
    if(mode=="DERIVS")
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

        return 0;
    }
    
    //===============================================================================
    // SZ signal for two-temperature case
    //===============================================================================
    if(mode=="TWOT")
    {
        double ftau=0.2;
        double Te1=5;
        double Delta=0.5;
        
        output_SZ_distortion_two_temp(outpath+"SZ_moments.two-temp.ftau_0.2.Te_5.Delta_0.5"+add, 
                                      xcmb, Dtau, ftau, Te1, Delta, 0*betac, muc);
         
        return 0;
    }
    
    //===============================================================================
    // compute derivatives of SZ signal 
    //===============================================================================
    if(mode=="DEGEN")
    {
        compute_degeneracy_functions(outpath+"SZ_moments_degeneracies.bc_0.0"+add, 
                                     0.1, 30.0, 200, 0.0, muc);
        
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
    add=".y_1rc.S1"+add; 
    
/*    output_TSZ_bestfit(outpath+"SZ_moments_", add, xcmb, CL, SZMoments, 
                       betac, muc, 10.0, 40);
    
    exit(0);
*/  
    //===============================================================================
    // compute SZ signal for different cluster profile according to 
    // fits of Vikhlinin et al 2006
    //===============================================================================
    if(mode=="VIK")
    {
        output_SZ_distortion_moments_cluster_VIK(outpath+"SZ_moments_clusters"+add, 
                                                 outpath, add, xcmb,
                                                 CL, xs, ys,
                                                 betac, muc, betao, muo, SZMoments,
                                                 Te_order, beta_order, 
                                                 accuracy_goal, Te_max);
    }
        
    else if(mode=="VIKOPT")
    {
        add=".opt"+add;
        output_SZ_distortion_moments_cluster_VIK(outpath+"SZ_moments_clusters"+add, 
                                                 outpath, add, xcmb, 
                                                 CL, xs, ys,
                                                 betac, muc, betao, muo, SZMomentsopt,
                                                 Te_order, beta_order, 
                                                 accuracy_goal, Te_max, " (opt)");
    }

    //===============================================================================
    // computing SZ signal by explicitly integrating along the line of sight for each
    // frequency bin. This routine is just for comparisons & rather slow.
    //===============================================================================
    else if(mode=="VIKEXP")
    {
        add=".explicit"+add;
        output_SZ_distortion_cluster_VIK_explicit(outpath+"SZ_moments_clusters"+add, 
                                                  xcmb, CL, xs, ys,
                                                  betac, muc, betao, muo);
    }
    
    return 0;
}

//==================================================================================================
//==================================================================================================
