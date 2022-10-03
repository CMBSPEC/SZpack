//==================================================================================================
//
// This program allows computing the thermal SZ and kinematic effect using the temperature-velocity 
// moment formalism described by CSNN 2012. 
//
//==================================================================================================
//
// setup the moment matrix. This function has to be called prior to the computations. The SZ 
// signal can then be computed in the form S = M * m. The accuracy levels are defines in XXX.  
// For a fixed maximal temperature the optimal kmax is determined internally.
//
//==================================================================================================
//
// Author: Jens Chluba  (CITA, University of Toronto)
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
#include "routines.h"
#include "Integration_routines.h"
#include "nPl_derivatives.h"
#include "SZ_moment_method.h"
#include "SZ_asymptotic.h"
#include "SZ_CNSN_basis.opt.h"
#include "SZpack.h"

//==================================================================================================
//
// namespaces
//
//==================================================================================================
using namespace std;

//==================================================================================================
double G_func(double x)
{ return x*exp(-x)/(1.0-exp(-x))/(1.0-exp(-x)); }

double Q_func(double x)
{ return x*exp(-x)/(1.0-exp(-x))/(1.0-exp(-x))*x*(1.0+exp(-x))/(1.0-exp(-x)); }

double M_func(double x)
{ return x*exp(-x)/(1.0-exp(-x))/(1.0-exp(-x))*(x*(1.0+exp(-x))/(1.0-exp(-x))-3.0); }

//==================================================================================================
void sanity_check(int &order_T_low, int &order_T_high, int &order_b)
{
    order_T_low = (int)min(order_T_low, 10);
    order_T_high = (int)min(order_T_high, 20);
    order_b = (int)min(order_b, 2);
    
    return;
}

//==================================================================================================
//
// setup the moment matrix. This function has to be called prior to the computations. The SZ signal
// can then be computed in the form S = M * m 
//
//==================================================================================================
void computeY(double x, vector<double> Y, int region = 0, bool CNSN = false){ //TODO: These shouldn't really be here!
    if (CNSN){
        return compute_Y_asymptotic(x, Y);
    }
    return compute_Y_CNSNopt(x, region, Y);
}

void computeD(double x, vector<double> Y, int region = 0, bool CNSN = false){
    if (CNSN){
        return compute_D_asymptotic(x, Y);
    }
    return compute_D_CNSNopt(x, region, Y);
}

void computeM(double x, vector<double> Y, int region = 0, bool CNSN = false){
    if (CNSN){
        return compute_Mcorr_asymptotic(x, Y);
    }
    return compute_Mcorr_CNSNopt(x, region, Y);
}

void computeQ(double x, vector<double> Y, int region = 0, bool CNSN = false){
    if (CNSN){
        return compute_Q_asymptotic(x, Y);
    }
    return compute_Q_CNSNopt(x, region, Y);
}

void SZ_moment_method::calculateGasMatrix(Parameters fp, int n, vector<double> &Y, int order_T, double x3, int &k, int region, bool CNSN){
    computeY(fp.xcmb[n], Y, region, CNSN);
    for(int l = 0; l <= order_T; l++){
        M.setValue(n, l + k, x3*Y[l]);
    }

    if(fp.beta_order > 0)
    {
        k += 1 + order_T;
        M.setValue(n, k, x3*G_func(fp.xcmb[n]));
        computeD(fp.xcmb[n], Y, region, CNSN);
        k++;
        for(int l = 0; l <= order_T; l++){
            M.setValue(n, l + k, x3*Y[l]);
        }

        if(fp.beta_order>1)
        {
            k += 1+order_T;
            M.setValue(n, k, x3*M_func(fp.xcmb[n])/3.0);
            computeM(fp.xcmb[n], Y, region, CNSN);
            k++;
            for(int l = 0; l <= order_T; l++){
                M.setValue(n, l + k, x3*Y[l]);
            }
            k += 1+order_T;
            
            M.setValue(n, k, x3*11.0/30.0*Q_func(fp.xcmb[n]));
            computeQ(fp.xcmb[n], Y, region, CNSN);
            k++;
            for(int l = 0; l <= order_T; l++){
                M.setValue(n, l + k, x3*Y[l]);
            }
        }
    }
}

void SZ_moment_method::setup_SZ_moment_matrix(int order_T_low, int order_T_high, Parameters fp)
{
    sanity_check(order_T_low, order_T_high, fp.beta_order);
    //Replaces low with 10 and high with 20 and beta_order with 2 if any value is bigger
    
    //==============================================================================================
    // save info
    //==============================================================================================
    index_high = 2 + order_T_low;
    nmom_high_one = 1 + order_T_high;
    
    if(fp.beta_order > 0)
    { 
        index_high    *= 2;                      // asymptotic expansion part
        nmom_high_one *= 2;                      // CNSN expansion part

        if(fp.beta_order > 1)
        {
            index_high    *= 2;                  // asymptotic expansion part
            nmom_high_one *= 2;                  // CNSN expansion part
        }
    }
    index_high -= 1;                            // asymptotic expansion part
    tot_moments = index_high+(nTregs-1)*nmom_high_one;

    //==============================================================================================
    // fill matrix
    //==============================================================================================
    vector<double> Y(order_T_low+1);
    M.setSize(nx, tot_moments);
    int ind;
    
    for(int n = 0; n < nx; n++)
    {
        int k = 0;
        //==========================================================================================
        // description of low temperature gas
        //==========================================================================================
        double x3 = pow(fp.xcmb[n], 3);

        calculateGasMatrix(fp, n, Y, order_T_low, x3, k);

        //==========================================================================================
        // description of high temperature gas
        //==========================================================================================
        if(nTregs > 1)
        {
            vector<double> Yb(order_T_high+1);

            for(int r = 1; r < nTregs; r++)
            {
                int ri = region_indices[r];
                calculateGasMatrix(fp, n, Yb, order_T_high, x3, k, ri, true);
            }
        }
    }
        

    //==============================================================================================
    // fill matrix for transformation from N(The) (The-The0)^k --> N(The) The^k
    //==============================================================================================
    T.setSize(tot_moments,tot_moments);
    
    for(int m = 0; m < index_high; m++) {
        T.setValue(m, m, 1.0);
    }
    
    for(int m = 0; m <= order_T_high && nTregs>1; m++) 
    {
        for(int k = 0; k <= m; k++) 
        {
            ind = index_high;
            
            for(int r = 1; r < nTregs; r++)
            {
                double Tij = Binomial_coeff(m, k)*pow(-Get_The_ref(r), m-k);

                T.setValue(ind+k, ind+m, Tij);
                ind += order_T_high+1;
                
                if(fp.beta_order > 0)
                {
                    T.setValue(ind+k, ind+m, Tij);
                    ind += order_T_high+1;
                    
                    if(fp.beta_order > 1)
                    {
                        T.setValue(ind+k, ind+m, Tij);
                        ind += order_T_high+1;
                        T.setValue(ind+k, ind+m, Tij);
                        ind += order_T_high+1;
                    }
                }
            }
        }
    }
    
    //==============================================================================================
    // compute transformed matrix MT = M*T
    //==============================================================================================
    MT.setSize(nx, tot_moments);
    for(int n=0; n<nx; n++)
    {
        for(int m = 0; m < tot_moments; m++)
        {
            double r=0.0;
            for(int k = 0; k < tot_moments; k++) {
                r += M.getValue(n,k)*T.getValue(m,k);
            }
            MT.setValue(n, m, r);
        }
    }
//    show_moment_matrix_M();
//    show_moment_matrix_MT();
//    show_conversion_matrix_T();     
    
    return;
}



//==================================================================================================
//
// show data of matrices
//
//==================================================================================================
SZ_moment_method::SZ_moment_method(){}

//==================================================================================================
SZ_moment_method::SZ_moment_method(Parameters fp, bool usekmax)
{
    nx = fp.gridpoints;
    kmax = fp.kmax;
    if (!usekmax){
        determine_optimal_kmax(fp.accuracy_level, fp.Te_max, kmax, nTregs);
    }
    Te_pivots = Get_temperature_pivots (kmax, fp.accuracy_level);
    Te_limits = Get_temperature_regions(kmax, fp.accuracy_level);
    region_indices = Get_region_indices(kmax, fp.accuracy_level);

    if (usekmax){
        nTregs = Te_limits.size();
        fp.Te_max = Te_limits.back();
    }

    setup_SZ_moment_matrix(kmax, kmax, fp);
    
    show_Temperature_limits();
    show_Temperature_pivots();
}

//==================================================================================================
SZ_moment_method::~SZ_moment_method(){}

//==================================================================================================
//
// simple functions 
//
//==================================================================================================
int SZ_moment_method::Get_index_of_Te_region(double Te)  // if Te > Temax --> -10 is returned
{
    if(Te>Te_limits[nTregs-1]) return -10;
    if(Te<Te_limits[0])        return 0;

    int reg=0;
    do reg++;
    while(Te>Te_limits[reg]);
    
    return reg;        
}

int SZ_moment_method::Get_start_index_of_moment_set(int reg)
{
    if(reg==0) return 0;
    return index_high+(reg-1)*Get_number_of_moments_for_one_high_region();
}

double SZ_moment_method::Get_Te_ref(int reg){ return Te_pivots[reg]*const_me; }
double SZ_moment_method::Get_The_ref(int reg){ return Te_pivots[reg]; }

//==================================================================================================
//
// show moment data
//
//==================================================================================================
void SZ_moment_method::show_Temperature_limits()
{
    cout << "\n SZ_moment_method::show_Temperature_limits: The temperature limits are:" << endl;
    cout << " " << 0 << "\t - \t" << Te_limits[0] << "\tkeV (Asymptotic expansion) " << endl;    

    for(int n=1; n<nTregs; n++) 
        cout << " " << Te_limits[n-1] << "\t - \t" << Te_limits[n] << "\tkeV (CNSN basis)" << endl; 
    
    cout << " Total number of Te-regions " << nTregs << endl;
    
    return;
}

//==================================================================================================
void SZ_moment_method::show_Temperature_pivots()
{
    cout << "\n SZ_moment_method::show_Temperature_pivots: The pivot temperatures are:" << endl;
    cout << " Region " << 0 << ": " << Te_pivots[0]*const_me << "\tkeV (Asymptotic expansion)" << endl;    

    for(int n=1; n<nTregs; n++) 
        cout << " Region " << n << ": " << Te_pivots[n]*const_me << "\tkeV (CNSN basis)" << endl;    

    return;
}

//==================================================================================================
//
// show data of matrices
//
//==================================================================================================
//==================================================================================================
void SZ_moment_method::show_matrix(simpleMatrix &B, Parameters fp, bool show_xcmb)
{
    for(int n=0; n < B.width; n++)
    {
        if (show_xcmb){
            cout << fp.xcmb[n] << " || ";
        }
        
        for(int m=0; m < B.height; m++)
        {
            if(m==index_high) cout << "|| ";
            cout << B.getValue(n, m) << " ";
        }
        
        cout << endl;
    }
    cout << endl;
    return;
}

void SZ_moment_method::show_moment_matrix_M(Parameters fp)
{
    show_matrix(M, fp);
    return;
}

void SZ_moment_method::show_moment_matrix_MT(Parameters fp)
{
    show_matrix(MT, fp);
    return;
}

//==================================================================================================
void SZ_moment_method::show_conversion_matrix_T(Parameters fp)
{
    show_matrix(T, fp, false);
    return;
}

//==================================================================================================

void SZ_moment_method::determine_Te_structure(double lres) //Uses 1d cluster functions
{
    Te_zeros.clear();
    Te_zeros.resize(nTregs+1);
    Te_region dum;
    
    int nres=floor(1.0/lres);
    
    double dx=2.0/(nres-1);
    
    double lmin=-1.0;

    double Te = clusterF_1D.Te(lmin, PSP);
    int ir = Get_index_of_Te_region(Te), irp;
    dum.lmin = lmin;
    
    for(int m=1; m<nres; m++)
    {
        lmin+=dx;
        Te=clusterF_1D.Te(lmin, PSP);
        irp=Get_index_of_Te_region(Te);
        
        if(ir!=irp)
        {
            dum.lmax=lmin;
            Te_zeros[ir].push_back(dum);
            dum.lmin=lmin;
            ir=irp;
        }        
    }
    dum.lmax = lmin;
    Te_zeros[ir].push_back(dum);

    dum.lmin = dum.lmax=0.0;
    for(int m=0; m<nTregs; m++) {
        if(Te_zeros[m].size()==0) Te_zeros[m].push_back(dum);
    }
    
/*    for(int m=0; m<nTregs; m++)
    {
        for(int z=0; z<(int)Te_zeros[m].size(); z++)
            cout << " Region: " << m << " " << z 
                 << " l-range: " << Te_zeros[m][z].lmin 
                 << " " << Te_zeros[m][z].lmax << endl;
    }
    wait_f_r();
*/    
    return;
}

//==================================================================================================
//==================================================================================================
