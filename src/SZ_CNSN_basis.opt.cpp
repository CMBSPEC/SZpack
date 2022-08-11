//==================================================================================================
//
// This program allows computing the thermal SZ and kinematic effect using the improved basis of 
// Chluba, Nagai, Sazonov & Nelson, 2012. Everything is similar to the routines provided in 
// SZ_CNSN_basis.h, but here the computations are carried out in the CMB frame only. Since only terms
// up to O(betac^2) are included, this leads to some tiny differences of the results. However, for 
// the temperature-velocity moment method explained by Chluba, Switzer, Nelson & Nagai, 2012, this
// choice is beneficial. 
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
#include "Relativistic_MB.h"
#include "nPl_derivatives.h"
#include "SZ_asymptotic.h"
#include "SZ_CNSN_basis.opt.h"

//==================================================================================================
//
// namespaces
//
//==================================================================================================
using namespace std;

//==================================================================================================
bool basis_loaded_CMB_opt=0;      

//==================================================================================================
// these are the reference temperatures that were used to compute basis;
//==================================================================================================
const int nregions=14;
const double Tref_basis[]={9.4/const_me, 11.0/const_me, 14.0/const_me, 18.5/const_me, 
                           22.0/const_me, 25.0/const_me, 30.0/const_me, 35.0/const_me, 
                           40.0/const_me, 45.0/const_me, 50.0/const_me, 55.0/const_me, 
                           60.0/const_me, 80.0/const_me}; 

#include "./src/database.opt/SZ_basis.Te_9.4keV.CMB.h" 
#include "./src/database.opt/SZ_basis.Te_11keV.CMB.h" 
#include "./src/database.opt/SZ_basis.Te_14keV.CMB.h"  
#include "./src/database.opt/SZ_basis.Te_18.5keV.CMB.h"
#include "./src/database.opt/SZ_basis.Te_22keV.CMB.h" 
#include "./src/database.opt/SZ_basis.Te_25keV.CMB.h" 
#include "./src/database.opt/SZ_basis.Te_30keV.CMB.h"
#include "./src/database.opt/SZ_basis.Te_35keV.CMB.h" 
#include "./src/database.opt/SZ_basis.Te_40keV.CMB.h" 
#include "./src/database.opt/SZ_basis.Te_45keV.CMB.h" 
#include "./src/database.opt/SZ_basis.Te_50keV.CMB.h" 
#include "./src/database.opt/SZ_basis.Te_55keV.CMB.h"  
#include "./src/database.opt/SZ_basis.Te_60keV.CMB.h"  
#include "./src/database.opt/SZ_basis.Te_80keV.CMB.h"  

//==================================================================================================
// mapping of basis functions
//==================================================================================================
int Get_pivot_index(string reg)
{
    if(reg=="9.4keV")  return 0;
    if(reg=="11keV")   return 1;
    if(reg=="14keV")   return 2;
    if(reg=="18.5keV") return 3;
    if(reg=="22keV")   return 4;
    if(reg=="25keV")   return 5;
    if(reg=="30keV")   return 6;
    if(reg=="35keV")   return 7;
    if(reg=="40keV")   return 8;
    if(reg=="45keV")   return 9;
    if(reg=="50keV")   return 10;
    if(reg=="55keV")   return 11;
    if(reg=="60keV")   return 12;
    if(reg=="80keV")   return 13;

    cerr << " This pivot temperature is not loaded. " << reg << endl; 
    exit(0);
    return -10;
}

const double (*Get_Data_thSZ(int i))[8]
{
    if(i==0) return Data_thSZ_0;
    else if(i==1) return Data_thSZ_1;
    else if(i==2) return Data_thSZ_2;
    else if(i==3) return Data_thSZ_3;
    else if(i==4) return Data_thSZ_4;
    else if(i==5) return Data_thSZ_5;
    else if(i==6) return Data_thSZ_6;
    else if(i==7) return Data_thSZ_7;
    else if(i==8) return Data_thSZ_8;
    else if(i==9) return Data_thSZ_9;
    else if(i==10) return Data_thSZ_10;
    else if(i==11) return Data_thSZ_11;
    else if(i==12) return Data_thSZ_12;
    
    return Data_thSZ_13;
}

const double (*Get_Data_kSZ(int i))[8]
{
    if(i==0) return Data_kSZ_0;
    else if(i==1) return Data_kSZ_1;
    else if(i==2) return Data_kSZ_2;
    else if(i==3) return Data_kSZ_3;
    else if(i==4) return Data_kSZ_4;
    else if(i==5) return Data_kSZ_5;
    else if(i==6) return Data_kSZ_6;
    else if(i==7) return Data_kSZ_7;
    else if(i==8) return Data_kSZ_8;
    else if(i==9) return Data_kSZ_9;
    else if(i==10) return Data_kSZ_10;
    else if(i==11) return Data_kSZ_11;
    else if(i==12) return Data_kSZ_12;
    
    return Data_kSZ_13;
}

const double (*Get_Data_qkSZ(int i))[8]
{
    if(i==0) return Data_qkSZ_0;
    else if(i==1) return Data_qkSZ_1;
    else if(i==2) return Data_qkSZ_2;
    else if(i==3) return Data_qkSZ_3;
    else if(i==4) return Data_qkSZ_4;
    else if(i==5) return Data_qkSZ_5;
    else if(i==6) return Data_qkSZ_6;
    else if(i==7) return Data_qkSZ_7;
    else if(i==8) return Data_qkSZ_8;
    else if(i==9) return Data_qkSZ_9;
    else if(i==10) return Data_qkSZ_10;
    else if(i==11) return Data_qkSZ_11;
    else if(i==12) return Data_qkSZ_12;
    
    return Data_qkSZ_13;
}

const double (*Get_Data_mkSZ(int i))[8]
{
    if(i==0) return Data_mkSZ_0;
    else if(i==1) return Data_mkSZ_1;
    else if(i==2) return Data_mkSZ_2;
    else if(i==3) return Data_mkSZ_3;
    else if(i==4) return Data_mkSZ_4;
    else if(i==5) return Data_mkSZ_5;
    else if(i==6) return Data_mkSZ_6;
    else if(i==7) return Data_mkSZ_7;
    else if(i==8) return Data_mkSZ_8;
    else if(i==9) return Data_mkSZ_9;
    else if(i==10) return Data_mkSZ_10;
    else if(i==11) return Data_mkSZ_11;
    else if(i==12) return Data_mkSZ_12;
    
    return Data_mkSZ_13;
}

//==================================================================================================
//
// optimal settings for kmax and temperature pivots according to Chluba et al. 2012
//
//==================================================================================================
struct Precision_settings
{
    int kmax;
    vector<double> Te_regions;
    vector<string> Te_pivots;
};

vector<Precision_settings> accurracy_I(3);    // absolute precision DIs
vector<Precision_settings> accurracy_II(3);   // absolute precision 0.1 DIs
vector<Precision_settings> accurracy_III(3);  // absolute precision 0.01 DIs
vector<Precision_settings> accurracy_IV(1);   // absolute precision 0.001 DIs

void init_precisions()
{
    //==============================================================================================
    accurracy_I[0].kmax=2;
    accurracy_I[0].Te_regions.resize(4, -10);
    accurracy_I[0].Te_pivots. resize(4, "");
    accurracy_I[0].Te_regions[0]=9.1; // asymptotic expansion part
    accurracy_I[0].Te_regions[1]=21.6;
    accurracy_I[0].Te_regions[2]=42.0;
    accurracy_I[0].Te_regions[3]=75.0;
    accurracy_I[0].Te_pivots[1]="14keV";
    accurracy_I[0].Te_pivots[2]="30keV";
    accurracy_I[0].Te_pivots[3]="55keV";
    
    accurracy_I[1].kmax=3;
    accurracy_I[1].Te_regions.resize(3, -10);
    accurracy_I[1].Te_pivots. resize(3, "");
    accurracy_I[1].Te_regions[0]=12.5; // asymptotic expansion part
    accurracy_I[1].Te_regions[1]=48.5; 
    accurracy_I[1].Te_regions[2]=90.0;
    accurracy_I[1].Te_pivots[1]="25keV";
    accurracy_I[1].Te_pivots[2]="80keV";

    accurracy_I[2].kmax=4;
    accurracy_I[2].Te_regions.resize(2, -10);
    accurracy_I[2].Te_pivots. resize(2, "");
    accurracy_I[2].Te_regions[0]=14.3; // asymptotic expansion part
    accurracy_I[2].Te_regions[1]=75.0; 
    accurracy_I[2].Te_pivots[1]="35keV";
    //==============================================================================================
    
    
    //==============================================================================================
    accurracy_II[0].kmax=3;
    accurracy_II[0].Te_regions.resize(4, -10);
    accurracy_II[0].Te_pivots. resize(4, "");
    accurracy_II[0].Te_regions[0]=6.8; // asymptotic expansion part
    accurracy_II[0].Te_regions[1]=17.25; 
    accurracy_II[0].Te_regions[2]=36.4;  
    accurracy_II[0].Te_regions[3]=68.0;
    accurracy_II[0].Te_pivots[1]="11keV";
    accurracy_II[0].Te_pivots[2]="25keV";
    accurracy_II[0].Te_pivots[3]="50keV";
    
    accurracy_II[1].kmax=4;
    accurracy_II[1].Te_regions.resize(3, -10);
    accurracy_II[1].Te_pivots. resize(3, "");
    accurracy_II[1].Te_regions[0]=9.3; // asymptotic expansion part
    accurracy_II[1].Te_regions[1]=33.1; 
    accurracy_II[1].Te_regions[2]=80.0;
    accurracy_II[1].Te_pivots[1]="18.5keV";
    accurracy_II[1].Te_pivots[2]="55keV";

    accurracy_II[2].kmax=5;
    accurracy_II[2].Te_regions.resize(3, -10);
    accurracy_II[2].Te_pivots. resize(3, "");
    accurracy_II[2].Te_regions[0]=10.76; // asymptotic expansion part
    accurracy_II[2].Te_regions[1]=48.5; 
    accurracy_II[2].Te_regions[2]=90.0;
    accurracy_II[2].Te_pivots[1]="25keV";
    accurracy_II[2].Te_pivots[2]="80keV";
    //==============================================================================================

    
    //==============================================================================================
    accurracy_III[0].kmax=4;
    accurracy_III[0].Te_regions.resize(4, -10);
    accurracy_III[0].Te_pivots. resize(4, "");
    accurracy_III[0].Te_regions[0]=5.76; // asymptotic expansion part
    accurracy_III[0].Te_regions[1]=14.7;
    accurracy_III[0].Te_regions[2]=32.0;
    accurracy_III[0].Te_regions[3]=61.0;
    accurracy_III[0].Te_pivots[1]="9.4keV";
    accurracy_III[0].Te_pivots[2]="22keV";
    accurracy_III[0].Te_pivots[3]="45keV";
    
    accurracy_III[1].kmax=5;
    accurracy_III[1].Te_regions.resize(3, -10);
    accurracy_III[1].Te_pivots. resize(3, "");
    accurracy_III[1].Te_regions[0]=7.3; // asymptotic expansion part
    accurracy_III[1].Te_regions[1]=23.86;
    accurracy_III[1].Te_regions[2]=62.0;
    accurracy_III[1].Te_pivots[1]="14keV";
    accurracy_III[1].Te_pivots[2]="40keV";
    
    accurracy_III[2].kmax=6;
    accurracy_III[2].Te_regions.resize(3, -10);
    accurracy_III[2].Te_pivots. resize(3, "");
    accurracy_III[2].Te_regions[0]=8.55; // asymptotic expansion part
    accurracy_III[2].Te_regions[1]=33.3; 
    accurracy_III[2].Te_regions[2]=80.0;
    accurracy_III[2].Te_pivots[1]="18.5keV";
    accurracy_III[2].Te_pivots[2]="60keV";
    //==============================================================================================

    //==============================================================================================
    accurracy_IV[0].kmax=6;
    accurracy_IV[0].Te_regions.resize(4, -10);
    accurracy_IV[0].Te_pivots. resize(4, "");
    accurracy_IV[0].Te_regions[0]=5.76; // asymptotic expansion part
    accurracy_IV[0].Te_regions[1]=14.7;
    accurracy_IV[0].Te_regions[2]=32.0;
    accurracy_IV[0].Te_regions[3]=63.5;
    accurracy_IV[0].Te_pivots[1]="9.4keV";
    accurracy_IV[0].Te_pivots[2]="22keV";
    accurracy_IV[0].Te_pivots[3]="45keV";
    //==============================================================================================

    return;
}


//==================================================================================================
//
// SZ effect using interpolation of derivative terms
//
//==================================================================================================
vector<vector<int> > spline_mem_indices_th_CMB_opt(nregions); 
vector<vector<int> > spline_mem_indices_d_CMB_opt(nregions);  
vector<vector<int> > spline_mem_indices_q_CMB_opt(nregions);  
vector<vector<int> > spline_mem_indices_m_CMB_opt(nregions);  

//==================================================================================================
void setup_expansion_splines_CMB(vector<int> &spline_mem_indices, const double D[][8])
{
    int cols=8;
    int np=1000;

    spline_mem_indices.resize(cols, -1);
    vector<double> xarr(np);
    vector<double> yarr(np);

    for(int i=0; i<np; i++) xarr[i] = D[i][0];
    
    for(int c=1; c<cols; c++)
    {
        for(int i=0; i<np; i++) yarr[i] = D[i][c];
        
        spline_mem_indices[c-1]=calc_spline_coeffies_JC(np, &xarr[0], &yarr[0], 
                                                        " setup_expansion_splines_CMB ");
    }
    
    return;
}

//==================================================================================================
// generic function for interpolation (x^3 Dn was used)
//==================================================================================================
double DI_SZ_splines_CMB(double x, double The, double The_ref, vector<int> &spline_mem, int kmax)
{
    double lx=log(x);
    double r=calc_spline_JC(lx, spline_mem[0]);

    for(int i=1; i<=kmax; i++) 
        r+=pow(The/The_ref-1.0, i)*calc_spline_JC(lx, spline_mem[i])/factorial(i);
    
    r*=norm_df_RM_dTheta(The)/norm_df_RM_dTheta(The_ref);

    return r;     
}

//==================================================================================================
//
// access basis functions and Normalization, etc. These functions are needed by the SZ moment 
// method but are otherwise not too useful for general purpose applications.
//
//==================================================================================================
double N_func_CNSN(double The){ return norm_df_RM_dTheta(The); }

void load_basis_CMB_opt()
{
    if(!basis_loaded_CMB_opt)
    {
        cout << " Loading CNSN 2012 basis functions (CMB frame/opt). This will only be done once. " 
             << endl;
        
        for(int r=0; r<nregions; r++)
        {
            setup_expansion_splines_CMB(spline_mem_indices_th_CMB_opt[r], Get_Data_thSZ(r));
            setup_expansion_splines_CMB(spline_mem_indices_d_CMB_opt [r], Get_Data_kSZ (r));
            setup_expansion_splines_CMB(spline_mem_indices_q_CMB_opt [r], Get_Data_qkSZ(r));
            setup_expansion_splines_CMB(spline_mem_indices_m_CMB_opt [r], Get_Data_mkSZ(r));
        }
        
        init_precisions();
        
        basis_loaded_CMB_opt=1;
    }
    
    return;
}

//==================================================================================================
//
// access basis functions
//
//==================================================================================================
void compute_XX_CNSN(double x, int r, vector<double> &XX, vector<int> &spline_mem)
{
    load_basis_CMB_opt();
    
    int Te_order=XX.size();
    double Tr=Tref_basis[r];
    double fac=1.0/N_func_CNSN(Tr);
    double lx=log(x), x3=pow(x, 3);
    
    for(int i=0; i<Te_order; i++) 
        XX[i]=fac/pow(Tr, i)*calc_spline_JC(lx, spline_mem[i])/factorial(i)/x3;
    
    return;
}

void compute_Y_CNSN(double x, int r, vector<double> &Y)
{
    compute_XX_CNSN(x, r, Y, spline_mem_indices_th_CMB_opt[r]); 
    return;
}

void compute_M_CNSN_CMB(double x, int r, vector<double> &Y)
{
    compute_XX_CNSN(x, r, Y, spline_mem_indices_m_CMB_opt[r]); 
    return;
}

void compute_D_CNSN_CMB(double x, int r, vector<double> &Y)
{
    compute_XX_CNSN(x, r, Y, spline_mem_indices_d_CMB_opt[r]); 
    return;
}

void compute_Q_CNSN_CMB(double x, int r, vector<double> &Y)
{
    compute_XX_CNSN(x, r, Y, spline_mem_indices_q_CMB_opt[r]); 
    return;
}

//==================================================================================================
// Get reference to required accuracy setting
//==================================================================================================
const Precision_settings& Get_accuracy_setting(int kmax, int accuracy_level)
{
    load_basis_CMB_opt();
    
    if(accuracy_level==0)
    { 
        if(kmax<2 || kmax>4) 
        {
            cerr << " compute_SZ_distortion_CNSN_basis_opt :: kmax = 2 - 4 only " 
                 << " (was " << kmax << ")" << endl; 
            exit(0);
        }
        else if(kmax==2) return accurracy_I[0];
        else if(kmax==3) return accurracy_I[1];
        else if(kmax==4) return accurracy_I[2];
    }
    else if(accuracy_level==1) 
    {
        if(kmax<3 || kmax>5) 
        {
            cerr << " compute_SZ_distortion_CNSN_basis_opt :: kmax = 3 - 5 only" 
                 << " (was " << kmax << ")" << endl;  
            exit(0);
        }
        else if(kmax==3) return accurracy_II[0];
        else if(kmax==4) return accurracy_II[1];
        else if(kmax==5) return accurracy_II[2];
    }
    else if(accuracy_level==2) 
    {
        if(kmax<4 || kmax>6) 
        {
            cerr << " compute_SZ_distortion_CNSN_basis_opt :: kmax = 4 - 6 only " 
                 << " (was " << kmax << ")" << endl;  
            exit(0);
        }
        else if(kmax==4) return accurracy_III[0];
        else if(kmax==5) return accurracy_III[1];
        else if(kmax==6) return accurracy_III[2];
    }
    else if(accuracy_level==3) 
    {
        if(kmax!=6) 
        {
            cerr << " compute_SZ_distortion_CNSN_basis_opt :: kmax = 6 only " 
                 << " (was " << kmax << ")" << endl;  
            exit(0);
        }
        else if(kmax==6) return accurracy_IV[0];
    }
    else
    { 
        cerr << " compute_SZ_distortion_CNSN_basis_opt :: accuracy level not available." 
             << " Should be 0-3 " << endl;
        exit(0); 
    }
    
    return accurracy_I[0];
}
    
//==================================================================================================
int determine_Te_region(Precision_settings &acc, double Te_max)
{
    int region=0, nmax=acc.Te_regions.size();
    
    if(Te_max>acc.Te_regions[nmax-1])
    {   
        cerr << " determine_Te_region:: temperature too high! Te_max was " << Te_max << endl; 
        exit(0);
    }
    else if(Te_max<acc.Te_regions[0]) return 1;

    do region++;
    while(Te_max>acc.Te_regions[region]);     

    return region+1;
}

//==================================================================================================
// computes the optimal value for kmax given the accuracy goal and required maximal temperature
//==================================================================================================
void determine_optimal_kmax(int accuracy_level, double Te_max, int &kmax, int &iregmax)
{
    load_basis_CMB_opt();

    if(accuracy_level==0)
    {
        int nreg0=determine_Te_region(accurracy_I[0], Te_max);
        int nreg1=determine_Te_region(accurracy_I[1], Te_max);
        int nreg2=determine_Te_region(accurracy_I[2], Te_max);
        
        if(nreg0*accurracy_I[0].kmax<nreg1*accurracy_I[1].kmax)
        { iregmax=nreg0; kmax=accurracy_I[0].kmax; }
        else{ iregmax=nreg1; kmax=accurracy_I[1].kmax; }
        if(iregmax>nreg2*accurracy_I[2].kmax)
        { iregmax=nreg2; kmax=accurracy_I[2].kmax; }
    }
    else if(accuracy_level==1)
    {
        int nreg0=determine_Te_region(accurracy_II[0], Te_max);
        int nreg1=determine_Te_region(accurracy_II[1], Te_max);
        int nreg2=determine_Te_region(accurracy_II[2], Te_max);
        
        if(nreg0*accurracy_II[0].kmax<nreg1*accurracy_II[1].kmax)
        { iregmax=nreg0; kmax=accurracy_II[0].kmax; }
        else{ iregmax=nreg1; kmax=accurracy_II[1].kmax; }
        if(iregmax>nreg2*accurracy_II[2].kmax)
        { iregmax=nreg2; kmax=accurracy_II[2].kmax; }
    }
    else if(accuracy_level==2)
    {
        int nreg0=determine_Te_region(accurracy_III[0], Te_max);
        int nreg1=determine_Te_region(accurracy_III[1], Te_max);
        int nreg2=determine_Te_region(accurracy_III[2], Te_max);
        
        if(nreg0*accurracy_III[0].kmax<nreg1*accurracy_III[1].kmax)
        { iregmax=nreg0; kmax=accurracy_III[0].kmax; }
        else{ iregmax=nreg1; kmax=accurracy_III[1].kmax; }
        if(iregmax>nreg2*accurracy_III[2].kmax)
        { iregmax=nreg2; kmax=accurracy_III[2].kmax; }
    }
    else if(accuracy_level==3)
    {
        int nreg0=determine_Te_region(accurracy_IV[0], Te_max);
        iregmax=nreg0; kmax=accurracy_IV[0].kmax;
    }
    else
    { 
        cerr << " determine_optimal_kmax :: accuracy level not available." 
             << " Should be 0-3 " << endl;
        exit(0); 
    }
    
    return;
}
    
//==================================================================================================
vector<double> Get_temperature_regions(int kmax, int accuracy_level)
{
    const Precision_settings *acc_p=&Get_accuracy_setting(kmax, accuracy_level);
    
    return acc_p->Te_regions;
}

vector<double> Get_temperature_pivots (int kmax, int accuracy_level)
{
    const Precision_settings *acc_p=&Get_accuracy_setting(kmax, accuracy_level);
    
    vector<double> Top(1, 0);
    
    for(int r=1; r<(int)acc_p->Te_pivots.size(); r++)
        Top.push_back(Tref_basis[Get_pivot_index(acc_p->Te_pivots[r])]);
    
    return Top;
}    

vector<int> Get_region_indices (int kmax, int accuracy_level)
{
    const Precision_settings *acc_p=&Get_accuracy_setting(kmax, accuracy_level);
    
    vector<int> dum(1, -10);

    for(int r=1; r<(int)acc_p->Te_pivots.size(); r++)
        dum.push_back(Get_pivot_index(acc_p->Te_pivots[r]));
    
    return dum;
}

//==================================================================================================
//
// compute Dn using improved expansion in CMB rest frame
//
// mode == "monopole"      --> only monopole part without second order kinematic corr
// mode == "dipole"        --> only dipolar part     (first order kinematic correction)
// mode == "quadrupole"    --> only quadrupolar part (second order kinematic correction)
// mode == "monopole_corr" --> only second order kinematic correction to monopole part
// mode == "all"           --> all terms added
// mode == "kin"           --> only kinematic terms
//
//==================================================================================================
double compute_SZ_distortion_CNSN_basis_opt(double x, 
                                            double The, double betac, double muc, 
                                            int kmax, int betac_order,
                                            string mode, int accuracy_level)
{
    //==============================================================================================
    // select precision of approximation
    //==============================================================================================
    const Precision_settings *acc_p=&Get_accuracy_setting(kmax, accuracy_level);
    
    //==============================================================================================
    // simplest diversion
    //==============================================================================================
    int nregs=acc_p->Te_regions.size();
    int piv, region=0;
    
    if(The<acc_p->Te_regions[0]/const_me) 
    {
        return compute_SZ_distortion_asymptotic_CMB(x, The, betac, muc, 
                                                    kmax-1, betac_order,
                                                    mode);
    }
    else if(The>acc_p->Te_regions[nregs-1]/const_me) 
    {
        cerr << " compute_SZ_distortion_CNSN_basis_opt :: Temperature too high! Should be < "  
             << acc_p->Te_regions[nregs-1] << " keV" << endl;
        exit(0); 
    }
    //==============================================================================================
    // selection of region and pivot
    //==============================================================================================
    else 
    {
        do region++;
        while(The>acc_p->Te_regions[region]/const_me);
        
        piv=Get_pivot_index(acc_p->Te_pivots[region]);
    }

    double r=0.0, rb, The_ref=Tref_basis[piv];
    double x3=pow(x, 3);
    
    //==============================================================================================
    // computation of SZ signal
    //==============================================================================================
    if(betac_order>0)
    {
        double G=x*exp(-x)/pow(one_minus_exp_mx(x), 2);
        
        if(mode=="dipole" || mode=="all" || mode=="kin") 
        {
            rb=DI_SZ_splines_CMB(x, The, The_ref, spline_mem_indices_d_CMB_opt[piv], kmax);
            rb+=x3*G;
            r+=betac*muc*rb;
        }
        
        if(betac_order>1)
        {
            double Q=x*x*exp(-x)*(exp(-x)+1.0)/pow(one_minus_exp_mx(x), 3);
            
            if(mode=="monopole_corr" || mode=="all" || mode=="kin") 
            {
                rb=DI_SZ_splines_CMB(x, The, The_ref, spline_mem_indices_m_CMB_opt[piv], kmax);
                rb+=x3*(Q-3.0*G)/3.0;
                r+=betac*betac*rb;
            }  
            
            if(mode=="quadrupole" || mode=="all" || mode=="kin") 
            {
                rb=DI_SZ_splines_CMB(x, The, The_ref, spline_mem_indices_q_CMB_opt[piv], kmax);
                rb+=x3*11.0/30.0*Q;
                r+=betac*betac*(1.5*muc*muc-0.5)*rb;
            }  
        }
    }
    
    if(mode=="monopole" || mode=="all") 
    {
        r+=DI_SZ_splines_CMB(x, The, The_ref, spline_mem_indices_th_CMB_opt[piv], kmax);
    }
    
    return r/x3;     
}

//==================================================================================================
//==================================================================================================
