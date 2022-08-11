//==================================================================================================
//
// This program allows computing the thermal and kinematic SZ effect using the basis functions of 
// Chluba, Nagai, Sazonov & Nelson, 2012. Up to 20th order temperature corrections are allowed. 
// Kinematic corrections can be included up to second order in the clusters peculiar velocity. 
// For 0.01 < x < 30, 2keV < Te < 75keV and beta < 0.01 the obtained SZ signal should be accurate 
// at the level of 0.001% when 20 temperature order are included. This precision is acheived using 
// three temperature pivots, The=0.01, 0.03, and 0.1. For Te < 2keV the asymptotic expansions 
// should be used.
//
//==================================================================================================
//
// Author: Jens Chluba  (CITA, University of Toronto)
//
// first implementation: April 2012
// last modification   : Sept  2013
//
//==================================================================================================
// 12th Sept 2013: fixed bug because of betac^2 P_2(muc) == beta_para^2 - beta_perp^2/2
//  4st  Aug: added derivatives of basis functions in the CMB rest frame
// 22th July: added low and high temperature expansions. Now 2keV < Te < 75keV is covered
// 21th July: added headers with required data; this avoids time taken for reading the files
// 10th July: added basis functions in the CMB rest frame
//  8th July: changed definition of S^kin; temperature-independent terms are fully canceled

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
#include "SZ_CNSN_basis.h"

//==================================================================================================
//
// namespaces
//
//==================================================================================================
using namespace std;

//==================================================================================================
// these are the reference temperatures that were used to compute basis;
//==================================================================================================
const double Tref_basis_low=0.01; 
const double Tref_basis=0.03; 
const double Tref_basis_high=0.1; 

bool basis_loaded=0, basis_loaded_CMB=0;

#include "./src/database/SZ_basis.The_0.01.h"
#include "./src/database/SZ_basis.The_0.03.h"
#include "./src/database/SZ_basis.The_0.1.h" 
// CMB frame data
#include "./src/database/SZ_basis.The_0.01.CMB.h"
#include "./src/database/SZ_basis.The_0.03.CMB.h"
#include "./src/database/SZ_basis.The_0.1.CMB.h" 

//==================================================================================================
//
// compute SZ effect using interpolation of derivative terms
//
//==================================================================================================
vector<int> spline_mem_indices_th_low, spline_mem_indices_th_low_CMB; 
vector<int> spline_mem_indices_d_low, spline_mem_indices_d_low_CMB;  
vector<int> spline_mem_indices_q_low, spline_mem_indices_q_low_CMB;  
vector<int> spline_mem_indices_m_low, spline_mem_indices_m_low_CMB;  

vector<int> spline_mem_indices_th, spline_mem_indices_th_CMB;
vector<int> spline_mem_indices_d, spline_mem_indices_d_CMB;
vector<int> spline_mem_indices_q, spline_mem_indices_q_CMB;
vector<int> spline_mem_indices_m, spline_mem_indices_m_CMB;

vector<int> spline_mem_indices_th_high, spline_mem_indices_th_high_CMB; 
vector<int> spline_mem_indices_d_high, spline_mem_indices_d_high_CMB;  
vector<int> spline_mem_indices_q_high, spline_mem_indices_q_high_CMB;  
vector<int> spline_mem_indices_m_high, spline_mem_indices_m_high_CMB;  

//==================================================================================================
void setup_expansion_splines(vector<int> &spline_mem_indices, const double D[][22])
{
    int cols=22;
    int np=1000;

    spline_mem_indices.resize(cols, -1);
    vector<double> xarr(np);
    vector<double> yarr(np);

    for(int i=0; i<np; i++) xarr[i] = D[i][0];
    
    for(int c=1; c<cols; c++)
    {
        for(int i=0; i<np; i++) yarr[i] = D[i][c];
        
        spline_mem_indices[c-1]=calc_spline_coeffies_JC(np, &xarr[0], &yarr[0], 
                                                        " setup_expansion_splines ");
    }
    
    return;
}

//==================================================================================================
// generic function for interpolation (x^3 Dn was used)
//==================================================================================================
double DI_SZ_splines(double x, double The, double The_ref, vector<int> &spline_mem, int Torder)
{
    double lx=log(x);
    double r=calc_spline_JC(lx, spline_mem[0]);

    for(int i=1; i<Torder; i++) 
        r+=pow(The/The_ref-1.0, i)*calc_spline_JC(lx, spline_mem[i])/factorial(i);
    
    r*=norm_df_RM_dTheta(The)/norm_df_RM_dTheta(The_ref);

    return r;     
}

//==================================================================================================
void load_basis()
{
    if(!basis_loaded)
    {
        cout << " Loading CNSN 2012 basis functions. This will only be done once. " << endl;
        
        setup_expansion_splines(spline_mem_indices_th_low, Data_thSZ_low);
        setup_expansion_splines(spline_mem_indices_d_low, Data_kSZ_low);
        setup_expansion_splines(spline_mem_indices_q_low, Data_qkSZ_low);
        setup_expansion_splines(spline_mem_indices_m_low, Data_mkSZ_low);
        
        setup_expansion_splines(spline_mem_indices_th, Data_thSZ);
        setup_expansion_splines(spline_mem_indices_d, Data_kSZ);
        setup_expansion_splines(spline_mem_indices_q, Data_qkSZ);
        setup_expansion_splines(spline_mem_indices_m, Data_mkSZ);
        
        setup_expansion_splines(spline_mem_indices_th_high, Data_thSZ_high);
        setup_expansion_splines(spline_mem_indices_d_high, Data_kSZ_high);
        setup_expansion_splines(spline_mem_indices_q_high, Data_qkSZ_high);
        setup_expansion_splines(spline_mem_indices_m_high, Data_mkSZ_high);
        
        basis_loaded=1;
    }
    
    return;
}

//==================================================================================================
void load_basis_CMB()
{
    if(!basis_loaded_CMB)
    {
        cout << " Loading CNSN 2012 basis functions (CMB frame). This will only be done once. " 
             << endl;
        
        setup_expansion_splines(spline_mem_indices_th_low_CMB, Data_thSZ_I);
        setup_expansion_splines(spline_mem_indices_d_low_CMB, Data_kSZ_I);
        setup_expansion_splines(spline_mem_indices_q_low_CMB, Data_qkSZ_I);
        setup_expansion_splines(spline_mem_indices_m_low_CMB, Data_mkSZ_I);
        
        setup_expansion_splines(spline_mem_indices_th_CMB, Data_thSZ_II);
        setup_expansion_splines(spline_mem_indices_d_CMB, Data_kSZ_II);
        setup_expansion_splines(spline_mem_indices_q_CMB, Data_qkSZ_II);
        setup_expansion_splines(spline_mem_indices_m_CMB, Data_mkSZ_II);
        
        setup_expansion_splines(spline_mem_indices_th_high_CMB, Data_thSZ_III);
        setup_expansion_splines(spline_mem_indices_d_high_CMB, Data_kSZ_III);
        setup_expansion_splines(spline_mem_indices_q_high_CMB, Data_qkSZ_III);
        setup_expansion_splines(spline_mem_indices_m_high_CMB, Data_mkSZ_III);
        
        basis_loaded_CMB=1;
    }
    
    return;
}

//==================================================================================================
//
// compute Dn using improved expansion in cluster frame (Chluba, Nagai, Sazonov, Nelson, 2012)
//
// mode == "monopole"      --> only scattering of monopole without second order kinematic corr
// mode == "dipole"        --> only scattering of dipole     (first order kinematic correction)
// mode == "quadrupole"    --> only scattering of quadrupole (second order kinematic correction)
// mode == "monopole_corr" --> only scattering of second order kinematic correction to monopole
// mode == "all"           --> all terms added
// mode == "kin"           --> only kinematic terms
//
//==================================================================================================
double compute_SZ_distortion_CNSN_basis(double x, 
                                        double The, double betac, double muc, 
                                        int Te_order, int betac_order,
                                        string mode)
{
    //==============================================================================================
    double r=0.0, rb, The_ref=Tref_basis;
    vector<int> *spline_th=&spline_mem_indices_th; 
    vector<int> *spline_d =&spline_mem_indices_d;  
    vector<int> *spline_q =&spline_mem_indices_q;  
    vector<int> *spline_m =&spline_mem_indices_m;  
    
    //==============================================================================================
    // load basis if needed
    //==============================================================================================
    load_basis();
    
    //==============================================================================================
    // check which temperature range applies
    //==============================================================================================
    if(The<=0.015)
    {
        The_ref=Tref_basis_low;
        spline_th=&spline_mem_indices_th_low; 
        spline_d =&spline_mem_indices_d_low;  
        spline_q =&spline_mem_indices_q_low;  
        spline_m =&spline_mem_indices_m_low;  
    }
    else if(The>=0.045)
    {
        The_ref=Tref_basis_high;
        spline_th=&spline_mem_indices_th_high; 
        spline_d =&spline_mem_indices_d_high;  
        spline_q =&spline_mem_indices_q_high;  
        spline_m =&spline_mem_indices_m_high;  
    }
    
    //==============================================================================================
    // computation of SZ signal
    //==============================================================================================
    Te_order++;
    
    double x3=pow(x, 3);
    if(betac_order>0)
    {
        if(mode=="dipole" || mode=="all" || mode=="kin") 
        {
            rb=DI_SZ_splines(x, The, The_ref, *spline_d, Te_order);
            rb+=x3*x*exp(-x)/pow(one_minus_exp_mx(x), 2);
            r+=betac*muc*rb;
        }

        if(betac_order>1)
        {
            if(mode=="monopole_corr" || mode=="all" || mode=="kin") 
            {
                rb=DI_SZ_splines(x, The, The_ref, *spline_m, Te_order);
                r+=betac*betac*rb;
            }  
            
            if(mode=="quadrupole" || mode=="all" || mode=="kin") 
            {
                rb=DI_SZ_splines(x, The, The_ref, *spline_q, Te_order);
                //rb-=x3*x*x*exp(-x)*(exp(-x)+1.0)/pow(one_minus_exp_mx(x), 3)/3.0; //JC08.07
                rb-=x3*0.3*x*x*exp(-x)*(exp(-x)+1.0)/pow(one_minus_exp_mx(x), 3);
                r+=betac*betac*(1.5*muc*muc-0.5)*rb;
            }  
        }
    }
    
    if(mode=="monopole" || mode=="all") 
    {
        r+=DI_SZ_splines(x, The, The_ref, *spline_th, Te_order);
    }
        
    return r/x3;     
}

//==================================================================================================
//
// compute Dn using expansion of CNSN2012 but in the CMB frame (added 10.07 by JC)
//
// mode == "monopole"      --> only scattering of monopole without second order kinematic corr
// mode == "dipole"        --> only scattering of dipole     (first order kinematic correction)
// mode == "quadrupole"    --> only scattering of quadrupole (second order kinematic correction)
// mode == "monopole_corr" --> only scattering of second order kinematic correction to monopole
// mode == "all"           --> all terms added
// mode == "kin"           --> only kinematic terms
//
// Note that since only terms up to O(betac^2) are included, the results obtained with this routine  
// are slightly different from those obtained by Lorentz transformation of the cluster-frame signal 
// into the CMB frame.
//
//==================================================================================================
double compute_SZ_distortion_CNSN_basis_CMB(double x, 
                                            double The, double betac, double muc, 
                                            int Te_order, int betac_order,
                                            string mode)
{
    //==============================================================================================
    double r=0.0, rb, The_ref=Tref_basis;
    vector<int> *spline_th=&spline_mem_indices_th_CMB; 
    vector<int> *spline_d =&spline_mem_indices_d_CMB;  
    vector<int> *spline_q =&spline_mem_indices_q_CMB;  
    vector<int> *spline_m =&spline_mem_indices_m_CMB;  
    
    //==============================================================================================
    // load basis if needed
    //==============================================================================================
    load_basis_CMB();
    
    //==============================================================================================
    // check which temperature range applies
    //==============================================================================================
    if(The<=0.015)
    {
        The_ref=Tref_basis_low;
        spline_th=&spline_mem_indices_th_low_CMB; 
        spline_d =&spline_mem_indices_d_low_CMB;  
        spline_q =&spline_mem_indices_q_low_CMB;  
        spline_m =&spline_mem_indices_m_low_CMB;  
    }
    else if(The>=0.045)
    {
        The_ref=Tref_basis_high;
        spline_th=&spline_mem_indices_th_high_CMB; 
        spline_d =&spline_mem_indices_d_high_CMB;  
        spline_q =&spline_mem_indices_q_high_CMB;  
        spline_m =&spline_mem_indices_m_high_CMB;  
    }
    
    //==============================================================================================
    // computation of SZ signal
    //==============================================================================================
    Te_order++;
    
    double x3=pow(x, 3);
    if(betac_order>0)
    {
        double G=x*exp(-x)/pow(one_minus_exp_mx(x), 2);

        if(mode=="dipole" || mode=="all" || mode=="kin") 
        {
            rb=DI_SZ_splines(x, The, The_ref, *spline_d, Te_order);
            rb+=x3*G;
            r+=betac*muc*rb;
        }
        
        if(betac_order>1)
        {
            double Q=x*x*exp(-x)*(1.0+exp(-x))/pow(one_minus_exp_mx(x), 3);
            
            if(mode=="monopole_corr" || mode=="all" || mode=="kin") 
            {
                rb=DI_SZ_splines(x, The, The_ref, *spline_m, Te_order);
                rb+=x3*(Q-3.0*G)/3.0;
                r+=betac*betac*rb;
            }  
            
            if(mode=="quadrupole" || mode=="all" || mode=="kin") 
            {
                rb=DI_SZ_splines(x, The, The_ref, *spline_q, Te_order);
                rb+=x3*11.0/30.0*Q;
                r+=betac*betac*(1.5*muc*muc-0.5)*rb;
            }  
        }
    }
    
    if(mode=="monopole" || mode=="all") 
    {
        r+=DI_SZ_splines(x, The, The_ref, *spline_th, Te_order);
    }
    
    return r/x3;         
}

//==================================================================================================
// analytic derivatives in betac_parallel and betac_perp
//==================================================================================================
double Dn_CNSN_CMB_for_The(double x, double The, double beta_para, double beta_perp,
                           double The_ref, 
                           vector<int> &spline_th, vector<int> &spline_d,  
                           vector<int> &spline_q, vector<int> &spline_m)
{ 
    double r=0.0, rb;
    int Te_order=21;
    double x3=pow(x, 3);
    double G=x*exp(-x)/pow(one_minus_exp_mx(x), 2);
    double Q=x*x*exp(-x)*(1.0+exp(-x))/pow(one_minus_exp_mx(x), 3);

    rb=DI_SZ_splines(x, The, The_ref, spline_d, Te_order);
    rb+=x3*G;
    r+=beta_para*rb;
        
    rb=DI_SZ_splines(x, The, The_ref, spline_m, Te_order);
    rb+=x3*(Q-3.0*G)/3.0;
    r+=(beta_para*beta_para+beta_perp*beta_perp)*rb;
            
    rb=DI_SZ_splines(x, The, The_ref, spline_q, Te_order);
    rb+=x3*11.0/30.0*Q;
    // 12th Sept 2013: betac^2 P_2(muc) == beta_para^2 - 0.5*beta_perp^2 !
    r+=(beta_para*beta_para-0.5*beta_perp*beta_perp)*rb;

    r+=DI_SZ_splines(x, The, The_ref, spline_th, Te_order);
    
    return r/x3;         
}

double Dn_dbeta_para_CNSN_CMB(double x, double The, double beta_para, double beta_perp,
                              double The_ref, 
                              vector<int> &spline_th, vector<int> &spline_d,  
                              vector<int> &spline_q, vector<int> &spline_m)
{ 
    double r=0.0, rb;
    int Te_order=21;
    double x3=pow(x, 3);
    double G=x*exp(-x)/pow(one_minus_exp_mx(x), 2);
    double Q=x*x*exp(-x)*(1.0+exp(-x))/pow(one_minus_exp_mx(x), 3);
    
    rb=DI_SZ_splines(x, The, The_ref, spline_d, Te_order);
    rb+=x3*G;
    r+=rb;
    
    rb=DI_SZ_splines(x, The, The_ref, spline_m, Te_order);
    rb+=x3*(Q-3.0*G)/3.0;
    r+=2.0*beta_para*rb;
    
    rb=DI_SZ_splines(x, The, The_ref, spline_q, Te_order);
    rb+=x3*11.0/30.0*Q;
    r+=2.0*beta_para*rb;
    
    return r/x3;         
}

double Dn_d2beta_para_CNSN_CMB(double x, double The, double beta_para, double beta_perp,
                               double The_ref, 
                               vector<int> &spline_th, vector<int> &spline_d,  
                               vector<int> &spline_q, vector<int> &spline_m)
{ 
    double r=0.0, rb;
    int Te_order=21;
    double x3=pow(x, 3);
    double G=x*exp(-x)/pow(one_minus_exp_mx(x), 2);
    double Q=x*x*exp(-x)*(1.0+exp(-x))/pow(one_minus_exp_mx(x), 3);
    
    rb=DI_SZ_splines(x, The, The_ref, spline_m, Te_order);
    rb+=x3*(Q-3.0*G)/3.0;
    r+=2.0*rb;
    
    rb=DI_SZ_splines(x, The, The_ref, spline_q, Te_order);
    rb+=x3*11.0/30.0*Q;
    r+=2.0*rb;
    
    return r/x3 / 2.0;         
}

double Dn_dbeta2_perp_CNSN_CMB(double x, double The, double beta_para, double beta_perp,
                               double The_ref, 
                               vector<int> &spline_th, vector<int> &spline_d,  
                               vector<int> &spline_q, vector<int> &spline_m)
{ 
    double r=0.0, rb;
    int Te_order=21;
    double x3=pow(x, 3);
    double G=x*exp(-x)/pow(one_minus_exp_mx(x), 2);
    double Q=x*x*exp(-x)*(1.0+exp(-x))/pow(one_minus_exp_mx(x), 3);
    
    rb=DI_SZ_splines(x, The, The_ref, spline_m, Te_order);
    rb+=x3*(Q-3.0*G)/3.0;
    r+= rb;
    
    rb=DI_SZ_splines(x, The, The_ref, spline_q, Te_order);
    rb+=x3*11.0/30.0*Q;
    // 12th Sept 2013: betac^2 P_2(muc) == beta_para^2 - 0.5*beta_perp^2 !
    r+=-0.5*rb;
    
    return r/x3;         
}

//==================================================================================================
void compute_all_Te_derivatives_upto_dThe_CNSN(double x, double The, 
                                               double beta_para, double beta_perp, double The_ref, 
                                               vector<int> &sp_th, vector<int> &sp_d,  
                                               vector<int> &sp_q, vector<int> &sp_m,
                                               int dThe, 
                                               double (*f)(double, double, double, double, double, 
                                                           vector<int> &, vector<int> &,  
                                                           vector<int> &, vector<int> &),
                                               vector<double> &dDn_dThe)
{
    double d0, dp, dm, dp2, dm2;
    double eps=0.01;
    
    dDn_dThe.resize(dThe+1);
    d0=f(x, The, beta_para, beta_perp, The_ref, sp_th, sp_d, sp_q, sp_m);
    dDn_dThe[0]=d0;
    
    if(dThe>0)
    {
        dp =f(x, The*(1.0+eps),   beta_para, beta_perp, The_ref, sp_th, sp_d, sp_q, sp_m);
        dm =f(x, The*(1.0-eps),   beta_para, beta_perp, The_ref, sp_th, sp_d, sp_q, sp_m);
        dp2=f(x, The*(1.0+2*eps), beta_para, beta_perp, The_ref, sp_th, sp_d, sp_q, sp_m);
        dm2=f(x, The*(1.0-2*eps), beta_para, beta_perp, The_ref, sp_th, sp_d, sp_q, sp_m);
        
        dDn_dThe[1]=(-dp2+8.0*dp-8.0*dm+dm2)/12.0/eps;
        
        if(dThe>1)
        {
            dDn_dThe[2]=(dp-2.0*d0+dm)/pow(eps, 2) / (2.0);
            
            if(dThe>2)
            {
                dDn_dThe[3]=(dp2-2.0*dp+2.0*dm-dm2)/2.0/pow(eps, 3) / (6.0);
                
                if(dThe>3)
                    dDn_dThe[4]=(dp2-4.0*dp+6.0*d0-4.0*dm+dm2)/pow(eps, 4) / (24.0);
            }
        }
    }
    
    return;
}

//==================================================================================================
// Derivatives (The^k d^k_dThe /k!) (betapara^m d^m_dbetapara /m!) (beta2perp^l d^l_beta2perp /l!) S
// in the CMB frame for a resting observer. Maximal orders in The and betac are used to compute 
// the derivatives.
//
// constraints: dThe<=4; dbeta_para<=2; dbeta2_perp<=1;
//==================================================================================================
void Dcompute_SZ_distortion_CNSN_basis_CMB(double x, 
                                           int dThe, int dbeta_para, int dbeta2_perp,
                                           double The, double betac, double muc,
                                           vector<double> &dDn_dThe)
{
    dDn_dThe.resize(dThe+1);
    for(int k=0; k<dThe+1; k++) dDn_dThe[k]=0.0;
    
    if(dThe>4) return;
    if(dbeta_para>2) return;
    if(dbeta2_perp>1) return;
    if(dbeta2_perp==1 && dbeta_para!=0) return;
    
    //==============================================================================================
    // select correct pivot point
    //==============================================================================================
    double The_ref=Tref_basis;
    vector<int> *spline_th=&spline_mem_indices_th_CMB; 
    vector<int> *spline_d =&spline_mem_indices_d_CMB;  
    vector<int> *spline_q =&spline_mem_indices_q_CMB;  
    vector<int> *spline_m =&spline_mem_indices_m_CMB;  
    
    //==============================================================================================
    // load basis if needed
    //==============================================================================================
    load_basis_CMB();
    
    //==============================================================================================
    // check which temperature range applies
    //==============================================================================================
    if(The<=0.015)
    {
        The_ref=Tref_basis_low;
        spline_th=&spline_mem_indices_th_low_CMB; 
        spline_d =&spline_mem_indices_d_low_CMB;  
        spline_q =&spline_mem_indices_q_low_CMB;  
        spline_m =&spline_mem_indices_m_low_CMB;  
    }
    else if(The>=0.045)
    {
        The_ref=Tref_basis_high;
        spline_th=&spline_mem_indices_th_high_CMB; 
        spline_d =&spline_mem_indices_d_high_CMB;  
        spline_q =&spline_mem_indices_q_high_CMB;  
        spline_m =&spline_mem_indices_m_high_CMB;  
    }

    double betac_perp=betac*sqrt(1.0-muc*muc);
    double betac_para=betac*muc;

    //==============================================================================================
    // analytic and numerical derivatives
    //==============================================================================================
    if(dbeta2_perp==1)
    {
        compute_all_Te_derivatives_upto_dThe_CNSN(x, The, betac_para, betac_perp, The_ref, 
                                                  *spline_th, *spline_d, *spline_q, *spline_m, dThe, 
                                                  Dn_dbeta2_perp_CNSN_CMB,
                                                  dDn_dThe);
    }
 
    else if(dbeta_para==1)
    {
        compute_all_Te_derivatives_upto_dThe_CNSN(x, The, betac_para, betac_perp, The_ref, 
                                                  *spline_th, *spline_d, *spline_q, *spline_m, dThe, 
                                                  Dn_dbeta_para_CNSN_CMB,
                                                  dDn_dThe);
    }

    else if(dbeta_para==2)
    {
        compute_all_Te_derivatives_upto_dThe_CNSN(x, The, betac_para, betac_perp, The_ref, 
                                                  *spline_th, *spline_d, *spline_q, *spline_m, dThe, 
                                                  Dn_d2beta_para_CNSN_CMB,
                                                  dDn_dThe);
    }
    
    //==============================================================================================
    // The derivatives only
    //==============================================================================================
    else
    {

        compute_all_Te_derivatives_upto_dThe_CNSN(x, The, betac_para, betac_perp, The_ref, 
                                                  *spline_th, *spline_d, *spline_q, *spline_m, dThe, 
                                                  Dn_CNSN_CMB_for_The,
                                                  dDn_dThe);
    }
    
    return;
}

//==================================================================================================
//==================================================================================================
