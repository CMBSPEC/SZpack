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
// Author: Jens Chluba & Elizabeth Lee
//
// first implementation: July 2012
// last modification   : February 2020
//
//==================================================================================================
// 02/2020: Added Class structure

#include "SZ_CNSN_basis.opt.h"

//==================================================================================================
bool basis_loaded_CMB_opt=0;      

//==================================================================================================
// these are the reference temperatures that were used to compute basis;
//==================================================================================================
const int nregions=14;
const vector<double> Tref_basis={9.4/const_me, 11.0/const_me, 14.0/const_me, 18.5/const_me, 
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
CNSNoptSplineMembers::CNSNoptSplineMembers(){
    loaded = false;
}

void CNSNoptSplineMembers::setup_expansion_splines(vector<int> &spline_mem_indices, const double D[][8])
{
    int cols=8;
    int np=1000;

    spline_mem_indices.resize(cols, -1);
    vector<double> xarr(np), yarr(np);

    for(int i=0; i<np; i++) xarr[i] = D[i][0];
    
    for(int c=1; c<cols; c++)
    {
        for(int i=0; i<np; i++) yarr[i] = D[i][c];
        
        spline_mem_indices[c-1]=calc_spline_coeffies_JC(np, &xarr[0], &yarr[0], 
                                                        " setup_expansion_splines ");
    }
}

void CNSNoptSplineMembers::loadSplines(double refThe, const double dY[][8],const double dD[][8],const double dQ[][8],const double dMcorr[][8]){
    if(!loaded){
        setup_expansion_splines(Y, dY);
        setup_expansion_splines(D, dD);
        setup_expansion_splines(Q, dQ);
        setup_expansion_splines(Mcorr, dMcorr);
        loaded = true;
        The_ref = refThe;
    }
}
//==================================================================================================
// mapping of basis functions
//==================================================================================================
namespace CNSNoptSplines{
    vector<CNSNoptSplineMembers> spline_vec(nregions);
    bool loadCNSNbasis_opt;

    //loading expansion splines
    void load_allSplines(){
        spline_vec[0].loadSplines(Tref_basis[0],Data_thSZ_0,Data_kSZ_0,Data_qkSZ_0,Data_mkSZ_0);
        spline_vec[1].loadSplines(Tref_basis[1],Data_thSZ_1,Data_kSZ_1,Data_qkSZ_1,Data_mkSZ_1);
        spline_vec[2].loadSplines(Tref_basis[2],Data_thSZ_2,Data_kSZ_2,Data_qkSZ_2,Data_mkSZ_2);
        spline_vec[3].loadSplines(Tref_basis[3],Data_thSZ_3,Data_kSZ_3,Data_qkSZ_3,Data_mkSZ_3);
        spline_vec[4].loadSplines(Tref_basis[4],Data_thSZ_4,Data_kSZ_4,Data_qkSZ_4,Data_mkSZ_4);
        spline_vec[5].loadSplines(Tref_basis[5],Data_thSZ_5,Data_kSZ_5,Data_qkSZ_5,Data_mkSZ_5);
        spline_vec[6].loadSplines(Tref_basis[6],Data_thSZ_6,Data_kSZ_6,Data_qkSZ_6,Data_mkSZ_6);
        spline_vec[7].loadSplines(Tref_basis[7],Data_thSZ_7,Data_kSZ_7,Data_qkSZ_7,Data_mkSZ_7);
        spline_vec[8].loadSplines(Tref_basis[8],Data_thSZ_8,Data_kSZ_8,Data_qkSZ_8,Data_mkSZ_8);
        spline_vec[9].loadSplines(Tref_basis[9],Data_thSZ_9,Data_kSZ_9,Data_qkSZ_9,Data_mkSZ_9);
        spline_vec[10].loadSplines(Tref_basis[10],Data_thSZ_10,Data_kSZ_10,Data_qkSZ_10,Data_mkSZ_10);
        spline_vec[11].loadSplines(Tref_basis[11],Data_thSZ_11,Data_kSZ_11,Data_qkSZ_11,Data_mkSZ_11);
        spline_vec[12].loadSplines(Tref_basis[12],Data_thSZ_12,Data_kSZ_12,Data_qkSZ_12,Data_mkSZ_12);
        spline_vec[13].loadSplines(Tref_basis[13],Data_thSZ_13,Data_kSZ_13,Data_qkSZ_13,Data_mkSZ_13);
    }
}

//==================================================================================================

void precision_settings::initialise_settings(int maxk, int len, const vector<double> Tregions, const vector<double> Tpivots){
    kmax = maxk;
    length = len;
    Te_regions = Tregions;
    Te_pivots = Tpivots;

    pivot_indices.resize(length-1);
    for (int i = 0; i<length-1; i++){
        pivot_indices[i] = get_pivot_index(i);
    }
}

int precision_settings::determine_Te_region(double Te_max){
    region = 0; //i've made region 0 indexed, idk, if that's an issue anywhere
    
    if(Te_max>Te_regions[length-1])
    {   
        exit_error("determine_Te_region:: temperature too high! Te_max was "+to_string(Te_max));
    }

    while(Te_max>Te_regions[region]) { region++; } 

    return region;
}

int precision_settings::get_pivot_index(int reg){
    // Check if pivot temperature exists as loaded temperature.
    for (int i = 0; i < nregions; i++){
        if (Tref_basis[i]==Te_pivots[reg]/const_me){
            return i;
        }
    }

    exit_error("This pivot temperature is not loaded. " +to_string(Te_pivots[region]) + " keV"); 
    return -10;
}


//==================================================================================================
// optimal settings for kmax and temperature pivots according to CSNN (Chluba et al. 2013)
//==================================================================================================
namespace CNSNopt_accuracies{
    vector<precision_settings> acc_I(3);    // absolute precision DIs
    vector<precision_settings> acc_II(3);   // absolute precision 0.1 DIs
    vector<precision_settings> acc_III(3);  // absolute precision 0.01 DIs
    vector<precision_settings> acc_IV(1);   // absolute precision 0.001 DIs

    void initialise_precisions(){ //below the first region, calculated through asymptotic regime
        acc_I[0].initialise_settings(2, 4,   {9.1,21.6,42.0,75.0},  {14.0,30.0, 55.0});
        acc_I[1].initialise_settings(3, 3,   {12.5,48.5,90.0},      {25.0,80.0});
        acc_I[2].initialise_settings(4, 2,   {14.3,75.0},           {35.0});

        acc_II[0].initialise_settings(3, 4,  {16.8,17.25,36.4,68.0},{11.0,25.0,50.0});
        acc_II[1].initialise_settings(4, 3,  {9.3,33.1,80.0},       {18.5,55.0});
        acc_II[2].initialise_settings(5, 3,  {10.76,48.5,90.0},     {25.0,80.0});

        acc_III[0].initialise_settings(4, 4, {5.76,14.7,32.0,61.0}, {9.4,22.0,45.0});
        acc_III[1].initialise_settings(5, 3, {7.3,23.86,62.0},      {14.0,40.0});
        acc_III[2].initialise_settings(6, 3, {8.55,33.3,80.0},      {18.5,60.0});

        acc_IV[0].initialise_settings(6, 4,  {5.76,14.7,32.0,63.5}, {9.4,22.0,45.0});
    }
}

//==================================================================================================
//
// SZ effect using interpolation of derivative terms
//
//==================================================================================================
IntegralCNSNopt::IntegralCNSNopt(double xfac_i, double x_i, double The_i, double betac_i, double muc_i, int kmax_i, int accuracy_level, int betac_order_i, bool CMB){
    xfac = xfac_i;
    x = xfac*x_i;
    The = The_i;
    betac = betac_i;
    muc = muc_i;
    CMBframe = CMB;
    kmax = kmax_i;
    betac_order = betac_order_i;

    x3=pow(x, 3);
    exp_mx=exp(-x);
    dex=one_minus_exp_mx(x, exp_mx);

    CNSNoptSplines::load_allSplines();

    getAccuracySettings(kmax, accuracy_level);

    region = accuracy.determine_Te_region(The*const_me);
    if (region != 0){
        piv = accuracy.pivot_indices[region-1];
        spline = CNSNoptSplines::spline_vec[piv];
    }
}

IntegralCNSNopt::IntegralCNSNopt(double x_i, double The_i, double betac_i, double muc_i, int kmax_i, int accuracy_level, int betac_order_i, bool CMB)
    : IntegralCNSNopt(1.0, x_i, The_i, betac_i, muc_i, kmax_i, accuracy_level, betac_order_i, CMB) {}

IntegralCNSNopt::IntegralCNSNopt()
    : IntegralCNSNopt(1.0, 0.0, 0.0, 0.0, 0.0, 6, 3, 2, true) {}

IntegralCNSNopt::IntegralCNSNopt(int k, Parameters fp, bool CMB)
    : IntegralCNSNopt(CMB ? fp.calc.xfacCMB : fp.calc.xfac, fp.xcmb[k], fp.calc.The, fp.betac, CMB ? fp.muc : fp.calc.mucc, 
                        fp.kmax, fp.accuracy_level, fp.beta_order, CMB) {}

void IntegralCNSNopt::Update_x(double x_i){
    x = xfac*x_i;
    x3=pow(x, 3);
    exp_mx=exp(-x);
    dex=one_minus_exp_mx(x, exp_mx);
}

IntegralCNSNopt::IntegralCNSNopt(int kmax, int accuracy_level){
    getAccuracySettings(kmax, accuracy_level);
} //For Accuracy settings only

IntegralCNSNopt::IntegralCNSNopt(double x_i, int region, int Te_order_i){
    x = x_i;
    Te_order = Te_order_i;
    x3 = pow(x,3);
    CNSNoptSplines::load_allSplines();

    spline = CNSNoptSplines::spline_vec[region];
    Compute_Y();
    Compute_D();
    Compute_Q();
    Compute_Mcorr();
} // for calculating the basis functions

// Get reference to required accuracy setting
void IntegralCNSNopt::getAccuracySettings(int kmax, int accuracy_level){
    CNSNopt_accuracies::initialise_precisions();
    
    if(accuracy_level==0){ 
        if(kmax<2 || kmax>4){
            exit_error(" compute_SZ_distortion_CNSN_basis_opt :: kmax = 2 - 4 only (was " + to_string(kmax) + ")");
        }
        accuracy = CNSNopt_accuracies::acc_I[kmax-2];
    }
    else if(accuracy_level==1){
        if(kmax<3 || kmax>5){
            exit_error("compute_SZ_distortion_CNSN_basis_opt :: kmax = 3 - 5 only (was " + to_string(kmax) + ")");
        }
        accuracy = CNSNopt_accuracies::acc_I[kmax-3];
    }
    else if(accuracy_level==2){
        if(kmax<4 || kmax>6){
            exit_error("compute_SZ_distortion_CNSN_basis_opt :: kmax = 4 - 6 only (was " + to_string(kmax) + ")");
        }
        accuracy = CNSNopt_accuracies::acc_I[kmax-4];
    }
    else if(accuracy_level==3) 
    {
        if(kmax!=6){
            exit_error("compute_SZ_distortion_CNSN_basis_opt :: kmax = 6 only (was " + to_string(kmax) + ")");
        }
        accuracy = CNSNopt_accuracies::acc_I[kmax-6];
    }
    else { 
        exit_error("compute_SZ_distortion_CNSN_basis_opt :: accuracy level not available. Should be 0-3 ");
    }
}

// generic function for interpolation (x^3 Dn was used). Outputs Dn
double IntegralCNSNopt::Dn_SZ_splines(double The_ref, vector<int> &spline_mem){
    double lx=log(x);
    double r=calc_spline_JC(lx, spline_mem[0]);

    for(int i=1; i<kmax; i++) 
        r+=pow(The/The_ref-1.0, i)*calc_spline_JC(lx, spline_mem[i])/factorial(i);
    
    r*=norm_df_RM_dTheta(The)/norm_df_RM_dTheta(The_ref)*The_ref/The;

    return r/x3;    
}

// access basis functions and Normalization, etc. These functions are needed by the SZ moment 
// function to derive the Y, D, Q and Mcorr vectors as in asymptotic
void IntegralCNSNopt::Compute_XX(double The_ref, vector<double> &XX, vector<int> &spline_mem) {
    double fac=1.0/norm_df_RM_dTheta(The_ref);
    double lx=log(x);
    
    for(int i=0; i<Te_order; i++) {
        XX[i]=fac/pow(The_ref, i)*calc_spline_JC(lx, spline_mem[i])/factorial(i)/x3;
    }
}

void IntegralCNSNopt::Compute_Y(){
    Compute_XX(spline.The_ref, Y, spline.Y); 
}

void IntegralCNSNopt::Compute_D(){
    Compute_XX(spline.The_ref, D, spline.D); 
}

void IntegralCNSNopt::Compute_Q(){
    Compute_XX(spline.The_ref, Q, spline.Q); 
}

void IntegralCNSNopt::Compute_Mcorr(){
    Compute_XX(spline.The_ref, Mcorr, spline.Mcorr); 
}

double IntegralCNSNopt::Calculate_monopole(){
    return Dn_SZ_splines(spline.The_ref, spline.Y);
}

double IntegralCNSNopt::Calculate_dipole(){
    if(betac_order < 1) return 0.0;
    double rnorm = x*exp_mx/dex/dex;
    double r = Dn_SZ_splines(spline.The_ref, spline.D);
    r+=rnorm;
    return r*betac*muc;
}

double IntegralCNSNopt::Calculate_quadrupole(){
    if(betac_order < 2) return 0.0;
    double rnorm = CMBframe ? 11.0/30.0 : -3.0/10.0;
    rnorm *= x*exp_mx/dex/dex*x*(1.0+exp_mx)/dex;

    if(Te_order==0) return betac*betac*(1.5*muc*muc-0.5)*rnorm;
    
    double r = Dn_SZ_splines(spline.The_ref, spline.Q);
    r+=rnorm;
    return r*betac*betac*(1.5*muc*muc-0.5);
}

double IntegralCNSNopt::Calculate_monopole_correction(){
    if(betac_order < 2) return 0.0;
    double rnorm = CMBframe ? x*exp_mx/dex/dex*(x*(1.0+exp_mx)/dex-3.0)/3.0 : 0.0;
    
    double r = Dn_SZ_splines(spline.The_ref, spline.Mcorr);
    r+=rnorm;
    return r*betac*betac;
}

double IntegralCNSNopt::Calculate_kinetic_correction(){
    return Calculate_dipole()+Calculate_quadrupole()+Calculate_monopole_correction();
}
double IntegralCNSNopt::Calculate_All(){
    return Calculate_kinetic_correction()+Calculate_monopole();
}

double IntegralCNSNopt::compute_distortion(string mode){
    //==============================================================================================
    // simplest diversion
    //==============================================================================================
    if (region == 0) {
        return compute_SZ_distortion_asymptotic(x, The, betac, muc, kmax-1, betac_order, mode, CMBframe);
    }

    run_mode="all";
    if(mode=="monopole" || mode=="dipole" || mode=="quadrupole" || mode=="monopole_corr" ||mode=="kin"){
        run_mode = mode;
    }

    if(run_mode=="monopole"){
        return Calculate_monopole();
    } 
    else if(run_mode=="dipole"){ 
        return Calculate_dipole();
    }
    else if(run_mode=="quadrupole"){
        return Calculate_quadrupole();
    }
    else if(run_mode=="monopole_corr"){ 
        return Calculate_monopole_correction();
    }
    else if(run_mode=="kin"){
        return Calculate_kinetic_correction();
    }
    else{
        return Calculate_All();
    }  
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
                                            string mode, int accuracy_level, bool CMBframe){
    IntegralCNSNopt szDistortion = IntegralCNSNopt(x,The,betac,muc,kmax,accuracy_level,betac_order,CMBframe);
    return szDistortion.compute_distortion(mode);   
}

void compute_SZ_distortion_CNSN_basis_opt(vector<double> &Dn, vector<double> x,
                                          double The, double betac, double muc, int kmax, int betac_order,
                                          string mode, int accuracy_level, bool DI, bool CMBframe){                
    int gridpoints = x.size();
    Dn.resize(gridpoints);
    Parameters fp = Parameters(); //This is just to get a value for the Dn_DI conversion 
    IntegralCNSNopt szDistortion = IntegralCNSNopt(x[0],The,betac,muc,kmax,accuracy_level,betac_order,CMBframe);
    Dn[0] = (DI ? pow(x[0],3.0)*fp.rare.Dn_DI_conversion() : 1.0)*szDistortion.compute_distortion(mode);
    for(int k = 1; k < gridpoints; k++){
        szDistortion.Update_x(x[k]);
        Dn[k] = szDistortion.compute_distortion(mode);
        if (DI) { Dn[k] *= pow(x[k],3.0)*fp.rare.Dn_DI_conversion(); }
    }
}

double compute_SZ_distortion_CNSN_basis_opt(int k, Parameters fp, bool CMBframe){
    IntegralCNSNopt szDistortion = IntegralCNSNopt(fp.k_inRange(k), fp, CMBframe);
    return szDistortion.compute_distortion(fp.rare.RunMode);
}

void compute_SZ_distortion_CNSN_basis_opt(vector<double> &Dn, Parameters fp, bool DI, bool CMBframe){
    Dn.resize(fp.gridpoints);
    IntegralCNSNopt szDistortion = IntegralCNSNopt(0, fp, CMBframe);
    Dn[0] = (DI ? pow(fp.xcmb[0],3.0)*fp.rare.Dn_DI_conversion() : 1.0)*fp.Dtau*szDistortion.compute_distortion(fp.rare.RunMode);
    for(int k = 1; k < fp.gridpoints; k++){
        szDistortion.Update_x(fp.xcmb[k]);
        Dn[k] = fp.Dtau*szDistortion.compute_distortion(fp.rare.RunMode);
        if (DI) { Dn[k] *= pow(fp.xcmb[k],3.0)*fp.rare.Dn_DI_conversion(); }
    }
}

//==================================================================================================
// computes the optimal value for kmax given the accuracy goal and required maximal temperature
//==================================================================================================
void determine_kmax_from_accuracy3vec(vector<precision_settings> acc, double Te_max, int &kmax, int &iregmax)
{
    int nreg0=acc[0].determine_Te_region(Te_max);
    int nreg1=acc[1].determine_Te_region(Te_max);
    int nreg2=acc[2].determine_Te_region(Te_max);
    
    //Determine the smallest number of variables i.e., (regions with T<T,max)*kmax
    if(nreg0*acc[0].kmax<nreg1*acc[1].kmax){ 
        iregmax=nreg0; kmax=acc[0].kmax; 
    }
    else{ 
        iregmax=nreg1; kmax=acc[1].kmax; 
    }

    if(iregmax*kmax>nreg2*acc[2].kmax){
        iregmax=nreg2; kmax=acc[2].kmax;
    }
}

void determine_optimal_kmax(int accuracy_level, double Te_max, int &kmax, int &iregmax)
{
    CNSNopt_accuracies::initialise_precisions();

    if(accuracy_level==0) {
        determine_kmax_from_accuracy3vec(CNSNopt_accuracies::acc_I, Te_max, kmax, iregmax);
    }
    else if(accuracy_level==1) {
        determine_kmax_from_accuracy3vec(CNSNopt_accuracies::acc_II, Te_max, kmax, iregmax);
    }
    else if(accuracy_level==2) {
        determine_kmax_from_accuracy3vec(CNSNopt_accuracies::acc_III, Te_max, kmax, iregmax);
    }
    else if(accuracy_level==3) { //on account of there only being the one set of precision_settings here.
        iregmax = CNSNopt_accuracies::acc_IV[0].determine_Te_region(Te_max);
        kmax = CNSNopt_accuracies::acc_IV[0].kmax;
    }
    else { 
        exit_error("determine_optimal_kmax :: accuracy level not available. Should be 0-3 "); 
    }
}

//==================================================================================================
// access basis functions (always in CMB frame)
//==================================================================================================
void compute_Y_CNSNopt(double x, int region, vector<double> &Y){
    int Te_order = Y.size();
    IntegralCNSNopt szDistortion = IntegralCNSNopt(x, region, Te_order);
    Y = szDistortion.Y;
}
void compute_D_CNSNopt(double x, int region, vector<double> &D){
    int Te_order = D.size();
    IntegralCNSNopt szDistortion = IntegralCNSNopt(x, region, Te_order);
    D = szDistortion.D;
}
void compute_Q_CNSNopt(double x, int region, vector<double> &Q){
    int Te_order = Q.size();
    IntegralCNSNopt szDistortion = IntegralCNSNopt(x, region, Te_order);
    Q = szDistortion.Q;
}
void compute_Mcorr_CNSNopt(double x, int region, vector<double> &Mcorr){
    int Te_order = Mcorr.size();
    IntegralCNSNopt szDistortion = IntegralCNSNopt(x, region, Te_order);
    Mcorr = szDistortion.Mcorr;
}
//==================================================================================================
//==================================================================================================
