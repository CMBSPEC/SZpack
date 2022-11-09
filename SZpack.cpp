//==================================================================================================
//
// SZpack functions
//
//==================================================================================================
//
// purpose: computation of the SZ signal according to Chluba, Nagai, Sazonov, Nelson, 2012
//          and Chluba, Switzer, Nagai, Nelson, 2012.
//
// comments: - computations are performed in the single-scattering approximation
//           - polarization effects are neglected
//           - the electron distribution function is assumed to be thermal
//==================================================================================================
//
// Author: Jens Chluba & Elizabeth Lee
//
// first implementation: May  2012
// last modification   : March 2020
//
//==================================================================================================
// 28th Aug 2017: added y-weighted moment method for temperature corrections
// 13th May 2015: improved performance of combo-means-methods && fixed bug for asymptotic derivative
// 20th  Dec: added optimized SZ signal function to minimize number of temperature terms
// 29th  Aug: added expansion of SZ signal around mean values of tau, TeSZ, and betac*muc with
//            higher moments included. 
//  8th  Aug: added function to compute SZ null.
//  6th  Aug: added expansion of SZ signal around mean values of tau, TeSZ, and betac*muc. 
//  4th  Aug: added derivatives of basis functions in the CMB rest frame
// 22th July: added combination of asymptotic expansion + CNSN basis functions. For Te < 2keV 
//            the asymptotic expansion is used while for 2keV < Te < 75keV the basis of CNSN is
//            applied. At 0.01 < x < 30 the precision should be similar to 0.001% for beta < 0.01. 
//            Also added simple sanity checks for approximation functions.
// 10th July: added functions to compute the SZ signal using temperature-velocity moments. These   
//            routines allow taking into account the detailed temperature and velocity structure   
//            of the cluster medium along the line-of-sight.
//--------------------------------------------------------------------------------------------------

//==================================================================================================
//
//  Unless stated otherwise the main parameters are:
//  
//  xo          : observer frame photon frequency xo == h nu / k T0 (T0 == todays CMB temperature)
//  Dtau        : line-of-sight scattering optical depth in the cluster frame
//  Te          : temperature of the (cluster frame) electron gas in keV
//  betac       : velocity beta=v/c of the cluster in the CMB frame
//  muc         : direction cosine of the cluster velocity with respect to the line-of-sight in  
//                the CMB rest frame
//  betao       : velocity beta=v/c of the observer in the CMB frame 
//                (from CMB dipole betao=1.241e-3)
//  muo         : direction cosine for the line-of-sight with respect to the observers velocity. 
//                The angle is measured in the observer frame
//
//==================================================================================================

#include "SZpack.h"

#include <string>
#include <cmath>
#include <vector>

#include "routines.h"
#include "Parameters.h"
#include "SZ_CNSN_basis.h"
#include "SZ_CNSN_basis.opt.h"
#include "SZ_Integral.3D.h"
#include "SZ_Integral.5D.h"
#include "SZ_asymptotic.h"
#include "SZ_nonrelativistic.h"
#include "global_functions.h"

using namespace std;

bool CNSN2012_convention = true;
void setConvention(bool UseNozawaConvention){
    CNSN2012_convention = !UseNozawaConvention;
}

// Note: If The or Te = 0.0, none of these will give the appropriate return. (Since y is calculated
// as Dtau*The). To calculate the nonrelativistic signal, use the specific function.
//==================================================================================================
double compute_signal_nonrelativistic(double x, Parameters &fp){
    return fp.Dtau*compute_SZ_distortion_nonrelativistic(x, fp);
}

double compute_signal_5D(double x, Parameters &fp){
    return fp.Dtau*compute_SZ_distortion_Patterson_5D(x, fp);
}

double compute_signal_3D(double x, Parameters &fp){
    return fp.Dtau*compute_SZ_distortion_Patterson_3D(x, fp);
}

//double compute_signal_Kernel(double x, Parameters &fp){
    //TODO: Add in Sanity Checks!
//    return compute_SZ_distortion_Kernel(x, fp);
//}

double compute_signal_asymptotic(double x, Parameters &fp, bool CMBframe){
    if(fp.Te>20.0) { 
        print_error("compute_SZ_signal_asymptotic :: Temperature really (!) high (Te = " + to_string(fp.Te) + " keV)");
    }
    
    if(fp.T_order>10){ fp.T_order=10; }
    if(fp.beta_order>2){ fp.beta_order=2; }

    return fp.Dtau*compute_SZ_distortion_asymptotic(x, fp, CMBframe);
}

double compute_signal_asymptotic(double x, Parameters &fp){
    return compute_signal_asymptotic(x,  fp, true);
}

double compute_signal_CNSN(double x, Parameters &fp, bool CMBframe){
    if(x<0.01 || x>50.0) { 
        exit_error("compute_SZ_signal_CNSN_basis :: You are outside of grid");
    }

    if(fp.Te<2.0 || fp.Te>75.0){ 
        exit_error("compute_SZ_signal_CNSN_basis :: Temperature too high/low: " + to_string(fp.Te) + "keV");
    }

    if(fp.T_order>20){ fp.T_order=20; }
    if(fp.beta_order>2){ fp.beta_order=2; }

    return fp.Dtau*compute_SZ_distortion_CNSN_basis(x, fp, CMBframe);
}

double compute_signal_CNSN(double x, Parameters &fp){
    return compute_signal_CNSN(x, fp, true);
}

double compute_signal_CNSN_opt(double x, Parameters &fp){
    //==============================================================================================
    // simple sanity checks
    //==============================================================================================
    if(x<0.01 || x>50.0){
        exit_error("compute_SZ_signal_CNSN_basis_opt :: You are outside of grid");
    }

    return fp.Dtau*compute_SZ_distortion_CNSN_basis_opt(x, fp, true);
}

double compute_signal_combo(double x, Parameters &fp, bool CMBframe){
    if(fp.Te<2.0){
        fp.T_order = 10;
        fp.beta_order = 2;
        return compute_signal_asymptotic(x, fp, CMBframe); 
    }
    fp.T_order = 20;
    fp.beta_order = 2;
    return compute_signal_CNSN(x, fp, CMBframe); 
}

double compute_signal_combo(double x, Parameters &fp){
    return compute_signal_combo(x, fp, true);
}

//==================================================================================================
// Derivatives (The^k d^k_dThe /k!) (d^m_dbetapara /m!) (d^l_beta2perp /l!) S(...)
// in the CMB frame for a resting observer. Maximal orders in The and betac are used to compute 
// the derivatives.
//
// constraints: dThe<=4; dbeta_para<=2; dbeta2_perp<=1;
//==================================================================================================
void Dcompute_signal_combo_CMB(double x, Parameters &fp, bool yw)
{
    if(fp.Te<2.0) Dcompute_SZ_distortion_asymptotic(x, fp, true); 
    else
    {
        //==========================================================================================
        // simple sanity checks
        //==========================================================================================
        if(x<0.01 || x>50.0){ 
            exit_error("Dcompute_SZ_signal_combo_CMB :: You are outside of grid");
        }
        
        if(fp.Te>75.0){ 
            exit_error("Dcompute_SZ_signal_combo_CMB :: Temperature too high: " + to_string(fp.Te) + "keV");
        }

        Dcompute_SZ_distortion_CNSN_basis(x, fp, true);
    }
    
    for(int g = 0; g <= fp.D.dThe; g++) {fp.D.dDn_dThe[g] *= fp.Dtau;}

    // nb, since Dcompute calculates assuming the tau is constant, and not y, as necessary for y weighting
    // we must rewrite our The derivatives noting, dDn_dThe[k] = The^k/k! d^k/dThe^k (<tau> The f(nu,The))
    // and we require, H[k] = The^k/k! d^k/dThe^k (<y> f(nu,The))
    if(yw) {
        for(int g = 1; g <= fp.D.dThe; g++) {
            fp.D.dDn_dThe[g] = fp.D.dDn_dThe[g]-fp.D.dDn_dThe[g-1];
        }
    }
}

//==================================================================================================
//
// compute x derivatives from the combo method
//
//==================================================================================================

double Dcompute_signal_combo_for_x(double x0, Parameters fp, int dx){
    double eps = 0.001;
    double Ix=pow(x0, 3)*compute_signal_combo(x0, fp, true);
    double xp=x0*(1.0+eps), xm=x0*(1.0-eps);
    double xp2=x0*(1.0+2*eps), xm2=x0*(1.0-2*eps);
    double Ixp=pow(xp, 3)*compute_signal_combo(xp, fp, true);
    double Ixm=pow(xm, 3)*compute_signal_combo(xm, fp, true);
    double Ixp2=pow(xp2, 3)*compute_signal_combo(xp2, fp, true);
    double Ixm2=pow(xm2, 3)*compute_signal_combo(xm2, fp, true);

    if (dx == 1) {
        return (-Ixp2+8.0*Ixp-8.0*Ixm+Ixm2)/12.0/eps;
    }
    if (dx == 2) {
        return (-Ixp2+16.0*Ixp-30.0*Ix+16.0*Ixm-Ixm2)/12.0/pow(eps, 2) / (2.0);
    }
    if (dx == 3) {
        return (Ixp2-2.0*Ixp+2.0*Ixm-Ixm2)/2.0/pow(eps, 3) / (6.0);
    }
    if (dx == 4) {
        return (Ixp2-4.0*Ixp+6.0*Ix-4.0*Ixm+Ixm2)/pow(eps, 4) / (24.0);
    }
    if (dx != 0) {
        print_error("Dcompute_signal_combo_for_x: dx not recognised. Must be an integer between 0 and 4 "
                    "(dx = " + to_string(dx) + ")");
    }
    return Ix;
}

//==================================================================================================

double compute_signal_means(double x, Parameters &fp, bool yw){
    // mean signal + temperature dispersion term
    fp.D.setValues(2,0,0);
    Dcompute_signal_combo_CMB(x, fp, yw);
    
    double r = fp.D.dDn_dThe[0]+fp.D.dDn_dThe[2]*fp.means.Omega;

    // velocity - temperature cross term
    if(fp.means.Sigma!=0.0)
    {
        fp.D.setValues(1,1,0);
        Dcompute_signal_combo_CMB(x, fp, yw);
        r += fp.D.dDn_dThe[1]*fp.means.Sigma;
    }
    
    // betac_parallel dispersion term
    if(fp.means.kappa!=0.0)
    {
        fp.D.setValues(0,2,0);
        Dcompute_signal_combo_CMB(x, fp, yw);
        r += fp.D.dDn_dThe[0]*fp.means.kappa;
    }

    // betac_perp dispersion term
    if(fp.calc.betac2_perp!=0.0)
    {
        fp.D.setValues(0,0,1);
        Dcompute_signal_combo_CMB(x, fp, yw);
        r += fp.D.dDn_dThe[0]*fp.calc.betac2_perp;
    }
    return r;
}

double compute_signal_means_tw(double x, Parameters &fp){
    //This uses the dimensionless tau-weighted terms
    return compute_signal_means(x, fp, false);
}

double compute_signal_means_yw(double x, Parameters &fp){
    //This uses the dimensionless y-weighted terms as described in Lee, Chluba, Kay & Barnes 2020
    return compute_signal_means(x, fp, true);
}

//==================================================================================================
// 
// Compute null of SZ signal using expansion around mean values of tau, TeSZ, and betac*muc. This  
// routine should reproduce the full numerical result for Te < 75keV, and beta < 0.01 with high 
// precision, assuming smooth cluster profile. (added 8th Aug, 2012, JC)
//
//==================================================================================================
double compute_null_of_SZ_signal(Parameters fp) {
    return find_root_brent([&fp](double x0) { return compute_signal_means(x0, fp, false);}, 0.1, 10.0, 1.0e-5);
}

//==================================================================================================
//
// Extended versions of expansions around average values. Up to The^4 terms can be included for
// the thSZ signal and up to The^3 for first order cross terms (temperature corrections to higher 
// order velocity terms can in principle be activated but have been omitted at this point). 
// The parameters are:
// 
// omega[0..2] == omega^(1..3)         [i.e., O(The), O(The^2), O(The^3) & O(The^4)]
// sigma[0..2] == sigma^(1..3, 1)      [i.e., O(betac The), O(betac The^2), O(betac The^3)]
// kappa       == kappa^(1)            [i.e., O(betac^2 muc^2)]
// betac2_perp == <beta_c^2 (1-muc^2)> 
//
// according to definitions of CSNN 2012. Values that are not need should be set to zero.
// (added 29th Aug, 2012, JC)
//
//==================================================================================================
double compute_signal_means_ex(double x, Parameters &fp, bool yw)
{
    // mean signal + temperature dispersion term
    fp.D.setValues(4,0,0);
    Dcompute_signal_combo_CMB(x, fp, yw);
    double r = fp.D.dDn_dThe[0]+fp.D.dDn_dThe[2]*fp.means.omegas[0]+fp.D.dDn_dThe[3]*fp.means.omegas[1]+fp.D.dDn_dThe[4]*fp.means.omegas[2];
    
    // velocity - temperature cross term
    if(fp.means.sigmas[0]!=0.0 || fp.means.sigmas[1]!=0.0 || fp.means.sigmas[2]!=0.0)
    {
        fp.D.setValues(3,1,0);
        Dcompute_signal_combo_CMB(x, fp, yw);
        r += fp.D.dDn_dThe[1]*fp.means.sigmas[0]+fp.D.dDn_dThe[2]*fp.means.sigmas[1]+fp.D.dDn_dThe[3]*fp.means.sigmas[2];
    }
    
    // betac_parallel dispersion term
    if(fp.means.kappa != 0.0)
    {
        fp.D.setValues(0,2,0);
        Dcompute_signal_combo_CMB(x, fp, yw);
        r += fp.D.dDn_dThe[0]*fp.means.kappa;
    }
    
    // betac_perp dispersion term
    if(fp.calc.betac2_perp != 0.0)
    {
        fp.D.setValues(0,0,1);
        Dcompute_signal_combo_CMB(x, fp, yw);
        r += fp.D.dDn_dThe[0]*fp.calc.betac2_perp;
    }
    return r;    
}  

//==================================================================================================

double compute_signal_RelCorrs(double x, Parameters &fp){
    double dum = compute_signal_combo(x, fp);
    dum -= compute_signal_nonrelativistic(x, fp);
    return dum;
}

double compute_signal_TDispersion(double x, Parameters &fp){
    double dum=compute_signal_means(x, fp, false);

    double Omega = fp.means.Omega, Sigma = fp.means.Sigma;
    fp.means.Omega = fp.means.Sigma = 0.0;
    dum -= compute_signal_means(x, fp, false);
    fp.means.Omega = Omega;
    fp.means.Sigma = Sigma;
    return dum;
}

void compute_signal_TwoTemperatures(vector<double> &Dn, double ftau, double DT_T, Parameters &fp, 
                                    bool DI, bool CMBframe){
    Dn.resize(fp.gridpoints);
    vector<double> temp1Dn, temp2Dn;
    
    Parameters temp1 = Parameters(), temp2 = Parameters();
    temp1.copyParameters(fp);
    temp2.copyParameters(fp);
    temp1.Dtau *= (1.0-ftau);
    temp2.Dtau *= ftau;
    temp2.updateT(temp2.Te*(1.0+DT_T));

    compute_signal_combo(temp1Dn, temp1, DI, CMBframe);
    compute_signal_combo(temp2Dn, temp2, DI, CMBframe);
    for (int k = 0; k < fp.gridpoints; k++){
        Dn[k] = temp1Dn[k]+temp2Dn[k];
    }
}

//==================================================================================================
// Vector forms of all of the signals
//==================================================================================================

void compute_signal_nonrelativistic(vector<double> &Dn, Parameters &fp, bool DI){
    compute_SZ_distortion_nonrelativistic(Dn, fp, DI);
}

void compute_signal_5D(vector<double> &Dn, Parameters &fp, bool DI){
    compute_SZ_distortion_Patterson_5D(Dn, fp, DI);
}

void compute_signal_3D(vector<double> &Dn, Parameters &fp, bool DI){
    compute_SZ_distortion_Patterson_3D(Dn, fp, DI);
}

//void compute_signal_Kernel(vector<double> &Dn, Parameters &fp, bool DI){
    //TODO: This!
//}

void compute_signal_asymptotic(vector<double> &Dn, Parameters &fp, bool DI, bool CMBframe){
    compute_SZ_distortion_asymptotic(Dn, fp, DI, CMBframe); 
}

void compute_signal_CNSN(vector<double> &Dn, Parameters &fp, bool DI, bool CMBframe){
    compute_SZ_distortion_CNSN_basis(Dn, fp, DI, CMBframe); 
}

void compute_signal_CNSN_opt(vector<double> &Dn, Parameters &fp, bool DI, bool CMBframe){
    compute_SZ_distortion_CNSN_basis_opt(Dn, fp, DI, CMBframe); 
}

void compute_signal_combo(vector<double> &Dn, Parameters &fp, bool DI, bool CMBframe){
    if(fp.Te<2.0){
        fp.T_order = 10;
        fp.beta_order = 2;
        compute_SZ_distortion_asymptotic(Dn, fp, DI, CMBframe); 
    }
    else{
        fp.T_order = 20;
        fp.beta_order = 2;
        compute_SZ_distortion_CNSN_basis(Dn, fp, DI, CMBframe); 
    }
}

void compute_signal_precise(vector<double> &Dn, Parameters &fp, bool DI){
    if (fp.Te > 75){ compute_signal_3D(Dn, fp, DI); }
    else{ compute_signal_combo(Dn, fp, DI, false); }
}

void compute_signal_means(vector<double> &Dn, Parameters &fp, bool DI, bool yw){
    Dn.resize(fp.gridpoints);
    for (int k = 0; k < fp.gridpoints; k++){
        Dn[k] = compute_signal_means(fp.xcmb[k], fp, yw);
        if (DI) { Dn[k] *= pow(fp.xcmb[k],3.0)*fp.rare.Dn_DI_conversion(); }
    }
}

void compute_signal_means_ex(vector<double> &Dn, Parameters &fp, bool DI, bool yw){
    Dn.resize(fp.gridpoints);
    for (int k = 0; k < fp.gridpoints; k++){
        Dn[k] = compute_signal_means_ex(fp.xcmb[k], fp, yw);
        if (DI) { Dn[k] *= pow(fp.xcmb[k],3.0)*fp.rare.Dn_DI_conversion(); }
    }
}

void compute_signal_RelCorrs(vector<double> &Dn, Parameters &fp, bool DI){
    Dn.resize(fp.gridpoints);
    for (int k = 0; k < fp.gridpoints; k++){
        Dn[k] = compute_signal_RelCorrs(fp.xcmb[k], fp);
        if (DI) { Dn[k] *= pow(fp.xcmb[k],3.0)*fp.rare.Dn_DI_conversion(); }
    }
}

void compute_signal_TDispersion(vector<double> &Dn, Parameters &fp, bool DI){
    Dn.resize(fp.gridpoints);
    for (int k = 0; k < fp.gridpoints; k++){
        Dn[k] = compute_signal_TDispersion(fp.xcmb[k], fp);
        if (DI) { Dn[k] *= pow(fp.xcmb[k],3.0)*fp.rare.Dn_DI_conversion(); }
    }
}

void Dcompute_signal_combo_CMB(vector<double> &Dn, Parameters &fp, 
                               int dThe, int dbeta_para, int dbeta2_perp, bool yw, bool DI){
    Dn.resize(fp.gridpoints);
    fp.D.dThe = dThe;
    fp.D.dbeta_para = dbeta_para;
    fp.D.dbeta2_perp = dbeta2_perp;
    for (int k=0; k< fp.gridpoints; k++){
        Dcompute_signal_combo_CMB(fp.xcmb[k], fp, yw);
        Dn[k] = fp.D.dDn_dThe[dThe];     
        if (DI) { Dn[k] *= pow(fp.xcmb[k],3.0)*fp.rare.Dn_DI_conversion(); }
    }
}

void Dcompute_signal_combo_for_x(vector<double> &Dn, Parameters fp, int dx){
    Dn.resize(fp.gridpoints);
    for (int k=0; k< fp.gridpoints; k++){
        Dn[k] = Dcompute_signal_combo_for_x(fp.xcmb[k], fp, dx)*fp.rare.Dn_DI_conversion();
    }
}

//==================================================================================================
// Vector methods to transform the signal from DI/Dn to DT/Tcmb
// i.e., the signal expressed as variations to the cmb temperature
// Note these take Dn and not DI.
//==================================================================================================

// This method converts through an inversion of the boltzmann distribution at each frequency.
void convert_signal_DT(vector<double> &DT, vector<double> xcmb, vector<double> Dn){
    int np = xcmb.size();
    DT.resize(np);
    for (int k = 0; k < np; k++){
        double expxMinusOne = exp(xcmb[k])-1;
        DT[k] = (xcmb[k]/log((expxMinusOne/(1+expxMinusOne*Dn[k]))+1))-1;
    }
}

// This method converts by using the first derivative and ignoring higher order terms.
void convert_signal_DT_approx(vector<double> &DT, vector<double> xcmb, vector<double> Dn){
    int np = xcmb.size();
    DT.resize(np);
    for (int k = 0; k < np; k++){
        DT[k] = Dn[k]*(exp(xcmb[k])-1)*(exp(xcmb[k])-1)/xcmb[k]/exp(xcmb[k]);
    }
}

//==================================================================================================
//==================================================================================================
