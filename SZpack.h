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

#ifndef SZPACK_H
#define SZPACK_H

#include <cmath>
#include <vector>
#include "Parameters.h"

#include "SZ_Integral.5D.h"
#include "SZ_Integral.3D.h"
#include "SZ_asymptotic.h"
#include "SZ_CNSN_basis.h"
#include "SZ_CNSN_basis.opt.h"
//#include "SZ_Integral.Kernel.h"
#include "SZ_nonrelativistic.h"

using namespace std;

const string SZpack_version="SZpack v2.0";

//==================================================================================================
//
// As discussed by Chluba et el. 2012, Nozawa 1998 and 2006 used another convention for the optical
// depth variable.
// 
//      Dtau --> gammac (1-betac muc) Dtau 
//
// for an observer at rest in the CMB frame. To use this convention call 'setConvention(true);' 
// before executing the other routines. To reset to the convention of Chluba et al. call 
// 'setConvention(false);', which is the default. These functions only have to be called once. 
// Subsequent calls of routines will use the convention that was set last.
// TODO: currently this is only used within the output_SZ_distortion() method.
//
//==================================================================================================
static bool CNSN2012_convention;

void setConvention(bool UseNozawaConvention);

double compute_signal_nonrelativistic(double x, Parameters &functionParameters);
double compute_signal_5D(double x, Parameters &functionParameters);
double compute_signal_3D(double x, Parameters &functionParameters);
//double compute_signal_Kernel(double x, Parameters &functionParameters);
double compute_signal_asymptotic(double x, Parameters &functionParameters, bool CMBframe);
double compute_signal_asymptotic(double x, Parameters &functionParameters);
double compute_signal_CNSN(double x, Parameters &functionParameters, bool CMBframe);
double compute_signal_CNSN(double x, Parameters &functionParameters);
double compute_signal_CNSN_opt(double x, Parameters &functionParameters);
double compute_signal_combo(double x, Parameters &functionParameters, bool CMBframe);
double compute_signal_combo(double x, Parameters &functionParameters);
double compute_signal_means(double x, Parameters &functionParameters, bool yw);
double compute_signal_means_tw(double x, Parameters &functionParameters);
double compute_signal_means_yw(double x, Parameters &functionParameters);
double compute_signal_RelCorrs(double x, Parameters &functionParameters);
double compute_signal_TDispersion(double x, Parameters &functionParameters);

//==================================================================================================
// Derivatives (The^k d^k_dThe /k!) (d^m_dbetapara /m!) (d^l_beta2perp /l!) S(...)
// in the CMB frame for a resting observer. Maximal orders in The and betac are used to compute 
// the derivatives.
//
// constraints: dThe<=4; dbeta_para<=2; dbeta2_perp<=1;
//==================================================================================================
void Dcompute_signal_combo_CMB(double x, Parameters &fp, bool yw=false);

//==================================================================================================
// compute x derivatives from the combo method
//==================================================================================================

double Dcompute_signal_combo_for_x(double x0, Parameters fp, int dx);

//==================================================================================================
// 
// Expansion of SZ signal around mean values of tau, TeSZ, and betac*muc. This routine should 
// reproduce the full numerical result for 0.01 < x < 30, Te < 75keV, and beta < 0.01 with 
// precision similar to 0.001%, assuming smooth cluster profile. (added 6th Aug, 2012, JC)
//
//==================================================================================================
double compute_null_of_SZ_signal(Parameters fp);

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
double compute_signal_means_ex(double x, Parameters &fp, bool yw=false);

// Calculates the combined signal to fp.Dn and the best fit omegas values for the two T signal
void compute_signal_TwoTemperatures(vector<double> &Dn, double ftau, double DT_T, Parameters &fp, 
                                    bool DI, bool CMBframe=true);

//==================================================================================================
// Vector forms of all of the signals
// These contain the Dtau factor already.
//==================================================================================================
void compute_signal_nonrelativistic(vector<double> &Dn, Parameters &fp, bool DI);
void compute_signal_5D(vector<double> &Dn, Parameters &fp, bool DI);
void compute_signal_3D(vector<double> &Dn, Parameters &fp, bool DI);
//void compute_signal_Kernel(vector<double> &Dn, Parameters &fp, bool DI);
void compute_signal_asymptotic(vector<double> &Dn, Parameters &fp, bool DI, bool CMBframe=false);
void compute_signal_CNSN(vector<double> &Dn, Parameters &fp, bool DI, bool CMBframe=false);
void compute_signal_CNSN_opt(vector<double> &Dn, Parameters &fp, bool DI, bool CMBframe=false);
void compute_signal_combo(vector<double> &Dn, Parameters &fp, bool DI, bool CMBframe=false);
void compute_signal_precise(vector<double> &Dn, Parameters &fp, bool DI);
void compute_signal_means(vector<double> &Dn, Parameters &fp, bool DI, bool yw=false);
void compute_signal_means_ex(vector<double> &Dn, Parameters &fp, bool DI, bool yw=false);
void compute_signal_RelCorrs(vector<double> &Dn, Parameters &fp, bool DI);
void compute_signal_TDispersion(vector<double> &Dn, Parameters &fp, bool DI);

// Vector forms of the Dcomputes. These return the derivative at each x in fp.xcmb.
void Dcompute_signal_combo_CMB(vector<double> &Dn, Parameters &fp, 
                               int dThe, int dbeta_para, int dbeta2_perp, bool yw, bool DI);
void Dcompute_signal_combo_for_x(vector<double> &Dn, Parameters fp, int dx);


//==================================================================================================
// Vector methods to transform the signal from DI/Dn to DT/Tcmb
// i.e., the signal expressed as variations to the cmb temperature
//==================================================================================================

// This method converts through an inversion of the boltzmann distribution at each frequency.
void convert_signal_DT(vector<double> &DT, vector<double> xcmb, vector<double> Dn);

// This method converts by using the first derivative and ignoring higher order terms.
void convert_signal_DT_approx(vector<double> &DT, vector<double> xcmb, vector<double> Dn);

#endif

//==================================================================================================
//==================================================================================================
