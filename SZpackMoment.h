//==================================================================================================
//
// SZpack functions unique to the Moment Method Run
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
// Author: Elizabeth Lee
// Based on Work by Jens Chluba
//
//==================================================================================================

#ifndef SZPack_Moment_H
#define SZPack_Moment_H

#include "output.h"
#include "SZ_moment_method.h"
#include "SZ_Integral_moment.h"
#include "IsothermalBestFit.h"

using namespace std;

extern Parameters parameters;
extern integrationValues IValues;

//==================================================================================================
//
// compute distortions using temperature-velocity moments
//
//==================================================================================================
void compute_isothermal_moments(bool Mtransformed, vector<double> &mv, SZ_moment_method &SZM, Parameters fp = parameters);

//==================================================================================================
//
// Compute the SZ signal using temperature-velocity moments. This routine is in particular 
// useful when simply trying to determine the different moments from observations, however, in 
// certain situations the user might want to compute the moment values independently from, e.g., 
// simulation data (especially if the cluster profiles are not very smooth, i.e. can be 
// approximated using interpolation). 
// 
// The moments have to be computed for different temperature ranges. For the  low temperature 
// gas the variable y^(k) = int The^(k+1) dtau is used, while the high temperature moments have 
// the form z^(k) = int G(The, k) dtau. The velocity moments are similar. The SZ signal is then 
// given by S = M * m. The moments have to be calculated for the different temperature ranges 
// setup by the code. 
// 
// input : mv(i) contains the ordered temperature-velocity moments. 
//         variable chooses the form of the high temperature moments. 
//
//         0: G(The, k) = N(The) (The-The0)^k
//         1: G(The, k) = N(The) The^k
//
// output: xv(i) frequencies of DIv(i); DIv(i): SZ signal in CMB frame.
//
//==================================================================================================
void compute_signal_moments(const vector<double> &mv, vector<double> &DIv, simpleMatrix &Mptr, Parameters fp, SZ_moment_method SZM);

void compute_cluster_moments(double zsc, vector<double> &mv, bool AdjustTemperature, double lres,
                             SZ_moment_method SZM, Parameters &fp);

void compute_cluster_moments_3D(double lx0, double ly0, double lsc, double lang, vector<double> &mv,
                                bool AdjustTemperature, double lres, SZ_moment_method SZM, Parameters &fp);

//==================================================================================================
//
// Compute the SZ signal using temperature-velocity moments. Here the moments are computed 
// internally according to the Ne_f(l), Te_f(l), betac_f(l), and muc_f(l) profile functions   
// provided by the user. In this l=z/zsc is so that z is a length along the line of sight and   
// 2*zsc is the total line-of-sight interval. It is assumed that the integration range is 
// -1 < l < 1, or -zsc < z <zsc (no particular symmetry assumed).  
// 
// input : Ne_f(l) and Te_f(l) functions for cluster profiles. [ Ne ] = 1/cm^3 and [ Te ]= keV
//       : betac_f(l) and muc_f(l) velocity of volume element and muc = ^beta . ^gamma  
//         If order_b==0 these functions will not be called.
//
// output : mv(i) contains the ordered temperature-velocity moments. 
//        : variable chooses the form of the high temperature moments (like above). 
//
//          0: G(The, k) = N(The) (The-The0)^k
//          1: G(The, k) = N(The) The^k
//
//        : xv(i) frequencies of DIv(i); DIv(i): SZ signal in CMB frame.
//        : lres determines the minimal size of temperature structures. Temperature  
//          structures that are smaller than ~zsc*lres will not be resolved.
//
//==================================================================================================
void compute_signal_moments_integrals(double zsc, vector<double> &mv, vector<double> &DIv,
                                      vector<double> &TSZetc, bool AdjustTemperature, double lres,
                                      SZ_moment_method SZM, Parameters &fp, simpleMatrix &Mptr);

//==================================================================================================
void compute_signal_moments_integrals(double lx0, double ly0, double lsc, double lang, vector<double> &mv, 
                               vector<double> &DIv, vector<double> &TSZetc, bool AdjustTemperature,
                               double lres, SZ_moment_method SZM, Parameters &fp, simpleMatrix &Mptr);

//==============================================================================================
//
// Compute the SZ signal using temperature-velocity moments, but assuming that the variance of
// the temperature and velocity field along the line-of-sight is small. In this case the SZ 
// signal is related to S_iso(tau, TeSZ) and its derivatives with respect to Te. The functions 
// and parameters are similar to those above.
//
//==============================================================================================
void compute_signal_moments_smooth(const vector<double> &omegav,
                                                        const vector<double> &TSZetc,
                                                        vector<double> &DIv, 
                                                        Parameters fp);

//==================================================================================================
void compute_signal_moments_smooth(double zsc, vector<double> &omegav, vector<double> &DIv,
                                   vector<double> &TSZetc, int kmax, Parameters fp);

//==================================================================================================
// For Computation of degeneracy functions
//==================================================================================================
void calculateSignalMatrices(vector<double> &M, vector<double> &A, Parameters &fp);
//==================================================================================================
// Functions for VIK
//==================================================================================================
void set_tau_ySZ_TSZ_VIK(double &tau, double &ySZ, double &TSZ);
//==================================================================================================
//==================================================================================================
#endif
