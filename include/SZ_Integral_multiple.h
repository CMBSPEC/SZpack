//==================================================================================================
//
// Computation multiple scattering correction to thermal SZ effect according to Chluba et al. 2013 
//
//==================================================================================================
//
// Author: Jens Chluba  (Johns Hopkins University)
//
// first implementation: Jan 2013
// last modification   : Jan 2013
//
//==================================================================================================

#ifndef SZ_MULTIPLE_H
#define SZ_MULTIPLE_H

using namespace std;

//==================================================================================================
//
// Compton scattering kernel ass defined in Chluba et al. 2013
//
// eps_Int : relative accuracy for numerical integration (lower than 10^-6 is hard to achieve)
//
//==================================================================================================
double P_l_Kernel(int l, double s, double The, double eps_Int=1.0e-4);
double P_l_Kernel_Norm(int l, double The, double eps_Int=1.0e-4);

void create_2D_Kernel_Lib(string filename, int lmax, int nps, int npT,
                          double The_min, double The_max, double
                          eps_Int=1.0e-4);

void read_2D_Kernel_Lib(string filename);

double Compute_Kernel_interpol(int l, double s, double The);


//==================================================================================================
//
// multiple scattering correction using the scattering kernel library
//
// eps_Int : relative accuracy for numerical integration (lower than 10^-6 is hard to achieve)
//
//==================================================================================================
double compute_SZ_distortion_Patterson_multiple_Kernel(double x, int l,
                                                       double The,
                                                       double eps_Int=1.0e-3);

double compute_SZ_distortion_Patterson_multiple_Kernel_2D(double x, int l,
                                                          double The,
                                                          double eps_Int=1.0e-3);


//==================================================================================================
//
// 3D integration for multiple scattering terms carried out using Patterson scheme
//
// eps_Int : relative accuracy for numerical integration (lower than 10^-6 is hard to achieve)
//
// mode == "SEC"      --> second scattering correction of thermal SZ effect (x is irrelevant)
// mode == "SIG"      --> compute scattering cross section correction
// mode == "CMB"      --> compute scattering of CMB multipoles
//
//==================================================================================================
double compute_SZ_distortion_Patterson_multiple(double x, int l,
                                                double The,
                                                string mode,
                                                double eps_Int=1.0e-3);

//==================================================================================================
//
// asymptotic expansion of multiple scattering terms according to Chluba et al. 2013. A maximum
// of 9 temperature correction terms, i.e., O(The^11) can be included (0<=Te_corr_order<=9)
//
//==================================================================================================
double compute_SZ_distortion_asym_multiple(double x, int l, double The, int Te_corr_order);

// access basis functions for different l
void compute_Yl0_k(double x, int l, vector<double> &Yk);

#endif

//==================================================================================================
//==================================================================================================
