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

#ifndef SZ_MOMENTS_H
#define SZ_MOMENTS_H

#include <string>
#include <vector>
#include "routines.h"
#include "ClusterFunctions.h"

using namespace std;

extern Parameters parameters;

struct Te_region
{
    double lmin, lmax;
};

class SZ_moment_method
{
private:
    
    //==============================================================================================
    //
    // data structures for storage
    //
    //==============================================================================================
    int nx;
    int nTregs;
    int index_high, nmom_high_one;
    int tot_moments;

    int kmax;
    
    vector<double> Te_pivots;
    vector<double> Te_limits;  
    vector<int> region_indices;

public: 
    simpleMatrix M;  // signal matrix
    simpleMatrix T;  // conversion matrix N(The) (The-The0)^k --> N(The) The^k
    simpleMatrix MT; // transformed signal matrix MT = M*T


    vector<vector<Te_region> > Te_zeros; //Contains the Te_structure
    void determine_Te_structure(double lres);

private:
    void show_matrix(simpleMatrix &B, Parameters fp, bool show_xcmb = true);

    void print_header_for_moment_matrix(ofstream &ofile, string mess, Parameters fp = parameters);
    void export_matrix(string fname, string form, simpleMatrix &M);
    
    void calculateGasMatrix(Parameters fp, int n, vector<double> &Y, int order_T, double x3, int &k, int region_indices = 0, bool CNSN = false);
    void setup_SZ_moment_matrix(int order_T_low, int order_T_high, Parameters fp = parameters);

    template<typename somestream>
    void output_profile_slice(somestream &output, double (*f)(double l), double lsc, int np);
  
public:
    
    //==============================================================================================
    //
    // Constructor & Destructor
    //
    //==============================================================================================
    SZ_moment_method();
    ~SZ_moment_method();
    
    SZ_moment_method(Parameters fp, bool usekmax = true);
    
    //==============================================================================================
    // This constructor allows to directly pick the accuracy level and kmax, with the 
    // definitions according to CSNN 2012. The maximal temperature range that is covered depends 
    // on kmax. Also, if the required Te_max is much smaller than 60keV, far fewer moments have
    // to be provided. In that case setting the moment method up  with the alternative constructor 
    // might be more practical.
    //
    // xcmb[i] contains the required frequency points for the SZ signal
    // kmax choses the maximal number of temperature terms within the given accuracy goal
    // accuracy_level =0, 1, 2, 3determines the required precision of the SZ approximation 
    // order_b == 0, 1, 2 define the order in betac
    //==============================================================================================
    //==============================================================================================
    // This constructor allows to pick the accuracy level but instead of kmax, defining Te_max. 
    // If the required maximal temperature is small compared to 60keV this setup might be beneficial
    // since far fewer moments are required. 
    //
    // xcmb[i] contains the required frequency points for the SZ signal
    // Te_max (in keV) defined the maximal required (expected) electron temperature
    // accuracy_level =0, 1, 2, 3determines the required precision of the SZ approximation 
    // according to the settings of CSNN 2012.
    // order_b == 0, 1, 2 define the order in betac
    //==============================================================================================
    
    
    //==============================================================================================
    //
    // display and access data of SZ moment formalism
    //
    //==============================================================================================
    void show_moment_matrix_M(Parameters fp = parameters);
    void show_moment_matrix_MT(Parameters fp = parameters);
    //
    void show_conversion_matrix_T(Parameters fp = parameters);
    
    void show_Temperature_limits();
    void show_Temperature_pivots();
    
    
    //==============================================================================================
    //
    // information needed to use moment formalism
    //
    //==============================================================================================
    vector<double> Get_Temperature_limits(){ return Te_limits; }
    double  Get_Te_ref(int reg);
    double  Get_The_ref(int reg);

    int Get_kmax(){ return kmax; }
    int Get_total_number_of_moments(){ return tot_moments; }
    int Get_total_number_of_Te_regions(){ return nTregs; }
    int Get_index_of_Te_region(double Te);  // if Te > Temax --> -10 is returned
    int Get_index_high(){ return index_high;}
    int Get_number_of_moments_low(){ return index_high;}
    int Get_number_of_moments_high(){ return nmom_high_one*(nTregs-1);}
    int Get_number_of_moments_for_one_high_region(){ return nmom_high_one;}
    int Get_start_index_of_moment_set(int reg);

    //==============================================================================================
    void show_profile_slice(double (*f)(double l), double lsc, int np);
    void export_profile_slice(string fname, double (*f)(double l), double lsc, int np);
};

#endif

//==================================================================================================
//==================================================================================================
