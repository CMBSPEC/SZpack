//==================================================================================================
//
// Cluster Ne and Te profiles according to definitions of Vikhlinin et al. 2006
//
//==================================================================================================
//==================================================================================================
//
// Author: Jens Chluba  (CITA, University of Toronto)
//
// first implementation: July 2012
// last modification   : Aug  2012
//
//==================================================================================================

#ifndef SZ_CLUSTER_PROFILES_H
#define SZ_CLUSTER_PROFILES_H

#include <string>
#include <vector>

using namespace std;

//==================================================================================================
struct Cluster_param_Te 
{
    // Te-profile
    double T0, rt;
    double a, b, c;
    double Tmin_T0;
    double rcool, acool;
};
    
struct Cluster_param_Ne 
{
    // Ne-profile
    double n0, rc, rs;
    double alpha, beta, eps;
    double n02, rc2, beta2;
};

//==================================================================================================
class SZ_cluster_profiles
{
    
private:
    
    //==============================================================================================
    //
    // parameters
    //
    //==============================================================================================
    Cluster_param_Ne SZpNe;
    Cluster_param_Te SZpTe;
    string label;
    
    template< typename somestream >
    void output_cluster_parameters(somestream &output);

public:
    
    //==============================================================================================
    //
    // Constructor & Destructor
    //
    //==============================================================================================
    SZ_cluster_profiles();
    ~SZ_cluster_profiles();
    
    //==============================================================================================
    // initialize Cluster profile
    //==============================================================================================
    SZ_cluster_profiles(Cluster_param_Ne &pNe, Cluster_param_Te &pTe, string label="");

    //==============================================================================================
    // access cluster parameters
    //==============================================================================================
    Cluster_param_Ne Get_Ne_parameters(){ return SZpNe; }
    Cluster_param_Te Get_Te_parameters(){ return SZpTe; }
    double Get_rc(){ return SZpNe.rc; }
    double Get_rc_cm();
    string Get_label(){ return label; }
    
    //==============================================================================================
    // output cluster parameters
    //==============================================================================================
    void show_cluster_parameters();
    void export_cluster_parameters(string fname);
    
    //==============================================================================================
    // computations (x, y, z all in units of the core-radius)
    //==============================================================================================
    double Ne(double x, double y, double z);  // Ne in cm^-3 
    double Te(double x, double y, double z);  // Te in keV
    
};

#endif

//==================================================================================================
//==================================================================================================
