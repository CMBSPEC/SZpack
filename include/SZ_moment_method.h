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

using namespace std;
/*
class SZ_moment_method
{
    
private:
    
    //==============================================================================================
    //
    // data structures for storage
    //
    //==============================================================================================
    int nx;
    int accuracy_level;
    int kmax;
    int nTregs;
    int order_b;
    int index_high, nmom_high_one;
    int tot_moments;
    double Te_max;
    
    vector<double> Te_pivots;
    vector<double> Te_limits;  
    vector<int> region_indices;
    
    vector<double> M;  // signal matrix
    vector<double> T;  // conversion matrix N(The) (The-The0)^k --> N(The) The^k
    vector<double> MT; // transformed signal matrix MT = M*T
    vector<double> x;  // frequencies  
    
    void show_moment_matrix(vector<double> &M);
    void show_conversion_matrix(vector<double> &C);

    void print_header_for_moment_matrix(ofstream &ofile, string mess);
    void export_moment_matrix(string fname, string form, vector<double> &M);
    
    void setup_SZ_moment_matrix(int order_T_low, int order_T_high, int order_b);
    
    template< typename somestream >
    void output_moment_vector(somestream &output, const vector<double> &mv);

    template<typename somestream>
    void output_profile_slice(somestream &output, double (*f)(double l), double lsc, int np);

    void determine_Te_structure(double lres);
    
    void compute_cluster_moments(double zsc, vector<double> &mv, int Te_order, 
                                 int variable, double lres);  
    
    void compute_cluster_moments_3D(double lx0, double ly0, double lsc, double lang, 
                                    vector<double> &mv, int Te_order, 
                                    int variable, double lres);    
public:
    
    //==============================================================================================
    //
    // Constructor & Destructor
    //
    //==============================================================================================
    SZ_moment_method();
    ~SZ_moment_method();
    
    
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
    SZ_moment_method(const vector<double> &xcmb, int kmax, int accuracy_level, int order_b);

    
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
    SZ_moment_method(const vector<double> &xcmb, double Te_max, int accuracy_level, int order_b);
    
    
    //==============================================================================================
    //
    // display and access data of SZ moment formalism
    //
    //==============================================================================================
    void show_x_vector();
    void show_moment_matrix_M();
    void show_moment_matrix_MT();
    //
    void show_conversion_matrix_T();
    void show_moment_vector(const vector<double> &mv);
   
    void export_x_vector(string fname);
    void export_moment_matrix_M(string fname);
    void export_moment_matrix_MT(string fname);
    void export_moment_vector(string fname, const vector<double> &mv);
    
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
    double Get_Te_max(){ return Te_max; }

    int Get_kmax(){ return kmax; }
    int Get_accuracy_level(){ return accuracy_level; }
    int Get_order_b(){ return order_b; }
    int Get_total_number_of_moments(){ return tot_moments; }
    int Get_total_number_of_Te_regions(){ return nTregs; }
    int Get_index_of_Te_region(double Te);  // if Te > Temax --> -10 is returned
    int Get_index_high(){ return index_high;}
    int Get_number_of_moments_low(){ return index_high;}
    int Get_number_of_moments_high(){ return nmom_high_one*(nTregs-1);}
    int Get_number_of_moments_for_one_high_region(){ return nmom_high_one;}
    int Get_start_index_of_moment_set(int reg);
    
    
    //==============================================================================================
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
    //==============================================================================================
    void compute_SZ_signal_moments(const vector<double> &mv, 
                                   vector<double> &xv, 
                                   vector<double> &DIv,
                                   int variable);

    //==============================================================================================
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
    //==============================================================================================
    void compute_SZ_signal_moments(double (*Ne_f)(double l), 
                                   double (*Te_f)(double l), 
                                   double (*betac_f)(double l), 
                                   double (*muc_f)(double l), 
                                   double zsc,
                                   vector<double> &mv, 
                                   vector<double> &xv, 
                                   vector<double> &DIv,
                                   vector<double> &TSZetc,
                                   int variable, double lres);
        
    //==============================================================================================
    void compute_SZ_signal_moments(double (*Ne_f)(double lx, double ly, double lz), 
                                   double (*Te_f)(double lx, double ly, double lz), 
                                   double (*betac_f)(double lx, double ly, double lz), 
                                   double (*muc_f)(double lx, double ly, double lz), 
                                   double lx0, double ly0, double lsc, double lang,
                                   vector<double> &mv, 
                                   vector<double> &xv, 
                                   vector<double> &DIv,
                                   vector<double> &TSZetc,
                                   int variable, double lres);

    //==============================================================================================
    //
    // Compute the SZ signal using temperature-velocity moments, but assuming that the variance of
    // the temperature and velocity field along the line-of-sight is small. In this case the SZ 
    // signal is related to S_iso(tau, TeSZ) and its derivatives with respect to Te. The functions 
    // and parameters are similar to those above. The value of kmax<=5 determines how many 
    // derivatives are used to approximate the SZ signal. 
    //
    //==============================================================================================
    void compute_SZ_signal_moments_smooth(double (*Ne_f)(double l), 
                                          double (*Te_f)(double l), 
                                          double (*betac_f)(double l), 
                                          double (*muc_f)(double l), 
                                          double zsc,
                                          vector<double> &omegav, 
                                          vector<double> &xv, 
                                          vector<double> &DIv,
                                          vector<double> &TSZetc, 
                                          int kmax);

    //==============================================================================================
    void compute_SZ_signal_moments_smooth(const vector<double> &omegav,
                                          const vector<double> &TSZetc,
                                          vector<double> &xv, 
                                          vector<double> &DIv, 
                                          int kmax);

    //==============================================================================================
    void show_profile_slice(double (*f)(double l), double lsc, int np);
    void export_profile_slice(string fname, double (*f)(double l), double lsc, int np);
};*/

#endif

//==================================================================================================
//==================================================================================================
