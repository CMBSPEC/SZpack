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
#include "Integration_routines.h"
#include "nPl_derivatives.h"
#include "SZ_moment_method.h"
#include "SZ_asymptotic.h"
#include "SZ_CNSN_basis.opt.h"
#include "SZpack.h"

//==================================================================================================
//
// namespaces
//
//==================================================================================================
using namespace std;

//==================================================================================================
/*double G_func(double x)
{ return x*exp(-x)/(1.0-exp(-x))/(1.0-exp(-x)); }

double Q_func(double x)
{ return x*exp(-x)/(1.0-exp(-x))/(1.0-exp(-x))*x*(1.0+exp(-x))/(1.0-exp(-x)); }

double M_func(double x)
{ return x*exp(-x)/(1.0-exp(-x))/(1.0-exp(-x))*(x*(1.0+exp(-x))/(1.0-exp(-x))-3.0); }

//==================================================================================================
void sanity_check(int &order_T_low, int &order_T_high, int &order_b)
{
    order_T_low=(int)min(order_T_low, 10);
    order_T_high=(int)min(order_T_high, 20);
    order_b=(int)min(order_b, 2);
    
    return;
}

//==================================================================================================
//
// setup the moment matrix. This function has to be called prior to the computations. The SZ signal
// can then be computed in the form S = M * m 
//
//==================================================================================================
void SZ_moment_method::setup_SZ_moment_matrix(int order_T_low, int order_T_high, int order_b)
{
    sanity_check(order_T_low, order_T_high, order_b);
    
    //==============================================================================================
    // save info
    //==============================================================================================
    index_high   =1+order_T_low;
    nmom_high_one=1+order_T_high;
    
    if(order_b>0)
    { 
        index_high   +=2+order_T_low;                     // asymptotic expansion part
        nmom_high_one+=1+order_T_high;                    // CNSN expansion part

        if(order_b>1)
        {
            index_high   +=2*(2+order_T_low);             // asymptotic expansion part
            nmom_high_one+=2*(1+order_T_high);            // CNSN expansion part
        }
    }
    tot_moments=index_high+(nTregs-1)*nmom_high_one;
    
    //==============================================================================================
    // fill matrix
    //==============================================================================================
    vector<double> Y(order_T_low+1);
    M.resize((tot_moments+1)*(nx+1), 0.0);
    int ind;
    
    //TODO: This function needs fixing
    print_error("This function is not working at the moment");
    for(int n=0; n<nx; n++)
    {
        //==========================================================================================
        // description of low temperature gas
        //==========================================================================================
        compute_Y(x[n], Y);
        double x3=pow(x[n], 3);
        ind=n*tot_moments;
        
        for(int k=0; k<=order_T_low; k++)
            M[k+ind]=x3*Y[k];
        ind+=1+order_T_low;
        
        if(order_b>0)
        {
            M[ind++]=x3*G_func(x[n]);

            compute_D_CMB(x[n], Y);
            
            for(int k=0; k<=order_T_low; k++)
                M[k+ind]=x3*Y[k];
            ind+=1+order_T_low;

            if(order_b>1)
            {
                M[ind++]=x3*M_func(x[n])/3.0;
                
                compute_M_CMB(x[n], Y);
                
                for(int k=0; k<=order_T_low; k++)
                    M[k+ind]=x3*Y[k];
                ind+=1+order_T_low;
                
                M[ind++]=x3*11.0/30.0*Q_func(x[n]);
                
                compute_Q_CMB(x[n], Y);
                
                for(int k=0; k<=order_T_low; k++)
                    M[k+ind]=x3*Y[k];
                ind+=1+order_T_low;
            }
        }

        //==========================================================================================
        // description of high temperature gas
        //==========================================================================================
        if(nTregs>1)
        {
            vector<double> Yb(order_T_high+1);

            for(int r=1; r<nTregs; r++)
            {
                int ri=region_indices[r];
                compute_Y_CNSN(x[n], ri, Yb);
            
                for(int k=0; k<=order_T_high; k++)
                    M[k+ind]=x3*Yb[k];
                ind+=order_T_high+1;
                
                if(order_b>0)
                {
                    compute_D_CNSN_CMB(x[n], ri, Yb);
                    
                    for(int k=0; k<=order_T_high; k++)
                        M[k+ind]=x3*Yb[k];
                    ind+=order_T_high+1;
                    
                    if(order_b>1)
                    {
                        compute_M_CNSN_CMB(x[n], ri, Yb);
                        
                        for(int k=0; k<=order_T_high; k++)
                            M[k+ind]=x3*Yb[k];
                        ind+=order_T_high+1;
                        
                        compute_Q_CNSN_CMB(x[n], ri, Yb);
                        
                        for(int k=0; k<=order_T_high; k++)
                            M[k+ind]=x3*Yb[k];
                        ind+=order_T_high+1;
                    }
                }
            }
        }
    }
        

    //==============================================================================================
    // fill matrix for transformation from N(The) (The-The0)^k --> N(The) The^k
    //==============================================================================================
    T.resize(tot_moments*tot_moments, 0.0);
    
    for(int m=0; m<index_high; m++) 
        T[m*(1+tot_moments)]=1.0;
    
    for(int m=0; m<=order_T_high && nTregs>1; m++) 
    {
        for(int k=0; k<=m; k++) 
        {
            ind=index_high;
            
            for(int r=1; r<nTregs; r++)
            {
                double Tij=Binomial_coeff(m, k)*pow(-Get_The_ref(r), m-k);

                T[ind+k+(ind+m)*tot_moments]=Tij;
                ind+=order_T_high+1;
                
                if(order_b>0)
                {
                    T[ind+k+(ind+m)*tot_moments]=Tij;
                    ind+=order_T_high+1;
                    
                    if(order_b>1)
                    {
                        T[ind+k+(ind+m)*tot_moments]=Tij;
                        ind+=order_T_high+1;
                        T[ind+k+(ind+m)*tot_moments]=Tij;
                        ind+=order_T_high+1;
                    }
                }
            }
        }
    }
    
    //==============================================================================================
    // compute transformed matrix MT = M*T
    //==============================================================================================
    MT.resize((tot_moments+1)*(nx+1), 0.0);

    for(int n=0; n<nx; n++)
    {
        for(int m=0; m<tot_moments; m++)
        {
            double r=0.0;
            for(int k=tot_moments-1; k>=0; k--) 
                r+=M[k+n*tot_moments]
                  *T[m+k*tot_moments];
                
            MT[m+n*tot_moments]=r;
        }
    }
        
//    show_moment_matrix_M();
//    show_moment_matrix_MT();
//    show_conversion_matrix_T();     
    
    return;
}



//==================================================================================================
//
// show data of matrices
//
//==================================================================================================
SZ_moment_method::SZ_moment_method(){}

//==================================================================================================
SZ_moment_method::SZ_moment_method(const vector<double> &xcmb, 
                                   int kmaxv, int accuracy_levelv, int order_bv)
{
    x=xcmb;
    nx=x.size();
    kmax=kmaxv;
    accuracy_level=accuracy_levelv;
    order_b=order_bv;
    
    Te_pivots=Get_temperature_pivots (kmax, accuracy_level);
    Te_limits=Get_temperature_regions(kmax, accuracy_level);
    region_indices=Get_region_indices(kmax, accuracy_level);

    nTregs=Te_limits.size();
    Te_max=Te_limits.back();
    
    setup_SZ_moment_matrix(kmax, kmax, order_b);
    
    show_Temperature_limits();
    show_Temperature_pivots();
}

//==================================================================================================
SZ_moment_method::SZ_moment_method(const vector<double> &xcmb, 
                                   double Te_maxv, int accuracy_levelv, int order_bv)
{
    x=xcmb;
    nx=x.size();
    accuracy_level=accuracy_levelv;
    order_b=order_bv;
    Te_max=Te_maxv;
    
    determine_optimal_kmax(accuracy_level, Te_max, kmax, nTregs);
    Te_pivots=Get_temperature_pivots (kmax, accuracy_level);
    Te_limits=Get_temperature_regions(kmax, accuracy_level);
    region_indices=Get_region_indices(kmax, accuracy_level);

    setup_SZ_moment_matrix(kmax, kmax, order_b);
    
    show_Temperature_limits();
    show_Temperature_pivots();
}

//==================================================================================================
SZ_moment_method::~SZ_moment_method(){}

//==================================================================================================
//
// simple functions 
//
//==================================================================================================
int SZ_moment_method::Get_index_of_Te_region(double Te)  // if Te > Temax --> -10 is returned
{
    if(Te>Te_limits[nTregs-1]) return -10;
    if(Te<Te_limits[0])        return 0;

    int reg=0;
    do reg++;
    while(Te>Te_limits[reg]);
    
    return reg;        
}

int SZ_moment_method::Get_start_index_of_moment_set(int reg)
{
    if(reg==0) return 0;
    return index_high+(reg-1)*Get_number_of_moments_for_one_high_region();
}

double SZ_moment_method::Get_Te_ref(int reg){ return Te_pivots[reg]*const_me; }
double SZ_moment_method::Get_The_ref(int reg){ return Te_pivots[reg]; }

//==================================================================================================
//
// show moment data
//
//==================================================================================================
void SZ_moment_method::show_Temperature_limits()
{
    cout << "\n SZ_moment_method::show_Temperature_limits: The temperature limits are:" << endl;
    cout << " " << 0 << "\t - \t" << Te_limits[0] << "\tkeV (Asymptotic expansion) " << endl;    

    for(int n=1; n<nTregs; n++) 
        cout << " " << Te_limits[n-1] << "\t - \t" << Te_limits[n] << "\tkeV (CNSN basis)" << endl; 
    
    cout << " Total number of Te-regions " << nTregs << endl;
    
    return;
}

//==================================================================================================
void SZ_moment_method::show_Temperature_pivots()
{
    cout << "\n SZ_moment_method::show_Temperature_pivots: The pivot temperatures are:" << endl;
    cout << " Region " << 0 << ": " << Te_pivots[0]*const_me << "\tkeV (Asymptotic expansion)" << endl;    

    for(int n=1; n<nTregs; n++) 
        cout << " Region " << n << ": " << Te_pivots[n]*const_me << "\tkeV (CNSN basis)" << endl;    

    return;
}

//==================================================================================================
template<typename somestream>
void SZ_moment_method::output_moment_vector(somestream &output, const vector<double> &mv)
{
    output.precision(6);

    output << "\n " << setfill('=') << setw(95) << "=" << endl;  

    //==============================================================================================
    output << "\n asymptotic expansion part" << endl;
   
    output << " " << setfill('-') << setw(42) << "-" 
           << " Region 0 " << setw(43) << "-" << endl;
    
    int dn=kmax+1;
    for(int m=0; m<kmax+1; m++)
    {
        output << scientific << " y^("+int_to_string(m)+ ") = " << mv[m];

        if(order_b>0)
        {
            output << scientific << "\t b_1^("+int_to_string(m)+ ") = " << mv[dn+m];
            
            if(order_b>1)
            {
                output << scientific << "\t b_0^("+int_to_string(m)+ ") = " << mv[2*dn+1+m];
                output << scientific << "\t b_2^("+int_to_string(m)+ ") = " << mv[3*dn+2+m];
            }
        }

        output << endl;
    }   

    //==============================================================================================
    // pure kinematic terms
    //==============================================================================================
    if(order_b>0)
    {
        output << "                     " << scientific 
               << "\t b_1^("+int_to_string(dn)+ ") = " << mv[dn+dn];
        
        if(order_b>1)
        {
            output << scientific << "\t b_0^("+int_to_string(dn)+ ") = " << mv[2*dn+1+dn];
            output << scientific << "\t b_2^("+int_to_string(dn)+ ") = " << mv[3*dn+2+dn];
        }

        output << endl;
    } 
    output << " " << setfill('-') << setw(95) << "-" << endl;
        
    //==============================================================================================
    output << "\n CNSN expansion part";
    
    for(int r=1; r<nTregs; r++)
    {
        int is=index_high+(r-1)*nmom_high_one;
        
        output << "\n " << setfill('-') << setw(42) << "-" 
               << " Region " << int_to_string(r) << " " << setw(43) << "-" << endl;
        
        for(int m=0; m<kmax+1; m++)
        {
            output << scientific << " z^("+int_to_string(m)+ ") = " << mv[is+m];
            
            if(order_b>0)
            {
                output << scientific << "\t c_1^("+int_to_string(m)+ ") = " << mv[is+dn+m];
                
                if(order_b>1)
                {
                    output << scientific << "\t c_0^("+int_to_string(m)+ ") = " << mv[is+2*dn+m];
                    output << scientific << "\t c_2^("+int_to_string(m)+ ") = " << mv[is+3*dn+m];
                }
            }
            
            output << endl;
        }   
    }        
    
    output << "\n " << setfill('=') << setw(95) << "=" << endl; 
    
    return;
}

//==================================================================================================
void SZ_moment_method::show_moment_vector(const vector<double> &mv)
{
    output_moment_vector(cout, mv);
    return;
}

//==================================================================================================
//
// show data of matrices
//
//==================================================================================================
void SZ_moment_method::show_x_vector()
{
    for(int n=0; n<nx; n++) cout << n << " " << x[n] << endl;
    return; 
}

//==================================================================================================
void SZ_moment_method::show_moment_matrix(vector<double> &B)
{
    for(int n=0; n<nx; n++)
    {
        cout << x[n] << " || ";
        
        for(int m=0; m<tot_moments; m++)
        {
            if(m==index_high) cout << "|| ";
            cout << B[m+n*tot_moments] << " ";
        }
        
        cout << endl;
    }
    cout << endl;
    return;
}

void SZ_moment_method::show_moment_matrix_M()
{
    show_moment_matrix(M);
    return;
}

void SZ_moment_method::show_moment_matrix_MT()
{
    show_moment_matrix(MT);
    return;
}


//==================================================================================================
void SZ_moment_method::show_conversion_matrix(vector<double> &C)
{
    for(int n=0; n<tot_moments; n++)
    {
        for(int m=0; m<tot_moments; m++)
        {
            if(m==index_high) cout << "|| ";
            cout << C[m+n*tot_moments] << " ";
        }
        
        cout << endl;
    }
    cout << endl;
    return;
}

void SZ_moment_method::show_conversion_matrix_T()
{
    show_conversion_matrix(T);
    return;
}

//==================================================================================================
//
// header for matrix info
//
//==================================================================================================
void SZ_moment_method::print_header_for_moment_matrix(ofstream &ofile, string mess)
{   
    ofile << setfill('#') << setw(90) << "#" << endl << "#" << endl;  
    ofile << "# SZ moment matrix " << mess << endl;
    ofile << "# Number of frequency points: " << nx << " (# of rows)\n";
    ofile << "# velocity order            : " << order_b << endl;
    ofile << "# accuracy level            : " << accuracy_level << endl;
    ofile << "# kmax+1                    : " << kmax+1 << " (# of Te (!) moments per region)\n";
    ofile << "# Number of Te regions      : " << nTregs << " \n";
    ofile << "# total number of moments   : " << tot_moments << " (# of cols)\n";
    ofile << "#\n" << setfill('#') << setw(90) << "#" << endl; 
    
    return;
}

//==================================================================================================
void SZ_moment_method::export_x_vector(string fname)
{
    ofstream ofile(fname.c_str());
    ofile.precision(16);
    
    for(int n=0; n<nx; n++) 
        ofile << n << " " << x[n] << endl;
    
    ofile.close();
    return; 
}

//==================================================================================================
void SZ_moment_method::export_moment_vector(string fname, const vector<double> &mv)
{
    ofstream ofile(fname.c_str());
    
    output_moment_vector(ofile, mv);
    
    ofile.close();
    return; 
}

//==================================================================================================
void SZ_moment_method::export_moment_matrix(string fname, string form, vector<double> &B)
{
    ofstream ofile(fname.c_str());
    ofile.precision(16);
    
    print_header_for_moment_matrix(ofile, form);
    for(int n=0; n<nx; n++)
    {
        for(int m=0; m<tot_moments; m++)
            ofile << B[m+n*tot_moments] << " ";
        
        ofile << endl;
    }    

    ofile.close();
    return; 
}

void SZ_moment_method::export_moment_matrix_M(string fname)
{
    export_moment_matrix(fname, "M", M);
    return; 
}

void SZ_moment_method::export_moment_matrix_MT(string fname)
{
    export_moment_matrix(fname, "MT", MT);
    return; 
}

//==================================================================================================
//
// methods to compute the required moments for given profile functions 
//
//==================================================================================================
double (*ptr_Ne_func)(double);
double (*ptr_Te_func)(double);
double (*ptr_betac_func)(double);
double (*ptr_muc_func)(double);

//==================================================================================================
struct Te_region
{
    double lmin, lmax;
};

vector<vector<Te_region> > Te_zeros;

void SZ_moment_method::determine_Te_structure(double lres)
{
    Te_zeros.clear();
    Te_zeros.resize(nTregs+1);
    Te_region dum;
    
    int nres=floor(1.0/lres);
    
    double dx=2.0/(nres-1);
    
    double lmin=-1.0;
    double Te=ptr_Te_func(lmin);
    int ir=Get_index_of_Te_region(Te), irp;
    dum.lmin=lmin;
    
    for(int m=1; m<nres; m++)
    {
        lmin+=dx;
        Te=ptr_Te_func(lmin);
        irp=Get_index_of_Te_region(Te);
        
        if(ir!=irp)
        {
            dum.lmax=lmin;
            Te_zeros[ir].push_back(dum);
            dum.lmin=lmin;
            ir=irp;
        }        
    }
    dum.lmax=lmin;
    Te_zeros[ir].push_back(dum);

    dum.lmin=dum.lmax=0.0;
    for(int m=0; m<nTregs; m++) 
        if(Te_zeros[m].size()==0) Te_zeros[m].push_back(dum);
    
/*    for(int m=0; m<nTregs; m++)
    {
        for(int z=0; z<(int)Te_zeros[m].size(); z++)
            cout << " Region: " << m << " " << z 
                 << " l-range: " << Te_zeros[m][z].lmin 
                 << " " << Te_zeros[m][z].lmax << endl;
    }
    wait_f_r();
*/   /*
    return;
}


//==================================================================================================
double compute_Integral(double lmin, double lmax, 
                        double (*func)(double l, void *userdata), 
                        void *userdata, double Int_eps)
{ 
    double r=Integrate_using_Patterson_adaptive(max(-1.0, lmin), min(-0.1, lmax), 
                                                Int_eps, 1.0e-300, func, userdata);
    r+=Integrate_using_Patterson_adaptive(max(-0.1, lmin), min(-0.01, lmax), 
                                          Int_eps, 1.0e-300, func, userdata);
    r+=Integrate_using_Patterson_adaptive(max(-0.01, lmin), min(-0.001, lmax), 
                                          Int_eps, 1.0e-300, func, userdata);
    r+=Integrate_using_Patterson_adaptive(max(-0.001, lmin), min(-0.0001, lmax), 
                                          Int_eps, 1.0e-300, func, userdata);
    r+=Integrate_using_Patterson_adaptive(max(-0.0001, lmin), min(-0.00001, lmax), 
                                          Int_eps, 1.0e-300, func, userdata);
    r+=Integrate_using_Patterson_adaptive(max(-0.00001, lmin), min(0.00001, lmax), 
                                          Int_eps, 1.0e-300, func, userdata);
    r+=Integrate_using_Patterson_adaptive(max(0.00001, lmin), min(0.0001, lmax), 
                                          Int_eps, 1.0e-300, func, userdata);
    r+=Integrate_using_Patterson_adaptive(max(0.0001, lmin), min(0.001, lmax), 
                                          Int_eps, 1.0e-300, func, userdata);
    r+=Integrate_using_Patterson_adaptive(max(0.001, lmin), min(0.01, lmax), 
                                          Int_eps, 1.0e-300, func, userdata);
    r+=Integrate_using_Patterson_adaptive(max(0.01, lmin), min(0.1, lmax), 
                                          Int_eps, 1.0e-300, func, userdata);
    r+=Integrate_using_Patterson_adaptive(max(0.1, lmin), min(1.0, lmax), 
                                          Int_eps, 1.0e-300, func, userdata);    
    return r;
}

//==================================================================================================
//
// functions to compute y, z, b, and c integrals according to CSNN 2012
//
//==================================================================================================
double (*Norm_func)(double The);
double Norm_func_1(double The){ return 1.0; }
double Norm_func_N(double The){ return N_func_CNSN(The); }
double The_ref_global;

double g_eps_Int=1.0e-6;

//==================================================================================================
double func_dtau(double l, void *userdata)
{ return ptr_Ne_func(l); }

double compute_tau(double zsc)
{ 
    void *userdata=NULL;
    return zsc*const_sigT*compute_Integral(-1.0, 1.0, func_dtau, userdata, g_eps_Int);
}

//==================================================================================================
double func_dySZ(double l, void *userdata)
{ return ptr_Ne_func(l)*ptr_Te_func(l)/const_me; }

double compute_ySZ(double zsc)
{ 
    void *userdata=NULL;
    return zsc*const_sigT*compute_Integral(-1.0, 1.0, func_dySZ, userdata, g_eps_Int);
}

double func_dyk(double l, void *userdata)
{ 
    double *d=(double *)userdata;
    double Te=ptr_Te_func(l);
    return ptr_Ne_func(l)*pow(Te/const_me, d[0]+1); 
}
                              
double compute_yk(double zsc, int k)
{ 
    double d[]={k};
    void *userdata=(void *)d;
    return zsc*const_sigT*compute_Integral(-1.0, 1.0, func_dyk, userdata, g_eps_Int);
}

//==================================================================================================
// low temperature moments
//==================================================================================================
double func_dy_k_cond(double l, void *userdata)
{ 
    double *d=(double *)userdata;
    double Te=ptr_Te_func(l);
    double D=Te/const_me-The_ref_global;
    return Norm_func(Te/const_me)*ptr_Ne_func(l)*pow(D, d[0]+1); 
}

double compute_y_k_cond(int k, int ir, double zsc)
{ 
    double d[]={k};
    void *userdata=(void *)d;
    double r=0.0;
    
    for(int z=0; z<(int)Te_zeros[ir].size(); z++)
        
        r+=compute_Integral(Te_zeros[ir][z].lmin, Te_zeros[ir][z].lmax, 
                            func_dy_k_cond, userdata, g_eps_Int);
    
    return zsc*const_sigT*r;
}

//==================================================================================================
double func_db_k_0_cond(double l, void *userdata)
{ 
    double *d=(double *)userdata;
    double Te=ptr_Te_func(l);
    double b=ptr_betac_func(l);
    double D=Te/const_me-The_ref_global;
    return Norm_func(Te/const_me)*ptr_Ne_func(l)*b*b*pow(D, d[0]); 
}

double compute_b_k_0_cond(int k, int ir, double zsc)
{ 
    double d[]={k};
    void *userdata=(void *)d;
    double r=0.0;
    
    for(int z=0; z<(int)Te_zeros[ir].size(); z++)
        
        r+=compute_Integral(Te_zeros[ir][z].lmin, Te_zeros[ir][z].lmax, 
                            func_db_k_0_cond, userdata, g_eps_Int);
    
    return zsc*const_sigT*r;
}

double compute_b_0(double zsc)
{ 
    double d[]={0};
    void *userdata=(void *)d;
    return zsc*const_sigT*compute_Integral(-1.0, 1.0, func_db_k_0_cond, userdata, g_eps_Int);
}

//==================================================================================================
double func_db_k_1_cond(double l, void *userdata)
{ 
    double *d=(double *)userdata;
    double Te=ptr_Te_func(l);
    double b=ptr_betac_func(l);
    double muc=ptr_muc_func(l);
    double D=Te/const_me-The_ref_global;
    return Norm_func(Te/const_me)*ptr_Ne_func(l)*b*muc*pow(D, d[0]); 
}

double compute_b_k_1_cond(int k, int ir, double zsc)
{ 
    double d[]={k};
    void *userdata=(void *)d;
    double r=0.0;
    
    for(int z=0; z<(int)Te_zeros[ir].size(); z++)
        
        r+=compute_Integral(Te_zeros[ir][z].lmin, Te_zeros[ir][z].lmax, 
                            func_db_k_1_cond, userdata, g_eps_Int);
    
    return zsc*const_sigT*r;
}

double compute_b_1(double zsc)
{ 
    double d[]={0};
    void *userdata=(void *)d;
    return zsc*const_sigT*compute_Integral(-1.0, 1.0, func_db_k_1_cond, userdata, g_eps_Int);
}

//==================================================================================================
double func_db_k_2_cond(double l, void *userdata)
{ 
    double *d=(double *)userdata;
    double Te=ptr_Te_func(l);
    double b=ptr_betac_func(l);
    double muc=ptr_muc_func(l);
    double D=Te/const_me-The_ref_global;
    return Norm_func(Te/const_me)*ptr_Ne_func(l)*b*b*(1.5*muc*muc-0.5)*pow(D, d[0]); 
}

double compute_b_k_2_cond(int k, int ir, double zsc)
{ 
    double d[]={k};
    void *userdata=(void *)d;
    double r=0.0;
    
    for(int z=0; z<(int)Te_zeros[ir].size(); z++)
        
        r+=compute_Integral(Te_zeros[ir][z].lmin, Te_zeros[ir][z].lmax, 
                            func_db_k_2_cond, userdata, g_eps_Int);
    
    return zsc*const_sigT*r;
}

double compute_b_2(double zsc)
{ 
    double d[]={0};
    void *userdata=(void *)d;
    return zsc*const_sigT*compute_Integral(-1.0, 1.0, func_db_k_2_cond, userdata, g_eps_Int);
}


//==================================================================================================
//
// additional average over x and y
// 
//==================================================================================================
double (*ptr_Ne_func_3D)(double, double, double);
double (*ptr_Te_func_3D)(double, double, double);
double (*ptr_betac_func_3D)(double, double, double);
double (*ptr_muc_func_3D)(double, double, double);
double (*ptr_Integrand_1D)(double lsc);
double (*ptr_Integrand_1D_yk)(double lsc, int k);


double lx_glob, ly_glob, lsc_glob, lang_glob;
double Ne_func_1D(double lz){ return ptr_Ne_func_3D(lx_glob, ly_glob, lz); }
double Te_func_1D(double lz){ return ptr_Te_func_3D(lx_glob, ly_glob, lz); }
double betac_func_1D(double lz){ return ptr_betac_func_3D(lx_glob, ly_glob, lz); }
double muc_func_1D(double lz){ return ptr_muc_func_3D(lx_glob, ly_glob, lz); }

int k_glob;
double Integrand_1D_yk(double lsc){ return ptr_Integrand_1D_yk(lsc, k_glob); }


//==================================================================================================
double beam_function(double Dlx, double Dly, double lang)
{
    double l2=Dlx*Dlx+Dly*Dly;
    return exp(-l2/2.0/lang/lang)/TWOPI/lang/lang;
}


//==================================================================================================
double compute_Int_xy(double ly, void *userdata)
{ 
    ly_glob=ly;
    
    double *d=(double *)userdata;
    double fac=beam_function(lx_glob-d[0], ly_glob-d[1], d[3]);
    
    return fac*ptr_Integrand_1D(d[2]); 
}

double compute_Int_x(double lx, void *userdata)
{ 
    lx_glob=lx;
    double *d=(double *)userdata;
    double a=d[1]-6*d[3], b=d[1]+6*d[3];
    
    return Integrate_using_Patterson_adaptive(a, b, 1.0e-4, 1.0e-300, compute_Int_xy, userdata);
}

//==================================================================================================
double compute_Int_3D(double lx0, double ly0, double lsc, double lang, double (*f)(double lsc))
{ 
    ptr_Integrand_1D=f;
    double d[]={lx0, ly0, lsc, lang};
    double a=d[0]-6*d[3], b=d[0]+6*d[3];
    void *userdata=(void *)d;
    
    return Integrate_using_Patterson_adaptive(a, b, 2.0e-4, 1.0e-300, compute_Int_x, userdata);
}

double compute_Int_3D_yk(double lx0, double ly0, double lsc, double lang, int k,
                         double (*f)(double lsc, int k))
{ 
    ptr_Integrand_1D_yk=f;
    k_glob=k;
    return compute_Int_3D(lx0, ly0, lsc, lang, Integrand_1D_yk);
}

//==================================================================================================
//
// some output functions
//
//==================================================================================================
template<typename somestream>
void SZ_moment_method::output_profile_slice(somestream &output, 
                                            double (*f)(double l), 
                                            double lsc, int np)
{
    output.precision(16);
    vector<double> lz(np); 
    init_xarr(1.0e-5, 1, &lz[0], np, 1, 0); 
    
    for(int k=np-1; k>=0; k--) output << " " << -lsc*lz[k] << " " << f(-lz[k]) << endl;
    output << " " << 0.0 << " " << f(0.0) << endl;
    for(int k=0; k<np; k++) output << " " << lsc*lz[k] << " " << f(lz[k]) << endl;
    
    return;
}

void SZ_moment_method::show_profile_slice(double (*f)(double l), double lsc, int np)
{
    output_profile_slice(cout, f, lsc, np);
    return;
}

void SZ_moment_method::export_profile_slice(string fname, double (*f)(double l), double lsc, int np)
{
    ofstream ofile(fname.c_str());
    output_profile_slice(ofile, f, lsc, np);
    ofile.close();
    return;
}

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
void SZ_moment_method::compute_SZ_signal_moments(const vector<double> &mv, 
                                                 vector<double> &xv, 
                                                 vector<double> &DIv,
                                                 int variable)
{
    if(tot_moments!=(int)mv.size())
    { cerr << "please check moment vector or moment matrix " << endl; exit(0); }
    
    DIv.resize(nx);
    xv=x;
    
    const double *Mptr=NULL;
    
    //==============================================================================================
    // check the variable
    //==============================================================================================
    if(variable==0) Mptr=&M[0];
    else if(variable==1) Mptr=&MT[0];
    else
    { 
        cerr << " compute_SZ_signal_moments :: This variable does not exist. Exiting." << endl; 
        exit(0); 
    }
    
    //==============================================================================================
    // S = F * m like in Chluba et. al, 2012
    //==============================================================================================
    for(int n=0; n<nx; n++)
    {
        DIv[n]=0.0;
        
        for(int m=0; m<tot_moments; m++) DIv[n]+=Mptr[m+n*tot_moments]*mv[m];
    }
    
    Mptr=NULL;
    
    return;
}

//==================================================================================================
// 
// compute moments over the cluster profiles.
//
//==================================================================================================
void SZ_moment_method::compute_cluster_moments(double zsc, vector<double> &mv, int Te_order, 
                                               int variable, double lres)
{
    double yk, bk;
    mv.clear();

    //==============================================================================================
    // brute force temperature structure computation; This has to be done to determine the intervals
    // in l=z/zsc that correspond to the different temperature regions. Temperature structures that 
    // are smaller than zsc/nres will not be resolved.
    //==============================================================================================
    determine_Te_structure(lres);
    
    //==============================================================================================
    // low temperature moments
    //==============================================================================================
    Norm_func=Norm_func_1;
    The_ref_global=0.0;
    
    for(int k=0; k<=Te_order; k++)
    {
        yk=compute_y_k_cond(k, 0, zsc);
        mv.push_back(yk);
    }
    
    //==============================================================================================
    if(order_b>0)
    {
        bk=compute_b_1(zsc);
        mv.push_back(bk);
        
        for(int k=0; k<=Te_order; k++)
        {
            bk=compute_b_k_1_cond(k+1, 0, zsc);
            mv.push_back(bk);
        }
        
        if(order_b>1)
        {
            bk=compute_b_0(zsc);
            mv.push_back(bk);
            
            for(int k=0; k<=Te_order; k++)
            {
                bk=compute_b_k_0_cond(k+1, 0, zsc);
                mv.push_back(bk);
            }
            
            bk=compute_b_2(zsc);
            mv.push_back(bk);
            
            for(int k=0; k<=Te_order; k++)
            {
                bk=compute_b_k_2_cond(k+1, 0, zsc);
                mv.push_back(bk);
            }
        }
    }
    
    //==============================================================================================
    // high temperature moments
    //==============================================================================================
    if(nTregs>1)
    {
        Norm_func=Norm_func_N;
        
        for(int r=1; r<nTregs; r++)
        {
            if(variable==0) The_ref_global=Get_The_ref(r);
            else The_ref_global=0.0;
            
            for(int k=0; k<=Te_order; k++)
            {
                // k-1 because of difference between y^(k)<-->z^(k) variable
                yk=compute_y_k_cond(k-1, r, zsc); 
                mv.push_back(yk);
            }
            
            //======================================================================================
            if(order_b>0)
            {
                for(int k=0; k<=Te_order; k++)
                {
                    bk=compute_b_k_1_cond(k, r, zsc);
                    mv.push_back(bk);
                }
                
                if(order_b>1)
                {
                    for(int k=0; k<=Te_order; k++)
                    {
                        bk=compute_b_k_0_cond(k, r, zsc);
                        mv.push_back(bk);
                    }
                    
                    for(int k=0; k<=Te_order; k++)
                    {
                        bk=compute_b_k_2_cond(k, r, zsc);
                        mv.push_back(bk);
                    }
                }
            }
        }
    }

    return;
}

//==================================================================================================
// 
// compute moments over the cluster profiles.
//
//==================================================================================================
void SZ_moment_method::compute_cluster_moments_3D(double lx0, double ly0, double lsc, double lang, 
                                                  vector<double> &mv, int Te_order, 
                                                  int variable, double lres)
{
//    double yk, bk;
    mv.clear();
    
    vector<double> dum_mv;
    
    int np=10;
    
    for(int c=0; c<np; c++)
    {
//        lx_glob=la[c];
        lx_glob=1.0/lsc;
        
        for(int r=0; r<np; r++)
        {
//            ly_glob=la[r];
            ly_glob=1.0/lsc;
            
            compute_cluster_moments(lsc, dum_mv, Te_order, variable, lres);
        }
    }
    
    return;
}

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
void SZ_moment_method::compute_SZ_signal_moments(double (*Ne_f)(double l), 
                                                 double (*Te_f)(double l), 
                                                 double (*betac_f)(double l), 
                                                 double (*muc_f)(double l), 
                                                 double zsc,
                                                 vector<double> &mv, 
                                                 vector<double> &xv, 
                                                 vector<double> &DIv,
                                                 vector<double> &TSZetc,
                                                 int variable, double lres)
{
    //==============================================================================================
    // clean plates & simple setup
    //==============================================================================================
    mv.clear();
    xv.clear();
    DIv.clear();
    TSZetc.clear();

    // set pointers to functions, so that profiles can be used
    ptr_Ne_func=Ne_f;
    ptr_Te_func=Te_f;
    ptr_betac_func=betac_f;
    ptr_muc_func=muc_f;

    // simple calculation of standard variables
    double tau=compute_tau(zsc);
    double ySZ=compute_ySZ(zsc);
    double TSZ=ySZ*const_me/tau;
    
    // return values
    TSZetc.push_back(tau);
    TSZetc.push_back(ySZ);
    TSZetc.push_back(TSZ);
    
    //==============================================================================================
    // compute the signal
    //==============================================================================================
    compute_cluster_moments(zsc, mv, kmax, variable, lres);
    compute_SZ_signal_moments(mv, xv, DIv, variable);
    
    return;
}

//==================================================================================================
void SZ_moment_method::compute_SZ_signal_moments(double (*Ne_f)(double lx, double ly, double lz), 
                                                 double (*Te_f)(double lx, double ly, double lz), 
                                                 double (*betac_f)(double lx, double ly, double lz), 
                                                 double (*muc_f)(double lx, double ly, double lz), 
                                                 double lx0, double ly0, double lsc, double lang,
                                                 vector<double> &mv, 
                                                 vector<double> &xv, 
                                                 vector<double> &DIv,
                                                 vector<double> &TSZetc,
                                                 int variable, double lres)
{
    
    //==============================================================================================
    // clean plates & simple setup
    //==============================================================================================
    mv.clear();
    xv.clear();
    DIv.clear();
    TSZetc.clear();
    
    // set pointers to functions, so that profiles can be used
    ptr_Ne_func_3D=Ne_f;
    ptr_Te_func_3D=Te_f;
    ptr_betac_func_3D=betac_f;
    ptr_muc_func_3D=muc_f;

    ptr_Ne_func=Ne_func_1D;
    ptr_Te_func=Te_func_1D;
    ptr_betac_func=betac_func_1D;
    ptr_muc_func=muc_func_1D;
    
    // simple calculation of standard variables
    double tau=compute_Int_3D(lx0, ly0, lsc, lang, compute_tau);
    double ySZ=compute_Int_3D(lx0, ly0, lsc, lang, compute_ySZ);
    double TSZ=ySZ*const_me/tau;
    
    // return values
    TSZetc.push_back(tau);
    TSZetc.push_back(ySZ);
    TSZetc.push_back(TSZ);
    
    //==============================================================================================
    // compute the signal
    //==============================================================================================
    compute_cluster_moments_3D(lx0, ly0, lsc, lang, mv, kmax, variable, lres);
    wait_f_r("STOP");
    compute_SZ_signal_moments(mv, xv, DIv, variable);
    
    return;
}

//==============================================================================================
//
// Compute the SZ signal using temperature-velocity moments, but assuming that the variance of
// the temperature and velocity field along the line-of-sight is small. In this case the SZ 
// signal is related to S_iso(tau, TeSZ) and its derivatives with respect to Te. The functions 
// and parameters are similar to those above.
//
//==============================================================================================
void SZ_moment_method::compute_SZ_signal_moments_smooth(double (*Ne_f)(double l), 
                                                        double (*Te_f)(double l), 
                                                        double (*betac_f)(double l), 
                                                        double (*muc_f)(double l), 
                                                        double zsc,
                                                        vector<double> &omegav, 
                                                        vector<double> &xv, 
                                                        vector<double> &DIv,
                                                        vector<double> &TSZetc, 
                                                        int kmax)
{
    //==============================================================================================
    // clean plates & simple setup
    //==============================================================================================
    omegav.clear();
    xv.clear();
    DIv.clear();
    TSZetc.clear();
    
    // set pointers to functions, so that profiles can be used
    ptr_Ne_func=Ne_f;
    ptr_Te_func=Te_f;
    ptr_betac_func=betac_f;
    ptr_muc_func=muc_f;
    
    // simple calculation of standard variables
    double tau=compute_tau(zsc);
    double ySZ=compute_ySZ(zsc);
    double TSZ=ySZ*const_me/tau;
    
    // return values
    TSZetc.push_back(tau);
    TSZetc.push_back(ySZ);
    TSZetc.push_back(TSZ);
    
    //==============================================================================================
    // compute omega(k)
    //============================================================================================== 
    vector<double> rhok(kmax+1);
    omegav.resize(kmax+1);
    omegav[0]=0.0;
    rhok[0]=1.0;

    for(int k=1; k<=kmax; k++) rhok[k]=compute_yk(zsc, k)/ySZ/pow(TSZ/const_me, k);
    
    for(int m=1; m<=kmax; m++)
    {
        omegav[m]=pow(-1.0, m+1);
        
        for(int l=1; l<=m+1; l++)
            omegav[m]+=Binomial_coeff(m+1, l)*pow(-1.0, m+1-l)*rhok[l-1];
    }
    
    //==============================================================================================
    // compute the signal
    //==============================================================================================    
    compute_SZ_signal_moments_smooth(omegav, TSZetc, xv, DIv, kmax);
    
    return;
}    

//==================================================================================================
void SZ_moment_method::compute_SZ_signal_moments_smooth(const vector<double> &omegav,
                                                        const vector<double> &TSZetc,
                                                        vector<double> &xv, 
                                                        vector<double> &DIv, 
                                                        int kmax)
{
    DIv.resize(nx);
    xv=x;
        
    double tau=TSZetc[0];
    double TeSZ=TSZetc[2];
    double betac=0.0;
    double muc=1.0;
    double betao=0.0;
    double muo=0.0;
    
    double eps=0.005;
    
    //==============================================================================================
    // S = F * m like in Chluba et. al, 2012
    //==============================================================================================
    double d0, dp, dm, dp2, dm2;
    double d2Sd2T, d3Sd3T, d4Sd4T;
    
    for(int n=0; n<nx; n++)
    {
        d0=compute_SZ_signal_combo(x[n], tau, TeSZ, betac, muc, betao, muo);
        
        DIv[n]=d0;
        
        if(kmax>0)
        {
            dp=compute_SZ_signal_combo(x[n], tau, TeSZ*(1.0+eps), betac, muc, betao, muo);
            dm=compute_SZ_signal_combo(x[n], tau, TeSZ*(1.0-eps), betac, muc, betao, muo);
            d2Sd2T=(dp-2.0*d0+dm)/pow(eps, 2) / (2.0);
            
            DIv[n]+=omegav[1]*d2Sd2T;
            
            if(kmax>1)
            {
                dp2=compute_SZ_signal_combo(x[n], tau, TeSZ*(1.0+2*eps), betac, muc, betao, muo);
                dm2=compute_SZ_signal_combo(x[n], tau, TeSZ*(1.0-2*eps), betac, muc, betao, muo);
                d3Sd3T=(dp2-2.0*dp+2.0*dm-dm2)/2.0/pow(eps, 3) / (6.0);
                
                DIv[n]+=omegav[2]*d3Sd3T;

                if(kmax>2)
                {
                    d4Sd4T=(dp2-4.0*dp+6.0*d0-4.0*dm+dm2)/pow(eps, 4) / (24.0);

                    DIv[n]+=omegav[3]*d4Sd4T;
                }
            }
        }
        
        DIv[n]*=pow(x[n], 3);
    }
            
    return;
}

*/

//==================================================================================================
//==================================================================================================
