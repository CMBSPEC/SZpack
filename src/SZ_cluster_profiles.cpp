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
#include "SZ_cluster_profiles.h"

//==================================================================================================
//
// namespaces
//
//==================================================================================================
using namespace std;

//==============================================================================================
//
// Constructor & Destructor
//
//==============================================================================================
SZ_cluster_profiles::SZ_cluster_profiles(){}
SZ_cluster_profiles::~SZ_cluster_profiles(){}

//==============================================================================================
// initialize Cluster profile
//==============================================================================================
SZ_cluster_profiles::SZ_cluster_profiles(Cluster_param_Ne &pNe, Cluster_param_Te &pTe, string lab)
{
    SZpNe=pNe;
    SZpTe=pTe;  
    label=lab;
}

//==============================================================================================
// output cluster parameters
//==============================================================================================
template< typename somestream >
void SZ_cluster_profiles::output_cluster_parameters(somestream &output)
{
    output.precision(6);
    
    output << "\n " << setfill('=') << setw(90) << "=" << endl;  
    
    //==============================================================================================
    output << "\n Ne parameters " << endl;
    output << " n0 = " << SZpNe.n0 << " x 10^-3 cm^-3  " 
           << " rc = " << SZpNe.rc << " kpc  "
           << " rs = " << SZpNe.rs << " kpc" << endl;
    output << " alpha = " << SZpNe.alpha << "  " 
           << " beta = " << SZpNe.beta << "  "
           << " eps = " << SZpNe.eps << endl;
    output << " n02 = " << SZpNe.n02 << " x 10^-1 cm^-3  " 
           << " rc2 = " << SZpNe.rc2 << " kpc  "
           << " beta2 = " << SZpNe.beta2 << endl;
    
    //==============================================================================================
    output << "\n Te parameters " << endl;
    output << " T0 = " << SZpTe.T0 << " keV  " 
           << " rt = " << SZpTe.rt << " Mpc" << endl;
    output << " a = " << SZpTe.a << "  " 
           << " b = " << SZpTe.b << "  "
           << " c = " << SZpTe.c << endl;
    output << " Tmin.T0 = " << SZpTe.Tmin_T0 << "  " 
           << " rcool = " << SZpTe.rcool << " kpc  "
           << " acool = " << SZpTe.acool << endl;
    
    output << "\n " << setfill('=') << setw(90) << "=" << endl; 

    return;
}

void SZ_cluster_profiles::show_cluster_parameters()
{
    output_cluster_parameters(cout);
    return;
}

void SZ_cluster_profiles::export_cluster_parameters(string fname)
{
    ofstream ofile(fname.c_str());
    
    output_cluster_parameters(ofile);
    
    ofile.close();    
}

double SZ_cluster_profiles::Get_rc_cm()
{ 
    return SZpNe.rc*const_Mpc*1.0e-3; 
}

//==================================================================================================
// computations (x, y, z all in units of the core-radius)
//==================================================================================================
double SZ_cluster_profiles::Ne(double x, double y, double z)  // Ne in cm^-3 
{
    double r=sqrt(x*x+y*y+z*z); // in units of rc
    double Ne0sq=SZpNe.n0*SZpNe.n0*pow(r, -SZpNe.alpha)
                  /pow(1.0+r*r, 3*SZpNe.beta-0.5*SZpNe.alpha)
                  /pow(1.0+pow(r*SZpNe.rc/SZpNe.rs, 3), SZpNe.eps/3.0);
    double Ne2sq=SZpNe.n02*SZpNe.n02/pow(1.0+pow(r*SZpNe.rc/SZpNe.rc2, 2), 3*SZpNe.beta2);
    
    return sqrt(Ne0sq*1.0e-6+Ne2sq*1.0e-2);
    
}

double SZ_cluster_profiles::Te(double x, double y, double z)  // Te in keV
{
    double r=sqrt(x*x+y*y+z*z); // in units of rc    
    double t=pow(r*SZpNe.rc/(SZpTe.rt*1000.0), -SZpTe.a)
              /pow(1.0+pow(r*SZpNe.rc/(SZpTe.rt*1000.0), SZpTe.b), SZpTe.c/SZpTe.b);
    double xx=pow(r*SZpNe.rc/SZpTe.rcool, SZpTe.acool), tcool=(xx+SZpTe.Tmin_T0)/(1.0+xx);
    
    return SZpTe.T0*tcool*t;
}

//==================================================================================================
//==================================================================================================
