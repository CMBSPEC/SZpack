//==================================================================================================
//
// Handling the output from SZpack, both to file and to screen
//
//==================================================================================================

#ifndef output_H
#define output_H

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>

#include "ModeFunctions.h"
#include "physical_consts.h"
#include "routines.h"

using namespace std;

extern Parameters parameters;

//==================================================================================================
//
// Output for files
//
//==================================================================================================
template <typename somestream>
void output_header(somestream &ofile, string methodDescription, Parameters fparams = parameters);

//--------------------------------------------------------------------------------------------------
void output_to_file(ofstream &f, double x, double Dn, double T0_CMB, double unit_conversion);

template <typename somestream>
void SetUpOutput(string description, method Method, Parameters fp, somestream &ofile);

template <typename somestream>
void SetUpOutput(string description, string fileAddition, Parameters fp, somestream &ofile,
                 bool printHeader = false);

//==================================================================================================
//
// compute distortion of a method.
//
//==================================================================================================

void output_SZ_distortion(method Method, Parameters functionParameters = parameters);

//==================================================================================================
//
// compute precision of method. This returns a fiducially weighted accuracy.
//
//==================================================================================================

template <typename somestream>
void compute_precision_of_basis(double Te, method Method, Parameters fp, somestream &ofile, double DI_s);

//This mainly exists to make the template for compute compiles properly.
void output_precision_of_basis(double Te, method Method, Parameters fp);

//==================================================================================================
//
// compute null of SZ signal
//
//==================================================================================================

//TODO: This shouldn't have to be declared here, but misbehaves otherwise.
template <typename somestream>
void IterateForNull(somestream &ofile, vector<double> &vec, double &iterated, Parameters &fp, double x=0.0){
    for (int i = 0; i < vec.size(); i++){
        iterated = vec[i];
        ofile << compute_null_of_SZ_signal(fp) - x << " ";
    }
    ofile << endl;
    iterated = 0.0;
}

//This mainly exists to make the template for Iterate compiles properly.
void output_SZ_null(string fileAddition, Parameters fp, vector<double> &it_vec, double &iterated);

//==================================================================================================
//
// output derivatives. Calculated using the combo method.
//
//==================================================================================================

void output_derivatives(double Te, int kmax, Parameters fp = parameters);

//==================================================================================================
//==================================================================================================
#endif
