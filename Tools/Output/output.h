//==================================================================================================
//
// Handling the output from SZpack, both to file and to screen
//
//==================================================================================================

#ifndef output_H
#define output_H

#include <fstream>
#include <string>
#include <vector>

#include "Parameters.h"
#include "SZpack.h"

struct method;

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

template <typename somestream, class fpGetter>
void IterateForNull(somestream &ofile, const vector<double> &vec, fpGetter &&getter, double x=0.0){
    for (const auto& val : vec){
        const Parameters &fp = getter(val);
        ofile << compute_null_of_SZ_signal(fp) - x << " ";
    }
    ofile << endl;
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
