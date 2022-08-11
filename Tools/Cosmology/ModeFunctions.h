//==================================================================================================
//
// Handling the various modes that SZ Pack can be run through
//
//==================================================================================================

#ifndef modeFunction_H
#define modeFunction_H

#include <string>
#include "Parameters.h"
#include "SZpack.h"

using namespace std;

typedef double (*modefunction)(int, Parameters &);

static double plainFunction(int a, Parameters &plainParameters){
    return (double) a;
}

struct method{
    string description;
    string fileAddition;
    modefunction modeFunction;

    method(modefunction ModeFunction, string FileAddition, string Description){
        description = Description;
        fileAddition = FileAddition;
        modeFunction = ModeFunction;
    }

    method(){
        description = "";
        fileAddition = "";
        modeFunction = &plainFunction;
    }
};

class distortionModes{
    public:
    method Int5D;
    method Int3D;
    //method Kernel;
    method Asymptotic;
    method CNSN;
    method CNSNopt;
    method Combo;
    method Means;
    method Means_Yweighted;
    method RelativisticCorrections;
    method TemperatureDispersion;

    distortionModes(){
        Int5D = method(&compute_signal_5D,"SZ_Integral.5D","5D Integral");
        Int3D = method(&compute_signal_3D,"SZ_Integral.3D","3D Integral");
        //Kernel = method(&compute_signal_Kernel,"SZ_Integral.Kernel","Kernel Integral");
        Asymptotic = method(&compute_signal_asymptotic,"SZ_asymptotic","Asymptotic Expansion of Collision Integral (Te<13 keV at most)");
        CNSN = method(&compute_signal_CNSN,"SZ_CNSN_basis","Integral using improved basis of CNSN 2012 (2keV < Te < 75 keV)");
        CNSNopt = method(&compute_signal_CNSN_opt,"SZ_CNSN_opt_basis","Integral using improved basis of CNSN 2012 with optimization of temperature terms");
        Combo = method(&compute_signal_combo,"SZ_combo","Integral using combo of asymptotic expansion and CNSN basis (Te < 75 keV at most)");
        Means = method(&compute_signal_means_tw,"SZ_means","Integral using expansion around mean values - tau weighted");
        Means_Yweighted = method(&compute_signal_means_yw,"SZ_means_yw","Integral using expansion around mean values - y weighted");
        RelativisticCorrections = method(&compute_signal_RelCorrs,"SZ_rel_corrs","Integral using expansion around mean values and computing relativistic correction");
        TemperatureDispersion = method(&compute_signal_TDispersion,"SZ_temp_dis","Expansion around mean values and computing temperature dispersion");
    }
};
//==================================================================================================
//==================================================================================================
#endif
