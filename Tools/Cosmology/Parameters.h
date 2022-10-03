//==================================================================================================
//
// Author: Elizabeth Lee
// first implementation: November 2018
// last modification: November 2018
//
//==================================================================================================

#ifndef Parameters_H
#define Parameters_H

#include <iostream>
#include <string>
#include <vector>
#include "parser.h"
#include "global_functions.h"
#include "routines.h"
#include "physical_consts.h"

using namespace std;

//=============================================================================================
// Parameters Class
//=============================================================================================

class MeansParameters{ //Contains the Parameters needed for running means
    public:
    double Omega, Sigma, kappa; // For temperature dispersion
    vector<double> omegas, sigmas; // For higher order temperature dispersion

    MeansParameters();

    void copyMeansParameters(MeansParameters copyMeansParameters);

    void assignOmegas(double omega0, double omega1, double omega2);
    void assignSigmas(double sigma0, double sigma1, double sigma2);
};

class KernelParameters{ //Contains the Parameters needed for running means
    public:
    double smin, smax;
    vector<double> srange;
    int l;

    KernelParameters();

    void copyKernelParameters(KernelParameters copyKernelParameters);
};

class TemperatureIterators{
    public:
    double Tmin, Tmax;
    int np;
    vector<double> Trange;

    TemperatureIterators();
    TemperatureIterators(double min, double max, int gridpoints, int mode = 0);
    TemperatureIterators(vector<double> Trange);

    void copyTemperatureIterators(TemperatureIterators original);
};

class DcomputeValues{
    public:
    int dThe, dbeta_para, dbeta2_perp;
    vector<double> dDn_dThe;

    DcomputeValues();
    void setValues(int DThe, int Dbeta_para, int Dbeta2_perp);
    void copyDcomputeValues(DcomputeValues original);
    void checkValues();
};

class RareParameters{
    private:
    double T0_CMB;
    double Dn_DI_conv;  
    double x_nu_conv;

    public:
    string RunMode; //Takes values such as "all", "kin", "monopole", "dipole", "quadrupole", "monopole_corr" or "full".
    //If left blank or unrecognised, defaults to containing "all" (or "full", for 5D)

    RareParameters();
    void copyRareParameters(RareParameters original);
    void setTCMB(double T0);
    double TCMB();
    double Dn_DI_conversion(); // conversion factor x^3 Dn --> DI in MJy / sr
    double x_nu_conversion();
};

class CalculatedParameters{
    public:
    double The, gammao, gammac, mucc, xfac, xfacCMB, CNSN_Dtau;
    double betac_para, betac2_perp;

    CalculatedParameters();
    void copyCalculatedParameters(CalculatedParameters original);
};

class Parameters{
    public:
    double xmin, xmax, Dtau, Te, betac, muc, betao, muo;
    vector<double> xcmb;
    int gridpoints;
    int T_order, beta_order; //only for Asymptotic & CNSN
    double relative_accuracy; //only for 5D & 3D
    int kmax, accuracy_level; //only for CNSNopt

    double Te_max; //only for the moment method running

    string outputPath, fileEnding; //for outputting data

    string mode;

    CalculatedParameters calc;
    KernelParameters kernel;
    MeansParameters means;
    RareParameters rare;

    TemperatureIterators T;
    DcomputeValues D;

    //default values set in constructor
    Parameters();
    void copyParameters(Parameters copyParameters);
    void PopulateRunningArrays();
    void SetParametersFromFile(string filename);

    void Set_x(double x_min, double x_max, int gridpoints_i);
    void Set_x_from_nu(double nu_min, double nu_max, int gridpoints_i);
    vector<double> get_nucmb();
    void Set_s(double s_min, double s_max, int gridpoints_i);
    void updateT(double T);
    void setCalcValues(); //needs to be called if muc, muo, betac or betao are changed!

    int k_inRange(int k);
    bool CheckValues();
    private:
    void CheckValuesNumericalIntegration();
    void CheckValuesCNSNopt();
    void CheckValues_orders();
    void CheckValues_Kernel();
    void CheckValues_means();
    void CheckValues_rare();
};

#endif
