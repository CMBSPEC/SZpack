//==================================================================================================
//
// Author: Elizabeth Lee
// first implementation: November 2018
// last modification: November 2018
//
//==================================================================================================

#include "Parameters.h"

#include <cmath>

#include "global_functions.h"
#include "parser.h"
#include "physical_consts.h"
#include "routines.h"

//=============================================================================================
// Parameters Class
//=============================================================================================
MeansParameters::MeansParameters(){
    omegas.resize(3);
    sigmas.resize(3);
    assignOmegas(0.2, 0.1, 0.05); // omega^(1), omega^(2), omega^(3)
    assignSigmas(0.0, 0.0, 0.0);

    Omega = 0.2;
    Sigma = 0.0;
    kappa = 0.0;
}

void MeansParameters::assignOmegas(double omega0, double omega1, double omega2){
    omegas[0] = omega0;
    omegas[1] = omega1;
    omegas[2] = omega2;
}
void MeansParameters::assignSigmas(double sigma0, double sigma1, double sigma2){
    sigmas[0] = sigma0;
    sigmas[1] = sigma1;
    sigmas[2] = sigma2;
}

void MeansParameters::copyMeansParameters(MeansParameters original){
    assignOmegas(original.omegas[0],original.omegas[1],original.omegas[2]);
    assignSigmas(original.sigmas[0],original.sigmas[1],original.sigmas[2]);

    Omega = original.Omega;
    Sigma = original.Sigma;
    kappa = original.kappa;
}

//=================================================================================================
KernelParameters::KernelParameters(){
    smin = -0.6; 
    smax = 0.6; 
    l = 2;  
}

void KernelParameters::copyKernelParameters(KernelParameters original){
    smin = original.smin;
    smax = original.smax;
    srange = original.srange;
    l = original.l;
}

//=================================================================================================
RareParameters::RareParameters(){
    setTCMB(2.726);
    RunMode = "";
}

void RareParameters::setTCMB(double T0){
    T0_CMB = T0;
    Dn_DI_conv = 13.33914078*pow(T0_CMB, 3);
    x_nu_conv = T0_CMB/const_h_kb/1.0e9;
}

void RareParameters::copyRareParameters(RareParameters original){
    setTCMB(original.T0_CMB);
    RunMode = original.RunMode;
}

double RareParameters::TCMB(){ return T0_CMB; }
double RareParameters::Dn_DI_conversion() { return Dn_DI_conv; };
double RareParameters::x_nu_conversion() {return x_nu_conv; }

//=================================================================================================
Parameters::Parameters(){
    xmin = 0.1;         //(>=0.1)
    xmax = 30.0;        //(<=50.0)
    gridpoints = 50;    //number of points for the log-frequency grid (>=2)
    
    Dtau = 0.01;        //optical depth > 0
    Te = 0.02*const_me; //electron temperature in keV

    betac = 0.01;       //peculiar velocity of cluster with respect to CMB frame (betac = v/c)
    muc = 1.0;          //direction cosine of the cluster velocity with respect to the line-of-sight in the CMB rest frame
    betao = 0.001241;   //peculiar velocity of observer with respect to CMB frame (betao = v/c)
    muo = 0.0;          //direction cosine for the line-of-sight with respect to the observers velocity. The angle is measured in the observer frame

    T_order = 5;        //Order of expansion with respect to Te. Used only for CNSN and aymptotic methods.
    beta_order = 2;     //Order of expansion with respect to beta. Used only for CNSN and aymptotic methods.
    relative_accuracy = 1.0e-4; //relative accuracy for numerical integral (only for "5D" and "3D" runmode)
    
    kmax = 5;           //kmax (only for 'CNSNopt' mode; depends on accuracy_level; see Table 1 of CSNN2012)
    accuracy_level = 1; //accuracy_level (0, 1, 2, 3; only for 'CNSNopt' runmode)

    Te_max = 30.0;      //Relevant only for the moment method running using CNSNopt mode.

    outputPath = "./outputs/";   //path for output
    fileEnding = ".dat";        //addition to name of files at the very end

    mode = "3D";        //the default mode used. The mode can still be chosen from the command line.

    kernel = KernelParameters();
    means = MeansParameters();
    rare = RareParameters();

    T = TemperatureIterators();
    D = DcomputeValues();

    PopulateRunningArrays();
    setCalcValues();
}

void Parameters::copyParameters(Parameters copyParameters){
    xmin = copyParameters.xmin;
    xmax = copyParameters.xmax;
    gridpoints = copyParameters.gridpoints;
    
    Dtau = copyParameters.Dtau;
    Te = copyParameters.Te;

    betac = copyParameters.betac;
    muc = copyParameters.muc;
    betao = copyParameters.betao;
    muo = copyParameters.muo;

    T_order = copyParameters.T_order;
    beta_order = copyParameters.beta_order;
    relative_accuracy = copyParameters.relative_accuracy;
    
    kmax = copyParameters.kmax;
    accuracy_level = copyParameters.accuracy_level;

    outputPath = copyParameters.outputPath;
    fileEnding = copyParameters.fileEnding;

    mode = copyParameters.mode;

    calc.copyCalculatedParameters(copyParameters.calc);
    kernel.copyKernelParameters(copyParameters.kernel);
    means.copyMeansParameters(copyParameters.means);
    rare.copyRareParameters(copyParameters.rare);
    T.copyTemperatureIterators(copyParameters.T);
    D.copyDcomputeValues(copyParameters.D);

    PopulateRunningArrays();
}

void Parameters::PopulateRunningArrays(){
    //set the arrays for the functions to run over.
    xcmb.resize(gridpoints);
    kernel.srange.resize(gridpoints);
    init_xarr(xmin, xmax, &xcmb[0], gridpoints, 1, 0); // change 1 --> 0 to have linear grid in x
    init_xarr(kernel.smin, kernel.smax, &kernel.srange[0], gridpoints, 0, 0);
}

void Parameters::SetParametersFromFile(string filename){
    if (filename == "") {
        print_error("No filename passed. Parameters unchanged.");
        return;
    }
    int ival;
    double dval;
    string sval;
    vector<double> vec3;
    vec3.resize(3);
    
    file_content pfc;
    parser_read_file(filename, pfc, 0);
    
    //==============================================================================================
    // Set Values
    //==============================================================================================
    if (parser_read(pfc, "outputPath", sval)){
        outputPath = sval;
    }
    if (parser_read(pfc, "fileEnding", sval)){
        fileEnding = sval;
    }

    if (parser_read(pfc, "Mode", sval)){
        mode = sval;
    }

    if (parser_read(pfc, "xmin", dval)){
        xmin = dval;
    }
    if (parser_read(pfc, "xmax", dval)){
        xmax = dval;
    }
    if (parser_read(pfc, "gridpoints", ival)){
        gridpoints = ival;
    }

    if (parser_read(pfc, "Dtau", dval)){
        Dtau = dval;
    }
    if (parser_read(pfc, "Te", dval)){
        updateT(dval);
    }

    if (parser_read(pfc, "betac", dval)){
        betac = dval;
    }
    if (parser_read(pfc, "muc", dval)){
        muc = dval;
    }
    if (parser_read(pfc, "betao", dval)){
        betao = dval;
    }
    if (parser_read(pfc, "muo", dval)){
        muo = dval;
    }

    if (parser_read(pfc, "T_order", ival)){
        T_order = ival;
    }
    if (parser_read(pfc, "beta_order", ival)){
        beta_order = ival;
    }
    if (parser_read(pfc, "relative_accuracy", dval)){
        relative_accuracy = dval;
    }

    if (parser_read(pfc, "kmax", dval)){
        kmax = dval;
    }
    if (parser_read(pfc, "accuracy_level", dval)){
        accuracy_level = dval;
    }

    if (parser_read(pfc, "Te_max", dval)){
        Te_max = dval;
    }

    if (parser_read(pfc, "smin", dval)){
        kernel.smin = dval;
    }
    if (parser_read(pfc, "smax", dval)){
        kernel.smax = dval;
    }
    if (parser_read(pfc, "l", ival)){
        kernel.l = ival;
    }

    if (parser_read_3vector(pfc, "omegas", vec3)){
        means.assignOmegas(vec3[0],vec3[1],vec3[2]);
    }
    if (parser_read_3vector(pfc, "sigmas", vec3)){
        means.assignSigmas(vec3[0],vec3[1],vec3[2]);
    }

    if (parser_read(pfc, "Omega", dval)){
        means.Omega = dval;
    }
    if (parser_read(pfc, "Sigma", dval)){
        means.Sigma = dval;
    }

    if (parser_read(pfc, "TCMB", dval)){
        rare.setTCMB(dval);
    }
    if (parser_read(pfc, "RunMode", sval)){
        //Takes values such as "all", "kin", "monopole", "dipole", "quadrupole", "monopole_corr" or "full".
        //If left blank or unrecognised, defaults to containing "all" (or "full", for 5D)
        rare.RunMode = sval;
    }

    //==============================================================================================
    // Clean up
    //==============================================================================================
    parser_free(pfc);

    PopulateRunningArrays();
    setCalcValues();
}

int Parameters::k_inRange(int k){
    if (k < 0 || k >= gridpoints){
        print_error("k out of range. k set to 0.");
        return 0;
    }
    return k;
}

bool Parameters::CheckValues(){
    bool valuesGood = true;
    if(xmin < 0.1){
        print_error("xmin must be >= 0.1.");
        valuesGood = false;
    }
    if (xmax > 50.0){
        print_error("xmax must be <= 50.");
        return false;
    }
    if (gridpoints < 2){
        print_error("Must be at least 2 points for log-frequency grid (Require gridpoints >= 2).");
        valuesGood = false;
    }
    if (Dtau < 0.0){
        print_error("Optical depth, Dtau, must be positive.");
        valuesGood = false;
    }
    if (Te <= 0.0){
        print_error("Electron Temperature, Te, must be greater than 0.");
        valuesGood = false;
    }
    if (betac < 0 || betac > 0.1){
        print_error("betac is out of bounds, must be between 0 and 0.1.");
        valuesGood = false;
    }
    if (betao < 0 || betao > 0.1){
        print_error("betao is out of bounds, must be between 0 and 0.1.");
        valuesGood = false;
    }
    if (muc < -1.0 || muc > 1.0){
        print_error("muc is a cosine, must take values between -1 and 1.");
        valuesGood = false;
    }
    if (muo < -1.0 || muo > 1.0){
        print_error("muo is a cosine, must take values between -1 and 1.");
        valuesGood = false;
    }
    if (Te_max < 0.0){
        print_error("Temperature maximum must be positive.");
        valuesGood = false;
    }
    //Check and fix the other values.
    CheckValuesNumericalIntegration();
    CheckValues_orders();
    CheckValuesCNSNopt();
    CheckValues_Kernel();
    CheckValues_means();
    CheckValues_rare();
    D.checkValues();

    return valuesGood;
}

void Parameters::CheckValuesNumericalIntegration(){
    if( 1e-8 >= relative_accuracy || relative_accuracy >= 0.001){
        print_error("The relative accuracy set is out of bounds. It must be >= 1e-8 and <= 0.001. Set = 1e-4.");
        relative_accuracy = 1.0e-4;
    }
}

void Parameters::CheckValuesCNSNopt(){
    if(0 >= accuracy_level || accuracy_level >= 4){
        print_error("The accuracy level set is out of bounds. It must take a value of 0, 1, 2 or 3. Set = 3.");
        accuracy_level = 3;
    }

    if(1 >= kmax || kmax >= 7){
        print_error("The kmax level set is out of bounds. It must take a value of 2, 3, 4 5 or 6. Set = 6.");
        accuracy_level = 3;
    }
}

void Parameters::CheckValues_orders(){
    if(T_order <0 || T_order > 21){
        print_error("The T_order set is out of bounds. It must be at least 0 and at most 21 for CNSN method. Set = 10.");
        T_order = 10;
    }

    if(beta_order < 0 || beta_order > 2){
        print_error("The beta_order set is out of bounds. It must be at least 0 and most 2. It is set = 2.");
        beta_order = 2;
    }
}

void Parameters::CheckValues_Kernel(){
    if(kernel.l<0){
        print_error("l cannot be less than 0. l set = 0.");
        kernel.l = 0;
    }
    //TODO: smin smax checks.
}

void Parameters::CheckValues_means(){
    if(means.Omega<0){
        print_error("Omega cannot be less than 0. l set = 0.");
        means.Omega = 0.0;
    }
    if(means.omegas[0]<0){
        print_error("omegas[0] cannot be less than 0. l set = 0.");
        means.omegas[0] = 0.0;
    }
    if(means.kappa<0){
        print_error("kappa cannot be less than 0. l set = 0.");
        means.kappa = 0.0;
    }
}

void Parameters::CheckValues_rare(){
    if(rare.TCMB()<0.0){
        print_error("TCMB should be positive. Set to 2.726 k.");
        rare.setTCMB(2.726);
    }
    if(rare.TCMB()>30.0){
        print_error("TCMB too large, SZ approximation will no longer work in this regime. Set to 2.726 k.");
        rare.setTCMB(2.726);
    }

    if (rare.RunMode != "full" && rare.RunMode != "all" && rare.RunMode != "kin" && rare.RunMode != "monopole" && rare.RunMode != "dipole"
        && rare.RunMode != "quadrupole" && rare.RunMode != "monopole_corr" && rare.RunMode != ""){
        print_error("The expansion runmode was not recognised. Set to full or all depending on the mode chosen.");
        rare.RunMode = "";
    }
}

void Parameters::Set_x(double x_min, double x_max, int gridpoints_i){
    if (x_min <= 0.0) {
        print_error("xmin must be > 0.0. xmin is set to 0.1.");
        x_min = 0.1;
    }
    xmin = x_min;
    xmax = x_max;
    gridpoints = gridpoints_i;
    PopulateRunningArrays();
}

void Parameters::Set_x_from_nu(double nu_min, double nu_max, int gridpoints_i){
    if (nu_min <= 0.0) {
        print_error("nu_min must be > 0.0. nu_min is set to 5.");
        nu_min = 5.0;
    }
    //Takes nu in GHZ
    xmin = nu_min/rare.x_nu_conversion();
    xmax = nu_max/rare.x_nu_conversion();
    gridpoints = gridpoints_i;
    PopulateRunningArrays();
}

vector<double> Parameters::get_nucmb(){
    vector<double> nucmb;
    nucmb.resize(gridpoints);
    for (int i=0; i < gridpoints; i++){
        nucmb[i] = xcmb[i]*rare.x_nu_conversion();
    }
    return nucmb;
}

void Parameters::Set_s(double s_min, double s_max, int gridpoints_i){
    kernel.smin = s_min;
    kernel.smax = s_max;
    gridpoints = gridpoints_i;
    PopulateRunningArrays();
}

void Parameters::updateT(double T){
    if (T<=0.0){
        print_error("Te must be greater than 0. For non relativistic method, please use specific function. Te = " 
                    + to_string(T) + "keV");
        return;
    }
    Te = T;
    calc.The = Te/const_me; // conversion keV --> The = k Te / me c^2
}

void Parameters::setCalcValues(){
    calc.The = Te/const_me; // conversion keV --> The = k Te / me c^2
    calc.gammao = 1.0/sqrt(1.0-betao*betao);
    calc.gammac = 1.0/sqrt(1.0-betac*betac);

    // transformation of muc to mucc
    calc.mucc = (muc - betac)/(1.0-betac*muc);
        
    // transformation of xo to xc
    calc.xfac = calc.gammac*(1.0-betac*muc)*calc.gammao*(1.0+betao*muo); // xfac*xo = xc;
    calc.xfacCMB = calc.gammao*(1.0+betao*muo);
    
    // change to Nozawa 2006 convention for scattering optical depth
    calc.CNSN_Dtau = Dtau*calc.gammac*(1.0-betac*calc.mucc);

    calc.betac_para = betac*muc;
    calc.betac2_perp = betac*betac*(1.0-muc*muc);
}

//=================================================================================================
CalculatedParameters::CalculatedParameters(){
    The = gammao = gammac = mucc = xfac = CNSN_Dtau = 0.0;
    betac_para = betac2_perp = 0.0;
}

void CalculatedParameters::copyCalculatedParameters(CalculatedParameters original){
    The = original.The;
    gammao = original.gammao;
    gammac = original.gammac;
    mucc = original.mucc;
    xfac = original.xfac;
    xfacCMB = original.xfacCMB;
    CNSN_Dtau = original.CNSN_Dtau;
    betac_para = original.betac_para;
    betac2_perp = original.betac2_perp;
}

//=================================================================================================
TemperatureIterators::TemperatureIterators(){
    Tmin = Tmax = 0;
    np = 0;
    Trange.resize(np);
}

TemperatureIterators::TemperatureIterators(double min, double max, int gridpoints, int mode){
    Tmin = min;
    Tmax = max;
    np = gridpoints;
    Trange.resize(np);
    init_xarr(Tmin, Tmax, &Trange[0], np, mode, 0); //linear! change mode to 1 for log
}

TemperatureIterators::TemperatureIterators(vector<double> range){
    Trange = range;
    np = range.size();
    Tmin = range[0];
    Tmax = range[np-1];
}

void TemperatureIterators::copyTemperatureIterators(TemperatureIterators original){
    Tmin = original.Tmin;
    Tmax = original.Tmax;
    np = original.np;
    Trange = original.Trange;
}

//=================================================================================================
DcomputeValues::DcomputeValues(){
    dThe = dbeta_para = dbeta2_perp = 0;
    dDn_dThe.resize(dThe+1);
}

void DcomputeValues::setValues(int DThe, int Dbeta_para, int Dbeta2_perp){
    dThe = DThe;
    dbeta_para = Dbeta_para;
    dbeta2_perp = Dbeta2_perp;
    checkValues();
    dDn_dThe.resize(dThe+1);
}

void DcomputeValues::copyDcomputeValues(DcomputeValues original){
    setValues(original.dThe, original.dbeta_para, original.dbeta2_perp);
    dDn_dThe = original.dDn_dThe;
}

void DcomputeValues::checkValues(){
    if(dThe < 0 || dThe > 4){
        print_error("dThe must be between 0 and 4. dThe set to 4.");
        dThe = 4;
    }

    if(dbeta_para < 0 || dbeta_para > 2){
        print_error("dbeta_para must be between 0 and 2. dbeta_para set to 0.");
        dbeta_para = 0;
    }

    if(dbeta2_perp < 0 || dbeta2_perp > 1){
        print_error("dbeta2_perp must be 0 or 1. dbeta2_perp set to 0.");
        dbeta2_perp = 0;
    }
}
