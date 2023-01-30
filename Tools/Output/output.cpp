//==================================================================================================
//
// Author: Elizabeth Lee
// first implementation: December 2018
// last modification: March 2020
//
//==================================================================================================

#include "output.h"

//==================================================================================================
//
// Output for files
//
//==================================================================================================
template <typename somestream>
void output_header(somestream &ofile, string methodDescription, Parameters fparams)
{
    ofile.precision(6);

    ofile << setfill('#') << setw(90) << "#" << endl << "#" << endl;  
    ofile << "# SZ signal computed with " << SZpack_version << endl;
    ofile << "# Runmode: " << methodDescription << endl;
    ofile << "#\n# General parameters:" << endl;
    ofile << "# xmin= " << fparams.xmin << " xmax= " << fparams.xmax << " np= " << fparams.gridpoints << endl;
    ofile << "# Dtau= " << fparams.Dtau << endl;
    ofile << "# Te= " << fparams.Te << "keV (The ~ " << fparams.calc.The << ")" << endl;
    ofile << "# betac= " << fparams.betac << " muc= " << fparams.muc << endl;
    ofile << "# betao= " << fparams.betao << " muo= " << fparams.muo << endl;
    ofile << "#\n# Runmode specific parameters:\n# eps= " << fparams.relative_accuracy << endl;
    ofile << "# Te_order= " << fparams.T_order << " beta_order= " << fparams.beta_order << endl;
    ofile << "# Accuracy_level= " << fparams.accuracy_level << " Te_max= " << fparams.Te_max << endl;
    ofile << "#\n" << setfill('#') << setw(90) << "#" << endl; 
    ofile << "# Here a CMB temperature of T0 = " << fparams.rare.TCMB() << " K is used." << endl;
    ofile << "#\n" << setfill('#') << setw(90) << "#" << endl; 
    
    return;
}

void output_to_file(ofstream &f, double x, double Dn, double T0_CMB, double unit_conversion)
{
    double x3=pow(x, 3), nu=x*T0_CMB/const_h_kb/1.0e+9;
    f << x << " " << nu << " " << x3*Dn << " " << unit_conversion*x3*Dn << endl;
    
    if(show_mess >= 1) cout << " x= " << x << " nu [GHz]" << nu << " "
                       << x3*Dn << " " << unit_conversion*x3*Dn << endl;
    
    return;
}

template <typename somestream>
void SetUpOutput(string description, string fileAddition, Parameters fp, somestream &ofile, bool printHeader){
    print_message(description);
    
    string filename = fp.outputPath+fileAddition+fp.fileEnding;

    ofile.open(filename.c_str());
    ofile.precision(16);

    if (printHeader) {
        output_header(ofile, description, fp);
    }
}

template <typename somestream>
void SetUpOutput(string description, method Method, Parameters fp, somestream &ofile){
    print_message(description + Method.description);
    
    string filename = fp.outputPath+Method.fileAddition+fp.fileEnding;

    ofile.open(filename.c_str());
    ofile.precision(16);

    output_header(ofile, Method.description, fp);
}

//==================================================================================================
//
// compute distortion
//
//==================================================================================================

void output_SZ_distortion(method Method, Parameters functionParameters)
{
    ofstream ofile;
    SetUpOutput("Carrying out ", Method, functionParameters, ofile);
    ofile << "#\n# Output format: x = (h nu/k T0) | nu [GHz] | x^3 Dn(x) | DI(x) in MJy/sr " << endl;

    if(!CNSN2012_convention) functionParameters.Dtau = functionParameters.calc.CNSN_Dtau;

    double TCMB = functionParameters.rare.TCMB(), unit_conversion = functionParameters.rare.Dn_DI_conversion();

    for(int k = 0; k < functionParameters.gridpoints; k++)
    {
        double dum = Method.modeFunction(k, functionParameters);
        output_to_file(ofile, functionParameters.xcmb[k], dum, TCMB, unit_conversion);
    }
    
    ofile.close();
}

void output_SZ_kernel(method Method, Parameters functionParameters)
{ //TODO: This is probably curently wrong. To be examined more closely.
    ofstream ofile;
    SetUpOutput("Carrying out ", Method, functionParameters, ofile);
    ofile << "#\n# Output format: x = (h nu/k T0) | nu [GHz] | x^3 Dn(x) | DI(x) in MJy/sr " << endl;

    double TCMB = functionParameters.rare.TCMB(), unit_conversion = functionParameters.rare.Dn_DI_conversion();

    for(int k = 0; k < functionParameters.gridpoints; k++)
    {
        double dum = Method.modeFunction(k, functionParameters);
        output_to_file(ofile, functionParameters.kernel.srange[k], dum, TCMB, unit_conversion);
    }
    
    ofile.close();
}

//==================================================================================================
//
// compute precision of method. This returns a fiducially weighted accuracy.
//
//==================================================================================================

template <typename somestream>
void compute_precision_of_basis(double Te, method Method, Parameters fp, somestream &ofile, double DI_s)
{
    vector<double> ya_3D(fp.gridpoints), ya(fp.gridpoints);

    fp.updateT(Te);
    ofile << fp.Te << " ";        
    if(show_mess >=1) cout << " Te= "<< fp.Te << endl;
    
    //==================================================================
    // compute precise solution
    //==================================================================
    compute_signal_precise(ya_3D, fp, true);

    //==================================================================
    // compute approximation
    //==================================================================
    for(int m = 0; m <= 10; m++)
    {
        fp.T_order = m;
        for(int k = 0; k < fp.gridpoints; k++) {
            ya[k] = pow(fp.xcmb[k], 3)*Method.modeFunction(fp.xcmb[k], fp);
        }
        //==============================================================
        // find maximal absolute deviation and save to file
        //==============================================================
        double max = 0.0;
        
        for(int k = 0; k < fp.gridpoints; k++){
            if(fabs(ya[k]-ya_3D[k])>=max){
                max = fabs(ya[k]-ya_3D[k]);
            }
        }
        ofile << max/DI_s << " ";
        if(show_mess >=1) cout << max/DI_s << endl;
    }
    
    ofile << endl;
}

//This mainly exists to make the template for compute compiles properly.
void output_precision_of_basis(double Te, method Method, Parameters fp){
    fp.fileEnding = "_Accuracy"+fp.fileEnding;

    ofstream ofile;
    SetUpOutput("Computing the precision of ", Method, fp, ofile);
    ofile << "#\n# Output format: Te [K] | Max error/fiducial accuracy " << endl;

    fp.gridpoints = 300; //This overwrites the user set value to ensure better accuracy
    fp.PopulateRunningArrays();

    double DI_s = 4.881e-5; // fiducial accuracy according to CSNN 2012
    compute_precision_of_basis(Te, Method, fp, ofile, DI_s);
    
    ofile.close();
}

//==================================================================================================
//
// output derivatives. Calculated using the combo method.
//
//==================================================================================================

void output_derivatives(double Te, int kmax, Parameters fp)
{
    if (kmax>4) {
        print_error("Too many derivatives require. kmax must be <= 4. kmax set = 4.");
        kmax = 4;
    }
    ofstream ofile;

    string fileAddition = "SZ_moments.derivs."+to_string(int(Te))+"keV";
    SetUpOutput("Computing derivatives of S in CMB rest frame", fileAddition, fp, ofile);
    ofile << "#\n# Output format: x = (h nu/k T0) | th-SZ terms in MJy/sr | k-SZ terms in MJy/sr "; 
    ofile << "| k-SZ second order terms in MJy/sr | k-SZ perpendicular terms in MJy/sr ";
    ofile << "| x derivatives in MJy/sr " << endl; 
    
    Parameters temp = Parameters();
    temp.copyParameters(fp);

    temp.Te = Te;
    temp.betao = temp.muo = 0;
    temp.setCalcValues();
 
    double norm=temp.rare.Dn_DI_conversion();   // in physical units; Set == 1 for dimensionless
    
    for(int k=0; k<temp.gridpoints; k++)
    {
        double x3=pow(temp.xcmb[k], 3);
        ofile << temp.xcmb[k] << " ";
        
        // th-SZ terms
        temp.D.setValues(kmax, 0, 0);
        Dcompute_signal_combo_CMB(temp.xcmb[k], temp);
        for(int m=0; m<=kmax; m++) ofile << norm*temp.D.dDn_dThe[m]*x3 << " ";
        
        // k-SZ terms
        temp.D.setValues(kmax, 1, 0);
        Dcompute_signal_combo_CMB(temp.xcmb[k], temp);
        for(int m=0; m<=kmax; m++) ofile << norm*temp.D.dDn_dThe[m]*x3 << " ";
        
        // k-SZ second order terms
        temp.D.setValues(kmax, 2, 0);
        Dcompute_signal_combo_CMB(temp.xcmb[k], temp);
        for(int m=0; m<=kmax; m++) ofile << norm*temp.D.dDn_dThe[m]*x3 << " ";
        
        // k-SZ perpendicular terms (beta_perp^2)
        temp.D.setValues(kmax, 0, 1);
        Dcompute_signal_combo_CMB(temp.xcmb[k], temp);
        for(int m=0; m<=kmax; m++) ofile << norm*temp.D.dDn_dThe[m]*x3 << " ";

        //------------------------------------------------------------
        // x-derivatives
        //------------------------------------------------------------
        double x0 = temp.xcmb[k];
        double dI = Dcompute_signal_combo_for_x(x0, temp, 1);
        double d2I = Dcompute_signal_combo_for_x(x0, temp, 2);
        
        ofile << norm*dI << " " << norm*d2I << endl;
        
        cout  << " x= " << x0 << " " << norm*temp.D.dDn_dThe[0]*x3 << endl;
    }
    
    ofile.close();
}