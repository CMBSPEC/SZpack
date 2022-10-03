//==================================================================================================
//
// Handling the output from SZpack, both to file and to screen
//
//==================================================================================================

#include "output_multiple_scattering.h"
using namespace std;

extern Parameters parameters;

//TODO: All of these functions cout the xcmb as they pass through, instead of anything to determine if it's done the correct thing
// Also could do with putting this output behind a show mess variable and such.

typedef double (*t_0)(double);
typedef double (*t_av)(double, int);

    //TODO: check output format -- overwriting SetUpOutput from outputs? Look into this more carefully
template <typename somestream>
void SetUpOutput(string description, string fileAddition, Parameters &fp, somestream &ofile, int precision, bool printHeader){
    print_message(description);
    
    string filename = fp.outputPath+fileAddition+fp.fileEnding;

    ofile.open(filename.c_str());
    ofile.precision(precision);

    if (printHeader) {
        output_header(ofile, description, fp);
    }
}


void output_tau_l(int lmax, int np, vector<double> &xc, t_0 tau0, t_av tau_av, ofstream &ofile){
    vector<double> taulc(lmax+2);
    vector<double> taulv(lmax+2);
    
    double r = 0.0;
    
    for(int k = 0; k <= lmax; k += 2)
    {
        taulc[k] = tau_av(0.0, k)/(0.5*pow(tau0(0), 2));
        r += taulc[k];
        
        cout << k << " " << taulc[k] << " " << r << endl;
    }
    
    for(int m = 0; m < np; m++)
    {
        cout << xc[m] << endl;
        ofile << scientific << xc[m] << " ";
        
        double tauscale = 0.5*pow(tau0(0), 2);
        
        r = 0;
        for(int k=0; k<=lmax; k += 2)
        {
            taulv[k] = tau_av(xc[m], k);
            
            ofile << taulv[k]/tauscale/taulc[k] << " ";
            
            r += taulv[k]/(0.5*pow(tau0(xc[m]), 2));
        }
        
        ofile << r << " ";
        
        double tauscale2 = tau0(0)*tau0(xc[m]);
        double kappa_E_iso = 0.5*pow(tau0(xc[m]), 2);
        double kappa_T = taulv[0]+taulv[2]*0.1-kappa_E_iso;
        double kappa_S = -0.6*taulv[2];
        double kappa_E = taulv[0]+taulv[2]*0.1;
        
        ofile << tau0(xc[m])/tau0(0) << " ";
        ofile << kappa_E_iso/tauscale << " " << kappa_T/tauscale << " "
              << kappa_S/tauscale << " " << kappa_E/tauscale << " ";
        
        ofile << kappa_E_iso/tauscale2 << " " << kappa_T/tauscale2 << " "
              << kappa_S/tauscale2 << " " << kappa_E/tauscale2 << " ";

        ofile << endl;
    }
}

//==================================================================================================
//
// plot <tau_l> using const density Sphere
//
//==================================================================================================
void output_tau_l_sphere(int lmax, int np)
{
    ofstream ofile;
    SetUpOutput("Calculating the monopoles of tau l using the constant density sphere.",
                "tau_l_sphere", parameters, ofile, 10);
    
    vector<double> xc(np);
    init_xarr(0.001, 1.0-1.0e-4, &xc[0], np, 0, 0);

    output_tau_l(lmax, np, xc, constDensitySphere::tau0, constDensitySphere::tau_av, ofile);
    ofile.close();
}

//==================================================================================================
//
// plot <tau_l> using the Isothermal beta model
//
//==================================================================================================
void output_tau_l_beta(int lmax, int np, double beta)
{
    isothermalBeta::setBeta(beta);

    ofstream ofile;
    SetUpOutput("Calculating the monopoles of tau l using the isothermal beta model.",
                "tau_l_beta", parameters, ofile, 10);


    vector<double> xc(np);
    init_xarr(1.0e-1, 1.0e+2, &xc[0], np, 1, 0);
    
    output_tau_l(lmax, np, xc, isothermalBeta::tau0, isothermalBeta::tau_av, ofile);
    ofile.close();
}

//==================================================================================================
//
// plot of Y_l0 functions for different temperatures
//
//==================================================================================================
void output_DI2_E_temperature(Parameters &fp)
{
    ofstream ofile;
    string temperature = to_string(fp.Te) + "keV";
    SetUpOutput("Calculating multiply, D12_E_: "+ temperature, "DI2_E."+temperature, fp, ofile, 8);

    double The2 = pow(fp.calc.The, 2);
    vector<double> a(3);
    a[0] =  The2;
    a[1] = -0.4*The2;
    a[2] =  0.1*The2;
    a[3] = -3.0/70.0*The2;
    //TODO: why do we have all these different weighting for each multipole?
    //Like this seems like it could be useful to actually properly understand.

    for(int m = 0; m < fp.gridpoints; m++)
    {
        cout << fp.xcmb[m] << endl;
        double x3 = pow(fp.xcmb[m], 3);
        
        ofile << fp.xcmb[m] << " ";
        //
        for (int l = 0; l <= 3; l++){
            ofile << x3*compute_SZ_distortion_asym_multiple(fp.xcmb[m], l, fp.calc.The, 10)/a[l] << " ";
        }
        //
        for (int l = 0; l <= 3; l++){
            ofile << x3*compute_SZ_distortion_Patterson_multiple_Kernel_2D(fp.xcmb[m], l, fp.calc.The)/a[l] << " ";
        }

        ofile << endl;
    }
    ofile.close();
}

//==================================================================================================
//
// plot of Y_l0 functions for different temperatures
//
//==================================================================================================
void output_SZ_signal_l(Parameters &fp)
{
    ofstream ofile;
    string temperature = to_string(fp.Te) + "keV";
    SetUpOutput("Calculating multiply, signal: "+ temperature, "DI_signal."+temperature, fp, ofile, 8);
    
    for(int m = 0; m < parameters.gridpoints; m++)
    {
        cout << parameters.xcmb[m] << endl;

        ofile << parameters.xcmb[m] << " ";
        
        double Dn0 = compute_signal_combo(m, parameters);
        double x3  = pow(parameters.xcmb[m], 3);
        double fac = 500.0*parameters.Dtau*parameters.Dtau/2.0;
        
        ofile << fp.rare.Dn_DI_conversion()*x3*Dn0 << " ";
        
        double S0 = compute_SZ_distortion_Patterson_multiple_Kernel_2D(parameters.xcmb[m], 0, fp.calc.The);
        double S2 = compute_SZ_distortion_Patterson_multiple_Kernel_2D(parameters.xcmb[m], 2, fp.calc.The);
        
        ofile << fp.rare.Dn_DI_conversion()*x3*fac*S0 << " ";

        double kappa_T = 0.69 + 0.1*0.14 - 1.0; //TODO: where do these numbers come from?
        //Presumably these form different weightings, but why these ones? what's the point?
        ofile << fp.rare.Dn_DI_conversion()*x3*fac*( kappa_T*Dn0/parameters.Dtau + 0.69*S0+0.1*0.14*S2 ) << " ";
        
        kappa_T = 0.87 + 0.1*0.15 -1.0; //TODO: where do these numbers come from?
        ofile << fp.rare.Dn_DI_conversion()*x3*fac*( kappa_T*Dn0/parameters.Dtau + 0.87*S0+0.1*0.15*S2 ) << " ";

        ofile << endl;
    } 
    ofile.close();
}

//==================================================================================================
//
// plot emission/absorption terms
//
//==================================================================================================
void output_emission_absorption(Parameters &fp)
{
    ofstream ofile;
    SetUpOutput("Calculating emission absorption", "spatial_data", fp, ofile, 10);
    
    double tau0v = isothermalBeta::tau0(0);
    
    for(int m = 0; m < fp.gridpoints; m++)
    {
        double tauv = isothermalBeta::tau0(fp.xcmb[m]);
        double dum = isothermalBeta::z0_Int(fp.xcmb[m], 0)+0.5*isothermalBeta::z0_Int(fp.xcmb[m], 2);
        
        cout << fp.xcmb[m] << endl;
        
        ofile << scientific << fp.xcmb[m] << " " << tauv/tau0v << " " << pow(tauv/tau0v, 2) << " ";
        ofile << (dum-0.5*tauv*tauv)/(0.5*tau0v*tau0v) << " " << dum/(0.5*tau0v*tau0v) << endl;
    }
    
    ofile.close();
}

//==================================================================================================
void output_emission_absorption_kinetic(Parameters &fp)
{
    ofstream ofile;
    SetUpOutput("Calculating kinetic emission absorption", "spatial_data.kinetic", fp, ofile, 10);
    
    double tau0v = isothermalBeta::tau0(0);
    
    for(int m = 0; m < fp.gridpoints; m++)
    {
        double tauv = isothermalBeta::tau0(fp.xcmb[m]);
        
        cout << fp.xcmb[m] << endl;
        
        ofile << scientific << fp.xcmb[m] << " " << tauv/tau0v << " " << pow(tauv/tau0v, 2) << " ";
        
        for(int l = 0; l <= 3; l++)
        {
            double dum  = isothermalBeta::lambda_function(fp.xcmb[m], l);
            double dum1 = isothermalBeta::kappa_function(fp.xcmb[m], l);
            
            cout << scientific << l << " " << dum/(0.5*tau0v*tau0v) << " " << dum1/(0.5*tau0v*tau0v) << endl;
            ofile << dum/(0.5*tau0v*tau0v) << " " << dum1/(0.5*tau0v*tau0v) << " ";
        }
        
        cout << defaultfloat << endl;
        ofile << defaultfloat << endl;
    }
    
    ofile.close();
}

//==================================================================================================
//
// CMB anisotropies signal
//
//==================================================================================================
void output_CMB_isotropy_signals(Parameters &fp)
{
    ofstream ofile;
    SetUpOutput("Calculating the CMB isotropy signals", "CMB_scattering.Te_"+to_string(fp.Te)+"keV",fp, ofile, 8);
    
    for(int m = 0; m < fp.gridpoints; m++)
    {
        cout << fp.xcmb[m] << endl;
        
        ofile << fp.xcmb[m] << " ";
        
       for(int l = 0; l <= 5; l++){
            ofile << pow(fp.xcmb[m], 3)*compute_SZ_distortion_Patterson_multiple(fp.xcmb[m], l, fp.Te/const_me, "CMB")
                         /(fp.Te/const_me) << " ";
        }
        ofile << endl;
    }
    
    ofile.close();
   
    return;
}

//==================================================================================================
//
// plot first few correction functions
//
//==================================================================================================
void output_CMB_isotropy_Y_functions(Parameters &fp)
{
    ofstream ofile;
    SetUpOutput("Calculating correction functions.", "CMB_scattering.Y", fp, ofile, 8);

    vector<double> Y(3);
    
    for(int m = 0; m < fp.gridpoints; m++)
    {
        cout << fp.xcmb[m] << endl;
        ofile << fp.xcmb[m] << " ";

        double x3 = pow(fp.xcmb[m],3);

        compute_Y_asymptotic(fp.xcmb[m], Y);
        for(int l = 0; l < Y.size(); l++) { ofile << x3*Y[l] << " "; }
        
        compute_Yl0_k(fp.xcmb[m], 0, Y);
        for(int l = 0; l < Y.size(); l++) { ofile << x3*Y[l] << " "; }

        compute_Yl0_k(fp.xcmb[m], 1, Y);
        for(int l = 0; l < Y.size(); l++) { ofile << x3*Y[l] << " "; }

        compute_Yl0_k(fp.xcmb[m], 2, Y);
        for(int l = 0; l < Y.size(); l++) { ofile << x3*Y[l] << " "; }

        compute_Yl0_k(fp.xcmb[m], 3, Y);
        for(int l = 0; l < Y.size(); l++) { ofile << x3*Y[l] << " "; }

        ofile << endl;
    }
    
    ofile.close();
}

void output_lowest_order_signal_model(double b, double ISA, double scale, double The, double x3,
                                      double Y0, double Yl0, t_0 tau0, t_av tau_av, ofstream &ofile)
{
    double tauscale = 0.5*tau0(b)*tau0(b);
    double tau_0 = tau_av(b, 0)/tauscale;
    double tau_2 = tau_av(b, 2)/tauscale;
    double kappa_T = tau_0+0.1*tau_2-1.0;
    double kappa_S = -0.6*tau_2;
    double kappa_E = tau_0+0.1*tau_2;
    
    cout << tau_0 << " " << tau_2 << endl;
    
    ofile << scale*ISA*x3*kappa_T*Y0 << " ";
    ofile << scale*ISA*x3*kappa_S*The*Y0 << " ";
    ofile << scale*ISA*x3*kappa_E*The*Yl0 << " ";
    
    ofile << scale*ISA*x3*( (kappa_T+kappa_S*The)*Y0 + kappa_E*The*Yl0 ) << " ";
}
//==================================================================================================
//
// plot lowest order signal
//
//==================================================================================================
void output_lowest_order_signal(Parameters &fp, double b)
{
    ofstream ofile;
    SetUpOutput("Calculating lowest order signal.", "lowest_order", fp, ofile, 8);

    vector<double> Y(3);

    double ISA = fp.Dtau*0.5;
    double scale = 500;
    
    for(int m = 0; m < fp.gridpoints; m++)
    {
        cout << fp.xcmb[m] << endl;
        
        ofile << fp.xcmb[m] << " ";
        
        compute_Y_asymptotic(fp.xcmb[m], Y);
        double x3 = pow(fp.xcmb[m], 3);
        double Y0 = Y[0];
        
        ofile << x3*Y0 << " ";
        
        compute_Yl0_k(fp.xcmb[m], 0, Y);
        ofile << scale*ISA*x3*fp.calc.The*Y[0] << " ";

        // constant density sphere
        output_lowest_order_signal_model(b, ISA, scale, fp.calc.The, x3, Y0, Y[0], constDensitySphere::tau0,
                                         constDensitySphere::tau_av, ofile);

        // isothermal beta model
        output_lowest_order_signal_model(b, ISA, scale, fp.calc.The, x3, Y0, Y[0], isothermalBeta::tau0,
                                         isothermalBeta::tau_av, ofile);

        ofile << endl;
    }

    ofile.close();
}
