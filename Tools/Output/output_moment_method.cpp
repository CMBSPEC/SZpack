//==================================================================================================
//
// Handling the output from SZpack, both to file and to screen
//
//==================================================================================================

#include "output_moment_method.h"

//==================================================================================================
template<typename somestream>
void output_moment_vector(somestream &output, const vector<double> &mv, Parameters fp, SZ_moment_method SZM)
{
    output.precision(6);
    output << scientific;

    output << "\n " << setfill('=') << setw(95) << "=" << endl;  

    //==============================================================================================
    output << "\n asymptotic expansion part" << endl;
   
    output << " " << setfill('-') << setw(42) << "-" << " Region 0 " << setw(43) << "-" << endl;
    
    int dn = fp.kmax + 1;
    for(int m = 0; m < dn; m++)
    {
        string mstring = to_string(m);
        output << " y^("+mstring+ ") = " << mv[m];

        if(fp.beta_order > 0)
        {
            output << "\t b_1^("+mstring+ ") = " << mv[dn+m];
            
            if(fp.beta_order > 1)
            {
                output << "\t b_0^("+mstring+ ") = " << mv[2*dn+1+m];
                output << "\t b_2^("+mstring+ ") = " << mv[3*dn+2+m];
            }
        }
        output << endl;
    }   

    //==============================================================================================
    // pure kinematic terms
    //==============================================================================================
    output << "\n kinematic part" << endl;
    if(fp.beta_order > 0)
    {
        string dn_string = to_string(dn);
        output << "                     " << "\t b_1^("+dn_string+ ") = " << mv[dn+dn];
        
        if(fp.beta_order > 1)
        {
            output << "\t b_0^("+dn_string+ ") = " << mv[2*dn+1+dn];
            output << "\t b_2^("+dn_string+ ") = " << mv[3*dn+2+dn];
        }

        output << endl;
    } 
    output << " " << setfill('-') << setw(95) << "-" << endl;
        
    //==============================================================================================
    output << "\n CNSN expansion part";
    
    for(int r = 1; r < SZM.Get_total_number_of_Te_regions(); r++)
    {
        int is = SZM.Get_index_high()+(r-1)*SZM.Get_number_of_moments_for_one_high_region();
        
        output << "\n " << setfill('-') << setw(42) << "-" 
               << " Region " << to_string(r) << " " << setw(43) << "-" << endl;
        
        for(int m = 0; m < fp.kmax+1; m++)
        {
            string mstring = to_string(m);
            output<< " z^("+mstring+ ") = " << mv[is+m];
            
            if(fp.beta_order > 0)
            {
                output << "\t c_1^("+mstring+ ") = " << mv[is+dn+m];
                
                if(fp.beta_order > 1)
                {
                    output << "\t c_0^("+mstring+ ") = " << mv[is+2*dn+m];
                    output << "\t c_2^("+mstring+ ") = " << mv[is+3*dn+m];
                }
            }
            
            output << endl;
        }   
    }        
    
    output << "\n " << setfill('=') << setw(95) << "=" << endl; 
    output << defaultfloat;
    
    return;
}

void saveAndPrintMomentVector(string fileAddition, vector<double> &mv, Parameters fp, SZ_moment_method SZM){
    output_moment_vector(cout, mv, fp, SZM); //TODO: Add show mess condition
    
    string name=fp.outputPath+"moments."+fileAddition+fp.fileEnding;
    ofstream newfile(name.c_str());

    output_moment_vector(newfile, mv, fp, SZM);
    newfile.close();
}

//==================================================================================================
//
// header for matrix info
//
//==================================================================================================
template<typename somestream>
void header_for_moment_matrix(somestream &output, string description, Parameters fp, SZ_moment_method SZM)
{   
    output << setfill('#') << setw(90) << "#" << endl << "#" << endl;  
    output << "# SZ moment matrix " << description << endl;
    output << "# Number of frequency points: " << fp.gridpoints << " (# of rows)\n";
    output << "# velocity order            : " << fp.beta_order << endl;
    output << "# accuracy level            : " << fp.accuracy_level << endl;
    output << "# kmax+1                    : " << fp.kmax+1 << " (# of Te (!) moments per region)\n";
    output << "# Number of Te regions      : " << SZM.Get_total_number_of_Te_regions() << " \n";
    output << "# total number of moments   : " << SZM.Get_total_number_of_moments() << " (# of cols)\n";
    output << "#\n" << setfill('#') << setw(90) << "#" << endl; 
    
    return;
}

void export_moment_matrix(string description, simpleMatrix &B, Parameters fp, SZ_moment_method SZM)
{
    string filename = fp.outputPath + "matrix." + description + fp.fileEnding;
    ofstream ofile(filename.c_str());
    ofile.precision(16);
    
    header_for_moment_matrix(ofile, description, fp, SZM);
    for(int n=0; n<B.width; n++)
    {
        for(int m=0; m<B.height; m++)
            ofile << B.getValue(m,n) << " ";
        
        ofile << endl;
    }    

    ofile.close();
    return; 
}

void export_vector(vector<double> A, string description, Parameters fp){
    string filename = fp.outputPath + "vector." + description + fp.fileEnding;
    ofstream ofile(filename.c_str());
    ofile.precision(16);
    ofile << "# Output of vector: " << description << endl;
    for(int n=0; n < A.size(); n++) 
        ofile << n << " " << A[n] << endl;
    
    ofile.close();
    return;
}

template<typename somestream>
void output_profile_slice(somestream &output, double (*f)(double l), double lsc, int np)
{
    output.precision(16);
    vector<double> lz(np); 
    init_xarr(1.0e-5, 1, &lz[0], np, 1, 0); 
    
    for(int k=np-1; k>=0; k--) output << " " << -lsc*lz[k] << " " << f(-lz[k]) << endl;
    output << " " << 0.0 << " " << f(0.0) << endl;
    for(int k=0; k<np; k++) output << " " << lsc*lz[k] << " " << f(lz[k]) << endl;
}

void show_profile_slice(double (*f)(double l), double lsc, int np)
{
    output_profile_slice(cout, f, lsc, np);
    return;
}

void export_profile_slice(string fname, double (*f)(double l), double lsc, int np)
{
    ofstream ofile(fname.c_str());
    output_profile_slice(ofile, f, lsc, np);
    ofile.close();
    return;
}

//==================================================================================================
//
// This function gives an example of how the moment method and the class `SZ_moment_method' are used
// to compute the SZ signal. Two variables type for the high temperature moments are available, 
// which should all give practically the same SZ signal. The different steps in the calculation are
// commented below and one can play with things a bit. //TODO: Sort comments!
//
//==================================================================================================
void output_distortion_moments_isothermal(SZ_moment_method &SZM, string mess, Parameters fp)
{
    ofstream ofile;

    SetUpOutput("Using moment method for isothermal cluster to obtain SZ signal"+mess, "SZ_moments_isothermal", fp, ofile);
    //TODO: check output format
    
    //===============================================================================
    //
    // setting up the moments in different form and computing SZ signal
    //
    //===============================================================================
    vector<double> mv;
    vector<vector<double> > DIv(2);
    //===============================================================================
    // Calling the functions of 'SZM' for M and MT;
    //===============================================================================
    compute_isothermal_moments(false, mv, SZM, fp);
    compute_signal_moments(mv, DIv[0], SZM.M, fp, SZM);
    saveAndPrintMomentVector("forM", mv, fp, SZM);

    compute_isothermal_moments(true, mv, SZM, fp);
    compute_signal_moments(mv, DIv[1], SZM.MT, fp, SZM);
    saveAndPrintMomentVector("forMT", mv, fp, SZM);
        
    //===============================================================================
    //
    // output SZ signal
    //
    //===============================================================================
    ofile << "# xcmb \t M distortion \t MT distortion \n";

    for(int n=0; n<(int)fp.xcmb.size(); n++){
        ofile << fp.xcmb[n] << "\t" << DIv[0][n] << "\t" << DIv[1][n] << endl;
    }
    ofile.close();

    mv.clear();
    DIv.clear();
    
    return;
}

//==================================================================================================
//
// In VIK the computation and output are very intertwined. This is all found below.
// A lot of these functions are not in the header.
//TODO: What is VIK exactly?
//TODO: Possibly separate these a little to clean up this file.
//
//==================================================================================================
template<typename somestream>
void outputMomentUpdates(somestream &output, double x_rc, double y_rc, double tau, double ySZ, double TSZ){
    output << " x_rc = " << x_rc << "\t\t y_rc = " << y_rc << endl
           << " tau = " << tau << "\t ySZ = " << ySZ << endl
           << " TSZ = " << TSZ << " keV" << endl;
}

template<typename somestream>
void outputBestFits(somestream &output, double tau, double tau_bestfit, double TSZ, double TSZ_bestfit){
    output << "\n tau = " << tau << "\t tau_bestfit = " << tau_bestfit << endl
           << " TSZ = " << TSZ << "\t\t TSZ_bestfit = " << TSZ_bestfit << endl;
}

void CalculateDIvs(simpleMatrix Mptr, vector<double> &mv, bool Mtransformed, SZ_moment_method SZM,
                   Parameters fp, vector<double> &TSZetc, vector<double> &DIv, vector<double> &DIviso)
{
    double x_rc = PSP.xs, y_rc = PSP.ys;
    compute_signal_moments_integrals(PSP.zsc*PSP.Cluster.Get_rc_cm(), mv, 
                            DIv, TSZetc, 
                            !Mtransformed, 2.0e-6, SZM, fp, Mptr);
                            //If using M we Adjust the temperature, otherwise, we do not
    outputMomentUpdates(cout, x_rc, y_rc, TSZetc[0], TSZetc[1], TSZetc[2]); //TODO: put this behind a show mess variable.

    compute_signal_moments_integrals(x_rc/PSP.zsc, y_rc/PSP.zsc, PSP.zsc*PSP.Cluster.Get_rc_cm(), 
                            1.0e-1/PSP.zsc, mv, DIv, TSZetc, 
                            !Mtransformed, 2.0e-6, SZM, fp, Mptr);
    outputMomentUpdates(cout, x_rc, y_rc, TSZetc[0], TSZetc[1], TSZetc[2]); //TODO: put this behind a show mess variable.
    
    //===========================================================================
    // output moment vector
    //===========================================================================
    fp.Dtau = TSZetc[0];
    fp.Te = TSZetc[2];

    string description = Mtransformed ? "var_MT" : "var_M";
    saveAndPrintMomentVector(description, mv, fp, SZM);
    //===========================================================================
    // compute SZ signal using isothermal moments
    //===========================================================================
    print_message("Now using isothermal approximation");

    compute_isothermal_moments(Mtransformed, mv, SZM, fp);
    compute_signal_moments(mv, DIviso, Mptr, fp, SZM);
    
    saveAndPrintMomentVector("iso."+description, mv, fp, SZM);
}

void CalculateDIvsNonRelAndSmooth(vector<double> &mv, SZ_moment_method SZM, Parameters fp,
                                  vector<double> &DInonrel, vector<double> &DIsmooth, double mv0){
    //=======================================================================
    // compute non-relativistic case
    //=======================================================================
    fill(mv.begin(), mv.end(), 0.0);
    mv[0] = mv0;

    compute_signal_moments(mv, DInonrel, SZM.M, fp, SZM);
    
    //=======================================================================
    // compute using smooth profile approximation
    //=======================================================================  
    vector<double> omegak, dumv;
    
    compute_signal_moments_smooth(PSP.zsc*PSP.Cluster.Get_rc_cm(), omegak, DIsmooth, dumv, 1, fp);    

    string name = fp.outputPath+"moments.omegak"+fp.fileEnding;
    ofstream f0(name.c_str());
    f0.precision(8);
    for(int k=0; k<(int) omegak.size(); k++) { f0 << k << " " << omegak[k] << endl; }
    f0.close();
}

void CalculateDIvsIsothermal(Parameters fp, vector<double> &DInonrel){
    //===========================================================================
    // compute high precision isothermal approximation
    //===========================================================================
    DInonrel.resize(fp.gridpoints);
    for(int k = 0; k < fp.gridpoints; k++){
        DInonrel[k] = pow(fp.xcmb[k], 3)*compute_signal_combo(fp.xcmb[k], fp);
    }
}

//--------------------------------------------------------------------------------------------------
void output_distortion_moments_cluster_VIK(string fileAddition, 
                                              SZ_cluster_profiles &CL, 
                                              double x_rc, double y_rc,
                                              SZ_moment_method &SZM,
                                              Parameters fp, string mess)
{
    //===============================================================================
    PSP.betac = fp.betac;
    PSP.muc = fp.muc;
    PSP.Cluster = CL;

    PSP.xs = x_rc;
    PSP.ys = y_rc;
    PSP.zsc = 1000.0;
    
    /*show_profile_slice(Ne_func, PSP.zsc, 200);
      export_profile_slice(outpath+"Ne_slice"+add, Ne_func, PSP.zsc, 200);
      show_profile_slice(Te_func, PSP.zsc, 200);
      export_profile_slice(outpath+"Te_slice"+add, Te_func, PSP.zsc, 200);*/

    //===============================================================================
    //
    // setting up the moments in different form and computing SZ signal
    //
    //===============================================================================
    vector<double> mv, DIsmooth;
    vector<vector<double> > DIv(2), DIviso(2), DInonrel(2);
    vector<vector<double> > TSZetc(2);

    //===============================================================================
    // Calling the various calculations for the two variable. The results should 
    // be identical to machine precision. For details about the variables see the 
    // moment setup routine above and CSNN 2012.
    //===============================================================================
    CalculateDIvs(SZM.M, mv, false, SZM, fp, TSZetc[0], DIv[0], DIviso[0]);
    CalculateDIvsNonRelAndSmooth(mv, SZM, fp, DInonrel[0], DIsmooth, TSZetc[0][1]);
    CalculateDIvs(SZM.MT, mv, true, SZM, fp, TSZetc[1], DIv[1], DIviso[1]);
    CalculateDIvsIsothermal(fp, DInonrel[1]);

    //===============================================================================
    // find best-fit isothermal model
    //===============================================================================
    print_message("Computing best-fit isothermal model");
    vector<double> solution, Isol;
    find_best_fit_isothermal_model(fp.xcmb, DIviso[1], solution, Isol, 3);
    //TODO: This was originally DIv rather than DIViso, but wasnt converging?
    
    outputBestFits(cout, TSZetc[0][0], solution[0], TSZetc[0][2], solution[1]);

    //===============================================================================
    //
    // output SZ signal
    //
    //===============================================================================
    fp.Dtau = TSZetc[0][0];
    fp.Te = TSZetc[0][1];
    
    ofstream ofile;
    SetUpOutput("Using moments VIK cluster profiles"+mess, fileAddition, fp, ofile, true);
    //TODO: check output format
        
    for(int n=0; n<(int)fp.xcmb.size(); n++)
        ofile << fp.xcmb[n] << " " 
              << DIv[0][n] << " " << DIv[1][n] << " " 
              << DIviso[0][n] << " " << DIviso[1][n] << " " 
              << DInonrel[0][n] << " " << DInonrel[1][n] << " " 
              << DIsmooth[n] << " " << Isol[n] << endl;
    
    ofile.close();

    cout << "Done!\n";
    
    //===============================================================================
    // clean up
    //===============================================================================
    mv.clear();
    DIv.clear();
    DIviso.clear();
    TSZetc.clear();
}

//==================================================================================================
//
// computing temperature moments
//
//==================================================================================================
void output_distortion_cluster_VIK_explicit(string fileAddition, SZ_cluster_profiles &CL,
                                               double x_rc, double y_rc, Parameters fp)
{
    ofstream ofile;
    SetUpOutput("Explicitly calculating line-of-sight-integral", fileAddition, fp, ofile, true);
    //TODO: check output format
    
    //===============================================================================
    PSP.setValues(x_rc, y_rc, 1000.0, fp.betac, fp.muc, CL);
    
    //===============================================================================
    fp.Dtau = fp.Te = 0;
    fp.accuracy_level = fp.Te_max = 0;
    fp.T_order = fp.beta_order = 0;

    vector<double> Ia(fp.gridpoints);

    for(int k = 0; k < fp.gridpoints; k++){
        Ia[k] = pow(fp.xcmb[k], 3)*integrate_SZ_VIK(fp.xcmb[k]);
    }
    
    print_message("Computing best-fit isothermal model");
    vector<double> solution, Isol;
    find_best_fit_isothermal_model(fp.xcmb, Ia, solution, Isol);
    
    double tau, ySZ, TSZ;
    set_tau_ySZ_TSZ_VIK(tau, ySZ, TSZ);
    
    //===============================================================================
    for(int k = 0; k < fp.gridpoints ; k++)
    {
        ofile << fp.xcmb[k] << " " << Ia[k] << " " << Isol[k] << endl;
        cout  << " x= " << fp.xcmb[k] << " " << Ia[k] << " " << Isol[k] << endl;
    }
    
    ofile.close();

    outputBestFits(cout, tau, solution[0], TSZ, solution[1]);
}
//==================================================================================================
template<typename somestream>
void write_yk_ykiso(somestream &f0, somestream &f1, somestream &f2, int kmax)
{
    double tau, ySZ, TSZ;
    set_tau_ySZ_TSZ_VIK(tau, ySZ, TSZ);

    f0 << PSP.xs << " " << tau << " " << TSZ << " " << ySZ << " ";
    f1 << PSP.xs << " " << tau << " " << TSZ << " " << ySZ << " ";
    f2 << PSP.xs << " " << tau << " " << TSZ << " " << ySZ << " ";

    vector<double> rhok(kmax+1, 1.0);
    
    for(int m = 1; m <= kmax; m++)
    {
        double yk = integrate_yk_VIK(m);
        double ykiso = ySZ*pow(TSZ/const_me, m);
        rhok[m] = yk/ykiso;
        
        f0 << yk << " ";
        f1 << rhok[m] << " ";
    }
    
    for(int m=0; m<=kmax; m++)
    {
        double omegak = pow(-1.0, m+1);
        for(int l = 1; l <= m+1; l++){
            omegak += Binomial_coeff(m+1, l)*pow(-1.0, m+1-l)*rhok[l-1];
        }
        f2 << omegak << " ";
    }

    f0 << endl;
    f1 << endl;
    f2 << endl;
}

void output_distortion_cluster_yk(double y_rc,
                                     SZ_cluster_profiles &CL, 
                                     Parameters fp)
{
    //===============================================================================
    PSP.Cluster = CL;
    PSP.ys=y_rc;
    PSP.zsc=1000.0;
    
    //===============================================================================
    vector<double> lz(fp.gridpoints); 
    init_xarr(1.0e-5, 1, &lz[0], fp.gridpoints, 1, 0); 

    fp.outputPath = fp.outputPath + "SZ_Moments_";
    fp.fileEnding = ".y_"+to_string((int) y_rc)+"kpc";

    ofstream ofile_y, ofile_rho, ofile_w;
    SetUpOutput("", "yk", fp, ofile_y);
    SetUpOutput("", "rhok", fp, ofile_rho);
    SetUpOutput("", "wk", fp, ofile_w);
    //TODO: check output format
    
    
/*    for(int k=np-1; k>=0; k--) 
    {
        xs=-PSP.zsc*lz[k];
        write_yk_ykiso(xs, ofile, ofile1, ofile2, kmax);
    }
*/    
    PSP.xs = 0.0;
    write_yk_ykiso(ofile_y, ofile_rho, ofile_w, fp.kmax);
    
    for(int k = 0; k < fp.gridpoints; k++)
    {
        PSP.xs = PSP.zsc*lz[k];
        write_yk_ykiso(ofile_y, ofile_rho, ofile_w, fp.kmax);
    }
    
    ofile_y.close();
    ofile_rho.close();
    ofile_w.close();
}

//==================================================================================================
//
// Compute degeneracy factors
//
//==================================================================================================
void compute_degeneracy_functions(Parameters fp)
{
    ofstream ofile;
    SetUpOutput("Computing the degeneracy coefficients", "SZ_moments_degeneracies.bc_"+ to_string(fp.betac), fp, ofile);
    //TODO: check output format
    
    fp.Dtau = 1.0;
    fp.T = TemperatureIterators(0.01, 50.0, 100, 1);
    
    vector<double> M(3*3, 0.0), A(4*3, 0.0);
        
    gsl_permutation * p = gsl_permutation_alloc (3);
    gsl_vector *b = gsl_vector_alloc (3);
    gsl_vector *x_GSL = gsl_vector_alloc (3);
    
    for(int m = 0; m < fp.T.np; m++)
    {
        fp.updateT(fp.T.Trange[m]);
        ofile << fp.Te << " ";
        
        calculateSignalMatrices(M, A, fp);
        //===============================================================================
        // no velocity vector present
        //===============================================================================
        double det2  =     M[0+3*1]*M[0+3*1] - M[0+3*0]*M[1+3*1];
        double alpha = - ( M[0+3*1]*A[0+4*1] - M[1+3*1]*A[0+4*0] ) / det2; 
        double beta  =   ( M[0+3*1]*A[0+4*0] - M[0+3*0]*A[0+4*1] ) / det2; 
        
        ofile << alpha << " " << beta << " ";

        //===============================================================================
        // velocity vector present and only tau & Te varied
        //===============================================================================
        //TODO: det2 unchanged, should it be?
        alpha = - ( M[0+3*1]*M[2+3*1] - M[1+3*1]*M[2+3*0] ) / det2; 
        beta  =   ( M[0+3*1]*M[2+3*0] - M[0+3*0]*M[2+3*1] ) / det2; 
                          
        ofile << alpha << " " << beta << " ";
        
        //===============================================================================
        // with velocity present
        //===============================================================================
        gsl_matrix_view mm = gsl_matrix_view_array (&M[0], 3, 3);
        int s;
        gsl_linalg_LU_decomp (&mm.matrix, p, &s);
        
        for(int c = 0; c < 4; c++) {
            for(int r = 0; r < 3; r++) { 
                b->data[r]=A[c+4*r];
            }
            gsl_linalg_LU_solve (&mm.matrix, p, b, x_GSL);
            
            for(int k = 0; k < 3; k++) {
                ofile << x_GSL->data[k] << " ";
            }
        }
        ofile << endl;
    }
    
    gsl_permutation_free (p);
    gsl_vector_free (x_GSL);    
    gsl_vector_free (b);    

    ofile.close();
    cout << "done!" << endl;
}

//==================================================================================================
void output_TSZ_bestfit(SZ_cluster_profiles &CL, SZ_moment_method &SZM, Parameters fp, double x_rc_max, int np)
{ //TODO: This is not currently used!
    //===============================================================================
    PSP.Cluster = CL;
    PSP.ys = 0.0;
    PSP.zsc = 1000.0;
    PSP.betac = fp.betac;
    PSP.muc = fp.muc;
    //===============================================================================
    vector<double> lz(np); 
    init_xarr(1.0e-5, x_rc_max/PSP.zsc, &lz[0], np, 1, 0); 
    
    ofstream ofile;
    SetUpOutput("Determining TSZ best fits", "SZ_moments_TSZ.bestfit", fp, ofile);
    //TODO: check output format
    
    double tau, ySZ, TSZ;

    vector<double> Ia(fp.gridpoints), mv, TSZetc;
    vector<double> solution, Isol;
    
    for(int k=0; k<np; k++)
    {
        PSP.xs = PSP.zsc*lz[k];
        set_tau_ySZ_TSZ_VIK(tau, ySZ, TSZ);
        
        ofile << PSP.xs << " " << tau << " " << TSZ << " ";

        //===========================================================================
        //for(int k=0; k<nx; k++)  Ia[k]=pow(fp.xcmb[k], 3)*integrate_SZ_VIK(fp.xcmb[k]);

        compute_signal_moments_integrals(PSP.zsc*PSP.Cluster.Get_rc_cm(), mv, 
                                      Ia, TSZetc, false, 2.0e-6, SZM, parameters, SZM.MT);
        find_best_fit_isothermal_model(fp.xcmb, Ia, solution, Isol);
        
        //===========================================================================
        // estimate best-fit case
        //===========================================================================
        double w1 = integrate_yk_VIK(1)/ySZ/(TSZ/const_me)-1.0;
        double est_tau = tau*( 1.0-w1/( 1.0+0.027*pow(TSZ, 0.86)) );
        double est_Te  = TSZ*( 1.0+w1*exp(-0.026*pow(TSZ, 0.86)) );
        
        ofile << " " << solution[0] << " " << solution[1] << " " << Te_func(lz[k])
              << " " << est_tau << " " << est_Te << endl;
    }
            
    ofile.close();
}
