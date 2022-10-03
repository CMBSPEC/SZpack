//==================================================================================================
//
// SZpack functions unique to the Moment Method Run
//
//==================================================================================================
//
// purpose: computation of the SZ signal according to Chluba, Nagai, Sazonov, Nelson, 2012
//          and Chluba, Switzer, Nagai, Nelson, 2012.
//
// comments: - computations are performed in the single-scattering approximation
//           - polarization effects are neglected
//           - the electron distribution function is assumed to be thermal
//==================================================================================================
//
// Author: Elizabeth Lee
// Based on Work by Jens Chluba
//
//==================================================================================================

#include "SZpackMoment.h"

using namespace std;

//==================================================================================================
//
// compute distortions using temperature-velocity moments
//
//==================================================================================================
void update_mv(vector<double> &mv, Parameters fp, int d, double The_eff, double Nf = 1.0, int i = 0){
    for(int m = 0; m <= fp.kmax; m++) {
        mv[i+m] = fp.Dtau*Nf*pow(The_eff, m);
    }
    
    if(fp.beta_order > 0)
    {
        for(int m = 0; m <= fp.kmax; m++) mv[i+d+m] = fp.Dtau*fp.betac*fp.muc*Nf*pow(The_eff, m);
        
        if(fp.beta_order > 1)
        {          
            double b2 = fp.betac*fp.betac, P2=1.5*fp.muc*fp.muc-0.5;
            for(int m = 0; m <= fp.kmax; m++) mv[i+2*d+m]=fp.Dtau*b2   *Nf*pow(The_eff, m);
            for(int m = 0; m <= fp.kmax; m++) mv[i+3*d+m]=fp.Dtau*b2*P2*Nf*pow(The_eff, m);
        }
    }
}

void compute_isothermal_moments(bool Mtransformed, vector<double> &mv, SZ_moment_method &SZM, Parameters fp)
{
    //===============================================================================
    // setting up the moments in different forms
    //===============================================================================
    mv.resize(SZM.Get_total_number_of_moments(), 0.0);
    
    //
    double The_eff = fp.calc.The;
    int i = 0;
    double Nf = 1.0;

    int d = 2 + fp.kmax;
    int ir = SZM.Get_index_of_Te_region(fp.Te);

    //===============================================================================
    // pure velocity moments
    //===============================================================================
    if(fp.beta_order > 0)
    {
        mv[d-1]=fp.Dtau*fp.betac*fp.muc;
        
        if(fp.beta_order > 1)
        {          
            mv[2*d-1] = fp.Dtau*fp.betac*fp.betac;
            mv[3*d-1] = fp.Dtau*fp.betac*fp.betac*(1.5*fp.muc*fp.muc-0.5);
        }
    }

    //===============================================================================
    // For low temperature moments (isothermal), ir == 0. This is all correct.
    // For high temperature moments, we must change our position in the matrix.
    //===============================================================================
    if(ir > 0)
    {
        // lookup where the moments have to be written in the moment vector
        i = SZM.Get_start_index_of_moment_set(ir);
        Nf = norm_df_RM_dTheta(The_eff)/The_eff;
        //TODO: This /The_eff scaling is because I changed this function many many commits before

        // (var 0)
        if(!Mtransformed)
        {
            The_eff -= SZM.Get_The_ref(ir);
        }
        
        d -= 1;
    }
    
    update_mv(mv, fp, d, The_eff, Nf, i);

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
void compute_signal_moments(const vector<double> &mv, vector<double> &DIv, simpleMatrix &Mptr,
                            Parameters fp, SZ_moment_method SZM)
{
    int tot_moments = SZM.Get_total_number_of_moments();
    if(tot_moments !=(int)mv.size() || Mptr.matrix.size() != tot_moments*fp.gridpoints)
    { exit_error("please check moment vector or moment matrix"); }

    DIv.resize(fp.gridpoints);
    
    //==============================================================================================
    // S = F * m like in Chluba et. al, 2012
    //==============================================================================================
    for(int n = 0; n < fp.gridpoints; n++)
    {
        DIv[n] = 0.0;
        
        for(int m = 0; m < tot_moments; m++){
            DIv[n] += Mptr.getValue(n,m)*mv[m];
        }
    }
}

//==================================================================================================
// 
// compute moments over the cluster profiles.
//
//==================================================================================================
void populate_mv(int r, double zsc, vector<double> &mv, Parameters fp, SZ_moment_method SZM, int k_offset = 0)
{
    double yk, bk;
    bool extraTerms = (r == 0);

    for(int k = k_offset; k <= fp.kmax+k_offset; k++)
    {
        // k-1 because of difference between y^(k)<-->z^(k) variable
        yk = compute_y_k_cond(k, r, zsc, SZM); 
        mv.push_back(yk);
    }
    
    //======================================================================================
    if(fp.beta_order > 0)
    {
        if (extraTerms){
            bk = compute_b_1(zsc);
            mv.push_back(bk);
        } 

        for(int k = k_offset; k <= fp.kmax+k_offset; k++)
        {
            bk=compute_b_k_1_cond(k+1, r, zsc, SZM);
            mv.push_back(bk);
        }
        
        if(fp.beta_order > 1)
        {
            if (extraTerms){
                bk = compute_b_0(zsc);
                mv.push_back(bk);
            } 
        
            for(int k = k_offset; k <= fp.kmax+k_offset; k++)
            {
                bk=compute_b_k_0_cond(k+1, r, zsc, SZM);
                mv.push_back(bk);
            }

            if (extraTerms){
                bk = compute_b_2(zsc);
                mv.push_back(bk);
            } 

            for(int k = k_offset; k <= fp.kmax+k_offset; k++)
            {
                bk=compute_b_k_2_cond(k+1, r, zsc, SZM);
                mv.push_back(bk);
            }
        }
    }
}

double N_func_1(double The){ return 1.0; }
void compute_cluster_moments(double zsc, vector<double> &mv, bool AdjustTemperature, double lres,
                             SZ_moment_method SZM, Parameters &fp)
{
    mv.clear();

    //==============================================================================================
    // brute force temperature structure computation; This has to be done to determine the intervals
    // in l=z/zsc that correspond to the different temperature regions. Temperature structures that 
    // are smaller than zsc/nres will not be resolved.
    //==============================================================================================
    SZM.determine_Te_structure(lres);
    
    //==============================================================================================
    // low temperature moments
    //==============================================================================================
    setNormFunc(N_func_1);
    The_ref_global = 0.0;

    fp.kmax = SZM.Get_kmax();

    populate_mv(0, zsc, mv, fp, SZM);
    //==============================================================================================
    // high temperature moments
    //==============================================================================================
    int nTregs = SZM.Get_total_number_of_Te_regions();
    if(nTregs > 1)
    {
        setNormFunc(norm_df_RM_dTheta);
        
        for(int r=1; r<nTregs; r++)
        {
            The_ref_global = AdjustTemperature ? SZM.Get_The_ref(r) : 0.0;

            populate_mv(r, zsc, mv, fp, SZM, -1); // k-1 because of difference between y^(k)<-->z^(k) variable
        }
    }
}


void compute_cluster_moments_3D(double lx0, double ly0, double lsc, double lang, vector<double> &mv,
                                bool AdjustTemperature, double lres, SZ_moment_method SZM, Parameters &fp)
{
/*    mv.clear();
    
    vector<double> dum_mv; //TODO: why is this?*/
    
    int np=10;
    
    for(int c=0; c<np; c++)
    {
//        IValues.lx = la[c];
        IValues.lx = 1.0/lsc;
        
        for(int r=0; r<np; r++)
        {
//            IValues.ly = la[r];
            IValues.ly =1.0/lsc;
            
            compute_cluster_moments(lsc, mv, AdjustTemperature, lres, SZM, fp);
        }
    }
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
void signal_moment_setUp(vector<double> &mv, vector<double> &DIv, vector<double> &TSZetc,
                         double zsc, double lx0 = 0.0, double ly0 = 0.0, double lang = 0.0)
{
    //==============================================================================================
    // clean plates & simple setup
    //==============================================================================================
    mv.clear();
    DIv.clear();
    TSZetc.clear();

    // simple calculation of standard variables
    double tau, ySZ;
    if (IValues.oneD) {
        tau = compute_tau(zsc);
        ySZ = compute_ySZ(zsc);
    }
    else {
        tau = compute_Int_3D(lx0, ly0, zsc, lang, compute_tau);
        ySZ = compute_Int_3D(lx0, ly0, zsc, lang, compute_ySZ);
    }
    double TSZ = ySZ*const_me/tau;
    
    // return values
    TSZetc.push_back(tau);
    TSZetc.push_back(ySZ);
    TSZetc.push_back(TSZ);
}

void compute_signal_moments_integrals(double zsc, vector<double> &mv, vector<double> &DIv, vector<double> &TSZetc,
                            bool AdjustTemperature, double lres, SZ_moment_method SZM, Parameters &fp, simpleMatrix &Mptr)
{
    //This inherits ClusterF_1D. Changes should be made to this global if they're needed.
    IValues = integrationValues(true);
    signal_moment_setUp(mv, DIv, TSZetc, zsc);
    
    //==============================================================================================
    // compute the signal
    //==============================================================================================
    compute_cluster_moments(zsc, mv, AdjustTemperature, lres, SZM, fp);
    compute_signal_moments(mv, DIv, Mptr, fp, SZM);
}

//==================================================================================================
void compute_signal_moments_integrals(double lx0, double ly0, double lsc, double lang, vector<double> &mv, 
                               vector<double> &DIv, vector<double> &TSZetc, bool AdjustTemperature,
                               double lres, SZ_moment_method SZM, Parameters &fp, simpleMatrix &Mptr)
{
    //This inherits ClusterF_3D. Changes should be made to this global if they're needed. 
    IValues = integrationValues(false);
    signal_moment_setUp(mv, DIv, TSZetc, lsc, lx0, ly0, lang);
    //==============================================================================================
    // compute the signal
    //==============================================================================================
    compute_cluster_moments_3D(lx0, ly0, lsc, lang, mv, AdjustTemperature, lres, SZM, fp);
    wait_f_r("STOP");
    compute_signal_moments(mv, DIv, Mptr, fp, SZM);
}

//==============================================================================================
//
// Compute the SZ signal using temperature-velocity moments, but assuming that the variance of
// the temperature and velocity field along the line-of-sight is small. In this case the SZ 
// signal is related to S_iso(tau, TeSZ) and its derivatives with respect to Te. The functions 
// and parameters are similar to those above.
//
//==============================================================================================
void compute_signal_moments_smooth(const vector<double> &omegav,
                                                        const vector<double> &TSZetc,
                                                        vector<double> &DIv, 
                                                        Parameters fp)
{
    DIv.resize(fp.gridpoints);
    
    Parameters temp = Parameters();
    temp.copyParameters(fp);
    temp.Te = TSZetc[2];
    temp.Dtau = TSZetc[0];
    temp.betac = temp.betao = 0.0;
    temp.muc = temp.muo = 1.0;
    temp.setCalcValues();
    
    double eps=0.005;
    double TeP = temp.Te*(1.0+eps);
    double TeL = temp.Te*(1.0-eps);
    double TePP = temp.Te*(1.0+2*eps);
    double TeLL = temp.Te*(1.0-2*eps);
    
    //==============================================================================================
    // S = F * m like in Chluba et. al, 2012
    //==============================================================================================
    double d0, dp, dm, dp2, dm2;
    double d2Sd2T, d3Sd3T, d4Sd4T;
    
    for(int n=0; n < temp.gridpoints; n++)
    {
        d0=compute_signal_combo(n, temp);
        
        DIv[n]=d0;
        
        if(fp.kmax > 0)
        {
            temp.Te = TeP;
            dp=compute_signal_combo(n, temp);
            temp.Te = TeL;
            dm=compute_signal_combo(n, temp);
            d2Sd2T=(dp-2.0*d0+dm)/pow(eps, 2) / (2.0);
            
            DIv[n]+=omegav[1]*d2Sd2T;
            
            if(fp.kmax > 1)
            {
                temp.Te = TePP;
                dp2=compute_signal_combo(n, temp);
                temp.Te = TeLL;
                dm2=compute_signal_combo(n, temp);
                d3Sd3T=(dp2-2.0*dp+2.0*dm-dm2)/2.0/pow(eps, 3) / (6.0);
                
                DIv[n]+=omegav[2]*d3Sd3T;

                if(fp.kmax > 2)
                {
                    d4Sd4T=(dp2-4.0*dp+6.0*d0-4.0*dm+dm2)/pow(eps, 4) / (24.0);

                    DIv[n]+=omegav[3]*d4Sd4T;
                }
            }
        }
        DIv[n]*=pow(fp.xcmb[n], 3);
    }
}

//==================================================================================================
void compute_signal_moments_smooth(double zsc, vector<double> &omegav, vector<double> &DIv,
                                   vector<double> &TSZetc, int kmax, Parameters fp)
{
    IValues.oneD = true;
    signal_moment_setUp(omegav, DIv, TSZetc, zsc);
    
    //==============================================================================================
    // compute omega(k)
    //============================================================================================== 
    vector<double> rhok(kmax+1);
    omegav.resize(kmax+1);
    omegav[0]=0.0;
    rhok[0]=1.0;

    for(int k=1; k<=kmax; k++){
        rhok[k] = compute_yk(zsc, k)/TSZetc[1]/pow(TSZetc[2]/const_me, k);
    }
    
    for(int m=1; m<=kmax; m++)
    {
        omegav[m]=pow(-1.0, m+1);
        
        for(int l=1; l<=m+1; l++){
            omegav[m]+= Binomial_coeff(m+1, l)*pow(-1.0, m+1-l)*rhok[l-1];
        }
    }
    
    //==============================================================================================
    // compute the signal
    //==============================================================================================    
    compute_signal_moments_smooth(omegav, TSZetc, DIv, fp);
}    

//==================================================================================================
// For Computation of degeneracy functions
//==================================================================================================
void calculateSignalMatrices(vector<double> &M, vector<double> &A, Parameters &fp){
    vector< vector<double> > a_v(3, vector<double>(fp.gridpoints, 0));
    vector< vector<double> > s_v(4, vector<double>(fp.gridpoints, 0));

    for(int k=0; k<fp.gridpoints; k++)
    {
        double x3 = pow(fp.xcmb[k], 3);
        
        fp.D.setValues(2,0,0);
        Dcompute_signal_combo_CMB(fp.xcmb[k], fp);
        a_v[0][k] = fp.D.dDn_dThe[0]*x3;
        a_v[1][k] = fp.D.dDn_dThe[1]*x3;
        //
        s_v[0][k] = fp.D.dDn_dThe[2]*x3; 
        
        fp.D.setValues(1,1,0);
        Dcompute_signal_combo_CMB(fp.xcmb[k], fp);
        a_v[2][k] = fp.D.dDn_dThe[0]*x3*0.01;
        //
        s_v[1][k] = fp.D.dDn_dThe[1]*x3*0.01;

        fp.D.setValues(0,2,0);
        Dcompute_signal_combo_CMB(fp.xcmb[k], fp);
        s_v[2][k] = fp.D.dDn_dThe[0]*x3*0.01*0.01;

        fp.D.setValues(0,0,1);
        Dcompute_signal_combo_CMB(fp.xcmb[k], fp);
        s_v[3][k] = fp.D.dDn_dThe[0]*x3*0.01*0.01;
    }
    
    //===============================================================================
    // compute signal matrices
    //===============================================================================
    for (int r = 0; r < 3; r++){
        for (int c = 0; c < 3; c++){
            M[c+3*r] = DC_sumprod(a_v[c], a_v[r]);
            A[c+4*r] = DC_sumprod(s_v[c], a_v[r]);
        }
        A[3+4*r] = DC_sumprod(s_v[3], a_v[r]);
    }

    a_v.clear();
    s_v.clear();
}

//==================================================================================================
// Functions for VIK
//==================================================================================================
void set_tau_ySZ_TSZ_VIK(double &tau, double &ySZ, double &TSZ){
    tau = integrate_yk_VIK(-1);
    ySZ = integrate_yk_VIK(0);
    TSZ = ySZ*const_me/tau;
}
//==================================================================================================
//==================================================================================================
