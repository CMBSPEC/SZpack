//==================================================================================================
//
// Computation of the thermal and k-SZ effect using explicit integration collision term.
// Using a precomputed kernel over the angular directions. 
// TODO: Update this!!!!
// The cluster is assumed to be isothermal and moving at a given speed betac with direction muc 
// relative to the line of sight. The integrals are carried out in the cluster frame.
//
//==================================================================================================
//
// Author: Elizabeth Lee
// Based on work by Jens Chluba (CITA, University of Toronto)
//
// first implementation: October 2018
// last modification: March 2021
//
//==================================================================================================

#include "SZ_multipoleKernel.h"
#include <cmath>
#include "global_functions.h"


//==================================================================================================
// Multipole Kernel
//==================================================================================================
MultipoleKernel::MultipoleKernel()
    : MultipoleKernel(0, 0.0, 1.0, 1.0e-4) {}

MultipoleKernel::MultipoleKernel(int l_i, double s_i, double eta_i)
    : MultipoleKernel(l_i, s_i, eta_i, 1.0e-4) {}

MultipoleKernel::MultipoleKernel(int l_i, double s_i, double eta_i, double Int_eps_i){
    Update_l(l_i);
    Update_s(s_i);
    eta = eta_i;
    Int_eps = Int_eps_i;
    photon_anisotropy = true;
    gamma0 = sqrt(1+eta*eta);
    s_low = log((gamma0-eta)/(gamma0+eta));
    s_high = log((gamma0+eta)/(gamma0-eta));
    Lpart = mus = mup = r = 0;
    P.resize(l+3);
    h0.resize(l+1);
    h2.resize(l+1);
    h4.resize(l+1);
    K1.resize(l+1);
    K2.resize(l+1);
}

vector<double> MultipoleKernel::s_limits(){
    vector<double> lims;
    lims.resize(2);
    lims[0] = s_low;
    lims[1] = s_high;
    return lims;
}

void MultipoleKernel::Update_s(double s_i){
    s = s_i;
    t = exp(s);
}

void MultipoleKernel::Update_l(int l_i){
    l = l_i;
    if (l < 0){
        l = 0;
        print_error("l is out of bounds, it must be >=0. l is now set to 0.");
    }
}

void MultipoleKernel::Calculate_integral_variables(){
    //Standard definitions
    //==============================================================================================
    // mu == cosine of angle of gamma  and beta 
    // mup== cosine of angle of gamma' and beta
    // mus== cosine of angle of gamma  and gamma'
    //==============================================================================================    
    double gamma2=1.0+eta*eta;
    //
    double kappap=(gamma0-eta*mup)/gamma0;
    double zeta=1.0/(t*pow(gamma0-eta*mup,2));
    double alpha_sc=1.0-mus;

    //d2sigma/dmu/dmup as defined equation 2 CNSN
    double dsig = 3.0/8.0/PI*t/kappap/gamma2*(1.0-zeta*alpha_sc*(1.0-0.5*zeta*alpha_sc));
    double dphidt = (gamma0-eta*mup)/sqrt(eta*eta*(-1+mup*mup)*(-1+mus*mus)-pow(gamma0-t*gamma0+eta*mup*(t-mus),2));
    double Plmu = 0;
    if (photon_anisotropy) {
        for (int k=0; k<l+1; k++){
            Plmu += pow(-1,k)*Binomial(l,k)*Binomial(l+k,k)*pow((1-mus)/2,k);
        }
    }
    else {
        for (int k=0; k<l+1; k++){
            Plmu += pow(-1,k)*Binomial(l,k)*Binomial(l+k,k)*pow((1-mup)/2,k);
        }
    }
    r = Plmu*dphidt*dsig;
}

double MultipoleKernel::sigl_Boltzmann_Compton(double int_mup){
    mup = int_mup;
    Calculate_integral_variables();
    return r;
}

double MultipoleKernel::mup_Int(double mus_int){
    mus = mus_int;
    double limitvar = sqrt((-1+mus*mus)*((-1+t)*(-1+t)*gamma0*gamma0-eta*eta*(1+t*t-2*t*mus)));
    double a=((-1 + t)*gamma0*(t-mus)-limitvar)/(eta*(1+t*t-2*t*mus));
    double b=((-1 + t)*gamma0*(t-mus)+limitvar)/(eta*(1+t*t-2*t*mus));
    double epsrel=Int_eps*0.8, epsabs=1.0e-300;
    
    return Integrate_using_Patterson_adaptive(a, b, epsrel, epsabs, [this](double int_var) { return this->sigl_Boltzmann_Compton(int_var);});
}

double MultipoleKernel::mus_Int(){
    double a=-1.0;
    double b=(eta*eta*(1+t*t)-pow(-1+t,2)*gamma0*gamma0)/(2*eta*eta*t);
    double epsrel=Int_eps*0.9, epsabs=1.0e-300;

    return Integrate_using_Patterson_adaptive(a, b, epsrel, epsabs, [this](double int_var) { return this->mup_Int(int_var);});
}

void MultipoleKernel::Calculate_formula_variables(){
    Lpart = 0.5*(fabs(s)-2.0*asinh(eta));
    for (int m = 0; m < l+3; m++){
        double ppart = pow(eta,2*m+1)/pow(1+eta*eta,2.5);
        double tpart = fabs(pow(t-1,2*m+1))/pow(4*t,m-2)/pow(1+t,5.0);
        P[m] = ppart-tpart;
    }

    h0[0] = -(2*P[2]+5*P[1])/15.0;
    h2[0] = -(P[2]+P[1])/3.0;
    h4[0] = Lpart+P[2]+2*P[1]+P[0];

    if (l >= 1){
        h0[1] = P[2]/5.0;
        h2[1] = -(3*Lpart+4*P[2]+7*P[1]+3*P[0])/3.0;
        h4[1] = (3*Lpart+P[3]+5*P[2]+7*P[1]+3*P[0])/2.0;

        if (l >= 2){
            h0[2] = Lpart+(23*P[2]+35*P[1]+15*P[0])/15.0;
            h2[2] = (-5*h0[2]-P[3])/2.0;
            h4[2] = (-3*h2[2]-P[4]-P[3])/4.0;
            
            if (l >= 3){
                for (int m = 3; m<l+1; m++){
                    h0[m] = ((2*m+1)*h0[m-1]-pow(-1,m)*P[m])/(2*(m-2));
                    h2[m] = (-5*h0[m]-pow(-1,m)*P[m+1])/(2*(m-1));
                    h4[m] = (15*h0[m]-pow(-1,m)*(2*m-5)*P[m+1]-pow(-1,m)*2*(m-1)*P[m+2])/(4*m*(m-1));
                }
            }
        }
    }

    for (int m = 0; m<l+1; m++){
        //Calculating K1
        double invpref = eta*gamma0*pow(4*t,m);
        double ksum = 0;
        for (int k = 0; k < m+1; k++){
            double pref = Binomial(m,k)*pow(-1,k)*pow(t-1,2*(m-k))/(2*k+1);
            ksum += pref*(pow(t+1,2*k+1)-(pow(gamma0*fabs(t-1)/eta,2*k+1)));
        }
        K1[m] = ksum/invpref;

        //Calculating K2
        K2[m] = 5*(1+eta*eta)*(1+t)*(1+t)*h0[m];
        K2[m] += -((5+2*eta*eta)*(1+t*t)+(22+16*eta*eta)*t)*h2[m];
        K2[m] += 4*(3+2*eta*eta)*t*h4[m];
    }
}

double MultipoleKernel::Calculate_electron_multipoles(){
    if (s<s_low|| s>s_high) { return 0.0; }

    Lpart = 0.5*(fabs(s)-2.0*asinh(eta));
    double beta0 = eta/gamma0;
    vector<double> K(l+5,0.0);
    for (int k = 0; k< l+5; k++){
        int m = k-3;
        double part1 = (pow(t,m)+1)*(pow(1+beta0,m)-pow(1-beta0,m));
        double part2 = (pow(t,m)-1)*(pow(1+beta0,m)+pow(1-beta0,m));
        double signt = t<1 ? -1.0 : 1.0;
        K[k] = (part1-signt*part2)/(2*pow(1-beta0,m)*pow(t,m));
    }

    vector<double> X0(l+1,0.0), X1(l+1,0.0), X2(l+1,0.0), X3(l+1,0.0), X4(l+1,0.0);
    for (int m = 0; m < l+1; m++){
        for (int k = 0; k < m+1; k++){
            double factor = Binomial(m,k)*pow(-1,k);
            X0[m] += factor*K[k+4]/(k+1.0);
            X1[m] += (k==0) ? -2.0*factor*Lpart : factor*K[k+3]/k;
            X2[m] += (k==1) ? -2.0*factor*Lpart : factor*K[k+2]/(k-1.0);
            X3[m] += (k==2) ? -2.0*factor*Lpart : factor*K[k+1]/(k-2.0);
            X4[m] += (k==3) ? -2.0*factor*Lpart : factor*K[k]/(k-3.0);
        }
        double pref = pow((gamma0-eta)/eta,m+1)/pow(2.0,m);
        X0[m] *= pref;
        X1[m] *= pref/(gamma0-eta);
        X2[m] *= pref/pow(gamma0-eta,2);
        X3[m] *= pref/pow(gamma0-eta,3);
        X4[m] *= pref/pow(gamma0-eta,4);
    }

    vector<double> Pn(l+1,0.0);
    for (int m = 0; m<l+1; m++){
        Pn[m] = 3*X4[m]/gamma0-6*(1+t)*X3[m]+(3*(1+4*t+t*t)+2*eta*eta*(1+6*t+t*t))*X2[m]/gamma0;
        Pn[m] += -2*(3+2*eta*eta)*t*(1+t)*X1[m]+(3+4*eta*eta+4*pow(eta,4))*t*t*X0[m]/gamma0;
    }

    double prefactor = 3/(32*pow(eta,5));
    double lsum = 0;
    for (int m = 0; m<l+1; m++){
        double pref = Binomial(l,m)*Binomial(l+m,m);
        lsum += pref*Pn[m];
    }
    return prefactor*lsum;
}

//TODO: The tlimits should probably be implemented here
double MultipoleKernel::Calculate_integrated(){
    if (s<s_low|| s>s_high) { return 0.0; }
    return mus_Int()*t;
}

double MultipoleKernel::Calculate_formula(){
    if (photon_anisotropy) {
        if (s<s_low|| s>s_high) { return 0.0; }
        Calculate_formula_variables();
        double prefactor = 3/(32*pow(eta,6));
        double lsum = 0;
        for (int n = 0; n<l+1; n++){
            double pref = Binomial(l,n)*Binomial(l+n,n);
            lsum += pref*(4*pow(eta,6)*K1[n]+(1+t)*K2[n]/t/pow(eta,2*n));
        }
        return prefactor*lsum*t;
    }
    return Calculate_electron_multipoles();
}

double MultipoleKernel::Calculate_stable(){
    if (s<s_low|| s>s_high) { return 0.0; }
    switch (l) {
    case 0:
        if (eta >= 0.0009690){return Calculate_formula();}
        return Calculate_integrated();
    case 1:
        if (eta >= 0.01749){return Calculate_formula();}
        return Calculate_integrated();
    case 2:
        if (eta >= 0.06710){return Calculate_formula();}
        return Calculate_integrated();
    case 3:
        if (eta >= 0.1467){return Calculate_formula();}
        return Calculate_integrated();
    case 4:
        if (eta >= 0.2426){return Calculate_formula();}
        return Calculate_integrated();
    case 5:
        if (eta >= 0.3499){return Calculate_formula();}
        return Calculate_integrated();
    case 6:
        if (eta >= 0.4566){return Calculate_formula();}
        return Calculate_integrated();
    case 7:
        if (eta >= 0.5680){return Calculate_formula();}
        return Calculate_integrated();
    case 8:
        if (eta >= 0.7208){return Calculate_formula();}
        return Calculate_integrated();
    case 9:
        if (eta >= 0.8667){return Calculate_formula();}
        return Calculate_integrated();
    default:
        return Calculate_integrated();
    }
}

//==================================================================================================
// Beam Kernel
//==================================================================================================
BeamKernel::BeamKernel()
    : BeamKernel(0, 0.0, 1.0, 0.0) {}

BeamKernel::BeamKernel(int l_i, double s_i, double eta_i, double mup_i){
    Update_l(l_i);
    Update_s(s_i);
    eta = eta_i;
    mup = mup_i;
    gamma0 = sqrt(1+eta*eta);
    f0 = f1 = f2 = 0.0;
}

vector<double> BeamKernel::s_limits(){
    vector<double> lims;
    lims.resize(2);
    lims[0] = log((gamma0-eta)/(gamma0-eta*mup));
    lims[1] = log((gamma0+eta)/(gamma0-eta*mup));
    return lims;
}

vector<double> BeamKernel::mup_limits(){
    vector<double> lims;
    lims.resize(2);
    double temp = (gamma0*(t-1)-eta)/(eta*t);
    lims[0] = (temp < -1.0) ? -1.0 : temp;
    temp = (gamma0*(t-1)+eta)/(eta*t);
    lims[1] = (temp > 1.0) ? 1.0 : temp;
    return lims;
}

void BeamKernel::Update_s(double s_i){
    s = s_i;
    t = exp(s);
}

void BeamKernel::Update_l(int l_i){
    l = l_i;
    if (l < 0){
        l = 0;
        print_error("l is out of bounds, it must be >=0. l is now set to 0.");
    }
}

void BeamKernel::Calculate_monopole(){
    double denominator = (32*(pow(eta,3.0))*gamma0*t*pow(gamma0 - eta*mup,4.0))/3.0;
    double part1a = (t-1)*(t-1)*(3*mup*mup - 1)+(eta*eta)*(2-2*t+3*t*t +(2-18*t+14*t*t)*mup*mup + 3*(t*t)*pow(mup,4.0));
    double part1b = 4*pow(eta,4.0)*t*(-1+2*t+(-3+9*t)*mup*mup +t*pow(mup,4.0))+4*pow(eta,6.0)*(t*t)*(1+6*mup*mup + pow(mup,4.0));
    double part2 = -2*eta*gamma0*mup*((t-1)*(-2+t+3*t*mup*mup)- 2*(eta*eta)*t*(3-5*t+(1-3*t)*mup*mup)+8*pow(eta,4.0)*(t*t)*(1+mup*mup));
    f0 = (part1a+part1b+part2)/denominator;
}

void BeamKernel::Calculate_dipole(){
    double denominator = (32*pow(eta,4.0)*gamma0*t*pow(gamma0 - eta*mup,4.0))/3.0;
    double part1a = eta*(t-1)*(t-1) *(2-2*t-3*(2+t)*mup*mup +15*t*pow(mup,4.0));
    double part1b = eta*(eta*eta)*(2*t*(-3+5*t-2*t*t)+(-4+24*t-26*t*t+11*t*t*t)*mup*mup +6*t*(3-12*t+8*t*t)*pow(mup,4.0) +5*(t*t*t)*pow(mup,6.0));
    double part1c = eta*2*pow(eta,4.0)*t*(t*(2-t)+(6-22*t+17*t*t)*mup*mup +(2-28*t+37*t*t)*pow(mup,4.0) + 3*t*t*pow(mup,6.0));
    double part1d = eta*4*pow(eta,6.0)*(t*t)*(mup*mup)*(-4+5*t+2*(-2+5*t)*mup*mup +t*pow(mup,4.0));
    double part2a = -gamma0*mup*(pow(t-1,3.0)*(5*mup*mup -3)+(eta*eta)*(t-1)*(2+10*t-7*(t*t)+(2-34*t+20*t*t)*mup*mup +15*(t*t)*pow(mup,4.0)));
    double part2b = -gamma0*mup*4*pow(eta,4.0)*t*(1+(-1+t)*(-3+14*t)*mup*mup + 3*t*(-1+2*t)*pow(mup,4.0));
    double part2c = -gamma0*mup*(-4*pow(eta,6.0)*(t*t)*(1-t+(6-10*t)*mup*mup+(1-5*t)*pow(mup,4.0)));
    f1 = (part1a+part1b+part1c+part1d+part2a+part2b+part2c)/denominator;
}

void BeamKernel::Calculate_quadrupole(){
    double denominator = (32*pow(eta,5.0)*gamma0*t*pow(gamma0 - eta*mup,4.0))/3.0; 
    double part1a = 3*pow(-1+t,4.0)*(3-30*mup*mup +35*pow(mup,4.0));
    double part1b = 2*(eta*eta)*(-1+t)*(-1+t) *((-4-6*t -3*t*t)-3*(8-108*t+69*t*t)*mup*mup +15*(4-34*t+9*t*t)*pow(mup,4.0) +315*(t*t)*pow(mup,6.0));
    double part1c = pow(eta,4.0)*((-8+56*t-132*t*t +156*t*t*t -63*pow(t,4.0))-4*(-4-84*t+437*t*t -561*t*t*t +225*pow(t,4.0))*mup*mup);
    double part1d = pow(eta,4.0)*(6*(4-156*t+574*t*t-586*t*t*t+189*pow(t,4.0))*pow(mup,4.0)+60*t*t*(15-47*t+29*t*t)*pow(mup,6.0)+105*pow(t,4.0)*pow(mup,8.0));
    double part1e = 8*pow(eta,6.0)*t*(-(-2+3*t)*(1-5*t+3*t*t)-3*t*(23-61*t+36*t*t)*mup*mup);
    double part1f = 8*pow(eta,6.0)*t*((-18+229*t-453*t*t +210*t*t*t)*pow(mup,4.0)+3*t*(15-93*t+92*t*t)*pow(mup,6.0)+15*(t*t*t)*pow(mup,8.0)) ;
    double part1g = 8*pow(eta,6.0)*(t*t)*(-1+3*mup*mup)*((2-6*t+3*t*t)+3*(4-20*t+15*t*t)*mup*mup +(2-30*t+45*t*t)*pow(mup,4.0) + 3*(t*t)*pow(mup,6.0));
    double part2a = -4*eta*gamma0*mup*(3*pow(t-1,3.0)*(6-3*t-10*(1+2*t)*mup*mup +35*t*pow(mup,4.0)));
    double part2b = -4*eta*gamma0*mup*(eta*eta *(t-1)*((4-74*t+120*t*t -63*t*t*t)-3*(4-32*t-28*t*t +31*t*t*t)*mup*mup));
    double part2c = -4*eta*gamma0*mup*(eta*eta *(t-1)*(45*t*(2-12*t +7*t*t)*pow(mup,4.0) + 105*(t*t*t)*pow(mup,6.0)));
    double part2d = -4*eta*gamma0*mup*(2*pow(eta,4.0)*t*((6-58*t+99*t*t -45*t*t*t)-(16-72*t+27*t*t +27*t*t*t)*mup*mup));
    double part2e = -4*eta*gamma0*mup*(2*pow(eta,4.0)*t*(3*(-2+54*t-153*t*t+95*t*t*t)*pow(mup,4.0) +15*t*t*(-3+5*t)*pow(mup,6.0))); 
    double part2f = -4*eta*gamma0*mup*(4*pow(eta,4.0)*(t*t)*(3*mup*mup-1)*((-4+3*t)*(-1+3*t)+2*(2-15*t+15*t*t)*mup*mup +3*t*(-1+3*t)*pow(mup,4.0)));
    f2 = (part1a+part1b+part1c+part1d+part1e+part1f+part1g+part2a+part2b+part2c+part2d+part2e+part2f)/denominator;
}

double BeamKernel::Calculate_formula(){
    if (l == 1){ 
        Calculate_dipole();
        return f1;
    }
    if (l == 2){ 
        Calculate_quadrupole();
        return f2;
    }
    Calculate_monopole();
    return f0;
}


//==================================================================================================
// Integration Class
//==================================================================================================
IntegralKernel::IntegralKernel()
    : IntegralKernel(0.0, 0.0, 0.0, 1.0e-4) {}

IntegralKernel::IntegralKernel(double x_i, double betac_i, double muc_i, double eps_Int_i){
    Int_eps = eps_Int_i;
    betac = betac_i;
    muc = muc_i;
    x = x_i;
    eta=s=xp=0.0;
    l = 0;
    run_mode = "";
    etaDistribution = plainDistribution;
    MK = MultipoleKernel(l, s, eta, Int_eps);
    BK = BeamKernel(l, s, eta, 0.0);
    xfac = 1.0;
    beam_kernel=fixed_eta=false;
    photon_anisotropy = true;
}

IntegralKernel::IntegralKernel(double x_i, Parameters fp)
    : IntegralKernel(fp.calc.xfac*x_i, fp.betac, fp.calc.mucc, fp.relative_accuracy) {
        xfac = fp.calc.xfac;
    }

void IntegralKernel::Update_x(double x_i){
    x = xfac*x_i;
}

void IntegralKernel::Calculate_shared_variables(){
    xp = x*exp(s);
    MK = MultipoleKernel(l, s, eta, Int_eps);
    BK = BeamKernel(l, s, eta, mup);
    MK.photon_anisotropy = photon_anisotropy;
}

double IntegralKernel::Calculate_kernel(int l_i){
    if (beam_kernel){
        BK.Update_l(l_i);
        return BK.Calculate_formula();
    }
    MK.Update_l(l_i);
    return MK.Calculate_stable();
}

double IntegralKernel::Calculate_monopole(int l_i){
    double Sx = xk_dk_nPl(0,x);
    double Sp= xk_dk_nPl(0,xp);
    double dist = etaDistribution(eta);
    MK.Update_l(l_i);
    double F = MK.Calculate_stable();
    return dist*F*(Sp-Sx);
}

double IntegralKernel::Calculate_dipole(int l_i){
    double Gx = xk_dk_nPl(1,x);
    double Gp = xk_dk_nPl(1,xp);
    double dist = etaDistribution(eta);
    MK.Update_l(l_i);
    double F = MK.Calculate_stable();
    return dist*F*(Gp-Gx);
}

double IntegralKernel::Calculate_quadrupole(int l_i){
    double Qx = xk_dk_nPl(2,x);
    double Qp = xk_dk_nPl(2,xp);
    double dist = etaDistribution(eta);
    MK.Update_l(l_i);
    double F = MK.Calculate_stable();

    return dist*F*(Qp-Qx);
}

double IntegralKernel::Calculate_monopole_correction(int l_i){
    double Gx = xk_dk_nPl(1,x);
    double Gp = xk_dk_nPl(1,xp);
    double Qx = xk_dk_nPl(2,x);
    double Qp = xk_dk_nPl(2,xp);
    double dist = etaDistribution(eta);
    MK.Update_l(l_i);
    double F = MK.Calculate_stable();
    
    return dist*F*((Qp-Qx)+3*(Gp-Gx));
}

double IntegralKernel::sig_Boltzmann_Compton(double int_eta){
    eta = int_eta;
    Calculate_shared_variables();

    double r = 0.0; 

    if(run_mode=="monopole"){
        r = Calculate_monopole(l);
    } 
    else if(run_mode=="dipole"){ 
        r = Calculate_dipole(l);
    }
    else if(run_mode=="quadrupole"){
        r = Calculate_quadrupole(l);
    }
    else if(run_mode=="monopole_corr"){ 
        r = Calculate_monopole_correction(l);
    }
    //TODO: Some kind of calculate Full
    else if(run_mode=="kernel"){
        double dist = etaDistribution(eta);
        MK.Update_l(l);
        double F = MK.Calculate_stable();
        r = dist*F;
    }
    else {
        r = Calculate_monopole(l);
    }
    return r;
}

double IntegralKernel::eta_Int(double int_s){
    s = int_s;
    if (s==0){ s = 1e-10; }
    if (fixed_eta){ return sig_Boltzmann_Compton(eta); }

    double a=sinh(fabs(s)/2.0), b = 30.0;//lim=30.0, b=lim*(1.0+0.5*lim*0.05);
    double epsrel=Int_eps, epsabs=1.0e-300;
    
    return Integrate_using_Patterson_adaptive(a, b, epsrel, epsabs, [this](double int_var) { return this->sig_Boltzmann_Compton(int_var);});
}

double IntegralKernel::s_Int(){
    double a=-2.0, b=2.0;
    double epsrel=Int_eps, epsabs=1.0e-300;
    double integral = Integrate_using_Patterson_adaptive(a, b, epsrel, epsabs, [this](double int_var) { return this->eta_Int(int_var);});
    return integral;
}

double IntegralKernel::compute_kernel(int l_i, double s_i, electronDistribution eDistribution){
    run_mode = "kernel";
    l = l_i;
    etaDistribution = eDistribution;
    return eta_Int(s_i);
}

double IntegralKernel::compute_electron_kernel(int l_i, double s_i, electronDistribution eDistribution){
    photon_anisotropy = false;
    double return_val = compute_kernel(l_i,s_i,eDistribution);
    photon_anisotropy = true;
    return return_val;
}

double IntegralKernel::compute_distortion(string mode, electronDistribution eDistribution, int l_i){
    run_mode="monopole";
    if(mode=="monopole" || mode=="dipole" || mode=="quadrupole" || mode=="monopole_corr"){
        run_mode = mode;
    }
    etaDistribution = eDistribution;
    l=l_i;
    return s_Int();
}

double IntegralKernel::compute_electron_distortion(string mode, electronDistribution eDistribution, int l_i){
    photon_anisotropy = false;
    double return_val = compute_distortion(mode,eDistribution,l_i);
    photon_anisotropy = true;
    return return_val;
}

double IntegralKernel::compute_distortion_fixed_eta(string mode, double eta_i, int l){
    run_mode="monopole";
    fixed_eta = true;
    eta = eta_i;
    double result = compute_distortion(mode, etaDistribution);
    fixed_eta = false;
    return result;
}

double IntegralKernel::compute_beam_distortion(double mup_i, string mode, electronDistribution eDistribution){
    beam_kernel = true;
    mup = mup_i;
    double result = compute_distortion(mode, eDistribution);
    beam_kernel = false;
    return result;
}

double IntegralKernel::compute_beam_kernel(double mup_i, int l_i, double s_i, electronDistribution eDistribution){
    beam_kernel = true;
    mup = mup_i;
    double result = compute_kernel(l_i, s_i, eDistribution);
    beam_kernel = false;
    return result;
}

double IntegralKernel::compute_beam_distortion_fixed_eta(double mup_i, string mode, double eta_i){
    fixed_eta = beam_kernel = true;
    eta = eta_i;
    mup = mup_i;
    double result = compute_distortion(mode, etaDistribution);
    fixed_eta = beam_kernel = false;
    return result;
}


//==================================================================================================
//
// 3D integration carried out using Patterson scheme
//
// eps_Int : relative accuracy for numerical integration (lower than 10^-6 is hard to achieve)
//
// mode == "monopole"      --> only scattering of monopole without second order kinematic corr
// mode == "dipole"        --> only scattering of dipole     (first order kinematic correction)
// mode == "quadrupole"    --> only scattering of quadrupole (second order kinematic correction)
// mode == "monopole_corr" --> only scattering of second order kinematic correction to monopole
// mode == "all"           --> all terms added
// mode == "kin"           --> only kinematic terms
//
//==================================================================================================

double compute_SZ_distortion_kernel(double x, double betac, double muc, double eps_Int, 
                                    std::function<double(double)> eDistribution, string mode, int l){
    IntegralKernel szDistortion = IntegralKernel(x, betac, muc, eps_Int);
    return szDistortion.compute_distortion(mode, eDistribution,l);
}

void compute_SZ_distortion_kernel(vector<double> &Dn, vector<double> x, double betac, double muc, double eps_Int, 
                                  bool DI, std::function<double(double)> eDistribution, string mode, int l){
    int gridpoints = x.size();
    Dn.resize(gridpoints);
    Parameters fp = Parameters(); //This is just to get a value for the Dn_DI conversion 
    IntegralKernel szDistortion = IntegralKernel(x[0], betac, muc, eps_Int);
    Dn[0] = (DI ? pow(x[0],3.0)*fp.rare.Dn_DI_conversion() : 1.0)*szDistortion.compute_distortion(mode, eDistribution,l);
    for(int k = 1; k < gridpoints; k++){
        szDistortion.Update_x(x[k]);
        Dn[k] = szDistortion.compute_distortion(mode, eDistribution,l);
        if (DI) { Dn[k] *= pow(x[k],3.0)*fp.rare.Dn_DI_conversion(); }
    }
}

double compute_SZ_distortion_kernel(double x, Parameters fp, std::function<double(double)> eDistribution, int l){
    IntegralKernel szDistortion = IntegralKernel(x, fp);
    return szDistortion.compute_distortion(fp.rare.RunMode, eDistribution,l);
}

void compute_SZ_distortion_kernel(vector<double> &Dn, Parameters &fp, bool DI, 
                                        std::function<double(double)> eDistribution, int l){
    Dn.resize(fp.gridpoints);
    IntegralKernel szDistortion = IntegralKernel(fp.xcmb[0], fp);
    Dn[0] = (DI ? pow(fp.xcmb[0],3.0)*fp.rare.Dn_DI_conversion() : 1.0)*fp.Dtau*szDistortion.compute_distortion(fp.rare.RunMode, eDistribution,l);
    for(int k = 1; k < fp.gridpoints; k++){
        szDistortion.Update_x(fp.xcmb[k]);
        Dn[k] = fp.Dtau*szDistortion.compute_distortion(fp.rare.RunMode, eDistribution,l);
        if (DI) { Dn[k] *= pow(fp.xcmb[k],3.0)*fp.rare.Dn_DI_conversion(); }
    }
}

void compute_SZ_distortion_electron_kernel(vector<double> &Dn, Parameters &fp, bool DI, std::function<double(double)> eDistribution, int l){
    Dn.resize(fp.gridpoints);
    IntegralKernel szDistortion = IntegralKernel(fp.xcmb[0], fp);
    Dn[0] = (DI ? pow(fp.xcmb[0],3.0)*fp.rare.Dn_DI_conversion() : 1.0)*fp.Dtau*szDistortion.compute_electron_distortion(fp.rare.RunMode, eDistribution,l);
    for(int k = 1; k < fp.gridpoints; k++){
        szDistortion.Update_x(fp.xcmb[k]);
        Dn[k] = fp.Dtau*szDistortion.compute_electron_distortion(fp.rare.RunMode, eDistribution,l);
        if (DI) { Dn[k] *= pow(fp.xcmb[k],3.0)*fp.rare.Dn_DI_conversion(); }
    }
}

double compute_averaged_kernel(int l, double s, double eps_Int, std::function<double(double)> eDistribution){
    IntegralKernel szDistortion = IntegralKernel(1.0, 0.0, 1.0, eps_Int);
    return szDistortion.compute_kernel(l, s, eDistribution);
}

void compute_averaged_kernel(vector<double> &Dn, int l, vector<double> s, double eps_Int, std::function<double(double)> eDistribution){
    int gridpoints = s.size();
    Dn.resize(gridpoints);
    IntegralKernel szDistortion = IntegralKernel(1.0, 0.0, 1.0, eps_Int);
    for(int k = 0; k < gridpoints; k++){
        Dn[k] = szDistortion.compute_kernel(l, s[k], eDistribution);
    }
}

double compute_averaged_kernel(double s, Parameters fp, std::function<double(double)> eDistribution){
    IntegralKernel szDistortion = IntegralKernel(0.1, fp);
    return szDistortion.compute_kernel(fp.kernel.l, s, eDistribution);
}

void compute_averaged_kernel(vector<double> &Dn, Parameters &fp, std::function<double(double)> eDistribution){
    Dn.resize(fp.gridpoints);
    IntegralKernel szDistortion = IntegralKernel(0.1, fp);
    for(int k = 0; k < fp.gridpoints; k++){
        Dn[k] = szDistortion.compute_kernel(fp.kernel.l, fp.kernel.srange[k], eDistribution);
    }
}

void compute_averaged_electron_kernel(vector<double> &Dn, Parameters &fp, std::function<double(double)> eDistribution){
    Dn.resize(fp.gridpoints);
    IntegralKernel szDistortion = IntegralKernel(0.1, fp);
    for(int k = 0; k < fp.gridpoints; k++){
        Dn[k] = szDistortion.compute_electron_kernel(fp.kernel.l, fp.kernel.srange[k], eDistribution);
    }
}


double compute_SZ_distortion_fixed_eta(double x, double eta, double betac, double muc, double eps_Int, string mode){
    IntegralKernel szDistortion = IntegralKernel(x, betac, muc, eps_Int);
    return szDistortion.compute_distortion_fixed_eta(mode, eta);
}

void compute_SZ_distortion_fixed_eta(vector<double> &Dn, vector<double> x, double eta, double betac, double muc,
                                     double eps_Int, bool DI, string mode){
    int gridpoints = x.size();
    Dn.resize(gridpoints);
    Parameters fp = Parameters(); //This is just to get a value for the Dn_DI conversion 
    IntegralKernel szDistortion = IntegralKernel(x[0], betac, muc, eps_Int);
    Dn[0] = (DI ? pow(x[0],3.0)*fp.rare.Dn_DI_conversion() : 1.0)*szDistortion.compute_distortion_fixed_eta(mode, eta);
    for(int k = 1; k < gridpoints; k++){
        szDistortion.Update_x(x[k]);
        Dn[k] = szDistortion.compute_distortion_fixed_eta(mode, eta);
        if (DI) { Dn[k] *= pow(x[k],3.0)*fp.rare.Dn_DI_conversion(); }
    }
}

double compute_SZ_distortion_fixed_eta(double x, Parameters fp, double eta){
    IntegralKernel szDistortion = IntegralKernel(x, fp);
    return szDistortion.compute_distortion_fixed_eta(fp.rare.RunMode, eta);
}

void compute_SZ_distortion_fixed_eta(vector<double> &Dn, Parameters &fp, bool DI, double eta){
    Dn.resize(fp.gridpoints);
    IntegralKernel szDistortion = IntegralKernel(fp.xcmb[0], fp);
    Dn[0] = (DI ? pow(fp.xcmb[0],3.0)*fp.rare.Dn_DI_conversion() : 1.0)*fp.Dtau*szDistortion.compute_distortion_fixed_eta(fp.rare.RunMode, eta);
    for(int k = 1; k < fp.gridpoints; k++){
        szDistortion.Update_x(fp.xcmb[k]);
        Dn[k] = fp.Dtau*szDistortion.compute_distortion_fixed_eta(fp.rare.RunMode, eta);
        if (DI) { Dn[k] *= pow(fp.xcmb[k],3.0)*fp.rare.Dn_DI_conversion(); }
    }
}


double compute_SZ_distortion_beam_kernel(double x, double mup, double betac, double muc, double eps_Int,
                                         std::function<double(double)> eDistribution, string mode){
    IntegralKernel szDistortion = IntegralKernel(x, betac, muc, eps_Int);
    return szDistortion.compute_beam_distortion(mup, mode, eDistribution);
}

void compute_SZ_distortion_beam_kernel(vector<double> &Dn, vector<double> x, double mup, double betac, double muc, double eps_Int, 
                                       bool DI, std::function<double(double)> eDistribution, string mode){
    int gridpoints = x.size();
    Dn.resize(gridpoints);
    Parameters fp = Parameters(); //This is just to get a value for the Dn_DI conversion 
    IntegralKernel szDistortion = IntegralKernel(x[0], betac, muc, eps_Int);
    Dn[0] = (DI ? pow(x[0],3.0)*fp.rare.Dn_DI_conversion() : 1.0)*szDistortion.compute_beam_distortion(mup, mode, eDistribution);
    for(int k = 1; k < gridpoints; k++){
        szDistortion.Update_x(x[k]);
        Dn[k] = szDistortion.compute_beam_distortion(mup, mode, eDistribution);
        if (DI) { Dn[k] *= pow(x[k],3.0)*fp.rare.Dn_DI_conversion(); }
    }
}

double compute_SZ_distortion_beam_kernel(double x, Parameters fp, double mup, std::function<double(double)> eDistribution){
    IntegralKernel szDistortion = IntegralKernel(x, fp);
    return szDistortion.compute_beam_distortion(mup, fp.rare.RunMode, eDistribution);
}

void compute_SZ_distortion_beam_kernel(vector<double> &Dn, Parameters &fp, bool DI, double mup, std::function<double(double)> eDistribution){
    Dn.resize(fp.gridpoints);
    IntegralKernel szDistortion = IntegralKernel(fp.xcmb[0], fp);
    Dn[0] = (DI ? pow(fp.xcmb[0],3.0)*fp.rare.Dn_DI_conversion() : 1.0)*fp.Dtau*szDistortion.compute_beam_distortion(mup, fp.rare.RunMode, eDistribution);
    for(int k = 1; k < fp.gridpoints; k++){
        szDistortion.Update_x(fp.xcmb[k]);
        Dn[k] = fp.Dtau*szDistortion.compute_beam_distortion(mup, fp.rare.RunMode, eDistribution);
        if (DI) { Dn[k] *= pow(fp.xcmb[k],3.0)*fp.rare.Dn_DI_conversion(); }
    }
}


double compute_averaged_beam_kernel(int l, double s, double mup, double eps_Int, std::function<double(double)> eDistribution){
    IntegralKernel szDistortion = IntegralKernel(1.0, 0.0, 1.0, eps_Int);
    return szDistortion.compute_beam_kernel(mup, l, s, eDistribution);
}

void compute_averaged_beam_kernel(vector<double> &Dn, int l, vector<double> s, double mup, double eps_Int,
                                  std::function<double(double)> eDistribution){
    int gridpoints = s.size();
    Dn.resize(gridpoints);
    IntegralKernel szDistortion = IntegralKernel(1.0, 0.0, 1.0, eps_Int);
    for(int k = 0; k < gridpoints; k++){
        Dn[k] = szDistortion.compute_beam_kernel(mup, l, s[k], eDistribution);
    }
}

double compute_averaged_beam_kernel(double s, Parameters fp, double mup, std::function<double(double)> eDistribution){
    IntegralKernel szDistortion = IntegralKernel(0.1, fp);
    return szDistortion.compute_beam_kernel(mup, fp.kernel.l, s, eDistribution);
}

void compute_averaged_beam_kernel(vector<double> &Dn, Parameters &fp, double mup, std::function<double(double)> eDistribution){
    Dn.resize(fp.gridpoints);
    IntegralKernel szDistortion = IntegralKernel(0,fp);
    for(int k = 0; k < fp.gridpoints; k++){
        Dn[k] = szDistortion.compute_beam_kernel(mup, fp.kernel.l, fp.kernel.srange[k], eDistribution);
    }
}


double compute_SZ_distortion_beam_kernel_fixed_eta(double x, double mup, double eta, double betac, double muc, 
                                                   double eps_Int, string mode){
    IntegralKernel szDistortion = IntegralKernel(x, betac, muc, eps_Int);
    return szDistortion.compute_beam_distortion_fixed_eta(mup, mode, eta);
}

void compute_SZ_distortion_beam_kernel_fixed_eta(vector<double> &Dn, vector<double> x, double mup, double eta, double betac, double muc, 
                                       double eps_Int, bool DI, string mode){
    int gridpoints = x.size();
    Dn.resize(gridpoints);
    Parameters fp = Parameters(); //This is just to get a value for the Dn_DI conversion 
    IntegralKernel szDistortion = IntegralKernel(x[0], betac, muc, eps_Int);
    Dn[0] = (DI ? pow(x[0],3.0)*fp.rare.Dn_DI_conversion() : 1.0)*szDistortion.compute_beam_distortion_fixed_eta(mup, mode, eta);
    for(int k = 1; k < gridpoints; k++){
        szDistortion.Update_x(x[k]);
        Dn[k] = szDistortion.compute_beam_distortion_fixed_eta(mup, mode, eta);
        if (DI) { Dn[k] *= pow(x[k],3.0)*fp.rare.Dn_DI_conversion(); }
    }
}

double compute_SZ_distortion_beam_kernel_fixed_eta(double x, Parameters fp, double mup, double eta){
    IntegralKernel szDistortion = IntegralKernel(x, fp);
    return szDistortion.compute_beam_distortion_fixed_eta(mup, fp.rare.RunMode, eta);
}

void compute_SZ_distortion_beam_kernel_fixed_eta(vector<double> &Dn, Parameters &fp, bool DI, double mup, double eta){
    Dn.resize(fp.gridpoints);
    IntegralKernel szDistortion = IntegralKernel(fp.xcmb[0], fp);
    Dn[0] = (DI ? pow(fp.xcmb[0],3.0)*fp.rare.Dn_DI_conversion() : 1.0)*fp.Dtau*szDistortion.compute_beam_distortion_fixed_eta(mup, fp.rare.RunMode, eta);
    for(int k = 1; k < fp.gridpoints; k++){
        szDistortion.Update_x(fp.xcmb[k]);
        Dn[k] = fp.Dtau*szDistortion.compute_beam_distortion_fixed_eta(mup, fp.rare.RunMode, eta);
        if (DI) { Dn[k] *= pow(fp.xcmb[k],3.0)*fp.rare.Dn_DI_conversion(); }
    }
}