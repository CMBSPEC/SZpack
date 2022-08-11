//==================================================================================================
//
// This program allows computing the thermal SZ and kinematic effect. Up to 10th order temperature
// corrections are allowed. The basis functions are computed using recursion relations, based on 
// Eulerian numbers to determine the derivatives x^k d^k nPl(x)/ dx^k.
//
// The expressions for the thermal SZ effect are equivalent to those of Itoh et al., 1998 and 
// Shimon & Rephaeli, 2004, but to higher order kinematic corrections are computed according to 
// Chluba, Nagai, Sazonov & Nelson, 2012. 
//
//==================================================================================================
//
// Author: Jens Chluba  (CITA, University of Toronto)
//
// first implementation: April 2012
// last modification   : Sept  2013
//
//==================================================================================================
// 12th Sept 2013: fixed bug for negative beta_para (thanks to John ZuHone & Tony Mroczkowski)
//  4st  Aug 2012: added derivatives of basis functions in the CMB rest frame
// 10th July 2012: added basis functions in the CMB rest frame

//==================================================================================================
// Standards
//==================================================================================================
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include <iomanip>
#include <cmath>
#include <limits.h>
#include <vector>

//==================================================================================================
// required libs
//==================================================================================================
#include "physical_consts.h"
#include "routines.h"
#include "nPl_derivatives.h"
#include "SZ_asymptotic.h"

//==================================================================================================
// namespace
//==================================================================================================
using namespace std;

//==================================================================================================
//
// analytic tables of moments for monopole, dipole and quadrupole scattering up to The^11
//
//==================================================================================================
const int order_Tmax_Y=11;

//==================================================================================================
// monopole scattering
//==================================================================================================
const double a00[11][23]={
//Y0
    {0, 4, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
//Y1
    {0, 10, 23.5, 8.4, 0.7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
//Y2
    {0, 7.5, 127.875, 173.6, 65.8, 8.8, 0.36666666666666664, 
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
//Y3
    {0, -7.5, 313.125, 1419.6, 1425.3, 531.2571428571429, 86.13571428571429, 6.095238095238095, 
     0.1523809523809524, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
//Y4
    {0, 4.21875, 237.3046875, 6239.1, 15368.175, 12438.9, 4446.2875, 788.952380952381, 
     71.58095238095238, 3.142857142857143, 0.05238095238095238, 
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
//Y5
    {0, 5.625, -234.84375, 14458.5, 99428.625, 165713.7857142857, 113100.74107142857, 
     38405.90476190476, 7103.333333333333, 740.936507936508, 43.067460317460316, 
     1.288888888888889, 0.015343915343915344, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
//Y6
    {0, -29.00390625, 125.6396484375, 11271.09375, 397293.3984375, 1.42874128125e6, 
     1.78985588671875e6, 1.0596171904761905e6, 340706.17619047617, 63964.445238095235, 
     7255.1140873015875, 499.1111111111111, 20.193650793650793, 0.43851851851851853, 
     0.003915343915343915, 0, 0, 0, 0, 0, 0, 0, 0},
//Y7
    {0, 84.375, 198.28125, -10843.875, 893108.53125, 8.269474660714285e6, 1.937579644419643e7, 
     1.9396852904761903e7, 1.0158613633333333e7, 3.0913358491341993e6, 579982.4908189033, 
     69245.74141414142, 5316.226118326118, 259.5399326599327, 7.733771043771044, 
     0.12744588744588745, 0.000885040885040885, 0, 0, 0, 0, 0, 0},
//Y8
    {0, -232.6080322265625, -971.7750549316406, 4978.23046875, 706801.3857421875, 
     3.1323163655691963e7, 1.4866270463148716e8, 2.5343754924107143e8, 2.1099909958125e8, 
     9.919729982728626e7, 2.8592066184899215e7, 5.309129966666667e6, 653653.1473665223, 
     54059.79919191919, 2997.7936616161614, 109.29316017316017, 2.499047619047619, 
     0.03232323232323232, 0.00017957351290684624, 0, 0, 0, 0},
//Y9
    {0, 711.9140625, 2835.791015625, 11977.875, -668374.03125, 6.906976891071428e7, 
     8.058455207879465e8, 2.440941352321429e9, 3.235491342175e9, 2.303830395858604e9, 
     9.833430327076434e8, 2.6885084103193474e8, 4.9103883772602394e7, 6.149341245345766e6, 
     535390.5228574204, 32501.263856143854, 1363.3057076257076, 38.54048174048174, 
     0.6976845376845376, 0.0072717406050739385, 0.00003305336638669972, 0, 0},
//Y10
    {0, -2582.6351165771484, -7972.811794281006, -52974.56359863281, 274433.91036987305, 
     5.515719944536482e7, 2.9496299591826196e9, 1.7426481072717632e10, 3.769307828481641e10, 
     4.060211026653589e10, 2.524367246341206e10, 9.871467951450539e9, 2.564478584196426e9, 
     4.5886954101723415e8, 5.7897765075480334e7, 5.2231502993006995e6, 338722.2291486291, 
     15739.505516705516, 516.9642912642913, 11.665205341395817, 0.1713468189658666, 
     0.001469852898424327, 5.5676246152436625e-6},
};

//==================================================================================================
// dipole scattering
//==================================================================================================
const double d10[11][23]={
//D0
    {-0.4, -1.6, -0.4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
//D1
    {-0.2, -4.8, -13.2, -4.8, -0.4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
//D2
    {2.907142857142857, -6.6571428571428575, -74.52142857142857, -110.74285714285715, 
     -43.22857142857143, -5.828571428571428, -0.24285714285714285, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
//D3
    {-12.964285714285714, 16.714285714285715, -195.10714285714286, -930.1142857142858, 
     -983.8428571428572, -376.1714285714286, -61.673809523809524, -4.380952380952381, 
     -0.10952380952380952, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
//D4
    {52.14955357142857, -54.25892857142857, -115.88616071428571, -4157.661904761905, 
     -10813.77738095238, -9066.62380952381, -3307.109325396825, -593.2698412698413, 
     -54.11746031746032, -2.380952380952381, -0.03968253968253968, 0, 0, 0, 0, 0, 
     0, 0, 0, 0, 0, 0, 0},
//D5
    {-222.27678571428572, 219.46428571428572, 30.330357142857142, -9612.547619047618, 
     -70650.18452380953, -122506.70952380952, -85787.77956349206, -29604.926984126985, 
     -5529.61746031746, -580.0857142857143, -33.81920634920635, -1.0133333333333334, 
     -0.012063492063492064, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
//D6
    {1049.978312702922, -1035.476359577922, 446.704976663961, -7682.669237012987, 
     -283443.83938717534, -1.0643261265422078e6, -1.3725463310301676e6, -828629.9151515152, 
     -269942.39145021647, -51116.84062049062, -5829.888084415585, -402.409696969697, 
     -16.311130351130352, -0.3544781144781145, -0.003164983164983165, 0, 0, 0, 0, 0, 0, 0, 0},
//D7
    {-5514.879058441558, 5472.691558441558, -2867.970779220779, 8154.86525974026, 
     -638448.7252435065, -6.18568069512987e6, -1.4954438377448592e7, -1.5302651941125542e7, 
     -8.138805677813853e6, -2.5035448205627706e6, -473232.2047871573, -56790.098585858585, 
     -4374.778326118326, -214.03420875420875, -6.385478595478595, -0.10528138528138528, 
     -0.0007311207311207312, 0, 0, 0, 0, 0, 0},
//D8
    {31915.806382643237, -31799.502366529956, 16517.707421057945, -8674.850044486764, 
     -503392.5539509319, -2.3479717113289054e7, -1.1518388765106362e8, -2.0106314515351316e8, 
     -1.7028461288606653e8, -8.105772537922598e7, -2.357550170187599e7, -4.406290473426573e6, 
     -545036.8997197247, -45225.05177933178, -2513.4934068709067, -91.76871128871129, 
     -2.1000865800865802, -0.027172827172827173, -0.00015096015096015096, 0, 0, 0, 0},
//D9
    {-201512.5000975587, 201156.5430663087, -102386.53325190434, 25706.426448551447, 
     469310.2594592907, -5.181291480457043e7, -6.257891822237918e8, -1.94348654784016e9, 
     -2.6240012870239573e9, -1.8942035577713912e9, -8.168251949195232e8, -2.2504155431751582e8, 
     -4.133810694754412e7, -5.19882679994228e6, -454043.6548696145, -27624.674966937822, 
     -1160.5508548594262, -32.842040499183355, -0.5948828596447644, -0.006201910963815726, 
     -0.00002819050438098057, 0, 0},
//D10
    {1.3770014727438318e6, -1.3757101551855432e6, 693026.1599368268, -194466.477051391, 
     -139990.75565867848, -4.1374198674687006e7, -2.2933736574899874e9, -1.3907025699113136e10, 
     -3.0671405859723873e10, -3.3528181988669834e10, -2.1080598242516373e10, -8.314546476846714e9, 
     -2.1742668777989497e9, -3.9101316848326534e8, -4.9526029225350976e7, -4.480926443446077e6, 
     -291220.874163297, -13553.74599559171, -445.67204592585546, -10.064008372579801, 
     -0.14789282674996962, -0.0012689109831966975, -4.806480996957188e-6}
};

//==================================================================================================
// quadrupole scattering
//==================================================================================================
const double q20[11][23]={
//Q0
    {-0.6, 0.4, 0.1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
//Q1
    {2.6142857142857143, -0.7142857142857143, 4.107142857142857, 1.7142857142857142, 
     0.14285714285714285, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
//Q2
    {-10.725, 10.35, 20.5875, 44, 19, 2.6285714285714286, 0.10952380952380952, 
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
//Q3
    {45.69642857142857, -45.32142857142857, 85.09821428571429, 380.1714285714286, 
     462.01428571428573, 188.97142857142856, 31.873809523809523, 2.2857142857142856, 
     0.05714285714285714, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
//Q4
    {-210.6328125, 210.421875, -57.64453125, 1804.1, 5240.425, 4772.928571428572, 
     1826.072023809524, 336, 31.02857142857143, 1.3714285714285714, 0.022857142857142857, 
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
//Q5
    {1063.182224025974, -1063.463474025974, 484.903612012987, 3967.1623376623374, 
     34982.81655844156, 66142.35454545454, 49045.89507575757, 17544.567965367965, 
     3348.8199134199135, 355.7183261183261, 20.873823953823955, 0.6270707070707071, 
     0.007465127465127465, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
//Q6
    {-5871.436751724837, 5872.886947037337, -2912.0406414874187, 4192.01663961039, 
     141206.77932224027, 583342.9385551948, 800866.4830458604, 504409.8787878788, 
     169022.6183982684, 32602.20432900433, 3762.170072150072, 261.5353535353535, 
     10.642135642135642, 0.23164983164983166, 0.002068302068302068, 
     0, 0, 0, 0, 0, 0, 0, 0},
//Q7
    {35292.037556193805, -35296.256306193805, 17689.893778096903, -8992.490634365635, 
     320867.1673014485, 3.4179723781468533e6, 8.835846921350524e6, 9.473327058408258e6, 
     5.205919926789877e6, 1.6383247194139194e6, 314605.851990232, 38161.21442113442, 
     2960.6207104007103, 145.49163429163428, 4.351446701446702, 0.07182151182151182, 
     0.0004987604987604987, 0, 0, 0, 0, 0, 0},
//Q8
    {-229542.7946249612, 229554.42502657254, -114977.38272507826, 39737.5062320492, 
     242648.61200225557, 1.3031721730250219e7, 6.858038253698304e7, 1.2583692569567932e8, 
     1.1047646109088099e8, 5.3976345275666e7, 1.5997834142372211e7, 3.0309264174159174e6, 
     378563.59738317237, 31625.883241203243, 1765.8162295244438, 64.66213405641977, 
     1.482292945150088, 0.01919350490779062, 0.000106630582821059, 0, 0, 0, 0},
//Q9
    {1.6067541794240915e6, -1.6067897751272165e6, 803979.8436182957, -264539.7860889111, 
     -172164.8838349151, 2.8786166000874124e7, 3.7429834605714595e8, 1.2250477237887113e9, 
     1.7188799543150475e9, 1.2767353532436564e9, 5.623919607973943e8, 1.5742119784871796e8, 
     2.9260416806568433e7, 3.712189417973138e6, 326285.12825202575, 19942.793846153847, 
     840.4897635697636, 23.834254634254634, 0.43224701224701223, 0.004508789270694032, 
     0.000020494496684972875, 0, 0},
//Q10
    {-1.2045816323789336e7, 1.2045945455545165e7, -6.024631856009353e6, 1.992983724701202e6, 
     -401994.0516848569, 2.3073725771758392e7, 1.3750719824278262e9, 8.806513463938858e9, 
     2.0224989460246635e10, 2.2793558864644966e10, 1.466735770630554e10, 5.88834861346897e9, 
     1.56079941696492e9, 2.836085274380134e8, 3.62061731549366e7, 3.295347646547738e6, 
     215120.70086976516, 10044.40915877773, 331.0326798333941, 7.4865996085043705, 
     0.1101162823067585, 0.0009451712308855166, 3.5801940563845327e-6}
};

//==================================================================================================
//
// spectral functions defined in Chluba, Nagai, Sazonov, Nelson, 2012
//
//==================================================================================================

//==================================================================================================
// th-SZ
//==================================================================================================
void compute_Y(double x, vector<double> &Yk)
{
    int ymax=Yk.size()-1;
    vector<double> Pk(2*(ymax+1)+1), beta(2*(ymax+1)+2);
    
    for(int k=0; k<=2*(ymax+1); k++) Pk[k]=Pfunc(k, x);
    
    double exp_mx=exp(-x);
    double fac=(-x)/(1.0-exp_mx);
    double nPl=exp_mx/(1.0-exp_mx);
    
    for(int k=0; k<=ymax; k++)
    {
        int kmax=2*(k+1);
        beta[kmax+1]=0.0;
        
        for(int l=kmax; l>=1; l--) beta[l]=fac*(a00[k][l]*Pk[l]+beta[l+1]);
        
        Yk[k]=nPl*beta[1];
    }
    
    return;
}

//==================================================================================================
// kinematic correction to monopole
//==================================================================================================
void compute_Ykin(double x, vector<double> &Ykin)
{
    int ymax=Ykin.size()-1;
    vector<double> Pk(2*(ymax+2)+1), beta(2*(ymax+2)+2);
    
    for(int k=0; k<=2*(ymax+2); k++) Pk[k]=Pfunc(k, x);
    
    double exp_mx=exp(-x);
    double fac=(-x)/(1.0-exp_mx);
    double nPl=exp_mx/(1.0-exp_mx);
    
    for(int k=0; k<=ymax; k++)
    {
        int kmax=2*(k+1);
        beta[kmax+1]=0.0;
        
        for(int l=kmax; l>=1; l--)
        {
            beta[l]=a00[k][l]*( Pk[l]*l*(l+2.0)+fac*(Pk[l+1]*(2.0*l+3.0)+fac*Pk[l+2]) );
            beta[l]=fac*(beta[l]+beta[l+1]);
        }
        
        Ykin[k]=nPl*beta[1]/6.0;
    }
    
    return;
}

//==================================================================================================
// kinematic correction to monopole in CMB frame (added 10.07 by JC)
//==================================================================================================
void compute_M_CMB(double x, vector<double> &Ykin)
{
    int ymax=Ykin.size()-1;
    vector<double> Pk(2*(ymax+2)+1), beta(2*(ymax+2)+2);
    
    for(int k=0; k<=2*(ymax+2); k++) Pk[k]=Pfunc(k, x);
    
    double exp_mx=exp(-x);
    double fac=(-x)/(1.0-exp_mx);
    double nPl=exp_mx/(1.0-exp_mx);
    
    for(int k=0; k<=ymax; k++)
    {
        int kmax=2*(k+1);
        beta[kmax+1]=0.0;
        
        for(int l=kmax; l>=1; l--)
        {
            beta[l]=(a00[k][l]-d10[k][l])*( Pk[l]*l*(l+2.0)+fac*(Pk[l+1]*(2.0*l+3.0)+fac*Pk[l+2]) );
            beta[l]=fac*(beta[l]+beta[l+1]);
        }
        
        Ykin[k]=nPl*( beta[1]-d10[k][0]*fac*(Pk[1]*3.0+fac*Pk[2]) )/3.0;
    }
    
    return;
}

//==================================================================================================
// kinematic correction to dipole
//==================================================================================================
void compute_Dkin(double x, vector<double> &Dkin)
{
    int ymax=Dkin.size()-1;
    vector<double> Pk(2*(ymax+2)+1), beta(2*(ymax+2)+2);
    
    for(int k=0; k<=2*(ymax+2); k++) Pk[k]=Pfunc(k, x);
    
    double exp_mx=exp(-x);
    double fac=(-x)/(1.0-exp_mx);
    double nPl=exp_mx/(1.0-exp_mx);
    
    for(int k=0; k<=ymax; k++)
    {
        int kmax=2*(k+1);
        beta[kmax+1]=0.0;
        
        for(int l=kmax; l>=1; l--)
        {
            beta[l]=d10[k][l]*( Pk[l]*l+fac*Pk[l+1] );
            beta[l]=fac*(beta[l]+beta[l+1]);
        }
        
        Dkin[k]=nPl*(beta[1]+d10[k][0]*fac*Pk[1]);
    }
    
    return;
}

//==================================================================================================
// kinematic correction to dipole in CMB frame (added 10.07 by JC)
//==================================================================================================
void compute_D_CMB(double x, vector<double> &Dkin)
{
    int ymax=Dkin.size()-1;
    vector<double> Pk(2*(ymax+2)+1), beta(2*(ymax+2)+2);
    
    for(int k=0; k<=2*(ymax+2); k++) Pk[k]=Pfunc(k, x);
    
    double exp_mx=exp(-x);
    double fac=(-x)/(1.0-exp_mx);
    double nPl=exp_mx/(1.0-exp_mx);
    
    for(int k=0; k<=ymax; k++)
    {
        int kmax=2*(k+1);
        beta[kmax+1]=0.0;
        
        for(int l=kmax; l>=1; l--)
        {
            beta[l]=(d10[k][l]-a00[k][l])*( Pk[l]*l+fac*Pk[l+1] );
            beta[l]=fac*(beta[l]+beta[l+1]);
        }
        
        Dkin[k]=nPl*(beta[1]+d10[k][0]*fac*Pk[1]);
    }
    
    return;
}

//==================================================================================================
// kinematic correction to quadrupole
//==================================================================================================
void compute_Qkin(double x, vector<double> &Qkin)
{
    int ymax=Qkin.size()-1;
    vector<double> Pk(2*(ymax+2)+1), beta(2*(ymax+2)+2);
    
    for(int k=0; k<=2*(ymax+2); k++) Pk[k]=Pfunc(k, x);
    
    double exp_mx=exp(-x);
    double fac=(-x)/(1.0-exp_mx);
    double nPl=exp_mx/(1.0-exp_mx);
    
    for(int k=0; k<=ymax; k++)
    {
        int kmax=2*(k+1);
        beta[kmax+1]=0.0;
        
        for(int l=kmax; l>=1; l--)
        {
            beta[l]=q20[k][l]*( Pk[l]*l*(l-1.0)+fac*(2.0*l*Pk[l+1]+fac*Pk[l+2]) );
            beta[l]=fac*(beta[l]+beta[l+1]);
        }
        
        Qkin[k]=nPl*(beta[1]+q20[k][0]*fac*fac*Pk[2])/3.0;
    }
    
    return;
}

//==================================================================================================
// kinematic correction to quadrupole in CMB frame (added 10.07 by JC)
//==================================================================================================
void compute_Q_CMB(double x, vector<double> &Qkin)
{
    int ymax=Qkin.size()-1;
    vector<double> Pk(2*(ymax+2)+1), beta(2*(ymax+2)+2);
    
    for(int k=0; k<=2*(ymax+2); k++) Pk[k]=Pfunc(k, x);
    
    double exp_mx=exp(-x);
    double fac=(-x)/(1.0-exp_mx);
    double nPl=exp_mx/(1.0-exp_mx);
    
    for(int k=0; k<=ymax; k++)
    {
        int kmax=2*(k+1);
        beta[kmax+1]=0.0;
        
        for(int l=kmax; l>=1; l--)
        {
            beta[l]=(q20[k][l]+a00[k][l]-2.0*d10[k][l])
                                 *( Pk[l]*l*(l-1.0)+fac*(2.0*l*Pk[l+1]+fac*Pk[l+2]) );
            beta[l]=fac*(beta[l]+beta[l+1]);
        }
        
        Qkin[k]=nPl*(beta[1]+(q20[k][0]-2.0*d10[k][0])*fac*fac*Pk[2])/3.0;
    }
    
    return;
}

//==================================================================================================
//
// approximation for relativistic SZ effect
//
//==================================================================================================
// const temperature case y^(k)= Theta^k tau_e
//==================================================================================================
double Dn_thSZ_appr(double x, double Theta, int max_order)
{
    vector<double> Yk(max_order);
    compute_Y(x, Yk);
    
    double r=Theta*Yk[max_order-1];
    for(int k=max_order-2; k>=0; k--) r=Theta*(r+Yk[k]);
    
    return r;
}

//==================================================================================================
double Dn_kSZ_appr1(double x, double Theta, double betac, double muc, int max_order)
{
    // normal k-SZ effect
    double rnorm=x*exp(-x)/(1.0-exp(-x))/(1.0-exp(-x));

    if(max_order==0) return betac*muc*rnorm;
    
    vector<double> Dkin(max_order);
    compute_Dkin(x, Dkin);
    
    double r=Theta*Dkin[max_order-1];
    for(int k=max_order-2; k>=0; k--) r=Theta*(r+Dkin[k]);
    
    // normal k-SZ effect
    r+=rnorm;
    
    return r*betac*muc;
}

//==================================================================================================
double Dn_kSZ_appr2(double x, double Theta, double betac, double muc, int max_order)
{
    // quadratic k-SZ effect
    double rnorm=-3.0/10.0*x*exp(-x)/(1.0-exp(-x))/(1.0-exp(-x))*x*(1.0+exp(-x))/(1.0-exp(-x));
    
    if(max_order==0) return betac*betac*(1.5*muc*muc-0.5)*rnorm;

    vector<double> Qkin(max_order);
    compute_Qkin(x, Qkin);
    
    double r=Theta*Qkin[max_order-1];
    for(int k=max_order-2; k>=0; k--) r=Theta*(r+Qkin[k]);
    
    // second order k-SZ effect
    r+=rnorm;

    return r*betac*betac*(1.5*muc*muc-0.5);
}

//==================================================================================================
double Dn_kSZ_appr2m(double x, double Theta, double betac, int max_order)
{
    if(max_order==0) return 0.0;
    
    vector<double> Ykin(max_order);
    compute_Ykin(x, Ykin);
    
    double r=Theta*Ykin[max_order-1];
    for(int k=max_order-2; k>=0; k--) r=Theta*(r+Ykin[k]);
    
    return r*betac*betac;
}

//==================================================================================================
//
// approximation for relativistic SZ effect for CMB rest frame terms (added 10.07 by JC)
//
//==================================================================================================
double Dn_kSZ_appr1_CMB(double x, double Theta, double betac, double muc, int max_order)
{
    // normal k-SZ effect
    double rnorm=x*exp(-x)/(1.0-exp(-x))/(1.0-exp(-x));

    if(max_order==0) return betac*muc*rnorm;
    
    vector<double> Dkin(max_order);
    compute_D_CMB(x, Dkin);
    
    double r=Theta*Dkin[max_order-1];
    for(int k=max_order-2; k>=0; k--) r=Theta*(r+Dkin[k]);
    
    // normal k-SZ effect
    r+=rnorm;
    
    return r*betac*muc;
}

//==================================================================================================
double Dn_kSZ_appr2_CMB(double x, double Theta, double betac, double muc, int max_order)
{
    // quadratic k-SZ effect
    double rnorm=11.0/30.0*x*exp(-x)/(1.0-exp(-x))/(1.0-exp(-x))*x*(1.0+exp(-x))/(1.0-exp(-x));

    if(max_order==0) return betac*betac*(1.5*muc*muc-0.5)*rnorm;
    
    vector<double> Qkin(max_order);
    compute_Q_CMB(x, Qkin);
    
    double r=Theta*Qkin[max_order-1];
    for(int k=max_order-2; k>=0; k--) r=Theta*(r+Qkin[k]);
    
    // second order k-SZ effect
    r+=rnorm;
    
    return r*betac*betac*(1.5*muc*muc-0.5);
}

//==================================================================================================
double Dn_kSZ_appr2m_CMB(double x, double Theta, double betac, int max_order)
{
    // quadratic k-SZ effect on monopole
    double rnorm=x*exp(-x)/(1.0-exp(-x))/(1.0-exp(-x))*(x*(1.0+exp(-x))/(1.0-exp(-x))-3.0)/3.0;

    if(max_order==0) return betac*betac*rnorm;
    
    vector<double> Ykin(max_order);
    compute_M_CMB(x, Ykin);
    
    double r=Theta*Ykin[max_order-1];
    for(int k=max_order-2; k>=0; k--) r=Theta*(r+Ykin[k]);

    // second order k-SZ effect on monopole
    r+=rnorm;

    return r*betac*betac;
}

//==================================================================================================
//
// compute Dn using asymptotic expansion in cluster frame
//
// mode == "monopole"      --> only scattering of monopole without second order kinematic corr
// mode == "dipole"        --> only scattering of dipole     (first order kinematic correction)
// mode == "quadrupole"    --> only scattering of quadrupole (second order kinematic correction)
// mode == "monopole_corr" --> only scattering of second order kinematic correction to monopole
// mode == "all"           --> all terms added
// mode == "kin"           --> only kinematic terms
//
//==================================================================================================
double compute_SZ_distortion_asymptotic(double x, 
                                        double The, double betac, double muc, 
                                        int Te_order, int betac_order, 
                                        string mode)
{
    double Dn_th=0.0, Dn_k1=0.0, Dn_k2=0.0, Dn_k2m=0.0; 
    
    Te_order++;
    
    if(betac_order>0)
    {
        if(mode=="dipole" || mode=="all" || mode=="kin") 
            Dn_k1=Dn_kSZ_appr1 (x, The, betac, muc, Te_order);
        
        if(betac_order>1)
        {
            if(mode=="monopole_corr" || mode=="all" || mode=="kin") 
                Dn_k2m=Dn_kSZ_appr2m(x, The, betac, Te_order);
           
            if(mode=="quadrupole" || mode=="all" || mode=="kin") 
                Dn_k2=Dn_kSZ_appr2 (x, The, betac, muc, Te_order);
        }
    }
    
    if(mode=="monopole" || mode=="all") 
        Dn_th=Dn_thSZ_appr (x, The, Te_order);
    
    return  Dn_th+Dn_k1+Dn_k2+Dn_k2m;
}

//==================================================================================================
//
// compute Dn using asymptotic expansion in CMB rest frame (added 10.07 by JC)
//
// mode == "monopole"      --> only monopole part without second order kinematic corr
// mode == "dipole"        --> only dipolar part     (first order kinematic correction)
// mode == "quadrupole"    --> only quadrupolar part (second order kinematic correction)
// mode == "monopole_corr" --> only second order kinematic correction to monopole part
// mode == "all"           --> all terms added
// mode == "kin"           --> only kinematic terms
//
//==================================================================================================
double compute_SZ_distortion_asymptotic_CMB(double x, 
                                            double The, double betac, double muc, 
                                            int Te_order, int betac_order,
                                            string mode)
{
    double Dn_th=0.0, Dn_k1=0.0, Dn_k2=0.0, Dn_k2m=0.0; 
    
    Te_order++;
    
    if(betac_order>0)
    {
        if(mode=="dipole" || mode=="all" || mode=="kin") 
            Dn_k1=Dn_kSZ_appr1_CMB(x, The, betac, muc, Te_order);
        
        if(betac_order>1)
        {
            if(mode=="monopole_corr" || mode=="all" || mode=="kin") 
                Dn_k2m=Dn_kSZ_appr2m_CMB(x, The, betac, Te_order);
            
            if(mode=="quadrupole" || mode=="all" || mode=="kin") 
                Dn_k2=Dn_kSZ_appr2_CMB(x, The, betac, muc, Te_order);
        }
    }
    
    if(mode=="monopole" || mode=="all") 
        Dn_th=Dn_thSZ_appr (x, The, Te_order);
    
    return  Dn_th+Dn_k1+Dn_k2+Dn_k2m;
}

//==================================================================================================
// analytic derivatives in betac_parallel and betac_perp
//==================================================================================================
double Dn_asymptotic_CMB_for_The(double x, double The, double beta_para, double beta_perp)
{ 
    double betac=sqrt(beta_para*beta_para+beta_perp*beta_perp);
    double muc= (betac==0 ? 0 : beta_para/betac);

    double Dn_th= Dn_thSZ_appr (x, The, 11);
    double Dn_k1= Dn_kSZ_appr1_CMB(x, The, betac, muc, 11);
    double Dn_k2 =Dn_kSZ_appr2_CMB(x, The, betac, muc, 11); 
    double Dn_k2m=Dn_kSZ_appr2m_CMB(x, The, betac, 11);
    
    return Dn_th+Dn_k1+Dn_k2+Dn_k2m;
}

double Dn_dbeta_para_asymptotic_CMB(double x, double The, double beta_para, double beta_perp)
{
    // 12th Sept 2013: fixed bug for negative beta_para (thanks to John ZuHone & Tony Mroczkowski)
    double sig=beta_para/fabs(beta_para);    
    double Dn_k1= Dn_kSZ_appr1_CMB(x, The, 1.0, 1.0, 11);
    double Dn_k2 =2.0*sig*Dn_kSZ_appr2_CMB(x, The, sqrt(fabs(beta_para)), 1.0, 11);
    double Dn_k2m=2.0*beta_para*Dn_kSZ_appr2m_CMB(x, The, 1.0, 11);

    return Dn_k1+Dn_k2+Dn_k2m;
}

double Dn_d2beta_para_asymptotic_CMB(double x, double The, double beta_para, double beta_perp)
{ 
    double Dn_k2 =2.0*Dn_kSZ_appr2_CMB(x, The, 1.0, 1.0, 11); 
    double Dn_k2m=2.0*Dn_kSZ_appr2m_CMB(x, The, 1.0, 11); 
    
    return  (Dn_k2+Dn_k2m) / 2.0;
}

double Dn_dbeta2_perp_asymptotic_CMB(double x, double The, double beta_para, double beta_perp)
{ 
    // 12th Sept 2013: betac^2 P_2(muc) == beta_para^2 - 0.5*beta_perp^2
    double Dn_k2 =-0.5*Dn_kSZ_appr2_CMB(x, The, 1.0, 1.0, 11);
    double Dn_k2m= Dn_kSZ_appr2m_CMB(x, The, 1.0, 11); 
    
    return  Dn_k2+Dn_k2m;
}

//==================================================================================================
void compute_all_Te_derivatives_upto_dThe(double x, double The, 
                                          double beta_para, double beta_perp, 
                                          int dThe, 
                                          double (*f)(double, double, double, double),
                                          vector<double> &dDn_dThe)
{
    double d0, dp, dm, dp2, dm2;
    double eps=0.01;
    
    dDn_dThe.resize(dThe+1);
    d0=f(x, The, beta_para, beta_perp);
    dDn_dThe[0]=d0;
    
    if(dThe>0)
    {
        dp =f(x, The*(1.0+eps),   beta_para, beta_perp);
        dm =f(x, The*(1.0-eps),   beta_para, beta_perp);
        dp2=f(x, The*(1.0+2*eps), beta_para, beta_perp);
        dm2=f(x, The*(1.0-2*eps), beta_para, beta_perp);
        
        dDn_dThe[1]=(-dp2+8.0*dp-8.0*dm+dm2)/12.0/eps;

        if(dThe>1)
        {
            dDn_dThe[2]=(dp-2.0*d0+dm)/pow(eps, 2) / (2.0);
        
            if(dThe>2)
            {
                dDn_dThe[3]=(dp2-2.0*dp+2.0*dm-dm2)/2.0/pow(eps, 3) / (6.0);
            
                if(dThe>3)
                    dDn_dThe[4]=(dp2-4.0*dp+6.0*d0-4.0*dm+dm2)/pow(eps, 4) / (24.0);
            }
        }
    }
    
    return;
}

//==================================================================================================
// Derivatives (The^k d^k_dThe /k!) (betapara^m d^m_dbetapara /m!) (beta2perp^l d^l_beta2perp /l!) S
// in the CMB frame for a resting observer. Maximal orders in The and betac are used to compute 
// the derivatives.
//
// constraints: dThe<=4; dbeta_para<=2; dbeta2_perp<=1;
//==================================================================================================
void Dcompute_SZ_distortion_asymptotic_CMB(double x, 
                                           int dThe, int dbeta_para, int dbeta2_perp,
                                           double The, double betac, double muc,
                                           vector<double> &dDn_dThe)
{
    dDn_dThe.resize(dThe+1);
    for(int k=0; k<dThe+1; k++) dDn_dThe[k]=0.0;
    
    if(dThe>4) return;
    if(dbeta_para>2) return;
    if(dbeta2_perp>1) return;
    if(dbeta2_perp==1 && dbeta_para!=0) return;
    
    double betac_perp=betac*sqrt(1.0-muc*muc);
    double betac_para=betac*muc;
    
    if(dbeta2_perp==1)
    {
        compute_all_Te_derivatives_upto_dThe(x, The, betac_para, betac_perp, dThe, 
                                             Dn_dbeta2_perp_asymptotic_CMB,
                                             dDn_dThe);
    }
    
    else if(dbeta_para==1)
    {
        compute_all_Te_derivatives_upto_dThe(x, The, betac_para, betac_perp, dThe, 
                                             Dn_dbeta_para_asymptotic_CMB,
                                             dDn_dThe);  
    }

    else if(dbeta_para==2)
    {
        compute_all_Te_derivatives_upto_dThe(x, The, betac_para, betac_perp, dThe, 
                                             Dn_d2beta_para_asymptotic_CMB,
                                             dDn_dThe);
    }

    //==============================================================================================
    // The derivatives only
    //==============================================================================================
    else
    {
        compute_all_Te_derivatives_upto_dThe(x, The, betac_para, betac_perp, dThe, 
                                             Dn_asymptotic_CMB_for_The,
                                             dDn_dThe);        
    }
    
    return;
}

//==================================================================================================
//==================================================================================================
