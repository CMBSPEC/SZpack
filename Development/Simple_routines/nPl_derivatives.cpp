//==================================================================================================
//
// Author: Jens Chluba 
// first implementation: April 2012
// last modification   : July  2012
//
//==================================================================================================

#include "physical_consts.h"
#include "routines.h"

//==================================================================================================
// namespaces
//==================================================================================================
using namespace std;

//==================================================================================================
vector<vector<double> > E_coeffies;
vector<vector<double> > Es_coeffies;

//==================================================================================================
//
// derivatives of planck spectrum
//
//==================================================================================================
double Binomial_coeff(int n, int k)
{
    double a=log10factorial(n)-log10factorial(k)-log10factorial(n-k);
    return pow(10.0, a);
}

double Eulerian_number(int n, int k)
{
    double r=0.0;
    
    for(int s=0; s<=k; s++)
        r+=pow(-1.0, s)*Binomial_coeff(n+1, s)*pow(1.0+k-s, n);
    
    return r;
}

//==================================================================================================
void update_Eulerian_numbers(int k)
{
    int kmax=E_coeffies.size();
    vector<double> nums;
    
    for(int l=kmax; l<=k; l++)
    {
        nums.clear();
        nums.push_back(1.0);
        
        for(int m=1; m<l-1; m++)
        {
            double dum=(1.0+m)*E_coeffies[l-1][m]+(l-m)*E_coeffies[l-1][m-1];
            nums.push_back(dum);
        }
        
        nums.push_back(1.0);
        
        E_coeffies.push_back(nums);
    }
    
    return;
}

//==================================================================================================
double Pfunc(int k, double x)
{
    if(k+1>(int)E_coeffies.size()) update_Eulerian_numbers(k);
    
    if(k==0) return 1.0;
    
    double exp_mx=exp(-x), P=0.0;
    
    for(int m=k-1; m>0; m--)
    {
        P+=E_coeffies[k][m];
        P*=exp_mx;
    }
    
    P+=E_coeffies[k][0];
    
    return P;        
}

//==================================================================================================
double one_minus_exp_mx(double x, double exp_mx)
{
    //==============================================================================================
    // for small x use series expansion for 1-exp(-x)
    //==============================================================================================
    if(fabs(x)<=0.04) 
    {
        double r=1.0+x/8;
        for(int k=7; k>1; k--) r=pow(-1.0, k)+r*x/k;
        return r*x;
    }
    else return 1.0-exp_mx;
}

double one_minus_exp_mx(double x){ return one_minus_exp_mx(x, exp(-x)); }

//==================================================================================================
double xk_dk_nPl(int k, double x)
{
    double exp_mx=exp(-x);
    
    if(k==0) return exp_mx/one_minus_exp_mx(x, exp_mx);
    
    double Hk=pow(-x/one_minus_exp_mx(x, exp_mx), k)*exp_mx/one_minus_exp_mx(x, exp_mx);
    double P=Pfunc(k, x);
    
    return Hk*P;        
}

//==================================================================================================
void update_s_coeffies(int k)
{
    int kmax=Es_coeffies.size();
    vector<double> nums;
    
    if(kmax==0)
    {
        nums.clear();
        nums.push_back(0.0);
        Es_coeffies.push_back(nums);
        nums.push_back(1.0);
        Es_coeffies.push_back(nums);
        kmax=2;
    }
    
    for(int l=kmax; l<=k; l++)
    {
        nums.clear();
        nums.push_back(0.0);
        nums.push_back(1.0);
        
        for(int m=2; m<l; m++)
        {
            double dum=Es_coeffies[l-1][m-1]+m*Es_coeffies[l-1][m];
            nums.push_back(dum);
        }
        
        nums.push_back(1.0);
        
        Es_coeffies.push_back(nums);
    }
    
    return;
}

//==================================================================================================
void dk_nPl_dks(int kmax, double x, vector<double> &results)
{
    if(kmax+1>(int)Es_coeffies.size()) update_s_coeffies(kmax);
    
    vector<double> derivs(kmax+1);
    for(int k=1; k<=kmax; k++)
        derivs[k]=xk_dk_nPl(k, x);
    
    results.resize(kmax+1);
    
    for(int k=1; k<=kmax; k++)
    {
        results[k]=derivs[1];
        for(int m=2; m<=k; m++) results[k]+=Es_coeffies[k][m]*derivs[m];
    }
    
    return;
}

//==================================================================================================
//==================================================================================================
