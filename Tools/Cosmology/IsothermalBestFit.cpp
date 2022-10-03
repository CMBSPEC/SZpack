//==================================================================================================
//
// routine to find best-fit solutions using isothermal model
//
//==================================================================================================
//
// Author: Elizabeth Lee
// Based off work by Jens Chluba
//
//==================================================================================================

#include "IsothermalBestFit.h"

//--------------------------------------------------------------------------------------------------
double f_chisq_betac(const gsl_vector *v, void *params)
{
    Parameters temp = Parameters();
    temp.copyParameters(parameters);

    temp.Dtau = gsl_vector_get(v, 0);
    temp.Te = gsl_vector_get(v, 1);
    temp.betac=gsl_vector_get(v, 2);
    temp.muc=1.0, temp.betao=0.0, temp.muo=0.0;
    
    minimizer_Data *p = (minimizer_Data *)params;
    temp.xcmb = p->xa;
    temp.gridpoints = temp.xcmb.size();
    temp.setCalcValues();
    
    double residual=0.0;
    for(int k=0; k < temp.gridpoints; k++)
    {
        p->Iapprox[k]=pow(temp.xcmb[k], 3)*compute_signal_combo(k, temp);
        
        residual+=pow(p->Iapprox[k]-p->Ia[k], 2);
    }
    //cout << residual << " " << temp.Dtau << " " << temp.Te << endl;
    
    return residual;
}

double f_chisq(const gsl_vector *v, void *params)
{
    Parameters temp = Parameters();
    temp.copyParameters(parameters);

    temp.betac=0.0, temp.muc=0.0, temp.betao=0.0, temp.muo=0.0;
    
    temp.Dtau = gsl_vector_get(v, 0);
    temp.Te = gsl_vector_get(v, 1);
    
    minimizer_Data p = *(minimizer_Data *)params;
    temp.xcmb = p.xa;
    temp.gridpoints = temp.xcmb.size();
    temp.setCalcValues();
    
    double residual=0.0;
    for(int k=0; k < temp.gridpoints; k++)
    {
        p.Iapprox[k]=pow(temp.xcmb[k], 3)*compute_signal_combo(k, temp);
        
        residual+=pow(p.Iapprox[k]-p.Ia[k], 2);
    }
    //cout << residual << " " << temp.Dtau << " " << temp.Te << endl;
    
    return residual;
}
//--------------------------------------------------------------------------------------------------
void find_best_fit_isothermal_model(vector<double> &xa, vector<double> &Ia,
                                   vector<double> &solution, vector<double> &Isol, 
                                   int numpar)
{
    minimizer_Data MinData;
    
    MinData.xa = xa;
    MinData.Ia = Ia;
    MinData.Iapprox.resize(Ia.size());
    
    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
    gsl_multimin_fminimizer *s = NULL;
    gsl_vector *ss, *x;
    gsl_multimin_function minex_func;
    
    size_t iter = 0;
    int status;
    double size;
    
    /* Starting point */
    x = gsl_vector_alloc (numpar);
    gsl_vector_set (x, 0, 0.01);
    gsl_vector_set (x, 1, 5.0);
    if(numpar>2) gsl_vector_set (x, 2, 0.01);
    
    /* Set initial step sizes to 1 */
    ss = gsl_vector_alloc (numpar);
    gsl_vector_set (ss, 0, 0.001);
    gsl_vector_set (ss, 1, 0.5);
    if(numpar>2) gsl_vector_set (ss, 2, 0.001);

    /* Initialize method and iterate */
    minex_func.n = numpar;
    
    if(numpar==2) minex_func.f = f_chisq;
    else if(numpar==3) minex_func.f = f_chisq_betac;
    else{ exit_error("find_best_fit_isothermal_model :: N/A "); }
    
    minex_func.params = (void *) &MinData;
    
    s = gsl_multimin_fminimizer_alloc (T, minex_func.n);
    gsl_multimin_fminimizer_set (s, &minex_func, x, ss);
    
    do
    {
        iter++;
        status = gsl_multimin_fminimizer_iterate(s);
        
        if(status) break;
        
        size = gsl_multimin_fminimizer_size (s);
        status = gsl_multimin_test_size (size, 1.0e-3);
        
        if (status == GSL_SUCCESS)
        {
            cout << " converged to minimum\n";
        }
    }
    while (status == GSL_CONTINUE && iter < 500);

    solution.clear();
    solution.push_back(gsl_vector_get (s->x, 0));
    solution.push_back(gsl_vector_get (s->x, 1));
    if(numpar>2) solution.push_back(gsl_vector_get (s->x, 2));
    Isol = MinData.Iapprox;
    
    for(int k = 0; k < (int)minex_func.n; k++){
        cout << "Parameter " << k << " : " << solution[k] << endl; //TODO: make this better output because what?
    }
    
    gsl_vector_free(x);
    gsl_vector_free(ss);
    gsl_multimin_fminimizer_free (s);
}
