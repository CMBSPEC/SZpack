#ifndef GlobalVariables
#define GlobalVariables

#include <string>
#include "Parameters.h"
using namespace std;

Parameters parameters = Parameters(); //Set initial values of parameters

//==================================================================
// change this to use the convention of Nozawa et al. 2006 for 
// the optical depth. Interpretation of the SZ signal is different
// in this case. The convention of Chluba et al., 2012 is used  
// otherwise, which allow a clean separation of kinematic and 
// scattering effect.
//==================================================================
static bool UseNozawaConvention = false;

//==================================================================================================
//
// As discussed by Chluba et el. 2012, Nozawa 1998 and 2006 used another convention for the optical
// depth variable.
// 
//      Dtau --> gammac (1-betac muc) Dtau 
//
// for an observer at rest in the CMB frame. To use this convention call 'setConvention(true);' 
// before executing the other routines. To reset to the convention of Chluba et al. call 
// 'setConvention(false);', which is the default. These functions only have to be called once. 
// Subsequent calls of routines will use the convention that was set last.
//
//==================================================================================================
#endif
