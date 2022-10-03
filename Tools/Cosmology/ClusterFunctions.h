//==================================================================================================
//
// Handling the passage between cluster profiles and the way they are used.
//
//==================================================================================================

#ifndef clusterFunctions_H
#define clusterFunctions_H

#include <string>
#include "Parameters.h"
#include "SZpack.h"
#include "SZ_cluster_profiles.h"

using namespace std;

//==================================================================================================
//
// Computation of SZ signal using the simple cluster profile functions of Vikhlinin et al. 2006.
// Here the moments are computed internally by the SZ_moment_method class functions. The profiles
// for Ne and Te are supplied at different slices through the cluster medium.
//
//==================================================================================================
struct profile_slice_params
{
    double xs, ys;
    double zsc;
    double betac, muc;
    double xcmb;
    int k;
    SZ_cluster_profiles Cluster;

    void setValues(double x, double y, double z, double beta, double mu, SZ_cluster_profiles &CL);
};

extern profile_slice_params PSP;

//--------------------------------------------------------------------------------------------------
double betac_func_3D(double lx, double ly, double lz, profile_slice_params psp = PSP);
double muc_func_3D(double lx, double ly, double lz, profile_slice_params psp = PSP);
double Ne_func_3D(double lx, double ly, double lz, profile_slice_params psp = PSP);
double Te_func_3D(double lx, double ly, double lz, profile_slice_params psp = PSP);

//--------------------------------------------------------------------------------------------------
double betac_func(double l, profile_slice_params psp = PSP);
double muc_func(double l, profile_slice_params psp = PSP);
double Ne_func(double l, profile_slice_params psp = PSP);
double Te_func(double l, profile_slice_params psp = PSP);

typedef double (*clusterFun1D)(double, profile_slice_params);
typedef double (*clusterFun3D)(double, double, double, profile_slice_params);

struct clusterFunctions1D{
    clusterFun1D Ne;
    clusterFun1D Te;
    clusterFun1D betac;
    clusterFun1D muc;

    clusterFunctions1D();
};

struct clusterFunctions3D{
    clusterFun3D Ne;
    clusterFun3D Te;
    clusterFun3D betac;
    clusterFun3D muc;

    clusterFunctions3D();
};

static clusterFunctions1D clusterF_1D = clusterFunctions1D();
static clusterFunctions3D clusterF_3D = clusterFunctions3D();


#endif
