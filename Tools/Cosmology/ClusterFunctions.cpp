//==================================================================================================
//
// Handling the passage between cluster profiles and the way they are used.
//
//==================================================================================================

#include "ClusterFunctions.h"

using namespace std;

profile_slice_params PSP;

void profile_slice_params::setValues(double x, double y, double z, double beta, double mu, SZ_cluster_profiles &CL){
    xs = x;
    ys = y;
    zsc = z;
    betac = beta;
    muc = mu;
    Cluster = CL;
}

//--------------------------------------------------------------------------------------------------
double betac_func_3D(double lx, double ly, double lz, profile_slice_params psp){
    return psp.betac;
}

double muc_func_3D(double lx, double ly, double lz, profile_slice_params psp){
    return psp.muc;
}

double Ne_func_3D(double lx, double ly, double lz, profile_slice_params psp){
    return psp.Cluster.Ne(lx*psp.zsc, ly*psp.zsc, lz*psp.zsc);
}

double Te_func_3D(double lx, double ly, double lz, profile_slice_params psp){
    return psp.Cluster.Te(lx*psp.zsc, ly*psp.zsc, lz*psp.zsc);
}

//--------------------------------------------------------------------------------------------------
double betac_func(double l, profile_slice_params psp){ 
    return psp.betac;
}

double muc_func(double l, profile_slice_params psp){
    return psp.muc;
}

double Ne_func(double l, profile_slice_params psp){ 
    return psp.Cluster.Ne(psp.xs, psp.ys, l*psp.zsc);
}

double Te_func(double l, profile_slice_params psp){
    return psp.Cluster.Te(psp.xs, psp.ys, l*psp.zsc);
}

clusterFunctions1D::clusterFunctions1D(){
    Ne = &Ne_func;
    Te = &Te_func;
    betac = &betac_func;
    muc = &muc_func;
}

clusterFunctions3D::clusterFunctions3D(){
    Ne = &Ne_func_3D;
    Te = &Te_func_3D;
    betac = &betac_func_3D;
    muc = &muc_func_3D;
}
