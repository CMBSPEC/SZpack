//======================================================================================
// Author: Jens Chluba 
// first implementation: Jul 2012
//======================================================================================

#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <vector>

using namespace std;

string int_to_string(int i)
{
    if(i==0) return "0";
    
    char   buf[16];
    int dig=10, num=0;
    string name="";
    
    if(i<0)
    {
        name="-";
        i=-i;
    }
    
    while(i>=dig) dig*=10;
    
    while(dig!=1)
    {
        dig/=10;
        num = (int)(i/dig);
        sprintf(buf, "%d", num); 
        name += *buf;
        i-=num*dig;
    }
    
    return name;
}

//======================================================================================
// load tables of data and convert them into header files
//======================================================================================
void load_data_from_file(string ifname, int cols, string name, ofstream &ofile)
{
    cout << " load_data_from_file :: loading data from file " << ifname << endl;
    
    vector<double> vdum(cols);
    vector<vector<double> > Xi_Data;

    ifstream ifile(ifname.c_str());
    while(ifile)
    {
        for(int i=0; i<cols && ifile; i++) ifile >> vdum[i]; 
        if(ifile) Xi_Data.push_back(vdum);
    }
    ifile.close();
    
    cout << " load_data_from_file :: setting up splines " << endl;
    
    vector<double> xarr(Xi_Data.size());
    vector<double> yarr(Xi_Data.size());
    
    for(int i=0; i<(int)Xi_Data.size(); i++) xarr[i] = log(Xi_Data[i][0]);
    
    ofile << "const double Data_"+name+"[" << int_to_string(Xi_Data.size()) << "][" 
                   << cols << "]={" << endl;
    for(int i=0; i<(int)Xi_Data.size(); i+=1) 
    {   
        ofile << "{" << scientific << xarr[i] << ", ";
        
        for(int c=1; c<cols-1; c++) 
        {
            ofile << scientific << Xi_Data[i][c] << ", ";
            if(!((c+1)%4)) ofile << endl;
        }
        
        ofile << scientific << Xi_Data[i][cols-1] << "}," << endl;
    }
    ofile << "};" << endl << endl;
    
    vdum.clear();
    Xi_Data.clear();
    xarr.clear();
    yarr.clear();
    
    return;
}

int main(int narg, char *args[])
{
    const string SZ_functions_path="./src/outputs/";
    //const string SZ_functions_path="./src/database/";
    //const string SZ_functions_path=SZPACKDIR+"./src/database.optimized/";
    
    const int mmax=20; // maximal temperature order

    string ofname;
    ofstream ofile;

/*
    //string ofname=SZ_functions_path+"SZ_basis.The_0.03.CMB.h";
    //string ofname=SZ_functions_path+"SZ_basis.The_0.1.h";
    string ofname=SZ_functions_path+"SZ_basis.Te_80keV.h";
    ofstream ofile(ofname.c_str());
    ofile.precision(16);
    
//    load_data_from_file(SZ_functions_path+"SZ_basis.3D.The_0.03.monopole.dat", mmax+2, "thSZ_CMB", ofile);
//    load_data_from_file(SZ_functions_path+"SZ_basis.3D.The_0.03.dipole.CMB.dat", mmax+2, "kSZ_CMB", ofile);
//    load_data_from_file(SZ_functions_path+"SZ_basis.3D.The_0.03.quadrupole.CMB.dat", mmax+2, "qkSZ_CMB", ofile);
//    load_data_from_file(SZ_functions_path+"SZ_basis.3D.The_0.03.monopole_corr.CMB.dat", mmax+2, "mkSZ_CMB", ofile);  
//    load_data_from_file(SZ_functions_path+"SZ_basis.3D.The_0.1.monopole.dat", mmax+2, "thSZ_high", ofile);
//    load_data_from_file(SZ_functions_path+"SZ_basis.3D.The_0.1.dipole.dat", mmax+2, "kSZ_high", ofile);
//    load_data_from_file(SZ_functions_path+"SZ_basis.3D.The_0.1.quadrupole.dat", mmax+2, "qkSZ_high", ofile);
//    load_data_from_file(SZ_functions_path+"SZ_basis.3D.The_0.1.monopole_corr.dat", mmax+2, "mkSZ_high", ofile);  
    load_data_from_file(SZ_functions_path+"SZ_basis.3D.Te_80keV.monopole.CMB.dat", mmax+2, "thSZ_c", ofile);
    load_data_from_file(SZ_functions_path+"SZ_basis.3D.Te_80keV.monopole.CMB.dat", mmax+2, "kSZ_c", ofile);
    load_data_from_file(SZ_functions_path+"SZ_basis.3D.Te_80keV.monopole.CMB.dat", mmax+2, "qkSZ_c", ofile);
    load_data_from_file(SZ_functions_path+"SZ_basis.3D.Te_80keV.monopole.CMB.dat", mmax+2, "mkSZ_c", ofile);  
*/    

    //==============================================================================================
    //==============================================================================================
    
    string addb="The_0.01";
    ofname=SZ_functions_path+"SZ_basis."+addb+".CMB.h";
    ofile.open(ofname.c_str());
    ofile.precision(16);
    
    load_data_from_file(SZ_functions_path+"SZ_basis.3D."+addb+".monopole.CMB.dat", mmax+2, "thSZ_I", ofile);
    load_data_from_file(SZ_functions_path+"SZ_basis.3D."+addb+".dipole.CMB.dat", mmax+2, "kSZ_I", ofile);
    load_data_from_file(SZ_functions_path+"SZ_basis.3D."+addb+".quadrupole.CMB.dat", mmax+2, "qkSZ_I", ofile);
    load_data_from_file(SZ_functions_path+"SZ_basis.3D."+addb+".monopole_corr.CMB.dat", mmax+2, "mkSZ_I", ofile);  

    ofile.close();

    return 0;
 
    //==============================================================================================
    //==============================================================================================

    int nT=14;
    string add[14]={"Te_9.4keV", "Te_11keV", "Te_14keV", "Te_18.5keV", "Te_22keV", 
                    "Te_25keV", "Te_30keV", "Te_35keV", "Te_40keV", "Te_45keV",
                    "Te_50keV", "Te_55keV", "Te_60keV", "Te_80keV"};
    
    for(int k=0; k<nT; k++)
    {
        ofname=SZ_functions_path+"SZ_basis."+add[k]+".CMB.h";
        ofile.open(ofname.c_str());
        ofile.precision(16);
        
        string ed=int_to_string(k);
        
        load_data_from_file(SZ_functions_path+"SZ_basis.3D."+add[k]+".monopole.CMB.dat", mmax+2, "thSZ_"+ed, ofile);
        load_data_from_file(SZ_functions_path+"SZ_basis.3D."+add[k]+".dipole.CMB.dat", mmax+2, "kSZ_"+ed, ofile);
        load_data_from_file(SZ_functions_path+"SZ_basis.3D."+add[k]+".quadrupole.CMB.dat", mmax+2, "qkSZ_"+ed, ofile);
        load_data_from_file(SZ_functions_path+"SZ_basis.3D."+add[k]+".monopole_corr.CMB.dat", mmax+2, "mkSZ_"+ed, ofile);  
        
        ofile.close();
    }    
    
    return 0;
}

